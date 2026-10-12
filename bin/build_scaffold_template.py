#!/usr/bin/env python3
"""
build_scaffold_template.py

Builds the polishing template of one scaffolding group (SCAFFOLD_MINIMAP2) from
the contig consensus along the reference, and records the reference position of
every template base.

Template case marks the source of each base: uppercase = contig, lowercase =
reference fill (replaced by read calls or N in POLISH_SCAFFOLD).

Reference positions (one token per template base, comma-separated in the .pos
file): 1..L for reference positions, "<p>i" for a contig insertion after
reference position p, <= 0 for 5' overhang bases and > L for 3' overhang bases.

Contig overhangs: contig bases beyond the reference ends (references often lack
UTR ends) are soft-clipped by minimap2 and so missing from the consensus. They
are added as contig bases when
  - the contig's alignment reaches within --end-tol bp of that reference end
    (otherwise the clip is a divergent end, not an overhang),
  - no host BLAST hit (>= --host-min-identity) covers >= 20 bp or half of the
    overhang (chimeric rRNA ends),
  - after removing a terminal homopolymer of >= --min-homopolymer bases (an
    assembler artefact, e.g. a poly-U tract assembled as a k+1 loop), at most
    --max-overhang bp remain (longer ones are rejected as likely chimeras).
The longest passing overhang per end is used.

Usage:
  build_scaffold_template.py --ref ref.fa --cons contig_cons.fa --sam group.sam
      --group GROUP --out-fa GROUP.fa --out-pos GROUP.pos
      [--host-blast host_blast.tsv] [--end-tol 10] [--max-overhang 300]
      [--min-homopolymer 10] [--host-min-identity 90]
Prints to stdout: contig_bp <tab> overhang_5p <tab> overhang_3p
"""

import argparse
import re
import sys


def read_fasta(path):
    seqs, name = {}, None
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if line.startswith('>'):
                name = line[1:].split()[0]
                seqs[name] = []
            elif name is not None:
                seqs[name].append(line)
    return {k: ''.join(v) for k, v in seqs.items()}


def log(msg):
    print(msg, file=sys.stderr)


def walk_consensus(cons, ref):
    """samtools consensus (--show-ins yes --mark-ins --show-del yes) along the
    reference -> (template bases, position tokens, reference bases covered).
    A marker ('_' in samtools 1.23, '+' in its help text) flags the next base as
    a contig insertion; every other character is reference position i:
    N (no contig) -> lowercase reference base, * (contig deletion) -> dropped,
    otherwise the contig base in uppercase."""
    bases, pos, covered, i, ins = [], [], 0, 0, False
    for c in cons.upper():
        if c in '_+':
            ins = True
            continue
        if ins:
            ins = False
            if c != 'N':
                bases.append(c)
                pos.append(f'{i}i')
            continue
        i += 1
        if c == 'N':
            bases.append(ref[i - 1].lower())
            pos.append(str(i))
        else:
            covered += 1
            if c != '*':
                bases.append(c)
                pos.append(str(i))
    if i != len(ref):
        log(f'WARNING: consensus covers {i} reference positions, reference has {len(ref)}')
    return bases, pos, covered


def host_hits(path, min_identity):
    hits = {}
    if not path:
        return hits
    with open(path) as fh:
        for line in fh:
            f = line.rstrip('\n').split('\t')
            if len(f) < 8 or float(f[2]) < min_identity:
                continue
            s, e = sorted((int(f[6]), int(f[7])))
            hits.setdefault(f[0], []).append((s, e))
    return hits


def strip_homopolymer(seq, outer_end, min_run):
    """Remove a run of >= min_run identical bases at the outer end."""
    pattern = r'^(A+|C+|G+|T+)' if outer_end == 'start' else r'(A+|C+|G+|T+)$'
    m = re.search(pattern, seq)
    if m and len(m.group(0)) >= min_run:
        return seq[len(m.group(0)):] if outer_end == 'start' else seq[:-len(m.group(0))]
    return seq


def overhang_candidates(sam_path, ref_len, args, hits):
    """Yield (end, seq, contig) for each passing overhang."""
    with open(sam_path) as fh:
        for line in fh:
            f = line.rstrip('\n').split('\t')
            if len(f) < 11 or f[9] == '*':
                continue
            contig, flag, pos, cigar, seq = f[0], int(f[1]), int(f[3]), f[5], f[9].upper()
            ops = [(int(n), o) for n, o in re.findall(r'(\d+)([MIDNSHP=X])', cigar)]
            aln_ref = sum(n for n, o in ops if o in 'MDN=X')
            end = pos + aln_ref - 1
            hard_l = ops[0][0] if ops[0][1] == 'H' else 0
            hard_r = ops[-1][0] if ops[-1][1] == 'H' else 0
            soft_l = ops[0][0] if ops[0][1] == 'S' else (ops[1][0] if len(ops) > 1 and ops[0][1] == 'H' and ops[1][1] == 'S' else 0)
            soft_r = ops[-1][0] if ops[-1][1] == 'S' else (ops[-2][0] if len(ops) > 1 and ops[-1][1] == 'H' and ops[-2][1] == 'S' else 0)
            qlen = len(seq) + hard_l + hard_r
            minus = bool(flag & 16)

            def contig_span(a, b):  # 0-based [a, b) in alignment orientation (incl. hard clips)
                return (a + 1, b) if not minus else (qlen - b + 1, qlen - a)

            sides = []
            if pos - 1 <= args.end_tol and soft_l > pos - 1:
                n = soft_l - (pos - 1)
                sides.append(('5p', seq[:n], contig_span(hard_l, hard_l + n)))
            if ref_len - end <= args.end_tol and soft_r > ref_len - end:
                n = soft_r - (ref_len - end)
                start = len(seq) - soft_r + (ref_len - end)
                sides.append(('3p', seq[start:], contig_span(hard_l + start, hard_l + len(seq))))
            for side, ov, (cs, ce) in sides:
                overlap = max((min(ce, e) - max(cs, s) + 1 for s, e in hits.get(contig, [])), default=0)
                if overlap >= 20 or overlap >= len(ov) / 2:
                    log(f'overhang {side} of {contig} ({len(ov)} bp) rejected: host hit over {overlap} bp')
                    continue
                trimmed = strip_homopolymer(ov, 'start' if side == '5p' else 'end', args.min_homopolymer)
                trimmed = trimmed.strip('N')
                if len(trimmed) != len(ov):
                    log(f'overhang {side} of {contig}: terminal homopolymer/N removed ({len(ov)} -> {len(trimmed)} bp)')
                if len(trimmed) > args.max_overhang:
                    log(f'overhang {side} of {contig} ({len(trimmed)} bp) rejected: longer than {args.max_overhang} bp')
                    continue
                if trimmed:
                    yield side, trimmed, contig


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--ref', required=True)
    ap.add_argument('--cons', required=True)
    ap.add_argument('--sam', required=True)
    ap.add_argument('--group', required=True)
    ap.add_argument('--out-fa', required=True)
    ap.add_argument('--out-pos', required=True)
    ap.add_argument('--host-blast', default='')
    ap.add_argument('--end-tol', type=int, default=10)
    ap.add_argument('--max-overhang', type=int, default=300)
    ap.add_argument('--min-homopolymer', type=int, default=10)
    ap.add_argument('--host-min-identity', type=float, default=90)
    args = ap.parse_args()

    ref = next(iter(read_fasta(args.ref).values())).upper()
    cons = next(iter(read_fasta(args.cons).values()))
    bases, pos, covered = walk_consensus(cons, ref)

    best = {'5p': ('', None), '3p': ('', None)}
    for side, ov, contig in overhang_candidates(args.sam, len(ref), args, host_hits(args.host_blast, args.host_min_identity)):
        if len(ov) > len(best[side][0]):
            best[side] = (ov, contig)
    ov5, ov3 = best['5p'][0], best['3p'][0]
    for side, (ov, contig) in best.items():
        if ov:
            log(f'{args.group}: {side} overhang of {len(ov)} bp from {contig} added')

    bases = list(ov5) + bases + list(ov3)
    pos = [str(k - len(ov5) + 1) for k in range(len(ov5))] + pos + [str(len(ref) + k + 1) for k in range(len(ov3))]
    with open(args.out_fa, 'w') as fh:
        fh.write(f'>{args.group}\n{"".join(bases)}\n')
    with open(args.out_pos, 'w') as fh:
        fh.write(','.join(pos) + '\n')
    print(f'{covered}\t{len(ov5)}\t{len(ov3)}')


if __name__ == '__main__':
    main()
