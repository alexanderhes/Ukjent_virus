#!/usr/bin/env python3
"""
polish_merge.py -- template bookkeeping for POLISH_SCAFFOLD.

Templates are tab files: name <tab> cased sequence <tab> reference positions
(comma-separated, one token per base; see build_scaffold_template.py:
1..L reference positions, "<p>i" insertions after p, <= 0 / > L overhangs).
Case marks the source: uppercase = contig base, lowercase = reference fill or
read-called base.

Subcommands
  merge      merge one round's reads into the template. Where the reads decide,
             the base is the most frequent read base among the base counts
             (--counts; a deletion if deletions are most frequent, the template
             base wins ties) when at least --min-depth reads cover the position,
             otherwise N. No threshold and no IUPAC codes: mixtures get their
             majority base. Read insertions come from the consensus (--cons,
             samtools consensus -c 0.5 --show-ins --mark-ins: inserted bases
             carried by the majority of the reads) and are kept unless both
             neighbours are kept contig bases. --contig-bases sets where the
             reads decide:
               fallback : everywhere except contig positions with fewer than
                          --min-depth reads, which keep the contig base (default;
                          SPAdes contigs can carry a minority allele);
               keep     : everywhere except contig positions (contig base kept);
               reads    : everywhere (N below --min-depth, also in contigs).
             A template with no reads gets N at every read-decided position, so
             reference fill never reaches the genome.
  finalize   write per-group genomes and statistics:
               {sample}_{group}.fasta          full length (end Ns kept): the
                                               reference length plus overhangs
               {sample}_{group}_trimmed.fasta  leading/trailing Ns removed
               {sample}_genome_stats.tsv       one row per group (see manual)
             and a final tab file (untrimmed, cased) for the positions step.
  counts     samtools mpileup (stdin) -> per-position base counts
             (name pos depth A C G T del; depth = A+C+G+T+del), used by merge
             --counts.
  positions  per genome base: group_id pos depth base refpos in_contig, for the
             plots (depth = reads that decided the base, from the last round's
             base counts; base case = source: uppercase = contig base kept;
             refpos = reference position token; in_contig = 1 on a contig
             position of the unpolished template).
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


def read_tab(path):
    out = {}
    with open(path) as fh:
        for line in fh:
            f = line.rstrip('\n').split('\t')
            if len(f) >= 3:
                out[f[0]] = (f[1], f[2].split(',') if f[2] else [])
    return out


def write_tab(path, tab):
    with open(path, 'w') as fh:
        for name, (seq, pos) in tab.items():
            fh.write(f'{name}\t{seq}\t{",".join(pos)}\n')


def ins_token(tok):
    return tok if tok.endswith('i') else tok + 'i'


def counts(args):
    """samtools mpileup text on stdin -> name pos depth A C G T del."""
    out = sys.stdout
    for line in sys.stdin:
        f = line.rstrip('\n').split('\t')
        if len(f) < 5:
            continue
        ref, s = f[2].upper(), f[4]
        cnt = {'A': 0, 'C': 0, 'G': 0, 'T': 0}
        dl, i, n = 0, 0, len(s)
        while i < n:
            c = s[i]
            if c == '^':                 # read start + mapping quality char
                i += 2
                continue
            if c in '+-':                # indel after this position: skip its bases
                m = re.match(r'[+-](\d+)', s[i:])
                i += len(m.group(0)) + int(m.group(1))
                continue
            if c in '.,':
                if ref in cnt:
                    cnt[ref] += 1
            elif c in '*#':
                dl += 1
            elif c.upper() in cnt:
                cnt[c.upper()] += 1
            i += 1
        out.write(f"{f[0]}\t{f[1]}\t{sum(cnt.values()) + dl}\t{cnt['A']}\t{cnt['C']}\t{cnt['G']}\t{cnt['T']}\t{dl}\n")


def read_counts(path):
    """counts file -> {name: {pos: (depth, {'A':..,'C':..,'G':..,'T':..,'*':..})}}"""
    out = {}
    if path:
        with open(path) as fh:
            for line in fh:
                f = line.rstrip('\n').split('\t')
                if len(f) >= 8:
                    out.setdefault(f[0], {})[int(f[1])] = (
                        int(f[2]), dict(zip('ACGT*', map(int, f[3:8]))))
    return out


def majority(cnt, contig_base):
    """Most frequent read base ('*' = deletion); the contig base wins ties."""
    best = max(cnt.values())
    if cnt.get(contig_base, -1) == best:
        return contig_base
    return next(b for b in 'ACGT*' if cnt[b] == best)


def merge(args):
    """Merge one round's reads into the template. Base calls come from the base
    counts (--counts); the consensus (--cons) only supplies read insertions."""
    tpl = read_tab(args.template)
    cons = read_fasta(args.cons)
    cnts = read_counts(args.counts)
    out = {}
    for name, (seq, pos) in tpl.items():
        d = cnts.get(name, {})
        depth = [d.get(k, (0, None))[0] for k in range(1, len(seq) + 1)]
        # Positions where the contig base is kept: all contig positions (keep),
        # none (reads), or contig positions below the min read depth (fallback).
        # Everywhere else the reads decide.
        is_contig = [b.isupper() and (args.contig_bases == 'keep' or
                                      (args.contig_bases == 'fallback' and dp < args.min_depth))
                     for b, dp in zip(seq, depth)]

        def read_call(k):
            """Read-decided base at template position k (1-based): the most
            frequent read base at >= min depth ('' for a deletion), else N."""
            if depth[k - 1] < args.min_depth:
                return 'n'
            b = majority(d[k][1], seq[k - 1].upper())
            return '' if b == '*' else b.lower()

        if name not in cons:
            # No read maps to this template: kept contig bases, N elsewhere.
            print(f'WARNING: no reads on {name}; read-decided positions set to N', file=sys.stderr)
            out[name] = (''.join(b if keep else 'n' for b, keep in zip(seq, is_contig)), pos)
            continue
        new_seq, new_pos, j, ins = [], [], 0, False
        for c in cons[name]:
            if c in '_+':
                ins = True
                continue
            if ins:
                # Read insertion (majority of the reads, samtools consensus -c 0.5),
                # kept unless both neighbours are kept contig bases.
                ins = False
                left = j >= 1 and is_contig[j - 1]
                right = j < len(seq) and is_contig[j]
                if not (left and right) and c.upper() in 'ACGT':
                    new_seq.append(c.lower())
                    new_pos.append(ins_token(pos[j - 1]) if j >= 1 else ins_token(str(int(pos[0].rstrip('i')) - 1)))
                continue
            j += 1
            if j > len(seq):
                break
            b = seq[j - 1] if is_contig[j - 1] else read_call(j)
            if b:
                new_seq.append(b)
                new_pos.append(pos[j - 1])
        if j != len(seq):
            print(f'WARNING: consensus of {name} has {j} positions, template {len(seq)}', file=sys.stderr)
        out[name] = (''.join(new_seq), new_pos)
    write_tab(args.out, out)


def finalize(args):
    tab = read_tab(args.template)
    depth_sum, depth_n = {}, {}
    with open(args.depth) as fh:
        for line in fh:
            f = line.rstrip('\n').split('\t')
            if len(f) >= 3 and f[2].isdigit() and int(f[2]) >= args.min_depth:
                depth_sum[f[0]] = depth_sum.get(f[0], 0) + int(f[2])
                depth_n[f[0]] = depth_n.get(f[0], 0) + 1

    cols = ['sample', 'safe_taxon', 'segment', 'group_id', 'ref_accession', 'ref_length',
            'n_contigs', 'n_contigs_used', 'contig_bp', 'contig_cov_pct', 'genome_len',
            'n_count', 'N_pct', 'ref_cov_pct', 'mean_depth', 'trimmed_len',
            'n_end_5p', 'n_end_3p', 'n_internal', 'ref_called_bp', 'overhang_5p', 'overhang_3p']
    rows, final = [], {}
    with open(args.groups) as fh:
        header = fh.readline().rstrip('\n').split('\t')
        for line in fh:
            g = dict(zip(header, line.rstrip('\n').split('\t')))
            gid, ref_len = g['group_id'], int(g['ref_length'] or 0)
            seq, pos = tab.get(gid, ('', []))
            up = seq.upper()
            lead = len(up) - len(up.lstrip('N'))
            trail = len(up) - len(up.rstrip('N'))
            trimmed = up[lead:len(up) - trail] if lead < len(up) else ''
            if not trimmed:                      # nothing called: no genome
                up, seq, pos, lead, trail = '', '', [], 0, 0
            called_ref = {int(t) for b, t in zip(up, pos)
                          if b != 'N' and not t.endswith('i') and 1 <= int(t) <= ref_len}
            n_count = up.count('N')
            contig_bp = int(g.get('contig_bp') or 0)
            row = {
                'sample': args.sample, 'safe_taxon': g['safe_taxon'], 'segment': g['segment'],
                'group_id': gid, 'ref_accession': g['ref_accession'], 'ref_length': ref_len,
                'n_contigs': g['n_contigs'], 'n_contigs_used': g.get('n_contigs_used', 0),
                'contig_bp': contig_bp,
                'contig_cov_pct': f'{contig_bp / ref_len * 100:.1f}' if ref_len else 'NA',
                'genome_len': len(up), 'n_count': n_count,
                'N_pct': f'{n_count / len(up) * 100:.1f}' if up else 'NA',
                'ref_cov_pct': f'{min(100, len(called_ref) / ref_len * 100):.1f}' if ref_len else '0.0',
                'mean_depth': f'{depth_sum[gid] / depth_n[gid]:.1f}' if depth_n.get(gid) else 0,
                'trimmed_len': len(trimmed), 'n_end_5p': lead, 'n_end_3p': trail,
                'n_internal': trimmed.count('N'), 'ref_called_bp': len(called_ref),
                'overhang_5p': g.get('overhang_5p', 0) if up else 0,
                'overhang_3p': g.get('overhang_3p', 0) if up else 0,
            }
            rows.append(row)
            if up:
                header_line = f'>{args.sample}|{g["safe_taxon"]}|{g["segment"]}|{g["ref_accession"]}'
                for suffix, s in (('', up), ('_trimmed', trimmed)):
                    with open(f'{args.outdir}/{args.sample}_{gid}{suffix}.fasta', 'w') as out:
                        out.write(header_line + '\n')
                        out.write('\n'.join(s[i:i + 60] for i in range(0, len(s), 60)) + '\n')
                final[gid] = (seq, pos)
    with open(f'{args.outdir}/{args.sample}_genome_stats.tsv', 'w') as fh:
        fh.write('\t'.join(cols) + '\n')
        for r in rows:
            fh.write('\t'.join(str(r[c]) for c in cols) + '\n')
    write_tab(args.final_tab, final)


def keyed_tokens(pos):
    """Position tokens made unique: repeated insertion tokens get an index."""
    seen, out = {}, []
    for t in pos:
        if t.endswith('i'):
            seen[t] = seen.get(t, 0) + 1
            t = f'{t}{seen[t]}'
        out.append(t)
    return out


def positions(args):
    """Per genome base: group_id pos depth base refpos in_contig.
    depth = the read depth that decided the base: the base counts of the last
    polishing round (--counts, on the template of that round, --template-called),
    matched by reference position; a read insertion takes the depth of the base
    before it. in_contig = 1 when the base sits on a contig position of the
    unpolished template (--template0; an inserted base counts when the position
    before it does)."""
    tab = read_tab(args.final_tab)
    contig_pos = {}
    for name, (seq, pos) in read_tab(args.template0).items():
        contig_pos[name] = {t for b, t in zip(seq, pos) if b.isupper() and not t.endswith('i')}
    cnts = read_counts(args.counts)
    depth_at = {}
    for name, (seq, pos) in read_tab(args.template_called).items():
        d = cnts.get(name, {})
        depth_at[name] = {t: d.get(k, (0, None))[0] for k, t in enumerate(keyed_tokens(pos), start=1)}

    with open(args.out, 'w') as out:
        out.write('group_id\tpos\tdepth\tbase\trefpos\tin_contig\n')
        for name, (seq, pos) in tab.items():
            dep, prev = depth_at.get(name, {}), 0
            for k, (b, t, key) in enumerate(zip(seq, pos, keyed_tokens(pos)), start=1):
                prev = dep.get(key, prev)
                ic = int(t.rstrip('i') in contig_pos.get(name, set()))
                out.write(f'{name}\t{k}\t{prev}\t{b}\t{t}\t{ic}\n')


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest='cmd', required=True)
    m = sub.add_parser('merge')
    m.add_argument('--template', required=True)
    m.add_argument('--cons', required=True)
    m.add_argument('--out', required=True)
    m.add_argument('--contig-bases', choices=['keep', 'reads', 'fallback'], default='fallback',
                   help='contig positions: keep the contig base, let the reads decide, '
                        'or let the reads decide where depth >= --min-depth')
    m.add_argument('--counts', required=True, help='per-position base counts of this round (counts subcommand)')
    m.add_argument('--min-depth', type=int, default=5)
    f = sub.add_parser('finalize')
    f.add_argument('--template', required=True)
    f.add_argument('--groups', required=True, help='SCAFFOLD_MINIMAP2 templates TSV')
    f.add_argument('--depth', required=True, help='base counts of the last round (name pos depth ...)')
    f.add_argument('--sample', required=True)
    f.add_argument('--outdir', required=True)
    f.add_argument('--final-tab', required=True)
    f.add_argument('--min-depth', type=int, default=5)
    sub.add_parser('counts')
    p = sub.add_parser('positions')
    p.add_argument('--final-tab', required=True)
    p.add_argument('--template0', required=True, help='template tab before polishing (contig positions)')
    p.add_argument('--counts', required=True, help='base counts of the last polishing round (counts subcommand)')
    p.add_argument('--template-called', required=True, help='template tab of the last polishing round (before its merge)')
    p.add_argument('--out', required=True)
    args = ap.parse_args()
    {'merge': merge, 'finalize': finalize, 'counts': counts, 'positions': positions}[args.cmd](args)


if __name__ == '__main__':
    main()
