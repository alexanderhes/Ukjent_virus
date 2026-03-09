FROM community.wave.seqera.io/library/bowtie2_esviritu_samtools_seqkit_r-tidyverse:3ee4a52f7d6ae7d9

# Install dataui R package (required for EsViritu sparkline coverage plots)
RUN R -e "remotes::install_github('timelyportfolio/dataui', upgrade = 'never')"
