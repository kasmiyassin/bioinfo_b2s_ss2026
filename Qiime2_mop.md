

conda activate qiime2-moshpit-2026.1
qiime info

wget -O reads.qza \
    https://polybox.ethz.ch/index.php/s/rgXpDtCMgRgyeKB/download

mosh assembly assemble-megahit \
    --i-reads reads.qza \
    --p-presets meta-sensitive \
    --p-num-cpu-threads 1 \
#    --p-min-contig 500 \
    --o-contigs contigs.qza \
    --verbose