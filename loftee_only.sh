#!/bin.bash
source ~/.bash_profile
path=/root/.vep/Plugins/
input=$1
output=$2

/opt/program/ensembl-vep/vep --assembly GRCh38 --offline --vcf \
-i $1 \
-o $2 \
--fork 8 --canonical \
--buffer_size 25000 \
--pick --pick_order canonical \
--force_overwrite \
--plugin LoF,loftee_path:${path},conservation_file:${path}loftee.sql\
,human_ancestor_fa:${path}human_ancestor.fa.gz\
,gerp_bigwig:${path}gerp_conservation_scores.homo_sapiens.GRCh38.bw \
--no_stats
