DeepBGC
================================================================================
/usr/local/bin/deepbgc pipeline --score 0.5 --prodigal-meta-mode --merge-max-protein-gap 0 --merge-max-nucl-gap 0 --min-nucl 1 --min-proteins 1 --min-domains 1 --min-bio-domains 0 --classifier-score 0.5 ECO010-megahit.fasta
================================================================================
LOG.txt	Log output of DeepBGC
ECO010-megahit.antismash.json 	AntiSMASH JSON file for sideloading.
ECO010-megahit.bgc.gbk 	Sequences and features of all detected BGCs in GenBank format
ECO010-megahit.bgc.tsv 	Table of detected BGCs and their properties
ECO010-megahit.full.gbk 	Fully annotated input sequence with proteins, Pfam domains (PFAM_domain features) and BGCs (cluster features)
ECO010-megahit.pfam.tsv 	Table of Pfam domains (pfam_id) from given sequence (sequence_id) in genomic order, with BGC detection scores
evaluation/ECO010-megahit.bgc.png 	Detected BGCs plotted by their nucleotide coordinates
evaluation/ECO010-megahit.pr.png 	Precision-Recall curve based on predicted per-Pfam BGC scores
evaluation/ECO010-megahit.roc.png 	ROC curve based on predicted per-Pfam BGC scores
evaluation/ECO010-megahit.score.png 	BGC detection scores of each Pfam domain in genomic order

