#!/bin/bash

## PrecessGnomadSV.sh

## Download genomAD SVs BED file, select 9 columns and del-dup entries.

## Download the BED file.
# https://discuss.gnomad.broadinstitute.org/t/gnomad-structural-variants-bed-file/259 or gnomAD download page 
wget https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/genome_sv/gnomad.v4.1.sv.sites.bed.gz 

## After inspection, select columns to use.
#chrom, start, end, name, SVTYPE, AN, AC, N_HOMALT, AC_XX
echo "Working on the file..."
echo ""
zless gnomad.v4.1.sv.sites.bed.gz | grep -v "#" |\
	cut -f1-4,45,47,48,53,250 |\
	# Remove some redundant words.
	sed "s/gnomAD-SV_v3_//g" |\
	# Del-dup only
	awk '$5 == "DEL"  || $5 == "DUP"' |\
	# Add 1 to start position to be in GenBank style.
	awk -v OFS="\t" '$2 = $2 + 1' \
	> sv_tmp
	
echo  -e "chrom\tstart\tend\tname\ttype\tan\tac\thom\txx" > sv_header
cat sv_header sv_tmp > gnomADsv

## The combined del-dup file is ~88MB, too large for GitHub, so separate the data by chromosome.
mkdir SV
for i in chr{1..22} chrX chrY; do
	awk -v chr="$i" -v OFS="\t" '$1 == chr {print $0}' sv_tmp > sv_"$i"
	cat sv_header sv_"$i" > SV/gnomADsv_"$i"
	rm sv_"$i"
done

## Clean up. 
rm *bed.gz sv_header sv_tmp

echo "Done."

