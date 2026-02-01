!#/bin/bash 

## processGnomADcnv.sh

## Download and process gnomAD CNVs BED file, 
# preserving 10 columns and change the coordinate into GenBank style (1-based). 
# https://discuss.gnomad.broadinstitute.org/t/gnomad-copy-number-variants-bed-files/422
wget https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/exome_cnv/gnomad.v4.1.cnv.all.bed

## columns to be preserved. Strand was not included because of being "NA". 
grep chrom gnomad.v4.1.cnv.all.bed | cut -f1-4,13,14,19,22,32,52 > cols
echo "Selected columns of the genomAD CNVs file: "
cat cols | sed "s/#//" 
# chrom	chromStart	chromEnd	name	SVLEN	SVTYPE	Genes	SC	SN	SC_XX

grep -v "#" gnomad.v4.1.cnv.all.bed |\
	cut -f1-4,13,14,19,22,32,52 |\
	# Remove "variant_is_" from the name.
	sed "s/variant_is_//g" |\
	# Some entries have "None" in gene symbol column. MANE gene names can be provided for some, 
	# but they seem to not be included in gnomAD CNVs display, therefore they were removed.
	awk '$7 != "None"' |\
	# Add 1 to start to be in GenBank style.
	awk '$2 = $2 + 1' |\
	# Shift the position of SN and SC to be the same as in SV file.
	awk -v OFS="\t" '{print $1, $2, $3, $4, $5, $6, $7, $9, $10, $8}' \
	> cnv_tmp

echo -e "chrom\tstart\tend\tname\tlength\ttype\tgene\tsn\tsc\txx" > header
cat header cnv_tmp > gnomADcnv

# Clean up.
rm cols header cnv_tmp *.bed.gz
echo ""
echo "Done."
	

