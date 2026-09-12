## Introduction

This is a R shiny app that allows a user to search gnomAD stucture variant data for matching variants using gene name and exon range or an interval. To keep it simple while potentially useful, we only include processed gnomAD SVs and CNVs for deletion and duplication and limit the transcripts to those from MANE Select Plus Clinical and a subset of HGMD transcript not in MANE Select. This program can also be used as an exon lookup table and for displaying the methionie positions of a protein. If the result is of clitical importance, please verify it with offical sources such as UCSC Genome Browser. The users can refer to the respective scripts for how the data were processed, and reproduce the data by running the scripts. 

## Files
**data** - include MANE Select and a subset of HGMD transcripts, summary, SV, and CNV tables.  
**scripts** - one app, others for processing MANE, HGMD and gnomAD data.

## App link
On Posit Connect Cloud:  
[gnomADsvSeek](https://019c265e-dfa0-22ff-f9fb-af131e0838e0.share.connect.posit.cloud)  

