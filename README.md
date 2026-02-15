## Introduction

This is a R shiny app that allows a user to search gnomAD stucture variant data for matching variants using gene name and exon range. To keep it simple while potentially useful, we only include processed gnomAD SVs and CNVs data for deletion and duplication, and limit the transcript to those from MANE and MANE Plus Clinical. Only variants *containing* the query del/dup will be displayed. No option for variants with partial exon change; these variants can be looked up by adjusting the exon range and assessed by comparing the result and exon data. The users can refer to the respective scripts for how the gnomAD data were processed, and reproduce the data by running the scripts. 

## Files
**data** - includes MANE exon, SV, and CNV tables.  
**scripts** - one app, three for processing MANE and gnomAD data.

## App link
On Posit Connect Cloud:  
[gnomADsvSeek](https://019c265e-dfa0-22ff-f9fb-af131e0838e0.share.connect.posit.cloud)  

