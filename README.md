# pangenomes_panel_study
Scripts used in "Pangenomes aid accurate detection of large insertion and deletions from gene panel data: the case of cardiomyopathies" by Mazzarotto et al. to post-process GATKhc/Manta/GRAF variant calls, process all samples with ExomeDepth, and merge/harmonize results.


**alignment_HG002.sh**: this script was used to perform the alignment of raw data from Genome in a Bottle sample HG002.

**exomedepth_analysis_HG002.R**: this R script was used to process the HG002 alignment obtained by the script above against a set of 22 whole-exome sequenced unaffected controls.

**exomedepth_analysis.R**: this is an R script that was used to process DCM/HCM samples with ExomeDepth (vs UVOLs) and half of the UVOL samples (n=902) vs the other half (n=903). 

**gatk_graf_manta_output_filtering.sh**: this script was used to post-process variant calls made by Manta, GATK Haplotype Caller and GRAF, and produce final VCF files as well as tab-separated tables with variants called by these 3 callers.

**results_processing_pipeline.R**: this R script was used to merge, harmonize and add extra information to the final variant calls of interest made by the 4 tools. **NOTE:** A minor **bug** in this script causes the assignment of a wrong exon numbering to exons of genes located on the reverse strand. Numbering has been corrected in the submitted Supplementary Tables.
