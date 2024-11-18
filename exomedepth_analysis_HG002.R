library("ExomeDepth")

#define reference genome
fasta="/rds/general/user/fmazzaro/home/WORK/Large_variants/resources/reference/Homo_sapiens_assembly38.fasta"

#define target regions
target_bed = read.table("/rds/general/user/fmazzaro/home/WORK/Large_variants/resources/target/hvol.WES.S02972011.CDS.annotate.hg38.bed")
target_bed = target_bed[,1:4]
colnames(target_bed) = c("chromosome","start","end","name")

if(!"ExomeDepth_HG002.RData" %in% list.files("/rds/general/user/fmazzaro/home/WORK/Large_variants/analysis/scripts")){
  
  bam_hg002 = "/rds/general/user/fmazzaro/home/WORK/Large_variants/data/hg002_wes/HG002.bam"
  controlfiles = subset(list.files("/rds/general/user/fmazzaro/ephemeral/HVOLWES_forHG002"),!grepl("bai",list.files("/rds/general/user/fmazzaro/ephemeral/HVOLWES_forHG002")))
  controlfiles = paste0("/rds/general/user/fmazzaro/ephemeral/HVOLWES_forHG002/",controlfiles)
  bamfiles = c(bam_hg002,controlfiles)
  
  my.counts = getBamCounts(bed.frame=target_bed,bam.files=bamfiles,include.chr=FALSE,referenceFasta=fasta)
  
  #select the most appropriate reference set from HVOL for HG002 and call CNVs
  test_sample = as.vector(my.counts[,6])
  controls = as.matrix(my.counts[,7:ncol(my.counts)])
  my.choice = select.reference.set(test.counts=test_sample,reference.counts=controls)
  
  #clean the workspace by keeping only the my.counts object
  rm(list=setdiff(ls(),c("my.choice","my.counts","bamfiles","bamfiles_cases","bamfiles_controls")))
  
  save.image("/rds/general/user/fmazzaro/home/WORK/Large_variants/analysis/scripts/ExomeDepth_HG002.RData")
} else{
    load("/rds/general/user/fmazzaro/home/WORK/Large_variants/analysis/scripts/ExomeDepth_HG002.RData")
}

if(length(my.choice$reference.choice)==1){
  my.matrix = as.matrix(my.counts[,my.choice$reference.choice,drop=FALSE])
} else{
  my.matrix = as.matrix(my.counts[,my.choice$reference.choice])
}
my.reference.selected = apply(X=my.matrix,MAR=1,FUN=sum)

#apply beta-binomial distribution to the full set of exons
all.exons = new("ExomeDepth",test=test_sample,reference=my.reference.selected,formula="cbind(test, reference) ~ 1")

#call CNVs and sort them by confidence level
all.exons = CallCNVs(x=all.exons,transition.probability=10^-4,chromosome=substr(as.character(my.counts$chromosome),4,5),start=my.counts$start,end=my.counts$end,name=my.counts$exon)

#save results
write.table(all.exons@CNV.calls,"/rds/general/user/fmazzaro/home/WORK/Large_variants/analysis/results/HG002_ED_calls.tsv",col.names=T,row.names=F,sep="\t",quote=F)

#create and save BED file for intersection with truth set
bed = all.exons@CNV.calls[,c(7,5,6)]
bed$chromosome = paste0("chr",bed$chromosome)
write.table(bed,"/rds/general/user/fmazzaro/home/WORK/Large_variants/analysis/results/HG002_ED_calls.bed",col.names=F,row.names=F,sep="\t",quote=F)
