# SNPsearcher
# Setup code (copy and paste into R console):

download.file("ftp://ftp.jax.org/sanger/current_snps/mgp.v5.merged.snps_all.dbSNP142.vcf.gz", "mgp.v5.merged.snps_all.dbSNP142.vcf.gz", method = "auto", quiet=TRUE)
source("https://bioconductor.org/biocLite.R")
biocLite("VariantAnnotation")
source("https://bioconductor.org/biocLite.R")
biocLite("AnnotationHub")

# Then download R script and run!
