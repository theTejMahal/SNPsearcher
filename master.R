# SNP searcher 
# Tejal Patwardhan Jackson Laboratory Summer Student Program 2016


### SETUP ###

cat("\014")
rm(list = ls())

# Set working directory to source file location 
setwd("~/file/path/master.R")
old.dir <- getwd()

# Set location of SNP file
snp.file = "~/file/path/mgp.v5.merged.snps_all.dbSNP142.vcf.gz"

# Set location to query
chr = 19
start = 14.88e6
end = 16.22e6

# Set strain to query ("A_J", "C57BL_6J", "129S1_SvImJ", "NOD_ShiLtJ", "NZO_HlLt", "CAST_EiJ", "PWK_PhJ", "WSB_EiJ", "C57BL_6J_WSB_EiJ")
strain = "WSB_EiJ"

# Load libraries
library(VariantAnnotation)
library(AnnotationHub)

# Get the Ensembl genes (ver. 80)
hub = AnnotationHub()
hub = query(hub, c("ensembl", "gtf", "mus musculus"))
ensembl = hub[[names(hub)[grep("80", hub$title)]]]

# Get the VCF header.
hdr = scanVcfHeader(snp.file)

# Make a GRanges object for the range you want to query. 
gr = GRanges(seqnames = chr, range = IRanges(start = start, end = end))

# Make a ScanVcfParam.
samples = samples(hdr)[c(5,2,26,28,16,30,35)] # DO founders
param = ScanVcfParam(samples = samples, geno = c("GT", "FI"),
                     which = gr)

# Get the data.
vcf = readVcf(file = snp.file, genome = "mm10", param = param)

# Convert SNPs to allele calls.
vcf = genotypeCodesToNucleotides(vcf)

# Get the SNPs (right now, C57BL/6J is used as a reference).
snps = geno(vcf)$GT

# Add in C57BL/6J.
ref = paste(as.character(fixed(vcf)$REF), as.character(fixed(vcf)$REF), sep = "/")
snps = cbind(snps[,1,drop = FALSE], C57BL_6J = ref,
             snps[,2:7])

### GET PRIVATE SNPS ###

if(strain="A_J"){
# Get A/J-private SNPs
keep = (rowSums(snps[,"A_J"] != snps[,c(2:8)]) == 7)
}

if(strain="C57BL_6J"){
# Get C57BL/6J-private SNPs
keep = (rowSums(snps[,"C57BL_6J"] != snps[,c(1,3:8)]) == 7)
}

if(strain="129S1_SvImJ"){
# Get 129S1/SvImJ-private SNPs
keep = (rowSums(snps[,"129S1_SvImJ"] != snps[,c(1:2,4:8)]) == 7)
}

if(strain="NOD_ShiLtJ"){
# Get NOD_ShiLtJ-private SNPs
keep = (rowSums(snps[,"NOD_ShiLtJ"] != snps[,c(1:3,5:8)]) == 7)
}

if(strain="NZO_HlLt"){
# Get NZO_HlLtJ-private SNPs
keep = (rowSums(snps[,"NZO_HlLtJ"] != snps[,c(1:4,6:8)]) == 7)
}

if(strain="CAST_EiJ"){
# Get CAST_EiJ-private SNPs
keep = (rowSums(snps[,"CAST_EiJ"] != snps[,c(1:5,7:8)]) == 7)
}

if(strain="PWK_PhJ"){
# Get PWK_PhJ-private SNPs
keep = (rowSums(snps[,"PWK_PhJ"] != snps[,c(1:6,8)]) == 7)
}

if(strain="WSB_EiJ"){
# Get WSB-private SNPs
keep = (rowSums(snps[,"WSB_EiJ"] != snps[,c(1:7)]) == 7)
}

if(strain="C57BL_6J_WSB_EiJ"){
# Get SNPs with C57BL/6J & WSB/EiJ having the same base and all the other strains different.
keep = ((snps[,"C57BL_6J"] == snps[,"WSB_EiJ"]) & (rowSums(snps[,"C57BL_6J"] != snps[,c(1,3:7)]) == 6))
}


### GETTING THE SNPS AND CONSEQUENCES ###

# vcf2keep and snps2keep
vcf2keep = vcf[keep]
snps2keep = snps[keep,]

# Get the consequences.
csq = as.list(info(vcf)$CSQ)
csq = lapply(csq, strsplit, split = "\\|")

# Convert each consequence (csq) to a matrix.
wh = which(sapply(csq, length) > 0)
for(i in wh) {
  csq[[i]] = matrix(unlist(csq[[i]]), ncol = length(csq[[i]][[1]]),
                    byrow = TRUE)
}

# Search for splice, stop or missense SNPs: "bad.snps"
bad.snps = sapply(csq, function(z) {
  if(length(z) > 0) {
    z[grep("stop|missense|splice", z[,5]),,drop = FALSE]
  } else {
    ""
  }
})
names(bad.snps) = names(rowRanges(vcf))
bad.snps = bad.snps[sapply(bad.snps, length) > 0]
allbadsnpnames = names(bad.snps)
bad.snps = bad.snps[sapply(bad.snps, length) < 14] 
newbadsnpnames = names(bad.snps)
multiples <- allbadsnpnames[!allbadsnpnames %in% newbadsnpnames]
rm(newbadsnpnames)

# Clean up bad SNPs data
df <- data.frame(cbind(names(bad.snps[1]),matrix(unlist(bad.snps[1]), ncol=13, byrow = T)))
for(i in 2:length(bad.snps)){
  df <- rbind(df,data.frame(cbind(names(bad.snps[i]),matrix(unlist(bad.snps[i], recursive = F), ncol = 13, byrow = T))))
}
invisible(unlist(bad.snps,recursive = FALSE));

cleandf <- data.frame(df$X1,df$X3,df$X6,df$X10,df$X11,df$X4,df$X7,df$X8,df$X9)
colnames(cleandf) <- c("SNP name", "Gene", "Consequence", "Amino Acid Change", "Codon Change", "Transcript", "Position in Transcript", "Position in CDS", "Position in Protein")


### OUTPUTS ###

# Create folder for outputs
if (dir.exists("SNPoutputs") == FALSE) {dir.create("SNPoutputs")}
outputs.dir <- file.path(old.dir, "SNPoutputs")
setwd(outputs.dir)

# Outputs all SNP names to a text file
if (file.exists("allSNPs.txt") == TRUE) {unlink("allSNPs.txt")}
write(names(vcf2keep), file = "allSNPs.txt", append = FALSE, sep = " ")

# Outputs all suspect SNP names to a text file
if (file.exists("suspectSNPs.txt") == TRUE) {unlink("suspectSNPs.txt")}
write(allbadsnpnames, file = "suspectSNPs.txt", append = FALSE, sep = " ")

# Outputs all multiples names to a text file
if (file.exists("multiples.txt") == TRUE) {unlink("multiples.txt")}
write(multiples, file = "multiples.txt", append = FALSE, sep = " ")

# Output Bad SNPs spreadsheet
if (file.exists("suspectSNPs.csv") == TRUE) {unlink("suspectSNPs.csv")}
write.csv(cleandf, "suspectSNPs.csv")

# Go back to old directory
setwd(old.dir)
