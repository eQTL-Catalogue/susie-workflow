suppressPackageStartupMessages(library("SNPRelate"))
suppressPackageStartupMessages(library("optparse"))

option_list <- list(
  optparse::make_option(c("--vcf"), type="character", default=NULL,
                        help="Input VCF file.", metavar = "type"),
  optparse::make_option(c("--gds"), type="character", default=NULL,
                        help="Output GDS file.", metavar = "type"))
opt <- optparse::parse_args(optparse::OptionParser(option_list=option_list))

#Perform conversion
SNPRelate::snpgdsVCF2GDS(opt$vcf, opt$gds, method="biallelic.only")