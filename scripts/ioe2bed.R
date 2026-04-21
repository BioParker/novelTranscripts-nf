#reconstruct bedfile of cassette exons from suppa ioe file of skipped exon (SE) events
library(optparse);library(tidyverse)
option_list = list(
  make_option(c("-e", "--ioefile"), type="character", default="files/suppa_outputs/transcript_model_as_events_SE_strict.ioe", 
              help="Path to .ioe file output from suppa", metavar="character"))
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)
ioe <- read.delim(opt$ioefile)
colnames(ioe)[colnames(ioe) == "alternative_transcripts"] <- "name"
bed <- ioe %>% 
       mutate(score = 0) %>% 
       separate_wider_delim(event_id,":",names = c("event","chrom","chromStart","chromEnd","strand")) %>% 
       select(chrom, chromStart,chromEnd,name,score,strand) %>% 
       mutate(chromStart = as.numeric(gsub("[[:digit:]]+-","",chromStart))-1, chromEnd=gsub("-[[:digit:]]+","",chromEnd))
write.table(bed, file = "cassetteExon.bed", sep = "\t", quote = F, row.names = F, col.names = F)
