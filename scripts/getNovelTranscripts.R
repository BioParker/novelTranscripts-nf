library(optparse);library(tidyverse)
option_list = list(
  make_option(c("-c", "--ce_vs_ref"), type="character", default=NULL, 
              help="Path output file from -wao bedtools intersect of cassette exon bed with gencode reference gtf", metavar="character"),
  make_option(c("-m", "--models"), type="character", default="NULL", 
              help="transcript models gtf from IsoQuant", metavar="character"),
  make_option(c("-p", "--psi"), type="character", default="NULL",    
              help="transcript level psi file from suppa psiPerIsoform on models gtf", metavar="character"),
  make_option(c("-t", "--tpms"), type="character", default="NULL",    
              help="model transcript tpms from IsoQuant", metavar="character")
)
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

#strand is already handled by the bedtools intersect

ce_vs_gen_wao <- read.delim(opt$ce_vs_ref,header = F)
model_gtf <- rtracklayer::readGFF(opt$models)
colnames(ce_vs_gen_wao) <- c(paste(c("seqid","start","end","transcript_id","score","strand"),"bed",sep = "_"),paste(colnames(model_gtf)[1:9],"gtf",sep = "_"),"overlap")

#unpacking transcripts

ce_vs_gen_wao <- ce_vs_gen_wao %>% separate_longer_delim(transcript_id_bed,",")
ce_vs_gen_wao <- ce_vs_gen_wao %>% mutate(transcript_id_gtf = str_extract(gene_id_gtf,"ENST[[:digit:]]+\\.[[:digit:]]+"),gene_id_gtf = str_extract(gene_id_gtf,"ENSG[[:digit:]]+\\.[[:digit:]]+"))
model_gtf_nc <- model_gtf %>% filter(type == "transcript", transcript_id %in% ce_vs_gen_wao$transcript_id_bed)

#removing decimal + numbers to allow for different annotation versions

ce_vs_gen_wao$gene_id_bed <- gsub("\\.[[:digit:]]+","",model_gtf_nc$gene_id[match(ce_vs_gen_wao$transcript_id_bed, model_gtf_nc$transcript_id)])
ce_vs_gen_wao$gene_id_gtf <- gsub("\\.[[:digit:]]+","",ce_vs_gen_wao$gene_id_gtf)

#get feature lengths and remove likely artifact micro-exons

ce_vs_gen_wao_u <- ce_vs_gen_wao %>% 
                   mutate(len = end_bed-start_bed,
                   len_gen = (end_gtf-start_gtf)+1) %>%
                   filter(len > 10) %>%
                   select(-transcript_id_gtf) %>%
                   unique()


#exon matching based on 5bp windows around model exon length + gene-aware overlap annotation (i.e overlapping exon needs to be from the same gene to qualify)

novel_cassettes <- ce_vs_gen_wao_u %>% 
                        group_by(start_bed,end_bed,transcript_id_bed) %>% 
                        mutate(exon_match = between(start_gtf,start_bed-4,start_bed+6) & between(end_gtf,end_bed-5,end_bed+5)) %>% filter(!any(exon_match)) %>% 
                        mutate(overlapping = case_when(overlap <= 5 | gene_id_bed != gene_id_gtf ~ F,
                                                       overlap > 5 & gene_id_bed == gene_id_gtf ~ T)) %>%
                        ungroup()

#adding in ranking criteria: cassette type, gene type and transcript psi

#summarising to single line per exon by overlap status

unovel_cassettes <- novel_cassettes %>% 
                    group_by(seqid_bed,start_bed,end_bed,strand_bed,transcript_id_bed, gene_id_bed) %>% 
                    summarise(overlapping=any(overlapping)) %>% 
                    mutate(cassette_type = case_when(overlapping == T ~ "overlapping",
                                                         overlapping == F ~ "classical")) %>%
                    select(-overlapping) %>%
                    ungroup()

colnames(unovel_cassettes) <- gsub("_bed","",colnames(unovel_cassettes))

model_gtf_ngenes <- model_gtf %>% 
                    mutate(gene_id = gsub("\\.[[:digit:]]+","",gene_id)) %>% 
                    filter(type == "gene", gene_id %in% unovel_cassettes$gene_id)

unovel_cassettes$gene_type <- model_gtf_ngenes$gene_type[match(unovel_cassettes$gene_id,model_gtf_ngenes$gene_id)]

#generating transcript level table

unovel_transcripts <- unovel_cassettes %>% 
                      group_by(transcript_id,gene_id, gene_type) %>% 
                      summarise(cassette_type = paste(cassette_type,collapse = ";"), n_cassettes = n()) %>%
                      ungroup()

# reading in suppa transcript psi results + appending relevant values to novel transcript table

suppa_tpsi <- read.delim(opt$psi)
colnames(suppa_tpsi) <- c("transcript_id","psi")
suppa_tpsi <- suppa_tpsi %>% separate_wider_delim(cols=transcript_id, delim=";", names=c("gene_id","transcript_id"))
unovel_transcripts$transcript_psi <- suppa_tpsi$psi[match(unovel_transcripts$transcript_id, suppa_tpsi$transcript_id)]

# reading in model transcript tpms and appending relevant values to novel transcript table

tpms <- read.delim(opt$tpms)
colnames(tpms) <- c("transcript_id","tpm")
unovel_transcripts$tpm <- tpms$tpm[match(unovel_transcripts$transcript_id, tpms$transcript_id)]

#writing exon-centric and transcript-centrix outputs

write.table(unovel_transcripts, file = "novelTranscripts.tsv", sep = "\t", quote = F, row.names = F, col.names = T)
write.table(unovel_cassettes, file = "novelCassetteExons.tsv", sep = "\t", quote = F, row.names = F, col.names = T)

