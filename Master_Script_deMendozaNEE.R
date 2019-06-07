library(data.table)
library(stringr)
library(ggplot2)
library(dplyr)
library(tidyr)
library(GenomicRanges)
library(R.utils)
library(cowplot)
library(reshape2)

#### Mapping WGBS or TAB-seq

# bs_seeker2-align.py --aligner=bowtie2 --bt2--end-to-end --bt2-p 10 -i sample.fastq -o sample.bam -g Amphimedon_queenslandica.genome.fasta 
# sambamba markdup -r -t 10 -p sample.bam sample.dedup.bam 
# bs_seeker2-call_methylation.py -i sample.dedup.bam -d ~/indexes/Amphimedon_queenslandica.genome.fasta_bowtie2/ -o sample

######################################################################
## functions to read, export and filter genomic ranges objects
######################################################################

read_bed_to_GRobject <- function(bedfile){
        dat <- fread(bedfile)
        gr <- GRanges(seqnames = Rle(dat$V1),
                      ranges = IRanges(start = dat$V2, end = dat$V3),
                      seqlengths = sLengths)
        return(gr)
}

gr_to_bed <- function(gr, out_path){
        dat <- as.data.frame(gr)[ ,1:3]
        dat$type <- gr$type
        write.table(x = dat, file = out_path, quote = FALSE, sep = "\t",
                    row.names = FALSE, col.names = FALSE)
}


filter_10kb <- function(gr){
        seqs <- sLengths[sLengths >= 10000]
        gr <- gr[seqnames(gr) %in% names(seqs)]
        return(gr)
}

filter_out_bacterial_scaffolds <- function (gr, bact_ids = bacterial_scaffolds ) {
        gr <- gr[!(seqnames(gr) %in% bact_ids)]
        return(gr)
}

add_loci <- function(gr){
        
        loci <- str_c(seqnames(gr), start(gr), sep=":") %>%
                str_c(end(gr), sep = "-")
        
        gr$loci <- loci
        
        return(gr)
        
}


######################################################################
## Reading CGmap into R and export RDS with bsseq objects (CG)
######################################################################

# Download genome from ftp://ftp.ensemblgenomes.org/pub/metazoa/release-40/fasta/amphimedon_queenslandica/dna/Amphimedon_queenslandica.Aqu1.dna.toplevel.fa.gz
# decompress and obtain faidx index "samtools faidx Amphimedon_queenslandica.Aqu1.dna.toplevel.fa"
# Import scaffold lengths
genome_fai <- fread(file = "Amphimedon_queenslandica.Aqu1.dna.toplevel.fa.fai")[,c(1,2)]
sLengths <- genome_fai$V2
names(sLengths) <- genome_fai$V1


# Read CGmap files, output from bs_seeker2-call_methylation.py
read_CGmap_into_CpG_granges <- function(CGmap, name, sLengths = sLengths){
        
        # Read the file
        dat <- fread(input = CGmap, sep = "\t", select = c(1,2,3,4,5,7,8),
                     col.names = c("chr", "base", "position", "context","dinucleotide",
                                   "C_reads", "CT_reads"))
        # Subset to CG context only
        message(paste0("processing CG..."))
        datCG <- dat[dat$context == "CG" & !(dat$chr %in% c("chrL","chrP")), ]
        
        # Set up the locus id
        datCG$locus <- NA
        
        # Get a locus id relative to forward strand
        datCG$locus <- ifelse(test = datCG$base == "G",
                            yes = paste0(datCG$chr, ":", datCG$position - 1, "-", datCG$position),
                            no = paste0(datCG$chr, ":", datCG$position, "-", datCG$position +1))
        
        # Drop the unused columns
        datCG <- datCG[ ,c("chr", "base", "position", "context") := NULL]
        
        #merge strands
        message(paste("Collapsing strands"))
        datCG <- datCG[, lapply(.SD, sum), by=.(locus), .SDcols=c("C_reads", "CT_reads")]
        
        #make GRanges object
        message(paste("Making GRanges"))
        gr_CG <- GRanges(datCG$locus, C_reads = datCG$C_reads, CT_reads = datCG$CT_reads)
        
        rm(datCG)
        
        message(paste("Making GRanges"))
        saveRDS(object = gr_CG, file = paste0(name,".CG_granges.rds"))
        
        ########### Export lambda for non-conversion rates
        message(paste0("processing Lambda..."))
        # Subset chrL
        dat_Lambda <- dat[dat$chr == "chrL", ]
        
        #add_sp_name
        dat_Lambda$sample <- name
        
        #save lambda genome
        saveRDS(object = dat_Lambda, file = paste0(name,".lambda.rds"))
        
        
        gc()
        
        # return the aggregated data object
        message(paste0("Objects created for ",name))
        return(gr_CG)
}


###### Download CGmap files from GEO

cleavage_gr <- read_CGmap_into_CpG_granges( CGmap = "GSM3518725_mC_Aqueenslandica_cleavage.CGmap.gz", name = "Aqu_cleavage")
precomp_gr <- read_CGmap_into_CpG_granges( CGmap = "GSM3518726_mC_Aqueenslandica_precompetent_larva.CGmap.gz", name = "Aqu_precompetent")
comp_gr <- read_CGmap_into_CpG_granges( CGmap = "GSM3518727_mC_Aqueenslandica_competent_larva.CGmap.gz", name = "Aqu_competent")
juvenile_gr <- read_CGmap_into_CpG_granges( CGmap = "GSM3518728_mC_Aqueenslandica_juvenile.CGmap.gz", name = "Aqu_juvenile")
adult_gr <- read_CGmap_into_CpG_granges( CGmap = "GSM3518729_mC_Aqueenslandica_adult.CGmap.gz", name = "Aqu_adult")


CGdat <- list(cleavage_gr, precomp_gr, comp_gr, juvenile_gr, adult_gr)
names(CGdat) <- c("Aqu_cleavage", "Aqu_precompetent", "Aqu_competent", "Aqu_juvenile", "Aqu_adult")


######################################################################
## Load all CpG positions in the genome
######################################################################

#first obtain all CG positions in the genome using coverage2cytosine script in Bismark
#then read into R and subset one strand
allCpG <- fread(input = "Amphimedon.allCpG.bismark.gz", sep = "\t", select = c(1,2,3),
                     col.names = c("chr", "position", "strand")) %>% dplyr::filter(strand == "+")
#turn into GRanges
allCpG <- GRanges(seqnames = allCpG$chr, ranges = IRanges(start = allCpG$position, end = allCpG$position+1))
allCpG$CpG <- 1



######################################################################
## Methylation levels per CpG position
######################################################################


bacterial_scaffolds <- fread("bacterial_contigs")

global_mC_filter_out_bacterial <- function(x, contigs_to_filter = bacterial_scaffolds)
{
        
        gr <- CGdat[[x]][!(CGdat[[x]]@seqnames %in% contigs_to_filter)]
        dat <- gr %>% as.data.frame()
        dat$pc <- dat$C_reads / dat$CT_reads
        
        # Number of covered positions
        covPos <- nrow(dat)
        
        # Global mCG level
        global <- sum(dat$C_reads) / sum(dat$CT_reads)
        
        # Summaries
        none <- sum(dat$pc == 0) / covPos
        low <- sum(dat$pc >0 & dat$pc <=0.2) / covPos
        inter <- sum(dat$pc >0.2 & dat$pc <=0.8) / covPos
        high <- sum(dat$pc >0.8) / covPos
        
        # Check
        sum(none, low, inter, high) == 1
        
        return(list(global=global, none=none, low=low, inter=inter, high=high))
}

globalDat <- lapply(1:length(CGdat), global_mC_filter_out_bacterial)
names(globalDat) <- names(CGdat)
globalDat <- data.frame(do.call(rbind, globalDat))
globalDat <- melt(data.matrix(globalDat))
colnames(globalDat) <- c("Time", "Class", "value")
globalDat$Time <- factor(globalDat$Time, levels = names(CGdat))
mC_class <- globalDat[globalDat$Class != "global", ]
global <- globalDat[globalDat$Class == "global", ]


# plot the data
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")


gg <- ggplot(data = mC_class, aes(x = as.factor(Time),y = value, fill = Class)) + 
        geom_bar(position = "fill", stat = "identity") +
        scale_y_continuous() + 
        
        # X axis label
        xlab(label = "") + 
        
        # Y axis label
        ylab(label = "mCG/CG") +
        
        # Add global line
        geom_line(data=global, aes(x = as.numeric(Time), y=value),
                  show.legend = FALSE, inherit.aes = FALSE) +
        
        # Add points for line
        geom_point(data=global, aes(x = as.numeric(Time), y=value),
                   show.legend = FALSE, inherit.aes = FALSE) +
        
        #theme with white background
        theme_bw() +
        
        #eliminates background, gridlines, and chart border
        theme( axis.text.x = element_text(angle = 45, hjust = 1)
               ,plot.background = element_blank()
               ,panel.grid.major = element_blank()
               ,panel.grid.minor = element_blank()
               ,panel.border = element_blank()
        ) +
        
        #draws x and y axis line
        theme(axis.line = element_line(color = 'black')) +
        
        # Edit legend
        scale_fill_manual(values = cbPalette, name="mCG/CG level",
                          labels=c("No (mCG/CG = 0) ",
                                   "Low (mCG/CG >0 and <0.2)",
                                   "Intermediate (mCG/CG >0.2 and <0.8)",
                                   "High (mCG/CG > 0.8)"))

ggsave(plot = gg, filename = "Amphimedon_methylation_levels_plot.pdf")

######################################################################
## Calculate non-conversion rates from lambda
######################################################################

Lambda_to_nonConversion <- function(bsobj){
        a <- readRDS(bsobj)
        name <- str_split(bsobj, pattern = "\\.", simplify = T)[1] 
        b <- a %>% data.frame() %>% group_by(.,dinucleotide) %>% summarise(sum(C_reads),sum(CT_reads),n())
        colnames(b) <- c("context","mC","C","positions")
        b <- mutate(b, mC_C = 100*mC/C)
        
        c <- data.frame(context = "CH", mC = sum(b$mC), C = sum(b$C), positions = sum(b$positions)) %>% mutate(non_conversion_rate = 100*mC/C)
        
        c$sample <- name
        
        return(c)
}

list_of_files <- list.files(pattern = ".lambda.rds$", full.names = TRUE) %>% str_replace(pattern = "\\./", replacement = "")

non_conversion_rates <- lapply(list_of_files, FUN = Lambda_to_nonConversion) %>% do.call(rbind,.) %>% data.frame()

ggplot(non_conversion_rates, aes(x = sample, y = non_conversion_rate)) + geom_bar( stat = "identity") + ylim(c(0,2)) + ylab("Non conversion rate %")

write.table(x = non_conversion_rates, file = "non_conversion_rates.tsv", quote = F, sep = "\t")


######################################################################
## Call Unmethylated Regions
######################################################################

library(MethylSeekR)
#build BSgenome from Amphimedon queenslandica methylome as instructed in http://bioconductor.org/packages/release/bioc/vignettes/BSgenome/inst/doc/BSgenomeForge.pdf
library("BSgenome.Aqueenslandica.JGI.rn4")

#this function takes a filename and a genome id (that has to be loaded previously) and obtains the UMRs and LMRs using MethylSeekR (see: https://bioconductor.org/packages/release/bioc/vignettes/MethylSeekR/inst/doc/MethylSeekR.pdf)
# it can take a long time...
run_UMR_finder <- function(cg_file){
        name_bed <- str_c(CGfile, "_UMRLMR.bed")
        
        meth.gr <- GRanges(seqnames = cg_file@seqnames, ranges = cg_file@ranges@start, `T` = as.numeric(cg_file$CT_reads), M = as.numeric(cg_file$C_reads) )
        
        UMRLMRsegments.gr <- segmentUMRsLMRs(m=meth.gr, meth.cutoff=0.5,
                                             nCpG.cutoff=4,
                                             num.cores=1, myGenomeSeq=Aqueenslandica,
                                             seqLengths=seqlengths(Aqueenslandica))
        saveUMRLMRSegments(segs=UMRLMRsegments.gr,
                           GRangesFilename=name_gr, TableFilename=name_bed)
}

UMRs_by_stage <- lapply(CGdat, FUN = run_UMR_finder) 


#####reading back the "tab" files spitted by the previous function, building GR objects of those .tab files
read_UMR_LMR_tab_to_GRobject <- function(UMRfile){
        dat <- fread(UMRfile)
        gr <- GRanges(seqnames = Rle(dat$chr),
                      ranges = IRanges(start = dat$start, end = dat$end),
                      type = dat$type, nCG.seq = dat$nCG.seq,
                      median.meth = dat$median.meth,
                      seqlengths = sLengths)
        return(gr)
}


list_of_UMR <- list.files(pattern = "UMRsLMRs.tab$",path = "./", full.names = TRUE) %>% str_replace(pattern = "\\./", replacement = "")


#import UMRs and LMRs and collapse
all_LMR_UMRs <- lapply(list_of_UMR, FUN = read_UMR_LMR_tab_to_GRobject) %>% GRangesList() %>% unlist() %>% reduce() 


# obtain mCG on CG levels for UMRs

get_mCG_on_CG_and_coverage <- function(gr, bedgraph_C_CT, allCpGfile = allCpG){
        hits <- findOverlaps(gr, bedgraph_C_CT, ignore.strand = TRUE)
        mcols(gr) <- aggregate(bedgraph_C_CT, hits, C_reads = sum(C_reads), 
                               CT_reads = sum(CT_reads), mCG_CG = sum(C_reads)/sum(CT_reads), 
                               drop = FALSE)
        hits_b <- findOverlaps(gr, allCpGfile, ignore.strand = TRUE)
        ans <- gr
        mcols(ans) <- aggregate(allCpGfile, hits_b, CpG_all=sum(CpG),drop = FALSE)
        gr$CpG_all <- ans$CpG_all
        gr$mean_coverage <- gr$CT_reads/gr$CpG_all
        gr$CpG_density <- (gr$CpG_all*100)/width(gr)
        return(gr)
}

# Filter UMRs for those in scaffolds longer than 10 kb and not in bacterial scaffolds, plus minimal coverage of 2 in the sample with highest coverage and 3 as minimal number of CpGs
all_LMR_UMRs_mCG <- get_mCG_on_CG_and_coverage(gr = all_LMR_UMRs, bedgraph_C_CT = juvenile_gr) 
all_LMR_UMRs_mCG <- all_UMRs_mCG[all_UMRs_mCG$CpG_all >= 3 & all_UMRs_mCG$mean_coverage >= 2] %>% filter_10kb() %>% filter_out_bacterial_scaffolds()

# Import N position stretches and split UMRs/LMRs that encompass a 40 bp long N stretch
Nstretch_gr <- read_bed_to_GRobject("Amphimedon_queenslandica.NNNNstrech.bed")
Nstretch_gr <- Nstretch_gr[width(Nstretch_gr) >= 40]

##get the section of the LMR that does not overlap with the NNNN stretch
all_LMR_UMR_filtered <- GenomicRanges::setdiff(all_LMR_UMRs_mCG,Nstretch_gr)

##export bed file:
gr_to_bed( gr = all_LMR_UMR_filtered, out_path = "Aqu_all_UMRsLMRs_merged.bed")

#### Intersect with features using BEDTools (bedtools intersect) and obtain motifs using HOMER: 
#findMotifsGenome.pl Aqu_all_UMRsLMRs_merged.bed Amphimedon_queenslandica.Aqu1.25.dna.genome.nobacterialscaffolds.fasta HOMER_output -p 10 -size given -len 6,8,10,12


######################################################################
## DAP-seq/ampDAP-seq analysis and plots
######################################################################

#Functions to import peaks and motif scans

read_narrowPeak <- function(bedfile){
        dat <- fread(bedfile)
        gr <- GRanges(seqnames = Rle(dat$V1),
                      ranges = IRanges(start = dat$V2, end = dat$V3),
                      seqlengths = sLengths, strand = "*", 
                      score = dat$V5, name = dat$V4, signalValue = dat$V7,
                      qValue = dat$V9)
        gr <- trim(gr)
        return(gr)
}

read_motifscan_to_GRobject <- function(bedfile){
        dat <- fread(bedfile)
        if (sum(grepl(unique(dat$V6), pattern = "C")) == 1){
                dat$V6 <- ifelse(dat$V6 == "C", yes = "-","+")
        }
        gr <- GRanges(seqnames = Rle(dat$V1),
                      ranges = IRanges(start = dat$V2, end = dat$V3),
                      seqlengths = sLengths, strand = dat$V6, 
                      motif_id = dat$V4, motif_score = dat$V5)
        return(gr)
}
# aligned reads with "bowtie -q --threads 10 -m 1 -v 1 -q -S -1 ${fastq1} -2 ${fastq2} 
#peaks called using: 
# macs2 callpeak -B -t TF_DAPseq.merge.bam -c plX_HALO_empty_DAPseq.merge.bam -n ${i}_B -f BAMPE -q 0.05 --down-sample -g 166679601
# macs2 callpeak -B -t TF_ampDAPseq.merge.bam -c plX_HALO_empty_ampDAPseq.merge.bam -n ${j}_B -f BAMPE -q 0.05 --down-sample -g 166679601

#or download peak files from GEO

#then read peaks into R using
narrowpeak_fls <- list.files(path = "./", pattern = "_peaks.txt$", full.names = TRUE)

narrowpeak_list <- lapply(X = narrowpeak_fls ,FUN = read_narrowPeak) %>% GRangesList()
names(narrowpeak_list) <- narrowpeak_fls %>% str_split(., pattern = "/", simplify = T) %>% .[,ncol(.)] %>% 
        str_remove(., pattern = "\\w+?_plX_HALO_") %>% str_remove(., pattern = "_B_peaks.txt") %>% 
        str_remove(., pattern = "narrowPeaks")

# merge same TF peak files (bedtools merge) and obtain best motif for each TF bed with HOMER:
# findMotifsGenome.pl TF_peak_file.bed Amphimedon_queenslandica.Aqu1.25.dna.genome.nobacteria.fa TF_sample -p 10 -size given -len 6,8,10,12
# then get motif 1 for each TF in homerResults/motif1.motif and scan the genome:
# scanMotifGenomeWide.pl motif1.motif Amphimedon_queenslandica.Aqu1.25.dna.genome.cleannames.fa -p 6 -bed > Aqu_TF_scan.bed

##import Motif scans
motifscans_fls <- list.files(path = "motif_hits/", pattern = "scan.bed$", full.names = TRUE)

motifscans_list <- lapply(X = motifscans_fls ,FUN = read_motifscan_to_GRobject) %>% GRangesList()
names(motifscans_list) <- motifscans_fls %>% str_split(., pattern = "/", simplify = T) %>% .[,ncol(.)] %>% 
        str_remove(., pattern = "Aqu_") %>% str_remove(., pattern = "_DAP_scan.bed")


# function to intersect with mCG data 
intersect_and_sumarise <- function(gr1, gr_CGs = precomp_gr, mincoverage = 10){
        hits <- findOverlaps(gr1, gr_CGs, ignore.strand = TRUE)
        df <- gr1[hits@from] %>% as.data.frame()
        df2 <- gr_CGs[hits@to] %>% as.data.frame() %>% dplyr::rename(seqnames_2 = seqnames, start_2 = start, end_2 = end) %>% dplyr::select(start_2, end_2,C_reads,CT_reads)
        df <- cbind(df, df2)
        df <- df[df$CT_reads >= mincoverage,]
        df <- df %>% mutate(position_within_motif = ifelse(strand == "+", start_2 - start, end - end_2), mCG = C_reads/CT_reads)
        df <- unique(df)
        return(df)
}


compare_DAP_and_ampDAP <- function(TF_id, peaklist = narrowpeak_list, Sample_id, threshold = 0.1){
        
        gr_motif <- motifscans_list[TF_id] %>% unlist() %>% add_loci()
        names(gr_motif) <- 1:length(gr_motif)
        gr_peaks_DAP <- peaklist[[paste0(Sample_id,"_DAPseq")]]
        gr_peaks_ampDAP <- peaklist[[paste0(Sample_id,"_ampDAPseq")]]
        
        gr_peaks_ampDAP_unique <- gr_peaks_ampDAP[!overlapsAny(gr_peaks_ampDAP,gr_peaks_DAP)]
        
        best_motif_peak <- function(gr_peaks){
                motif_vs_peaks_hits <- findOverlaps( gr_motif, gr_peaks, ignore.strand = TRUE)
                df <- gr_motif[motif_vs_peaks_hits@from] %>% as.data.frame() %>% dplyr::select(loci, strand, motif_score) %>% dplyr::rename(strand_o = strand)
                df2 <- gr_peaks[motif_vs_peaks_hits@to]@elementMetadata %>% as.data.frame() %>% dplyr::select(name, signalValue, qValue)
                df <- cbind(df, df2)
                df <- df %>% group_by(name) %>% dplyr::arrange(desc(motif_score)) %>% dplyr::slice(1)
                gr <- GRanges(df$loci, strand = df$strand_o)
                values(gr) <- df
                gr$strand_o <- NULL
                return(gr)
        }
        
        dap_gr <- best_motif_peak(gr_peaks = gr_peaks_DAP)
        ampdap_gr <- best_motif_peak(gr_peaks = gr_peaks_ampDAP_unique)
        
        dap_gr_cgs <- intersect_and_sumarise(gr1 = dap_gr) %>% dplyr::mutate(library = "DAP-seq", sample = Sample_id)
        ampdap_gr_cgs <- intersect_and_sumarise(gr1 = ampdap_gr) %>% dplyr::mutate(library = "ampDAP-seq", sample = Sample_id)
        
        dap_and_ampdap <- rbind(dap_gr_cgs, ampdap_gr_cgs)
        
        
        return(dap_and_ampdap)
}


NRF <- compare_DAP_and_ampDAP(TF_id = "NRF", Sample_id = "NRF" )
NRF_DBD <- compare_DAP_and_ampDAP(TF_id = "NRF", Sample_id = "NRF_DBD" )
YY1 <- compare_DAP_and_ampDAP(TF_id = "YY1", Sample_id = "YY1" )
YY1_DBD <- compare_DAP_and_ampDAP(TF_id = "YY1", Sample_id = "YY1_DBD" )
SP4 <- compare_DAP_and_ampDAP(TF_id = "SP41", Sample_id = "SP41" )
SP4_DBD <- compare_DAP_and_ampDAP(TF_id = "SP41", Sample_id = "SP41_DBD" )
EGR <- compare_DAP_and_ampDAP(TF_id = "EGR", Sample_id = "EGR_DBD" )
GLI <- compare_DAP_and_ampDAP(TF_id = "GLI", Sample_id = "GLI_DBD" )

#plot mCG on CG on all position across the motif

test_CG_proportions_on_positions <- function(df, threshold = 0.05, pval = 0.05){
        
        good_positions <- df$position_within_motif %>% unique()
        
        #test the position hyper and hypo between DAP and ampDAP
        get_enrichment_per_position <- function(x){
                pos <- good_positions[x]
                df_pos <- df[df$position_within_motif == pos,]
                
                df_dap <- df_pos[df_pos$library == "DAP-seq",]
                df_amp <- df_pos[df_pos$library == "ampDAP-seq",]
                
                to_test <- cbind(c(sum(df_dap$mCG <= 0.2), sum(df_dap$mCG >= 0.8)), c(sum(df_amp$mCG <= 0.2), sum(df_amp$mCG >= 0.8)))
                enrichment <- fisher.test(to_test, alternative = "two.sided")
                enrichment <- c(enrichment$p.value, enrichment$estimate)
                return(enrichment)
        }
        enrichments_all <- lapply(X = 1:length(good_positions), FUN = get_enrichment_per_position )
        enrichments_all <- do.call(rbind, enrichments_all) %>% as.data.frame()
        colnames(enrichments_all) <- c("p.value", "odds.ratio")
        enrichments_all$position <- good_positions
        enrichments_all <- enrichments_all %>% dplyr::mutate(significance = ifelse(p.value < pval & odds.ratio > 1, "hypomethylated", 
                                                                                   ifelse(p.value < pval & odds.ratio < 1, "hypermethylated", "NA")),
                                                             test_type = "hypo_vs_hyper")
        
        #test the position C vs T calls between DAP and ampDAP
        get_enrichment_per_position_Cs <- function(x){
                pos <- good_positions[x]
                df_pos <- df[df$position_within_motif == pos,]
                
                df_dap <- df_pos[df_pos$library == "DAP-seq",]
                df_amp <- df_pos[df_pos$library == "ampDAP-seq",]
                
                to_test <- cbind(c(sum(df_dap$CT_reads) - sum(df_dap$C_reads),sum(df_dap$C_reads)), c(sum(df_amp$CT_reads) - sum(df_amp$C_reads), sum(df_amp$C_reads)))
                enrichment <- fisher.test(to_test, alternative = "two.sided")
                enrichment <- c(enrichment$p.value, enrichment$estimate)
                return(enrichment)
        }
        
        enrichments_all_Cs <- lapply(X = 1:length(good_positions), FUN = get_enrichment_per_position_Cs )
        enrichments_all_Cs <- do.call(rbind, enrichments_all_Cs) %>% as.data.frame()
        colnames(enrichments_all_Cs) <- c("p.value_CT", "odds.ratio_CT")
        enrichments_all$position <- good_positions
        enrichments_all_Cs <- enrichments_all_Cs %>% dplyr::mutate(significance_CT = ifelse(p.value_CT < pval & odds.ratio_CT > 1, "hypomethylated", 
                                                                                            ifelse(p.value_CT < pval & odds.ratio_CT < 1, "hypermethylated", "NA")),
                                                                   test_type_CT = "C_vs_T_reads")
        
        
        both <- cbind(enrichments_all, enrichments_all_Cs)
        
        return(both)
}


plot_violins <- function(df, threshold = 0.05){
        #filter CG positions that are rare in the motif
        position_freqs <- df$position_within_motif %>% table() %>% as.data.frame() %>% dplyr::rename(position = '.')
        lower_thres <- sum(position_freqs$Freq)*threshold
        good_positions <- position_freqs$position[position_freqs$Freq >= lower_thres]
        
        #only look at positions with data
        df <- df[df$position %in% good_positions, ]
        
        #set DAP before ampDAP
        df$library <-  factor(df$library, 
                              levels = c("DAP-seq","ampDAP-seq") )
        
        #get the significance values
        significance <- test_CG_proportions_on_positions(df = df ) %>% dplyr::select(position, significance_CT) %>% dplyr::rename(position_within_motif = position)
        
        df <- left_join(df, significance)
        
        df$significance_CT <- factor(df$significance_CT, levels = c("hypomethylated", "hypermethylated", "NA")) 
        
        myColors <- c("#0072B2", "#D55E00", "grey")
        names(myColors) <- c("hypomethylated", "hypermethylated", "NA")
        #colScale <- scale_colour_manual(name = "grp",values = myColors)
        
        
        gg_violin <- ggplot(df, aes(x = as.factor(position_within_motif), y = mCG, fill = significance_CT)) + facet_grid(~library) + ggtitle(label = unique(df$sample)) + 
                theme_bw() + geom_violin() + xlab("position within motif") + theme(legend.position="none") +
                scale_fill_manual(name = "significance_CT", values= myColors)
        # + geom_jitter(alpha = 0.2)
        return(gg_violin)
}


gg_violin_mCG_on_motifs <- plot_grid( plot_violins(NRF), plot_violins(NRF_DBD),  plot_violins(YY1), plot_violins(YY1_DBD),
           plot_violins(SP4), plot_violins(SP4_DBD), plot_violins(EGR), plot_violins(GLI), ncol = 2)
ggsave( gg_violin_mCG_on_motifs, filename = "Per_position_CG_levels_on_motifs.pdf")


# plot C vs T calls on DAP-seq peaks versus ampDAP-seq specific peaks

plot_and_test_CG_proportions <- function(df, threshold = 0.05, TF_id){
        
        #filter CG positions that are rare in the motif
        position_freqs <- df$position_within_motif %>% table() %>% as.data.frame() %>% dplyr::rename(position = '.')
        lower_thres <- sum(position_freqs$Freq)*threshold
        good_positions <- position_freqs$position[position_freqs$Freq >= lower_thres]
        
        #only look at positions with data
        df <- df[df$position %in% good_positions, ]
        
        message(df$library %>% table())
        
        df_dap <- df[df$library == "DAP-seq",]
        df_amp <- df[df$library == "ampDAP-seq",]
        Ts_dap <- sum(df_dap$CT_reads) - sum(df_dap$C_reads)
        Cs_dap <- sum(df_dap$C_reads)
        Ts_amp <- sum(df_amp$CT_reads) - sum(df_amp$C_reads)
        Cs_amp <- sum(df_amp$C_reads)        
        to_test <- cbind(c(Ts_dap, Cs_dap), c(Ts_amp,Cs_amp))
        enrichment <- fisher.test(to_test, alternative = "two.sided")
        enrichment <- c(enrichment$p.value, enrichment$estimate)
        
        to_plot <- data.frame(sample = c(rep("DAP-seq", times = 2), rep("ampDAP-seq", times = 2)),
                              values = c(Ts_dap, Cs_dap, Ts_amp, Cs_amp),
                              Freq = c(Ts_dap/(Ts_dap+Cs_dap),Cs_dap/(Ts_dap+Cs_dap),Ts_amp/(Ts_amp+Cs_amp),Cs_amp/(Ts_amp+Cs_amp)),
                              basecall = c("T","C","T","C"), 
                              TF =  rep(TF_id, times = 4))
        
        
        
        
        
        return(to_plot)
}



proportions_plot <- rbind(plot_and_test_CG_proportions(SP4, TF_id = "SP4"),
                    plot_and_test_CG_proportions(SP4_DBD, TF_id = "SP4_DBD"),
                    plot_and_test_CG_proportions(YY1, TF_id = "YY1"),
                    plot_and_test_CG_proportions(YY1_DBD, TF_id = "YY1_DBD"),
                    plot_and_test_CG_proportions(NRF, TF_id = "NRF"),
                    plot_and_test_CG_proportions(NRF_DBD, TF_id = "NRF_DBD"),
                    plot_and_test_CG_proportions(EGR, TF_id = "EGR"),
                    plot_and_test_CG_proportions(GLI, TF_id = "GLI"))

final_plot$sample <-  factor(final_plot$sample, 
                             levels = c("DAP-seq","ampDAP-seq") )
final_plot$TF <-  factor(final_plot$TF, 
                         levels = c("NRF","NRF_DBD", "YY1", "YY1_DBD","EGR","GLI","SP4","SP4_DBD") )

gg_proportions_plot <- ggplot(data = final_plot, 
             aes(x = sample, y = Freq*100, fill = basecall, label = round(Freq*100,digits = 1))) + 
        geom_bar( stat = "identity" ) +
        geom_text(size = 5, color = "white", position = position_stack(vjust = 0.5), angle = 90) +
        scale_y_continuous() + 
        
        # X axis label
        xlab(label = "") + 
        
        # Y axis label
        ylab(label = "percentage (%)") +
        
        
        #theme with white background
        theme_bw() + 
        #eliminates background, gridlines, and chart border
        theme( axis.text.x = element_text(angle = 90, hjust = 1)
               ,plot.background = element_blank()
               ,panel.grid.major = element_blank()
               ,panel.grid.minor = element_blank()
               ,panel.border = element_blank()
        ) +
        
        #draws x and y axis line
        theme(axis.line = element_line(color = 'black')) +
        
        # Edit legend
        scale_fill_brewer(palette="Paired") +
        facet_grid(~TF)

ggsave(filename = "percentage_CG_on_motif.pdf", plot = gg_proportions_plot, height = 5, width = 6)


######################################################################
## DMR calling and analysis
######################################################################
library(DSS)
library(bsseq)
library(caTools)


## get all sample combinations without repeats
pairwiseCombinations <- function(files){
        
        # Compute number of pairwise combinations
        combinations <- combs(1:length(files), k = 2)
        
        # Function for getting each file pair
        myFunction <- function(x){
                f1 <- combinations[x, 1]
                f2 <- combinations[x, 2]
                file1 <- files[f1]
                file2 <- files[f2]
                return(c(file1, file2))
        }
        
        # Apply function over all combinations
        out <- lapply(X = 1:nrow(combinations), FUN = myFunction)
        
        # Format and return the data
        out <- do.call(what = rbind, out)
        return(data.frame(out))
}

pairwise_comparisons <- pairwiseCombinations(names(CGdat))


# Function for reading and formatting aggregatedCG files into DSS format
gr_to_DSS_format <- function(gr){
        dat <- as.data.frame(gr)
        dat <- dat[,c(1,2,7,6)]
        colnames(dat) <- c("chr", "pos","N","X")
        dat$chr <- as.character(dat$chr)
        dat <- arrange(dat, chr, pos)
        return(dat)
}

CGdat_DSSformat <- lapply(CGdat, gr_to_DSS_format)


### this is very slow, it will take a long while to run
DSS_caller <- function(x){
        
        #get the files in order
        s1_name <- pairwise_comparisons[x,1] %>% as.character()
        s2_name <- pairwise_comparisons[x,2] %>% as.character()
        
        sample_names <- as.vector(as.character(c(s1_name,s2_name)))
        
        s1 <- CGdat_DSSformat[[s1_name]] %>% data.frame()
        s2 <- CGdat_DSSformat[[s2_name]] %>% data.frame()
        
        # Create the BSobject
        BSobj <- makeBSseqData(dat = list(s1, s2), sampleNames = sample_names)
        
        #message(BSobj)
        # Perform the DML test
        dmlTest <- DMLtest(BSobj, group1=s1_name, group2=s2_name, smoothing=TRUE, smoothing.span=100)
        # Call DMRs
        dmrs <- callDMR(dmlTest, delta=0.3, p.threshold=0.05, minCG=5, dis.merge=100)
        
        # Convert dmrs into gr
        gr <- GRanges(seqnames = dmrs$chr,
                      ranges = IRanges(start = as.numeric(dmrs$start),
                                       end = as.numeric(dmrs$end)),
                      nCG=as.numeric(dmrs$nCG),
                      areaStat=as.numeric(dmrs$areaStat),
                      length=as.numeric(dmrs$length)
        )
        add_loci(gr)              
        return(gr)
}

dmrs <- lapply(X = 1:nrow(pairwise_comparisons),FUN = DSS_caller)

# Import DMRs
dmr_fls <- list.files(path = "./", pattern = "DSSoutDMR$", full.names = TRUE)

dmrs <- lapply(X = dmr_fls ,FUN = read_bed_to_GRobject) %>% GRangesList()


#collapse all DMRs
dmrs_collapsed <- dmrs %>% unlist() %>% reduce()


dmr_filter <- function(gr, CGlist, mincoverage = 4, minCpG = 3, min_delta = 0.4, allCpGfile = allCpG){
        
        get_mCG_on_CG_and_filter_by_coverage <- function(x){
                bedgraph_C_CT <- CGlist[[x]]
                hits <- findOverlaps(gr, bedgraph_C_CT, ignore.strand = TRUE)
                mcols(gr) <- aggregate(bedgraph_C_CT, hits, C_reads = sum(C_reads), 
                                       CT_reads = sum(CT_reads), mCG_CG = sum(C_reads)/sum(CT_reads), 
                                       drop = FALSE)
                hits_b <- findOverlaps(gr, allCpGfile, ignore.strand = TRUE)
                ans <- gr
                mcols(ans) <- aggregate(allCpGfile, hits_b, CpG_all=sum(CpG),drop = FALSE)
                gr$CpG_all <- ans$CpG_all
                gr$mean_coverage <- gr$CT_reads/gr$CpG_all
                gr_filtered <- gr[(gr$CpG_all >= minCpG) & (gr$mean_coverage >= mincoverage)]
                gr_filtered <- add_loci(gr_filtered)
                return(gr_filtered)
        }
        # Apply function to all CG files
        message("Calculating CG levels...")
        cg_dat <- lapply(1:length(CGlist), get_mCG_on_CG_and_filter_by_coverage)
        
        # Reduce the data to covered loci in all samples
        message("Subselecting for data covered in all samples...")
        get_all_loci <- function(x){
                loci_list <- cg_dat[[x]]$loci %>% unique()
                return(loci_list)
        }
        covered_in_all <- lapply(X=1:length(cg_dat), FUN = get_all_loci) %>% 
                unlist() %>% table() %>% 
                as.data.frame() %>% filter(Freq == length(cg_dat)) %>%
                dplyr::rename(loci = ".")
        
        reduce_to_covered_in_all <- function(x){
                gr_1 <- cg_dat[[x]]
                gr_1 <- gr_1[gr_1$loci %in% covered_in_all$loci]
                names(gr_1) <- gr_1$loci
                gr_1 <- gr_1[covered_in_all$loci]
                return(gr_1$mCG_CG)
        }
        
        message("Reducing to common loci...")
        df <- lapply(1:length(cg_dat), reduce_to_covered_in_all)
        df <- do.call(cbind, df) %>% as.data.frame()
        rownames(df) <- covered_in_all$loci
        colnames(df) <- names(CGlist)
        max <- apply(df,MARGIN = 1,FUN = max)
        min <- apply(df,MARGIN = 1,FUN = min)
        max_delta <- max - min
        #message of percentatge of discarded DMRs
        total_dmrs <- length(gr)
        to_retain_logical <- max_delta >= min_delta
        retained_dmrs <- sum(to_retain_logical)
        str_c("Of ",total_dmrs," initial DMRs you retain ",retained_dmrs," (",round((retained_dmrs/total_dmrs)*100,digits = 2),"%)") %>% message()
        df <- df[to_retain_logical,]
        return(df)
}

dss_dmrs_filtered <- dmr_filter(gr = dmrs_collapsed, CGlist = CGdat, min_delta = 0.4 ) 

dmrs_in_10kb_scaffolds <- GRanges(rownames(dss_dmrs_filtered)) %>% filter_10kb() %>% filter_out_bacterial_scaffolds() %>% add_loci()

dss_dmrs_filtered <- dss_dmrs_filtered[rownames(dss_dmrs_filtered) %in% dmrs_in_10kb_scaffolds$loci,] 

#cluster using k-means
kmeans_ordered_clusters <- function(df, knum = 5){
        set.seed(20)
        df <- kmeans(df, centers = knum, nstart = 20)
        df <- df$cluster %>% as.data.frame() %>% dplyr::rename(cluster = ".") 
        df$loci <- rownames(df)
        df <- arrange(df, cluster)
        rownames(df) <- df$loci
        return(df)
}

dss_dmrs_filtered_5_clusters <- kmeans_ordered_clusters( df = dss_dmrs_filtered, knum = 5 )

#merge cluster IDs with mCG values

add_cluster_id_to_DF_to_ggplot <- function(df,cluster_id_dic){
        sample_count <- ncol(df)
        df$cluster <- cluster_id_dic[rownames(df),]$cluster %>% as.factor()
        df <- melt(df)
        counts <- table(df$cluster)/sample_count
        counts <- as.data.frame(counts) %>% dplyr::rename(cluster = Var1) %>% melt()
        counts$value <- paste('Cluster ',counts$cluster,' (n = ', counts$value,')', sep = "")
        df$cluster2 <- factor(df$cluster, labels = as.character(counts$value))
        
        return(df)
}

dss_dmrs_filtered_5_clusters <- add_cluster_id_to_DF_to_ggplot(df = dss_dmrs_filtered, cluster_id_dic = dss_dmrs_filtered_5_clusters )

#plot boxplot by DMR cluster

boxplot_DMRs <- ggplot(dss_dmrs_filtered_5_clusters, aes(y=value, x=variable )) +
        #geom_violin(alpha=0.1) +
        geom_boxplot(width=0.8, notch = TRUE, 
                     outlier.size = 1, outlier.shape = 1, outlier.stroke = NA) +
        ylab("mCG/CG") + theme_bw() +
        xlab("stages") +
        theme( axis.text.x = element_text(angle = 45, hjust = 1)) + ylim(c(0,1))
boxplot_DMRs <- boxplot_DMRs + facet_grid(. ~ cluster2) + ggtitle(label = "DSS DMRs ∆0.4 by kmeans clusters")

ggsave(boxplot_DMRs, file = "Boxplot_mCG_DSS_DMRs.pdf")

# obtain closest gene ID 

aqu_genes <- fread("Aqu2.1.genes.bed")
aqu_genes <- GRanges(seqnames = Rle(aqu_genes$V1),
                ranges = IRanges(start = aqu_genes$V2, end = aqu_genes$V3),
                gene_id = aqu_genes$V4, strand = aqu_genes$V6,
                seqlengths = sLengths)

aqu_genes.tss <- resize(aqu_genes, width=1, fix='start')
aqu_promoters <- promoters(aqu_genes, upstream = 1000, downstream = 200) %>% trim()

find_nearest_and_add_metadata <- function (gr1, gr2 = aqu_genes.tss){
        gr1 <- add_loci(gr1)
        hits <- distanceToNearest(gr1, gr2, ignore.strand = TRUE, select = "all")
        gr_out <- gr1[hits@from]
        gr_out@elementMetadata <- cbind(gr_out@elementMetadata,gr2[hits@to]@elementMetadata )
        gr_out@elementMetadata$distance <- hits@elementMetadata$distance
        gr_out <- unique(gr_out)
        return(gr_out)
}

dss_dmrs_filtered_nearest_tss <- find_nearest_and_add_metadata(gr1 = GRanges(rownames(dss_dmrs_filtered)) ) %>% .[.$distance < 10000,]
dss_dmrs_filtered_nearest_promoter <- find_nearest_and_add_metadata(gr1 = GRanges(rownames(dss_dmrs_filtered)), gr2 = aqu_promoters ) %>% .[.$distance == 0,]

dmrs_to_gene_id <- c(dss_dmrs_filtered_nearest_tss,dss_dmrs_filtered_nearest_promoter ) %>% data.frame() %>% dplyr::select(loci, gene_id) %>% unique()


# Load gene expression data
celseq_tpm <- read.table("CELSEQ_Aqu2.summarised.tsv")

#filter non-expressed/detected genes
celseq_tpm <- celseq_tpm[rowSums2(celseq_tpm >= 5) >= 1,]

gene_expr_to_DMR_mC_correlation_matcher <- function(geneExpr, mCG_levels, gene_to_loci){
        #clust counts
        #gene_to_loci <- gene_to_loci %>% as.data.frame() %>% dplyr::select(gene_id, loci)
        total_lines <- nrow(gene_to_loci)
        
        get_corr_gene <- function(x){
                expr_level <- geneExpr[gene_to_loci[x,]$gene_id,] %>% as.vector() %>% as.numeric()
                mCG_level <- mCG_levels[gene_to_loci[x,]$loci,] %>% as.vector() %>% as.numeric()
                corval <- cor(expr_level, mCG_level, method ="pearson")
                return(corval)
        }
        corrvals_all <- lapply(X = 1:total_lines, FUN = get_corr_gene )
        corrvals_all <- do.call(rbind, corrvals_all) %>% as.data.frame()
        colnames(corrvals_all) <- "corr"
        corrvals_all <- cbind(corrvals_all, gene_to_loci)
        return(corrvals_all)
}

gene_dmr_correlations <- gene_expr_to_DMR_mC_correlation_matcher(geneExpr = celseq_tpm, mCG_levels = dss_dmrs_filtered, gene_to_loci = dmrs_to_gene_id)
gene_dmr_correlations$corr %>% hist(, breaks = 20, main = "Correlation distribution between gene TPM and all DMRs mCG/CG levels")

promoter_dmr <- dss_dmrs_filtered_nearest_promoter %>% data.frame() %>% dplyr::select(gene_id, loci) %>% unique()
promoter_dmr_correlations <- gene_expr_to_DMR_mC_correlation_matcher(geneExpr = celseq_tpm, mCG_levels = dss_dmrs_filtered, gene_to_loci = promoter_dmr)
promoter_dmr_correlations$corr %>% hist(, breaks = 20, main = "Correlation distribution between gene TPM and promoter DMR mCG/CG levels")



######################################################################
## hmCG loading
######################################################################

# Read CGmap files, output from bs_seeker2-call_methylation.py
read_hmC_CGmap_into_CpG_granges <- function(CGmap, name, sLengths = sLengths){
        
        # Read the file
        dat <- fread(input = CGmap, sep = "\t", select = c(1,2,3,4,5,7,8),
                     col.names = c("chr", "base", "position", "context","dinucleotide",
                                   "C_reads", "CT_reads"))
        # Subset to CG context only
        message(paste0("processing CG..."))
        datCG <- dat[dat$context == "CG" & !(dat$chr %in% c("chrL","chrP")), ]
        
        # Set up the locus id
        datCG$locus <- NA
        
        # Get a locus id relative to forward strand
        datCG$locus <- ifelse(test = datCG$base == "G",
                              yes = paste0(datCG$chr, ":", datCG$position - 1, "-", datCG$position),
                              no = paste0(datCG$chr, ":", datCG$position, "-", datCG$position +1))
        
        # Drop the unused columns
        datCG <- datCG[ ,c("chr", "base", "position", "context") := NULL]
        
        #merge strands
        message(paste("Collapsing strands"))
        datCG <- datCG[, lapply(.SD, sum), by=.(locus), .SDcols=c("C_reads", "CT_reads")]
        
        #make GRanges object
        message(paste("Making GRanges"))
        gr_CG <- GRanges(datCG$locus, C_reads = datCG$C_reads, CT_reads = datCG$CT_reads)
        
        rm(datCG)
        
        message(paste("Making GRanges"))
        saveRDS(object = gr_CG, file = paste0(name,".CG_granges.rds"))
        
        ########### Export lambda for non-conversion rates
        message(paste0("processing Lambda..."))
        # Subset chrL
        dat_Lambda <- dat[dat$chr == "chrL", ]
        
        #add_sp_name
        dat_Lambda$sample <- name
        
        #save lambda genome
        saveRDS(object = dat_Lambda, file = paste0(name,".lambda.rds"))
        
        ########### Export PUC19 for hmCG protection rates
        message(paste0("processing PUC19..."))
        # Subset chrL
        dat_PUC19 <- dat[dat$chr == "chrP", ]
        
        #add_sp_name
        dat_PUC19$sample <- name
        
        #save lambda genome
        saveRDS(object = dat_PUC19, file = paste0(name,".puc19.rds"))
        
        gc()
        
        # return the aggregated data object
        message(paste0("Objects created for ",name))
        return(gr_CG)
}


precomp_hmC_gr <- read_hmC_CGmap_into_CpG_granges( CGmap = "GSM3518740_hmC_Aqueenslandica_precompetent_larva.CGmap.gz", name = "Aqu_precompetent_TAB")
juvenile_hmC_gr <- read_hmC_CGmap_into_CpG_granges( CGmap = "GSM3518741_hmC_Aqueenslandica_juvenile.CGmap.gz", name = "Aqu_juvenile_TAB")


hmCGdat <- list(precomp_hmC_gr, juvenile_hmC_gr)
names(hmCGdat) <- c("Aqu_precompetent_TAB", "Aqu_juvenile_TAB")

######################################################################
## hmCG global levels
######################################################################



get_mCG_and_hmCG_values_as_percent <- function(gr_hmCG, gr_mCG, TAB_seq_stats, stats_df, contigs_to_filter = bacterial_scaffolds){
        #clean the bad scaffolds
        gr_mCG <- gr_mCG[!seqnames(gr_mCG) %in% contigs_to_filter]
        gr_hmCG <- gr_hmCG[!seqnames(gr_hmCG) %in% contigs_to_filter]
        
        #get the percents
        mCG_percent <- (sum(gr_mCG$C_reads)/sum(gr_mCG$CT_reads))*100 - stats_df["MethylC-seq non-conversion rate (Lambda mC/C)",1]
        C_percent <- 100 - mCG_percent
        hmCG_percent <- sum(gr_hmCG$C_reads)/sum(gr_hmCG$CT_reads)*100
        hmCG_percent_corrected <- 100*(hmCG_percent - stats_df["Non-conversion rate (Lambda mCH/CH)",1] )/(100 - stats_df["Non-conversion rate (Lambda mCH/CH)",1] - (100 - stats_df["Protection rate (PUC19 mCG/CG)",1]))
        mCG_percent <- mCG_percent - hmCG_percent_corrected
        df <- data.frame(context = c("CG","hmCG","mCG"), percent = c( C_percent, hmCG_percent_corrected, mCG_percent))
        return(df)
}

#include the non conversion rates and other stuff from TAB-seq and get the global levels
juvenile_TABseq_stats <- data.frame(row.names = c("Non-oxidation rate (Lambda mCG/CG)","Non-conversion rate (Lambda mCH/CH)", "Protection rate (PUC19 mCG/CG)","MethylC-seq non-conversion rate (Lambda mC/C)" ),
                                    percent = as.numeric(c( "1.5329", "0.372557", "50.6192", "0.366823")))
juvenil_global <- get_mCG_and_hmCG_values_as_percent(gr_hmCG = hmCGdat[["Aqu_juvenile_TAB"]], gr_mCG = CGdat[["Aqu_juvenile"]], stats_df = juvenile_TABseq_stats)
juvenil_global$stage <- "Juvenile"

precompetent_TABseq_stats <- data.frame(row.names = c("Non-oxidation rate (Lambda mCG/CG)","Non-conversion rate (Lambda mCH/CH)", "Protection rate (PUC19 mCG/CG)","MethylC-seq non-conversion rate (Lambda mC/C)" ),
                                        percent = as.numeric(c( "1.87414", "0.396234", "52.3326", "0.434742")))
precompetent_global <- get_mCG_and_hmCG_values_as_percent(gr_hmCG = hmCGdat[["Aqu_precompetent_TAB"]], gr_mCG = CGdat[["Aqu_precompetent"]], stats_df = precompetent_TABseq_stats)
precompetent_global$stage <- "Precompetent Larvae"

greensPalette <- c("#EBEB9C","#9ACC4F","#4C9100")

gg_global_hmCG <- ggplot(data = rbind(juvenil_global,precompetent_global), aes(x = as.factor(stage),y = percent, fill = context, label = round(percent, digits = 1))) + 
        geom_bar(stat = "identity", width = 0.5) +
        scale_y_continuous() + 
        geom_text(size = 3, color = "white", position = position_stack(vjust = 0.5)) +
        # X axis label
        xlab(label = "") + 
        
        # Y axis label
        ylab(label = "Cytosine modification (%)") +
        
        #theme with white background
        theme_bw() + 
        #eliminates background, gridlines, and chart border
        theme( axis.text.x = element_text(angle = 45, hjust = 1)
               ,plot.background = element_blank()
               ,panel.grid.major = element_blank()
               ,panel.grid.minor = element_blank()
               ,panel.border = element_blank()
        ) +
        
        #draws x and y axis line
        theme(axis.line = element_line(color = 'black')) +
        
        # Edit legend
        scale_fill_manual(values = greensPalette, name="CG type")

ggsave( filename = "hmCG_CG_global_levels_corrected.pdf", plot = gg_global_hmCG)

######################################################################
## hmCG sliding windows and motifs
######################################################################


#get sliding windows 150 bp long with a 50 bp overlap

contig_gr <- GRanges(seqnames = names(sLengths), ranges = IRanges(start = 1, end = sLengths), seqlengths = sLengths)

sliding_windows_150bp <- slidingWindows(contig_gr, width = 150, step = 50L) %>% unlist() %>% filter_out_bacterial_scaffolds()


# obtain hmC level for each window, filtering for windows that have at least 3 CpGs and mean coverage of 4

sliding_windows_150bp_hmCG <- dmr_filter(gr = sliding_windows_150bp, CGlist = hmCGdat, min_delta = 0)

#plot scatterplot of correlations
plotAnnotatedScatter <- function(x, y, title = "", x_name = "", y_name = "",
                                 pointCol=colorRampPalette(c('lemonchiffon','lightblue','darkblue')),
                                 legendPos="topleft", legendCex=1, ... ){
        # Generate a linear model summary
        correl <- cor(x,y)
        
        # Format the legend for r and p values
        rp <- vector('expression',1)
        rp[1] <- substitute(expression(italic(r) == valueA),
                            list(valueA = format(correl,dig=3)))[2]
        #rp[2] <- substitute(expression(italic(p) == valueB),
        #                    list(valueB = format(pVal)))[2]
        
        # Plot the data
        smoothScatter(x, y,colramp=pointCol, main=title, xlab=x_name, ylab=y_name,...)
        
        # Add the legend
        legend(legendPos, legend = rp, bty = 'n', cex=legendCex)
}

non_zero_in_both <- sliding_windows_150bp_hmCG$Aqu_juvenile_TAB > 0 & sliding_windows_150bp_hmCG$Aqu_precompetent_TAB > 0

pdf("Aque_hmCG_window_comparison.pdf")
plotAnnotatedScatter(sliding_windows_150bp_hmCG[non_zero_in_both, "Aqu_juvenile_TAB"], sliding_windows_150bp_hmCG[non_zero_in_both,"Aqu_precompetent_TAB"], 
                     title = "hmCG 150 bp sliding windows", x_name = "hmCG/CG\nJuvenile TAB-seq", y_name = "hmCG/CG\nPrecompetent larvae TAB-seq")
abline(h = 0.1, v = 0.1, col = "darkgray", lty = 3)
legend("topright", legend = "~ 1%", bty = 'n', cex=1)
dev.off()

#filter windows for those that have hmCG/CG > 0.1
over01_hmC_juv_gr <- sliding_windows_150bp_hmCG[sliding_windows_150bp_hmCG$Aqu_juvenile_TAB >= 0.1,] %>% rownames() %>% GRanges() %>% reduce() %>% add_loci()
over01_hmC_lar_gr <- sliding_windows_150bp_hmCG[sliding_windows_150bp_hmCG$Aqu_precompetent_TAB >= 0.1,] %>% rownames() %>% GRanges() %>% reduce() %>% add_loci()

#common in both samples
over01_hmC_both_gr <- over01_hmC_juv_gr[overlapsAny(over01_hmC_juv_gr,over01_hmC_lar_gr )]

#export 
gr_to_bed(gr = over01_hmc_both_gr, out_path = "hmC_150windows_both_stages_01.bed")

#find motifs using HOMER: 
#findMotifsGenome.pl hmC_150windows_both_stages_01.bed Amphimedon_queenslandica.Aqu1.25.dna.genome.nobacterialscaffolds.fasta hmC_150windows_HOMER_output -p 10 -size given -len 6,8,10,12

#overlaps with DSS DMRs
hmc_DSS_overlap <- over01_hmC_both_gr[overlapsAny(over01_hmC_both_gr,GRanges(rownames(dss_dmrs_filtered)) )]

#Obtain mCG values on the hmCG rich regions

#get mCG values
over01_hmc_both_mCG <- dmr_filter(gr = over01_hmC_both_gr, CGlist = CGdat, min_delta = 0)
#get hmC values
over01_hmc_both_hmCG <- dmr_filter(gr = over01_hmC_both_gr, CGlist = hmCGdat, min_delta = 0)

#plot both boxplots (mCG and hmCG)
over01_hmc_both_mCG_melt <- melt(over01_hmc_both_mCG)
over01_hmc_both_mCG_melt$variable <- factor(over01_hmc_both_mCG_melt$variable, levels = names(CGdat))
over01_hmc_both_hmCG_melt <- melt(over01_hmc_both_hmCG)
over01_hmc_both_hmCG_melt$variable <- factor(over01_hmc_both_hmCG_melt$variable, levels = names(hmCGdat), labels = c("precomp larvae", "juvenile"))
pmCG <- ggplot(over01_hmc_both_mCG_melt, aes(y=value, x=variable)) +
        #geom_violin(alpha=0.1) +
        geom_boxplot(width=0.8, notch = TRUE, 
                     outlier.size = 1, outlier.shape = 1, outlier.stroke = 0.4) +
        ylab("mCG/CG") + theme_bw() +
        ggtitle(">0.1 hmCG windows") + theme( axis.text.x = element_text(angle = 45, hjust = 1)) +
        xlab("")
phmCG <- ggplot(over01_hmc_both_hmCG_melt, aes(y=value, x=variable)) +
        #geom_violin(alpha=0.1) +
        geom_boxplot(width=0.8, notch = TRUE, 
                     outlier.size = 1, outlier.shape = 1, outlier.stroke = 0.4) +
        ylab("hmCG/CG") + theme_bw() + ylim( c(0,1) ) + ggtitle(">0.1 hmCG windows") +
        xlab("") + theme( axis.text.x = element_text(angle = 45, hjust = 1))



mCG_hmCG_boxplot <- plot_grid(pmCG, phmCG, ncol=2)

ggsave(plot = mCG_hmCG_boxplot, file = "High_hmCG_windows_vs_mCG_boxplot.pdf")

#test for significance against random set of regions
random_windows <- sample(sliding_windows_150bp, size = 10000, replace = FALSE)
random_mCG <- dmr_filter(gr = random_windows, CGlist = CGdat, min_delta = 0)
random_mCG <- melt(random_mCG)
random_mCG$variable <- factor(random_mCG$variable, levels = names(CGdat))

test_all_stages <- function(df1, df2){
        stages <- df1$variable %>% unique()
        tester <- function(x){
                return(wilcox.test(df1$value[df1$variable == x],
                                   df2$value[df2$variable == x], alternative = "less")$p.value)
        }
        pvals <- lapply(X = stages, tester) %>% unlist() %>% p.adjust(., method = "BH")
        names(pvals) <- stages
        return(pvals)
}


test_all_stages(df1 = over01_hmc_both_mCG_melt, df2 = random_mCG)


