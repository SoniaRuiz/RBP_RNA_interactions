library(BiocManager)
options(repos = BiocManager::repositories())
library(tidyverse)
library(ggplot2)
library(ggtranscript)
library(DBI)
library(GenomicRanges)
library(rsconnect)

## 1. Declare main functions -----------------------------------------------------------------

# saveRDS(object = MANE_reference %>%
#           as_tibble() %>%
#           dplyr::filter(type %in% c("exon", "CDS", "UTR"),
#                         gene_name %in% gene_names),
#         file = "dependencies/MANE_genes_CLIP_sites.rds")

visualiseCLIP <- function(target_RBP, target_gene, allRBPs) {
  
  # target_RBP <- "TARDBP"
  # target_gene <- "UNC13A"
  
  message("Loading data...")
  
  ## LOAD FILE DEPENDENCIES
  
  if (allRBPs) {
    target_RBP <- CLIP_RBPs 
  }
  
  MANE_local <- MANE %>%
    filter(gene_name == target_gene)
  
  ## GET EXON COORDINATES
  
  MANE_exons <- MANE_local %>% dplyr::filter(type == "exon")
  MANE_CDS <- MANE_local %>% dplyr::filter(type == "CDS")
  MANE_UTR <- MANE_local %>% dplyr::filter(type == "UTR")
  
  
  ## GET RBP BINDING INFO FROM THE GENE OF INTEREST
  
  RBP_CLIP_sites <- introns_CLIP_sites %>%
    filter(RBP %in% target_RBP,
           gene_name == target_gene) %>%
    group_by(iCLIP_seqnames, iCLIP_start, iCLIP_end) %>%
    distinct(RBP, .keep_all = T) %>%
    ungroup() %>%
    dplyr::select(seqnames = iCLIP_seqnames,
                  start = iCLIP_start,
                  end = iCLIP_end,
                  width = iCLIP_width,
                  strand = iCLIP_strand,
                  RBP) %>%
    drop_na()
  
  
  if ( RBP_CLIP_sites %>% nrow() == 0 ) {
    ggplot() +
      theme_void(base_size = 14) +
      geom_text(aes(0,0,
                    label=paste0("No data for the gene and RBP selected."))) %>%
      return()
    
  } else if ( RBP_CLIP_sites %>% nrow() > 100 ) {
      ggplot() +
        theme_void(base_size = 14) +
        geom_text(aes(0,0,
                      label=paste0("Too many RBPs to plot (>100 labels). Please, consider reducing the number of RBPs selected."))) %>%
        return()
      
    
  } else {
    
    message("Plotting data...")
    
    MANE_exons %>%
      dplyr::filter(type == "exon") %>%
      ggplot(aes(
        xstart = start,
        xend = end,
        y = MANE_local$transcript_id %>% unique() 
      )) +
      ggtranscript::geom_intron(
        data = to_intron(MANE_exons, "transcript_name"),
        aes(strand = MANE_exons$strand %>% unique)
      ) +   
      # geom_half_range(
      #   range.orientation = "top",
      #   data = RBP_CLIP_sites
      # ) +
      # geom_half_range(
      #   range.orientation = "top",
      #   data = ref_introns_MSRA,
      #   mapping = aes(height = MSR_A, 
      #                 fill = "MSR Acceptor")
      # ) + 
      geom_range(
        data = MANE_CDS,
        fill = "#999999",
        height = 0.2
      ) +
      geom_range(
        data = MANE_UTR,
        fill = "#333333",
        height = 0.05
      ) +
      ggrepel::geom_label_repel(
        data = RBP_CLIP_sites,
        aes(x = start, 
            label = paste0(RBP),#, " \n ", seqnames, ":", start, "-", end),#),
            color = RBP),
        #nudge_y = 0.5,
        #nudge_x = 3,
        box.padding = 5, 
        #point.padding = unit(x = 0.01, "lines"),
        #segment.color = 'grey50',
        direction = c("both"),
        max.overlaps = 100
      ) +
      theme_light() +
      theme(axis.text.y = element_text(angle = 90, hjust = 0.5),
            axis.text = element_text(size = "14"),
            axis.title = element_text(size = "14"),
            plot.title = element_text(size = "22"),
            plot.subtitle = element_text(size = "18"),
            plot.caption = element_text(size = "11"),
            legend.position = "none",
            legend.text = element_text(size = "13"),
            legend.title = element_text(size = "13")) +
      xlab(paste0("Genomic position (", MANE_exons$seqnames %>% unique() ,")")) + 
      ylab("MANE Transcript") +
      guides(fill = guide_legend(element_blank())) %>%
      return()
  }
  
  
}



tableCLIP <- function(target_RBP, target_gene, allRBPs) {
  
  # target_RBP <- "TARDBP"
  # target_gene <- "GBA"
  
  message("Loading data...")
  
  ## LOAD FILE DEPENDENCIES
  
  if (allRBPs) {
    target_RBP <- CLIP_RBPs 
  }
  
  MANE_local <- MANE %>%
    filter(gene_name == target_gene)
  
  ## GET EXON COORDINATES
  
  MANE_exons <- MANE_local %>% dplyr::filter(type == "exon")
  MANE_CDS <- MANE_local %>% dplyr::filter(type == "CDS")
  MANE_UTR <- MANE_local %>% dplyr::filter(type == "UTR")
  
  
  ## GET RBP BINDING INFO FROM THE GENE OF INTEREST
  
  RBP_CLIP_sites <- introns_CLIP_sites %>%
    filter(RBP %in% target_RBP,
           gene_name == target_gene) %>%
    group_by(iCLIP_seqnames, iCLIP_start, iCLIP_end) %>%
    distinct(RBP, .keep_all = T) %>%
    ungroup() %>%
    dplyr::select(seqnames = iCLIP_seqnames,
                  start = iCLIP_start,
                  end = iCLIP_end,
                  binding_width = iCLIP_width,
                  strand = iCLIP_strand,
                  RBP) %>%
    drop_na()
  
  
  if ( RBP_CLIP_sites %>% nrow() == 0 ) {
    
      return(NULL)
    
  } else {
    
    message("Table data...")
    
    RBP_CLIP_sites %>%
      return()
  }
  
  
}


## 2. Fill in select dropdown objects --------------------------------------------------------


CLIP_RBPs <- readRDS(file = "dependencies/list_CLIP_RBPs.rds")

RBP_choices <- c(CLIP_RBPs) %>% as.list()
names(RBP_choices) <- c(CLIP_RBPs) %>% as.list()


gene_names <- readRDS(file = "dependencies/all_genes_CLIP_sites.rds")
gene_names <- gene_names[!is.na(gene_names)] %>% sort()

gene_choices <- c(gene_names) %>% as.list()
names(gene_choices) <- c(gene_names) %>% as.list()


## 3. Load dependecies ----------------------------------------------------------------------

MANE <- readRDS(file = "dependencies/MANE_genes_CLIP_sites.rds")
introns_CLIP_sites <- readRDS(file = "dependencies/all_introns_CLIP_sites.rds")
