# Separation of phyloseq objects and generating plots

library(phyloseq)
library(here)
library(microViz)
library(ggplot2)


DNA_fixed <- readRDS("DNA_fixed_ps.rds")
ps_DNA <- DNA_fixed
ps_DNA <-subset_taxa(ps_DNA, order !="Chloroplast")
ps_DNA <-subset_samples(ps_DNA, Names != "IS19_Culture")

ps_DNA_Euk <- subset_taxa(ps_DNA, domain =="Eukaryota") #For comparison to 18S
ps_DNA_BacArch <- subset_taxa(ps_DNA, domain !="Eukaryota") #For comparison to 16S




ps_RNA <- TotalRNA_fixed
ps_RNA <-subset_taxa(ps_RNA, order !="Chloroplast")

ps_RNA_Euk <- subset_taxa(ps_RNA, domain =="Eukaryota") #For comparison to 18S
ps_RNA_BacArch <- subset_taxa(ps_RNA, domain !="Eukaryota") #For comparison to 16S


sample_data(ps_DNA_Euk)$DNA_or_RNA <- "DNA"
sample_data(ps_DNA_BacArch)$DNA_or_RNA <- "DNA"

sample_data(ps_RNA_Euk)$DNA_or_RNA <- "RNA"
sample_data(ps_RNA_BacArch)$DNA_or_RNA <- "RNA"



# David Barnett's code for specific functions to modify the color palettes to use with microViz. This code addresses this issue: https://github.com/david-barnett/microViz/issues/16

# this makes a named palette vector from your combined dataset (considering overall abundance to assign colours)
tax_palette <- function(data, # phyloseq or ps_extra 
                        rank, # e.g. "Genus"
                        n, # n colours / taxa not including other
                        by = sum, # method for tax_sort
                        pal = "brewerPlus", # palette name from distinct_palette
                        add = c(other = "lightgrey"), # name = value pairs appended to end of output
                        ... # other args passed to tax_sort
) {
  taxa <- tax_top(data = data, rank = rank, n = n, by = by, ...)
  taxColours <- distinct_palette(n = n, pal = pal, add = NA)
  
  names(taxColours) <- taxa
  taxColours <- c(taxColours, add)
  return(taxColours)
}

# palette viewer function (unnecessary but maybe handy)
tax_palette_plot <- function(
  named_pal_vec # named vector of colours
) {
  stopifnot(!anyNA(names(named_pal_vec))) # all colours need names
  p <- 
    named_pal_vec %>% 
    as.data.frame() %>% 
    tibble::rownames_to_column("taxon") %>% 
    dplyr::mutate(taxon = factor(taxon, levels = rev(taxon))) %>% 
    dplyr::rename("hex" = ".") %>% 
    ggplot2::ggplot(ggplot2::aes(y = taxon, fill = hex, x = "")) +
    ggplot2::geom_raster() +
    ggplot2::scale_fill_identity() + 
    ggplot2::labs(x = NULL, y = NULL) + 
    ggplot2::theme_minimal()
  return(p)
}


# Generating palette's based on functions above
myPal <- tax_palette(ps_merged, rank = "class", pal = "brewerPlus", n = 40)
myPal %>% tax_palette_plot() # just to check the palette

EukPal <- tax_palette(ps_merged_Euk, rank = "class", pal = "brewerPlus", n = 40)
EukPal %>% tax_palette_plot() # just to check the palette

BacPal <- tax_palette(ps_merged_BacArch, rank = "class", pal = "kelly", n = 20)
BacPal %>% tax_palette_plot()


#####################
##### 18S Plots #####
#####################


# 18S barplot from TotalRNA
ps_RNA_Euk %>%
  comp_barplot(tax_level = "class", 
               n_taxa = 15,
               label = "Preservation",
               sample_order = "default",
               bar_outline_colour = NA, 
               tax_transform_for_plot = "compositional",
               bar_width = 0.9,
               merge_other = TRUE, 
               group_by = "class",
               palette = myPal 
               ) + 
 facet_grid(~factor(Site, levels=c('IS19_10', 'IS19_11', 'IS19_12','IS19_13','IS19_14','Control')), scales = "free", space = "free" # these options are critically important!
  ) +
  scale_y_continuous(expand = c(0.005, 0.005)) + 
  ylab("Relative Abundance") + 
  theme(text = element_text(size = 15), axis.text.x = element_text(size = 8))

# Save the created plot
ggsave(filename = here("plots", "composition_RNA_18S_class.png") ,plot = last_plot(),width = 30,height = 15,units = "cm")


# 18S barplot from DNA
ps_DNA_Euk %>%
  comp_barplot(tax_level = "class", 
               n_taxa = 15,
               label = "Preservation",
               sample_order = "default",
               bar_outline_colour = NA, 
               tax_transform_for_plot = "compositional",
               bar_width = 0.9,
               merge_other = TRUE, 
               group_by = "class",
               palette = myPal 
               ) + 
  facet_grid(~factor(Site, levels=c('IS19_10', 'IS19_11', 'IS19_12','IS19_13','IS19_14','Control')), scales = "free", space = "free" # these options are critically important!
  ) +
  scale_y_continuous(expand = c(0.005, 0.005)) + 
  ylab("Relative Abundance") + 
  theme(text = element_text(size = 15), axis.text.x = element_text(size = 8))

# Save the created plot
ggsave(filename = here("plots", "composition_DNA_18S_class.png") ,plot = last_plot(),width = 30,height = 15,units = "cm")



#####################
##### 16S Plots #####
#####################

# 16S barplot from RNA
ps_RNA_BacArch %>%
  comp_barplot(tax_level = "class", 
               n_taxa = 15,
               label = "Preservation",
               sample_order = "default",
               bar_outline_colour = NA, 
               tax_transform_for_plot = "compositional",
               bar_width = 0.9,
               merge_other = TRUE, 
               group_by = "class",
               palette = BacPal 
               ) + 
  facet_grid(~factor(Site, levels=c('IS19_10', 'IS19_11', 'IS19_12','IS19_13','IS19_14','Control')), scales = "free", space = "free" # these options are critically important!
  ) +
  scale_y_continuous(expand = c(0.005, 0.005)) + 
  ylab("Relative Abundance") + 
  theme(text = element_text(size = 15), axis.text.x = element_text(size = 8))

# Save the created plot
ggsave(filename = here("plots", "composition_RNA_16S_class.png") ,plot = last_plot(),width = 30,height = 15,units = "cm")



# 16S barplot from DNA
ps_DNA_BacArch %>%
  comp_barplot(tax_level = "class", 
               n_taxa = 15,
               label = "Preservation",
               sample_order = "default",
               bar_outline_colour = NA, 
               tax_transform_for_plot = "compositional",
               bar_width = 0.9,
               merge_other = TRUE, 
               group_by = "class",
               palette = BacPal 
               ) + 
  facet_grid(~factor(Site, levels=c('IS19_10', 'IS19_11', 'IS19_12','IS19_13','IS19_14','Control')), scales = "free", space = "free" # these options are critically important!
  ) +
  scale_y_continuous(expand = c(0.005, 0.005)) + 
  ylab("Relative Abundance") + 
  theme(text = element_text(size = 15), axis.text.x = element_text(size = 8))

# Save the created plot
ggsave(filename = here("plots", "composition_DNA_16S_class.png") ,plot = last_plot(),width = 30,height = 15,units = "cm")


#####################
##### RDA Plots #####
#####################


# TotalRNA constrained ordination RDA
ps<-TotalRNA_fixed
ps<-subset_taxa(ps,order !="Chloroplast")
ps<-subset_samples(ps,Sample != "IS19_12")
ps<-scale_reads(ps,100000)
ps <- ps %>%
  ps_mutate(
    Lib.conc_scaled = c(scale(Library.conc, center = TRUE, scale = TRUE)),
    TOC.scaled = c(scale(TOC, center =TRUE, scale = TRUE)),
    F = if_else(Preservation == "F", 1, 0),
    Z = if_else(Preservation == "Z", 1, 0),
    L = if_else(Preservation == "L", 1, 0)
  )
constr_ord <- ps %>%
  tax_filter(min_prevalence = 1/100, tax_level = "family") %>%
  tax_agg("family") %>%
  tax_transform("total") %>%
  ord_calc(method = "RDA", constraints = c("TOC.scaled","Lib.conc_scaled", "F", "L", "Z"))

ord_plot(constr_ord, plot_taxa = 1:8,tax_vec_length = 1.5,constraint_vec_length = 0.8,constraint_lab_style = list(alpha=0.2),tax_lab_style = list(alpha=0.2),size = 3, colour = "Preservation", shape = "Sample") +
  scale_color_viridis(discrete=TRUE)

ggsave(filename = "constrained_family.png",plot = last_plot(),width = 20,height = 15,units = "cm")




# DNA constrained ordination RDA
ps<-DNA_fixed
ps<-subset_taxa(ps,order !="Chloroplast")
ps<- subset_samples(ps, Names != "IS19_Culture" & Site != "IS19_12")

ps<-scale_reads(ps,100000)
ps <- ps %>%
  ps_mutate(
    Lib.conc_scaled = c(scale(Library.conc, center = TRUE, scale = TRUE)),
    TOC.scaled = c(scale(TOC, center =TRUE, scale = TRUE)),
    F = if_else(Preservation == "F", 1, 0),
    Z = if_else(Preservation == "Z", 1, 0),
    L = if_else(Preservation == "L", 1, 0)
  )
constr_ord <- ps %>%
  tax_filter(min_prevalence = 1/100, tax_level = "family") %>%
  tax_agg("family") %>%
  tax_transform("total") %>%
  ord_calc(method = "RDA", constraints = c("TOC.scaled","Lib.conc_scaled", "F", "L", "Z"))

ord_plot(constr_ord, plot_taxa = 1:8,tax_vec_length = 1.5,constraint_vec_length = 0.8,constraint_lab_style = list(alpha=0.2),tax_lab_style = list(alpha=0.2),size = 3, colour = "Preservation", shape = "Site") +
  scale_color_viridis(discrete=TRUE)

ggsave(filename = "constrained_family.png",plot = last_plot(),width = 20,height = 15,units = "cm")
