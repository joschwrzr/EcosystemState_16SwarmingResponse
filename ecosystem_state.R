# set up -----
# load libraries

library(phyloseq)
library(microViz)
library(patchwork)
library(ggplot2)
library(tidyverse)
library(viridis)
library(magrittr)

# load data 
setwd("data/")
ps <- readRDS("16S_phyloseq.RDS")



# PCA country and depth -----------------
country_col <- c("Sweden" = "#005293", "Denmark" = "#C8102E")

ps |>
  tax_transform(rank = "unique", trans = "clr") |>
  ord_calc("PCA") |>
  ord_plot(axes = c(1, 2), color = "Country", shape = "Depth", size = 2) +
  
  # Draw ellipses based on Depth with greyscale colors
  ggplot2::stat_ellipse(ggplot2::aes(group = interaction(Depth, Country), linetype = Depth), color = "grey50") +
  #ggplot2::stat_ellipse(ggplot2::aes(group = Depth, linetype = Depth), color = "grey50") +
  # Define specific colors for Country
  scale_colour_manual(values = country_col, name = "Origin") +
  scale_fill_manual(values = country_col, name = "Origin") +
  # Define specific shapes for Depth
  scale_shape_manual(values = c("0-5 cm" = 16, "5-10 cm" = 24, "10-20 cm" = 15), name = "Depth") +
  # Define linetype for Depth to differentiate ellipses
  scale_linetype_manual(values = 
                          c("0-5 cm" = "solid", 
                            "5-10 cm" = "dashed",
                            "10-20 cm" = "dotted"), name = "Depth") +
  theme_bw() +
  #ggtitle(label = "16S rDNA: Sample Origin and Depth Cluster")+
  labs(caption="") +
  
  # Add side density plots with Country colors
  ggside::geom_xsidedensity(aes(fill = Country), alpha = 0.5, show.legend = FALSE) +
  ggside::geom_ysidedensity(aes(fill = Depth), alpha = 0.5, show.legend = F) +
  ggside::theme_ggside_void()+
  # Adjust text size for axis titles, labels, and legend
  theme(
    text = element_text(size = 16), 
    axis.title = element_text(size = 18), 
    axis.text = element_text(size = 14), 
    legend.title = element_text(size = 16), 
    legend.text = element_text(size = 14)
  )


# PCA warming ---------------------
waring_col<-  c("ambient" = "#4575b4", "+3.0°C" = "#fdae61", "+6.0°C" = "#d73027" ) 

ps_D5 <- ps %>% subset_samples(Country == "Denmark") %>%
  subset_samples(Depth == "0-5 cm")
ps_D10 <- ps %>% subset_samples(Country == "Denmark") %>%
  subset_samples(Depth == "5-10 cm")
ps_D20 <- ps %>% subset_samples(Country == "Denmark") %>%
  subset_samples(Depth == "10-20 cm")

ps_S5 <- ps %>% subset_samples(Country == "Sweden") %>%
  subset_samples(Depth == "0-5 cm")
ps_S10 <- ps %>% subset_samples(Country == "Sweden") %>%
  subset_samples(Depth == "5-10 cm")
ps_S20 <- ps %>% subset_samples(Country == "Sweden") %>%
  subset_samples(Depth == "10-20 cm")



# Create a list of your phyloseq objects
ps_list <- list(ps_D5, ps_D10, ps_D20, 
                ps_S5, ps_S10, ps_S20)
# Assign names to the list elements
names(ps_list) <- c("D5", "D10", "D20", 
                    "S5", "S10", "S20")
#remove single ds
rm(ps_D5, ps_D10, ps_D20, 
   ps_S5, ps_S10, ps_S20)

clrPCA_list <- list() # for the clr matrixes
clrPCA_plot <- list() # for the plots

new_labels <- c(
  "D5"  = "Denmark 0-5cm",
  "D10" = "Denmark 5-10cm",
  "D20" = "Denmark 10-20cm",
  "S5"  = "Sweden 0-5cm",
  "S10" = "Sweden 5-10cm",
  "S20" = "Sweden 10-20cm"
)

for (name in names(ps_list)) {
  
  clrPCA_list[[name]] <- ps_list[[name]]%>%
    tax_fix(unknowns = c("Incertae Sedis", "Unknown Family")) |>  # clean unknowns
    tax_transform("clr", rank = "unique") |>  
    ord_calc(method = "PCA") 
  
  clrPCA_plot[[name]]  <- clrPCA_list[[name]]  |> 
    ord_plot(axes = c(1, 2), color = "Warming", shape = "Warming") +
    scale_color_manual(values = waring_col, name = "Warming") +
    scale_fill_manual(values = waring_col, name = "Warming") +
    theme_bw() +
    stat_chull(aes(colour = Warming)) +
    ggside::geom_xsidedensity(aes(fill = Warming), alpha = 0.8, show.legend = FALSE) +  
    ggside::geom_ysidedensity(aes(fill = Warming), alpha = 0.8, show.legend = FALSE) +
    ggside::theme_ggside_void() +
    ggtitle(label = new_labels[name]) +
    labs(caption="") 
}


clrPCA_plot |> 
  wrap_plots(nrow = 2, ncol = 3) +
  plot_layout(guides = 'collect') 


## PCA axes 3 and 4 ------------------
for (name in names(ps_list)) {
  
  clrPCA_list[[name]] <- ps_list[[name]]%>%
    tax_fix(unknowns = c("Incertae Sedis", "Unknown Family")) |>  # clean unknowns
    tax_transform("clr", rank = "unique") |>  
    ord_calc(method = "PCA") 
  
  clrPCA_plot[[name]]  <- clrPCA_list[[name]]  |> 
    ord_plot(axes = c(3,4), color = "Warming", shape = "Warming") +
    scale_color_manual(values = waring_col, name = "Warming") +
    scale_fill_manual(values = waring_col, name = "Warming") +
    theme_bw() +
    stat_chull(aes(colour = Warming)) +
    ggside::geom_xsidedensity(aes(fill = Warming), alpha = 0.8, show.legend = FALSE) +  
    ggside::geom_ysidedensity(aes(fill = Warming), alpha = 0.8, show.legend = FALSE) +
    ggside::theme_ggside_void() +
    ggtitle(label = new_labels[name]) +
    labs(caption="") 
}

clrPCA_plot |> 
  wrap_plots(nrow = 2, ncol = 3) +
  plot_layout(guides = 'collect') 


## PERMANOVA using Aitchison distance ----
Ai_dist <- ps %>%
  tax_transform("clr") %>%   
  dist_calc("euclidean")  # Aitchison distance (euclidean on CLR)

# Test influence of Country, Depth, Warming
Ai_dist %>%
  dist_permanova(
    seed = 1,
    variables = " Country + Depth + WarmingT", 
    n_perms = 9999)

# for subsets
dist_list <- lapply(ps_list, function(ps) {
  ps %>%
    tax_transform("clr") %>%
    dist_calc("euclidean")
})


permanova_list <- lapply(dist_list, function(x) {
  dist_permanova(
    x,
    seed = 1,
    variables = "WarmingT",
    n_perms = 9999
  )
})

permanova_list$D5
permanova_list$D10



# Differential Abundance ---------------------------
library(DESeq2)
# Apply phyloseq_to_deseq2 on each object in the list

dds_list <- lapply(ps_list, function(ps_obj) {
  phyloseq_to_deseq2(ps_obj, ~ WarmingT)
})
# Assign names to the list elements
names(dds_list) <-  c("D5", "D10", "D20", 
                      "S5", "S10", "S20")




# Apply DESeq() to each DESeqDataSet in the list
dds_list <- lapply(dds_list, function(dds_obj) {
  DESeq(dds_obj)
})

# Apply results() to each DESeqDataSet in the list
res_list <- lapply(dds_list, function(dds_obj) {
  results(dds_obj)
})

##  Vulcano Plot --------------
# Initialize an empty list to store the results
res_df_list <- list()

# Loop through the res_list
for (name in names(res_list)) {
  # Create a data frame
  res_df <- as.data.frame(res_list[[name]])
  res_df$ASV <- rownames(res_df)
  
  # Add significance label
  res_df$Significant <- ifelse(res_df$padj < 0.05, "padj < 0.05", "ns")
  
  
  # Save the data frame to the new list, using the same name
  res_df_list[[name]] <- res_df  # Or res_df_list[[name]] <- sign_ASVs if you filtered
  rm(res_df)
}


# Initialize list to store the volcano plots
Volcano_list <- list()
sign_list <- list()
# Loop through the res_df_list
for (name in names(res_df_list)) {
  
  # Filter for significant ASVs
  sign_list[[name]] <- res_df_list[[name]] %>%
    filter(padj < 0.05)
  # Create the volcano plot
  Volcano_list[[name]] <- res_df_list[[name]] %>%
    filter(!is.na(padj)) %>%
    ggplot(aes(x = log2FoldChange, y = -log10(padj), color = Significant)) +
    geom_point(alpha = 0.8) +
    scale_color_manual(values = c("padj < 0.05" = "red", "ns" = "grey")) +
    theme_minimal() +
    ggtitle(label = new_labels[name]) +
    labs(
      # title = paste("Volcano Plot -", name),
      x = "Log2 fold change",
      y = "-Log10 adjusted p-value"
    ) +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
    ggrepel::geom_text_repel(data = sign_list[[name]], 
                             aes(label = ASV), size = 3, max.overlaps = 20) +
    coord_cartesian(xlim = c(-2, 2), ylim = c(0, 4))
}



Volcano_list |> 
  wrap_plots(nrow = 2, ncol = 3) +
  
  plot_layout(guides = 'collect') &
  theme(legend.position = 'bottom')


## legend --------------------
# Ensure consistent legend levels
Volcano_list <- list()
sign_list <- list()

for (name in names(res_df_list)) {
  
  df <- res_df_list[[name]] %>%
    filter(!is.na(padj)) %>%
    mutate(Significant = ifelse(padj < 0.05, "padj < 0.05", "ns"),
           Significant = factor(Significant, levels = c("padj < 0.05", "ns")))  # <- Force both levels
  
  sign_list[[name]] <- df %>% filter(Significant == "padj < 0.05")
  
  Volcano_list[[name]] <- ggplot(df, aes(x = log2FoldChange, y = -log10(padj), color = Significant)) +
    geom_point(alpha = 0.8) +
    scale_color_manual(values = c("padj < 0.05" = "red", "ns" = "grey")) +
    theme_minimal() +
    ggtitle(label = new_labels[name]) +
    labs(
      x = "Log2 fold change",
      y = "p-value"
    ) +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
    ggrepel::geom_text_repel(data = sign_list[[name]], 
                             aes(label = ASV), size = 3, max.overlaps = 20) +
    coord_cartesian(xlim = c(-2, 2), ylim = c(0, 4))
  
  
  p <- Volcano_list[[name]]
  
  # Only first plot keeps the legend
  if (name != "D5") {
    p <- p + theme(legend.position = "none")
  }
  
  # Remove y-axis unless it's the first column
  if (name == "D10" |name ==  "D20"| name ==  "S10" |name ==  "S20") {
    p <- p + theme(axis.title.y = element_blank(),
                   axis.text.y = element_blank(),
                   axis.ticks.y = element_blank())
  }
  
  # Remove x-axis unless it's bottom row
  if (name == "D5" |name ==  "D10" | name== "D20") {
    p <- p + theme(axis.title.x = element_blank(),
                   axis.text.x = element_blank(),
                   axis.ticks.x = element_blank())
  }
  
  Volcano_list[[name]] <- p
  
  
}

# Final plot: share guides and put legend at the bottom
Volcano_list |> 
  wrap_plots(nrow = 2, ncol = 3, guides = "collect") 


## summary table -----------------------------------------
# Get taxonomy table
Comp_list <- list()
# Loop through the res_list
for (name in names(ps_list)) {
  tax_df <- as.data.frame(tax_table(ps_list[[name]])) %>%
    rownames_to_column(var = "ASV")
  
  Comp_list[[name]] <- sign_list[[name]] %>% 
    select(ASV, padj, log2FoldChange) %>%
    left_join(tax_df, by = "ASV") %>%
    mutate(label = paste(ASV, Genus, sep = " | "),
           Method = "DESeq2",
           Data = name)
  rm(tax_df)
}


Comp_df <- Comp_list %>%
  map(~ select(.x, -label)) %>%  # remove 'label' from each data frame
  bind_rows()
Comp_df <-Comp_df %>%
  mutate(Direction = ifelse(log2FoldChange > 0, "Up", "Down"))%>%
  arrange(Direction, Phylum, Class, Order, Family)

Comp_df


# tree ----------------------

library(Biostrings)
library(seqinr)

seqs <- read.table('table_sequences_asvNames.csv', sep = ";", stringsAsFactors = F, header = T)
DESeq_df<- Comp_df 

# Subset of interest
seqs <- subset(seqs, Name %in%  DESeq_df$ASV)

# als DNAstrings
library(msa)
seqs_string <- DNAStringSet(seqs[,1])
names(seqs_string) <- seqs$Name # keep ASV names!


# multiple sequence alignment
seq_msa <- msa(seqs_string)
# distance
library(seqinr)
d <- dist.alignment(msaConvert(seq_msa, type="seqinr::alignment"))
# tree
library(ape)
tree <- nj(d)

library(ggtree)
ggtree(tree, branch.length='none', layout='circular') + 
  geom_tiplab()

## root tree ---------
tree_rooted <- phytools::midpoint.root(tree)

## label ----
taxonomy_df<-Comp_df %>%
  mutate(Family = ifelse(Family == "Rhizobiales Incertae Sedis", "Rhizobiales", Family))
# Combine taxonomic ranks into a path string for tree conversion
taxonomy_df$pathString <- taxonomy_df %>%
  unite("path", Kingdom, Phylum, Class, Order, Family, Genus, ASV, sep = "/",
        na.rm = T, remove = FALSE) %>%
  pull(path)


ggtree(tree, branch.length="none") %<+% taxonomy_df + 
  geom_tippoint(aes(color = log2FoldChange), size = 3) +
  scale_color_gradient2(low = "blue", mid = "gray90", high = "red", midpoint = 0) +
  xlim(NA,15)+
  geom_tiplab(aes(label = Phylum), size = 3)


# Family
p <- ggtree(tree_rooted, 
            layout = "circular", 
            branch.length = "none") 
# Add  annotations
p %<+% taxonomy_df +
  geom_tippoint(aes(color = log2FoldChange), size = 3) +
  scale_color_gradient2(low = "blue", mid = "gray90", high = "red", midpoint = 0) +
  geom_tiplab(aes(label = Family), size = 3, hjust = -0.1) +
  ggtitle("midpoint rooted, phyla labelled, tipponits are Family")+
  
  # Add node numbers to the plot
  #geom_text2(aes(label = node), color = "black", size = 2, hjust = -0.3) + 
  
  # Cladelabels
  geom_cladelab(node = 80, label = "Planctomycetota", angle = 0, 
                fontsize = 3.5, offset = 12, hjust = 0.7,  vjust = -9) +
  geom_cladelab(node = 76, label = "Bacteroidota", angle = 0, 
                fontsize = 3.5, offset = 12, hjust = -0.1, vjust = 0)+
  geom_cladelab(node = 75, label = "Acidobacteriota", angle = 0, 
                fontsize = 3.5, offset = 12, hjust = 0.7, vjust = -1.2)+
  geom_cladelab(node = 29, label = "Acidobacteriota", angle = 0, 
                fontsize = 3.5, offset = 12, hjust = 1, vjust = 1.5)+
  geom_cladelab(node = 73, label = "Cyanobacteria", angle = 0, 
                fontsize = 3.5, offset = 12, hjust = 1.2, vjust = 1.5)+
  geom_cladelab(node = 74, label = "Actinobacteriota", angle = 0, 
                fontsize = 3.5, offset = 12, hjust = 0.8, vjust = 1.8)+
  #geom_cladelab(node = 47, label = "Actinobacteriota", angle = 0, 
  #             fontsize = 3.5, offset = 12, hjust = -0, vjust = 0.5)+
  geom_cladelab(node = 77, label = "Proteobacteria", angle = 0, 
                fontsize = 3.5, offset = 12, hjust = -0.1, vjust = 0) +
  # geom_cladelab(node = 11, label = "Proteobacteria", angle = 0, 
  #        fontsize = 3.5, offset = 12, hjust = -0.1, vjust = 0) +
  geom_cladelab(node = 113, label = "Firmicutes", angle = 0, 
                fontsize = 3.5, offset = 12, hjust = 0.5,  vjust = -1)+
  geom_strip(3, 10,   label = "Chloroflexi", angle = 0, 
             fontsize = 3.5, offset = 12, hjust = -0.2)






# Vegetation barplot -----------------------------------
plants <- readRDS("Vegetation_cover.RDS")

plants_long <- plants |>
  pivot_longer(cols = `bare soil`:Unspecified, 
               names_to = "Species", 
               values_to = "Cover")


plants_long$Species <- factor(plants_long$Species,
                              levels = c("bare soil",
                                         # tree
                                         "Alnus glutinosa","Betula pubescens" ,
                                         # graminos
                                         "Agrostis stolonifera", "Agrostis capilaris" , 
                                         "Carex sp.","Elymus repens" ,
                                         "Festuca rubra","Phleum pratense",
                                         "Phragmites australis",
                                         # herbs
                                         "Aster tripodium","Atriplex prostata" ,
                                         "Cirsium sp.","Erigeron acris" ,
                                         "Epilobium hirsutum","Epilobium sp.",
                                         "Eupatorium cannabinum","Galium palustre",
                                         "Lycopus europaeus" ,"Lythrum salicaria" ,
                                         "Matricaria maritima","Potentilla anserina",
                                         "Senecio sp." ,"Solidago sp.",
                                         "Sonchus sp.","Solanum dulcamare",
                                         "Taraxacum sp." ,"Viccia cracca", 
                                         "Unspecified"))


custom_colors <-  c("bare soil"= "#806E65",
                    # tree
                    "Alnus glutinosa" = "#E39D47",
                    "Betula pubescens"  ="#CDA70A",
                    # graminos
                    "Agrostis stolonifera"= "#82EBD6", 
                    "Agrostis capilaris" ="#00C89B",
                    "Carex sp." ="#00C685",
                    "Elymus repens" ="#00CD79",
                    "Festuca rubra" = "#00C57C",
                    "Phleum pratense" = "#74BE3F",
                    "Phragmites australis" = "#34A58D",
                    # herbs
                    "Aster tripodium"="#DA8DF7",
                    "Atriplex prostata" ="#BF98FF",
                    "Cirsium sp."="#9CA3FF", #kratzdistel 
                    "Erigeron acris" = "#E6D7FF",
                    "Epilobium hirsutum" ="#F582DE",
                    "Epilobium sp." ="#ED85EA",
                    "Eupatorium cannabinum" = "#DEA7C8", #wasserdost
                    "Galium palustre" ="#EDC6FC", # sumpflabkraut
                    "Lycopus europaeus" ="#FD80CA" , 
                    "Lythrum salicaria" = "#CE69E4", #blutweiderich
                    "Matricaria maritima" = "#EAADBE" , #strandkamille
                    "Potentilla anserina" = "#FFE7A7", #gänsefingerkraut
                    "Senecio sp." = "#B69BFF",
                    "Solidago sp."= "#FFE78D"  , #goldrute
                    "Sonchus sp." = "#EBE97C", #gänsedistel
                    "Solanum dulcamare" = "#5E23A4", #bittersüßer nachtschatten
                    "Taraxacum sp." = "#F0EC13",
                    "Vicia cracca" = "#D390FA" ,# Vogelwicke #D9678B
                    "Unspecified" = "#137F41")


plants_long |>
  group_by(`Warming treatment`, Species, Country) |>
  summarise(mean_cover = mean(Cover, na.rm = TRUE), .groups = "keep") |>
  #mutate(Country = factor(Country, levels = c("Sweden", "Denmark"), ordered = TRUE)) |>
  arrange(Country, `Warming treatment`, Species) |> 
  ggplot(aes(x = factor(`Warming treatment`), y = mean_cover, fill = Species)) +
  geom_bar(stat = "identity") +
  facet_grid(~ Country, scales = "free_x", space = "free") +
  labs(x = "Warming treatment", 
       y = "Cover (%)", 
       # title = "Vegetation Cover by Plot and Country", 
       fill = "Species") +
  theme_minimal() +
  scale_fill_manual(values = custom_colors)+ 
  
  theme(#axis.text.x = element_text(angle = 45, hjust = 1),
    #legend.position = "right",
    legend.text = element_text(size = 12), 
    legend.title = element_text(size = 14),
    legend.direction = "vertical" ,
    legend.key.size = unit(0.35, "cm"))+
  guides(color = guide_legend(ncol = 1), fill = guide_legend(ncol = 1), shape = guide_legend(ncol = 1))





# EEA ------

eea <- sample_data(ps)
eea <- data.frame(eea)
eea %<>% rownames_to_column("SampleID")

eea %<>% select(SampleID, Country, Warming, Depth, OM_perc,
                GLU = "OM_Activity_GLU",
                CEL = "OM_Activity_CEL",
                CHI = "OM_Activity_CHI",
                LEU = "OM_Activity_LEU",
                PHO = "OM_Activity_PHO"
                ) 


## Outlier Removal----
# Remove outliers higher than 10,000 for each enzyme --> set to NA
cat("GLU outliers:", sum(eea$GLU > 10000, na.rm = TRUE))
cat("CEL outliers:", sum(eea$CEL > 10000, na.rm = TRUE)) 
cat("CHI outliers:", sum(eea$CHI > 10000, na.rm = TRUE)) #2 outliers are removed here
cat("LEU outliers:", sum(eea$LEU > 10000, na.rm = TRUE))
cat("PHO outliers:", sum(eea$PHO > 10000, na.rm = TRUE))

eea_filtr <- eea %>% mutate(
  GLU = ifelse(GLU > 10000, NA, GLU),
  CEL = ifelse(CEL > 10000, NA, CEL),
  CHI = ifelse(CHI > 10000, NA, CHI),
  LEU = ifelse(LEU > 10000, NA, LEU),
  PHO = ifelse(PHO > 10000, NA, PHO)
)

## Data Transformation for Enzyme groups (Acquisition Types)----
# Create enzyme group data for acquisition types
eea_grouped <- eea_filtr %>%
  mutate(
    Carbon = GLU + CEL,    # Carbon acquisition enzymes
    Nitrogen = CHI + LEU,    # Nitrogen acquisition enzymes  
    Phosphorous = PHO           # Phosphorous acquisition enzyme
  ) %>%
  # Remove rows with NA values for PCA
  na.omit() %>%
  dplyr::select(SampleID, Country, Warming, Depth, OM_perc,
                Carbon, Nitrogen, Phosphorous) %>%
  as.data.frame()

## Boxplots aquisition Types---------------

# Plot enzyme groups (Carbon, Nitrogen, Phosphorous acquisition)
eea_grouped_aqtype <-
  eea_grouped %>% select(!OM_perc) %>%
  pivot_longer(cols = c(Carbon, 
                        Nitrogen, 
                        Phosphorous),
               names_to = "Enzyme_Group", values_to = "Activity") %>%
  mutate(Enzyme_Group = factor(Enzyme_Group, 
                               levels = c("Carbon", 
                                          "Nitrogen", 
                                          "Phosphorous"))
  ) 

plot_enzyme_groups <- eea_grouped_aqtype %>%
  ggplot(aes(x = Warming, y = Activity, fill = Warming)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA, width = 0.6, size = 0.3) +
  geom_jitter(width = 0.2, alpha = 0.6, size = 0.5) +
  # stat_summary(fun = mean, geom = "point", shape = 4, size = 2, color = "black", stroke = 1) +
  scale_fill_manual(values = waring_col,
                    name = "Warming") +
  labs(y = expression(paste("Activity [", 
                            "nmol ", g^{-1}, " OM ", h^{-1}, 
                            "]")
  ),
  x = "Temperature treatment") +
  theme_minimal() +
  facet_grid(Country ~ Enzyme_Group) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none") 


plot_enzyme_groups


## Principal component analysis -----


library(FactoMineR)
pca_eea <- eea_grouped %>% select(!OM_perc) %>% 
  PCA(scale.unit = TRUE,
      quali.sup = 1:4,    
      graph = F)

# PCA diagnostics
# Scree plot as diagnostic tool
library(factoextra)
fviz_eig(pca_eea, addlabels = TRUE, ylim = c(0, 70)) +
  ggtitle("Scree Plot - Variance Explained by Principal Components") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))



# Quality of representation plot
library(ggcorrplot)
quality_plot <- ggcorrplot(pca_eea$var$cos2, method = "circle") +
  scale_fill_gradient2(low = "white", high = "#8B4513",
                       limit = c(0, 1)) + 
  labs(fill = "Cos²") +
  ggtitle("Quality of Representation of Variables") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

print(quality_plot)

## PCA Biplots
pca_coords <- data.frame(
  PC1 = pca_eea$ind$coord[,1],
  PC2 = pca_eea$ind$coord[,2],
  Country = eea_grouped$Country,
  Warming = eea_grouped$Warming
)

# Extract variable coordinates for arrows
var_coords <- data.frame(
  PC1 = pca_eea$var$coord[,1] * 4,  
  PC2 = pca_eea$var$coord[,2] * 4,
  Variable = rownames(pca_eea$var$coord)
)

# Plot PCA
plot_biplot_main <-
  ggplot() +
  geom_point(data = pca_coords, 
             aes(x = PC1, y = PC2, color = Country, shape = Warming), 
             size = 3, alpha = 0.7) +
  stat_ellipse(data = pca_coords, 
               aes(x = PC1, y = PC2, color = Country, fill = Country), 
               alpha = 0.2, level = 0.95, geom = "polygon") +
  geom_segment(data = var_coords, 
               aes(x = 0, y = 0, xend = PC1, yend = PC2),
               arrow = arrow(length = unit(0.2, "cm")), 
               color = "black", linewidth  = 0.5) +
  geom_text(data = var_coords, 
            aes(x = PC1 * 1.1, y = PC2 * 1.1, label = Variable),
            color = "black", size = 4) +
  scale_color_manual(values = country_col,
                     name = "Country") +
  scale_fill_manual(values = country_col,
                    name = "Country") +
  scale_shape_manual(values = c("ambient" = 16, "+3.0°C" = 17, "+6.0°C" = 15),
                     name = "Treatment") +
  labs(x = paste0("PC1 (", round(pca_eea$eig[1,2], 1), "%)"),
       y = paste0("PC2 (", round(pca_eea$eig[2,2], 1), "%)")) +
  theme_minimal() +
  geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.5)

plot_biplot_main

library(ggpubr)
plot_panel_eea <- ggarrange(plot_enzyme_groups, plot_biplot_main, 
                            ncol = 2, labels = c("A", "B"), 
                            font.label = list(size = 14))
plot_panel_eea



#  Interpretation

cat("- Together, PC1 and PC2 explain", round(sum(pca_eea$eig[1:2,2]), 1), "% of the total variance\n")


## Correlation analysis -----

eea_grouped_corr <- eea_grouped %>%
  #filter(!is.na(OM) & !is.na(Carbon) & !is.na(Nitrogen) & !is.na(Phosphorous)) %>%
  mutate(Activity_Type = "OM-based")

# OM-based correlations by country
sweden_data_OM <- eea_grouped_corr[eea_grouped_corr$Country == "Sweden", ]
denmark_data_OM <- eea_grouped_corr[eea_grouped_corr$Country == "Denmark", ]

sweden_OM_C <- cor.test(sweden_data_OM$Carbon, sweden_data_OM$OM)
sweden_OM_N <- cor.test(sweden_data_OM$Nitrogen, sweden_data_OM$OM)
sweden_OM_P <- cor.test(sweden_data_OM$Phosphorous, sweden_data_OM$OM)

denmark_OM_C <- cor.test(denmark_data_OM$Carbon, denmark_data_OM$OM)
denmark_OM_N <- cor.test(denmark_data_OM$Nitrogen, denmark_data_OM$OM)
denmark_OM_P <- cor.test(denmark_data_OM$Phosphorous, denmark_data_OM$OM)



# OM-based activity plots
plot_C_OM_based <- ggplot(eea_grouped_corr, aes(x = OM_perc, y = Carbon)) +
  geom_point(aes(color = Country, shape = Warming), size = 2.5, alpha = 0.7) +
  geom_smooth(aes(color = Country), method = "lm", se = TRUE, linewidth = 1) +
  scale_color_manual(values = country_col) +
  scale_shape_manual(values = c("ambient" = 16, "+3.0°C" = 17, "+6.0°C" = 15)) +
  labs(x = "Organic Matter (%)", 
       y = expression(paste("Activity [", 
                            "nmol ", g^{-1}, " OM ", h^{-1}, 
                            "]"))) +
  theme_minimal() +
  theme(legend.position = "none") +
  annotate("text", x = Inf, y = Inf, hjust = 1.1, vjust = 2,
           label = paste0("Sweden: r = ", round(sweden_OM_C$estimate, 3), 
                          ", R² = ", round(sweden_OM_C$estimate^2, 3),
                          ", p = ", round(sweden_OM_C$p.value, 3)),
           color = "#2166ac", size = 3) +
  annotate("text", x = Inf, y = Inf, hjust = 1.1, vjust = 3.5,
           label = paste0("Denmark: r = ", round(denmark_OM_C$estimate, 3), 
                          ", R² = ", round(denmark_OM_C$estimate^2, 3),
                          ", p = ", round(denmark_OM_C$p.value, 3)),
           color = "#C8102E", size = 3)+
  coord_cartesian(ylim=c(0,5000))

plot_N_OM_based <- ggplot(eea_grouped_corr, aes(x = OM_perc, y = Nitrogen)) +
  geom_point(aes(color = Country, shape = Warming), size = 2.5, alpha = 0.7) +
  geom_smooth(aes(color = Country), method = "lm", se = TRUE, linewidth = 1) +
  scale_color_manual(values = country_col) +
  scale_shape_manual(values = c("ambient" = 16, "+3.0°C" = 17, "+6.0°C" = 15)) +
  labs(x = "Organic Matter (%)", 
       y = expression(paste("Activity [", 
                            "nmol ", g^{-1}, " OM ", h^{-1}, 
                            "]"))) +
  theme_minimal() +
  theme(legend.position = "none") +
  annotate("text", x = Inf, y = Inf, hjust = 1.1, vjust = 2,
           label = paste0("Sweden: r = ", round(sweden_OM_N$estimate, 3), 
                          ", R² = ", round(sweden_OM_N$estimate^2, 3),
                          ", p = ", round(sweden_OM_N$p.value, 3)),
           color = "#005293", size = 3) +
  annotate("text", x = Inf, y = Inf, hjust = 1.1, vjust = 3.5,
           label = paste0("Denmark: r = ", round(denmark_OM_N$estimate, 3), 
                          ", R² = ", round(denmark_OM_N$estimate^2, 3),
                          ", p = ", round(denmark_OM_N$p.value, 3)),
           color = "#C8102E", size = 3)+
  coord_cartesian(ylim=c(0,5000))

plot_P_OM_based <- ggplot(eea_grouped_corr, aes(x = OM_perc, y = Phosphorous)) +
  geom_point(aes(color = Country, shape = Warming), size = 2.5, alpha = 0.7) +
  geom_smooth(aes(color = Country), method = "lm", se = TRUE, linewidth = 1) +
  scale_color_manual(values = country_col) +
  scale_shape_manual(values = c("ambient" = 16, "+3.0°C" = 17, "+6.0°C" = 15)) +
  labs(x = "Organic Matter (%)", 
       y = expression(paste("Activity [", 
                            "nmol ", g^{-1}, " OM ", h^{-1}, 
                            "]"))) +
  theme_minimal() +
  theme(legend.position = "bottom") +
  guides(color = guide_legend(title = "Country"),
         shape = guide_legend(title = "Temperature")) +
  annotate("text", x = Inf, y = Inf, hjust = 1.1, vjust = 2,
           label = paste0("Sweden: r = ", round(sweden_OM_P$estimate, 3), 
                          ", R² = ", round(sweden_OM_P$estimate^2, 3),
                          ", p = ", round(sweden_OM_P$p.value, 3)),
           color = "#005293", size = 3) +
  annotate("text", x = Inf, y = Inf, hjust = 1.1, vjust = 3.5,
           label = paste0("Denmark: r = ", round(denmark_OM_P$estimate, 3), 
                          ", R² = ", round(denmark_OM_P$estimate^2, 3),
                          ", p = ", round(denmark_OM_P$p.value, 3)),
           color = "#C8102E", size = 3)+
  coord_cartesian(ylim=c(0,5000))


plot_C_OM_based <- plot_C_OM_based + labs(title = "Carbon acquisition")
plot_N_OM_based <- plot_N_OM_based + labs(title = "Nitrogen acquisition")
plot_P_OM_based <- plot_P_OM_based + labs(title = "Phosphorous acquisition")

correlation_panel_OM <- ggarrange(
  plot_C_OM_based, plot_N_OM_based, plot_P_OM_based,
  ncol = 3, common.legend = TRUE, legend = "bottom",
  labels = c("A", "B", "C"), font.label = list(size = 14)
)

correlation_panel_OM

## Statistical Analysis ----
library(lme4)
library(emmeans)
# Linear Mixed-Effects Models (main effects)
eea_grouped_aqtype %<>% mutate(
  Plot_ID = substr(SampleID, start = 1, stop = 2)
)

# Calculate statistical models lmer
# C-acquisition 
stat_model_C_acq <- lmer(Activity ~ Country * Warming * Depth + (1|Plot_ID),  
                         data = eea_grouped_aqtype %>% filter(Enzyme_Group == "Carbon"))

print(car::Anova(stat_model_C_acq, type = "II"))
# N-acquisition   
stat_model_N_acq <- lmer(Activity ~ Country * Warming * Depth + (1|Plot_ID), 
                         data = eea_grouped_aqtype %>% filter(Enzyme_Group == "Nitrogen"))

print(car::Anova(stat_model_N_acq, type = "II"))

# P-acquisition 
stat_model_P_acq <- lmer(Activity ~ Country * Warming * Depth + (1|Plot_ID), 
                         data = eea_grouped_aqtype %>% filter(Enzyme_Group == "Phosphorous"))

print(car::Anova(stat_model_P_acq, type = "II"))



## Post-hoc test (pairwise comparison) ---
# Carbon acquisition
print(pairs(emmeans(stat_model_C_acq, ~ Country)))
print(pairs(emmeans(stat_model_C_acq, ~ Warming))) # remember interaction effect on C-aq
pairs(emmeans(stat_model_C_acq, ~ Warming | Country)) # Sweden: +3.0°C significantly reduces carbon acquisition compared to ambient.
emmip(stat_model_C_acq, ~ Warming | Country)
print(pairs(emmeans(stat_model_C_acq, ~ Depth)))



# Nitrogen acquisition
print(pairs(emmeans(stat_model_N_acq, ~ Country)))
print(pairs(emmeans(stat_model_N_acq, ~ Warming)))
print(pairs(emmeans(stat_model_N_acq, ~ Depth)))

# Phosphorous acquisition
print(pairs(emmeans(stat_model_P_acq, ~ Country)))
print(pairs(emmeans(stat_model_P_acq, ~ Warming)))
print(pairs(emmeans(stat_model_P_acq, ~ Depth)))



