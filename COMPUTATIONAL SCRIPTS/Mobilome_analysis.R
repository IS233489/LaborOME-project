#Retrieve deduplicated AMR matrix (.CSV/.TXT) from AMR++ pipeline
library("installr")
install.packages("rlang")
library("rlang")
library("openxlsx")
library("phyloseq")
library("ggplot2")
library("dplyr")
library("ggthemes")
library("lme4")
library("lmerTest")
library("lsmeans")
library("ggrepel")
library("pbkrtest")
install.packages("ggpubr")
library("ggpubr")
install.packages("vioplot")
library(vioplot)
install.packages("patchwork")
library(patchwork)
library(metagenomeSeq)
library(scales)
library(vegan)
library(DESeq2)
library(pairwiseAdonis)
devtools::install_github("vmikk/metagMisc")
library(metagMisc)
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("ALDEx2")
library(ALDEx2)
install.packages("randomcoloR")
library(randomcoloR)
install.packages("zCompositions")
library(zCompositions)
pd = position_jitterdodge(dodge.width = 0.75, jitter.width = 0.3)


#setwd("")



########################################################
##                                                     #
###                                                    #
####Mobilome  analysis                                 #
###                                                    # 
##                                                     #
########################################################


## phyloseq object:
mob_mat<-read.csv(file.choose())      ## AMR_analytic_matrix_updatedID
mobgen_mat<- read.csv(file.choose()) ##WORKER_MOBILOME_ANNOTATION_FINAL
mobilomesamples_df <- read.csv(file.choose()) ##FINAL_SHOTGUN_MCOHS_WORKER_METADATA

##need to have row.names
row.names(mob_mat) <- mob_mat$Genes
mob_mat <- mob_mat %>% dplyr::select(-Genes) #remove the column Genes since it is now used as a row name

row.names(mobgen_mat) <- mobgen_mat$Genes
mobgen_mat <- mobgen_mat %>%dplyr::select(-Genes) 

row.names(mobilomesamples_df) <- mobilomesamples_df$NovSeq_SS_Lib_ID
mobilomesamples_df <- mobilomesamples_df %>%dplyr::select (-NovSeq_SS_Lib_ID)

#Transform into matrixes otu and tax tables (sample table can be left as data frame)

mob_mat <- as.matrix(mob_mat)
mobgen_mat <- as.matrix(mobgen_mat)

gplots::venn(list(mobilomesamples_df=rownames(mobilomesamples_df), featuretable=colnames(mob_mat))) #checking to make sure that metadata (mobilomesamples_df) sample names match mob_mat counts matrix sample names

##phyloseq:
MOBILOME = otu_table(mob_mat, taxa_are_rows = T)
MOBGENES = tax_table(mobgen_mat)
MOB_SAMPLES = sample_data(mobilomesamples_df)
MOBILOME_physeq <- phyloseq(MOBILOME, MOBGENES, MOB_SAMPLES)
saveRDS(MOBILOME_physeq, "initial.mobilome.phyloseq")
MOBILOME_physeq <- readRDS("initial.mobilome.phyloseq")

sample_data(MOBILOME_physeq)
##No environment option
noenvMOBILOME_physeq<-subset_samples(MOBILOME_physeq, CollectionPhase== "POST_SHOWER" | CollectionPhase== "SWINE" | CollectionPhase== "WORKDAY_END" | CollectionPhase== "WORKDAY_START")


##Without environment option metadata
noenvMOBILOME_physeq
NOENV_META<-subset_samples(MOB_SAMPLES, CollectionPhase)


mobtype.ps <-  tax_glom(noenvMOBILOME_physeq, "Type")
mobtype_VIRUS.ps <- subset_taxa(noenvMOBILOME_physeq, Type=="Virus")
mobtype_PLASMID.ps <- subset_taxa(noenvMOBILOME_physeq, Type=="Plasmid" | Type=="Plasmid replicon" | Type=="UNCLASSIFIED" | Type=="Plasmid ARG")
mobtype_PROPHAGE.ps <- subset_taxa(noenvMOBILOME_physeq, Type=="Prophage")
mobtype_TE.ps <- subset_taxa(noenvMOBILOME_physeq, Type=="TE")
mobtype_IS.ps <- subset_taxa(noenvMOBILOME_physeq, Type=="IS")
mobtype_ICE.ps <- subset_taxa(noenvMOBILOME_physeq, Type=="ICE")
mobtype_PLASMIDREP.ps <- subset_taxa(noenvMOBILOME_physeq, Type=="Plasmid replicon")
mobtype_PLASMIDARG.ps <- subset_taxa(noenvMOBILOME_physeq, Type=="Plasmid ARG")
mobtype_UNCLASSIFIED.ps <- subset_taxa(noenvMOBILOME_physeq, Type=="UNCLASSIFIED")



##With environment option

mobtype.ps <-  tax_glom(MOBILOME_physeq, "Type")
mobtype_VIRUS.ps <- subset_taxa(MOBILOME_physeq, Type=="Virus")
mobtype_PLASMID.ps <- subset_taxa(MOBILOME_physeq, Type=="Plasmid" | Type=="Plasmid replicon" | Type=="UNCLASSIFIED" | Type=="Plasmid ARG")
mobtype_PLASMIDonly.ps <- subset_taxa(MOBILOME_physeq, Type=="Plasmid")
mobtype_PROPHAGE.ps <- subset_taxa(MOBILOME_physeq, Type=="Prophage")
mobtype_TE.ps <- subset_taxa(MOBILOME_physeq, Type=="TE")
mobtype_IS.ps <- subset_taxa(MOBILOME_physeq, Type=="IS")
mobtype_ICE.ps <- subset_taxa(MOBILOME_physeq, Type=="ICE")
mobtype_PLASMIDREP.ps <- subset_taxa(MOBILOME_physeq, Type=="Plasmid replicon")
mobtype_PLASMIDARG.ps <- subset_taxa(MOBILOME_physeq, Type=="Plasmid ARG")
mobtype_UNCLASSIFIED.ps <- subset_taxa(MOBILOME_physeq, Type=="UNCLASSIFIED")

################################################################################
#                           ALPHA DIVERSITY                                     #      
################################################################################

### Relative abudance-- Mobilome type: Summarized across all samples per collectionphase
mobtype_relative.css<- phyloseq_transform_css(MOBILOME_physeq)
mobtype_relative.css
mobtype_relative <- transform_sample_counts(MOBILOME_physeq, function(x) {x/sum(x)}*100)
mobtype_relative_aglom <- tax_glom(mobtype_relative, taxrank = "Type") # 9 MOB types

mobtype_relative_long <- psmelt(mobtype_relative_aglom)
mobtype_palette <- distinctColorPalette(9)

ggplot(mobtype_relative_long, aes(x= CollectionPhase, y= Abundance, fill= factor(Type, levels = c("Virus", "Prophage","TE", "Plasmid replicon", "Plasmid ARG", "ICE", "IS", "Plasmid", "UNCLASSIFIED")))) +
  theme_bw() + 
  labs(y= "Relative abundance (%)") +
  geom_bar(stat = "summary", colour = "black", position = "fill") +
  scale_fill_manual(values = mobtype_palette) +
  scale_y_continuous(expand = c(0.0005,0,0.0005,0)) +
  scale_x_discrete(limits = c("WORKDAY_START", "WORKDAY_END","POST_SHOWER","SWINE", "ENVIRONMENT"), labels = c("Workday start","Workday end","Post-shower","Swine", "Environment")) +
  theme(#legend.position = "none",
    panel.border = element_rect(colour = "black", size = 1),
    plot.margin = unit(c(0.1,0.5,0.5,0.5), "cm"),
    strip.background = element_blank(),
    strip.text = element_text(size =24, colour = "black"),
    axis.text.y = element_text(size = 20, colour = "black"),
    axis.text.x = element_text(size =32, colour = "black",angle = 45, hjust = 0.95, vjust = 0.95),
    axis.title.y = element_text(size = 32, vjust = 1.75),
    axis.title.x = element_blank(),
    plot.title = element_text(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.y = element_blank())

### Relative abudance-- Mobilome type: Displaying all samples per collectionphase
MOBILOME_physeq <- readRDS("initial.mobilome.phyloseq")
mobtype_palette <- distinctColorPalette(9)
mob_relative <- transform_sample_counts(MOBILOME_physeq, function(x) x / sum(x) )
mob_relative_long <- psmelt(mob_relative)
mobtype_relative_long <- mob_relative_long %>%
  group_by(Type) %>%
  mutate(mean_relative_abund = mean(Abundance))

mobtype_relative_long$Type <- as.character(mobtype_relative_long$Type)
mobtype_relative_long$mean_relative_abund <- as.numeric(mobtype_relative_long$mean_relative_abund)

mobtype_relative_long$CollectionPhase <- factor(mobtype_relative_long$CollectionPhase,
                                                levels = c("WORKDAY_START", "WORKDAY_END", "POST_SHOWER", 'SWINE', "ENVIRONMENT"))

mobtype_relative_long %>%
  ggplot(aes(x = Sample, y = Abundance * 100, fill = factor(Type, levels = c("Virus", "Prophage", "TE", "Plasmid replicon", "Plasmid ARG", "ICE", "IS", "Plasmid", "UNCLASSIFIED")))) +
  geom_bar(stat = "identity") +
  facet_wrap(~CollectionPhase, scales = "free_x", nrow = 1, labeller = labeller(CollectionPhase = c("WORKDAY_START" = "Workday start", "WORKDAY_END" = "Workday end", "POST_SHOWER" = "Post-shower", "SWINE" = "Swine", "ENVIRONMENT" = "Environment"))) +
  geom_bar(stat = "identity", alpha = 0.9, position = "fill") +
  scale_fill_manual(values = mobtype_palette) +
  theme(plot.margin = unit(c(0.1, 0.5, 0.5, 0.5), "cm"),
        axis.title.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(size = 14),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 14),
        plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
        panel.background = element_rect(fill = "white"),
        strip.background = element_rect(fill = "white"),
        strip.text = element_text(face = "bold"),
        text = element_text(family = "Times New Roman", size = 14)) +
  ylab("Relative abundance (%)") +
  xlab("Sample")

# Define a custom color palette
type_palette <- c(
  "Virus" = "#FF7F00",
  "Prophage" = "#984EA3",
  "TE" = "#4DAF4A",
  "Plasmid replicon" = "#377EB8",
  "Plasmid ARG" = "#E41A1C",
  "ICE" = "#00A08A",
  "IS" = "#1B9E77",
  "Plasmid" = "#D95F02",
  "UNCLASSIFIED" = "#7570B3"
)

# Customize the theme
visual_theme <- theme_minimal(base_size = 12) +
  theme(
    axis.title = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 12),
    axis.text.x = element_blank(),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 12),
    plot.title = element_text(size = 16, face = "bold"),
    plot.margin = margin(1, 1, 1, 1, "cm"),
    panel.background = element_rect(fill = "white"),
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "white"),
    strip.text = element_text(size = 14, face = "bold"),
    text = element_text(family = "serif")
  )+theme(legend.position = "none")

# Create the plot
ggplot(mobtype_relative_long, aes(x = Sample, y = Abundance * 100, fill = Type)) +
  geom_bar(stat = "identity", alpha = 0.9) +
  facet_wrap(~ CollectionPhase, scales = "free_x", nrow = 1, labeller = labeller(CollectionPhase = c("WORKDAY_START" = "Workday start", "WORKDAY_END" = "Workday end", "POST_SHOWER" = "Post-shower", "SWINE" = "Swine", "ENVIRONMENT" = "Environment"))) +
  scale_fill_manual(values = type_palette) +
  labs(x = "Sample", y = "Relative abundance (%)") +
  visual_theme


### Relative abudance-- Plasmid type
mobtype_relative.css<- phyloseq_transform_css(mobtype_PLASMIDonly.ps)
mobtype_relative.css
mobtype_relative <- transform_sample_counts(mobtype_PLASMIDonly.ps, function(x) {x/sum(x)}*100)
mobtype_relative
mobtype_relative_aglom <- tax_glom(mobtype_relative, taxrank = "Plasmid_structural_accessory_genes") # 12 Plasmid types
otu_table(mobtype_relative)
mobtype_relative_long <- psmelt(mobtype_relative_aglom)
plasmid_palette <- distinctColorPalette(k=12)

ggplot(mobtype_relative_long, aes(x= CollectionPhase, y= Abundance, fill= factor(Plasmid_structural_accessory_genes, levels = c("Plasmid virulence and pathogenicity domain","Plasmid secretion or conjugation machinery","Plasmid non_conjugative efflux_transport","Plasmid stress_SOS_tox_antitox_partitioning","Plasmid relaxase and mobilization machinery","Plasmid transposition_recombination_integration","Plasmid biosynthesis_regulatory_enzymes","PTTREG","PREPIM","Hypothetical plasmid protein","Unclassified plasmid gene","Plasmid ARG")))) +
  theme_bw() + 
  labs(y= "Relative abundance (%)") +
  geom_bar(stat = "summary", colour = "black") +
  scale_fill_manual(values = plasmid_palette) +
  scale_y_continuous(expand = c(0.0005,0,0.0005,0)) +
  scale_x_discrete(limits = c("WORKDAY_START", "WORKDAY_END","POST_SHOWER","SWINE", "ENVIRONMENT"), labels = c("Workday start","Workday end","Post-shower","Swine", "Environment")) +
  theme(#legend.position = "none",
    panel.border = element_rect(colour = "black", size = 1),
    plot.margin = unit(c(0.1,0.5,0.5,0.5), "cm"),
    strip.background = element_blank(),
    strip.text = element_text(size =24, colour = "black"),
    axis.text.y = element_text(size = 20, colour = "black"),
    axis.text.x = element_text(size =32, colour = "black",angle = 45, hjust = 0.95, vjust = 0.95),
    axis.title.y = element_text(size = 32, vjust = 1.75),
    axis.title.x = element_blank(),
    plot.title = element_text(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.y = element_blank())

### Relative abudance-- Plasmidome genes: Displaying all samples per collectionphase
mobtype_PLASMIDonly.ps

PLASMIDonly_relative <- transform_sample_counts(mobtype_PLASMIDonly.ps, function(x) x / sum(x) )
PLASMIDonly_relative_long <- psmelt(PLASMIDonly_relative)
PLASMIDonly_relative_long <- PLASMIDonly_relative_long %>%
  group_by("Plasmid_structural_accessory_genes")  %>%  #12 Plasmid types)
  mutate(mean_relative_abund = mean(Abundance))

PLASMIDonly_relative_long$Plasmid_structural_accessory_genes <- as.character(PLASMIDonly_relative_long$Plasmid_structural_accessory_genes)
PLASMIDonly_relative_long$mean_relative_abund <- as.numeric(PLASMIDonly_relative_long$mean_relative_abund)

PLASMIDonly_relative_long$CollectionPhase <- factor(PLASMIDonly_relative_long$CollectionPhase, levels = c("WORKDAY_START", "WORKDAY_END", "POST_SHOWER", 'SWINE', "ENVIRONMENT"))

plasmid_palette <-c( "#1F78B4",
                     "#33A02C",
                     "#B15928",
                     "#E31A1C",
                     "#A6CEE3",
                     "#B2DF8A",
                     "#FDBF6F",
                     "#FF7F00",
                     "#CAB2D6",
                     "#6A3D9A",
                     "#FFFF99",
                     "#FB9A99"
)

# Customize the theme
visual_theme <- theme_minimal(base_size = 12) +
  theme(
    axis.title = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 12),
    axis.text.x = element_blank(),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 12),
    plot.title = element_text(size = 16, face = "bold"),
    plot.margin = margin(1, 1, 1, 1, "cm"),
    panel.background = element_rect(fill = "white"),
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "white"),
    strip.text = element_text(size = 14, face = "bold"),
    text = element_text(family = "serif")
  )+theme(legend.position = "none")

# Create the plot
ggplot(PLASMIDonly_relative_long, aes(x = Sample, y = Abundance * 100, fill = Plasmid_structural_accessory_genes)) +
  geom_bar(stat = "identity", alpha = 0.9) +
  facet_wrap(~ CollectionPhase, scales = "free_x", nrow = 1, labeller = labeller(CollectionPhase = c("WORKDAY_START" = "Workday start", "WORKDAY_END" = "Workday end", "POST_SHOWER" = "Post-shower", "SWINE" = "Swine", "ENVIRONMENT" = "Environment"))) +
  scale_fill_manual(values = plasmid_palette) +
  labs(x = "Sample", y = "Relative abundance (%)") +
  visual_theme


## Relative abudance--IS type
mobtype_IS.ps.pruned<-subset_samples(mobtype_IS.ps, sample_names(mobtype_IS.ps) != "R1_19")
mobtype_relative.css<- phyloseq_transform_css(mobtype_IS.ps.pruned)
mobtype_relative.css
mobtype_relative <- transform_sample_counts(mobtype_IS.ps.pruned, function(x){x/sum(x)}*100)
mobtype_relative
mobtype_relative_aglom <- tax_glom(mobtype_relative, taxrank = "IS_Family") # 24 IS families
mobtype_relative_aglom
mobtype_relative_MELT <- psmelt(mobtype_relative_aglom)
mobtype_relative_MELT
is_palette <- distinctColorPalette(k=24)

ggplot(mobtype_relative_MELT, aes(x= CollectionPhase, y= Abundance, fill= factor(IS_Family, levels = c("None","ISAs1","IS21","IS3", "IS110","IS66", "IS1595",  "IS1182", "IS30", "IS630", "IS982", "ISNCY","IS5","IS200/IS605","Tn3","IS1","IS1380","ISL3","IS256","IS91","IS481","IS4", "IS6","ISLre2")))) +
  theme_bw() + 
  geom_bar(stat = "summary", colour = "black")+
  labs(y= "Normalized abundance (%)") +
  geom_bar(stat = "summary", colour = "black") +
 # geom_errorbar()+
  scale_fill_manual(values = is_palette) +
  scale_y_continuous(expand = c(0.0005,0,0.0005,0)) +
  scale_x_discrete(limits = c("WORKDAY_START", "WORKDAY_END","POST_SHOWER","SWINE", "ENVIRONMENT"), labels = c("Workday start","Workday end","Post-shower","Swine", "Environment")) +
  theme(#legend.position = "none",
    panel.border = element_rect(colour = "black", size = 1),
    plot.margin = unit(c(0.1,0.5,0.5,0.5), "cm"),
    strip.background = element_blank(),
    strip.text = element_text(size =24, colour = "black"),
    axis.text.y = element_text(size = 20, colour = "black"),
    axis.text.x = element_text(size =32, colour = "black",angle = 45, hjust = 0.95, vjust = 0.95),
    axis.title.y = element_text(size = 30, vjust = 1.75),
    axis.title.x = element_blank(),
    plot.title = element_text(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.y = element_blank())




### Relative abudance-- IS genes: Displaying all samples per collectionphase
mobtype_IS.ps

IS_relative <- transform_sample_counts(mobtype_IS.ps, function(x) x / sum(x) )
IS_relative_long <- psmelt(IS_relative)
IS_relative_long <- IS_relative_long %>%
  group_by("IS_Family")  %>%  #4 IS families)
  mutate(mean_relative_abund = mean(Abundance))

IS_relative_long$IS_Family <- as.character(IS_relative_long$IS_Family)
IS_relative_long$mean_relative_abund <- as.numeric(IS_relative_long$mean_relative_abund)

IS_relative_long$CollectionPhase <- factor(IS_relative_long$CollectionPhase, levels = c("WORKDAY_START", "WORKDAY_END", "POST_SHOWER", 'SWINE', "ENVIRONMENT"))


IS_palette <-      c(
  "#B15928", "#FF7F0E", "#2CA02C", "#D62728", "#9467BD", "#8C564B",
  "#E377C2", "#7F7F7F", "#BCBD22", "#17BECF", "#525252", "#6A3D9A",
  "#33A02C", "#1F77B4", "#A6CEE3", "#B2DF8A", "#FB9A99", "#FDBF6F",
  "#CAB2D6", "#FEE08B", "#E31A1C", "#FDBF6F", "#1B9E77", "#D95F02"
)

# Customize the theme
visual_theme <- theme_minimal(base_size = 12) +
  theme(
    axis.title = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 12),
    axis.text.x = element_blank(),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 12),
    plot.title = element_text(size = 16, face = "bold"),
    plot.margin = margin(1, 1, 1, 1, "cm"),
    panel.background = element_rect(fill = "white"),
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "white"),
    strip.text = element_text(size = 14, face = "bold"),
    text = element_text(family = "serif")
  )+theme(legend.position = "none")

# Create the plot
ggplot(IS_relative_long, aes(x = Sample, y = Abundance * 100, fill = IS_Family)) +
  geom_bar(stat = "identity", alpha = 0.9) +
  facet_wrap(~ CollectionPhase, scales = "free_x", nrow = 1, labeller = labeller(CollectionPhase = c("WORKDAY_START" = "Workday start", "WORKDAY_END" = "Workday end", "POST_SHOWER" = "Post-shower", "SWINE" = "Swine", "ENVIRONMENT" = "Environment"))) +
  scale_fill_manual(values = IS_palette) +
  labs(x = "Sample", y = "Relative abundance (%)") +
  visual_theme

###############################################################################
###ALPHA DIVERSITY#############################################################
################################################################################

alpha_div_mob <- estimate_richness(MOBILOME_physeq, measures = c("Observed","Shannon","Simpson","InvSimpson"))

alpha_div_mob
alpha_div.mobilome.df <- as(sample_data(MOBILOME_physeq), "data.frame")
alpha_div_mobilome.meta<- cbind(alpha_div_mob, alpha_div.mobilome.df)
alpha_div_mobilome.meta

mobrichplot<-
  alpha_div_mobilome.meta %>%
mutate(CollectionPhase= fct_relevel(CollectionPhase,"WORKDAY_START", "WORKDAY_END", "POST_SHOWER", "SWINE", "ENVIRONMENT"))%>%
ggplot(aes(x=CollectionPhase, y = Observed, fill = CollectionPhase)) + theme_bw() + labs(title= "", y= "Mobilome richness", x= "Collection phase") + geom_jitter(width = 0.2, shape=21, color= "black", size= 3) +
  geom_boxplot(alpha=0.3, shape=31, color= "black", outlier.colour = 'white', width=0.6, coef=5000, fatten=1, lwd=1) + ylim(8,NA) + 
  theme(plot.margin = unit(c(0.3,0.06,0.3,0.06), "cm"),
        strip.background = element_blank(),
        strip.text = element_text(size =24, colour = "black"),
        axis.text = element_text(size = 14, colour = "black"),
        axis.text.x = element_blank(),
        axis.title.y = element_text(size = 20, vjust = 1.75),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.title = element_text(),
        panel.border = element_rect(colour = "black", size = 1.0),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank()) +  scale_fill_manual(values=c("#789ED0", "#F79D40", "#72BF56", "#E4655A", "#AB8BCB"))+ scale_color_manual(values=c("#789ED0", "#F79D40", "#72BF56", "#E4655A", "#AB8BCB"))+ scale_fill_discrete("Collection phase", labels= c("WORKDAY_START"= "Workday start", "WORKDAY_END"="Workday end", "POST_SHOWER"="Post-shower", "SWINE"="Swine", "ENVIRONMENT"="Environment"))

mobrichplot

#richness_stats <- pairwise.wilcox.test(alpha_div_meta$Observed, alpha_div_meta$CollectionPhase, p.adjust.method = "BH")
#richness_stats$p.value
#write.csv(richness_stats[["p.value"]],"richness_stats.csv")

mob_rich_model <-subset(alpha_div.mobilome.df, !CollectionPhase %in% c("MockComm", "NegCtrl"))
mobmodel_richness<-
  lmer(alpha_div_mobilome.meta$Observed ~ CollectionPhase + Shotgun_Post_Host_Rem + Post_captureSampleConc_ng_ul + (1|Worker), data = mob_rich_model)
summary(mobmodel_richness)
anova(mobmodel_richness, type='III',test="F")
lsmeans(mobmodel_richness, pairwise~CollectionPhase, adjust="Tukey")
VarCorr(mob_rich_model)
#$contrasts
#contrast                    estimate  SE   df t.ratio p.value
#ENVIRONMENT - POST_SHOWER     -652.1 335 27.1  -1.944  0.3193
#ENVIRONMENT - SWINE           -593.7 332 26.8  -1.791  0.3996
#ENVIRONMENT - WORKDAY_END     -630.1 384 32.8  -1.642  0.4823
#ENVIRONMENT - WORKDAY_START   -985.9 354 29.1  -2.785  0.0653
#POST_SHOWER - SWINE             58.4 201 28.4   0.290  0.9984
#POST_SHOWER - WORKDAY_END       22.0 171 19.4   0.129  0.9999
#POST_SHOWER - WORKDAY_START   -333.8 142 19.5  -2.356  0.1699
#SWINE - WORKDAY_END            -36.3 246 34.4  -0.148  0.9999
#SWINE - WORKDAY_START         -392.2 239 33.4  -1.640  0.4834
#WORKDAY_END - WORKDAY_START   -355.8 175 18.9  -2.037  0.2871


mobshannonplot<-
  alpha_div_mobilome.meta %>%
  mutate(CollectionPhase= fct_relevel(CollectionPhase,"WORKDAY_START", "WORKDAY_END", "POST_SHOWER", "SWINE", "ENVIRONMENT"))%>%

ggplot(aes(x=CollectionPhase, y = Shannon, fill = CollectionPhase)) + theme_bw() + labs(title= "", y= "Shannon's diversity of MGEs", x= "Collection phase") + geom_jitter(width = 0.2, shape=21, color= "black", size= 3) +
  geom_boxplot(alpha=0.3, shape=31, color= "black", outlier.colour = 'white', width=0.6, coef=5000, fatten=1, lwd=1) + ylim(0,NA) + scale_fill_tableau(palette = "Classic 10 Medium", type = "regular",direction = 1) + scale_color_tableau(palette = "Classic 10 Medium", type = "regular",direction = 1) + 
  scale_fill_discrete(name= "CollectionPhase", labels= c("WORKDAY_START"= "Workday start", "WORKDAY_END"="Workday end", "POST_SHOWER"="Post-shower", "SWINE"="Swine", "ENVIRONMENT" = "Environment")) + 
  theme(plot.margin = unit(c(0.1,0.5,0.5,0.5), "cm"),
        strip.background = element_blank(),
        strip.text = element_text(size =24, colour = "black"),
        axis.text = element_text(size = 14, colour = "black"),
        axis.text.x = element_blank(),
        axis.title.y = element_text(size = 20, vjust = 1.75),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.title = element_text(),
        panel.border = element_rect(colour = "black", size = 1.0),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank())+ scale_fill_manual(values=c("#789ED0", "#F79D40", "#72BF56", "#E4655A", "#AB8BCB"))+ scale_color_manual(values=c("#789ED0", "#F79D40", "#72BF56", "#E4655A", "#AB8BCB"))+ scale_fill_discrete("Collection phase", labels= c("WORKDAY_START"= "Workday start", "WORKDAY_END"="Workday end", "POST_SHOWER"="Post-shower", "SWINE"="Swine", "ENVIRONMENT"="Environment"))
mobshannonplot


mob_shanon_model <-subset(alpha_div.mobilome.df, !CollectionPhase %in% c("MockComm", "NegCtrl"))
mobmodel_shannon<-
  lmer(alpha_div_mobilome.meta$Shannon ~ CollectionPhase + Shotgun_Post_Host_Rem + Post_captureSampleConc_ng_ul + (1|Worker), data = mob_shanon_model)
summary(mobmodel_shannon)
anova(mobmodel_shannon, type='III',test="F")
lsmeans(mobmodel_shannon, pairwise~CollectionPhase, adjust="Tukey")
VarCorr(mobmodel_shannon)

#$contrasts
#contrast                    estimate    SE   df t.ratio p.value
#ENVIRONMENT - POST_SHOWER    -0.8379 0.432 32.9  -1.941  0.3169
#ENVIRONMENT - SWINE          -0.5244 0.426 32.8  -1.230  0.7341
#ENVIRONMENT - WORKDAY_END    -0.7607 0.519 35.0  -1.466  0.5909
#ENVIRONMENT - WORKDAY_START  -1.3110 0.461 33.7  -2.841  0.0548
#POST_SHOWER - SWINE           0.3135 0.261 33.4   1.201  0.7506
#POST_SHOWER - WORKDAY_END     0.0771 0.274 22.1   0.282  0.9985
#POST_SHOWER - WORKDAY_START  -0.4731 0.227 20.8  -2.086  0.2630
#SWINE - WORKDAY_END          -0.2364 0.337 35.0  -0.701  0.9549
#SWINE - WORKDAY_START        -0.7866 0.322 34.9  -2.442  0.1281
#WORKDAY_END - WORKDAY_START  -0.5502 0.282 21.8  -1.952  0.3214

############################RAREFACTION ANALYSIS

subsetted_biome_asv.ps<- subset_samples(biome_asv.ps, CollectionPhase=="WORKDAY_START" | CollectionPhase=="WORKDAY_END" | CollectionPhase=="POST_SHOWER" |CollectionPhase=="SWINE" | CollectionPhase=="ENVIRONMENT")

################################################################################
#                           BETA DIVERSITY                                     #      
################################################################################
#Virus/Virulence factors--> VIRUS ordination, even with imputation was too sparse to yield a convergent NMDS solution, and when one was produced, none of the collection phase groups were statistically significantly different

#sum(taxa_sums(mobtype_VIRUS.ps)==0) # 0 taxa not present across all samples
#VIRUS_physeq.ps.pruned <- prune_taxa(taxa_sums(mobtype_VIRUS.ps) > 30, mobtype_VIRUS.ps)#prunning low-abundance taxa #of gene groups

#VIRUS_physeq.ps.pruned.css <- phyloseq_transform_css(VIRUS_physeq.ps.pruned, log = F)
#IRUS_physeq

#VIRUS_physeq.ps.pruned.dist <- vegdist(decostand(t(otu_table(VIRUS_physeq.ps.pruned)),"rclr"), method = "euclidean")
#set.seed(1999)
#VIRUS_physeq.ps.pruned.ord <- vegan::metaMDS(comm = dist((VIRUS_physeq.ps.pruned.dist)), distance = "euclidean", k=5, try = 50, trymax = 1000, autotransform = F)

#VIRUS_plot<-plot_ordination(VIRUS_physeq.ps.pruned.css,VIRUS_physeq.ps.pruned.ord, type = "samples", color="CollectionPhase")
#VIRUS_plot

#VIRUScentroid<- VIRUS_plot$data %>% group_by(CollectionPhase) %>% summarise(NMDS1=mean(NMDS1), NMDS2=mean(NMDS2), .groups = "drop")

#stressplot(plot_ord_all)
#str(plot_ord_all)

#centroid
#VIRUS_plot$mapping
#VIRUS_ordplot <- ggplot(VIRUS_plot$data,VIRUS_plot$mapping) + stat_ellipse(geom = "polygon", type= "t", level = 0.85, lty =1, lwd = 1.3, alpha= 0.07, show.legend =FALSE) + theme_bw() +geom_point(aes(fill=CollectionPhase), size=6, alpha=0.3, stroke=1.0,) + geom_point(data=VIRUScentroid, size=13, shape=16, fill="black",mapping= aes(colour=CollectionPhase), stroke=4, alpha= 0.8, show.legend = FALSE) + scale_color_manual(name="Collection phase", breaks = c("WORKDAY_START", "WORKDAY_END", "POST_SHOWER", "SWINE", "ENVIRONMENT"), values = c("#789ED0", "#F79D40", "#72BF56", "#E4655A", "#AB8BCB"), labels= c("WORKDAY_START"= "Workday start", "WORKDAY_END"="Workday end", "POST_SHOWER"="Post-shower", "SWINE"="Swine", "ENVIRONMENT"="Environment")) +
 # theme(axis.line = element_line(), 
  #      text=element_text(family="Times",  size=10),
  #      panel.background = element_blank(),
  #      plot.margin = unit(c(0.3,0.06,0.3,0.06), "cm"),
  #      strip.background = element_blank(),
  #      strip.text = element_blank(),
  #      legend.key.size = unit(0.02, 'cm'),
  #      legend.position = c(0.8,0.8),
  #      legend.background = element_rect(fill = NA),
  #      legend.key.height = unit(.08, 'cm'), 
  #      legend.key.width = unit(.08, 'cm'),
  #      panel.border = element_rect(colour = "black", size = 1.0),
  #      axis.text = element_text(size = 28, colour = "black"),
  #      axis.title.y = element_text(size = 18, vjust = 1.75),
  #      axis.title.x = element_text(size = 18, vjust = -1.5)) + coord_fixed(xlim=c(-65,125))  

#VIRUS_ordplot



stressplot(group.ps.pruned.dist.ord)


#PLASMIDS

IMPUTEDmobtype_PLASMID.ps <- subset_taxa(MOBILOME_physeq, Type=="Plasmid")
IMPUTEDmobtype_PLASMID.ps
PlasmidCounts <- as.data.frame(phyloseq::otu_table(IMPUTEDmobtype_PLASMID.ps))
if(phyloseq::taxa_are_rows(IMPUTEDmobtype_PLASMID.ps) == FALSE){ otus <- t(PlasmidCounts) }
PlasmidCounts

##ALTERNATIVE FOR VOLCANO PLOTS

PlasmidCounts <- as.data.frame(phyloseq::otu_table(mobtype_PLASMID.ps))
if(phyloseq::taxa_are_rows(mobtype_PLASMID.ps) == FALSE){ otus <- t(PlasmidCounts) }
PlasmidCounts
mobtype_PLASMID.ps <-phyloseq::otu_table(mobtype_PLASMID.ps) <- phyloseq::otu_table(PLASMID_cmultRepl, taxa_are_rows = T)

return(IMPUTEDmobtype_PLASMID.ps)
mobtype_PLASMID.ps


#zCompositions component of the anlaysis, using bayesian approach:

PLASMID_cmultRepl<-cmultRepl(PlasmidCounts,  label = 0, method = c("CZM"), output = c("p-counts"), frac = 0.2, threshold = 0.5, adjust = TRUE, t= NULL, s= NULL, z.warning = 0.8, suppress.print = FALSE)

PLASMID_cmultRepl

#Now we must reconvert this OTU table to the phyloseq object

IMPUTEDmobtype_PLASMID.ps <-phyloseq::otu_table(IMPUTEDmobtype_PLASMID.ps) <- phyloseq::otu_table(PLASMID_cmultRepl, taxa_are_rows = T)

return(IMPUTEDmobtype_PLASMID.ps)
IMPUTEDmobtype_PLASMID.ps

#Create a new phyloseq object to include the imputed virus OTU table

IMPPLASMID_MOBILOME = otu_table(mobtype_PLASMID.ps, taxa_are_rows = T)
MOBGENES = phyloseq::tax_table(mobgen_mat)
MOB_SAMPLES = sample_data(mobilomesamples_df)
IMPUTED_PLASMID_physeq <- phyloseq(IMPUTEDmobtype_PLASMID.ps, MOBGENES, MOB_SAMPLES)
IMPUTED_PLASMID_physeq

#Now let's retry the ordination w/ this new phyloseq object

sum(taxa_sums(mobtype_PLASMID.ps)==0) # 0 taxa not present across all samples
PLASMID_physeq.ps.pruned <- prune_taxa(taxa_sums(IMPUTED_PLASMID_physeq) > 0.000001, IMPUTED_PLASMID_physeq)#prunning low-abundance taxa of gene groups

PLASMID_physeq.ps.pruned.css <- phyloseq_transform_css(IMPUTED_PLASMID_physeq)

PLASMID_physeq.ps.pruned.dist <- vegdist(decostand(t(otu_table(IMPUTED_PLASMID_physeq)),"rclr"), method = "euclidean")
set.seed(1999)
#PROPHAGE_physeq.ps.pruned.ord <- vegan::rda(comm = dist((PROPHAGE_physeq.ps.pruned.dist)), distance = "euclidean", k=3, try = 50, trymax = 1000, autotransform = F)

PLASMID_physeq.ps.pruned.ord<- ordinate(PLASMID_physeq.ps.pruned.css, method = "PCoA", distance = "euclidean")

PLASMID_plot<-plot_ordination(PLASMID_physeq.ps.pruned.css,PLASMID_physeq.ps.pruned.ord, type = "samples", color="CollectionPhase")
PLASMID_plot

#PROPHAGEcentroid<- PROPHAGE_plot$data %>% group_by(CollectionPhase) %>% summarise(NMDS1=mean(NMDS1), NMDS2=mean(NMDS2), .groups = "drop")


PLASMIDcentroid<- PLASMID_plot$data %>% group_by(CollectionPhase) %>% summarise(Axis.1=mean(Axis.1), Axis.2=mean(Axis.2), .groups = "drop")

#stressplot(plot_ord_all)
#str(plot_ord_all)

PLASMID_ordplot <- ggplot(PLASMID_plot$data,PLASMID_plot$mapping) + stat_ellipse(geom = "polygon", type= "t", level = 0.70, lty =1, lwd = 1.3, alpha= 0.07, show.legend =FALSE) + theme_bw() +geom_point(aes(fill=CollectionPhase), size=6, alpha=0.3, stroke=1.0,) + geom_point(data=PLASMIDcentroid, size=13, shape=16, fill="black",mapping= aes(colour=CollectionPhase), stroke=4, alpha= 0.8, show.legend = FALSE) + scale_color_manual(name="Collection phase", breaks = c("WORKDAY_START", "WORKDAY_END", "POST_SHOWER", "SWINE", "ENVIRONMENT"), values = c("#789ED0", "#F79D40", "#72BF56", "#E4655A", "#AB8BCB"), labels= c("WORKDAY_START"= "Workday start", "WORKDAY_END"="Workday end", "POST_SHOWER"="Post-shower", "SWINE"="Swine", "ENVIRONMENT"="Environment")) +
  theme(axis.line = element_line(), 
        text=element_text(family="Times",  size=10),
        panel.background = element_blank(),
        plot.margin = unit(c(0.3,0.06,0.3,0.06), "cm"),
        strip.background = element_blank(),
        strip.text = element_blank(),
        legend.key.size = unit(0.02, 'cm'),
        legend.position = c(0.8,0.8),
        legend.background = element_rect(fill = NA),
        legend.key.height = unit(.08, 'cm'), 
        legend.key.width = unit(.08, 'cm'),
        panel.border = element_rect(colour = "black", size = 1.0),
        axis.text = element_text(size = 28, colour = "black"),
        axis.title.y = element_text(size = 18, vjust = 1.75),
        axis.title.x = element_text(size = 18, vjust = -1.5)) + coord_fixed(xlim=c(-100,100), ylim = c(-50,50))  
dev.off()
PLASMID_ordplot



plasmid_metadata<- sample_data(IMPUTED_PLASMID_physeq)
pmeta<-as.data.frame(plasmid_metadata)
pmeta

plasmid.anosim<- anosim(PLASMID_physeq.ps.pruned.dist, pmeta$CollectionPhase, perm = 10000)
plasmid.anosim
#ANOSIM statistic R: 0.2398 
#Significance: 9.999e-05 
#
plasmid.adonis <- pairwise.adonis(PLASMID_physeq.ps.pruned.dist, pmeta$CollectionPhase, perm = 10000, p.adjust.m = "BH")
plasmid.adonis


#pairs Df SumsOfSqs   F.Model         R2    p.value p.adjusted sig
#1   WORKDAY_START vs POST_SHOWER  1  1278.266 0.9380651 0.04953331 0.59784022 0.59784022    
#2   WORKDAY_START vs ENVIRONMENT  1  2166.205 1.9293579 0.16173192 0.01369863 0.02283105   .
#3   WORKDAY_START vs SWINE        1  7425.544 5.6427556 0.23866743 0.00009999 0.00033330  **
#4   WORKDAY_START vs WORKDAY_END  1  2759.365 2.3733631 0.11649344 0.00119988 0.00299970   *
#5   POST_SHOWER vs ENVIRONMENT    1  2318.914 1.5447073 0.13380221 0.07539246 0.10770352    
#6   POST_SHOWER vs SWINE          1  7607.862 4.9848810 0.21687652 0.00009999 0.00033330  **
#7   POST_SHOWER vs WORKDAY_END    1  2726.597 1.9860367 0.09937121 0.00799920 0.01599840   .
#8   ENVIRONMENT vs SWINE          1  1596.413 1.1265322 0.10124738 0.30156984 0.33507760    
#9   ENVIRONMENT vs WORKDAY_END    1  1400.912 1.2276239 0.10933960 0.13968603 0.17460754    
#10  SWINE vs WORKDAY_END          1  3401.387 2.5648320 0.12471933 0.00009999 0.00033330  **

#Subset for procrustes analysis
IMPUTED_PLASMID_physeq_SUB<- subset_samples(IMPUTED_PLASMID_physeq, CollectionPhase== "POST_SHOWER" | CollectionPhase== "SWINE" | CollectionPhase== "WORKDAY_END" | CollectionPhase== "WORKDAY_START")

PLASMID_physeq.ps.pruned.dist.SUB <- vegdist(decostand(t(otu_table(IMPUTED_PLASMID_physeq_SUB)),"rclr"), method = "euclidean")

#Phages

IMPUTEDmobtype_PROPHAGE.ps <- subset_taxa(MOBILOME_physeq, Type=="Prophage")

PhageCounts <- as.data.frame(phyloseq::otu_table(IMPUTEDmobtype_PROPHAGE.ps))
if(phyloseq::taxa_are_rows(IMPUTEDmobtype_PROPHAGE.ps) == FALSE){ otus <- t(PhageCounts) }
PhageCounts

#zCompositions component of the anlaysis, using bayesian approach:


PHAGE_cmultRepl<-cmultRepl(PhageCounts,  label = 0, method = c("GBM"), output = c("p-counts"), frac = 0.2, threshold = 0.5, adjust = TRUE, t= NULL, s= NULL, z.warning = 0.8, suppress.print = FALSE)

PHAGE_cmultRepl

#Now we must reconvert this OTU table to the phyloseq object

IMPUTEDmobtype_PROPHAGE.ps <-phyloseq::otu_table(IMPUTEDmobtype_PROPHAGE.ps) <- phyloseq::otu_table(PHAGE_cmultRepl, taxa_are_rows = T)

IMPUTEDmobtype_PROPHAGE.ps

#Create a new phyloseq object to include the imputed virus OTU table

IMPPROPHAGE_MOBILOME = otu_table(IMPUTEDmobtype_PROPHAGE.ps, taxa_are_rows = T)
MOBGENES = tax_table(mobgen_mat)
MOB_SAMPLES = sample_data(mobilomesamples_df)
IMPUTED_PROPHAGE_physeq <- phyloseq(IMPUTEDmobtype_PROPHAGE.ps, MOBGENES, MOB_SAMPLES)


#Now let's retry the ordination w/ this new phyloseq object

sum(taxa_sums(IMPUTED_PROPHAGE_physeq)==0) # 0 taxa not present across all samples
PROPHAGE_physeq.ps.pruned <- prune_taxa(taxa_sums(IMPUTED_PROPHAGE_physeq) > 0.000001, IMPUTED_PROPHAGE_physeq)#prunning low-abundance taxa of gene groups

PROPHAGE_physeq.ps.pruned.css <- phyloseq_transform_css(IMPUTED_PROPHAGE_physeq)

PROPHAGE_physeq.ps.pruned.dist <- vegdist(decostand(t(otu_table(IMPUTED_PROPHAGE_physeq)),"rclr"), method = "euclidean")
set.seed(1999)
#PROPHAGE_physeq.ps.pruned.ord <- vegan::rda(comm = dist((PROPHAGE_physeq.ps.pruned.dist)), distance = "euclidean", k=3, try = 50, trymax = 1000, autotransform = F)

PROPHAGE_physeq.ps.pruned.ord<- ordinate(PROPHAGE_physeq.ps.pruned.css, method = "PCoA", distance = "euclidean")

PROPHAGE_plot<-plot_ordination(PROPHAGE_physeq.ps.pruned.css,PROPHAGE_physeq.ps.pruned.ord, type = "samples", color="CollectionPhase")
PROPHAGE_plot

#PROPHAGEcentroid<- PROPHAGE_plot$data %>% group_by(CollectionPhase) %>% summarise(NMDS1=mean(NMDS1), NMDS2=mean(NMDS2), .groups = "drop")


PROPHAGEcentroid<- PROPHAGE_plot$data %>% group_by(CollectionPhase) %>% summarise(Axis.1=mean(Axis.1), Axis.2=mean(Axis.2), .groups = "drop")

#stressplot(plot_ord_all)
#str(plot_ord_all)

PROPHAGE_ordplot <- ggplot(PROPHAGE_plot$data,PROPHAGE_plot$mapping) + stat_ellipse(geom = "polygon", type= "t", level = 0.7, lty =1, lwd = 1.3, alpha= 0.07, show.legend =FALSE) + theme_bw() +geom_point(aes(fill=CollectionPhase), size=6, alpha=0.3, stroke=1.0,) + geom_point(data=PROPHAGEcentroid, size=13, shape=16, fill="black",mapping= aes(colour=CollectionPhase), stroke=4, alpha= 0.8, show.legend = FALSE) + scale_color_manual(name="Collection phase", breaks = c("WORKDAY_START", "WORKDAY_END", "POST_SHOWER", "SWINE", "ENVIRONMENT"), values = c("#789ED0", "#F79D40", "#72BF56", "#E4655A", "#AB8BCB"), labels= c("WORKDAY_START"= "Workday start", "WORKDAY_END"="Workday end", "POST_SHOWER"="Post-shower", "SWINE"="Swine", "ENVIRONMENT"="Environment")) +
  theme(axis.line = element_line(), 
        text=element_text(family="Times",  size=10),
        panel.background = element_blank(),
        plot.margin = unit(c(0.3,0.06,0.3,0.06), "cm"),
        strip.background = element_blank(),
        strip.text = element_blank(),
        legend.key.size = unit(0.02, 'cm'),
        legend.position = c(0.8,0.8),
        legend.background = element_rect(fill = NA),
        legend.key.height = unit(.08, 'cm'), 
        legend.key.width = unit(.08, 'cm'),
        panel.border = element_rect(colour = "black", size = 1.0),
        axis.text = element_text(size = 28, colour = "black"),
        axis.title.y = element_text(size = 18, vjust = 1.75),
        axis.title.x = element_text(size = 18, vjust = -1.5)) + coord_fixed(xlim=c(-100,100), ylim = c(-50,50))  
PROPHAGE_physeq.ps.pruned.dist
PROPHAGE_ordplot

prophage_metadata<- sample_data(IMPUTED_PROPHAGE_physeq)
pmeta<-as.data.frame(prophage_metadata)
pmeta

prophage.anosim<- anosim(PROPHAGE_physeq.ps.pruned.dist, pmeta$CollectionPhase, perm = 10000)
prophage.anosim
#ANOSIM statistic R: 0.0898 
#Significance: 0.020798 

prophage.adonis <- pairwise.adonis(PROPHAGE_physeq.ps.pruned.dist, pmeta$CollectionPhase, perm = 10000, p.adjust.m = "BH")
prophage.adonis

#> prophage.adonis
#pairs Df SumsOfSqs   F.Model         R2    p.value p.adjusted sig
#1  WORKDAY_START vs POST_SHOWER  1  1649.776 0.9662134 0.05094393 0.46005399 0.57260941    
#2  WORKDAY_START vs ENVIRONMENT  1  1628.473 1.1643435 0.10429126 0.26807319 0.45812085    
#3        WORKDAY_START vs SWINE  1  8038.588 2.2997553 0.11328980 0.00019998 0.00199980   *
#4  WORKDAY_START vs WORKDAY_END  1  3378.753 2.0441223 0.10198114 0.00559944 0.01866480   .
#5    POST_SHOWER vs ENVIRONMENT  1  1901.331 1.0257512 0.09303232 0.46305369 0.57260941    
#6          POST_SHOWER vs SWINE  1  7782.356 2.0763043 0.10342064 0.00049995 0.00249975   *
#7    POST_SHOWER vs WORKDAY_END  1  2753.445 1.4448638 0.07430568 0.09929007 0.24822518    
#8          ENVIRONMENT vs SWINE  1  2498.817 0.4926787 0.04695452 0.64743526 0.64743526    
#9    ENVIRONMENT vs WORKDAY_END  1  1518.876 0.8652582 0.07963531 0.51534847 0.57260941    
#10         SWINE vs WORKDAY_END  1  4167.731 1.1283585 0.05898878 0.27487251 0.45812085    



#Subset for procrustes analysis
IMPUTED_PROPHAGE_physeq_SUB<- subset_samples(IMPUTED_PROPHAGE_physeq, CollectionPhase== "POST_SHOWER" | CollectionPhase== "SWINE" | CollectionPhase== "WORKDAY_END" | CollectionPhase== "WORKDAY_START")

PROPHAGE_physeq.ps.pruned.dist.SUB <- vegdist(decostand(t(otu_table(IMPUTED_PROPHAGE_physeq_SUB)),"rclr"), method = "euclidean")

#ICE

IMPUTEDmobtype_ICE.ps <- subset_taxa(MOBILOME_physeq, Type=="ICE")
IMPUTEDmobtype_ICE.ps.pruned<-subset_samples(IMPUTEDmobtype_ICE.ps, sample_names(IMPUTEDmobtype_ICE.ps) != "R1_19")

ICECounts <- as.data.frame(phyloseq::otu_table(IMPUTEDmobtype_ICE.ps.pruned))
if(phyloseq::taxa_are_rows(IMPUTEDmobtype_ICE.ps.pruned) == FALSE){ otus <- t(ICECounts) }
ICECounts

#zCompositions component of the anlaysis, using bayesian approach:


ICE_cmultRepl<-cmultRepl(ICECounts,  label = 0, method = c("GBM"), output = c("p-counts"), frac = 0.4, threshold = 0.5, adjust = TRUE, t= NULL, s= NULL, z.warning = 0.8, suppress.print = FALSE)

ICE_cmultRepl

#Now we must reconvert this OTU table to the phyloseq object

IMPUTEDmobtype_ICE.ps.pruned <-phyloseq::otu_table(IMPUTEDmobtype_ICE.ps.pruned) <- phyloseq::otu_table(ICE_cmultRepl, taxa_are_rows = T)

IMPUTEDmobtype_ICE.ps.pruned

#Create a new phyloseq object to include the imputed virus OTU table

IMPICE_MOBILOME = otu_table(IMPUTEDmobtype_ICE.ps.pruned, taxa_are_rows = T)
MOBGENES = tax_table(mobgen_mat)
MOB_SAMPLES = sample_data(mobilomesamples_df)
IMPUTED_ICE_physeq <- phyloseq(IMPICE_MOBILOME, MOBGENES, MOB_SAMPLES)

IMPUTED_ICE_physeq
#Now let's retry the ordination w/ this new phyloseq object

sum(taxa_sums(IMPUTED_ICE_physeq)==0) # 0 taxa not present across all samples
#PROPHAGE_physeq.ps.pruned <- prune_taxa(taxa_sums(IMPUTED_ICE_physeq) > 0.000001, IMPUTED_PROPHAGE_physeq)#prunning low-abundance taxa of gene groups

ICE_physeq.ps.pruned.css <- phyloseq_transform_css(IMPUTED_ICE_physeq)

ICE_physeq.ps.pruned.dist <- vegdist(decostand(t(otu_table(IMPUTED_ICE_physeq)),"rclr"), method = "euclidean")
set.seed(1999)
#PROPHAGE_physeq.ps.pruned.ord <- vegan::rda(comm = dist((PROPHAGE_physeq.ps.pruned.dist)), distance = "euclidean", k=3, try = 50, trymax = 1000, autotransform = F)

ICE_physeq.ps.pruned.ord<- ordinate(ICE_physeq.ps.pruned.css, method = "PCoA", distance = "euclidean")

ICE_plot<-plot_ordination(ICE_physeq.ps.pruned.css,ICE_physeq.ps.pruned.ord, type = "samples", color="CollectionPhase")
ICE_plot

#PROPHAGEcentroid<- PROPHAGE_plot$data %>% group_by(CollectionPhase) %>% summarise(NMDS1=mean(NMDS1), NMDS2=mean(NMDS2), .groups = "drop")


ICEcentroid<- ICE_plot$data %>% group_by(CollectionPhase) %>% summarise(Axis.1=mean(Axis.1), Axis.2=mean(Axis.2), .groups = "drop")

#stressplot(plot_ord_all)
#str(plot_ord_all)

ICE_ordplot <- ggplot(ICE_plot$data,ICE_plot$mapping) + stat_ellipse(geom = "polygon", type= "t", level = 0.7, lty =1, lwd = 1.3, alpha= 0.07, show.legend =FALSE) + theme_bw() +geom_point(aes(fill=CollectionPhase), size=6, alpha=0.3, stroke=1.0,) + geom_point(data=ICEcentroid, size=13, shape=16, fill="black",mapping= aes(colour=CollectionPhase), stroke=4, alpha= 0.8, show.legend = FALSE) + scale_color_manual(name="Collection phase", breaks = c("WORKDAY_START", "WORKDAY_END", "POST_SHOWER", "SWINE", "ENVIRONMENT"), values = c("#789ED0", "#F79D40", "#72BF56", "#E4655A", "#AB8BCB"), labels= c("WORKDAY_START"= "Workday start", "WORKDAY_END"="Workday end", "POST_SHOWER"="Post-shower", "SWINE"="Swine", "ENVIRONMENT"="Environment")) +
  theme(axis.line = element_line(), 
        text=element_text(family="Times",  size=10),
        panel.background = element_blank(),
        plot.margin = unit(c(0.3,0.06,0.3,0.06), "cm"),
        strip.background = element_blank(),
        strip.text = element_blank(),
        legend.key.size = unit(0.02, 'cm'),
        legend.position = c(0.8,0.8),
        legend.background = element_rect(fill = NA),
        legend.key.height = unit(.08, 'cm'), 
        legend.key.width = unit(.08, 'cm'),
        panel.border = element_rect(colour = "black", size = 1.0),
        axis.text = element_text(size = 28, colour = "black"),
        axis.title.y = element_text(size = 18, vjust = 1.75),
        axis.title.x = element_text(size = 18, vjust = -1.5)) + coord_fixed(xlim=c(-100,100), ylim = c(-50,50))  
ICE_ordplot

ICE_metadata<- sample_data(IMPUTED_ICE_physeq)
pmeta<-as.data.frame(ICE_metadata)
pmeta

ice.anosim<- anosim(ICE_physeq.ps.pruned.dist, pmeta$CollectionPhase, perm = 10000)
ice.anosim
#ANOSIM statistic R: 0.4434 
#Significance: 9.999e-05 

ice.adonis <- pairwise.adonis(ICE_physeq.ps.pruned.dist, pmeta$CollectionPhase, perm = 10000, p.adjust.m = "BH")
ice.adonis

#> ice.adonis
#pairs Df SumsOfSqs  F.Model         R2    p.value p.adjusted sig
#1  WORKDAY_START vs POST_SHOWER  1  1392.063 1.472582 0.07562337 0.07669233 0.07669233    
#2  WORKDAY_START vs ENVIRONMENT  1  2363.576 2.041933 0.16956852 0.01549845 0.01937306   .
#3        WORKDAY_START vs SWINE  1  6454.564 6.176066 0.25546199 0.00009999 0.00033330  **
#4  WORKDAY_START vs WORKDAY_END  1  4302.941 4.186186 0.19759035 0.00019998 0.00049995  **
#5    POST_SHOWER vs ENVIRONMENT  1  2357.899 3.421228 0.25491171 0.01389861 0.01937306   .
#6          POST_SHOWER vs SWINE  1  4736.633 6.034588 0.25107932 0.00009999 0.00033330  **
#7    POST_SHOWER vs WORKDAY_END  1  3189.419 4.238954 0.19958394 0.00009999 0.00033330  **
#8          ENVIRONMENT vs SWINE  1  1862.764 2.144098 0.17655471 0.01469853 0.01937306   .
#9    ENVIRONMENT vs WORKDAY_END  1  2956.385 3.619900 0.28684062 0.01809819 0.02010910   .
#10         SWINE vs WORKDAY_END  1  2808.130 3.272695 0.16143363 0.00129987 0.00259974   *



#Transposable elements (TEs)

mobtype_TE.ps <- subset_taxa(MOBILOME_physeq, Type=="TE")
mobtype_TE.ps.pruned<-subset_samples(mobtype_TE.ps, sample_names(mobtype_TE.ps) != "R1_19")
mobtype_TE.ps.pruned


TECounts <- as.data.frame(phyloseq::otu_table(mobtype_TE.ps.pruned))
if(phyloseq::taxa_are_rows(mobtype_TE.ps.pruned) == FALSE){ otus <- t(TECounts) }
TECounts

#zCompositions component of the anlaysis, using bayesian approach:


TE_cmultRepl<-cmultRepl(TECounts, label = 0, method = c("GBM"), output = c("p-counts"), frac = 0.3, threshold = 0.5, adjust = TRUE, t= NULL, s= NULL, z.warning = 0.8, suppress.print = FALSE)

TE_cmultRepl

#Now we must reconvert this OTU table to the phyloseq object

mobtype_TE.ps.pruned <-phyloseq::otu_table(mobtype_TE.ps.pruned) <- phyloseq::otu_table(TE_cmultRepl, taxa_are_rows = T)
return(IMPUTEDmobtype_TE.ps.pruned)
mobtype_TE.ps.pruned

#Create a new phyloseq object to include the imputed virus OTU table

IMPTE_MOBILOME = otu_table(mobtype_TE.ps.pruned, taxa_are_rows = T)
MOBGENES = tax_table(mobgen_mat)
MOB_SAMPLES = sample_data(mobilomesamples_df)
IMPUTED_TE_physeq <- phyloseq(IMPTE_MOBILOME, MOBGENES, MOB_SAMPLES)

IMPUTED_TE_physeq
#Now let's retry the ordination w/ this new phyloseq object

sum(taxa_sums(IMPUTED_TE_physeq)==0) # 0 taxa not present across all samples
#PROPHAGE_physeq.ps.pruned <- prune_taxa(taxa_sums(IMPUTED_ICE_physeq) > 0.000001, IMPUTED_PROPHAGE_physeq)#prunning low-abundance taxa of gene groups

TE_physeq.ps.pruned.css <- phyloseq_transform_css(IMPUTED_TE_physeq)

TE_physeq.ps.pruned.dist <- vegdist(decostand(t(otu_table(IMPUTED_TE_physeq)),"rclr"), method = "euclidean")
set.seed(1999)
#PROPHAGE_physeq.ps.pruned.ord <- vegan::rda(comm = dist((PROPHAGE_physeq.ps.pruned.dist)), distance = "euclidean", k=3, try = 50, trymax = 1000, autotransform = F)

TE_physeq.ps.pruned.ord<- ordinate(TE_physeq.ps.pruned.css, method = "PCoA", distance = "euclidean")

TE_plot<-plot_ordination(TE_physeq.ps.pruned.css,TE_physeq.ps.pruned.ord, type = "samples", color="CollectionPhase")
TE_plot

#PROPHAGEcentroid<- PROPHAGE_plot$data %>% group_by(CollectionPhase) %>% summarise(NMDS1=mean(NMDS1), NMDS2=mean(NMDS2), .groups = "drop")


TEcentroid<- TE_plot$data %>% group_by(CollectionPhase) %>% summarise(Axis.1=mean(Axis.1), Axis.2=mean(Axis.2), .groups = "drop")

#stressplot(plot_ord_all)
#str(plot_ord_all)

TE_ordplot <- ggplot(TE_plot$data,TE_plot$mapping) + stat_ellipse(geom = "polygon", type= "t", level = 0.85, lty =1, lwd = 1.3, alpha= 0.07, show.legend =FALSE) + theme_bw() +geom_point(aes(fill=CollectionPhase), size=6, alpha=0.3, stroke=1.0,) + geom_point(data=TEcentroid, size=13, shape=16, fill="black",mapping= aes(colour=CollectionPhase), stroke=4, alpha= 0.8, show.legend = FALSE) + scale_color_manual(name="Collection phase", breaks = c("WORKDAY_START", "WORKDAY_END", "POST_SHOWER", "SWINE", "ENVIRONMENT"), values = c("#789ED0", "#F79D40", "#72BF56", "#E4655A", "#AB8BCB"), labels= c("WORKDAY_START"= "Workday start", "WORKDAY_END"="Workday end", "POST_SHOWER"="Post-shower", "SWINE"="Swine", "ENVIRONMENT"="Environment")) +
  theme(axis.line = element_line(), 
        text=element_text(family="Times",  size=10),
        panel.background = element_blank(),
        plot.margin = unit(c(0.3,0.06,0.3,0.06), "cm"),
        strip.background = element_blank(),
        strip.text = element_blank(),
        legend.key.size = unit(0.02, 'cm'),
        legend.position = c(0.8,0.8),
        legend.background = element_rect(fill = NA),
        legend.key.height = unit(.08, 'cm'), 
        legend.key.width = unit(.08, 'cm'),
        panel.border = element_rect(colour = "black", size = 1.0),
        axis.text = element_text(size = 28, colour = "black"),
        axis.title.y = element_text(size = 18, vjust = 1.75),
        axis.title.x = element_text(size = 18, vjust = -1.5)) + coord_fixed(xlim=c(-100,100), ylim = c(-50,50))  
TE_ordplot

TE_metadata<- sample_data(IMPUTED_TE_physeq)
pmeta<-as.data.frame(TE_metadata)
pmeta


TE.anosim<- anosim(TE_physeq.ps.pruned.dist, pmeta$CollectionPhase, perm = 1000)
TE.anosim
#ANOSIM statistic R: 0.4024 
#Significance: 9.999e-05 

te.adonis <- pairwise.adonis(TE_physeq.ps.pruned.dist, pmeta$CollectionPhase, perm = 1000, p.adjust.m = "BH")
te.adonis

#pairs Df SumsOfSqs   F.Model         R2     p.value  p.adjusted sig
#1  WORKDAY_START vs POST_SHOWER  1  1310.674 1.6313133 0.08309751 0.027972028 0.037462537   .
#2  WORKDAY_START vs ENVIRONMENT  1   742.273 0.9936652 0.09038526 0.366633367 0.366633367    
#3        WORKDAY_START vs SWINE  1  4248.060 4.9858079 0.21690810 0.000999001 0.002497502   *
#4  WORKDAY_START vs WORKDAY_END  1  1683.098 2.1348131 0.11156697 0.000999001 0.002497502   *
#5    POST_SHOWER vs ENVIRONMENT  1  1133.879 1.4789861 0.12884292 0.083916084 0.093240093    
#6          POST_SHOWER vs SWINE  1  4235.173 4.9077860 0.21424096 0.000999001 0.002497502   *
#7    POST_SHOWER vs WORKDAY_END  1  1425.044 1.7813783 0.09484811 0.002997003 0.005994006   *
#8          ENVIRONMENT vs SWINE  1  1505.068 1.7621490 0.14981522 0.029970030 0.037462537   .
#9    ENVIRONMENT vs WORKDAY_END  1  1235.704 1.6831535 0.15755212 0.015984016 0.026640027   .
#10         SWINE vs WORKDAY_END  1  2501.970 2.9386268 0.14738361 0.000999001 0.002497502   *



#Insertional sequences (ISs)

mobtype_IS.ps <- subset_taxa(MOBILOME_physeq, Type=="IS")
mobtype_IS.ps.pruned<-subset_samples(mobtype_IS.ps, sample_names(mobtype_IS.ps) != "R1_19")
mobtype_IS.ps.pruned


ISCounts <- as.data.frame(phyloseq::otu_table(mobtype_IS.ps.pruned))
if(phyloseq::taxa_are_rows(mobtype_IS.ps.pruned) == FALSE){ otus <- t(ISCounts) }
ISCounts
ISCounts<-ISCounts[rowSums(ISCounts[])>0,]#remove any rows (taxa) with all zeros
ISCounts
#zCompositions component of the anlaysis, using bayesian approach:


IS_cmultRepl<-cmultRepl(ISCounts, label = 0, method = c("GBM"), output = c("p-counts"), frac = 0.3, threshold = 0.5, adjust = TRUE, t= NULL, s= NULL, z.warning = 0.8, suppress.print = FALSE)

IS_cmultRepl

#Now we must reconvert this OTU table to the phyloseq object

mobtype_IS.ps.pruned <-phyloseq::otu_table(mobtype_IS.ps.pruned) <- phyloseq::otu_table(IS_cmultRepl, taxa_are_rows = T)
return(mobtype_IS.ps.pruned)
mobtype_IS.ps.pruned

#Create a new phyloseq object to include the imputed virus OTU table

IMPUTED_IS_MOBILOME = otu_table(mobtype_IS.ps.pruned, taxa_are_rows = T)
MOBGENES = tax_table(mobgen_mat)
MOB_SAMPLES = sample_data(mobilomesamples_df, CollectionPhase== "POST_SHOWER" | CollectionPhase== "SWINE" | CollectionPhase== "WORKDAY_END" | CollectionPhase== "WORKDAY_START"| CollectionPhase== "ENVIRONMENT")
IMPUTED_IS_physeq <- phyloseq(IMPUTED_IS_MOBILOME, MOBGENES, MOB_SAMPLES)

IMPUTED_IS_physeq
#Now let's retry the ordination w/ this new phyloseq object

sum(taxa_sums(IMPUTED_IS_physeq)==0) # 0 taxa not present across all samples
#PROPHAGE_physeq.ps.pruned <- prune_taxa(taxa_sums(IMPUTED_ICE_physeq) > 0.000001, IMPUTED_PROPHAGE_physeq)#prunning low-abundance taxa of gene groups

IS_physeq.ps.pruned.css <- phyloseq_transform_css(IMPUTED_IS_physeq)

IS_physeq.ps.pruned.dist <- vegdist(decostand(t(otu_table(IMPUTED_IS_physeq)),"rclr"), method = "euclidean")
set.seed(1999)
#PROPHAGE_physeq.ps.pruned.ord <- vegan::rda(comm = dist((PROPHAGE_physeq.ps.pruned.dist)), distance = "euclidean", k=3, try = 50, trymax = 1000, autotransform = F)

IS_physeq.ps.pruned.ord<- ordinate(IS_physeq.ps.pruned.css, method = "PCoA", distance = "euclidean")

IS_plot<-plot_ordination(IS_physeq.ps.pruned.css,IS_physeq.ps.pruned.ord, type = "samples", color="CollectionPhase")
IS_plot

#PROPHAGEcentroid<- PROPHAGE_plot$data %>% group_by(CollectionPhase) %>% summarise(NMDS1=mean(NMDS1), NMDS2=mean(NMDS2), .groups = "drop")


IScentroid<- IS_plot$data %>% group_by(CollectionPhase) %>% summarise(Axis.1=mean(Axis.1), Axis.2=mean(Axis.2), .groups = "drop")

#stressplot(plot_ord_all)
#str(plot_ord_all)

IS_ordplot <- ggplot(IS_plot$data,IS_plot$mapping) + stat_ellipse(geom = "polygon", type= "t", level = 0.85, lty =1, lwd = 1.3, alpha= 0.07, show.legend =FALSE) + theme_bw() +geom_point(aes(fill=CollectionPhase), size=6, alpha=0.3, stroke=1.0,) + geom_point(data=IScentroid, size=13, shape=16, fill="black",mapping= aes(colour=CollectionPhase), stroke=4, alpha= 0.8, show.legend = FALSE) + scale_color_manual(name="Collection phase", breaks = c("WORKDAY_START", "WORKDAY_END", "POST_SHOWER", "SWINE", "ENVIRONMENT"), values = c("#789ED0", "#F79D40", "#72BF56", "#E4655A", "#AB8BCB"), labels= c("WORKDAY_START"= "Workday start", "WORKDAY_END"="Workday end", "POST_SHOWER"="Post-shower", "SWINE"="Swine", "ENVIRONMENT"="Environment")) +
  theme(axis.line = element_line(), 
        text=element_text(family="Times",  size=10),
        panel.background = element_blank(),
        plot.margin = unit(c(0.3,0.06,0.3,0.06), "cm"),
        strip.background = element_blank(),
        strip.text = element_blank(),
        legend.key.size = unit(0.02, 'cm'),
        legend.position = c(0.8,0.8),
        legend.background = element_rect(fill = NA),
        legend.key.height = unit(.08, 'cm'), 
        legend.key.width = unit(.08, 'cm'),
        panel.border = element_rect(colour = "black", size = 1.0),
        axis.text = element_text(size = 28, colour = "black"),
        axis.title.y = element_text(size = 18, vjust = 1.75),
        axis.title.x = element_text(size = 18, vjust = -1.5)) + coord_fixed(xlim=c(-100,100), ylim = c(-50,50))  
IS_ordplot

IS_metadata<- sample_data(IMPUTED_IS_physeq)
pmeta<-as.data.frame(IS_metadata)
pmeta


is.anosim<- anosim(IS_physeq.ps.pruned.dist, pmeta$CollectionPhase, perm = 1000)
is.anosim
#ANOSIM statistic R: 0.2393 
#Significance: 0.000999 

is.adonis <- pairwise.adonis(IS_physeq.ps.pruned.dist, pmeta$CollectionPhase, perm = 1000, p.adjust.m = "BH")
is.adonis
#pairs Df SumsOfSqs  F.Model         R2     p.value  p.adjusted sig
#1 WORKDAY_START vs POST_SHOWER  1  4604.617 1.449184 0.07451128 0.047952048 0.047952048   .
#2       WORKDAY_START vs SWINE  1 18768.346 4.933650 0.21512712 0.000999001 0.002997003   *
#  3 WORKDAY_START vs WORKDAY_END  1  8935.919 2.748624 0.13918055 0.001998002 0.002997003   *
#  4         POST_SHOWER vs SWINE  1 18044.307 4.858849 0.21255877 0.000999001 0.002997003   *
#  5   POST_SHOWER vs WORKDAY_END  1  7062.841 2.238420 0.11635154 0.003996004 0.004795205   *
#  6         SWINE vs WORKDAY_END  1  8976.207 2.350462 0.12146799 0.001998002 0.002997003   *

#############################################################################################
#DIFFERENTIAL ABUNDANCE TESTING##############################################################
#############################################################################################
##Performing initial imputation via zCompositions using bayesian approach for entire mobilome:


ALLMOBCounts <- as.data.frame(phyloseq::otu_table(MOBILOME_physeq))
if(phyloseq::taxa_are_rows(MOBILOME_physeq) == FALSE){ otus <- t(ALLMOBCounts) }
ALLMOBCounts

##There are rows in the count matrix where, after filtering and subsetting, row sums are all zeros, and thus must be removed. To identify which these are:
# Compute row sums
#row_sums <- rowSums(IMPUTED_ALLMOBCounts)
# Find row indices with sum 0
#rows_with_zero_sum <- which(row_sums == 0)
# Extract row names for rows with sum 0
#rows_with_zero_sum_names <- rownames(IMPUTED_ALLMOBCounts)[rows_with_zero_sum]
#print(rows_with_zero_sum_names)
# Print row names with sum 0
#cat("Rows with sum 0:", paste(rows_with_zero_sum_names, collapse = ", "))

ALLMOB_cmultRepl<-cmultRepl(ALLMOBCounts,  label = 0, method = c("GBM"), output = c("p-counts"), frac = 0.65, threshold = 0.5, adjust = TRUE, t= NULL, s= NULL, z.warning = 0.8, suppress.print = FALSE)
ALLMOB_cmultRepl

#Now we must reconvert this OTU table to the phyloseq object

MOBILOME_physeq <-phyloseq::otu_table(MOBILOME_physeq) <- phyloseq::otu_table(ALLMOB_cmultRepl, taxa_are_rows = T)

return(MOBILOME_physeq)

#Create a new phyloseq object to include the imputed 

IMPALL_MOBILOME = otu_table(MOBILOME_physeq, taxa_are_rows = T)
MOBGENES = tax_table(mobgen_mat)
MOB_SAMPLES = sample_data(mobilomesamples_df)
IMPUTED_ALLMOBILOME_physeq <- phyloseq(IMPALL_MOBILOME, MOBGENES, MOB_SAMPLES)
IMPUTED_ALLMOBILOME_physeq

sample_data(IMPUTED_ALLMOBILOME_physeq)$CollectionPhase <- as.factor(sample_data(IMPUTED_ALLMOBILOME_physeq)$CollectionPhase) # factorize for DESeq2

sample_data(IMPUTED_ALLMOBILOME_physeq)$CollectionPhase <- as.factor(sample_data(IMPUTED_ALLMOBILOME_physeq)$CollectionPhase) # factorize for DESeq2
sample_data(IMPUTED_ALLMOBILOME_physeq)$EVSMOK <- as.factor(sample_data(IMPUTED_ALLMOBILOME_physeq)$EVSMOK) # factorize for DESeq2
sample_data(IMPUTED_ALLMOBILOME_physeq)$Shotgun_Post_Host_Rem  <- as.numeric(sample_data(IMPUTED_ALLMOBILOME_physeq)$Shotgun_Post_Host_Rem ) # factorize for DESeq2
sample_data(IMPUTED_ALLMOBILOME_physeq)$Shotgun_Post_Host_Rem  <- scale(sample_data(IMPUTED_ALLMOBILOME_physeq)$Shotgun_Post_Host_Rem , center = TRUE, scale=TRUE) # factorize for DESeq2
sample_data(IMPUTED_ALLMOBILOME_physeq)$Post_captureSampleConc_ng_ul    <- as.numeric(sample_data(IMPUTED_ALLMOBILOME_physeq)$Post_captureSampleConc_ng_ul) # factorize for DESeq2
sample_data(IMPUTED_ALLMOBILOME_physeq)$Post_captureSampleConc_ng_ul  <- scale(sample_data(IMPUTED_ALLMOBILOME_physeq)$Post_captureSampleConc_ng_ul , center = TRUE, scale=TRUE) # factorize for DESeq2
sample_data(IMPUTED_ALLMOBILOME_physeq)$GEND   <- as.factor(sample_data(IMPUTED_ALLMOBILOME_physeq)$GEND) # factorize for DESeq2
sample_data(IMPUTED_ALLMOBILOME_physeq)$AGE   <- scale(sample_data(IMPUTED_ALLMOBILOME_physeq)$AGE  , center = TRUE, scale=TRUE) # factorize for DESeq2

sample_data(IMPUTED_ALLMOBILOME_physeq)$BMI   <- scale(sample_data(IMPUTED_ALLMOBILOME_physeq)$BMI  , center = TRUE, scale=TRUE) # factorize for DESeq2
sample_data(IMPUTED_ALLMOBILOME_physeq)$EATPRK   <- scale(sample_data(IMPUTED_ALLMOBILOME_physeq)$EATPRK  , center = TRUE, scale=TRUE) # factorize for DESeq2
sample_data(IMPUTED_ALLMOBILOME_physeq)$HNDPRK   <- scale(sample_data(IMPUTED_ALLMOBILOME_physeq)$HNDPRK  , center = TRUE, scale=TRUE) # factorize for DESeq2


##############SUBSET FOR WORKSTART_WORK_END


filtered_IMPUTED_ALLMOBILOME_physeq<-taxa_filter(IMPUTED_ALLMOBILOME_physeq, frequency = 0.8, treatment = 'CollectionPhase')
filtered_IMPUTED_ALLMOBILOME_physeq

IMPUTED_ALL_MOBILOME <- subset_samples(filtered_IMPUTED_ALLMOBILOME_physeq, CollectionPhase=="WORKDAY_START" | CollectionPhase=="WORKDAY_END")
IMPUTED_ALL_MOBILOME
#Filter sparse features of IMPUTED_ALLMOBILOME_physeq object so that we only consider those features with less than 90% zeros
#IMPUTED_ALL_MOBILOME
#nosparse_ALLMOBILOME<- prune_taxa(rowSums(otu_table(IMPUTED_ALL_MOBILOME)==0) < ncol(otu_table(IMPUTED_ALL_MOBILOME))*0.9, IMPUTED_ALL_MOBILOME)

nosparse_ALLMOBILOME<-filter_taxa(IMPUTED_ALL_MOBILOME, function(x){(sum(x > 20) > nsamples(IMPUTED_ALL_MOBILOME)*0.1)}, prune = TRUE)
nosparse_ALLMOBILOME

#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 2208 taxa and 20 samples ]
#sample_data() Sample Data:       [ 20 samples by 93 sample variables ]
#tax_table()   Taxonomy Table:    [ 2208 taxa by 10 taxonomic ranks ]

DS_ALLMOBILOME= phyloseq_to_deseq2(nosparse_ALLMOBILOME, ~CollectionPhase+ Post_captureSampleConc_ng_ul + Shotgun_Post_Host_Rem + GEND + BMI + EVSMOK + EATPRK)
DS_ALLMGES<- estimateSizeFactors(DS_ALLMOBILOME, type="poscounts")
DS_ALLMGES$CollectionPhase <- relevel(DS_ALLMGES$CollectionPhase, ref = "WORKDAY_START")
DS_ALLMGES= DESeq(DS_ALLMGES, test= "Wald", fitType= "parametric") 
alpha= 0.01
DS_ALLMGES_RESULT= results(DS_ALLMGES, alpha = alpha)
DS_ALLMGES <- lfcShrink(DS_ALLMGES, res= DS_ALLMGES_RESULT, type="ashr")
DS_ALLMGES_RESULT= DS_ALLMGES_RESULT[order(DS_ALLMGES_RESULT$padj, na.last = NA), ]
DS_ALLMGES_RESULT
summary(DS_ALLMGES_RESULT)
#out of 2463 with nonzero total read count
#adjusted p-value < 0.01
#LFC > 0 (up)       : 115, 4.7%
#LFC < 0 (down)     : 61, 2.5%
#out of 1996 with nonzero total read count
#adjusted p-value < 0.01
#LFC > 0 (up)       : 65, 3.3%
#LFC < 0 (down)     : 12, 0.6%

R<-as.data.frame(DS_ALLMGES_RESULT)
R_SIG_ALLMGES<- head(R[, 1:6], 77)
R_SIG_ALLMGES<- R_SIG_ALLMGES[order(R_SIG_ALLMGES$baseMean, na.last = NA), ]
write.csv(R_SIG_ALLMGES, "C:/Users/is233/OneDrive/Desktop/AG_WORKER_METAGENOMICS_STUDY/MANUSCRIPT/ANALYSIS/WORKER_MOBILOME/DESEQ_ALLMGE_WSWE.csv", row.names = TRUE)
DS_ALLARGS <- estimateSizeFactors(DS_ALLARGS)
DS_ALLARGS <- estimateDispersions(DS_ALLARGS)
plotDispEsts(DS_ALLARGS)
par(mfrow=c(1,1))
# Make a basic volcano plot
with(R, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot", xlim=c(-50,50)))
with(subset(R, padj<.01 ), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
with(subset(R, padj<.01 & abs(log2FoldChange)>1.5), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))
EnhancedVolcano(R,
                lab = rownames(R),
                x = 'log2FoldChange',
                y = 'padj',  pCutoff = 10e-3,
                FCcutoff = 1.5,
                pointSize = log10(R$baseMean)+2,
                labSize = 5,
                shape=16,
                colAlpha = 0.4, legendLabels=c('Not sig.','Log (base 2) FC','FDRadj p <0.01',
                                               'FDRadj p <0.01 & Log (base 2) FC'),
                legendPosition = 'right',
                legendLabSize = 12,
                legendIconSize = 4.0)

##ggplot2 Volcanoplot

R$diffexpressed<- 'NO'
R$diffexpressed[R$log2FoldChange>1.5 & R$padj<0.01]<- 'Workday end_T2'
R$diffexpressed[R$log2FoldChange<(-1.5) & R$padj<0.01]<- 'Workday start_T1'

ggplot(data=R, aes(x=log2FoldChange, y=-log10(padj), col=diffexpressed))+ geom_vline(xintercept = c(-1.5, 1.5), col='gray', linetype='dashed') + geom_hline(yintercept= c(-log10(10e-3)), col='grey', linetype= 'dashed') + geom_point(size=log10(R$baseMean), shape=16, stroke=2)+
  
  theme_classic(base_size = 20)+ theme(
    axis.title.y = element_text(face = "bold", margin = margin(0,20,0,0), size = rel(1.1)), 
    axis.title.x = element_text(hjust = 0.5, face = "bold", margin = margin(20,0,0,0), size = rel(1.1)), 
    plot.title = element_text(hjust = 0.5))

sum( rowMeans( counts(DS_ALLMGES, normalized=TRUE)) > 5 )#check number of rows for nsub
vsdata <- vst(DS_ALLMGES, blind=TRUE, nsub=1949)
DESEQ_PCA_MGE_GENES_WSWE<-plotPCA(vsdata, intgroup="CollectionPhase", ntop=1000) + theme_bw()
DESEQ_PCA_MGE_GENES_WSWE

DS_SIG_ALLMGES<- rownames(DS_ALLMGES_RESULT[1:200, ]) #Selecting the bottom 100 with the lowest adjusted p values
DS_SIG_ALLMGES
#total=median(sample_sums(IMPUTED_ALLRESISTOME_physeq))
#standf=function(x, t=total) round(t * (x / sum(x)))
DS_REL_ALLMGES<-phyloseq_transform_css(IMPUTED_ALL_MOBILOME)
#DS_REL_ALLARGS<- transform_sample_counts(IMPUTED_ALLRESISTOME_physeq, function(x) x/sum(x)*10)
DS_SIG_REL_ALLMGES<- prune_taxa(DS_SIG_ALLMGES, DS_REL_ALLMGES)

DS_SIG_REL_ALLMGES


##############SUBSET FOR WORK_END_POST_SHOWER

IMPALL_MOBILOME = otu_table(MOBILOME_physeq, taxa_are_rows = T)
MOBGENES = tax_table(mobgen_mat)
MOB_SAMPLES = sample_data(mobilomesamples_df)
IMPUTED_ALLMOBILOME_physeq <- phyloseq(IMPALL_MOBILOME, MOBGENES, MOB_SAMPLES)
IMPUTED_ALLMOBILOME_physeq


sample_data(IMPUTED_ALLMOBILOME_physeq)$CollectionPhase <- as.factor(sample_data(IMPUTED_ALLMOBILOME_physeq)$CollectionPhase) # factorize for DESeq2

sample_data(IMPUTED_ALLMOBILOME_physeq)$CollectionPhase <- as.factor(sample_data(IMPUTED_ALLMOBILOME_physeq)$CollectionPhase) # factorize for DESeq2
sample_data(IMPUTED_ALLMOBILOME_physeq)$EVSMOK <- as.factor(sample_data(IMPUTED_ALLMOBILOME_physeq)$EVSMOK) # factorize for DESeq2
sample_data(IMPUTED_ALLMOBILOME_physeq)$Shotgun_Post_Host_Rem  <- as.numeric(sample_data(IMPUTED_ALLMOBILOME_physeq)$Shotgun_Post_Host_Rem ) # factorize for DESeq2
sample_data(IMPUTED_ALLMOBILOME_physeq)$Shotgun_Post_Host_Rem  <- scale(sample_data(IMPUTED_ALLMOBILOME_physeq)$Shotgun_Post_Host_Rem , center = TRUE, scale=TRUE) # factorize for DESeq2
sample_data(IMPUTED_ALLMOBILOME_physeq)$Post_captureSampleConc_ng_ul    <- as.numeric(sample_data(IMPUTED_ALLMOBILOME_physeq)$Post_captureSampleConc_ng_ul) # factorize for DESeq2
sample_data(IMPUTED_ALLMOBILOME_physeq)$Post_captureSampleConc_ng_ul  <- scale(sample_data(IMPUTED_ALLMOBILOME_physeq)$Post_captureSampleConc_ng_ul , center = TRUE, scale=TRUE) # factorize for DESeq2
sample_data(IMPUTED_ALLMOBILOME_physeq)$GEND   <- as.factor(sample_data(IMPUTED_ALLMOBILOME_physeq)$GEND) # factorize for DESeq2
sample_data(IMPUTED_ALLMOBILOME_physeq)$AGE   <- scale(sample_data(IMPUTED_ALLMOBILOME_physeq)$AGE  , center = TRUE, scale=TRUE) # factorize for DESeq2

sample_data(IMPUTED_ALLMOBILOME_physeq)$BMI   <- scale(sample_data(IMPUTED_ALLMOBILOME_physeq)$BMI  , center = TRUE, scale=TRUE) # factorize for DESeq2

sample_data(IMPUTED_ALLMOBILOME_physeq)$EATPRK   <- scale(sample_data(IMPUTED_ALLMOBILOME_physeq)$EATPRK  , center = TRUE, scale=TRUE) # factorize for DESeq2
sample_data(IMPUTED_ALLMOBILOME_physeq)$HNDPRK   <- scale(sample_data(IMPUTED_ALLMOBILOME_physeq)$HNDPRK  , center = TRUE, scale=TRUE) # factorize for DESeq2


filtered_IMPUTED_ALLMOBILOME_physeq<-taxa_filter(IMPUTED_ALLMOBILOME_physeq, frequency = 0.8, treatment = 'CollectionPhase')
filtered_IMPUTED_ALLMOBILOME_physeq

IMPUTED_ALL_MOBILOME <- subset_samples(filtered_IMPUTED_ALLMOBILOME_physeq, CollectionPhase=="POST_SHOWER" | CollectionPhase=="WORKDAY_END")
IMPUTED_ALL_MOBILOME
#Filter sparse features of IMPUTED_ALLMOBILOME_physeq object so that we only consider those features with less than 90% zeros
#nosparse_ALLMOBILOME<- prune_taxa(rowSums(otu_table(IMPUTED_ALL_MOBILOME)==0) < ncol(otu_table(IMPUTED_ALL_MOBILOME))*0.9, IMPUTED_ALL_MOBILOME)

nosparse_ALLMOBILOME<-filter_taxa(IMPUTED_ALL_MOBILOME, function(x){(sum(x > 20) > nsamples(IMPUTED_ALL_MOBILOME)*0.1)}, prune = TRUE)
nosparse_ALLMOBILOME

#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 2052 taxa and 20 samples ]
#sample_data() Sample Data:       [ 20 samples by 93 sample variables ]
#tax_table()   Taxonomy Table:    [ 2052 taxa by 10 taxonomic ranks ]

DS_ALLMOBILOME= phyloseq_to_deseq2(nosparse_ALLMOBILOME, ~CollectionPhase+ Post_captureSampleConc_ng_ul + Shotgun_Post_Host_Rem + GEND + BMI + EVSMOK + EATPRK)
DS_ALLMGES<- estimateSizeFactors(DS_ALLMOBILOME, type="poscounts")
DS_ALLMGES$CollectionPhase <- relevel(DS_ALLMGES$CollectionPhase, ref = "WORKDAY_END")
DS_ALLMGES= DESeq(DS_ALLMGES, test= "Wald", fitType= "parametric") 
alpha= 0.01
DS_ALLMGES_RESULT= results(DS_ALLMGES, alpha = alpha)
DS_ALLMGES <- lfcShrink(DS_ALLMGES, res= DS_ALLMGES_RESULT, type="ashr")
DS_ALLMGES_RESULT= DS_ALLMGES_RESULT[order(DS_ALLMGES_RESULT$padj, na.last = NA), ]
DS_ALLMGES_RESULT
summary(DS_ALLMGES_RESULT)
#out of 2497 with nonzero total read count
#adjusted p-value < 0.01
#LFC > 0 (up)       : 182, 7.3%
#LFC < 0 (down)     : 112, 4.5%
#out of 2052 with nonzero total read count
#adjusted p-value < 0.01
#LFC > 0 (up)       : 103, 5%
#LFC < 0 (down)     : 51, 2.5%

R<-as.data.frame(DS_ALLMGES_RESULT)
R_SIG_ALLMGES<- head(R[, 1:6], 154)
R_SIG_ALLMGES<- R_SIG_ALLMGES[order(R_SIG_ALLMGES$baseMean, na.last = NA), ]
write.csv(R_SIG_ALLMGES, "C:/Users/is233/OneDrive/Desktop/AG_WORKER_METAGENOMICS_STUDY/MANUSCRIPT/ANALYSIS/WORKER_MOBILOME/DESEQ_ALLMGE_WEPS.csv", row.names = TRUE)
DS_ALLARGS <- estimateSizeFactors(DS_ALLARGS)
DS_ALLARGS <- estimateDispersions(DS_ALLARGS)
plotDispEsts(DS_ALLARGS)
par(mfrow=c(1,1))
# Make a basic volcano plot
with(R, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot", xlim=c(-50,50)))
with(subset(R, padj<.01 ), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
with(subset(R, padj<.01 & abs(log2FoldChange)>1.5), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))
EnhancedVolcano(R,
                lab = rownames(R),
                x = 'log2FoldChange',
                y = 'padj',  pCutoff = 10e-3,
                FCcutoff = 1.5,
                pointSize = log10(R$baseMean)+2,
                labSize = 5,
                shape=16,
                colAlpha = 0.4, legendLabels=c('Not sig.','Log (base 2) FC','FDRadj p <0.01',
                                               'FDRadj p <0.01 & Log (base 2) FC'),
                legendPosition = 'right',
                legendLabSize = 12,
                legendIconSize = 4.0)

##ggplot2 Volcanoplot

R$diffexpressed<- 'NO'
R$diffexpressed[R$log2FoldChange>1.5 & R$padj<0.01]<- 'Post-shower_T3'
R$diffexpressed[R$log2FoldChange<(-1.5) & R$padj<0.01]<- 'Workday end_T2'

ggplot(data=R, aes(x=log2FoldChange, y=-log10(padj), col=diffexpressed))+ geom_vline(xintercept = c(-1.5, 1.5), col='gray', linetype='dashed') + geom_hline(yintercept= c(-log10(10e-3)), col='grey', linetype= 'dashed') + geom_point(size=log10(R$baseMean), shape=16, stroke=2)+
  
  theme_classic(base_size = 20)+ theme(
    axis.title.y = element_text(face = "bold", margin = margin(0,20,0,0), size = rel(1.1)), 
    axis.title.x = element_text(hjust = 0.5, face = "bold", margin = margin(20,0,0,0), size = rel(1.1)), 
    plot.title = element_text(hjust = 0.5))

sum( rowMeans( counts(DS_ALLMGES, normalized=TRUE)) > 5 )#check number of rows for nsub
vsdata <- vst(DS_ALLMGES, blind=TRUE, nsub=2015)
DESEQ_PCA_MGE_GENES_WEPS<-plotPCA(vsdata, intgroup="CollectionPhase", ntop=2000) + theme_bw()
DESEQ_PCA_MGE_GENES_WEPS


DS_SIG_ALLMGES<- rownames(DS_ALLMGES_RESULT[1:90, ]) #Selecting the bottom 100 with the lowest adjusted p values
DS_SIG_ALLMGES
#total=median(sample_sums(IMPUTED_ALLRESISTOME_physeq))
#standf=function(x, t=total) round(t * (x / sum(x)))
DS_REL_ALLMGES<-phyloseq_transform_css(IMPUTED_ALL_MOBILOME)
#DS_REL_ALLARGS<- transform_sample_counts(IMPUTED_ALLRESISTOME_physeq, function(x) x/sum(x)*10)
DS_SIG_REL_ALLMGES<- prune_taxa(DS_SIG_ALLMGES, DS_REL_ALLMGES)

DS_SIG_REL_ALLMGES


##############SUBSET FOR WORK_START_POST_SHOWER

IMPALL_MOBILOME = otu_table(MOBILOME_physeq, taxa_are_rows = T)
MOBGENES = tax_table(mobgen_mat)
MOB_SAMPLES = sample_data(mobilomesamples_df)
IMPUTED_ALLMOBILOME_physeq <- phyloseq(IMPALL_MOBILOME, MOBGENES, MOB_SAMPLES)
IMPUTED_ALLMOBILOME_physeq

sample_data(IMPUTED_ALLMOBILOME_physeq)$CollectionPhase <- as.factor(sample_data(IMPUTED_ALLMOBILOME_physeq)$CollectionPhase) # factorize for DESeq2

sample_data(IMPUTED_ALLMOBILOME_physeq)$CollectionPhase <- as.factor(sample_data(IMPUTED_ALLMOBILOME_physeq)$CollectionPhase) # factorize for DESeq2
sample_data(IMPUTED_ALLMOBILOME_physeq)$EVSMOK <- as.factor(sample_data(IMPUTED_ALLMOBILOME_physeq)$EVSMOK) # factorize for DESeq2
sample_data(IMPUTED_ALLMOBILOME_physeq)$Shotgun_Post_Host_Rem  <- as.numeric(sample_data(IMPUTED_ALLMOBILOME_physeq)$Shotgun_Post_Host_Rem ) # factorize for DESeq2
sample_data(IMPUTED_ALLMOBILOME_physeq)$Shotgun_Post_Host_Rem  <- scale(sample_data(IMPUTED_ALLMOBILOME_physeq)$Shotgun_Post_Host_Rem , center = TRUE, scale=TRUE) # factorize for DESeq2
sample_data(IMPUTED_ALLMOBILOME_physeq)$Post_captureSampleConc_ng_ul    <- as.numeric(sample_data(IMPUTED_ALLMOBILOME_physeq)$Post_captureSampleConc_ng_ul) # factorize for DESeq2
sample_data(IMPUTED_ALLMOBILOME_physeq)$Post_captureSampleConc_ng_ul  <- scale(sample_data(IMPUTED_ALLMOBILOME_physeq)$Post_captureSampleConc_ng_ul , center = TRUE, scale=TRUE) # factorize for DESeq2
sample_data(IMPUTED_ALLMOBILOME_physeq)$GEND   <- as.factor(sample_data(IMPUTED_ALLMOBILOME_physeq)$GEND) # factorize for DESeq2
sample_data(IMPUTED_ALLMOBILOME_physeq)$AGE   <- scale(sample_data(IMPUTED_ALLMOBILOME_physeq)$AGE  , center = TRUE, scale=TRUE) # factorize for DESeq2

sample_data(IMPUTED_ALLMOBILOME_physeq)$BMI   <- scale(sample_data(IMPUTED_ALLMOBILOME_physeq)$BMI  , center = TRUE, scale=TRUE) # factorize for DESeq2
sample_data(IMPUTED_ALLMOBILOME_physeq)$EATPK   <- scale(sample_data(IMPUTED_ALLMOBILOME_physeq)$EATPRK  , center = TRUE, scale=TRUE) # factorize for DESeq2
sample_data(IMPUTED_ALLMOBILOME_physeq)$HNDPRK   <- scale(sample_data(IMPUTED_ALLMOBILOME_physeq)$HNDPRK  , center = TRUE, scale=TRUE) # factorize for DESeq2
sample_data(IMPUTED_ALLMOBILOME_physeq)$EATPRK   <- scale(sample_data(IMPUTED_ALLMOBILOME_physeq)$EATPRK  , center = TRUE, scale=TRUE) # factorize for DESeq2
sample_data(IMPUTED_ALLMOBILOME_physeq)$HNDPRK   <- scale(sample_data(IMPUTED_ALLMOBILOME_physeq)$HNDPRK  , center = TRUE, scale=TRUE) # factorize for DESeq2

filtered_IMPUTED_ALLMOBILOME_physeq<-taxa_filter(IMPUTED_ALLMOBILOME_physeq, frequency = 0.8, treatment = 'CollectionPhase')
filtered_IMPUTED_ALLMOBILOME_physeq

#Filter sparse features of IMPUTED_ALLMOBILOME_physeq object so that we only consider those features with less than 90% zeros
#nosparse_ALLMOBILOME<- prune_taxa(rowSums(otu_table(IMPUTED_ALL_MOBILOME)==0) < ncol(otu_table(IMPUTED_ALL_MOBILOME))*0.9, IMPUTED_ALL_MOBILOME)

IMPUTED_ALL_MOBILOME <- subset_samples(filtered_IMPUTED_ALLMOBILOME_physeq, CollectionPhase=="POST_SHOWER" | CollectionPhase=="WORKDAY_START")
IMPUTED_ALL_MOBILOME
#Filter sparse features of IMPUTED_ALLMOBILOME_physeq object so that we only consider those features with less than 90% zeros
#IMPUTED_ALL_MOBILOME
#nosparse_ALLMOBILOME<- prune_taxa(rowSums(otu_table(IMPUTED_ALL_MOBILOME)==0) < ncol(otu_table(IMPUTED_ALL_MOBILOME))*0.9, IMPUTED_ALL_MOBILOME)

nosparse_ALLMOBILOME<-filter_taxa(IMPUTED_ALL_MOBILOME, function(x){(sum(x > 20) > nsamples(IMPUTED_ALL_MOBILOME)*0.1)}, prune = TRUE)
nosparse_ALLMOBILOME

#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 1989 taxa and 20 samples ]
#sample_data() Sample Data:       [ 20 samples by 94 sample variables ]
#tax_table()   Taxonomy Table:    [ 1989 taxa by 10 taxonomic ranks ]

DS_ALLMOBILOME= phyloseq_to_deseq2(nosparse_ALLMOBILOME, ~CollectionPhase+ Post_captureSampleConc_ng_ul + Shotgun_Post_Host_Rem + GEND + BMI + EVSMOK + EATPRK)
DS_ALLMGES<- estimateSizeFactors(DS_ALLMOBILOME, type="poscounts")
DS_ALLMGES$CollectionPhase <- relevel(DS_ALLMGES$CollectionPhase, ref = "WORKDAY_START")
DS_ALLMGES= DESeq(DS_ALLMGES, test= "Wald", fitType= "parametric") 
alpha= 0.01
DS_ALLMGES_RESULT= results(DS_ALLMGES, alpha = alpha)
DS_ALLMGES <- lfcShrink(DS_ALLMGES, res= DS_ALLMGES_RESULT, type="ashr")
DS_ALLMGES_RESULT= DS_ALLMGES_RESULT[order(DS_ALLMGES_RESULT$padj, na.last = NA), ]
DS_ALLMGES_RESULT
summary(DS_ALLMGES_RESULT)
#out of 2400 with nonzero total read count
#adjusted p-value < 0.01
#LFC > 0 (up)       : 150, 6.2%
#LFC < 0 (down)     : 101, 4.2%
#out of 1989 with nonzero total read count
#adjusted p-value < 0.01
#LFC > 0 (up)       : 74, 3.7%
#LFC < 0 (down)     : 56, 2.8%
R<-as.data.frame(DS_ALLMGES_RESULT)
R_SIG_ALLMGES<- head(R[, 1:6], 130)
R_SIG_ALLMGES<- R_SIG_ALLMGES[order(R_SIG_ALLMGES$baseMean, na.last = NA), ]
write.csv(R_SIG_ALLMGES, "C:/Users/is233/OneDrive/Desktop/AG_WORKER_METAGENOMICS_STUDY/MANUSCRIPT/ANALYSIS/WORKER_MOBILOME/DESEQ_ALLMGE_WSPS.csv", row.names = TRUE)
DS_ALLARGS <- estimateSizeFactors(DS_ALLARGS)
DS_ALLARGS <- estimateDispersions(DS_ALLARGS)
plotDispEsts(DS_ALLARGS)
par(mfrow=c(1,1))
# Make a basic volcano plot
with(R, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot", xlim=c(-50,50)))
with(subset(R, padj<.01 ), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
with(subset(R, padj<.01 & abs(log2FoldChange)>1.5), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))
EnhancedVolcano(R,
                lab = rownames(R),
                x = 'log2FoldChange',
                y = 'padj',  pCutoff = 10e-3,
                FCcutoff = 1.5,
                pointSize = log10(R$baseMean)+2,
                labSize = 5,
                shape=16,
                colAlpha = 0.4, legendLabels=c('Not sig.','Log (base 2) FC','FDRadj p <0.01',
                                               'FDRadj p <0.01 & Log (base 2) FC'),
                legendPosition = 'right',
                legendLabSize = 12,
                legendIconSize = 4.0)

##ggplot2 Volcanoplot

R$diffexpressed<- 'NO'
R$diffexpressed[R$log2FoldChange>1.5 & R$padj<0.01]<- 'Post-shower_T3'
R$diffexpressed[R$log2FoldChange<(-1.5) & R$padj<0.01]<- 'Workday start_T1'

ggplot(data=R, aes(x=log2FoldChange, y=-log10(padj), col=diffexpressed))+ geom_vline(xintercept = c(-1.5, 1.5), col='gray', linetype='dashed') + geom_hline(yintercept= c(-log10(10e-3)), col='grey', linetype= 'dashed') + geom_point(size=log10(R$baseMean), shape=16, stroke=2)+
  
  theme_classic(base_size = 20)+ theme(
    axis.title.y = element_text(face = "bold", margin = margin(0,20,0,0), size = rel(1.1)), 
    axis.title.x = element_text(hjust = 0.5, face = "bold", margin = margin(20,0,0,0), size = rel(1.1)), 
    plot.title = element_text(hjust = 0.5))


sum( rowMeans( counts(DS_ALLMGES, normalized=TRUE)) > 5 )#check number of rows for nsub
vsdata <- vst(DS_ALLMGES, blind=TRUE, nsub=1964)
DESEQ_PCA_MGE_GENES_PSWS<-plotPCA(vsdata, intgroup="CollectionPhase", ntop=2000) + theme_bw()
DESEQ_PCA_MGE_GENES_PSWS

DS_SIG_ALLMGES<- rownames(DS_ALLMGES_RESULT[1:250, ]) #Selecting the bottom 100 with the lowest adjusted p values
DS_SIG_ALLMGES
#total=median(sample_sums(IMPUTED_ALLRESISTOME_physeq))
#standf=function(x, t=total) round(t * (x / sum(x)))
DS_REL_ALLMGES<-phyloseq_transform_css(IMPUTED_ALL_MOBILOME)
DS_REL_ALLMGES<- transform_sample_counts(IMPUTED_ALLMOBILOME_physeq, function(x) x/sum(x)*10)
DS_SIG_REL_ALLMGES<- prune_taxa(DS_SIG_ALLMGES, DS_REL_ALLMGES)

DS_SIG_REL_ALLMGES

##############SUBSET FOR WORK_END_SWINE

IMPALL_MOBILOME = otu_table(MOBILOME_physeq, taxa_are_rows = T)
MOBGENES = tax_table(mobgen_mat)
MOB_SAMPLES = sample_data(mobilomesamples_df)
IMPUTED_ALLMOBILOME_physeq <- phyloseq(IMPALL_MOBILOME, MOBGENES, MOB_SAMPLES)
IMPUTED_ALLMOBILOME_physeq

sample_data(IMPUTED_ALLMOBILOME_physeq)$CollectionPhase <- as.factor(sample_data(IMPUTED_ALLMOBILOME_physeq)$CollectionPhase) # factorize for DESeq2

sample_data(IMPUTED_ALLMOBILOME_physeq)$CollectionPhase <- as.factor(sample_data(IMPUTED_ALLMOBILOME_physeq)$CollectionPhase) # factorize for DESeq2
sample_data(IMPUTED_ALLMOBILOME_physeq)$EVSMOK <- as.factor(sample_data(IMPUTED_ALLMOBILOME_physeq)$EVSMOK) # factorize for DESeq2
sample_data(IMPUTED_ALLMOBILOME_physeq)$Shotgun_Post_Host_Rem  <- as.numeric(sample_data(IMPUTED_ALLMOBILOME_physeq)$Shotgun_Post_Host_Rem ) # factorize for DESeq2
sample_data(IMPUTED_ALLMOBILOME_physeq)$Shotgun_Post_Host_Rem  <- scale(sample_data(IMPUTED_ALLMOBILOME_physeq)$Shotgun_Post_Host_Rem , center = TRUE, scale=TRUE) # factorize for DESeq2
sample_data(IMPUTED_ALLMOBILOME_physeq)$Post_captureSampleConc_ng_ul    <- as.numeric(sample_data(IMPUTED_ALLMOBILOME_physeq)$Post_captureSampleConc_ng_ul) # factorize for DESeq2
sample_data(IMPUTED_ALLMOBILOME_physeq)$Post_captureSampleConc_ng_ul  <- scale(sample_data(IMPUTED_ALLMOBILOME_physeq)$Post_captureSampleConc_ng_ul , center = TRUE, scale=TRUE) # factorize for DESeq2
sample_data(IMPUTED_ALLMOBILOME_physeq)$GEND   <- as.factor(sample_data(IMPUTED_ALLMOBILOME_physeq)$GEND) # factorize for DESeq2
sample_data(IMPUTED_ALLMOBILOME_physeq)$AGE   <- scale(sample_data(IMPUTED_ALLMOBILOME_physeq)$AGE  , center = TRUE, scale=TRUE) # factorize for DESeq2

sample_data(IMPUTED_ALLMOBILOME_physeq)$BMI   <- scale(sample_data(IMPUTED_ALLMOBILOME_physeq)$BMI  , center = TRUE, scale=TRUE) # factorize for DESeq2

filtered_IMPUTED_ALLMOBILOME_physeq<-taxa_filter(IMPUTED_ALLMOBILOME_physeq, frequency = 0.8, treatment = 'CollectionPhase')
filtered_IMPUTED_ALLMOBILOME_physeq

IMPUTED_ALL_MOBILOME <- subset_samples(filtered_IMPUTED_ALLMOBILOME_physeq, CollectionPhase=="SWINE" | CollectionPhase=="WORKDAY_END")
IMPUTED_ALL_MOBILOME

#Filter sparse features of IMPUTED_ALLMOBILOME_physeq object so that we only consider those features with less than 90% zeros
#IMPUTED_ALL_MOBILOME
#nosparse_ALLMOBILOME<- prune_taxa(rowSums(otu_table(IMPUTED_ALL_MOBILOME)==0) < ncol(otu_table(IMPUTED_ALL_MOBILOME))*0.9, IMPUTED_ALL_MOBILOME)


nosparse_ALLMOBILOME<-filter_taxa(IMPUTED_ALL_MOBILOME, function(x){(sum(x > 20) > nsamples(IMPUTED_ALL_MOBILOME)*0.1)}, prune = TRUE)
nosparse_ALLMOBILOME
#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 2186 taxa and 20 samples ]
#sample_data() Sample Data:       [ 20 samples by 93 sample variables ]
#tax_table()   Taxonomy Table:    [ 2186 taxa by 10 taxonomic ranks ]

DS_ALLMOBILOME= phyloseq_to_deseq2(nosparse_ALLMOBILOME, ~CollectionPhase+ Post_captureSampleConc_ng_ul + Shotgun_Post_Host_Rem)
DS_ALLMGES<- estimateSizeFactors(DS_ALLMOBILOME, type="poscounts")
DS_ALLMGES$CollectionPhase <- relevel(DS_ALLMGES$CollectionPhase, ref = "SWINE")
DS_ALLMGES= DESeq(DS_ALLMGES, test= "Wald", fitType= "parametric") 
alpha= 0.01
DS_ALLMGES_RESULT= results(DS_ALLMGES, alpha = alpha)
DS_ALLMGES <- lfcShrink(DS_ALLMGES, res= DS_ALLMGES_RESULT, type="ashr")
DS_ALLMGES_RESULT= DS_ALLMGES_RESULT[order(DS_ALLMGES_RESULT$padj, na.last = NA), ]
DS_ALLMGES_RESULT
summary(DS_ALLMGES_RESULT)
#out of 2372 with nonzero total read count
#adjusted p-value < 0.01
#LFC > 0 (up)       : 89, 3.8%
#LFC < 0 (down)     : 26, 1.1%

#out of 2186 with nonzero total read count
#adjusted p-value < 0.01
#LFC > 0 (up)       : 64, 2.9%
#LFC < 0 (down)     : 7, 0.32%

R<-as.data.frame(DS_ALLMGES_RESULT)
R_SIG_ALLMGES<- head(R[, 1:6], 71)
R_SIG_ALLMGES<- R_SIG_ALLMGES[order(R_SIG_ALLMGES$baseMean, na.last = NA), ]
write.csv(R_SIG_ALLMGES, "C:/Users/is233/OneDrive/Desktop/AG_WORKER_METAGENOMICS_STUDY/MANUSCRIPT/ANALYSIS/WORKER_MOBILOME/DESEQ_ALLMGE_SWINE_WE.csv", row.names = TRUE)
DS_ALLARGS <- estimateSizeFactors(DS_ALLARGS)
DS_ALLARGS <- estimateDispersions(DS_ALLARGS)
plotDispEsts(DS_ALLARGS)
par(mfrow=c(1,1))
# Make a basic volcano plot
with(R, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot", xlim=c(-50,50)))
with(subset(R, padj<.01 ), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
with(subset(R, padj<.01 & abs(log2FoldChange)>1.5), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))
EnhancedVolcano(R,
                lab = rownames(R),
                x = 'log2FoldChange',
                y = 'padj',  pCutoff = 10e-3,
                FCcutoff = 1.5,
                pointSize = log10(R$baseMean)+2,
                labSize = 5,
                shape=16,
                colAlpha = 0.4, legendLabels=c('Not sig.','Log (base 2) FC','FDRadj p <0.01',
                                               'FDRadj p <0.01 & Log (base 2) FC'),
                legendPosition = 'right',
                legendLabSize = 12,
                legendIconSize = 4.0)

##ggplot2 Volcanoplot

R$diffexpressed<- 'NO'
R$diffexpressed[R$log2FoldChange>1.5 & R$padj<0.01]<- 'Workday end_T2'
R$diffexpressed[R$log2FoldChange<(-1.5) & R$padj<0.01]<- 'Swine'

ggplot(data=R, aes(x=log2FoldChange, y=-log10(padj), col=diffexpressed))+ geom_vline(xintercept = c(-1.5, 1.5), col='gray', linetype='dashed') + geom_hline(yintercept= c(-log10(10e-3)), col='grey', linetype= 'dashed') + geom_point(size=log10(R$baseMean), shape=16, stroke=2)+
  
  theme_classic(base_size = 20)+ theme(
    axis.title.y = element_text(face = "bold", margin = margin(0,20,0,0), size = rel(1.1)), 
    axis.title.x = element_text(hjust = 0.5, face = "bold", margin = margin(20,0,0,0), size = rel(1.1)), 
    plot.title = element_text(hjust = 0.5))

sum( rowMeans( counts(DS_ALLMGES, normalized=TRUE)) > 5 )#check number of rows for nsub
vsdata <- vst(DS_ALLMGES, blind=TRUE, nsub=2142)
DESEQ_PCA_MGE_GENES_SWINE_WE<-plotPCA(vsdata, intgroup="CollectionPhase", ntop=2000) + theme_bw()
DESEQ_PCA_MGE_GENES_SWINE_WE

#########################################
#PLASMID DIFFERENTIAL abundance
#########################################

IMPPLASMID_MOBILOME = otu_table(IMPUTEDmobtype_PLASMID.ps, taxa_are_rows = T)
MOBGENES = tax_table(mobgen_mat)
MOB_SAMPLES = sample_data(mobilomesamples_df)
IMPUTED_PLASMID_physeq <- phyloseq(IMPUTEDmobtype_PLASMID.ps, MOBGENES, MOB_SAMPLES)



sample_data(IMPUTED_PLASMID_physeq)$CollectionPhase <- as.factor(sample_data(IMPUTED_PLASMID_physeq)$CollectionPhase) # factorize for DESeq2

sample_data(IMPUTED_PLASMID_physeq)$CollectionPhase <- as.factor(sample_data(IMPUTED_PLASMID_physeq)$CollectionPhase) # factorize for DESeq2
sample_data(IMPUTED_PLASMID_physeq)$EVSMOK <- as.factor(sample_data(IMPUTED_ALLMOBILOME_physeq)$EVSMOK) # factorize for DESeq2
sample_data(IMPUTED_PLASMID_physeq)$Shotgun_Post_Host_Rem  <- as.numeric(sample_data(IMPUTED_PLASMID_physeq)$Shotgun_Post_Host_Rem ) # factorize for DESeq2
sample_data(IMPUTED_PLASMID_physeq)$Shotgun_Post_Host_Rem  <- scale(sample_data(IMPUTED_PLASMID_physeq)$Shotgun_Post_Host_Rem , center = TRUE, scale=TRUE) # factorize for DESeq2
sample_data(IMPUTED_PLASMID_physeq)$Post_captureSampleConc_ng_ul    <- as.numeric(sample_data(IMPUTED_PLASMID_physeq)$Post_captureSampleConc_ng_ul) # factorize for DESeq2
sample_data(IMPUTED_PLASMID_physeq)$Post_captureSampleConc_ng_ul  <- scale(sample_data(IMPUTED_PLASMID_physeq)$Post_captureSampleConc_ng_ul , center = TRUE, scale=TRUE) # factorize for DESeq2
sample_data(IMPUTED_PLASMID_physeq)$GEND   <- as.factor(sample_data(IMPUTED_PLASMID_physeq)$GEND) # factorize for DESeq2
sample_data(IMPUTED_PLASMID_physeq)$AGE   <- scale(sample_data(IMPUTED_PLASMID_physeq)$AGE  , center = TRUE, scale=TRUE) # factorize for DESeq2

sample_data(IMPUTED_PLASMID_physeq)$BMI   <- scale(sample_data(IMPUTED_ALLMOBILOME_physeq)$BMI  , center = TRUE, scale=TRUE) # factorize for DESeq2

##############SUBSET FOR WORKSTART_WORK_END

IMPUTED_ALL_MOBILOME <- subset_samples(IMPUTED_PLASMID_physeq, CollectionPhase=="WORKDAY_START" | CollectionPhase=="WORKDAY_END")
IMPUTED_ALL_MOBILOME
#Filter sparse features of IMPUTED_ALLMOBILOME_physeq object so that we only consider those features with less than 90% zeros
IMPUTED_PLASMID_physeq
nosparse_ALLMOBILOME<- prune_taxa(rowSums(otu_table(IMPUTED_PLASMID_physeq)==0) < ncol(otu_table(IMPUTED_ALL_MOBILOME))*0.9, IMPUTED_ALL_MOBILOME)




DS_ALLMOBILOME= phyloseq_to_deseq2(nosparse_ALLMOBILOME, ~CollectionPhase+ Post_captureSampleConc_ng_ul + Shotgun_Post_Host_Rem + GEND + BMI + EVSMOK)
DS_ALLMGES<- estimateSizeFactors(DS_ALLMOBILOME, type="poscounts")
DS_ALLMGES$CollectionPhase <- relevel(DS_ALLMGES$CollectionPhase, ref = "WORKDAY_START")
DS_ALLMGES= DESeq(DS_ALLMGES, test= "Wald", fitType= "parametric") 
alpha= 0.01
DS_ALLMGES_RESULT= results(DS_ALLMGES, alpha = alpha)
DS_ALLMGES <- lfcShrink(DS_ALLMGES, res= DS_ALLMGES_RESULT, type="ashr")
DS_ALLMGES_RESULT= DS_ALLMGES_RESULT[order(DS_ALLMGES_RESULT$padj, na.last = NA), ]
DS_ALLMGES_RESULT
summary(DS_ALLMGES_RESULT)
R<-as.data.frame(DS_ALLMGES_RESULT)
R_SIG_ALLMGES<- head(R[, 1:6], 125)
R_SIG_ALLMGES<- R_SIG_ALLMGES[order(R_SIG_ALLMGES$baseMean, na.last = NA), ]
DS_ALLARGS <- estimateSizeFactors(DS_ALLARGS)
DS_ALLARGS <- estimateDispersions(DS_ALLARGS)
plotDispEsts(DS_ALLARGS)
par(mfrow=c(1,1))
# Make a basic volcano plot
with(R, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot", xlim=c(-50,50)))
with(subset(R, padj<.05 ), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
with(subset(R, padj<.01 & abs(log2FoldChange)>2), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))
EnhancedVolcano(R,
                lab = rownames(R),
                x = 'log2FoldChange',
                y = 'padj',  pCutoff = 10e-3,
                FCcutoff = 4,
                pointSize = 4,
                labSize = 2.0,
                shape=16,
                colAlpha = 0.4, legendLabels=c('Not sig.','Log (base 2) FC','FDRadj p <0.01',
                                               'FDRadj p <0.01 & Log (base 2) FC'),
                legendPosition = 'right',
                legendLabSize = 12,
                legendIconSize = 4.0)

##############SUBSET FOR WORK_END_POST_SHOWER

IMPPLASMID_MOBILOME = otu_table(IMPUTEDmobtype_PLASMID.ps, taxa_are_rows = T)
MOBGENES = tax_table(mobgen_mat)
MOB_SAMPLES = sample_data(mobilomesamples_df)
IMPUTED_PLASMID_physeq <- phyloseq(IMPUTEDmobtype_PLASMID.ps, MOBGENES, MOB_SAMPLES)



sample_data(IMPUTED_PLASMID_physeq)$CollectionPhase <- as.factor(sample_data(IMPUTED_PLASMID_physeq)$CollectionPhase) # factorize for DESeq2

sample_data(IMPUTED_PLASMID_physeq)$CollectionPhase <- as.factor(sample_data(IMPUTED_PLASMID_physeq)$CollectionPhase) # factorize for DESeq2
sample_data(IMPUTED_PLASMID_physeq)$EVSMOK <- as.factor(sample_data(IMPUTED_ALLMOBILOME_physeq)$EVSMOK) # factorize for DESeq2
sample_data(IMPUTED_PLASMID_physeq)$Shotgun_Post_Host_Rem  <- as.numeric(sample_data(IMPUTED_PLASMID_physeq)$Shotgun_Post_Host_Rem ) # factorize for DESeq2
sample_data(IMPUTED_PLASMID_physeq)$Shotgun_Post_Host_Rem  <- scale(sample_data(IMPUTED_PLASMID_physeq)$Shotgun_Post_Host_Rem , center = TRUE, scale=TRUE) # factorize for DESeq2
sample_data(IMPUTED_PLASMID_physeq)$Post_captureSampleConc_ng_ul    <- as.numeric(sample_data(IMPUTED_PLASMID_physeq)$Post_captureSampleConc_ng_ul) # factorize for DESeq2
sample_data(IMPUTED_PLASMID_physeq)$Post_captureSampleConc_ng_ul  <- scale(sample_data(IMPUTED_PLASMID_physeq)$Post_captureSampleConc_ng_ul , center = TRUE, scale=TRUE) # factorize for DESeq2
sample_data(IMPUTED_PLASMID_physeq)$GEND   <- as.factor(sample_data(IMPUTED_PLASMID_physeq)$GEND) # factorize for DESeq2
sample_data(IMPUTED_PLASMID_physeq)$AGE   <- scale(sample_data(IMPUTED_PLASMID_physeq)$AGE  , center = TRUE, scale=TRUE) # factorize for DESeq2

sample_data(IMPUTED_PLASMID_physeq)$BMI   <- scale(sample_data(IMPUTED_ALLMOBILOME_physeq)$BMI  , center = TRUE, scale=TRUE) # factorize for DESeq2


IMPUTED_ALL_MOBILOME <- subset_samples(IMPUTED_PLASMID_physeq, CollectionPhase=="POST_SHOWER" | CollectionPhase=="WORKDAY_END")
IMPUTED_ALL_MOBILOME
#Filter sparse features of IMPUTED_ALLMOBILOME_physeq object so that we only consider those features with less than 90% zeros
IMPUTED_PLASMID_physeq
nosparse_ALLMOBILOME<- prune_taxa(rowSums(otu_table(IMPUTED_PLASMID_physeq)==0) < ncol(otu_table(IMPUTED_PLASMID_physeq))*0.9, IMPUTED_ALL_MOBILOME)




DS_ALLMOBILOME= phyloseq_to_deseq2(nosparse_ALLMOBILOME, ~CollectionPhase+ Post_captureSampleConc_ng_ul + Shotgun_Post_Host_Rem + GEND + BMI + EVSMOK)
DS_ALLMGES<- estimateSizeFactors(DS_ALLMOBILOME, type="poscounts")
DS_ALLMGES$CollectionPhase <- relevel(DS_ALLMGES$CollectionPhase, ref = "WORKDAY_END")
DS_ALLMGES= DESeq(DS_ALLMGES, test= "Wald", fitType= "parametric") 
alpha= 0.01
DS_ALLMGES_RESULT= results(DS_ALLMGES, alpha = alpha)
DS_ALLMGES <- lfcShrink(DS_ALLMGES, res= DS_ALLMGES_RESULT, type="ashr")
DS_ALLMGES_RESULT= DS_ALLMGES_RESULT[order(DS_ALLMGES_RESULT$padj, na.last = NA), ]
DS_ALLMGES_RESULT
summary(DS_ALLMGES_RESULT)
R<-as.data.frame(DS_ALLMGES_RESULT)
R_SIG_ALLMGES<- head(R[, 1:6], 166)
R_SIG_ALLMGES<- R_SIG_ALLMGES[order(R_SIG_ALLMGES$baseMean, na.last = NA), ]
DS_ALLARGS <- estimateSizeFactors(DS_ALLARGS)
DS_ALLARGS <- estimateDispersions(DS_ALLARGS)
plotDispEsts(DS_ALLARGS)
par(mfrow=c(1,1))
# Make a basic volcano plot
with(R, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot", xlim=c(-50,50)))
with(subset(R, padj<.05 ), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
with(subset(R, padj<.01 & abs(log2FoldChange)>2), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))
EnhancedVolcano(R,
                lab = rownames(R),
                x = 'log2FoldChange',
                y = 'padj',  pCutoff = 10e-3,
                FCcutoff = 4,
                pointSize = 4,
                labSize = 2.0,
                shape=16,
                colAlpha = 0.4, legendLabels=c('Not sig.','Log (base 2) FC','FDRadj p <0.01',
                                               'FDRadj p <0.01 & Log (base 2) FC'),
                legendPosition = 'right',
                legendLabSize = 12,
                legendIconSize = 4.0)


##############SUBSET FOR WORK_START_POST_SHOWER

IMPPLASMID_MOBILOME = otu_table(IMPUTEDmobtype_PLASMID.ps, taxa_are_rows = T)
MOBGENES = tax_table(mobgen_mat)
MOB_SAMPLES = sample_data(mobilomesamples_df)
IMPUTED_PLASMID_physeq <- phyloseq(IMPUTEDmobtype_PLASMID.ps, MOBGENES, MOB_SAMPLES)


sample_data(IMPUTED_PLASMID_physeq)$CollectionPhase <- as.factor(sample_data(IMPUTED_PLASMID_physeq)$CollectionPhase) # factorize for DESeq2

sample_data(IMPUTED_PLASMID_physeq)$CollectionPhase <- as.factor(sample_data(IMPUTED_PLASMID_physeq)$CollectionPhase) # factorize for DESeq2
sample_data(IMPUTED_PLASMID_physeq)$EVSMOK <- as.factor(sample_data(IMPUTED_ALLMOBILOME_physeq)$EVSMOK) # factorize for DESeq2
sample_data(IMPUTED_PLASMID_physeq)$Shotgun_Post_Host_Rem  <- as.numeric(sample_data(IMPUTED_PLASMID_physeq)$Shotgun_Post_Host_Rem ) # factorize for DESeq2
sample_data(IMPUTED_PLASMID_physeq)$Shotgun_Post_Host_Rem  <- scale(sample_data(IMPUTED_PLASMID_physeq)$Shotgun_Post_Host_Rem , center = TRUE, scale=TRUE) # factorize for DESeq2
sample_data(IMPUTED_PLASMID_physeq)$Post_captureSampleConc_ng_ul    <- as.numeric(sample_data(IMPUTED_PLASMID_physeq)$Post_captureSampleConc_ng_ul) # factorize for DESeq2
sample_data(IMPUTED_PLASMID_physeq)$Post_captureSampleConc_ng_ul  <- scale(sample_data(IMPUTED_PLASMID_physeq)$Post_captureSampleConc_ng_ul , center = TRUE, scale=TRUE) # factorize for DESeq2
sample_data(IMPUTED_PLASMID_physeq)$GEND   <- as.factor(sample_data(IMPUTED_PLASMID_physeq)$GEND) # factorize for DESeq2
sample_data(IMPUTED_PLASMID_physeq)$AGE   <- scale(sample_data(IMPUTED_PLASMID_physeq)$AGE  , center = TRUE, scale=TRUE) # factorize for DESeq2

sample_data(IMPUTED_PLASMID_physeq)$BMI   <- scale(sample_data(IMPUTED_ALLMOBILOME_physeq)$BMI  , center = TRUE, scale=TRUE) # factorize for DESeq2

IMPUTED_PLASMID_physeq <- subset_samples(IMPUTED_PLASMID_physeq, CollectionPhase=="POST_SHOWER" | CollectionPhase=="WORKDAY_START")
IMPUTED_PLASMID_physeq
#Filter sparse features of IMPUTED_ALLMOBILOME_physeq object so that we only consider those features with less than 90% zeros
IMPUTED_PLASMID_physeq
nosparse_ALLMOBILOME<- prune_taxa(rowSums(otu_table(IMPUTED_PLASMID_physeq)==0) < ncol(otu_table(IMPUTED_PLASMID_physeq))*0.9, IMPUTED_PLASMID_physeq)




DS_ALLMOBILOME= phyloseq_to_deseq2(nosparse_ALLMOBILOME, ~CollectionPhase+ Post_captureSampleConc_ng_ul + Shotgun_Post_Host_Rem + GEND + BMI + EVSMOK)
DS_ALLMGES<- estimateSizeFactors(DS_ALLMOBILOME, type="poscounts")
DS_ALLMGES$CollectionPhase <- relevel(DS_ALLMGES$CollectionPhase, ref = "WORKDAY_START")
DS_ALLMGES= DESeq(DS_ALLMGES, test= "Wald", fitType= "parametric") 
alpha= 0.01
DS_ALLMGES_RESULT= results(DS_ALLMGES, alpha = alpha)
DS_ALLMGES <- lfcShrink(DS_ALLMGES, res= DS_ALLMGES_RESULT, type="ashr")
DS_ALLMGES_RESULT= DS_ALLMGES_RESULT[order(DS_ALLMGES_RESULT$padj, na.last = NA), ]
DS_ALLMGES_RESULT
summary(DS_ALLMGES_RESULT)
R<-as.data.frame(DS_ALLMGES_RESULT)
R_SIG_ALLMGES<- head(R[, 1:6], 481)
R_SIG_ALLMGES<- R_SIG_ALLMGES[order(R_SIG_ALLMGES$baseMean, na.last = NA), ]
DS_ALLARGS <- estimateSizeFactors(DS_ALLARGS)
DS_ALLARGS <- estimateDispersions(DS_ALLARGS)
plotDispEsts(DS_ALLARGS)
par(mfrow=c(1,1))
# Make a basic volcano plot
with(R, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot", xlim=c(-50,50)))
with(subset(R, padj<.05 ), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
with(subset(R, padj<.01 & abs(log2FoldChange)>2), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))
EnhancedVolcano(R,
                lab = rownames(R),
                x = 'log2FoldChange',
                y = 'padj',  pCutoff = 10e-3,
                FCcutoff = 4,
                pointSize = 4,
                labSize = 2.0,
                shape=16,
                colAlpha = 0.4, legendLabels=c('Not sig.','Log (base 2) FC','FDRadj p <0.01',
                                               'FDRadj p <0.01 & Log (base 2) FC'),
                legendPosition = 'right',
                legendLabSize = 12,
                legendIconSize = 4.0)


##############SUBSET FOR WORK_END_SWINE

IMPPLASMID_MOBILOME = otu_table(IMPUTEDmobtype_PLASMID.ps, taxa_are_rows = T)
MOBGENES = tax_table(mobgen_mat)
MOB_SAMPLES = sample_data(mobilomesamples_df)
IMPUTED_PLASMID_physeq <- phyloseq(IMPUTEDmobtype_PLASMID.ps, MOBGENES, MOB_SAMPLES)


sample_data(IMPUTED_PLASMID_physeq)$CollectionPhase <- as.factor(sample_data(IMPUTED_PLASMID_physeq)$CollectionPhase) # factorize for DESeq2

sample_data(IMPUTED_PLASMID_physeq)$CollectionPhase <- as.factor(sample_data(IMPUTED_PLASMID_physeq)$CollectionPhase) # factorize for DESeq2
sample_data(IMPUTED_PLASMID_physeq)$EVSMOK <- as.factor(sample_data(IMPUTED_ALLMOBILOME_physeq)$EVSMOK) # factorize for DESeq2
sample_data(IMPUTED_PLASMID_physeq)$Shotgun_Post_Host_Rem  <- as.numeric(sample_data(IMPUTED_PLASMID_physeq)$Shotgun_Post_Host_Rem ) # factorize for DESeq2
sample_data(IMPUTED_PLASMID_physeq)$Shotgun_Post_Host_Rem  <- scale(sample_data(IMPUTED_PLASMID_physeq)$Shotgun_Post_Host_Rem , center = TRUE, scale=TRUE) # factorize for DESeq2
sample_data(IMPUTED_PLASMID_physeq)$Post_captureSampleConc_ng_ul    <- as.numeric(sample_data(IMPUTED_PLASMID_physeq)$Post_captureSampleConc_ng_ul) # factorize for DESeq2
sample_data(IMPUTED_PLASMID_physeq)$Post_captureSampleConc_ng_ul  <- scale(sample_data(IMPUTED_PLASMID_physeq)$Post_captureSampleConc_ng_ul , center = TRUE, scale=TRUE) # factorize for DESeq2
sample_data(IMPUTED_PLASMID_physeq)$GEND   <- as.factor(sample_data(IMPUTED_PLASMID_physeq)$GEND) # factorize for DESeq2
sample_data(IMPUTED_PLASMID_physeq)$AGE   <- scale(sample_data(IMPUTED_PLASMID_physeq)$AGE  , center = TRUE, scale=TRUE) # factorize for DESeq2

sample_data(IMPUTED_PLASMID_physeq)$BMI   <- scale(sample_data(IMPUTED_ALLMOBILOME_physeq)$BMI  , center = TRUE, scale=TRUE) # factorize for DESeq2

IMPUTED_PLASMID_physeq <- subset_samples(IMPUTED_PLASMID_physeq, CollectionPhase=="SWINE" | CollectionPhase=="WORKDAY_END")
IMPUTED_PLASMID_physeq
#Filter sparse features of IMPUTED_ALLMOBILOME_physeq object so that we only consider those features with less than 90% zeros
IMPUTED_PLASMID_physeq
nosparse_ALLMOBILOME<- prune_taxa(rowSums(otu_table(IMPUTED_PLASMID_physeq)==0) < ncol(otu_table(IMPUTED_PLASMID_physeq))*0.9, IMPUTED_PLASMID_physeq)

DS_ALLMOBILOME= phyloseq_to_deseq2(nosparse_ALLMOBILOME, ~CollectionPhase+ Post_captureSampleConc_ng_ul + Shotgun_Post_Host_Rem)
DS_ALLMGES<- estimateSizeFactors(DS_ALLMOBILOME, type="poscounts")
DS_ALLMGES$CollectionPhase <- relevel(DS_ALLMGES$CollectionPhase, ref = "SWINE")
DS_ALLMGES= DESeq(DS_ALLMGES, test= "Wald", fitType= "parametric") 
alpha= 0.01
DS_ALLMGES_RESULT= results(DS_ALLMGES, alpha = alpha)
DS_ALLMGES <- lfcShrink(DS_ALLMGES, res= DS_ALLMGES_RESULT, type="ashr")
DS_ALLMGES_RESULT= DS_ALLMGES_RESULT[order(DS_ALLMGES_RESULT$padj, na.last = NA), ]
DS_ALLMGES_RESULT
summary(DS_ALLMGES_RESULT)
R<-as.data.frame(DS_ALLMGES_RESULT)
R_SIG_ALLMGES<- head(R[, 1:6], 307)
R_SIG_ALLMGES<- R_SIG_ALLMGES[order(R_SIG_ALLMGES$baseMean, na.last = NA), ]
DS_ALLARGS <- estimateSizeFactors(DS_ALLARGS)
DS_ALLARGS <- estimateDispersions(DS_ALLARGS)
plotDispEsts(DS_ALLARGS)
par(mfrow=c(1,1))
# Make a basic volcano plot
with(R, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot", xlim=c(-50,50)))
with(subset(R, padj<.05 ), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
with(subset(R, padj<.01 & abs(log2FoldChange)>2), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))
EnhancedVolcano(R,
                lab = rownames(R),
                x = 'log2FoldChange',
                y = 'padj',  pCutoff = 10e-3,
                FCcutoff = 4,
                pointSize = 4,
                labSize = 2.0,
                shape=16,
                colAlpha = 0.4, legendLabels=c('Not sig.','Log (base 2) FC','FDRadj p <0.01',
                                               'FDRadj p <0.01 & Log (base 2) FC'),
                legendPosition = 'right',
                legendLabSize = 12,
                legendIconSize = 4.0)



################################################################################
#                           cmultRepl in zCompositions                         #      
################################################################################
#Trialing a bayesian approach for zero imputations, and using the virus phyloseq object for analysis. Principled methods for the imputation of zeros, left-censored and missing data in compositional data sets (Palarea-Albaladejo and MartinFernandez (2015) <doi:10.1016/j.chemolab.2015.02.019>).

IMPUTEDmobtype_VIRUS.ps <- subset_taxa(MOBILOME_physeq, Type=="Virus")

VirusCounts <- as.data.frame(phyloseq::otu_table(IMPUTEDmobtype_VIRUS.ps))
if(phyloseq::taxa_are_rows(IMPUTEDmobtype_VIRUS.ps) == FALSE){ otus <- t(VirusCounts) }
VirusCounts

#zCompositions component of the anlaysis, using bayesian approach:


VIRUS_cmultRepl<-cmultRepl(VirusCounts,  label = 0, method = c("GBM"), output = c("p-counts"), frac = 0.2,
          threshold = 0.5, adjust = TRUE, t= NULL, s= NULL, z.warning = 0.8, suppress.print = FALSE)

VIRUS_cmultRepl

#Now we must reconvert this OTU table to the phyloseq object

IMPUTEDmobtype_VIRUS.ps <-phyloseq::otu_table(IMPUTEDmobtype_VIRUS.ps) <- phyloseq::otu_table(VIRUS_cmultRepl, taxa_are_rows = T)

return(IMPUTEDmobtype_VIRUS.ps)
IMPUTEDmobtype_VIRUS.ps

#Create a new phyloseq object to include the imputed virus OTU table

IMPVIRUS_MOBILOME = otu_table(IMPUTEDmobtype_VIRUS.ps, taxa_are_rows = T)
MOBGENES = tax_table(mobgen_mat)
MOB_SAMPLES = sample_data(mobilomesamples_df)
IMPUTED_VIRUS_physeq <- phyloseq(IMPVIRUS_MOBILOME, MOBGENES, MOB_SAMPLES)


#Now let's retry the ordination w/ this new phyloseq object

sum(taxa_sums(mobtype_PLASMID.ps)==0) # 0 taxa not present across all samples
PLASMID_physeq.ps.pruned <- prune_taxa(taxa_sums(IMPUTED_PLASMID_physeq) > 0.000001, IMPUTED_PROPHAGE_physeq)#prunning low-abundance taxa of gene groups

VIRUS_physeq.ps.pruned.css <- phyloseq_transform_css(IMPUTED_VIRUS_physeq)

VIRUS_physeq.ps.pruned.dist <- vegdist(decostand(t(otu_table(IMPUTED_VIRUS_physeq)),"rclr"), method = "euclidean")
set.seed(1999)
#PROPHAGE_physeq.ps.pruned.ord <- vegan::rda(comm = dist((PROPHAGE_physeq.ps.pruned.dist)), distance = "euclidean", k=3, try = 50, trymax = 1000, autotransform = F)

VIRUS_physeq.ps.pruned.ord<- ordinate(VIRUS_physeq.ps.pruned.css, method = "PCoA", distance = "euclidean")

VIRUS_plot<-plot_ordination(VIRUS_physeq.ps.pruned.css,VIRUS_physeq.ps.pruned.ord, type = "samples", color="CollectionPhase")
VIRUS_plot

#PROPHAGEcentroid<- PROPHAGE_plot$data %>% group_by(CollectionPhase) %>% summarise(NMDS1=mean(NMDS1), NMDS2=mean(NMDS2), .groups = "drop")


VIRUScentroid<- VIRUS_plot$data %>% group_by(CollectionPhase) %>% summarise(Axis.1=mean(Axis.1), Axis.2=mean(Axis.2), .groups = "drop")

#stressplot(plot_ord_all)
#str(plot_ord_all)

VIRUS_ordplot <- ggplot(VIRUS_plot$data,VIRUS_plot$mapping) + stat_ellipse(geom = "polygon", type= "t", level = 0.70, lty =1, lwd = 1.3, alpha= 0.07, show.legend =FALSE) + theme_bw() +geom_point(aes(fill=CollectionPhase), size=6, alpha=0.3, stroke=1.0,) + geom_point(data=VIRUScentroid, size=13, shape=16, fill="black",mapping= aes(colour=CollectionPhase), stroke=4, alpha= 0.8, show.legend = FALSE) + scale_color_manual(name="Collection phase", breaks = c("WORKDAY_START", "WORKDAY_END", "POST_SHOWER", "SWINE", "ENVIRONMENT"), values = c("#789ED0", "#F79D40", "#72BF56", "#E4655A", "#AB8BCB"), labels= c("WORKDAY_START"= "Workday start", "WORKDAY_END"="Workday end", "POST_SHOWER"="Post-shower", "SWINE"="Swine", "ENVIRONMENT"="Environment")) +
  theme(axis.line = element_line(), 
        text=element_text(family="Times",  size=10),
        panel.background = element_blank(),
        plot.margin = unit(c(0.3,0.06,0.3,0.06), "cm"),
        strip.background = element_blank(),
        strip.text = element_blank(),
        legend.key.size = unit(0.02, 'cm'),
        legend.position = c(0.8,0.8),
        legend.background = element_rect(fill = NA),
        legend.key.height = unit(.08, 'cm'), 
        legend.key.width = unit(.08, 'cm'),
        panel.border = element_rect(colour = "black", size = 1.0),
        axis.text = element_text(size = 28, colour = "black"),
        axis.title.y = element_text(size = 18, vjust = 1.75),
        axis.title.x = element_text(size = 18, vjust = -1.5)) + coord_fixed(xlim=c(-100,100), ylim = c(-50,50))  
VIRUS_ordplot

virus.anosim<- anosim(VIRUS_physeq.ps.pruned.dist, pmeta$CollectionPhase, perm = 10000)
virus.anosim
#ANOSIM statistic R: 0.08793 
#Significance: 0.024398 

virus.adonis <- pairwise.adonis(VIRUS_physeq.ps.pruned.dist, pmeta$CollectionPhase, perm = 10000, p.adjust.m = "BH")
virus.adonis

#pairs Df SumsOfSqs   F.Model         R2    p.value p.adjusted sig
#1  WORKDAY_START vs POST_SHOWER  1  927.6498 0.4947159 0.02674904 0.98230177 0.98230177    
#2  WORKDAY_START vs ENVIRONMENT  1 3192.3893 2.0267821 0.16852240 0.04619538 0.10958904    
#3        WORKDAY_START vs SWINE  1 9699.7543 3.9827270 0.18117529 0.00009999 0.00099990  **
#4  WORKDAY_START vs WORKDAY_END  1 3594.0975 2.2940468 0.11304038 0.02719728 0.09065760    
#5    POST_SHOWER vs ENVIRONMENT  1 2957.1154 1.5774124 0.13624913 0.13888611 0.19840873    
#6          POST_SHOWER vs SWINE  1 8989.5088 3.4550083 0.16103505 0.00049995 0.00249975   *
#7    POST_SHOWER vs WORKDAY_END  1 3271.2539 1.8874851 0.09490818 0.05479452 0.10958904    
#8          ENVIRONMENT vs SWINE  1  993.9725 0.3447376 0.03332492 0.83471653 0.92746281    
#9    ENVIRONMENT vs WORKDAY_END  1 1034.4802 0.7839798 0.07269856 0.74062594 0.92578242    
#10         SWINE vs WORKDAY_END  1 3223.1359 1.4053550 0.07242099 0.07629237 0.12715395   


