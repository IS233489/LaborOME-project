#Retrieve deduplicated AMR matrix (.CSV/.TXT) from AMR++ pipeline
#Censor out row-level ARG accessions requires SNP confirmation ('RequiresSNPConfirmation') using a linux machine: sed '/RequiresSNPConfirmation/d' AMR_analytic_matrix.csv > noSNPConfirmation_resistome_notransposition.csv (this one has updated IDs w/o the sample number added from AMR++)


library("installr")
install.packages("rlang")
library("rlang")
remove.packages("rlang")
library("openxlsx")
BiocManager::install("phyloseq")
library(phyloseq)
detach(package:phyloseq, unload=TRUE)
install.packages("vctrs")
remove.packages("vctrs")
install.packages("ggplot2")
library("ggplot2")
install.packages("dplyr")
remove.packages("dplyr")
library("dplyr")
detach(package:dplyr, unload=TRUE)
update.packages("dplyr")
library("ggthemes")
install.packages("lme4")
library("lme4")
install.packages("glmm")
library(glmm)
library("lmerTest")
library("lsmeans")
library("ggrepel")
install.packages("purrr")
library(purrr)
library("pbkrtest")
uninstall.packages("package_name")
install.packages("ggpubr")
library("ggpubr")
detach(package:ggpubr, unload=TRUE)
install.packages("patchwork")
library(patchwork)
library(cowplot)
library(metagenomeSeq)
library(scales)
library(vegan)
library(pairwiseAdonis)
devtools::install_github("vmikk/metagMisc")
library(metagMisc)
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("ALDEx2")
library(ALDEx2)
install.packages("randomcoloR")
library(randomcoloR)
install.packages("htmltools") 
library(htmltools)
BiocManager::install("microbiome")
library(microbiome)
devtools::install_github("microsud/microbiomeutilities", force=TRUE)
library(microbiome)
install.packages("BiocManager")
if (!requireNamespace("remotes", quietly=TRUE))
  install.packages("remotes")
remotes::install_github("YuLab-SMU/MicrobiotaProcess")
install.packages("DESeq2")
library(DESeq2)
pd = position_jitterdodge(dodge.width = 0.75, jitter.width = 0.3)
BiocManager::install("apeglm")
library(apeglm)
BiocManager::install("ashr", force = TRUE)
library(ashr)
remove.packages("remotes")
install.packages("remotes")
remotes::install_github("gauravsk/ranacapa")
library(ranacapa)
remotes::install_github("YuLab-SMU/MicrobiotaProcess")
library(MicrobiotaProcess)
install.packages(c("devtools", "RcppEigen", "RcppParallel", "Rtsne", "ggforce", "units"))
library(units)
remotes::install_github('schuyler-smith/phylosmith')
library(phylosmith)
library(matrixStats)
BiocManager::install(sparseMatrixStats)

install.packages()

install.packages(
  "microViz",
  repos = c(davidbarnett = "https://david-barnett.r-universe.dev", getOption("repos"))
)
library(microViz)


#setwd(""") set working directory here



#Metadata file name:: FINAL_SHOTGUN_MCOHS_WORKER_METADATA

resistome_metadata <- read.csv("./FINAL_SHOTGUN_MCOHS_WORKER_METADATA.csv")   ##Metadata from Microbiome files with added shotgun seq data
head(resistome_metadata)
tail(resistome_metadata)

resistome_metadata$CollectionPhase <- factor (resistome_metadata$CollectionPhase, levels= c("WORKDAY_START",
                                                                                    "WORKDAY_END",
                                                                                    "POST_SHOWER",
                                                                                    'SWINE',
                                                                                    "ENVIRONMENT",
                                                                                    "MockComm",
                                                                                    "NegCtrl"))

resistome_metadata$OccupationalTask <- factor (resistome_metadata$OccupationalTask, levels= c("FARROWING",
                                                                                      "GESTATION",
                                                                                      "SOW_GESTATION",
                                                                                      'SWINE_PEN',
                                                                                      "ENV_FARROWING",
                                                                                      "ENV_GESTATION",
                                                                                      "MOVING_FOSTERING",
                                                                                      "VISITING VETERINARIAN",
                                                                                      "MockComm",
                                                                                      "NegCtrl"))
resistome_metadata$TASKEXPOSURE <- factor (resistome_metadata$TASKEXPOSURE, levels= c("DIRECT", "INDIRECT","SWINE", "ENVIRONMENT", "MockComm", "NegCtrl"))

resistome_metadata$qpcr_16s_copies_ul<- as.numeric(as.character(resistome_metadata$qpcr_16s_copies_ul))

resistome_metadata$DIRCOHRRATE<- as.numeric(as.character(resistome_metadata$DIRCOHRRATE))

resistome_metadata$DIRCOWKRATE <- as.numeric(as.character(resistome_metadata$DIRCOWKRATE))

resistome_metadata$Shotgun_PE_Total <- as.numeric(as.character(resistome_metadata$Shotgun_PE_Total))

resistome_metadata$Shotgun_PE_Post_Trim <- as.numeric(as.character(resistome_metadata$Shotgun_PE_Post_Trim))

resistome_metadata$Shotgun_Post_Host_Rem <- as.numeric(as.character(resistome_metadata$Shotgun_Post_Host_Rem))

resistome_metadata<- resistome_metadata %>%
  mutate(Host_ReadCount= Shotgun_PE_Post_Trim - Shotgun_Post_Host_Rem)

resistome_metadata$Host_ReadCount <- as.numeric(as.character(resistome_metadata$Host_ReadCount))

###### ANALYSIS OF SHOTGUN RESISTOME READ DEPTH (TRIMMED/HOST-REMOVED) AS A FUNCTION OF COLLECTION PHASE, OCCUPATIONAL TASK, AND TASK EXPOSURE (DIRCT VS. INDIRECT)  ######


#Raw untrimmed read depth analysis

Untrimmed_raw_reads <-
  resistome_metadata %>%
  filter(!CollectionPhase %in% c("MockComm","NegCtrl")) %>%
  ggplot(aes(CollectionPhase, y= Shotgun_PE_Total, fill=CollectionPhase))+
  geom_boxplot(alpha=0.3, outlier.colour = 'white')+ylim(0,NA)+
  geom_point(aes(color = CollectionPhase), size=2, alpha=0.6, position = pd)+
  scale_fill_tableau(palette = "Classic 10 Medium", type = "regular",direction = 1)+  #Tableau 10
  scale_color_tableau(palette = "Classic 10 Medium", type = "regular",direction = 1)+  #Tableau 10
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  theme(legend.position = "none")
plot(Untrimmed_raw_reads)
topptx(filename = "rawreadbyCollectionPhase_plot.pptx", width=5, height=4) 

model_censored_Untrimmed_raw_reads <-subset (resistome_metadata, !CollectionPhase %in% c("MockComm", "NegCtrl", "ENVIRONMENT"))
model_Untrimmed_raw_reads<-
  lmer((Shotgun_PE_Total) ~ CollectionPhase +(1|Worker), data = model_censored_Untrimmed_raw_reads, REML=F)
summary(model_Untrimmed_raw_reads)

#Fixed effects:
#  Estimate Std. Error         df t value Pr(>|t|)    
#(Intercept)                   6.470e+07  5.448e+06  5.979e+43  11.876   <2e-16 ***
#  CollectionPhaseSWINE          1.878e+07  7.704e+06  5.979e+43   2.438   0.0148 *  
#  CollectionPhaseWORKDAY_END    4.353e+06  7.704e+06  5.979e+43   0.565   0.5721    
#  CollectionPhaseWORKDAY_START -2.156e+06  7.704e+06  5.979e+43  -0.280   0.7796 

lsmeans(model_Untrimmed_raw_reads, pairwise~CollectionPhase, adjust="tukey")
VarCorr(model_Untrimmed_raw_reads)
#$contrasts
#contrast                     estimate      SE       df t.ratio p.value
#POST_SHOWER - SWINE         -18783227 7704110 5.98e+43  -2.438  0.0701
#POST_SHOWER - WORKDAY_END    -4352593 7704110 5.98e+43  -0.565  0.9424
#POST_SHOWER - WORKDAY_START   2155977 7704110 5.98e+43   0.280  0.9924
#SWINE - WORKDAY_END          14430634 7704110 5.98e+43   1.873  0.2396
#SWINE - WORKDAY_START        20939204 7704110 5.98e+43   2.718  0.0333
#WORKDAY_END - WORKDAY_START   6508570 7704110 5.98e+43   0.845  0.8330

#Raw read medians / IQR / SUMS, ETC
resistome_metadata %>%
  group_by(CollectionPhase) %>%
  summarize(across(starts_with('Perc_HostCount'),IQR, na.rm=TRUE)) #the function like median, sum, IQR can be supplanted to obtain calculations

censored_raw_reads1_nohost <-
  resistome_metadata %>%
  filter(!CollectionPhase %in% c("MockComm","NegCtrl")) %>%
  ggplot(aes(CollectionPhase, y= Shotgun_Post_Host_Rem, fill=CollectionPhase))+
  geom_boxplot(alpha=0.3, outlier.colour = 'white')+ylim(0,1.2e+08)+
  geom_point(aes(color = CollectionPhase), size=2, alpha=0.6, position = pd)+
  scale_fill_tableau(palette = "Classic 10 Medium", type = "regular",direction = 1)+  #Tableau 10
  scale_color_tableau(palette = "Classic 10 Medium", type = "regular",direction = 1)+  #Tableau 10
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  theme(legend.position = "none")
plot(censored_raw_reads1_nohost)
#topptx(filename = "rawreadbyCollectionPhase_plot.pptx", width=5, height=4) 

model_censored_raw_reads1_nohost <-subset (resistome_metadata, !CollectionPhase %in% c("MockComm", "NegCtrl", "ENVIRONMENT"))
model_raw_reads1<-
  lmer((Shotgun_Post_Host_Rem) ~ CollectionPhase +(1|Worker), data = model_censored_raw_reads1_nohost, REML=F)
summary(model_raw_reads1)
#Fixed effects:
# Estimate Std. Error        df t value Pr(>|t|)    
#(Intercept)                1.749e+07  2.767e+06 8.953e+42   6.322 2.59e-10 ***
#CollectionPhaseWORKDAY_END 3.459e+06  3.913e+06 8.953e+42   0.884  0.37672    
#CollectionPhasePOST_SHOWER 1.156e+07  3.913e+06 8.953e+42   2.955  0.00313 ** 
#CollectionPhaseSWINE       2.344e+07  3.913e+06 8.953e+42   5.989 2.11e-09 ***
  
  #Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

anova(model_raw_reads1, type='III',test="F")

lsmeans(model_raw_reads1, pairwise~CollectionPhase, adjust="Tukey")
VarCorr(model_raw_reads1)

#$contrasts
#contrast                     estimate      SE   df t.ratio p.value
#WORKDAY_START - WORKDAY_END  -3459321 4125114 26.7  -0.839  0.8356
#WORKDAY_START - POST_SHOWER -11563859 4125114 26.7  -2.803  0.0434
#WORKDAY_START - SWIN28.13.75  -23438848 4125114 44.4  -5.682  <.0001
#WORKDAY_END - POST_SHOWER    -8104538 4125114 26.7  -1.965  0.2264
#WORKDAY_END - SWINE         -19979527 4125114 44.4  -4.843  0.0001
#POST_SHOWER - SWINE         -11874989 4125114 44.4  -2.879  0.0300

#Displaying above plot for host and non-host reads

censored_raw_reads1_host <-
  resistome_metadata %>%
  filter(!CollectionPhase %in% c("MockComm","NegCtrl", "ENVIRONMENT")) %>%
  ggplot(aes(CollectionPhase, y= Shotgun_PE_Post_Trim, fill=CollectionPhase))+
  geom_boxplot(alpha=0.3, outlier.colour = 'white')+ylim(0,1.2e+08)+
  geom_point(aes(color = CollectionPhase), size=2, alpha=0.6, position = pd)+
  scale_fill_tableau(palette = "Classic 10 Medium", type = "regular",direction = 1)+  #Tableau 10
  scale_color_tableau(palette = "Classic 10 Medium", type = "regular",direction = 1)+  #Tableau 10
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  theme(legend.position = "none")
plot(censored_raw_reads1_host)

patchedPlot<- censored_raw_reads1_host + censored_raw_reads1_nohost +plot_annotation(tag_levels = 'A')
patchedPlot

#Assessing depth of host reads across libraries
Host_reads_only <-
  resistome_metadata %>%
 # filter(!CollectionPhase %in% c("MockComm","NegCtrl")) %>%
  ggplot(aes(CollectionPhase, y= Host_ReadCount, fill=CollectionPhase))+
  geom_boxplot(alpha=0.3, outlier.colour = 'white')+ylim(0,NA)+
  geom_point(aes(color = CollectionPhase), size=2, alpha=0.6, position = pd)+
  scale_fill_tableau(palette = "Classic 10 Medium", type = "regular",direction = 1)+  #Tableau 10
  scale_color_tableau(palette = "Classic 10 Medium", type = "regular",direction = 1)+  #Tableau 10
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  theme(legend.position = "none")
plot(Untrimmed_raw_reads)

model_censored_raw_reads1_host <-subset (resistome_metadata, !CollectionPhase %in% c("MockComm", "NegCtrl", "ENVIRONMENT"))
model_raw_reads1_host<-
  lmer((Host_ReadCount) ~ CollectionPhase +(1|Worker), data = model_censored_raw_reads1_host, REML=F)
summary(model_raw_reads1_host)

anova(model_raw_reads1_host, type='III',test="F")

lsmeans(model_raw_reads1_host, pairwise~CollectionPhase, adjust="Tukey")
VarCorr(model_raw_reads1_host) #Significant difference in host read content and CollectionPhase accounting for raw trimmed/filtered sequencing depth

#$contrasts
#contrast                     estimate      SE       df t.ratio p.value
#POST_SHOWER - SWINE           1967471 2985345 3.68e+01   0.659  0.9117
#POST_SHOWER - WORKDAY_END   -11731478 2325891 6.67e+04  -5.044  <.0001
#POST_SHOWER - WORKDAY_START -10159200 2283196 2.75e+06  -4.450  0.0001
#SWINE - WORKDAY_END         -13698949 2806957 3.28e+01  -4.880  0.0002
#SWINE - WORKDAY_START       -12126671 3071945 3.87e+01  -3.948  0.0018
#WORKDAY_END - WORKDAY_START   1572278 2371457 1.95e+04   0.663  0.9110

#$contrasts
#contrast                     estimate      SE       df t.ratio p.value
#POST_SHOWER - SWINE         -11052893 4867704 6.40e+01  -2.271  0.1158
#POST_SHOWER - WORKDAY_END   -16497966 4478321 6.95e+22  -3.684  0.0013
#POST_SHOWER - WORKDAY_START  -8313211 4478321 6.67e+21  -1.856  0.2471
#SWINE - WORKDAY_END          -5445073 4867704 6.40e+01  -1.119  0.6795
#SWINE - WORKDAY_START         2739682 4867704 6.40e+01   0.563  0.9427
#WORKDAY_END - WORKDAY_START   8184755 4478321 4.06e+28   1.828  0.2602#Significant difference in host read content only for WORKDAY_END samples and CollectionPhase when not accounting for raw trimmed/filtered sequencing depth


raw_reads2 <- 
  resistome_metadata %>%
  filter(!OccupationalTask %in% c("MockComm", "NegCtrl", "ENV_FARROWING", "ENV_GESTATION")) %>%
  ggplot(aes(OccupationalTask, y=Shotgun_Post_Host_Rem, fill=OccupationalTask))+
  geom_boxplot(alpha=0.45)+ylim(0,NA)+
  geom_point(aes(color = OccupationalTask),size=2, alpha=0.6, position = pd)+
  scale_fill_tableau(palette = "Superfishel Stone", type = "regular",direction = 1)+  
  scale_color_tableau(palette = "Superfishel Stone", type = "regular", direction = 1)+
  theme_bw()+
  theme(legend.position = "none")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  facet_wrap(~Exposure,scales='free_x', nrow=1)
#topptx(filename = "rawreadbyCollectionPhase_plot2.pptx", width=5, height=4) 

plot(raw_reads2)

model_raw_reads2 <-subset (resistome_metadata, !OccupationalTask %in% c("MockComm", "NegCtrl", "ENV_FARROWING", "ENV_GESTATION"))
model_rawreads2<-
  lmer(Shotgun_Post_Host_Rem ~ OccupationalTask +(1|Worker), data = model_raw_reads2, REML=F)
summary(model_rawreads2)
anova(model_rawreads2, type='III',test="F")
#Fixed effects:
# Estimate Std. Error         df t value Pr(>|t|)    
#(Intercept)                            2.355e+07  2.156e+06  6.309e+43  10.923  < 2e-16 ***
#OccupationalTaskGESTATION             -1.354e+06  4.313e+06  6.309e+43  -0.314 0.753613    
#OccupationalTaskSOW_GESTATION          2.921e+07  6.819e+06  6.309e+43   4.284 1.84e-05 *** 
#OccupationalTaskSWINE_PEN              1.442e+07  3.887e+06  6.309e+43   3.710 0.000207 ***
#OccupationalTaskMOVING_FOSTERING      -4.168e+06  5.705e+06  6.309e+43  -0.731 0.465027    
#OccupationalTaskVISITING VETERINARIAN -3.648e+06  5.705e+06  6.309e+43  -0.639 0.522530    

#Type III Analysis of Variance Table with Satterthwaite's method
#                     Sum Sq    Mean Sq NumDF DenDF F value   Pr(>F)    
#OccupationalTask 2.9675e+15 5.9351e+14     5   Inf  7.0913 1.22e-06 ***

lsmeans(model_rawreads2, pairwise~OccupationalTask, adjust="tukey")
VarCorr(model_raw_reads1)
#No statistically significant differences in Shotgun depth by Occupational Task, however, there were differences in 


########################################################
##                                                     #
###                                                    #
####Resistome  analysis                                #
###                                                    # 
##                                                     #
########################################################


## phyloseq object:
res_mat<-read.csv(file.choose())      ## 
gen_mat<- read.csv(file.choose())
resistomesamples_df <- read.csv(file.choose())



##need to have row.names
row.names(res_mat) <- res_mat$Genes
res_mat <- res_mat %>% dplyr::select(-Genes) #remove the column Genes since it is now used as a row name

row.names(gen_mat) <- gen_mat$Genes
gen_mat <- gen_mat %>% dplyr::select(-Genes) 

row.names(resistomesamples_df) <- resistomesamples_df$NovSeq_SS_Lib_ID
resistomesamples_df <- resistomesamples_df %>% dplyr::select(-NovSeq_SS_Lib_ID)
#Transform into matrixes otu and tax tables (sample table can be left as data frame)

res_mat <- as.matrix(res_mat)
res_mat
gen_mat <- as.matrix(gen_mat)

##phyloseq:
RESISTOME = otu_table(res_mat, taxa_are_rows = T)
RESISTOME
GENES = tax_table(gen_mat)
RES_SAMPLES = sample_data(resistomesamples_df)
RESISTOME_physeq <- phyloseq(RESISTOME, GENES, RES_SAMPLES)
saveRDS(RESISTOME_physeq, "initial.resistome.phyloseq")
RESISTOME_physeq <- readRDS("C:/Users/is233/OneDrive/Desktop/AG_WORKER_METAGENOMICS_STUDY/MANUSCRIPT/ANALYSIS/WORKER_RESISTOME/prefiltered/initial.resistome.phyloseq")
sample_data(RESISTOME_physeq)
type.ps <-  tax_glom(RESISTOME_physeq, "Type")
class.ps <-   tax_glom(RESISTOME_physeq, "Class")
mechanism.ps <-  tax_glom(RESISTOME_physeq, "Mechanism")
group.ps <-   tax_glom(RESISTOME_physeq, "Group")
medimp_group.ps <- subset_taxa(group.ps, Group=="CTX" |
                              Group=="GES" |
                              Group=="IMI" |
                              Group=="KPC" |
                              Group=="SHV" |
                              Group=="TEM" |
                              Group=="IMP" |
                              Group=="NDM" |
                              Group=="CMY" |
                              Group=="OXA" |
                              Group=="MEC" |
                              Group=="MCR" |
                              Group=="VAT" |
                              Group=="VGB" |
                              Group=="CFR" |
                              Group=="VGA" |
                              Group=="SME" |
                              Group=="AAC6-PRIME" |
                              Group=="BLAZ" |
                              Group=="NDM" |
                              Group=="VIM" |
                              Group=="MCR" |
                              Group=="ERMB" |
                              Group=="QNRA" |
                              Group=="QNRB" |
                              Group=="QNRA" |
                              Group=="TETM" |
                              Group=="DFRA" |
                              Group=="VANYB" |
                              Group=="VANYD" |
                              Group=="VANYA" |
                              Group=="SULI")
write.csv(otu_table(medimp_group.ps), "MedImpAMR_count.csv")
write.csv(tax_table(medimp_group.ps), "MedImpAMR_taxa.csv")
#SWK9 has zero counts across the entire sample for medically important ARGs

### Relative abudance-- Type

type_relative <- transform_sample_counts(type.ps, function(x) x / sum(x) )
type_relative_long <- psmelt(type_relative)
type_relative_long <- type_relative_long %>%
  group_by(Type) %>%
  mutate(mean_relative_abund = mean(Abundance))
medimp_group.ps
type_relative_long$Type <- as.character(type_relative_long$Type)
type_relative_long$mean_relative_abund <- as.numeric(type_relative_long$mean_relative_abund)
type_relative_long$Type[type_relative_long$mean_relative_abund < 0.001] <- "Type (< 0.5%)"  ## mean_relative_abund < 0.005


type_relative_long$CollectionPhase <- factor (type_relative_long$CollectionPhase, levels= c("WORKDAY_START",
                                                                                          "WORKDAY_END",
                                                                                          "POST_SHOWER",
                                                                                          'SWINE',
                                                                                          "ENVIRONMENT"))


type_relative_long$OccupationalTask <- factor (type_relative_long$OccupationalTask, levels= c("FARROWING",
                                                                                            "GESTATION",
                                                                                            "SOW_GESTATION",
                                                                                            'SWINE_PEN',
                                                                                            "ENV_FARROWING",
                                                                                            "ENV_GESTATION",
                                                                                            "MOVING_FOSTERING",
                                                                                            "VISITING VETERINARIAN"))

type_relative_long %>%
  #filter(!CollectionPhase %in% c("NegCtrl"))%>%
  ggplot(aes(x = Sample, y = Abundance*100, fill = Type)) +
  geom_bar(stat = "identity") +
  facet_wrap(~CollectionPhase, scales = "free_x", nrow = 1) +
  geom_bar(stat = "identity", position=0, alpha=0.9) +
  theme(plot.margin = unit(c(0.1,0.5,0.5,0.5), "cm"),
        axis.title.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y = element_text(size = 14),
        legend.title = element_text(size = 14),
        legend.text = element_text( size = 14),
        plot.title = element_text(size=12, face = "bold", hjust = 0.5),
        panel.background = element_rect(fill = "white"),
        strip.background = element_rect(fill="white"),
        strip.text = element_text(face = "bold"),
        text=element_text(family="Times New Roman",  size=14)) +
  ylab("Relative abundance (%)")+
  xlab('Sample') 
#scale_fill_div
scale_fill_tableau(palette = "Color Blind", type = "regular",direction = 1) 
topptx(filename = "phy_relative_plot1.pptx", width=8, height=4)


## class:

class_relative <- transform_sample_counts(class.ps, function(x) x / sum(x) )
class_relative_long <- psmelt(class_relative)
class_relative_long <- class_relative_long %>%
  group_by(Class) %>%
  mutate(mean_relative_abund = mean(Abundance))

class_relative_long$Class <- as.character(class_relative_long$Class)
class_relative_long$mean_relative_abund <- as.numeric(class_relative_long$mean_relative_abund)
class_relative_long$Class[class_relative_long$mean_relative_abund < 0.005] <- "Low abundance classes (< 1%)"  ## mean_relative_abund < 0.01

class_relative_long$CollectionPhase <- factor (class_relative_long$CollectionPhase, levels= c("WORKDAY_START",
                                                                                          "WORKDAY_END",
                                                                                          "POST_SHOWER",
                                                                                          'SWINE',
                                                                                          "ENVIRONMENT"))

class_relative_long %>%
  ggplot(aes(x = Sample, y = Abundance*100, fill = Class)) +
  geom_bar(stat = "identity") +
  facet_wrap(~CollectionPhase, scales = "free_x", nrow = 1) +
  geom_bar(stat = "identity", alpha=0.9) +
  theme(axis.title.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y = element_text(size = 14),
        legend.title = element_text(size = 14),
        legend.text = element_text( size = 14),
        plot.title = element_text(size=12, face = "bold", hjust = 0.5),
        panel.background = element_rect(fill = "white"),
        strip.background = element_rect(fill="white"),
        strip.text = element_text(face = "bold"),
        text=element_text(family="Times New Roman",  size=14)) +
  ylab("Relative abundance (%)")+
  xlab('Sample') 
#scale_fill_div
#scale_fill_tableau(palette = "Color Blind", type = "regular",direction = 1) 
topptx(filename = "family_relative_plot1.pptx", width=9, height=4)

## mechanism:

mech_relative <- transform_sample_counts(mechanism.ps, function(x) x / sum(x) )
mech_relative_long <- psmelt(mech_relative)
mech_relative_long <- mech_relative_long %>%
  group_by(Mechanism) %>%
  mutate(mean_relative_abund = mean(Abundance))

mech_relative_long$Mechanism <- as.character(mech_relative_long$Mechanism)
mech_relative_long$mean_relative_abund <- as.numeric(mech_relative_long$mean_relative_abund)
mech_relative_long$Mechanism[mech_relative_long$mean_relative_abund < 0.005] <- "Low abundant mechanisms (< 0.5%)"  ## mean_relative_abund < 0.005

mech_relative_long$CollectionPhase <- factor (mech_relative_long$CollectionPhase, levels= c("WORKDAY_START",
                                                                                          "WORKDAY_END",
                                                                                          "POST_SHOWER",
                                                                                          'SWINE',
                                                                                          "ENVIRONMENT"))

mech_relative_long %>%
  ggplot(aes(x = Sample, y = Abundance*100, fill = Mechanism)) +
  geom_bar(stat = "identity") +
  facet_wrap(~CollectionPhase, scales = "free_x", nrow = 1) +
  geom_bar(stat = "identity", alpha=0.9) +
  theme(axis.title.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y = element_text(size = 14),
        legend.title = element_text(size = 14),
        legend.text = element_text( size = 14),
        plot.title = element_text(size=12, face = "bold", hjust = 0.5),
        panel.background = element_rect(fill = "white"),
        strip.background = element_rect(fill="white"),
        strip.text = element_text(face = "bold"),
        text=element_text(family="Times New Roman",  size=14)) +
  ylab("Relative abundance (%)")+
  xlab('Sample') 
#scale_fill_div
#scale_fill_tableau(palette = "Color Blind", type = "regular",direction = 1) 
topptx(filename = "genus_relative_plot1.pptx", width=9, height=4)


## MEDICALLY IMPORTANT ARGS: RELATIVE ABUNDANCE
medimp_group_subset.ps = subset_samples(medimp_group.ps, sample_names(medimp_group.ps) != "R1_19") #mustremove this column as all counts are empty for SWK9

medimp_group_glomtogroup <- tax_glom(medimp_group_subset.ps, taxrank = "Group") # 19 groups
medimp_group_glomtogroup

medimp_group_relab<- transform_sample_counts(medimp_group_glomtogroup, function(x) {x*100/sum(x)})


medimp_group_relab 

medimpra_group_melt <- psmelt(medimp_group_relab)#%>%
medimpra_group_melt
#  mutate(mean_relative_abund = mean(Abundance))
#medimpra_group_melt
#medimpra_group_palette <- distinctColorPalette(19)

medimpra_group_melt$Group <- as.character(medimpra_group_melt$Group)
medimpra_group_melt_medians <- medimpra_group_melt %>%
  filter(!CollectionPhase %in% c("MockComm","NegCtrl"))%>%
  group_by(CollectionPhase, Group) %>%
  mutate(median=median(Abundance), IQR= IQR(Abundance))
medimpra_group_melt_medians

# select group median > 1
#keep <- unique(medimpra_group_melt$Group[medimpra_group_melt_medians$median > 1])
#medimpra_group_melt_medians$Group[!(medimpra_group_melt_medians$Group %in% keep)] <- "< 1%"

#to get the same rows together
medimpra_group_melt_sum <- medimpra_group_melt_medians %>%
  group_by(Sample, CollectionPhase, Group)%>%
  summarise(Abundance=sum(Abundance))
  #summarise(Abundance=sum(Abundance))
medimpra_group_melt_sum


ggplot(medimpra_group_melt_sum, aes(x = Sample, y = Abundance)) + 
  geom_bar(stat = "identity", aes(fill=Group)) + 
  labs(x="", y="%") +
  facet_wrap(~CollectionPhase, scales= "free_x", nrow=1) +
  theme_classic()  +
  theme(strip.background = element_blank(), axis.text.x.bottom = element_text(angle = -90))

medimpra_group_melt_allsamps <- medimpra_group_melt_medians %>%
  group_by(CollectionPhase, Group)%>%
  summarise(Sample, Abundance, median, IQR)
#summarise(Abundance=sum(Abundance))
medimpra_group_melt_allsamps

#Bubbleplot of medimparg median abundances

xx = ggplot(medimpra_group_melt_allsamps, aes(x = CollectionPhase, y = Group)) + 
  geom_point(aes(size = IQR, fill = median), alpha = 0.75, shape = 21) + 
  scale_size_continuous(limits = c(0.000001, 100), range = c(1,17), breaks = c(1,10,50,75)) + 
  labs( x= "", y = "", size = "IQR relative abundance (%)", fill = "median relative abundance (%)")  + 
  theme(legend.key=element_blank(), axis.text.x = element_text(colour = "black", size = 12, face = "bold", angle = 90, vjust = 0.3, hjust = 1), 
        axis.text.y = element_text(colour = "black", face = "bold", size = 11), 
        legend.text = element_text(size = 10, face ="bold", colour ="black"), 
        legend.title = element_text(size = 12, face = "bold"), 
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.2), 
        legend.position = "right") #+  
  #scale_fill_manual(values = colours, guide = FALSE) #+ 
  #scale_y_continuous(limits = rev(levels(medimpra_group_melt_allsamps$Group))) 

xx

# group dataframe by CollectionPhase, calculate median rel. abundance
medimp_medians <- mutate(medimp_group_relab_glomtogroup, ~CollectionPhase, function(x) c(median=median(x$Abundance)))

dfmedimp = data.frame(Group = phyloseq::tax_table(medimp_group_relab_glomtogroup)[,"Group"], Median = rowMedians(otu_table(medimp_group_relab_glomtogroup)), row.names = NULL)

df = df[order(-df$Mean),]
head(dfmedimp)


subsetted_resistome_group.ps<- subset_samples(group.ps, CollectionPhase=="WORKDAY_START" | CollectionPhase=="WORKDAY_END" | CollectionPhase=="POST_SHOWER" |CollectionPhase=="SWINE" | CollectionPhase=="ENVIRONMENT")

subsetted_resistome_all.ps<- subset_samples(RESISTOME_physeq, CollectionPhase=="WORKDAY_START" | CollectionPhase=="WORKDAY_END" | CollectionPhase=="POST_SHOWER" |CollectionPhase=="SWINE" | CollectionPhase=="ENVIRONMENT")


##MICROBIOTAPROCESS RAREFACTION

mpse<- subsetted_resistome_group.ps %>% as.MPSE()
mpse
data(mpse)
mpse %<>% 
mp_cal_rarecurve(
    .abundance = Abundance,
    chunks = 100, action = "add", force=TRUE
  )
mpse.data<-as.data.frame(mpse)

set.seed(1234)
plotmpse<- mpse %>% 
    mp_plot_rarecurve(
.rare = AbundanceRarecurve,   
      .alpha = "Observe", 
    .group = CollectionPhase,
    plot.group = TRUE
  ) + theme_bw() + coord_cartesian(xlim=c(0,NA), ylim = c(0,NA))+
  scale_color_manual(values=c("#789ED0", "#F79D40", "#72BF56", "#E4655A", "#AB8BCB")) +
  scale_fill_manual(values=c("#789ED0", "#F79D40", "#72BF56", "#E4655A", "#AB8BCB"))  
plotmpse


rare.level <- min(sample_sums(subsetted_resistome_group.ps))
rareres <- get_rarecurve(obj=subsetted_resistome_all.ps, chunks=50)
prare1 <- ggrarecurve(obj=rareres, factorNames="CollectionPhase", indexNames=c("Observe"), se=F, shadow = FALSE) + scale_fill_manual(values=c("#789ED0", "#F79D40", "#72BF56", "#E4655A", "#AB8BCB"))+scale_color_manual(values=c("#789ED0", "#F79D40", "#72BF56", "#E4655A", "#AB8BCB"))+
  theme_bw()+ theme(axis.text=element_text(size=8), panel.grid=element_blank(), strip.background = element_rect(colour=NA,fill="grey"), strip.text = element_blank()) + xlab("Read depth") + ylab("Unique number of ARG groups")+ theme(strip.text = element_blank()) + coord_cartesian(xlim=c(0,10000000))
prare1


p <- ggrare(subsetted_resistome_all.ps, sample, step = 10, color = "CollectionPhase", label = NULL, se = TRUE)
p <- p + facet_wrap(~CollectionPhase, nrow = 1) + theme_bw() + ylab("Genus richness")+  scale_fill_manual(values=c("#789ED0", "#F79D40", "#72BF56", "#E4655A", "#AB8BCB"))+ scale_color_manual(values=c("#789ED0", "#F79D40", "#72BF56", "#E4655A", "#AB8BCB")) 
plot(p)


################################################################################
#                             ALPHA DIVERSITY                                  #      
################################################################################
group.ps
###################All ARGs GROUP LEVEL###########
alpha_div <- estimate_richness(group.ps, measures = c("Observed","Shannon","Simpson","InvSimpson"))

alpha_div
alpha_div.df <- as(sample_data(group.ps), "data.frame")
alpha_div_meta<- cbind(alpha_div, alpha_div.df)
alpha_div_meta

ggplot(alpha_div_meta, aes(x=CollectionPhase, y = Observed, fill = CollectionPhase)) + theme_bw() + labs(title= "", y= "ARG group richness", x= "Collection phase") + geom_jitter(width = 0.2, shape=21, color= "black", size= 3) +
  geom_boxplot(alpha=0.3, shape=31, color= "black", outlier.colour = 'white', width=0.6, coef=5000, fatten=1, lwd=1) + ylim(8,NA) + scale_fill_tableau(palette = "Classic 10 Medium", type = "regular",direction = 1) + scale_color_tableau(palette = "Classic 10 Medium", type = "regular",direction = 1) + 
  scale_fill_discrete(name= "CollectionPhase", labels= c("WORKDAY_START"= "Workday start", "WORKDAY_END"="Workday end", "POST_SHOWER"="Post-shower", "SWINE"="Swine", "ENVIRONMENT" = "Environment")) + 
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
        panel.grid.minor.y = element_blank()) +   scale_x_discrete(limits = c("WORKDAY_START", "WORKDAY_END","POST_SHOWER","SWINE", "ENVIRONMENT"), labels = c("Workday start","Workday end","Post-shower","Swine", "Environment")) 


#richness_stats <- pairwise.wilcox.test(alpha_div_meta$Observed, alpha_div_meta$CollectionPhase, p.adjust.method = "BH")
#richness_stats$p.value
#write.csv(richness_stats[["p.value"]],"richness_stats.csv")

all_rich_model <-subset(alpha_div.df, !CollectionPhase %in% c("MockComm", "NegCtrl"))
allmodel_richness<-
  lmer(alpha_div_meta$Observed ~ CollectionPhase + Shotgun_Post_Host_Rem + Post_captureSampleConc_ng_ul + (1|Worker), data = all_rich_model)
summary(allmodel_richness)
anova(allmodel_richness, type='III',test="F")
lsmeans(allmodel_richness, pairwise~CollectionPhase, adjust="Tukey")
VarCorr(allmodel_richness)
#$contrasts
#contrast                    estimate   SE df t.ratio p.value
#ENVIRONMENT - POST_SHOWER     -143.9 57.7 42  -2.496  0.1109
#ENVIRONMENT - SWINE           -111.5 56.2 42  -1.985  0.2910
#ENVIRONMENT - WORKDAY_END     -195.9 60.9 42  -3.215  0.0201
#ENVIRONMENT - WORKDAY_START   -230.8 62.8 42  -3.674  0.0057
#POST_SHOWER - SWINE             32.4 35.9 42   0.902  0.8944
#POST_SHOWER - WORKDAY_END      -52.0 34.1 42  -1.526  0.5518
#POST_SHOWER - WORKDAY_START    -86.9 35.7 42  -2.433  0.1265
#SWINE - WORKDAY_END            -84.4 41.6 42  -2.031  0.2693
#SWINE - WORKDAY_START         -119.3 44.5 42  -2.681  0.0740
#WORKDAY_END - WORKDAY_START    -34.9 32.7 42  -1.067  0.8220


ggplot(alpha_div_meta, aes(x=CollectionPhase, y = Shannon, fill = CollectionPhase)) + theme_bw() + labs(title= "", y= "Shannon's diversity, all unique ARGs", x= "Collection phase") + geom_jitter(width = 0.2, shape=21, color= "black", size= 3) +
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
        panel.grid.minor.y = element_blank()) +  scale_x_discrete(limits = c("WORKDAY_START", "WORKDAY_END","POST_SHOWER","SWINE", "ENVIRONMENT"), labels = c("Workday start","Workday end","Post-shower","Swine", "Environment")) 

shannon_stats<- pairwise.wilcox.test(alpha_div_meta$Simpson, alpha_div_meta$CollectionPhase, p.adjust.method = "BH")
shannon_stats$p.value
write.csv(shannon_stats[["p.value"]],"shannon_stats.csv")

all_shannon_model <-subset(alpha_div.df, !CollectionPhase %in% c("MockComm", "NegCtrl"))
allmodel_shannon<-
  lmer(alpha_div_meta$Shannon ~ CollectionPhase + Shotgun_Post_Host_Rem +Post_captureSampleConc_ng_ul
       + (1|Worker), data = all_shannon_model, REML=F)
summary(allmodel_shannon)
anova(allmodel_shannon, type='III',test="F")
lsmeans(allmodel_shannon, pairwise~CollectionPhase, adjust="Tukey")
VarCorr(allmodel_shannon)

#$contrasts
#contrast                    estimate    SE   df t.ratio p.value
#ENVIRONMENT - POST_SHOWER    -0.3782 0.438 41.1  -0.864  0.9084
#ENVIRONMENT - SWINE          -0.5558 0.427 40.9  -1.301  0.6919
#ENVIRONMENT - WORKDAY_END    -0.3568 0.462 41.5  -0.773  0.9369
#ENVIRONMENT - WORKDAY_START  -0.8735 0.475 41.7  -1.838  0.3662
#POST_SHOWER - SWINE          -0.1776 0.272 41.6  -0.654  0.9650
#POST_SHOWER - WORKDAY_END     0.0213 0.241 25.4   0.089  1.0000
#POST_SHOWER - WORKDAY_START  -0.4953 0.254 26.9  -1.951  0.3162
#SWINE - WORKDAY_END           0.1990 0.313 42.0   0.636  0.9683
#SWINE - WORKDAY_START        -0.3177 0.334 42.0  -0.950  0.8755
#WORKDAY_END - WORKDAY_START  -0.5167 0.230 24.1  -2.242  0.1987

##################Medically important ARGs###########

medimp_alpha_div <- estimate_richness(medimp_group_subset.ps, measures = c("Observed","Shannon","Simpson","InvSimpson"))
medimp_alpha_div
medimp_alpha_div.df

medimp_alpha_div.df <- as(sample_data(medimp_group_subset.ps), "data.frame")
medimp_alpha_div_meta<- cbind(medimp_alpha_div, medimp_alpha_div.df)

medimp_alpha_div_meta


ggplot(medimp_alpha_div_meta, aes(x=CollectionPhase, y = Observed, fill = CollectionPhase)) + theme_bw() + labs(title= "", y= "Clinically important ARG group, richness", x= "Collection phase") + geom_jitter(width = 0.2, shape=21, color= "black", size= 3) +
  geom_boxplot(alpha=0.3, shape=31, color= "black", outlier.colour = 'white', width=0.6, coef=5000, fatten=1, lwd=1) + ylim(8,NA) + scale_fill_tableau(palette = "Classic 10 Medium", type = "regular",direction = 1) + scale_color_tableau(palette = "Classic 10 Medium", type = "regular",direction = 1) + 
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
panel.grid.minor.y = element_blank())  + scale_x_discrete(limits = c("WORKDAY_START", "WORKDAY_END","POST_SHOWER","SWINE", "ENVIRONMENT"), labels = c("Workday start","Workday end","Post-shower","Swine", "Environment")) + geom_text(data = tibble(x=1.75, y=8.5), size=3, aes(x=x, y=y, label= "Type III ANOVA P < 0.001", fontface="italic"),inherit.aes=FALSE) + geom_line(data=tibble(x=c(2,1), y=c(17.5,17.5)), aes(x=x, y=y), inherit.aes=FALSE)+ geom_line(data=tibble(x=c(2,3), y=c(17.75,17.75)), aes(x=x, y=y), inherit.aes=FALSE)+ geom_text(data = tibble(x=1.5, y=17.6), size=4, aes(x=x, y=y, label= "*"),inherit.aes=FALSE)+geom_text(data = tibble(x=2.5, y=17.8), size=4, aes(x=x, y=y, label= "*"),inherit.aes=FALSE)

richness_stats <- pairwise.wilcox.test(medimp_alpha_div_meta$Observed, medimp_alpha_div_meta$CollectionPhase, p.adjust.method = "BH")
richness_stats$p.value
write.csv(richness_stats[["p.value"]],"richness_stats.csv")

medimport_rich_model <-subset(medimp_alpha_div.df, !CollectionPhase %in% c("MockComm", "NegCtrl"))
medmodel_richness<-
  lmer(medimp_alpha_div_meta$Observed ~ CollectionPhase + (1|Worker), data = medimport_rich_model, REML=F)
summary(medmodel_richness)
anova(medmodel_richness, type='III',test="F")
lsmeans(medmodel_richness, pairwise~CollectionPhase, adjust="Tukey")
VarCorr(medmodel_richness)

#$contrasts
#contrast                    estimate    SE   df t.ratio p.value
#ENVIRONMENT - POST_SHOWER     -2.111 1.156 35.0  -1.827  0.3747
#ENVIRONMENT - SWINE           -1.111 1.125 34.3  -0.987  0.8593
#ENVIRONMENT - WORKDAY_END     -4.184 1.206 36.5  -3.469  0.0112
#ENVIRONMENT - WORKDAY_START   -2.497 1.262 36.9  -1.978  0.2966
#POST_SHOWER - SWINE            1.001 0.721 36.7   1.388  0.6392
#POST_SHOWER - WORKDAY_END     -2.073 0.557 20.0  -3.719  0.0105
#POST_SHOWER - WORKDAY_START   -0.385 0.603 24.0  -0.638  0.9672
#SWINE - WORKDAY_END           -3.074 0.808 39.4  -3.802  0.0042
#SWINE - WORKDAY_START         -1.386 0.898 39.6  -1.543  0.5418
#WORKDAY_END - WORKDAY_START    1.688 0.561 21.1   3.010  0.0467


ggplot(medimp_alpha_div_meta, aes(x=CollectionPhase, y = Shannon, fill = CollectionPhase)) + theme_bw() + labs(title= "", y= "Shannon's diversity, clinically important ARGs", x= "Collection phase") + geom_jitter(width = 0.2, shape=21, color= "black", size= 3) +
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
        panel.grid.minor.y = element_blank())+ scale_x_discrete(limits = c("WORKDAY_START", "WORKDAY_END","POST_SHOWER","SWINE", "ENVIRONMENT"), labels = c("Workday start","Workday end","Post-shower","Swine", "Environment"))

shannon_stats<- pairwise.wilcox.test(medimp_alpha_div_meta$Simpson, medimp_alpha_div_meta$CollectionPhase, p.adjust.method = "BH")
shannon_stats$p.value
write.csv(shannon_stats[["p.value"]],"shannon_stats.csv")

medimport_shannon_model <-subset(medimp_alpha_div.df, !CollectionPhase %in% c("MockComm", "NegCtrl"))
medmodel_shannon<-
  lmer(medimp_alpha_div_meta$Shannon ~ CollectionPhase + Shotgun_Post_Host_Rem +
         Post_captureSampleConc_ng_ul + (1|Worker), data = medimport_shannon_model, REML=F)
summary(medmodel_shannon)
anova(medmodel_shannon, type='III',test="F")
lsmeans(medmodel_shannon, pairwise~CollectionPhase, adjust="Tukey")
VarCorr(medmodel_shannon)

#$contrasts
#contrast                    estimate     SE   df t.ratio p.value
#ENVIRONMENT - POST_SHOWER   -0.29652 0.1941 31.7  -1.527  0.5530
#ENVIRONMENT - SWINE         -0.09181 0.1891 30.8  -0.486  0.9881
#ENVIRONMENT - WORKDAY_END   -0.35265 0.2025 33.8  -1.742  0.4233
#ENVIRONMENT - WORKDAY_START -0.34520 0.2119 34.5  -1.629  0.4895
#POST_SHOWER - SWINE          0.20471 0.1211 34.2   1.691  0.4527
#POST_SHOWER - WORKDAY_END   -0.05613 0.0919 15.8  -0.611  0.9713
#POST_SHOWER - WORKDAY_START -0.04868 0.0997 19.9  -0.488  0.9876
#SWINE - WORKDAY_END         -0.26084 0.1355 38.2  -1.926  0.3216
#SWINE - WORKDAY_START       -0.25339 0.1505 38.7  -1.684  0.4555
#WORKDAY_END - WORKDAY_START  0.00745 0.0925 16.9   0.080  1.0000

#ARGPrevalence_abundance:########################### WORKDAY_START
WORK_START.group.ps <- subset_samples(RESISTOME_physeq, CollectionPhase=="WORKDAY_START")
medimp_WORK_START.group.ps<- subset_samples(RESISTOME_physeq, CollectionPhase=="WORKDAY_START")
medimp_WORK_START.group.ps
medimp_WORK_START.group_subset.ps = subset_samples(medimp_WORK_START.group.ps, sample_names(medimp_WORK_START.group.ps) != "R1_19")
medimp_WORK_START.group_subset.ps

row

prevalence <- function(medimp_WORK_START.group_subset.ps, add_tax = TRUE){
  set.seed(87)                             
  
  ## Check if taxa are rows
  trows <- taxa_are_rows(medimp_WORK_START.group_subset.ps)
  
  ## Extract OTU table
  otutab <- as.data.frame(otu_table(medimp_WORK_START.group_subset.ps))
  otutab
  ## Transpose OTU table (species should be arranged by rows)
  if(trows == FALSE){
    otutab <- t(otutab)
  }
  ## Extract sequencing depth data
  
  WORK_START_DEPTH<- sample_data(medimp_WORK_START.group_subset.ps)   
  WORK_START_DEPTH.df<-as.data.frame(WORK_START_DEPTH)
  
  ## Estimate prevalence (number of samples with OTU present)
  prevdf <- apply(X = otutab,
                  #MARGIN = ifelse(trows, yes = 1, no = 2),  # for a non-transposed data
                  MARGIN = 1,
                  FUN = function(x){sum(((x > 0.01))+0.01)})
  ## Add total and average read counts per OTU
  prevdf <- data.frame(Prevalence = log10(prevdf)*10,
                       TotalAbundance = ((taxa_sums(phyloseq_transform_css(medimp_WORK_START.group_subset.ps))/(log10(WORK_START_DEPTH.df$Shotgun_Post_Host_Rem)))),
                       MeanAbundance = rowMeans(otutab),
                       MedianAbundance = apply(otutab, 1, median))
  
  ## Add taxonomy table
  if(add_tax == TRUE && !is.null(tax_table(medimp_WORK_START.group_subset.ps, errorIfNULL = F))){
    prevdf <- cbind(prevdf, tax_table(medimp_WORK_START.group_subset.ps))
  }
  return(prevdf)
}
phyloseq_prevalence_plot <- function(medimp_WORK_START.group_subset.ps, prev.trh = NULL, taxcolor = NULL, facet = FALSE, point_alpha = 0.3, showplot = T){
  
  require(ggplot2)
  
  ## Compute prevalence of each species
  prevdf <- prevalence(medimp_WORK_START.group_subset.ps)
  prevdf$PrevFrac <- with(prevdf, Prevalence / nsamples(medimp_WORK_START.group_subset.ps))
  
  ## Prepare a plot
  if(is.null(taxcolor)){ pp <- ggplot(prevdf, aes(x = TotalAbundance, y = PrevFrac)) }
  if(!is.null(taxcolor)){ pp <- ggplot(prevdf, aes_string(x = "TotalAbundance", y = "PrevFrac", color = taxcolor)) }
  
  ## Add prevalence threshold line
  if(!is.null(prev.trh)){ pp <- pp + geom_hline(yintercept = prev.trh, alpha = 0.7, linetype = 1, size=1, color= "#063970") }
  
  pp <- pp +
    geom_point(size = 3, alpha = point_alpha) +
    theme_bw() + theme(plot.margin = unit(c(0.1,0.4,0.4,0.4), "cm"))+
    xlab("Normalized abundance") + ylim(0,1.3)
  ylab("Prevalence [Frac. samples]") +
    theme(legend.position="none")
  
  if(facet == TRUE){ pp <- pp + facet_wrap(as.formula(paste("~", taxcolor))) }
  if(showplot == TRUE){ print(pp) }
  invisible(pp)
}

subWORK_START.group.ps<-subset_taxa(medimp_WORK_START.group_subset.ps, Group=="CTX" | 
              Group=="GES" |
              Group=="IMI" |
              Group=="KPC" |
              Group=="SHV" |
              Group=="TEM" |
              Group=="IMP" |
              Group=="NDM" |
              Group=="CMY" |
              Group=="OXA" |
              Group=="MEC" |
              Group=="MCR" |
              Group=="VAT" |
              Group=="VGB" |
              Group=="CFR" |
              Group=="VGA" |
              Group=="SME" |
              Group=="AAC6-PRIME" |
              Group=="BLAZ" |
              Group=="NDM" |
              Group=="VIM" |
              Group=="MCR" |
              Group=="ERMB" |
              Group=="QNRA" |
              Group=="QNRB" |
              Group=="QNRA" |
              Group=="TETM" |
              Group=="DFRA" |
              Group=="VANYB" |
              Group=="VANYD" |
              Group=="VANYA" |
              Group=="SULI")

phyloseq_prevalence_plot_medimp <- function(subWORK_START.group.ps, prev.trh = NULL, taxcolor = NULL, facet = TRUE, point_alpha = 0.3, showplot = T){
  
  require(ggplot2)
  
  ## Compute prevalence of each species
  prevdf <- prevalence(subWORK_START.group.ps)
  prevdf$PrevFrac <- with(prevdf, Prevalence / nsamples(subWORK_START.group.ps))
  
  ## Prepare a plot
  if(is.null(taxcolor)){ pp <- ggplot(prevdf, aes(x = TotalAbundance, y = PrevFrac)) }
  if(!is.null(taxcolor)){ pp <- ggplot(prevdf, aes_string(x = "TotalAbundance", y = "PrevFrac", color = taxcolor)) }
  
  ## Add prevalence threshold line
  if(!is.null(prev.trh)){ pp <- pp + geom_hline(yintercept = prev.trh, alpha = 0.7, linetype = 1, size=1, color= "#063970") }
  
  pp <- pp +
    geom_point(size = 3, alpha = point_alpha) + #stat_smooth(method = "lm",formula = y ~ x, geom = "smooth", se=FALSE, span=0.8, fullrange=TRUE, color="black", size=0.6)+ 
    geom_point(shape = 1,size = 3,colour = "black", stroke=0.01) + theme_bw() +
    theme_bw() + theme(plot.margin = unit(c(0.1,0.4,0.4,0.4), "cm"))+ ylim(0,1.05)+
    xlim(0,15)+
    xlab("CSS Normalized allele abundance") +
    ylab("Prevalence [Frac. samples]") +
    theme(legend.position="none")
  
  if(facet == TRUE){ pp <- pp + facet_wrap(as.formula(paste("~", taxcolor))) }
  if(showplot == TRUE){ print(pp) }
  invisible(pp)
}


Prevplot_WORKDAY_START<- #phyloseq_prevalence_plot(WORK_START.group.ps, taxcolor = "Type", facet =TRUE, point_alpha = 0.4, prev.trh = 0.05) + 
  phyloseq_prevalence_plot_medimp(subWORK_START.group.ps, taxcolor="Class",facet = TRUE, point_alpha = 0.3, prev.trh = 0.05)
Prevplot_WORKDAY_START

Prevplot_WORKDAY_START.df<- Prevplot_WORKDAY_START$data
Prevplot_WORKDAY_START.df


#ARGPrevalence_abundance:########################### WORKDAY_END
WORK_END.group.ps <- subset_samples(RESISTOME_physeq, CollectionPhase=="WORKDAY_END")
medimp_WORK_END.group.ps<- subset_samples(RESISTOME_physeq, CollectionPhase=="WORKDAY_END")
medimp_WORK_END.group.ps
medimp_WORK_END.group_subset.ps = subset_samples(medimp_WORK_END.group.ps, sample_names(medimp_WORK_END.group.ps) != "R1_19")

prevalence <- function(WORK_END.group.ps, add_tax = TRUE){
  set.seed(87)                             
  
  ## Check if taxa are rows
  trows <- taxa_are_rows(WORK_END.group.ps)
  
  ## Extract OTU table
  otutab <- as.data.frame(otu_table(WORK_END.group.ps))
  otutab
  ## Transpose OTU table (species should be arranged by rows)
  if(trows == FALSE){
    otutab <- t(otutab)
  }

  ## Extract sequencing depth data
  
 
  WORK_END_DEPTH<- sample_data(medimp_WORK_END.group_subset.ps)   
  WORK_END_DEPTH.df<-as.data.frame(WORK_END_DEPTH)
  
    ## Estimate prevalence (number of samples with OTU present)
  prevdf <- apply(X = otutab,
                  #MARGIN = ifelse(trows, yes = 1, no = 2),  # for a non-transposed data
                  MARGIN = 1,
                  FUN = function(x){sum(((x > 0.01))+0.01)})
## Add total and average read counts per OTU
  prevdf <- data.frame(Prevalence = log10(prevdf)*10,
                       TotalAbundance = ((taxa_sums(phyloseq_transform_css(WORK_END.group.ps))/(log10(WORK_END_DEPTH$Shotgun_Post_Host_Rem)))),
                       MeanAbundance = rowMeans(otutab),
                       MedianAbundance = apply(otutab, 1, median))
  
  ## Add taxonomy table
  if(add_tax == TRUE && !is.null(tax_table(WORK_END.group.ps, errorIfNULL = F))){
    prevdf <- cbind(prevdf, tax_table(WORK_END.group.ps))
  }
  return(prevdf)
}
phyloseq_prevalence_plot <- function(WORK_END.group.ps, prev.trh = NULL, taxcolor = NULL, facet = FALSE, point_alpha = 0.3, showplot = T){
  
  require(ggplot2)
  
  ## Compute prevalence of each species
  prevdf <- prevalence(WORK_END.group.ps)
  prevdf$PrevFrac <- with(prevdf, Prevalence / nsamples(WORK_END.group.ps))
  
  ## Prepare a plot
  if(is.null(taxcolor)){ pp <- ggplot(prevdf, aes(x = TotalAbundance, y = PrevFrac)) }
  if(!is.null(taxcolor)){ pp <- ggplot(prevdf, aes_string(x = "TotalAbundance", y = "PrevFrac", color = taxcolor)) }
  
  ## Add prevalence threshold line
  if(!is.null(prev.trh)){ pp <- pp + geom_hline(yintercept = prev.trh, alpha = 0.7, linetype = 1, size=1, color= "#063970") }
  
  pp <- pp +
    geom_point(size = 3, alpha = point_alpha) +
    theme_bw() + theme(plot.margin = unit(c(0.1,0.4,0.4,0.4), "cm"))+
    xlab("Normalized abundance") + ylim(0,1.3)
  ylab("Prevalence [Frac. samples]") +
    theme(legend.position="none")
  
  if(facet == TRUE){ pp <- pp + facet_wrap(as.formula(paste("~", taxcolor))) }
  if(showplot == TRUE){ print(pp) }
  invisible(pp)
}

subWORK_END.group.ps<-subset_taxa(medimp_WORK_END.group_subset.ps, Group=="CTX" | 
                                      Group=="GES" |
                                      Group=="IMI" |
                                      Group=="KPC" |
                                      Group=="SHV" |
                                      Group=="TEM" |
                                      Group=="IMP" |
                                      Group=="NDM" |
                                      Group=="CMY" |
                                      Group=="OXA" |
                                      Group=="MEC" |
                                      Group=="MCR" |
                                      Group=="VAT" |
                                      Group=="VGB" |
                                      Group=="CFR" |
                                      Group=="VGA" |
                                      Group=="SME" |
                                      Group=="AAC6-PRIME" |
                                      Group=="BLAZ" |
                                      Group=="NDM" |
                                      Group=="VIM" |
                                      Group=="MCR" |
                                      Group=="ERMB" |
                                      Group=="QNRA" |
                                      Group=="QNRB" |
                                      Group=="QNRA" |
                                      Group=="TETM" |
                                      Group=="DFRA" |
                                      Group=="VANYB" |
                                      Group=="VANYD" |
                                      Group=="VANYA" |
                                      Group=="SULI")
subWORK_END.group.ps
phyloseq_prevalence_plot_medimp <- function(subWORK_END.group.ps, prev.trh = NULL, taxcolor = NULL, facet = TRUE, point_alpha = 0.3, showplot = T){
  
  require(ggplot2)
  
  ## Compute prevalence of each species
  prevdf <- prevalence(subWORK_END.group.ps)
  prevdf$PrevFrac <- with(prevdf, Prevalence / nsamples(subWORK_END.group.ps))
  
  ## Prepare a plot
  if(is.null(taxcolor)){ pp <- ggplot(prevdf, aes(x = TotalAbundance, y = PrevFrac)) }
  if(!is.null(taxcolor)){ pp <- ggplot(prevdf, aes_string(x = "TotalAbundance", y = "PrevFrac", color = taxcolor)) }
  
  ## Add prevalence threshold line
  if(!is.null(prev.trh)){ pp <- pp + geom_hline(yintercept = prev.trh, alpha = 0.7, linetype = 1, size=1, color= "#063970") }
  
  pp <- pp +
    geom_point(size = 3, alpha = point_alpha) + #stat_smooth(method = "lm", formula = y ~ x, geom = "smooth", se=FALSE, span=0.8, fullrange=TRUE, color="black", size=0.6)+ 
    geom_point(shape = 1,size = 3,colour = "black", stroke=0.01) + theme_bw() + theme(plot.margin = unit(c(0.1,0.4,0.4,0.4), "cm"))+ ylim(0,1.05)+ xlim(0,15)+
    xlab("CSS Normalized allele abundance") +
    ylab("Prevalence [Frac. samples]") +
    theme(legend.position="none")
  
  if(facet == TRUE){ pp <- pp + facet_wrap(as.formula(paste("~", taxcolor))) }
  if(showplot == TRUE){ print(pp) }
  invisible(pp)
}


Prevplot_WORKDAY_END<- #phyloseq_prevalence_plot(WORK_END.group.ps, taxcolor = "Type", facet =TRUE, point_alpha = 0.4, prev.trh = 0.05) + 
  phyloseq_prevalence_plot_medimp(subWORK_END.group.ps, taxcolor="Class",facet = TRUE, point_alpha = 0.3, prev.trh = 0.05)

#ARGPrevalence_abundance:########################### POST_SHOWER
POST_SHOWER.group.ps <- subset_samples(RESISTOME_physeq, CollectionPhase=="POST_SHOWER")
medimp_POST_SHOWER.group.ps<- subset_samples(RESISTOME_physeq, CollectionPhase=="POST_SHOWER")
medimp_POST_SHOWER.group.ps
medimp_POST_SHOWER.group.ps.group_subset.ps = subset_samples(medimp_POST_SHOWER.group.ps, sample_names(medimp_POST_SHOWER.group.ps) != "R1_19")
medimp_POST_SHOWER.group.ps.group_subset.ps
prevalence <- function(medimp_POST_SHOWER.group.ps.group_subset.ps, add_tax = TRUE){
  set.seed(87)                             
  
  ## Check if taxa are rows
  trows <- taxa_are_rows(medimp_POST_SHOWER.group.ps.group_subset.ps)
  
  ## Extract OTU table
  otutab <- as.data.frame(otu_table(medimp_POST_SHOWER.group.ps.group_subset.ps))
  otutab
  ## Transpose OTU table (species should be arranged by rows)
  if(trows == FALSE){
    otutab <- t(otutab)
  }
  POST_SHOWER_DEPTH<- sample_data(medimp_POST_SHOWER.group.ps.group_subset.ps)   
  POST_SHOWER_DEPTH.df<-as.data.frame(POST_SHOWER_DEPTH)
  
  ## Estimate prevalence (number of samples with OTU present)
  prevdf <- apply(X = otutab,
                  #MARGIN = ifelse(trows, yes = 1, no = 2),  # for a non-transposed data
                  MARGIN = 1,
                  FUN = function(x){sum(((x > 0.01))+0.01)})
  ## Add total and average read counts per OTU
  prevdf <- data.frame(Prevalence = log10(prevdf)*10,
                       TotalAbundance = ((taxa_sums(phyloseq_transform_css(medimp_POST_SHOWER.group.ps.group_subset.ps))/(log10(POST_SHOWER_DEPTH.df$Shotgun_Post_Host_Rem)))),
                       MeanAbundance = rowMeans(otutab),
                       MedianAbundance = apply(otutab, 1, median))
  
  ## Add taxonomy table
  if(add_tax == TRUE && !is.null(tax_table(medimp_POST_SHOWER.group.ps.group_subset.ps, errorIfNULL = F))){
    prevdf <- cbind(prevdf, tax_table(medimp_POST_SHOWER.group.ps.group_subset.ps))
  }
  return(prevdf)
}
phyloseq_prevalence_plot <- function(medimp_POST_SHOWER.group.ps.group_subset.ps, prev.trh = NULL, taxcolor = NULL, facet = FALSE, point_alpha = 0.3, showplot = T){
  
  require(ggplot2)
  
  ## Compute prevalence of each species
  prevdf <- prevalence(medimp_POST_SHOWER.group.ps.group_subset.ps)
  prevdf$PrevFrac <- with(prevdf, Prevalence / nsamples(medimp_POST_SHOWER.group.ps.group_subset.ps))
  
  ## Prepare a plot
  if(is.null(taxcolor)){ pp <- ggplot(prevdf, aes(x = TotalAbundance, y = PrevFrac)) }
  if(!is.null(taxcolor)){ pp <- ggplot(prevdf, aes_string(x = "TotalAbundance", y = "PrevFrac", color = taxcolor)) }
  
  ## Add prevalence threshold line
  if(!is.null(prev.trh)){ pp <- pp + geom_hline(yintercept = prev.trh, alpha = 0.7, linetype = 1, size=1, color= "#063970") }
  
  pp <- pp +
    geom_point(size = 3, alpha = point_alpha) +
    theme_bw() + theme(plot.margin = unit(c(0.1,0.4,0.4,0.4), "cm"))+
    xlab("Normalized abundance") + ylim(0,1.3)
  ylab("Prevalence [Frac. samples]") +
    theme(legend.position="none")
  
  if(facet == TRUE){ pp <- pp + facet_wrap(as.formula(paste("~", taxcolor))) }
  if(showplot == TRUE){ print(pp) }
  invisible(pp)
}

subPOST_SHOWER.group.ps<-subset_taxa(medimp_POST_SHOWER.group.ps.group_subset.ps, Group=="CTX" | 
                                    Group=="GES" |
                                    Group=="IMI" |
                                    Group=="KPC" |
                                    Group=="SHV" |
                                    Group=="TEM" |
                                    Group=="IMP" |
                                    Group=="NDM" |
                                    Group=="CMY" |
                                    Group=="OXA" |
                                    Group=="MEC" |
                                    Group=="MCR" |
                                    Group=="VAT" |
                                    Group=="VGB" |
                                    Group=="CFR" |
                                    Group=="VGA" |
                                    Group=="SME" |
                                    Group=="AAC6-PRIME" |
                                    Group=="BLAZ" |
                                    Group=="NDM" |
                                    Group=="VIM" |
                                    Group=="MCR" |
                                    Group=="ERMB" |
                                    Group=="QNRA" |
                                    Group=="QNRB" |
                                    Group=="QNRA" |
                                    Group=="TETM" |
                                    Group=="DFRA" |
                                    Group=="VANYB" |
                                    Group=="VANYD" |
                                    Group=="VANYA" |
                                    Group=="SULI")
subPOST_SHOWER.group.ps
phyloseq_prevalence_plot_medimp <- function(subPOST_SHOWER.group.ps, prev.trh = NULL, taxcolor = NULL, facet = TRUE, point_alpha = 0.3, showplot = T){
  
  require(ggplot2)
  
  ## Compute prevalence of each species
  prevdf <- prevalence(subPOST_SHOWER.group.ps)
  prevdf$PrevFrac <- with(prevdf, Prevalence / nsamples(subPOST_SHOWER.group.ps))
  
  ## Prepare a plot
  if(is.null(taxcolor)){ pp <- ggplot(prevdf, aes(x = TotalAbundance, y = PrevFrac)) }
  if(!is.null(taxcolor)){ pp <- ggplot(prevdf, aes_string(x = "TotalAbundance", y = "PrevFrac", color = taxcolor)) }
  
  ## Add prevalence threshold line
  if(!is.null(prev.trh)){ pp <- pp + geom_hline(yintercept = prev.trh, alpha = 0.7, linetype = 1, size=1, color= "#063970") }
  
  pp <- pp +
    geom_point(size = 3, alpha = point_alpha) + 
    #stat_smooth(method = "lm", formula = y ~ x, geom = "smooth", se=FALSE, span=0.8, fullrange=TRUE, color="black", size=0.6) + 
    geom_point(shape = 1,size = 3,colour = "black", stroke=0.01) + theme_bw() +
    theme_bw() + theme(plot.margin = unit(c(0.1,0.4,0.4,0.4), "cm"))+ ylim(0,1.05)+xlim(0,15)+
    xlab("CSS Normalized allele abundance") +
    ylab("Prevalence [Frac. samples]") +
    theme(legend.position="none")
  
  if(facet == TRUE){ pp <- pp + facet_wrap(as.formula(paste("~", taxcolor))) }
  if(showplot == TRUE){ print(pp) }
  invisible(pp)
}


Prevplot_POST_SHOWER<- #phyloseq_prevalence_plot(WORK_END.group.ps, taxcolor = "Type", facet =TRUE, point_alpha = 0.4, prev.trh = 0.05) + 
  phyloseq_prevalence_plot_medimp(subPOST_SHOWER.group.ps, taxcolor="Class",facet = TRUE, point_alpha = 0.3, prev.trh = 0.05)

#ARGPrevalence_abundance:########################### SWINE
SWINE.group.ps <- subset_samples(RESISTOME_physeq, CollectionPhase=="SWINE")
medimp_SWINE.group.ps<- subset_samples(RESISTOME_physeq, CollectionPhase=="SWINE")
medimp_SWINE.group.ps


prevalence <- function(medimp_SWINE.group.ps, add_tax = TRUE){
  set.seed(87)                             
  
  ## Check if taxa are rows
  trows <- taxa_are_rows(medimp_SWINE.group.ps)
  
  ## Extract OTU table
  otutab <- as.data.frame(otu_table(medimp_SWINE.group.ps))
  otutab
  ## Transpose OTU table (species should be arranged by rows)
  if(trows == FALSE){
    otutab <- t(otutab)
  }
  
  SWINE_DEPTH<- sample_data(medimp_SWINE.group.ps)   
  SWINE_DEPTH.df<-as.data.frame(SWINE_DEPTH)
  
  ## Estimate prevalence (number of samples with OTU present)
  prevdf <- apply(X = otutab,
                  #MARGIN = ifelse(trows, yes = 1, no = 2),  # for a non-transposed data
                  MARGIN = 1,
                  FUN = function(x){sum(((x > 0.01))+0.01)})
  ## Add total and average read counts per OTU
  prevdf <- data.frame(Prevalence = log10(prevdf)*10,
                       TotalAbundance = ((taxa_sums(phyloseq_transform_css(medimp_SWINE.group.ps))/(log10(SWINE_DEPTH.df$Shotgun_Post_Host_Rem)))),
                       MeanAbundance = rowMeans(otutab),
                       MedianAbundance = apply(otutab, 1, median))
  
  ## Add taxonomy table
  if(add_tax == TRUE && !is.null(tax_table(medimp_SWINE.group.ps, errorIfNULL = F))){
    prevdf <- cbind(prevdf, tax_table(medimp_SWINE.group.ps))
  }
  return(prevdf)
}
phyloseq_prevalence_plot <- function(medimp_SWINE.group.ps, prev.trh = NULL, taxcolor = NULL, facet = FALSE, point_alpha = 0.3, showplot = T){
  
  require(ggplot2)
  
  ## Compute prevalence of each species
  prevdf <- prevalence(medimp_SWINE.group.ps)
  prevdf$PrevFrac <- with(prevdf, Prevalence / nsamples(medimp_SWINE.group.ps))
  
  ## Prepare a plot
  if(is.null(taxcolor)){ pp <- ggplot(prevdf, aes(x = TotalAbundance, y = PrevFrac)) }
  if(!is.null(taxcolor)){ pp <- ggplot(prevdf, aes_string(x = "TotalAbundance", y = "PrevFrac", color = taxcolor)) }
  
  ## Add prevalence threshold line
  if(!is.null(prev.trh)){ pp <- pp + geom_hline(yintercept = prev.trh, alpha = 0.7, linetype = 1, size=1, color= "#063970") }
  
  pp <- pp +
    geom_point(size = 3, alpha = point_alpha) +
    theme_bw() + theme(plot.margin = unit(c(0.1,0.4,0.4,0.4), "cm"))+
    xlab("Normalized abundance") + ylim(0,1.3)
  ylab("Prevalence [Frac. samples]") +
    theme(legend.position="none")
  
  if(facet == TRUE){ pp <- pp + facet_wrap(as.formula(paste("~", taxcolor))) }
  if(showplot == TRUE){ print(pp) }
  invisible(pp)
}

subSWINE.group.ps<-subset_taxa(medimp_SWINE.group.ps, Group=="CTX" | 
                                       Group=="GES" |
                                       Group=="IMI" |
                                       Group=="KPC" |
                                       Group=="SHV" |
                                       Group=="TEM" |
                                       Group=="IMP" |
                                       Group=="NDM" |
                                       Group=="CMY" |
                                       Group=="OXA" |
                                       Group=="MEC" |
                                       Group=="MCR" |
                                       Group=="VAT" |
                                       Group=="VGB" |
                                       Group=="CFR" |
                                       Group=="VGA" |
                                       Group=="SME" |
                                       Group=="AAC6-PRIME" |
                                       Group=="BLAZ" |
                                       Group=="NDM" |
                                       Group=="VIM" |
                                       Group=="MCR" |
                                       Group=="ERMB" |
                                       Group=="QNRA" |
                                       Group=="QNRB" |
                                       Group=="QNRA" |
                                       Group=="TETM" |
                                       Group=="DFRA" |
                                       Group=="VANYB" |
                                       Group=="VANYD" |
                                       Group=="VANYA" |
                                       Group=="SULI")
subSWINE.group.ps
phyloseq_prevalence_plot_medimp <- function(subSWINE.group.ps, prev.trh = NULL, taxcolor = NULL, facet = TRUE, point_alpha = 0.3, showplot = T){
  
  require(ggplot2)
  
  ## Compute prevalence of each species
  prevdf <- prevalence(subSWINE.group.ps)
  prevdf$PrevFrac <- with(prevdf, Prevalence / nsamples(subSWINE.group.ps))
  
  ## Prepare a plot
  if(is.null(taxcolor)){ pp <- ggplot(prevdf, aes(x = TotalAbundance, y = PrevFrac)) }
  if(!is.null(taxcolor)){ pp <- ggplot(prevdf, aes_string(x = "TotalAbundance", y = "PrevFrac", color = taxcolor)) }
  
  ## Add prevalence threshold line
  if(!is.null(prev.trh)){ pp <- pp + geom_hline(yintercept = prev.trh, alpha = 0.7, linetype = 1, size=1, color= "#063970") }
  
  pp <- pp +
    geom_point(size = 3, alpha = point_alpha) + #stat_smooth(method = "lm", formula = y ~ x,geom = "smooth", se=FALSE, span=0.8, fullrange=TRUE, color="black", size=0.6)+ 
    geom_point(shape = 1,size = 3,colour = "black", stroke=0.01) +
    theme_bw() + theme(plot.margin = unit(c(0.1,0.4,0.4,0.4), "cm"))+ ylim(0,1.05)+ xlim(0,15.1)+
xlab("CSS Normalized allele abundance") +
    ylab("Prevalence [Frac. samples]") +
    theme(legend.position="none")
  
  if(facet == TRUE){ pp <- pp + facet_wrap(as.formula(paste("~", taxcolor))) }
  if(showplot == TRUE){ print(pp) }
  invisible(pp)
}


Prevplot_SWINE<- #phyloseq_prevalence_plot(WORK_END.group.ps, taxcolor = "Type", facet =TRUE, point_alpha = 0.4, prev.trh = 0.05) + 
  phyloseq_prevalence_plot_medimp(subSWINE.group.ps, taxcolor="Class",facet = TRUE, point_alpha = 0.3, prev.trh = 0.05)


patchedPrevAbundPlot<-plot_grid(Prevplot_WORKDAY_START, Prevplot_WORKDAY_END, Prevplot_POST_SHOWER, Prevplot_SWINE,labels=c('A. Workday start', 'B. Workday end', 'C. Post-shower', 'D. Swine'), align= "h", label_size=14)
          
patchedPrevAbundPlot

################################################################################
#                           BETA DIVERSITY                                     #      
################################################################################

##All ARG groups

group.ps
ALLARGCounts <- as.data.frame(phyloseq::otu_table(group.ps))
if(phyloseq::taxa_are_rows(group.ps) == FALSE){ otus <- t(ALLARGCounts) }
ALLARGCounts

#zCompositions component of the anlaysis, using bayesian approach:


ALLARG_cmultRepl<-cmultRepl(ALLARGCounts,  label = 0, method = c("GBM"), output = c("p-counts"), frac = 0.2, threshold = 0.5, adjust = TRUE, t= NULL, s= NULL, z.warning = 0.8, suppress.print = FALSE)

ALLARG_cmultRepl

#Now we must reconvert this OTU table to the phyloseq object

group.ps <-phyloseq::otu_table(group.ps) <- phyloseq::otu_table(ALLARG_cmultRepl, taxa_are_rows = T)

return(IMPUTEDmobtype_PLASMID.ps)
group.ps

#Create a new phyloseq object to include the imputed virus OTU table

IMPALL_RESISTOME = otu_table(group.ps, taxa_are_rows = T)
GENES = tax_table(gen_mat)
ARG_SAMPLES = sample_data(resistomesamples_df)
IMPUTED_ALLRESISTOME_physeq <- phyloseq(IMPALL_RESISTOME, GENES, ARG_SAMPLES)
IMPUTED_ALLRESISTOME_physeq

#Now let's retry the ordination w/ this new phyloseq object

sum(taxa_sums(IMPUTED_ALLRESISTOME_physeq)==0) # 0 taxa not present across all samples
IMPUTED_ALLARGS<- prune_taxa(taxa_sums(IMPUTED_ALLRESISTOME_physeq) > 0, IMPUTED_ALLRESISTOME_physeq)#prunning low-abundance taxa of gene groups

ALLARGS_physeq.css <- phyloseq_transform_css(IMPUTED_ALLARGS)

ALLARGS_physeq.dist <- vegdist(decostand(t(otu_table(IMPUTED_ALLARGS)),"rclr"), method = "euclidean")
set.seed(1999)
ALLARGS_physeq.dist
ALLARGS_physeq.ord<- ordinate(ALLARGS_physeq.css, method = "PCoA", distance = "euclidean")

ALLARGS_plot<-plot_ordination(ALLARGS_physeq.css,ALLARGS_physeq.ord, type = "samples", color="CollectionPhase")
ALLARGS_plot

ALLARGcentroid<- ALLARGS_plot$data %>% group_by(CollectionPhase) %>% summarise(Axis.1=mean(Axis.1), Axis.2=mean(Axis.2), .groups = "drop")

#stressplot(plot_ord_all)
#str(plot_ord_all)

ALLARG_ordplot <- ggplot(ALLARGS_plot$data,ALLARGS_plot$mapping) + stat_ellipse(geom = "polygon", type= "t", level = 0.85, lty =1, lwd = 1.3, alpha= 0.07, show.legend =FALSE) + theme_bw() +geom_point(aes(fill=CollectionPhase), size=6, alpha=0.3, stroke=1.0,) + geom_point(data=ALLARGcentroid, size=13, shape=16, fill="black",mapping= aes(colour=CollectionPhase), stroke=4, alpha= 0.8, show.legend = FALSE) + scale_color_manual(name="Collection phase", breaks = c("WORKDAY_START", "WORKDAY_END", "POST_SHOWER", "SWINE", "ENVIRONMENT"), values = c("#789ED0", "#F79D40", "#72BF56", "#E4655A", "#AB8BCB"), labels= c("WORKDAY_START"= "Workday start", "WORKDAY_END"="Workday end", "POST_SHOWER"="Post-shower", "SWINE"="Swine", "ENVIRONMENT"="Environment")) +
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
        axis.title.x = element_text(size = 18, vjust = -1.5)) + coord_fixed(xlim=c(-50,50), ylim = c(-50,50))  
ALLARG_ordplot



ALLARG_metadata<- sample_data(IMPUTED_ALLARGS)
pmeta<-as.data.frame(ALLARG_metadata)
pmeta

allarg.anosim<- anosim(ALLARGS_physeq.dist, pmeta$CollectionPhase, perm = 1000)
allarg.anosim
#ANOSIM statistic R: 0.3562 
#Significance: 0.000999 
#

allarg.adonis <- pairwise.adonis(ALLARGS_physeq.dist, pmeta$CollectionPhase, perm = 1000, p.adjust.m = "BH")
allarg.adonis
summary(allarg.adonis)

#> allarg.adonis
#                           pairs Df SumsOfSqs  F.Model         R2    p.value p.adjusted sig
#1  WORKDAY_START vs POST_SHOWER  1  5271.285 1.804409 0.09111149 0.00079992 0.00199980   *
#2        WORKDAY_START vs SWINE  1 36787.206 9.401556 0.34310301 0.00009999 0.00033330  **
#3  WORKDAY_START vs ENVIRONMENT  1  7349.079 2.103744 0.17380939 0.01379862 0.01724828   .
#4  WORKDAY_START vs WORKDAY_END  1  9289.303 2.194412 0.10866433 0.00209979 0.00349965   *
#5          POST_SHOWER vs SWINE  1 27257.231 8.421394 0.31873391 0.00009999 0.00033330  **
#6    POST_SHOWER vs ENVIRONMENT  1  5491.456 2.412622 0.19436842 0.01319868 0.01724828   .
#7    POST_SHOWER vs WORKDAY_END  1  6977.389 1.961625 0.09826983 0.00179982 0.00349965   *
#8          SWINE vs ENVIRONMENT  1 10423.997 2.566902 0.20425895 0.01579842 0.01755380   .
#9          SWINE vs WORKDAY_END  1 23514.201 5.169670 0.22312230 0.00009999 0.00033330  **
#10   ENVIRONMENT vs WORKDAY_END  1  7653.776 1.650437 0.14166308 0.18548145 0.18548145   


allarg.disper <- betadisper(ALLARGS_physeq.dist, pmeta$CollectionPhase)
plot(allarg.disper)

allarg.permdisp<- permutest(allarg.disper, permutations = 1000, pairwise=T, p.adjust.m="BH")
allarg.permdisp

#Permutation test for homogeneity of multivariate dispersions
#Permutation: free
#Number of permutations: 1000

#Response: Distances
#             Df   Sum Sq         Mean Sq     F       N.Perm  Pr(>F)  
#Groups       4     2648.9        662.22    2.4681     1000   0.05794 .
#Residuals    37    9927.5        268.31                    


#Subsetting for procrustes analysis

IMPUTED_ALLARGS_SUBSET<- subset_samples(IMPUTED_ALLARGS, CollectionPhase== "POST_SHOWER" | CollectionPhase== "SWINE" | CollectionPhase== "WORKDAY_END" | CollectionPhase== "WORKDAY_START")
ALLARGS_physeq.dist.SUB <- vegdist(decostand(t(otu_table(IMPUTED_ALLARGS_SUBSET)),"rclr"), method = "euclidean")

###################medimp_group_subset.ps
medimp_group_subset.ps

MEDIMPARGCounts <- as.data.frame(phyloseq::otu_table(medimp_group_subset.ps))
if(phyloseq::taxa_are_rows(medimp_group_subset.ps) == FALSE){ otus <- t(MEDIMPARGCounts) }
MEDIMPARGCounts

#zCompositions component of the anlaysis, using bayesian approach:


MEDIMPARG_cmultRepl<-cmultRepl(MEDIMPARGCounts,  label = 0, method = c("GBM"), output = c("p-counts"), frac = 0.2, threshold = 0.5, adjust = TRUE, t= NULL, s= NULL, z.warning = 0.8, suppress.print = FALSE)

MEDIMPARG_cmultRepl

#Now we must reconvert this OTU table to the phyloseq object

medimp_group_subset.ps <-phyloseq::otu_table(medimp_group_subset.ps) <- phyloseq::otu_table(MEDIMPARG_cmultRepl, taxa_are_rows = T)
medimp_group_subset.ps

#Create a new phyloseq object to include the imputed OTU table

IMPMED_RESISTOME = otu_table(medimp_group_subset.ps, taxa_are_rows = T)
GENES = tax_table(gen_mat)
ARG_SAMPLES = sample_data(resistomesamples_df)
IMPMED_RESISTOME_physeq <- phyloseq(IMPMED_RESISTOME, GENES, ARG_SAMPLES)
IMPMED_RESISTOME_physeq

#Now let's retry the ordination w/ this new phyloseq object

sum(taxa_sums(IMPMED_RESISTOME_physeq)==0) # 0 taxa not present across all samples
IMPUTED_IMPARGS<- prune_taxa(taxa_sums(IMPMED_RESISTOME_physeq) > 0, IMPMED_RESISTOME_physeq)#prunning low-abundance taxa of gene groups

IMPARGS_physeq.css <- phyloseq_transform_css(IMPUTED_IMPARGS)

IMPARGS_physeq.dist <- vegdist(decostand(t(otu_table(IMPUTED_IMPARGS)),"rclr"), method = "euclidean")
set.seed(1999)

IMPARGS_physeq.ord<- ordinate(IMPARGS_physeq.css, method = "PCoA", distance = "euclidean")

IMPARGS_plot<-plot_ordination(IMPARGS_physeq.css,IMPARGS_physeq.ord, type = "samples", color="CollectionPhase")
IMPARGS_plot

IMPARGcentroid<- IMPARGS_plot$data %>% group_by(CollectionPhase) %>% summarise(Axis.1=mean(Axis.1), Axis.2=mean(Axis.2), .groups = "drop")

#stressplot(plot_ord_all)
#str(plot_ord_all)

IMPARG_ordplot <- ggplot(IMPARGS_plot$data,IMPARGS_plot$mapping) + stat_ellipse(geom = "polygon", type= "t", level = 0.95, lty =1, lwd = 1.3, alpha= 0.07, show.legend =FALSE) + theme_bw() +geom_point(aes(fill=CollectionPhase), size=6, alpha=0.3, stroke=1.0,) + geom_point(data=IMPARGcentroid, size=13, shape=16, fill="black",mapping= aes(colour=CollectionPhase), stroke=4, alpha= 0.8, show.legend = FALSE) + scale_color_manual(name="Collection phase", breaks = c("WORKDAY_START", "WORKDAY_END", "POST_SHOWER", "SWINE", "ENVIRONMENT"), values = c("#789ED0", "#F79D40", "#72BF56", "#E4655A", "#AB8BCB"), labels= c("WORKDAY_START"= "Workday start", "WORKDAY_END"="Workday end", "POST_SHOWER"="Post-shower", "SWINE"="Swine", "ENVIRONMENT"="Environment")) +
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
        axis.title.x = element_text(size = 18, vjust = -1.5)) + coord_fixed(xlim=c(-25,25), ylim = c(-15,15))  
IMPARG_ordplot



IMP_IMPARG_metadata<- sample_data(IMPMED_RESISTOME_physeq)
pmeta<-as.data.frame(IMP_IMPARG_metadata)
pmeta

impmedarg.anosim<- anosim(IMPARGS_physeq.dist, pmeta$CollectionPhase, perm = 1000)
impmedarg.anosim
#ANOSIM statistic R: 0.2059 
#Significance: 0.000999 
impmadonis <- adonis2(IMPARGS_physeq.dist~pmeta$CollectionPhase, perm = 1000, p.adjust.m = "BH")
impmadonis
impmedarg.adonis <- pairwise.adonis(IMPARGS_physeq.dist, pmeta$CollectionPhase, perm = 1000, p.adjust.m = "BH")
impmedarg.adonis


#                           pairs Df SumsOfSqs   F.Model         R2     p.value  p.adjusted sig
#1  WORKDAY_START vs POST_SHOWER  1  65.03746 0.9114664 0.04819649 0.547452547 0.547452547    
#2        WORKDAY_START vs SWINE  1 381.09941 4.3827729 0.19581010 0.000999001 0.004995005   *
#3  WORKDAY_START vs ENVIRONMENT  1 223.93513 2.9206085 0.22604264 0.023976024 0.064935065    
#4  WORKDAY_START vs WORKDAY_END  1 148.16270 1.7711499 0.09435490 0.086913087 0.101010101    
#5          POST_SHOWER vs SWINE  1 379.14020 4.4260784 0.19736301 0.000999001 0.004995005   *
#6    POST_SHOWER vs ENVIRONMENT  1 198.61872 2.6715466 0.21083035 0.031968032 0.064935065    
#7    POST_SHOWER vs WORKDAY_END  1 192.31873 2.3372587 0.12086815 0.032967033 0.064935065    
#8          SWINE vs ENVIRONMENT  1 222.24029 2.1697956 0.17829351 0.090909091 0.101010101    
#9          SWINE vs WORKDAY_END  1 208.31494 2.1084362 0.11034059 0.050949051 0.072784358    
#10   ENVIRONMENT vs WORKDAY_END  1 272.93332 2.7876196 0.23648707 0.038961039 0.064935065  


impmedarg.disper <- betadisper(IMPARGS_physeq.dist, pmeta$CollectionPhase)
plot(impmedarg.disper)
impmedarg.disper

impmedarg.permdisp<- permutest(impmedarg.disper, permutations = 1000, pairwise=T, p.adjust.m="BH")
impmedarg.permdisp

#Response: Distances
#Df  Sum Sq Mean Sq      F N.Perm Pr(>F)
#Groups     4  18.139  4.5347 0.5581   1000 0.7323
#Residuals 36 292.537  8.1260     

plot_grid(ALLARG_ordplot, IMPARG_ordplot,labels=c('A. Total resistome', 'B. Medically important resistome fraction'), align= "vh", rel_widths = c(10, 10), jhust=0, axis='bt', label_x=-0.04, label_size=14)

#########################################################################
#MEDIMPARG ANALYSIS USING >99% GF cut-off
#########################################################################
## phyloseq object:
imparg99_mat<-read.csv(file.choose()) ##MedImp_work_resistome_coverage99_gene_count_matrix 
imparg99_tax_mat<- read.csv(file.choose())#MedImpAMR_taxa
imparg99_resistomesamples_df <- read.csv(file.choose())#FINAL_SHOTGUN_MCOHS_WORKER_METADATA



##need to have row.names
row.names(imparg99_mat) <- imparg99_mat$gene
imparg99_mat <- imparg99_mat %>% dplyr::select(-gene) #remove the column Genes since it is now used as a row name

row.names(imparg99_tax_mat) <- imparg99_tax_mat$gene
imparg99_tax_mat <- imparg99_tax_mat %>% dplyr::select(-gene) 

row.names(imparg99_resistomesamples_df) <- imparg99_resistomesamples_df$NovSeq_SS_Lib_ID
imparg99_resistomesamples_df <- imparg99_resistomesamples_df %>% dplyr::select(-NovSeq_SS_Lib_ID)
#Transform into matrixes otu and tax tables (sample table can be left as data frame)

imparg99_mat <- as.matrix(imparg99_mat)
imparg99_mat
imparg99_tax_mat <- as.matrix(imparg99_tax_mat)
imparg99_tax_mat

##Creating IMP_ARG99 phyloseq object:
imparg99_resistome = otu_table(imparg99_mat, taxa_are_rows = T)
imparg99_resistome

imparg99_tax = tax_table(imparg99_tax_mat)

imparg99_SAMPLES = sample_data(imparg99_resistomesamples_df)

MEDIMP99_physeq <- phyloseq(imparg99_resistome, imparg99_tax, imparg99_SAMPLES)
MEDIMP99_physeq
saveRDS(MEDIMP99_physeq, "IMPARG_99.phyloseq.rds")
MEDIMP99_physeq <- readRDS("C:/Users/is233/OneDrive/Desktop/AG_WORKER_METAGENOMICS_STUDY/MANUSCRIPT/ANALYSIS/WORKER_RESISTOME/prefiltered/IMPARG_99.phyloseq.rds")

imparg99_type.ps <-  tax_glom(MEDIMP99_physeq, "Type")
imparg99_class.ps <-   tax_glom(MEDIMP99_physeq, "Class")
#10 drug classes
imparg99_mechanism.ps <-  tax_glom(MEDIMP99_physeq, "Mechanism")
#15 medimp mechanisms
imparg99_group.ps <-   tax_glom(MEDIMP99_physeq, "Group")
#19 medimp groups
#487 medimp genes

##Order phyloseq samples by medimp99 resistome similarity



ps_seriate <- function(imparg99_group.ps,
                       method = "OLO_Ward",
                       dist = "bray",
                       tax_transform = "compositional",
                       add_variable = FALSE,
                       rank = Mechanism) {
  imparg99_group.ps <- ps_get(imparg99_group.ps)
  if (phyloseq::nsamples(imparg99_group.ps) <= 1) {
    return(imparg99_group.ps) # return early if no ordering
  }
  
  # aggregate taxa for ordering (possibly, as NA is no aggregation)
  psX <- tax_agg(imparg99_group.ps, rank = Group)
  # transform taxa for ordering (facilitated primarily for clr for PCA methods)
  ps_transformed <- tax_transform(psX, trans = tax_transform) %>% ps_get()
  
  if (method %in% seriation::list_seriation_methods(kind = "matrix")) {
    # directly seriate the otu matrix
    ser <- seriation::seriate(x = otu_get(ps_transformed), method = method)
  } else if (method %in% seriation::list_seriation_methods(kind = "dist")) {
    # calculate distance between samples
    distMat <- dist_get(dist_calc(data = ps_transformed, dist = dist))
    ser <- seriation::seriate(x = distMat, method = method)
  } else {
    stop(
      method, " is not a valid method in seriation::seriate!\n",
      "Nearest match is: ", findNearestSeriationMethods(method)[[1]]
    )
  }
  s_order <- seriation::get_order(ser)
  if (isTRUE(add_variable)) add_variable <- ".seriation_order" # default name
  if (!isFALSE(add_variable)) phyloseq::sample_data(imparg99_group.ps)[[add_variable]] <- s_order
  
  imparg99_group.ps <- ps_reorder(imparg99_group.ps, sample_order = s_order)
  return(imparg99_group.ps)
}

cols <- distinct_palette(n = 5, add = NA)
names(cols) <- unique(samdat_tbl(imparg99_group.ps)$CollectionPhase)

imparg99_group.ps %>%
  tax_transform("compositional",rank = "Group") %>%
  comp_heatmap(tax_anno = taxAnnotation(Prev. = anno_tax_prev(bar_width = 0.6, size = grid::unit(1, "cm"))), sample_ser_dist = "bray", sample_ser_counts ="bray",
               colors = heat_palette(palette = "Rocket", rev = TRUE),
               sample_anno = sampleAnnotation(
                 "Collection phase" = anno_sample("CollectionPhase"),
                 col = list(CollectionPhase = cols), border = TRUE))


imparg99_group.ps %>%
  tax_transform("compositional", rank = "Group") %>%
  comp_heatmap(
    tax_anno = taxAnnotation(Prev. = anno_tax_prev(bar_width = 0.6, size = grid::unit(1, "cm")), Abund.= anno_tax_box(bar_width = 0.6, size = grid::unit(1, "cm"))),
    colors = heat_palette(palette = "Reds", rev = TRUE),
    sample_anno = sampleAnnotation(
      "Collection phase" = anno_sample("CollectionPhase"),
      col = list(CollectionPhase = cols),
      border = TRUE
    )
  )

imparg99_group.ps %>%
  tax_transform("compositional", rank="Group")%>%
  comp_heatmap(colors = heat_palette(sym = FALSE), name = "rCLR",
    row_title = NULL,
    tax_anno = taxAnnotation(Prev. = anno_tax_prev(bar_width = 0.6, size = grid::unit(1, "cm"))),
    sample_anno = sampleAnnotation(
      "Collection phase" = anno_sample("CollectionPhase"),
      col = list(Collection_phase = cols), border = TRUE))

#Subsetting for procrustes analysis
IMPUTED_IMPARGS_SUB<- subset_samples(IMPUTED_ALLARGS, CollectionPhase== "POST_SHOWER" | CollectionPhase== "SWINE" | CollectionPhase== "WORKDAY_END" | CollectionPhase== "WORKDAY_START")
IMPARGS_physeq.dist.SUB <- vegdist(decostand(t(otu_table(IMPUTED_ALLARGS_SUBSET)),"rclr"), method = "euclidean")

############################################################################################
#DIFFERENTIAL ABUNDANCE TESTING#############################################################
############################################################################################

ALLARGCounts <- as.data.frame(phyloseq::otu_table(RESISTOME_physeq))
if(phyloseq::taxa_are_rows(RESISTOME_physeq) == FALSE){ otus <- t(ALLARGCounts) }
ALLARGCounts

#zCompositions component of the anlaysis, using bayesian approach:


ALLARG_cmultRepl<-cmultRepl(ALLARGCounts,  label = 0, method = c("GBM"), output = c("p-counts"), frac = 0.65, threshold = 0.5, adjust = TRUE, t= NULL, s= NULL, z.warning = 0.8, suppress.print = FALSE)

ALLARG_cmultRepl

#Now we must reconvert this OTU table to the phyloseq object

RESISTOME_physeq <-phyloseq::otu_table(RESISTOME_physeq) <- phyloseq::otu_table(ALLARG_cmultRepl, taxa_are_rows = T)


#Create a new phyloseq object to include the imputed 

IMPALL_RESISTOME = otu_table(RESISTOME_physeq, taxa_are_rows = T)
GENES = tax_table(gen_mat)
ARG_SAMPLES = sample_data(resistomesamples_df)
IMPUTED_ALLRESISTOME_physeq <- phyloseq(IMPALL_RESISTOME, GENES, ARG_SAMPLES)
IMPUTED_ALLRESISTOME_physeq

#All ARG group genes

sample_data(IMPUTED_ALLRESISTOME_physeq)$CollectionPhase <- as.factor(sample_data(IMPUTED_ALLRESISTOME_physeq)$CollectionPhase) # factorize for DESeq2
sample_data(IMPUTED_ALLRESISTOME_physeq)$EVSMOK <- as.factor(sample_data(IMPUTED_ALLRESISTOME_physeq)$EVSMOK) # factorize for DESeq2
sample_data(IMPUTED_ALLRESISTOME_physeq)$Shotgun_Post_Host_Rem  <- as.numeric(sample_data(IMPUTED_ALLRESISTOME_physeq)$Shotgun_Post_Host_Rem ) # factorize for DESeq2
sample_data(IMPUTED_ALLRESISTOME_physeq)$Shotgun_Post_Host_Rem  <- scale(sample_data(IMPUTED_ALLRESISTOME_physeq)$Shotgun_Post_Host_Rem , center = TRUE, scale=TRUE) # factorize for DESeq2
sample_data(IMPUTED_ALLRESISTOME_physeq)$Post_captureSampleConc_ng_ul    <- as.numeric(sample_data(IMPUTED_ALLRESISTOME_physeq)$Post_captureSampleConc_ng_ul) # factorize for DESeq2
sample_data(IMPUTED_ALLRESISTOME_physeq)$Post_captureSampleConc_ng_ul  <- scale(sample_data(IMPUTED_ALLRESISTOME_physeq)$Post_captureSampleConc_ng_ul , center = TRUE, scale=TRUE) # factorize for DESeq2
sample_data(IMPUTED_ALLRESISTOME_physeq)$GEND   <- as.factor(sample_data(IMPUTED_ALLRESISTOME_physeq)$GEND) # factorize for DESeq2
sample_data(IMPUTED_ALLRESISTOME_physeq)$AGE   <- scale(sample_data(IMPUTED_ALLRESISTOME_physeq)$AGE  , center = TRUE, scale=TRUE) # factorize for DESeq2
sample_data(IMPUTED_ALLRESISTOME_physeq)$BMI   <- scale(sample_data(IMPUTED_ALLRESISTOME_physeq)$BMI  , center = TRUE, scale=TRUE) # factorize for DESeq2

sample_data(IMPUTED_ALLRESISTOME_physeq)$HNDPRK   <- scale(sample_data(IMPUTED_ALLRESISTOME_physeq)$HNDPRK  , center = TRUE, scale=TRUE) # factorize for DESeq2
sample_data(IMPUTED_ALLRESISTOME_physeq)$EATPRK   <- scale(sample_data(IMPUTED_ALLRESISTOME_physeq)$EATPRK  , center = TRUE, scale=TRUE) # factorize for DESeq2

filtered_IMPUTED_ALLRESISTOME_physeq<-taxa_filter(IMPUTED_ALLRESISTOME_physeq, frequency = 0.8, treatment = 'CollectionPhase')
filtered_IMPUTED_ALLRESISTOME_physeq

IMPUTED_ALLRESISTOME_physeq <- tax_glom(filtered_IMPUTED_ALLRESISTOME_physeq, taxrank = 'Group', NArm = FALSE)
IMPUTED_ALLRESISTOME_physeq
BiocManager::install('EnhancedVolcano')
library('EnhancedVolcano')
##############SUBSET FOR WORKSTART_WORK_END

IMPUTED_ALLRESISTOME_physeq <- subset_samples(IMPUTED_ALLRESISTOME_physeq, CollectionPhase=="WORKDAY_START" | CollectionPhase=="WORKDAY_END")

#Filter sparse features of IMPUTED_ALLRESISTOME_physeq object so that we only consider those features with less than 90% zeros
#sample_data(IMPUTED_ALLRESISTOME_physeq)$CollectionPhase<- as.factor(sample_data(IMPUTED_ALLRESISTOME_physeq))
#IMPUTED_ALLRESISTOME_physeq
#nosparse_ALLRESISTOME<- prune_taxa(rowSums(otu_table(IMPUTED_ALLRESISTOME_physeq)==0) < ncol(otu_table(IMPUTED_ALLRESISTOME_physeq))*0.9, IMPUTED_ALLRESISTOME_physeq)

###Filter by prevalence: ARG groups w/ >10% sample prevalence with >2 hits produces a reduction from 685 unique ARG groups to 395 unique ARG groups

nosparse_ALLRESISTOME<-filter_taxa(IMPUTED_ALLRESISTOME_physeq, function(x){(sum(x > 10) > nsamples(IMPUTED_ALLRESISTOME_physeq)*0.1)}, prune = TRUE)
nosparse_ALLRESISTOME

DS_ALLRESISTOME= phyloseq_to_deseq2(nosparse_ALLRESISTOME,design= ~CollectionPhase+ Post_captureSampleConc_ng_ul + Shotgun_Post_Host_Rem + GEND + BMI + EVSMOK + EATPRK)
DS_ALLARGS<- estimateSizeFactors(DS_ALLRESISTOME, type="poscounts")
DS_ALLARGS$CollectionPhase <- relevel(DS_ALLARGS$CollectionPhase, ref = "WORKDAY_START")
DS_ALLARGS= DESeq(DS_ALLARGS, test= "Wald", fitType= "parametric") 
resultsNames(DS_ALLARGS)
alpha= 0.01
DS_ALLARGS_RESULT= results(DS_ALLARGS, alpha = alpha)
DS_ALLARGS <- lfcShrink(DS_ALLARGS, res= DS_ALLARGS_RESULT, type="ashr")
DS_ALLARGS_RESULT= DS_ALLARGS_RESULT[order(DS_ALLARGS_RESULT$padj, na.last = NA), ]
DS_ALLARGS_RESULT
summary(DS_ALLARGS_RESULT)
#out of 668 with nonzero total read count
#adjusted p-value < 0.01
#LFC > 0 (up)       : 46, 6.9%
#LFC < 0 (down)     : 17, 2.5%
#out of 533 with nonzero total read count with filtering
#adjusted p-value < 0.01
#LFC > 0 (up)       : 27, 5.1%
#LFC < 0 (down)     : 2, 0.38%
R<-as.data.frame(DS_ALLARGS_RESULT)
R_SIG_ALLARGS<- head(R[, 1:6], 29)
R_SIG_ALLARGS<- R_SIG_ALLARGS[order(R_SIG_ALLARGS$baseMean, na.last = NA), ]
write.csv(R_SIG_ALLARGS, "C:/Users/is233/OneDrive/Desktop/AG_WORKER_METAGENOMICS_STUDY/MANUSCRIPT/ANALYSIS/WORKER_RESISTOME/DESEQ_GROUP_WSWE.csv", row.names = TRUE)
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

sum( rowMeans( counts(DS_ALLARGS, normalized=TRUE)) > 5 )#check number of rows for nsub
vsdata <- vst(DS_ALLARGS, blind=TRUE, nsub=516)
vsdata_WSWE.df<-assay(vsdata)
write.csv(vsdata_WSWE.df, "C:/Users/is233/OneDrive/Desktop/AG_WORKER_METAGENOMICS_STUDY/MANUSCRIPT/ANALYSIS/WORKER_RESISTOME/DESEQ_GROUP_WSWE_vsdata.csv", row.names = TRUE)
DESEQ_PCA_GENE_WSWE<-plotPCA(vsdata, intgroup="CollectionPhase", ntop=1000) + theme_bw()
DESEQ_PCA_GENE_WSWE

ggplot(R_SIG_ALLARGS, aes(x=log2FoldChange, y = log10(baseMean))) + geom_point(aes(size=baseMean), shape=21, stroke=1, color="black", fill="gray")+ theme_bw()

DS_SIG_ALLARGS<- rownames(DS_ALLARGS_RESULT[1:29, ]) #Selecting the bottom 100 with the lowest adjusted p values
DS_SIG_ALLARGS
#total=median(sample_sums(IMPUTED_ALLRESISTOME_physeq))
#standf=function(x, t=total) round(t * (x / sum(x)))
DS_REL_ALLARGS<-phyloseq_transform_css(nosparse_ALLRESISTOME)
#DS_REL_ALLARGS<- transform_sample_counts(IMPUTED_ALLRESISTOME_physeq, function(x) x/sum(x)*10)
DS_SIG_REL_ALLARGS<- prune_taxa(DS_SIG_ALLARGS, DS_REL_ALLARGS)

DS_SIG_REL_ALLARGS


##############SUBSET FOR WORK_END_POST_SHOWER
IMPALL_RESISTOME = otu_table(RESISTOME_physeq, taxa_are_rows = T)
GENES = tax_table(gen_mat)
ARG_SAMPLES = sample_data(resistomesamples_df)
IMPUTED_ALLRESISTOME_physeq <- phyloseq(IMPALL_RESISTOME, GENES, ARG_SAMPLES)
IMPUTED_ALLRESISTOME_physeq

sample_data(IMPUTED_ALLRESISTOME_physeq)$CollectionPhase <- as.factor(sample_data(IMPUTED_ALLRESISTOME_physeq)$CollectionPhase) # factorize for DESeq2
sample_data(IMPUTED_ALLRESISTOME_physeq)$EVSMOK <- as.factor(sample_data(IMPUTED_ALLRESISTOME_physeq)$EVSMOK) # factorize for DESeq2
sample_data(IMPUTED_ALLRESISTOME_physeq)$Shotgun_Post_Host_Rem  <- as.numeric(sample_data(IMPUTED_ALLRESISTOME_physeq)$Shotgun_Post_Host_Rem ) # factorize for DESeq2
sample_data(IMPUTED_ALLRESISTOME_physeq)$Shotgun_Post_Host_Rem  <- scale(sample_data(IMPUTED_ALLRESISTOME_physeq)$Shotgun_Post_Host_Rem , center = TRUE, scale=TRUE) # factorize for DESeq2
sample_data(IMPUTED_ALLRESISTOME_physeq)$Post_captureSampleConc_ng_ul    <- as.numeric(sample_data(IMPUTED_ALLRESISTOME_physeq)$Post_captureSampleConc_ng_ul) # factorize for DESeq2
sample_data(IMPUTED_ALLRESISTOME_physeq)$Post_captureSampleConc_ng_ul  <- scale(sample_data(IMPUTED_ALLRESISTOME_physeq)$Post_captureSampleConc_ng_ul , center = TRUE, scale=TRUE) # factorize for DESeq2
sample_data(IMPUTED_ALLRESISTOME_physeq)$GEND   <- as.factor(sample_data(IMPUTED_ALLRESISTOME_physeq)$GEND) # factorize for DESeq2
sample_data(IMPUTED_ALLRESISTOME_physeq)$AGE   <- scale(sample_data(IMPUTED_ALLRESISTOME_physeq)$AGE  , center = TRUE, scale=TRUE) # factorize for DESeq2
sample_data(IMPUTED_ALLRESISTOME_physeq)$BMI   <- scale(sample_data(IMPUTED_ALLRESISTOME_physeq)$BMI  , center = TRUE, scale=TRUE) # factorize for DESeq2

sample_data(IMPUTED_ALLRESISTOME_physeq)$HNDPRK   <- scale(sample_data(IMPUTED_ALLRESISTOME_physeq)$HNDPRK  , center = TRUE, scale=TRUE) # factorize for DESeq2
sample_data(IMPUTED_ALLRESISTOME_physeq)$EATPRK   <- scale(sample_data(IMPUTED_ALLRESISTOME_physeq)$HNDPRK  , center = TRUE, scale=TRUE) # factorize for DESeq2

filtered_IMPUTED_ALLRESISTOME_physeq<-taxa_filter(IMPUTED_ALLRESISTOME_physeq, frequency = 0.8, treatment = 'CollectionPhase')
filtered_IMPUTED_ALLRESISTOME_physeq

IMPUTED_ALLRESISTOME_physeq <- tax_glom(filtered_IMPUTED_ALLRESISTOME_physeq, taxrank = 'Group', NArm = FALSE)
IMPUTED_ALLRESISTOME_physeq

IMPUTED_ALLRESISTOME_physeq <- subset_samples(IMPUTED_ALLRESISTOME_physeq, CollectionPhase=="POST_SHOWER" | CollectionPhase=="WORKDAY_END")

#Filter sparse features of IMPUTED_ALLRESISTOME_physeq object so that we only consider those features with less than 90% zeros
#sample_data(IMPUTED_ALLRESISTOME_physeq)$CollectionPhase<- as.factor(sample_data(IMPUTED_ALLRESISTOME_physeq))

#IMPUTED_ALLRESISTOME_physeq
#nosparse_ALLRESISTOME<- prune_taxa(rowSums(otu_table(IMPUTED_ALLRESISTOME_physeq)==0) < ncol(otu_table(IMPUTED_ALLRESISTOME_physeq))*0.9, IMPUTED_ALLRESISTOME_physeq)

###Filter by prevalence: ARG groups w/ >50% sample prevalence with >2 hits produces a reduction from 685 unique ARG groups to 395 unique ARG groups

nosparse_ALLRESISTOME<-filter_taxa(IMPUTED_ALLRESISTOME_physeq, function(x){(sum(x > 10) > nsamples(IMPUTED_ALLRESISTOME_physeq)*0.1)}, prune = TRUE)
nosparse_ALLRESISTOME

DS_ALLRESISTOME= phyloseq_to_deseq2(nosparse_ALLRESISTOME,design= ~CollectionPhase+ Post_captureSampleConc_ng_ul + Shotgun_Post_Host_Rem + GEND + BMI + EVSMOK + EATPRK)
DS_ALLARGS<- estimateSizeFactors(DS_ALLRESISTOME, type="poscounts")
DS_ALLARGS$CollectionPhase <- relevel(DS_ALLARGS$CollectionPhase, ref = "WORKDAY_END")
DS_ALLARGS= DESeq(DS_ALLARGS, test= "Wald", fitType= "parametric") 
resultsNames(DS_ALLARGS)
alpha= 0.01
DS_ALLARGS_RESULT= results(DS_ALLARGS, alpha = alpha)
DS_ALLARGS <- lfcShrink(DS_ALLARGS, res= DS_ALLARGS_RESULT, type="ashr")
DS_ALLARGS_RESULT= DS_ALLARGS_RESULT[order(DS_ALLARGS_RESULT$padj, na.last = NA), ]
DS_ALLARGS_RESULT
summary(DS_ALLARGS_RESULT)
#out of 629 with nonzero total read count
#adjusted p-value < 0.01
#LFC > 0 (up)       : 8, 1.3%
#LFC < 0 (down)     : 17, 2.7%
#outliers [1]       : 0, 0%
#low counts [2]     : 0, 0%
#out of 490 with nonzero total read count
#adjusted p-value < 0.01
#LFC > 0 (up)       : 4, 0.82%
#LFC < 0 (down)     : 14, 2.9%

R<-as.data.frame(DS_ALLARGS_RESULT)
R_SIG_ALLARGS<- head(R[, 1:6], 18)
R_SIG_ALLARGS<- R_SIG_ALLARGS[order(R_SIG_ALLARGS$baseMean, na.last = NA), ]
write.csv(R_SIG_ALLARGS, "C:/Users/is233/OneDrive/Desktop/AG_WORKER_METAGENOMICS_STUDY/MANUSCRIPT/ANALYSIS/WORKER_RESISTOME/DESEQ_GROUP_WEPS.csv", row.names = TRUE)
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

sum( rowMeans( counts(DS_ALLARGS, normalized=TRUE)) > 5 )#check number of rows for nsub
vsdata <- vst(DS_ALLARGS, blind=TRUE, nsub=467)
DESEQ_PCA_GENE_WEPS<-plotPCA(vsdata, intgroup="CollectionPhase", ntop=1000) + theme_bw()
DESEQ_PCA_GENE_WEPS


DS_SIG_ALLARGS<- rownames(DS_ALLARGS_RESULT[1:80, ]) #Selecting the bottom 100 with the lowest adjusted p values
#total=median(sample_sums(IMPUTED_ALLRESISTOME_physeq))
#standf=function(x, t=total) round(t * (x / sum(x)))
DS_REL_ALLARGS<-phyloseq_transform_css(IMPUTED_ALLRESISTOME_physeq)
#DS_REL_ALLARGS<- transform_sample_counts(IMPUTED_ALLRESISTOME_physeq, function(x) x/sum(x)*10)
DS_SIG_REL_ALLARGS<- prune_taxa(DS_SIG_ALLARGS, DS_REL_ALLARGS)

DS_SIG_REL_ALLARGS


##############SUBSET FOR WORK_START_POST_SHOWER
IMPALL_RESISTOME = otu_table(RESISTOME_physeq, taxa_are_rows = T)
GENES = tax_table(gen_mat)
ARG_SAMPLES = sample_data(resistomesamples_df)
IMPUTED_ALLRESISTOME_physeq <- phyloseq(IMPALL_RESISTOME, GENES, ARG_SAMPLES)

sample_data(IMPUTED_ALLRESISTOME_physeq)$CollectionPhase <- as.factor(sample_data(IMPUTED_ALLRESISTOME_physeq)$CollectionPhase) # factorize for DESeq2
sample_data(IMPUTED_ALLRESISTOME_physeq)$EVSMOK <- as.factor(sample_data(IMPUTED_ALLRESISTOME_physeq)$EVSMOK) # factorize for DESeq2
sample_data(IMPUTED_ALLRESISTOME_physeq)$Shotgun_Post_Host_Rem  <- as.numeric(sample_data(IMPUTED_ALLRESISTOME_physeq)$Shotgun_Post_Host_Rem ) # factorize for DESeq2
sample_data(IMPUTED_ALLRESISTOME_physeq)$Shotgun_Post_Host_Rem  <- scale(sample_data(IMPUTED_ALLRESISTOME_physeq)$Shotgun_Post_Host_Rem , center = TRUE, scale=TRUE) # factorize for DESeq2
sample_data(IMPUTED_ALLRESISTOME_physeq)$Post_captureSampleConc_ng_ul    <- as.numeric(sample_data(IMPUTED_ALLRESISTOME_physeq)$Post_captureSampleConc_ng_ul) # factorize for DESeq2
sample_data(IMPUTED_ALLRESISTOME_physeq)$Post_captureSampleConc_ng_ul  <- scale(sample_data(IMPUTED_ALLRESISTOME_physeq)$Post_captureSampleConc_ng_ul , center = TRUE, scale=TRUE) # factorize for DESeq2
sample_data(IMPUTED_ALLRESISTOME_physeq)$GEND   <- as.factor(sample_data(IMPUTED_ALLRESISTOME_physeq)$GEND) # factorize for DESeq2
sample_data(IMPUTED_ALLRESISTOME_physeq)$AGE   <- scale(sample_data(IMPUTED_ALLRESISTOME_physeq)$AGE  , center = TRUE, scale=TRUE) # factorize for DESeq2
sample_data(IMPUTED_ALLRESISTOME_physeq)$BMI   <- scale(sample_data(IMPUTED_ALLRESISTOME_physeq)$BMI  , center = TRUE, scale=TRUE) # factorize for DESeq2

#HNDPRK	EATPRK

sample_data(IMPUTED_ALLRESISTOME_physeq)$HNDPRK   <- scale(sample_data(IMPUTED_ALLRESISTOME_physeq)$HNDPRK  , center = TRUE, scale=TRUE) # factorize for DESeq2
sample_data(IMPUTED_ALLRESISTOME_physeq)$EATPRK   <- scale(sample_data(IMPUTED_ALLRESISTOME_physeq)$HNDPRK  , center = TRUE, scale=TRUE) # factorize for DESeq2

filtered_IMPUTED_ALLRESISTOME_physeq<-taxa_filter(IMPUTED_ALLRESISTOME_physeq, frequency = 0.8, treatment = 'CollectionPhase')
filtered_IMPUTED_ALLRESISTOME_physeq

IMPUTED_ALLRESISTOME_physeq <- tax_glom(filtered_IMPUTED_ALLRESISTOME_physeq, taxrank = 'Group', NArm = FALSE)
IMPUTED_ALLRESISTOME_physeq

IMPUTED_ALLRESISTOME_physeq <- subset_samples(IMPUTED_ALLRESISTOME_physeq, CollectionPhase=="POST_SHOWER" | CollectionPhase=="WORKDAY_START")

#Filter sparse features of IMPUTED_ALLRESISTOME_physeq object so that we only consider those features with less than 90% zeros
#sample_data(IMPUTED_ALLRESISTOME_physeq)$CollectionPhase<- as.factor(sample_data(IMPUTED_ALLRESISTOME_physeq))

#IMPUTED_ALLRESISTOME_physeq
#nosparse_ALLRESISTOME<- prune_taxa(rowSums(otu_table(IMPUTED_ALLRESISTOME_physeq)==0) < ncol(otu_table(IMPUTED_ALLRESISTOME_physeq))*0.9, IMPUTED_ALLRESISTOME_physeq)


nosparse_ALLRESISTOME<-filter_taxa(IMPUTED_ALLRESISTOME_physeq, function(x){(sum(x > 10) > nsamples(IMPUTED_ALLRESISTOME_physeq)*0.1)}, prune = TRUE)

DS_ALLRESISTOME= phyloseq_to_deseq2(nosparse_ALLRESISTOME,design= ~CollectionPhase+ Post_captureSampleConc_ng_ul + Shotgun_Post_Host_Rem + GEND + BMI + EVSMOK +EATPRK)
DS_ALLARGS<- estimateSizeFactors(DS_ALLRESISTOME, type="poscounts")
DS_ALLARGS$CollectionPhase <- relevel(DS_ALLARGS$CollectionPhase, ref = "WORKDAY_START")
DS_ALLARGS= DESeq(DS_ALLARGS, test= "Wald", fitType= "parametric") 
resultsNames(DS_ALLARGS)
alpha= 0.01
DS_ALLARGS_RESULT= results(DS_ALLARGS, alpha = alpha)
DS_ALLARGS <- lfcShrink(DS_ALLARGS, res= DS_ALLARGS_RESULT, type="ashr")
DS_ALLARGS_RESULT= DS_ALLARGS_RESULT[order(DS_ALLARGS_RESULT$padj, na.last = NA), ]
DS_ALLARGS_RESULT
summary(DS_ALLARGS_RESULT)
#out of 647 with nonzero total read count
#adjusted p-value < 0.01
#LFC > 0 (up)       : 31, 4.8%
#LFC < 0 (down)     : 21, 3.2%
#out of 488 with nonzero total read count
#adjusted p-value < 0.01
#LFC > 0 (up)       : 3, 0.61%
#LFC < 0 (down)     : 10, 2%
R<-as.data.frame(DS_ALLARGS_RESULT)
R_SIG_ALLARGS<- head(R[, 1:6], 13)
R_SIG_ALLARGS<- R_SIG_ALLARGS[order(R_SIG_ALLARGS$baseMean, na.last = NA), ]
write.csv(R_SIG_ALLARGS, "C:/Users/is233/OneDrive/Desktop/AG_WORKER_METAGENOMICS_STUDY/MANUSCRIPT/ANALYSIS/WORKER_RESISTOME/DESEQ_GROUP_WSPS.csv", row.names = TRUE)
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


sum( rowMeans( counts(DS_ALLARGS, normalized=TRUE)) > 5 )#check number of rows for nsub
vsdata <- vst(DS_ALLARGS, blind=TRUE, nsub=481)
DESEQ_PCA_GENE_PSWS<-plotPCA(vsdata, intgroup="CollectionPhase", ntop=1000) + theme_bw()
DESEQ_PCA_GENE_PSWS


DS_SIG_ALLARGS<- rownames(DS_ALLARGS_RESULT[1:50, ]) #Selecting the bottom 100 with the lowest adjusted p values

#total=median(sample_sums(IMPUTED_ALLRESISTOME_physeq))
#standf=function(x, t=total) round(t * (x / sum(x)))
DS_REL_ALLARGS<-phyloseq_transform_css(IMPUTED_ALLRESISTOME_physeq)
#DS_REL_ALLARGS<- transform_sample_counts(IMPUTED_ALLRESISTOME_physeq, function(x) x/sum(x)*10)
DS_SIG_REL_ALLARGS<- prune_taxa(DS_SIG_ALLARGS, DS_REL_ALLARGS)

DS_SIG_REL_ALLARGS

#All ARG group genes
sample_data(IMPUTED_ALLRESISTOME_physeq)$CollectionPhase <- as.factor(sample_data(IMPUTED_ALLRESISTOME_physeq)$CollectionPhase) # factorize for DESeq2
IMPUTED_ALLRESISTOME_physeq <- tax_glom(IMPUTED_ALLRESISTOME_physeq, taxrank = 'Group', NArm = FALSE)

##############SUBSET FOR SWINE_WORK_END
IMPALL_RESISTOME = otu_table(RESISTOME_physeq, taxa_are_rows = T)
GENES = tax_table(gen_mat)
ARG_SAMPLES = sample_data(resistomesamples_df)
IMPUTED_ALLRESISTOME_physeq <- phyloseq(IMPALL_RESISTOME, GENES, ARG_SAMPLES)
IMPUTED_ALLRESISTOME_physeq

sample_data(IMPUTED_ALLRESISTOME_physeq)$CollectionPhase <- as.factor(sample_data(IMPUTED_ALLRESISTOME_physeq)$CollectionPhase) # factorize for DESeq2
sample_data(IMPUTED_ALLRESISTOME_physeq)$EVSMOK <- as.factor(sample_data(IMPUTED_ALLRESISTOME_physeq)$EVSMOK) # factorize for DESeq2
sample_data(IMPUTED_ALLRESISTOME_physeq)$Shotgun_Post_Host_Rem  <- as.numeric(sample_data(IMPUTED_ALLRESISTOME_physeq)$Shotgun_Post_Host_Rem ) # factorize for DESeq2
sample_data(IMPUTED_ALLRESISTOME_physeq)$Shotgun_Post_Host_Rem  <- scale(sample_data(IMPUTED_ALLRESISTOME_physeq)$Shotgun_Post_Host_Rem , center = TRUE, scale=TRUE) # factorize for DESeq2
sample_data(IMPUTED_ALLRESISTOME_physeq)$Post_captureSampleConc_ng_ul    <- as.numeric(sample_data(IMPUTED_ALLRESISTOME_physeq)$Post_captureSampleConc_ng_ul) # factorize for DESeq2
sample_data(IMPUTED_ALLRESISTOME_physeq)$Post_captureSampleConc_ng_ul  <- scale(sample_data(IMPUTED_ALLRESISTOME_physeq)$Post_captureSampleConc_ng_ul , center = TRUE, scale=TRUE) # factorize for DESeq2
sample_data(IMPUTED_ALLRESISTOME_physeq)$GEND   <- as.factor(sample_data(IMPUTED_ALLRESISTOME_physeq)$GEND) # factorize for DESeq2
sample_data(IMPUTED_ALLRESISTOME_physeq)$AGE   <- scale(sample_data(IMPUTED_ALLRESISTOME_physeq)$AGE  , center = TRUE, scale=TRUE) # factorize for DESeq2
sample_data(IMPUTED_ALLRESISTOME_physeq)$BMI   <- scale(sample_data(IMPUTED_ALLRESISTOME_physeq)$BMI  , center = TRUE, scale=TRUE) # factorize for DESeq2

filtered_IMPUTED_ALLRESISTOME_physeq<-taxa_filter(IMPUTED_ALLRESISTOME_physeq, frequency = 0.8, treatment = 'CollectionPhase')
filtered_IMPUTED_ALLRESISTOME_physeq

IMPUTED_ALLRESISTOME_physeq <- tax_glom(filtered_IMPUTED_ALLRESISTOME_physeq, taxrank = 'Group', NArm = FALSE)
IMPUTED_ALLRESISTOME_physeq

IMPUTED_ALLRESISTOME_physeq <- subset_samples(IMPUTED_ALLRESISTOME_physeq, CollectionPhase=="SWINE" | CollectionPhase=="WORKDAY_END")

#Filter sparse features of IMPUTED_ALLRESISTOME_physeq object so that we only consider those features with less than 90% zeros
#sample_data(IMPUTED_ALLRESISTOME_physeq)$CollectionPhase<- as.factor(sample_data(IMPUTED_ALLRESISTOME_physeq))

#IMPUTED_ALLRESISTOME_physeq
#nosparse_ALLRESISTOME<- prune_taxa(rowSums(otu_table(IMPUTED_ALLRESISTOME_physeq)==0) < ncol(otu_table(IMPUTED_ALLRESISTOME_physeq))*0.9, IMPUTED_ALLRESISTOME_physeq)


nosparse_ALLRESISTOME<-filter_taxa(IMPUTED_ALLRESISTOME_physeq, function(x){(sum(x > 10) > nsamples(IMPUTED_ALLRESISTOME_physeq)*0.1)}, prune = TRUE)

DS_ALLRESISTOME= phyloseq_to_deseq2(nosparse_ALLRESISTOME,design= ~CollectionPhase+ Post_captureSampleConc_ng_ul + Shotgun_Post_Host_Rem)
DS_ALLARGS<- estimateSizeFactors(DS_ALLRESISTOME, type="poscounts")
DS_ALLARGS$CollectionPhase <- relevel(DS_ALLARGS$CollectionPhase, ref = "SWINE")
DS_ALLARGS= DESeq(DS_ALLARGS, test= "Wald", fitType= "parametric") 
resultsNames(DS_ALLARGS)
alpha= 0.01
DS_ALLARGS_RESULT= results(DS_ALLARGS, alpha = alpha)
DS_ALLARGS <- lfcShrink(DS_ALLARGS, res= DS_ALLARGS_RESULT, type="ashr")
DS_ALLARGS_RESULT= DS_ALLARGS_RESULT[order(DS_ALLARGS_RESULT$padj, na.last = NA), ]
DS_ALLARGS_RESULT
summary(DS_ALLARGS_RESULT)
#out of 298 with nonzero total read count
#adjusted p-value < 0.01
#LFC > 0 (up)       : 102, 34%
#LFC < 0 (down)     : 0, 0%
#out of 505 with nonzero total read count
#adjusted p-value < 0.01
#LFC > 0 (up)       : 10, 2%
#LFC < 0 (down)     : 2, 0.4%
R<-as.data.frame(DS_ALLARGS_RESULT)
R_SIG_ALLARGS<- head(R[, 1:6], 12)
R_SIG_ALLARGS<- R_SIG_ALLARGS[order(R_SIG_ALLARGS$baseMean, na.last = NA), ]
write.csv(R_SIG_ALLARGS, "C:/Users/is233/OneDrive/Desktop/AG_WORKER_METAGENOMICS_STUDY/MANUSCRIPT/ANALYSIS/WORKER_RESISTOME/DESEQ_GROUP_SWINE_WE.csv", row.names = TRUE)
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

sum( rowMeans( counts(DS_ALLARGS, normalized=TRUE)) > 5 )#check number of rows for nsub
vsdata <- vst(DS_ALLARGS, blind=TRUE, nsub=490)
DESEQ_PCA_GENE_SWINE_WE<-plotPCA(vsdata, intgroup="CollectionPhase", ntop=1000) + theme_bw()
DESEQ_PCA_GENE_SWINE_WE


DS_SIG_ALLARGS<- rownames(DS_ALLARGS_RESULT[1:100, ]) #Selecting the bottom 100 with the lowest adjusted p values
DS_SIG_ALLARGS
#total=median(sample_sums(IMPUTED_ALLRESISTOME_physeq))
#standf=function(x, t=total) round(t * (x / sum(x)))
DS_REL_ALLARGS<-phyloseq_transform_css(IMPUTED_ALLRESISTOME_physeq)
#DS_REL_ALLARGS<- transform_sample_counts(IMPUTED_ALLRESISTOME_physeq, function(x) x/sum(x)*10)
DS_SIG_REL_ALLARGS<- prune_taxa(DS_SIG_ALLARGS, IMPUTED_ALLRESISTOME_physeq)

DS_SIG_REL_ALLARGS


###############################
#DIFFERENTIAL ABUNDANCE MECHANISM LEVEL
###############################
IMPALL_RESISTOME = otu_table(RESISTOME_physeq, taxa_are_rows = T)
GENES = tax_table(gen_mat)
ARG_SAMPLES = sample_data(resistomesamples_df)
IMPUTED_ALLRESISTOME_physeq <- phyloseq(IMPALL_RESISTOME, GENES, ARG_SAMPLES)
IMPUTED_ALLRESISTOME_physeq

sample_data(IMPUTED_ALLRESISTOME_physeq)$CollectionPhase <- as.factor(sample_data(IMPUTED_ALLRESISTOME_physeq)$CollectionPhase) # factorize for DESeq2
sample_data(IMPUTED_ALLRESISTOME_physeq)$EVSMOK <- as.factor(sample_data(IMPUTED_ALLRESISTOME_physeq)$EVSMOK) # factorize for DESeq2
sample_data(IMPUTED_ALLRESISTOME_physeq)$Shotgun_Post_Host_Rem  <- as.numeric(sample_data(IMPUTED_ALLRESISTOME_physeq)$Shotgun_Post_Host_Rem ) # factorize for DESeq2
sample_data(IMPUTED_ALLRESISTOME_physeq)$Shotgun_Post_Host_Rem  <- scale(sample_data(IMPUTED_ALLRESISTOME_physeq)$Shotgun_Post_Host_Rem , center = TRUE, scale=TRUE) # factorize for DESeq2
sample_data(IMPUTED_ALLRESISTOME_physeq)$Post_captureSampleConc_ng_ul    <- as.numeric(sample_data(IMPUTED_ALLRESISTOME_physeq)$Post_captureSampleConc_ng_ul) # factorize for DESeq2
sample_data(IMPUTED_ALLRESISTOME_physeq)$Post_captureSampleConc_ng_ul  <- scale(sample_data(IMPUTED_ALLRESISTOME_physeq)$Post_captureSampleConc_ng_ul , center = TRUE, scale=TRUE) # factorize for DESeq2
sample_data(IMPUTED_ALLRESISTOME_physeq)$GEND   <- as.factor(sample_data(IMPUTED_ALLRESISTOME_physeq)$GEND) # factorize for DESeq2
sample_data(IMPUTED_ALLRESISTOME_physeq)$AGE   <- scale(sample_data(IMPUTED_ALLRESISTOME_physeq)$AGE  , center = TRUE, scale=TRUE) # factorize for DESeq2
sample_data(IMPUTED_ALLRESISTOME_physeq)$BMI   <- scale(sample_data(IMPUTED_ALLRESISTOME_physeq)$BMI  , center = TRUE, scale=TRUE) # factorize for DESeq2

sample_data(IMPUTED_ALLRESISTOME_physeq)$HNDPRK   <- scale(sample_data(IMPUTED_ALLRESISTOME_physeq)$HNDPRK  , center = TRUE, scale=TRUE) # factorize for DESeq2
sample_data(IMPUTED_ALLRESISTOME_physeq)$EATPRK   <- scale(sample_data(IMPUTED_ALLRESISTOME_physeq)$HNDPRK  , center = TRUE, scale=TRUE) # factorize for DESeq2

filtered_IMPUTED_ALLRESISTOME_physeq<-taxa_filter(IMPUTED_ALLRESISTOME_physeq, frequency = 0.8, treatment = 'CollectionPhase')
filtered_IMPUTED_ALLRESISTOME_physeq

IMPUTED_ALLRESISTOME_physeq <- tax_glom(IMPUTED_ALLRESISTOME_physeq, taxrank = 'Mechanism', NArm = FALSE)

##############SUBSET FOR WORKSTART_WORK_END

IMPUTED_ALLRESISTOME_physeq <- subset_samples(IMPUTED_ALLRESISTOME_physeq, CollectionPhase=="WORKDAY_START" | CollectionPhase=="WORKDAY_END")

#Filter sparse features of IMPUTED_ALLRESISTOME_physeq object so that we only consider those features with less than 90% zeros
#sample_data(IMPUTED_ALLRESISTOME_physeq)$CollectionPhase<- as.factor(sample_data(IMPUTED_ALLRESISTOME_physeq))
IMPUTED_ALLRESISTOME_physeq
#nosparse_ALLRESISTOME<- prune_taxa(rowSums(otu_table(IMPUTED_ALLRESISTOME_physeq)==0) < ncol(otu_table(IMPUTED_ALLRESISTOME_physeq))*0.9, IMPUTED_ALLRESISTOME_physeq)

nosparse_ALLRESISTOME<-filter_taxa(IMPUTED_ALLRESISTOME_physeq, function(x){(sum(x > 10) > nsamples(IMPUTED_ALLRESISTOME_physeq)*0.1)}, prune = TRUE)

DS_ALLRESISTOME= phyloseq_to_deseq2(IMPUTED_ALLRESISTOME_physeq,design= ~CollectionPhase+ Post_captureSampleConc_ng_ul + Shotgun_Post_Host_Rem + GEND + BMI + EVSMOK + EATPRK)
DS_ALLARGS<- estimateSizeFactors(DS_ALLRESISTOME, type="poscounts")
DS_ALLARGS$CollectionPhase <- relevel(DS_ALLARGS$CollectionPhase, ref = "WORKDAY_START")
DS_ALLARGS= DESeq(DS_ALLARGS, test= "Wald", fitType= "parametric") 
resultsNames(DS_ALLARGS)
alpha= 0.05
DS_ALLARGS_RESULT= results(DS_ALLARGS, alpha = alpha)
DS_ALLARGS <- lfcShrink(DS_ALLARGS, res= DS_ALLARGS_RESULT, type="ashr")
DS_ALLARGS_RESULT= DS_ALLARGS_RESULT[order(DS_ALLARGS_RESULT$padj, na.last = NA), ]
DS_ALLARGS_RESULT
summary(DS_ALLARGS_RESULT)
#out of 131 with nonzero total read count
#adjusted p-value < 0.05
#LFC > 0 (up)       : 12, 9.2%
#LFC < 0 (down)     : 3, 2.3%
R<-as.data.frame(DS_ALLARGS_RESULT)
R_SIG_ALLARGS<- head(R[, 1:6], 15)
R_SIG_ALLARGS<- R_SIG_ALLARGS[order(R_SIG_ALLARGS$baseMean, na.last = NA), ]
write.csv(R_SIG_ALLARGS, "C:/Users/is233/OneDrive/Desktop/AG_WORKER_METAGENOMICS_STUDY/MANUSCRIPT/ANALYSIS/WORKER_RESISTOME/DESEQ_MECHANISM_WSWE.csv", row.names = TRUE)
DS_ALLARGS <- estimateSizeFactors(DS_ALLARGS)
DS_ALLARGS <- estimateDispersions(DS_ALLARGS)
plotDispEsts(DS_ALLARGS)
par(mfrow=c(1,1))
# Make a basic volcano plot
with(R, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot", xlim=c(-50,50)))
with(subset(R, padj<.05 ), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
with(subset(R, padj<.05 & abs(log2FoldChange)>1.5), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))
EnhancedVolcano(R,
                lab = rownames(R),
                x = 'log2FoldChange',
                y = 'padj',  pCutoff = 50e-3,
                FCcutoff = 1.5,
                pointSize = log10(R$baseMean)+2,
                labSize = 5,
                shape=16,
                colAlpha = 0.4, legendLabels=c('Not sig.','Log (base 2) FC','FDRadj p <0.05',
                                               'FDRadj p <0.05 & Log (base 2) FC'),
                legendPosition = 'right',
                legendLabSize = 12,
                legendIconSize = 4.0)

##ggplot2 Volcanoplot

R$diffexpressed<- 'NO'
R$diffexpressed[R$log2FoldChange>1.5 & R$padj<0.05]<- 'Workday end_T2'
R$diffexpressed[R$log2FoldChange<(-1.5) & R$padj<0.05]<- 'Workday start_T1'

ggplot(data=R, aes(x=log2FoldChange, y=-log10(padj), col=diffexpressed))+ geom_vline(xintercept = c(-1.5, 1.5), col='gray', linetype='dashed') + geom_hline(yintercept= c(-log10(50e-3)), col='grey', linetype= 'dashed') + geom_point(size=log10(R$baseMean), shape=16, stroke=2)+
  
  theme_classic(base_size = 20)+ theme(
    axis.title.y = element_text(face = "bold", margin = margin(0,20,0,0), size = rel(1.1)), 
    axis.title.x = element_text(hjust = 0.5, face = "bold", margin = margin(20,0,0,0), size = rel(1.1)), 
    plot.title = element_text(hjust = 0.5))

sum( rowMeans( counts(DS_ALLARGS, normalized=TRUE)) > 5 )#check number of rows for nsub
vsdata <- vst(DS_ALLARGS, blind=TRUE, nsub=120)
DESEQ_PCA_MECHANISM_WSWE<-plotPCA(vsdata, intgroup="CollectionPhase", ntop=1000) + theme_bw()
DESEQ_PCA_MECHANISM_WSWE


DS_SIG_ALLARGS<- rownames(DS_ALLARGS_RESULT[1:258, ]) #Selecting the bottom 100 with the lowest adjusted p values
DS_SIG_ALLARGS
#total=median(sample_sums(IMPUTED_ALLRESISTOME_physeq))
#standf=function(x, t=total) round(t * (x / sum(x)))
DS_REL_ALLARGS<-phyloseq_transform_css(IMPUTED_ALLRESISTOME_physeq)
#DS_REL_ALLARGS<- transform_sample_counts(IMPUTED_ALLRESISTOME_physeq, function(x) x/sum(x)*10)
DS_SIG_REL_ALLARGS<- prune_taxa(DS_SIG_ALLARGS, DS_REL_ALLARGS)

DS_SIG_REL_ALLARGS


##############SUBSET FOR WORK_END_POST_SHOWER
IMPALL_RESISTOME = otu_table(RESISTOME_physeq, taxa_are_rows = T)
GENES = tax_table(gen_mat)
ARG_SAMPLES = sample_data(resistomesamples_df)
IMPUTED_ALLRESISTOME_physeq <- phyloseq(IMPALL_RESISTOME, GENES, ARG_SAMPLES)
IMPUTED_ALLRESISTOME_physeq

sample_data(IMPUTED_ALLRESISTOME_physeq)$CollectionPhase <- as.factor(sample_data(IMPUTED_ALLRESISTOME_physeq)$CollectionPhase) # factorize for DESeq2
sample_data(IMPUTED_ALLRESISTOME_physeq)$EVSMOK <- as.factor(sample_data(IMPUTED_ALLRESISTOME_physeq)$EVSMOK) # factorize for DESeq2
sample_data(IMPUTED_ALLRESISTOME_physeq)$Shotgun_Post_Host_Rem  <- as.numeric(sample_data(IMPUTED_ALLRESISTOME_physeq)$Shotgun_Post_Host_Rem ) # factorize for DESeq2
sample_data(IMPUTED_ALLRESISTOME_physeq)$Shotgun_Post_Host_Rem  <- scale(sample_data(IMPUTED_ALLRESISTOME_physeq)$Shotgun_Post_Host_Rem , center = TRUE, scale=TRUE) # factorize for DESeq2
sample_data(IMPUTED_ALLRESISTOME_physeq)$Post_captureSampleConc_ng_ul    <- as.numeric(sample_data(IMPUTED_ALLRESISTOME_physeq)$Post_captureSampleConc_ng_ul) # factorize for DESeq2
sample_data(IMPUTED_ALLRESISTOME_physeq)$Post_captureSampleConc_ng_ul  <- scale(sample_data(IMPUTED_ALLRESISTOME_physeq)$Post_captureSampleConc_ng_ul , center = TRUE, scale=TRUE) # factorize for DESeq2
sample_data(IMPUTED_ALLRESISTOME_physeq)$GEND   <- as.factor(sample_data(IMPUTED_ALLRESISTOME_physeq)$GEND) # factorize for DESeq2
sample_data(IMPUTED_ALLRESISTOME_physeq)$AGE   <- scale(sample_data(IMPUTED_ALLRESISTOME_physeq)$AGE  , center = TRUE, scale=TRUE) # factorize for DESeq2
sample_data(IMPUTED_ALLRESISTOME_physeq)$BMI   <- scale(sample_data(IMPUTED_ALLRESISTOME_physeq)$BMI  , center = TRUE, scale=TRUE) # factorize for DESeq2

sample_data(IMPUTED_ALLRESISTOME_physeq)$HNDPRK   <- scale(sample_data(IMPUTED_ALLRESISTOME_physeq)$HNDPRK  , center = TRUE, scale=TRUE) # factorize for DESeq2
sample_data(IMPUTED_ALLRESISTOME_physeq)$EATPRK   <- scale(sample_data(IMPUTED_ALLRESISTOME_physeq)$HNDPRK  , center = TRUE, scale=TRUE) # factorize for DESeq2


filtered_IMPUTED_ALLRESISTOME_physeq<-taxa_filter(IMPUTED_ALLRESISTOME_physeq, frequency = 0.8, treatment = 'CollectionPhase')
filtered_IMPUTED_ALLRESISTOME_physeq

IMPUTED_ALLRESISTOME_physeq <- tax_glom(filtered_IMPUTED_ALLRESISTOME_physeq, taxrank = 'Mechanism', NArm = FALSE)

IMPUTED_ALLRESISTOME_physeq <- subset_samples(IMPUTED_ALLRESISTOME_physeq, CollectionPhase=="POST_SHOWER" | CollectionPhase=="WORKDAY_END")

#Filter sparse features of IMPUTED_ALLRESISTOME_physeq object so that we only consider those features with less than 90% zeros
#sample_data(IMPUTED_ALLRESISTOME_physeq)$CollectionPhase<- as.factor(sample_data(IMPUTED_ALLRESISTOME_physeq))

IMPUTED_ALLRESISTOME_physeq
#nosparse_ALLRESISTOME<- prune_taxa(rowSums(otu_table(IMPUTED_ALLRESISTOME_physeq)==0) < ncol(otu_table(IMPUTED_ALLRESISTOME_physeq))*0.9, IMPUTED_ALLRESISTOME_physeq)

nosparse_ALLRESISTOME<-filter_taxa(IMPUTED_ALLRESISTOME_physeq, function(x){(sum(x > 10) > nsamples(IMPUTED_ALLRESISTOME_physeq)*0.1)}, prune = TRUE)

DS_ALLRESISTOME= phyloseq_to_deseq2(IMPUTED_ALLRESISTOME_physeq,design= ~CollectionPhase+ Post_captureSampleConc_ng_ul + Shotgun_Post_Host_Rem + GEND + BMI + EVSMOK + EATPRK)
DS_ALLARGS<- estimateSizeFactors(DS_ALLRESISTOME, type="poscounts")
DS_ALLARGS$CollectionPhase <- relevel(DS_ALLARGS$CollectionPhase, ref = "WORKDAY_END")
DS_ALLARGS= DESeq(DS_ALLARGS, test= "Wald", fitType= "parametric") 
resultsNames(DS_ALLARGS)
alpha= 0.05
DS_ALLARGS_RESULT= results(DS_ALLARGS, alpha = alpha)
DS_ALLARGS <- lfcShrink(DS_ALLARGS, res= DS_ALLARGS_RESULT, type="ashr")
DS_ALLARGS_RESULT= DS_ALLARGS_RESULT[order(DS_ALLARGS_RESULT$padj, na.last = NA), ]
DS_ALLARGS_RESULT
summary(DS_ALLARGS_RESULT)
#out of 126 with nonzero total read count
#adjusted p-value < 0.05
#LFC > 0 (up)       : 8, 6.3%
#LFC < 0 (down)     : 2, 1.6%
R<-as.data.frame(DS_ALLARGS_RESULT)
R_SIG_ALLARGS<- head(R[, 1:6], 10)
R_SIG_ALLARGS<- R_SIG_ALLARGS[order(R_SIG_ALLARGS$baseMean, na.last = NA), ]
write.csv(R_SIG_ALLARGS, "C:/Users/is233/OneDrive/Desktop/AG_WORKER_METAGENOMICS_STUDY/MANUSCRIPT/ANALYSIS/WORKER_RESISTOME/DESEQ_MECHANISM_WEPS.csv", row.names = TRUE)
DS_ALLARGS <- estimateSizeFactors(DS_ALLARGS)
DS_ALLARGS <- estimateDispersions(DS_ALLARGS)
plotDispEsts(DS_ALLARGS)
par(mfrow=c(1,1))
# Make a basic volcano plot
with(R, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot", xlim=c(-50,50)))
with(subset(R, padj<.05 ), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
with(subset(R, padj<.05 & abs(log2FoldChange)>1.5), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))
EnhancedVolcano(R,
                lab = rownames(R),
                x = 'log2FoldChange',
                y = 'padj',  pCutoff = 50e-3,
                FCcutoff = 1.5,
                pointSize = log10(R$baseMean)+2,
                labSize = 5,
                shape=16,
                colAlpha = 0.4, legendLabels=c('Not sig.','Log (base 2) FC','FDRadj p <0.05',
                                               'FDRadj p <0.05 & Log (base 2) FC'),
                legendPosition = 'right',
                legendLabSize = 12,
                legendIconSize = 4.0)

##ggplot2 Volcanoplot

R$diffexpressed<- 'NO'
R$diffexpressed[R$log2FoldChange>1.5 & R$padj<0.05]<- 'Post-shower_T3'
R$diffexpressed[R$log2FoldChange<(-1.5) & R$padj<0.05]<- 'Workday end_T2'

ggplot(data=R, aes(x=log2FoldChange, y=-log10(padj), col=diffexpressed))+ geom_vline(xintercept = c(-1.5, 1.5), col='gray', linetype='dashed') + geom_hline(yintercept= c(-log10(50e-3)), col='grey', linetype= 'dashed') + geom_point(size=log10(R$baseMean), shape=16, stroke=2)+
  
  theme_classic(base_size = 20)+ theme(
    axis.title.y = element_text(face = "bold", margin = margin(0,20,0,0), size = rel(1.1)), 
    axis.title.x = element_text(hjust = 0.5, face = "bold", margin = margin(20,0,0,0), size = rel(1.1)), 
    plot.title = element_text(hjust = 0.5))

sum( rowMeans( counts(DS_ALLARGS, normalized=TRUE)) > 5 )#check number of rows for nsub
vsdata <- vst(DS_ALLARGS, blind=TRUE, nsub=111)
DESEQ_PCA_MECHANISM_WEPS<-plotPCA(vsdata, intgroup="CollectionPhase", ntop=1000) + theme_bw()
DESEQ_PCA_MECHANISM_WEPS


ComplexHeatmap::pheatmap(R_SIG_ALLARGS, cluster_rows=FALSE, show_rownames=TRUE,
         cluster_cols=FALSE)

DS_SIG_ALLARGS<- rownames(DS_ALLARGS_RESULT[1:80, ]) #Selecting the bottom 100 with the lowest adjusted p values
#total=median(sample_sums(IMPUTED_ALLRESISTOME_physeq))
#standf=function(x, t=total) round(t * (x / sum(x)))
DS_REL_ALLARGS<-phyloseq_transform_css(IMPUTED_ALLRESISTOME_physeq)
#DS_REL_ALLARGS<- transform_sample_counts(IMPUTED_ALLRESISTOME_physeq, function(x) x/sum(x)*10)
DS_SIG_REL_ALLARGS<- prune_taxa(DS_SIG_ALLARGS, DS_REL_ALLARGS)

DS_SIG_REL_ALLARGS



##############SUBSET FOR WORK_START_POST_SHOWER
IMPALL_RESISTOME = otu_table(RESISTOME_physeq, taxa_are_rows = T)
GENES = tax_table(gen_mat)
ARG_SAMPLES = sample_data(resistomesamples_df)
IMPUTED_ALLRESISTOME_physeq <- phyloseq(IMPALL_RESISTOME, GENES, ARG_SAMPLES)
IMPUTED_ALLRESISTOME_physeq

sample_data(IMPUTED_ALLRESISTOME_physeq)$CollectionPhase <- as.factor(sample_data(IMPUTED_ALLRESISTOME_physeq)$CollectionPhase) # factorize for DESeq2
sample_data(IMPUTED_ALLRESISTOME_physeq)$EVSMOK <- as.factor(sample_data(IMPUTED_ALLRESISTOME_physeq)$EVSMOK) # factorize for DESeq2
sample_data(IMPUTED_ALLRESISTOME_physeq)$Shotgun_Post_Host_Rem  <- as.numeric(sample_data(IMPUTED_ALLRESISTOME_physeq)$Shotgun_Post_Host_Rem ) # factorize for DESeq2
sample_data(IMPUTED_ALLRESISTOME_physeq)$Shotgun_Post_Host_Rem  <- scale(sample_data(IMPUTED_ALLRESISTOME_physeq)$Shotgun_Post_Host_Rem , center = TRUE, scale=TRUE) # factorize for DESeq2
sample_data(IMPUTED_ALLRESISTOME_physeq)$Post_captureSampleConc_ng_ul    <- as.numeric(sample_data(IMPUTED_ALLRESISTOME_physeq)$Post_captureSampleConc_ng_ul) # factorize for DESeq2
sample_data(IMPUTED_ALLRESISTOME_physeq)$Post_captureSampleConc_ng_ul  <- scale(sample_data(IMPUTED_ALLRESISTOME_physeq)$Post_captureSampleConc_ng_ul , center = TRUE, scale=TRUE) # factorize for DESeq2
sample_data(IMPUTED_ALLRESISTOME_physeq)$GEND   <- as.factor(sample_data(IMPUTED_ALLRESISTOME_physeq)$GEND) # factorize for DESeq2
sample_data(IMPUTED_ALLRESISTOME_physeq)$AGE   <- scale(sample_data(IMPUTED_ALLRESISTOME_physeq)$AGE  , center = TRUE, scale=TRUE) # factorize for DESeq2
sample_data(IMPUTED_ALLRESISTOME_physeq)$BMI   <- scale(sample_data(IMPUTED_ALLRESISTOME_physeq)$BMI  , center = TRUE, scale=TRUE) # factorize for DESeq2

sample_data(IMPUTED_ALLRESISTOME_physeq)$HNDPRK   <- scale(sample_data(IMPUTED_ALLRESISTOME_physeq)$HNDPRK  , center = TRUE, scale=TRUE) # factorize for DESeq2
sample_data(IMPUTED_ALLRESISTOME_physeq)$EATPRK   <- scale(sample_data(IMPUTED_ALLRESISTOME_physeq)$HNDPRK  , center = TRUE, scale=TRUE) # factorize for DESeq2

filtered_IMPUTED_ALLRESISTOME_physeq<-taxa_filter(IMPUTED_ALLRESISTOME_physeq, frequency = 0.8, treatment = 'CollectionPhase')
filtered_IMPUTED_ALLRESISTOME_physeq

IMPUTED_ALLRESISTOME_physeq <- tax_glom(filtered_IMPUTED_ALLRESISTOME_physeq, taxrank = 'Mechanism', NArm = FALSE)

IMPUTED_ALLRESISTOME_physeq <- subset_samples(IMPUTED_ALLRESISTOME_physeq, CollectionPhase=="POST_SHOWER" | CollectionPhase=="WORKDAY_START")

#Filter sparse features of IMPUTED_ALLRESISTOME_physeq object so that we only consider those features with less than 90% zeros
#sample_data(IMPUTED_ALLRESISTOME_physeq)$CollectionPhase<- as.factor(sample_data(IMPUTED_ALLRESISTOME_physeq))

IMPUTED_ALLRESISTOME_physeq
#nosparse_ALLRESISTOME<- prune_taxa(rowSums(otu_table(IMPUTED_ALLRESISTOME_physeq)==0) < ncol(otu_table(IMPUTED_ALLRESISTOME_physeq))*0.9, IMPUTED_ALLRESISTOME_physeq)


nosparse_ALLRESISTOME<-filter_taxa(IMPUTED_ALLRESISTOME_physeq, function(x){(sum(x > 10) > nsamples(IMPUTED_ALLRESISTOME_physeq)*0.1)}, prune = TRUE)


DS_ALLRESISTOME= phyloseq_to_deseq2(IMPUTED_ALLRESISTOME_physeq,design= ~CollectionPhase+ Post_captureSampleConc_ng_ul + Shotgun_Post_Host_Rem + GEND + BMI + EVSMOK + EATPRK)
DS_ALLARGS<- estimateSizeFactors(DS_ALLRESISTOME, type="poscounts")
DS_ALLARGS$CollectionPhase <- relevel(DS_ALLARGS$CollectionPhase, ref = "WORKDAY_START")
DS_ALLARGS= DESeq(DS_ALLARGS, test= "Wald", fitType= "parametric") 
resultsNames(DS_ALLARGS)
alpha= 0.05
DS_ALLARGS_RESULT= results(DS_ALLARGS, alpha = alpha)
DS_ALLARGS <- lfcShrink(DS_ALLARGS, res= DS_ALLARGS_RESULT, type="ashr")
DS_ALLARGS_RESULT= DS_ALLARGS_RESULT[order(DS_ALLARGS_RESULT$padj, na.last = NA), ]
DS_ALLARGS_RESULT
summary(DS_ALLARGS_RESULT)
#out of 128 with nonzero total read count
#adjusted p-value < 0.05
#LFC > 0 (up)       : 6, 4.7%
#LFC < 0 (down)     : 7, 5.5%
R<-as.data.frame(DS_ALLARGS_RESULT)
R_SIG_ALLARGS<- head(R[, 1:6], 13)
R_SIG_ALLARGS<- R_SIG_ALLARGS[order(R_SIG_ALLARGS$baseMean, na.last = NA), ]
write.csv(R_SIG_ALLARGS, "C:/Users/is233/OneDrive/Desktop/AG_WORKER_METAGENOMICS_STUDY/MANUSCRIPT/ANALYSIS/WORKER_RESISTOME/DESEQ_MECHANISM_WSPS.csv", row.names = TRUE)
DS_ALLARGS <- estimateSizeFactors(DS_ALLARGS)
DS_ALLARGS <- estimateDispersions(DS_ALLARGS)
plotDispEsts(DS_ALLARGS)
par(mfrow=c(1,1))
# Make a basic volcano plot
with(R, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot", xlim=c(-50,50)))
with(subset(R, padj<.05 ), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
with(subset(R, padj<.05 & abs(log2FoldChange)>1.5), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))
EnhancedVolcano(R,
                lab = rownames(R),
                x = 'log2FoldChange',
                y = 'padj',  pCutoff = 50e-3,
                FCcutoff = 1.5,
                pointSize = log10(R$baseMean)+2,
                labSize = 5,
                shape=16,
                colAlpha = 0.4, legendLabels=c('Not sig.','Log (base 2) FC','FDRadj p <0.05',
                                               'FDRadj p <0.05 & Log (base 2) FC'),
                legendPosition = 'right',
                legendLabSize = 12,
                legendIconSize = 4.0)

##ggplot2 Volcanoplot

R$diffexpressed<- 'NO'
R$diffexpressed[R$log2FoldChange>1.5 & R$padj<0.05]<- 'Post-shower_T3'
R$diffexpressed[R$log2FoldChange<(-1.5) & R$padj<0.05]<- 'Workday start_T1'

ggplot(data=R, aes(x=log2FoldChange, y=-log10(padj), col=diffexpressed))+ geom_vline(xintercept = c(-1.5, 1.5), col='gray', linetype='dashed') + geom_hline(yintercept= c(-log10(50e-3)), col='grey', linetype= 'dashed') + geom_point(size=log10(R$baseMean), shape=16, stroke=2)+
  
  theme_classic(base_size = 20)+ theme(
    axis.title.y = element_text(face = "bold", margin = margin(0,20,0,0), size = rel(1.1)), 
    axis.title.x = element_text(hjust = 0.5, face = "bold", margin = margin(20,0,0,0), size = rel(1.1)), 
    plot.title = element_text(hjust = 0.5))

sum( rowMeans( counts(DS_ALLARGS, normalized=TRUE)) > 5 )#check number of rows for nsub
vsdata <- vst(DS_ALLARGS, blind=TRUE, nsub=119)
DESEQ_PCA_MECHANISM_PSWS<-plotPCA(vsdata, intgroup="CollectionPhase", ntop=1000) + theme_bw()
DESEQ_PCA_MECHANISM_PSWS

DS_SIG_ALLARGS<- rownames(DS_ALLARGS_RESULT[1:50, ]) #Selecting the bottom 100 with the lowest adjusted p values

#total=median(sample_sums(IMPUTED_ALLRESISTOME_physeq))
#standf=function(x, t=total) round(t * (x / sum(x)))
DS_REL_ALLARGS<-phyloseq_transform_css(IMPUTED_ALLRESISTOME_physeq)
#DS_REL_ALLARGS<- transform_sample_counts(IMPUTED_ALLRESISTOME_physeq, function(x) x/sum(x)*10)
DS_SIG_REL_ALLARGS<- prune_taxa(DS_SIG_ALLARGS, DS_REL_ALLARGS)

DS_SIG_REL_ALLARGS

##############SUBSET FOR SWINE_WORK_END
IMPALL_RESISTOME = otu_table(RESISTOME_physeq, taxa_are_rows = T)
GENES = tax_table(gen_mat)
ARG_SAMPLES = sample_data(resistomesamples_df)
IMPUTED_ALLRESISTOME_physeq <- phyloseq(IMPALL_RESISTOME, GENES, ARG_SAMPLES)
IMPUTED_ALLRESISTOME_physeq

sample_data(IMPUTED_ALLRESISTOME_physeq)$CollectionPhase <- as.factor(sample_data(IMPUTED_ALLRESISTOME_physeq)$CollectionPhase) # factorize for DESeq2
sample_data(IMPUTED_ALLRESISTOME_physeq)$EVSMOK <- as.factor(sample_data(IMPUTED_ALLRESISTOME_physeq)$EVSMOK) # factorize for DESeq2
sample_data(IMPUTED_ALLRESISTOME_physeq)$Shotgun_Post_Host_Rem  <- as.numeric(sample_data(IMPUTED_ALLRESISTOME_physeq)$Shotgun_Post_Host_Rem ) # factorize for DESeq2
sample_data(IMPUTED_ALLRESISTOME_physeq)$Shotgun_Post_Host_Rem  <- scale(sample_data(IMPUTED_ALLRESISTOME_physeq)$Shotgun_Post_Host_Rem , center = TRUE, scale=TRUE) # factorize for DESeq2
sample_data(IMPUTED_ALLRESISTOME_physeq)$Post_captureSampleConc_ng_ul    <- as.numeric(sample_data(IMPUTED_ALLRESISTOME_physeq)$Post_captureSampleConc_ng_ul) # factorize for DESeq2
sample_data(IMPUTED_ALLRESISTOME_physeq)$Post_captureSampleConc_ng_ul  <- scale(sample_data(IMPUTED_ALLRESISTOME_physeq)$Post_captureSampleConc_ng_ul , center = TRUE, scale=TRUE) # factorize for DESeq2
sample_data(IMPUTED_ALLRESISTOME_physeq)$GEND   <- as.factor(sample_data(IMPUTED_ALLRESISTOME_physeq)$GEND) # factorize for DESeq2
sample_data(IMPUTED_ALLRESISTOME_physeq)$AGE   <- scale(sample_data(IMPUTED_ALLRESISTOME_physeq)$AGE  , center = TRUE, scale=TRUE) # factorize for DESeq2
sample_data(IMPUTED_ALLRESISTOME_physeq)$BMI   <- scale(sample_data(IMPUTED_ALLRESISTOME_physeq)$BMI  , center = TRUE, scale=TRUE) # factorize for DESeq2

filtered_IMPUTED_ALLRESISTOME_physeq<-taxa_filter(IMPUTED_ALLRESISTOME_physeq, frequency = 0.8, treatment = 'CollectionPhase')
filtered_IMPUTED_ALLRESISTOME_physeq

IMPUTED_ALLRESISTOME_physeq <- tax_glom(filtered_IMPUTED_ALLRESISTOME_physeq, taxrank = 'Mechanism', NArm = FALSE)

IMPUTED_ALLRESISTOME_physeq <- subset_samples(IMPUTED_ALLRESISTOME_physeq, CollectionPhase=="SWINE" | CollectionPhase=="WORKDAY_END")

#Filter sparse features of IMPUTED_ALLRESISTOME_physeq object so that we only consider those features with less than 90% zeros
#sample_data(IMPUTED_ALLRESISTOME_physeq)$CollectionPhase<- as.factor(sample_data(IMPUTED_ALLRESISTOME_physeq))

IMPUTED_ALLRESISTOME_physeq
#nosparse_ALLRESISTOME<- prune_taxa(rowSums(otu_table(IMPUTED_ALLRESISTOME_physeq)==0) < ncol(otu_table(IMPUTED_ALLRESISTOME_physeq))*0.9, IMPUTED_ALLRESISTOME_physeq)

nosparse_ALLRESISTOME<-filter_taxa(IMPUTED_ALLRESISTOME_physeq, function(x){(sum(x > 10) > nsamples(IMPUTED_ALLRESISTOME_physeq)*0.1)}, prune = TRUE)

DS_ALLRESISTOME= phyloseq_to_deseq2(IMPUTED_ALLRESISTOME_physeq,design= ~CollectionPhase+ Post_captureSampleConc_ng_ul + Shotgun_Post_Host_Rem)
DS_ALLARGS<- estimateSizeFactors(DS_ALLRESISTOME, type="poscounts")
DS_ALLARGS$CollectionPhase <- relevel(DS_ALLARGS$CollectionPhase, ref = "SWINE")
DS_ALLARGS= DESeq(DS_ALLARGS, test= "Wald", fitType= "parametric") 
resultsNames(DS_ALLARGS)
alpha= 0.05
DS_ALLARGS_RESULT= results(DS_ALLARGS, alpha = alpha)
DS_ALLARGS <- lfcShrink(DS_ALLARGS, res= DS_ALLARGS_RESULT, type="ashr")
DS_ALLARGS_RESULT= DS_ALLARGS_RESULT[order(DS_ALLARGS_RESULT$padj, na.last = NA), ]
DS_ALLARGS_RESULT
summary(DS_ALLARGS_RESULT)
#out of 78 with nonzero total read count
#adjusted p-value < 0.05
#LFC > 0 (up)       : 29, 37%
#LFC < 0 (down)     : 2, 2.6%
R<-as.data.frame(DS_ALLARGS_RESULT)
R_SIG_ALLARGS<- head(R[, 1:6], 31)
R_SIG_ALLARGS<- R_SIG_ALLARGS[order(R_SIG_ALLARGS$baseMean, na.last = NA), ]
write.csv(R_SIG_ALLARGS, "C:/Users/is233/OneDrive/Desktop/AG_WORKER_METAGENOMICS_STUDY/MANUSCRIPT/ANALYSIS/WORKER_RESISTOME/DESEQ_MECHANISM_SWINE_WE.csv", row.names = TRUE)
DS_ALLARGS <- estimateSizeFactors(DS_ALLARGS)
DS_ALLARGS <- estimateDispersions(DS_ALLARGS)
plotDispEsts(DS_ALLARGS)
par(mfrow=c(1,1))
# Make a basic volcano plot
with(R, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot", xlim=c(-50,50)))
with(subset(R, padj<.05 ), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
with(subset(R, padj<.05 & abs(log2FoldChange)>1.5), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))
EnhancedVolcano(R,
                lab = rownames(R),
                x = 'log2FoldChange',
                y = 'padj',  pCutoff = 50e-3,
                FCcutoff = 1.5,
                pointSize = log10(R$baseMean)+2,
                labSize = 5,
                shape=16,
                colAlpha = 0.4, legendLabels=c('Not sig.','Log (base 2) FC','FDRadj p <0.05',
                                               'FDRadj p <0.05 & Log (base 2) FC'),
                legendPosition = 'right',
                legendLabSize = 12,
                legendIconSize = 4.0)

##ggplot2 Volcanoplot

R$diffexpressed<- 'NO'
R$diffexpressed[R$log2FoldChange>1.5 & R$padj<0.05]<- 'Workday end_T2'
R$diffexpressed[R$log2FoldChange<(-1.5) & R$padj<0.05]<- 'Swine'

ggplot(data=R, aes(x=log2FoldChange, y=-log10(padj), col=diffexpressed))+ geom_vline(xintercept = c(-1.5, 1.5), col='gray', linetype='dashed') + geom_hline(yintercept= c(-log10(50e-3)), col='grey', linetype= 'dashed') + geom_point(size=log10(R$baseMean), shape=16, stroke=2)+
  
  theme_classic(base_size = 20)+ theme(
    axis.title.y = element_text(face = "bold", margin = margin(0,20,0,0), size = rel(1.1)), 
    axis.title.x = element_text(hjust = 0.5, face = "bold", margin = margin(20,0,0,0), size = rel(1.1)), 
    plot.title = element_text(hjust = 0.5))


sum( rowMeans( counts(DS_ALLARGS, normalized=TRUE)) > 5 )#check number of rows for nsub
vsdata <- vst(DS_ALLARGS, blind=TRUE, nsub=117)
DESEQ_PCA_MECHANISM_SWINE_WE<-plotPCA(vsdata, intgroup="CollectionPhase", ntop=1000) + theme_bw()
DESEQ_PCA_MECHANISM_SWINE_WE

DS_SIG_ALLARGS<- rownames(DS_ALLARGS_RESULT[1:100, ]) #Selecting the bottom 100 with the lowest adjusted p values
DS_SIG_ALLARGS
#total=median(sample_sums(IMPUTED_ALLRESISTOME_physeq))
#standf=function(x, t=total) round(t * (x / sum(x)))
DS_REL_ALLARGS<-phyloseq_transform_css(IMPUTED_ALLRESISTOME_physeq)
#DS_REL_ALLARGS<- transform_sample_counts(IMPUTED_ALLRESISTOME_physeq, function(x) x/sum(x)*10)
DS_SIG_REL_ALLARGS<- prune_taxa(DS_SIG_ALLARGS, IMPUTED_ALLRESISTOME_physeq)

DS_SIG_REL_ALLARGS


