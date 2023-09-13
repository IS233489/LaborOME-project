
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

BiocManager::install("decontam")
library(decontam)



install.packages("dplyr")
library("dplyr")
install.packages("forcats")
library("forcats")
#Metadata file name:: MCOHS_WORKER_METADATA

metadata <- read.csv("./MCOHS_WORKER_METADATA.csv")   ##  Swine_worker_16s_metadata.csv
head(metadata)
tail(metadata)

MiSeqMetadata <- read.csv ("./Noyes_Project_032_MiSeq_Summary.csv") ##Metadata output from 16S (MiSeq) run
head(MiSeqMetadata)
tail(MiSeqMetadata)

updated_metadata <- left_join(metadata, MiSeqMetadata, by= "SampleName")
tail(updated_metadata)

write.csv(updated_metadata, "C:/Users/is233/OneDrive/Desktop/AG_WORKER_METAGENOMICS_STUDY/MANUSCRIPT/ANALYSIS/WORKER_MICROBIOME/UPDATED_MCOHS_WORKER_METADATA.csv", row.names = FALSE)
  
track_file <- read.csv ("./track_asvtab_edited68samples.csv") ## track file from dada2
head(track_file)

final_updated_metadata <- left_join(updated_metadata,track_file, by="SampleName")
tail(final_updated_metadata)

write.csv(final_updated_metadata, "C:/Users/is233/OneDrive/Desktop/AG_WORKER_METAGENOMICS_STUDY/MANUSCRIPT/ANALYSIS/WORKER_MICROBIOME/FINAL_MCOHS_WORKER_METADATA.csv", row.names = FALSE)

final_metadata <- read.csv("FINAL_MCOHS_WORKER_METADATA.csv")

## Analysis of 16S rawreads:
install.packages("ggplot2")

library("ggplot2")
install.packages("ggthemes")
library("ggthemes")
install.packages("lme4")
library("lme4")
install.packages("lmerTest")
library("lmerTest")
install.packages("lsmeans")
library("lsmeans")
install.packages("ggrepel")
library(ggrepel)
install.packages("pbkrtest")
library(pbkrtest)
library(devtools)
devtools::install_github("hadley/devtools")

pd = position_jitterdodge(dodge.width = 0.75, jitter.width = 0.3)

final_metadata$CollectionPhase <- factor (final_metadata$CollectionPhase, levels= c("WORKDAY_START",
                                                                                        "WORKDAY_END",
                                                                                        "POST_SHOWER",
                                                                                        'SWINE',
                                                                                        "ENVIRONMENT",
                                                                                        "MockComm",
                                                                                        "NegCtrl"))

final_metadata$OccupationalTask <- factor (final_metadata$OccupationalTask, levels= c("FARROWING",
                                                                                          "GESTATION",
                                                                                          "SOW_GESTATION",
                                                                                          'SWINE_PEN',
                                                                                          "ENV_FARROWING",
                                                                                          "ENV_GESTATION",
                                                                                          "MOVING_FOSTERING",
                                                                                          "VISITING VETERINARIAN",
                                                                                          "MockComm",
                                                                                          "NegCtrl"))
final_metadata$TASKEXPOSURE <- factor (final_metadata$TASKEXPOSURE, levels= c("DIRECT", "INDIRECT","SWINE",
                                                                          "ENVIRONMENT", "MockComm", "NegCtrl"))

final_metadata$qpcr_16s_copies_ul <- gsub(",", "", final_metadata$qpcr_16s_copies_ul) ##Remove the commas in qpcr copy column

final_metadata$qpcr_16s_copies_ul<- as.numeric(as.character(final_metadata$qpcr_16s_copies_ul))

final_metadata$DIRCOHRRATE<- as.numeric(as.character(final_metadata$DIRCOHRRATE))

final_metadata$DIRCOWKRATE <- as.numeric(as.character(final_metadata$DIRCOWKRATE))

###### ANALYSIS OF READ DEPTH AS A FUNCTION OF COLLECTION PHASE, OCCUPATIONAL TASK, AND TASK EXPOSURE (DIRCT VS. INDIRECT)  ######

raw_reads1 <- 
  final_metadata %>%
  filter(!CollectionPhase %in% c("MockComm","NegCtrl", "ENVIRONMENT")) %>%
  ggplot(aes(CollectionPhase, y= DADAReadPairInput, fill=CollectionPhase))+
ylim(0,NA)+
  geom_jitter(width = 0.25, shape=21, color= "black", size= 3) +
  geom_boxplot(alpha=0.3, shape=31, color= "black", outlier.colour = 'white', width=0.6, coef=5000, fatten=1, lwd=1)+
  #geom_point(aes(color = CollectionPhase), size=2, alpha=0.6, position = pd)+
  scale_fill_tableau(palette = "Classic 10 Medium", type = "regular",direction = 1)+  #Tableau 10
  scale_color_tableau(palette = "Classic 10 Medium", type = "regular",direction = 1)+  #Tableau 10
  theme_classic()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.line = element_line(size = 1))+
  theme(legend.position = "none") + ylab("Raw reads") + scale_x_discrete("Collection phase", labels= c("WORKDAY_START"= "Workday start", "WORKDAY_END"="Workday end", "POST_SHOWER"="Post-shower", "SWINE"="Swine", "ENVIRONMENT"="Environment", "MockComm"="Mock community", "NegCtrl"="Negative control"))
plot(raw_reads1)
topptx(filename = "rawreadbyCollectionPhase_plot.pptx", width=5, height=4) 
model_rawreads1<-
  lmer(DADAReadPairInput ~ CollectionPhase +(1|Worker), data = final_metadata, REML=F)
summary(model_rawreads1)
anova(model_rawreads1)
summary(anova)

#anova(model_rawreads1)
#Type III Analysis of Variance Table with Satterthwaite's method
 #               Sum Sq    Mean Sq NumDF DenDF F value   Pr(>F)   
#CollectionPhase 1.1737e+10 1956169348     6    47  4.2679 0.001642 **  

lsmeans(model_rawreads1, pairwise~CollectionPhase, adjust="tukey")
VarCorr(model_rawreads1)

#$contrasts
#contrast                    estimate    SE df t.ratio p.value
#WORKDAY_END - ENVIRONMENT     -63167 16583 47  -3.809  0.0069
#POST_SHOWER - ENVIRONMENT     -55451 16583 47  -3.344  0.0254

####################################Just trying something
censored_raw_reads1 <-
  final_metadata %>%
  filter(!CollectionPhase %in% c("MockComm","NegCtrl", "ENVIRONMENT")) %>%
  ggplot(aes(CollectionPhase, y= DADAReadPairInput, fill=CollectionPhase))+
  geom_boxplot(alpha=0.3, outlier.colour = 'white')+ylim(0,NA)+
  geom_point(aes(color = CollectionPhase), size=2, alpha=0.6, position = pd)+
  scale_fill_tableau(palette = "Classic 10 Medium", type = "regular",direction = 1)+  #Tableau 10
  scale_color_tableau(palette = "Classic 10 Medium", type = "regular",direction = 1)+  #Tableau 10
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  theme(legend.position = "none")
plot(censored_raw_reads1)
topptx(filename = "rawreadbyCollectionPhase_plot.pptx", width=5, height=4) 
model_censored_raw_reads1 <-subset (final_metadata, !CollectionPhase %in% c("MockComm", "NegCtrl"))
model_raw_reads1<-
  lmer((DADAReadPairInput) ~ CollectionPhase +(1|Worker), data = model_censored_raw_reads1, REML=F)
summary(model_raw_reads1)
anova(model_raw_reads1, type='III',test="F")

lsmeans(model_raw_reads1, pairwise~CollectionPhase, adjust="tukey")
VarCorr(model_raw_reads1)
########################################Just trying something

raw_reads2 <- 
  final_metadata %>%
  filter(!OccupationalTask %in% c("MockComm", "NegCtrl", "ENV_FARROWING", "ENV_GESTATION")) %>%
  ggplot(aes(OccupationalTask, y=DADAReadPairInput, fill=OccupationalTask))+
  geom_boxplot(alpha=0.45)+ylim(0,NA)+
  geom_point(aes(color = OccupationalTask),size=2, alpha=0.6, position = pd)+
  scale_fill_tableau(palette = "Superfishel Stone", type = "regular",direction = 1)+  
  scale_color_tableau(palette = "Superfishel Stone", type = "regular", direction = 1)+
  theme_bw()+
  theme(legend.position = "none")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
facet_wrap(~TASKEXPOSURE,scales='free_x', nrow=1)
#topptx(filename = "rawreadbyCollectionPhase_plot2.pptx", width=5, height=4) 

plot(raw_reads2)

model_raw_reads2 <-subset (final_metadata, !OccupationalTask %in% c("MockComm", "NegCtrl", "ENV_FARROWING", "ENV_GESTATION"))
model_rawreads2<-
  lmer(DADAReadPairInput ~ OccupationalTask +(1|Worker), data = model_raw_reads2, REML=F)
summary(model_rawreads2)
anova(model_rawreads2, type='III',test="F")
#> anova(model_rawreads2, type='III',test="F")
#Type III Analysis of Variance Table with Satterthwaite's method
                    # Sum Sq    Mean Sq NumDF DenDF F value  Pr(>F)  
#OccupationalTask 8136058169 1162294024     7    42  2.2252 0.05112 
###                  

raw_reads3 <- 
  final_metadata %>%
  filter(!TASKEXPOSURE %in% c("MockComm", "NegCtrl", "ENVIRONMENT")) %>%
  ggplot(aes(TASKEXPOSURE, y=DADAReadPairInput, fill=TASKEXPOSURE))+
  geom_boxplot(alpha=0.45)+ylim(0,NA)+
  geom_point(aes(color = TASKEXPOSURE),size=2, alpha=0.6, position = pd)+
  scale_fill_tableau(palette = "Superfishel Stone", type = "regular",direction = 1)+  
  scale_color_tableau(palette = "Superfishel Stone", type = "regular", direction = 1)+
  theme_bw()+
  theme(legend.position = "none")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
 facet_wrap(~CollectionPhase,scales='free_x', nrow=1)
#topptx(filename = "rawreadbyCollectionPhase_plot2.pptx", width=5, height=4) 

plot(raw_reads3)

model_raw_reads3 <-subset (final_metadata, !TASKEXPOSURE %in% c("MockComm", "NegCtrl", "ENVIRONMENT"))
model_rawreads3<-
  lmer(DADAReadPairInput ~ TASKEXPOSURE +(1|Worker), data = model_raw_reads3, REML=F)
summary(model_rawreads3)
anova(model_rawreads3, type='III',test="F")

#Type III Analysis of Variance Table with Satterthwaite's method
 #                Sum Sq    Mean Sq NumDF DenDF F value Pr(>F)
#TASKEXPOSURE 2173108218 1086554109     2    40  1.8371 0.1725

lsmeans(model_rawreads3, pairwise~TASKEXPOSURE, adjust="tukey")

VarCorr(model_rawreads3)

####### qpcr values ######
#Unfiltered description of collection phase on qPCR copy number
qpcr_copies1 <- 
  ggplot(final_metadata, aes(CollectionPhase, y= log10(qpcr_16s_copies_ul), fill=CollectionPhase))+
#  stat_summary(fun= median, width=0.6, size=0.6)+
#geom_dotplot(binaxis ="y", binwidth = 0.15, stackdir= "center")+
  geom_jitter(width = 0.25, shape=21, color= "black", size= 3) +
  geom_boxplot(alpha=0.3, shape=31, color= "black", outlier.colour = 'white', width=0.6, coef=5000, fatten=1, lwd=1)+ylim(0,NA)+
 #geom_point(aes(color = CollectionPhase), size=2, alpha=0.6, position = pd)+
  scale_fill_tableau(palette = "Classic 10 Medium", type = "regular",direction = 1)+  #Tableau 10
  scale_color_tableau(palette = "Classic 10 Medium", type = "regular",direction = 1)+  #Tableau 10
  theme_classic()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.line = element_line(size = 1))+
  theme(legend.position = "none")
plot(qpcr_copies1)
#topptx(filename = "rawreadbyCollectionPhase_plot.pptx", width=5, height=4) 
model_qpcr_copies1<-
  lmer(log(qpcr_16s_copies_ul) ~ CollectionPhase +(1|Worker), data = final_metadata, REML=F)
summary(model_qpcr_copies1)
anova(model_qpcr_copies1)
#> anova(model_qpcr_copies1)
#Type III Analysis of Variance Table with Satterthwaite's method
#               Sum Sq Mean Sq NumDF  DenDF F value    Pr(>F)    
#CollectionPhase 199.04  33.173     6 25.797  23.767 2.359e-09 ***
summary(anova)
lsmeans(model_qpcr_copies1, pairwise~CollectionPhase, adjust="tukey")
#$contrasts
#contrast                    estimate    SE   df t.ratio p.value
#ENVIRONMENT - MockComm        -8.279 2.121 27.6  -3.903  0.0089
#ENVIRONMENT - SWINE           -9.229 1.442 34.4  -6.401  <.0001
#ENVIRONMENT - WORKDAY_END     -5.971 1.442 34.4  -4.141  0.0036
#MockComm - NegCtrl            10.292 2.303 21.9   4.470  0.0031
#NegCtrl - SWINE              -11.242 1.697 21.4  -6.623  <.0001
#NegCtrl - WORKDAY_END         -7.983 1.697 21.4  -4.703  0.0019
#POST_SHOWER - SWINE           -6.422 0.832 34.4  -7.715  <.0001
#POST_SHOWER - WORKDAY_END     -3.164 0.528 25.6  -5.989  0.0001
#SWINE - WORKDAY_END            3.258 0.832 34.4   3.914  0.0068
#SWINE - WORKDAY_START          6.505 0.832 34.4   7.815  <.0001
#WORKDAY_END - WORKDAY_START    3.247 0.528 25.6   6.146  <.0001

VarCorr(model_qpcr_copies1)

#Filtered description of collection phase on qPCR copy number
qpcr_copies1_censored <-  
  final_metadata %>%
  filter(!CollectionPhase %in% c("MockComm", "NegCtrl", "ENVIRONMENT"))%>%
  mutate(CollectionPhase= fct_relevel(CollectionPhase,"WORKDAY_START", "WORKDAY_END", "POST_SHOWER", "SWINE", "ENVIRONMENT"))%>% #performed to change x axis order 
  ggplot(aes(CollectionPhase, y=log10(qpcr_16s_copies_ul), fill=CollectionPhase))+
  geom_jitter(shape=21, color= "black", size= 3, alpha=0.6, position = pd) +
  geom_boxplot(alpha=0.45, shape=31, color= "black", outlier.colour = 'white', width=0.6, coef=5000, fatten=1, lwd=1)+ylim(0,NA)+ xlab("Collection phase")+ ylab("Skin 16S qPCR gene copies/ ul (log10)")+
  #geom_point(aes(color = CollectionPhase), size=2, alpha=0.6, position = pd)+
  scale_fill_tableau(palette = "Classic 10 Medium", type = "regular",direction = 1)+  #Tableau 10
  scale_color_tableau(palette = "Classic 10 Medium", type = "regular",direction = 1)+  #Tableau 10
  theme_classic()+ theme(axis.text.x = element_text(size = 16),
                           axis.text.y = element_text(size = 16),
                           text=element_text(family="Helvetica Neue Medium",  size=18))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.line = element_line(size = 1))+
  theme(legend.position = "none") +  scale_x_discrete("Collection phase", labels= c("WORKDAY_START"= "Workday start", "WORKDAY_END"="Workday end", "POST_SHOWER"="Post-shower", "SWINE"="Swine"))
  
qpcr_copies1_censored + geom_text(data = tibble(x=1, y=0.1), size=3, aes(x=x, y=y, label= "Type III ANOVA P < 0.0001", fontface="italic"),inherit.aes=FALSE) + 
  geom_line(data=tibble(x=c(4,3), y=c(7.5,7.5)), aes(x=x, y=y), inherit.aes=FALSE)+ geom_line(data=tibble(x=c(4,2), y=c(8.4,8.4)), aes(x=x, y=y), inherit.aes=FALSE)+ 
  geom_line(data=tibble(x=c(4,1), y=c(9.4,9.4)), aes(x=x, y=y), inherit.aes=FALSE) + geom_line(data=tibble(x=c(2,1), y=c(6.6,6.6)), aes(x=x, y=y), inherit.aes=FALSE) + geom_line(data=tibble(x=c(2,3), y=c(7,7)), aes(x=x, y=y), inherit.aes=FALSE) + geom_text(data = tibble(x=3.5, y=7.8), size=5, aes(x=x, y=y, label= "***"),inherit.aes=FALSE)+
  geom_text(data = tibble(x=3, y=8.8), size=5, aes(x=x, y=y, label= "***"),inherit.aes=FALSE)+
  geom_text(data = tibble(x=2.5, y=9.8), size=5, aes(x=x, y=y, label= "***"),inherit.aes=FALSE)+geom_text(data = tibble(x=1.5, y=7.0), size=5, aes(x=x, y=y, label= "***"),inherit.aes=FALSE) + geom_text(data = tibble(x=2.5, y=7.4), size=5, aes(x=x, y=y, label= "***"),inherit.aes=FALSE)


#topptx(filename = "qPCRbyCollectionPhase_plot.pptx", width=3.55, height=4) 

model_qpcr_copies1_censored <-subset (final_metadata, !CollectionPhase %in% c("MockComm", "NegCtrl","ENVIRONMENT"))
model_logqpcrcensored<-
  lmer(log10(qpcr_16s_copies_ul) ~ CollectionPhase +(1|Worker), data = model_qpcr_copies1_censored, REML=F)
summary(model_logqpcrcensored)
anova(model_logqpcrcensored, type='III',test="F")
#> anova(model_logqpcrcensored, type='III',test="F")
#Type III Analysis of Variance Table with Satterthwaite's method
#                Sum Sq Mean Sq NumDF  DenDF F value    Pr(>F)    
#CollectionPhase 24.525   8.175     3 22.038  28.991 7.813e-08 ***

lsmeans(model_logqpcrcensored, pairwise~CollectionPhase, adjust="tukey")
#$contrasts
#contrast                    estimate    SE   df t.ratio p.value
#POST_SHOWER - SWINE          -2.7892 0.390 27.9  -7.147  <.0001
#POST_SHOWER - WORKDAY_END    -1.3741 0.237 21.8  -5.786  <.0001
#SWINE - WORKDAY_END           1.4151 0.390 27.9   3.626  0.0059
#SWINE - WORKDAY_START         2.8253 0.390 27.9   7.240  <.0001
#WORKDAY_END - WORKDAY_START   1.4102 0.237 21.8   5.938  <.0001
VarCorr(model_logqpcrcensored)



qpcr_copies2 <-  
  final_metadata %>%
  filter(!OccupationalTask %in% c("MockComm", "NegCtrl", "SWINE", "ENV_FARROWING", "ENV_GESTATION"))%>%
  ggplot(aes(OccupationalTask, y=log10(qpcr_16s_copies_ul), fill=OccupationalTask))+
  geom_jitter(width = 0.25, shape=21, color= "black", size= 3) +
  geom_boxplot(alpha=0.3, shape=31, color= "black", outlier.colour = 'white', width=0.6, coef=5000, fatten=1, lwd=1)+ylim(0,NA)+
  #geom_point(aes(color = CollectionPhase), size=2, alpha=0.6, position = pd)+
  scale_fill_tableau(palette = "Classic 10 Medium", type = "regular",direction = 1)+  #Tableau 10
  scale_color_tableau(palette = "Classic 10 Medium", type = "regular",direction = 1)+  #Tableau 10
  theme_classic()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.line = element_line(size = 1))+
  theme(legend.position = "none")+
  
  facet_wrap(~TASKEXPOSURE,scales='free_x', nrow=1)
#topptx(filename = "qpCRreadbyCollectionPhase_plot2.pptx", width=5, height=4) 
plot(qpcr_copies2)

model_qpcr_copies2_censored <-subset (final_metadata, !OccupationalTask %in% c("MockComm", "NegCtrl","ENV_FARROWING", "ENV_GESTATION"))
model_logqpcrcensored2<-
  lmer(log10(qpcr_16s_copies_ul) ~ OccupationalTask +(1|Worker), data = model_qpcr_copies2_censored, REML=F)
summary(model_logqpcrcensored2)
anova(model_logqpcrcensored2, type='III',test="F")
#> anova(model_logqpcrcensored, type='III',test="F")
#Type III Analysis of Variance Table with Satterthwaite's method
# Sum Sq Mean Sq NumDF  DenDF F value    Pr(>F)    
#OccupationalTask 60.059  8.5798     7 22.291  10.904 6.354e-06 ***

lsmeans(model_logqpcrcensored2, pairwise~OccupationalTask, adjust="tukey")
#$contrasts
#contrast                    estimate    SE   df t.ratio p.value
#FARROWING - GESTATION                       0.333 0.551 17.6   0.604  0.9894
#FARROWING - SOW_GESTATION                  -2.312 0.769 40.4  -3.008  0.0480
#FARROWING - SWINE_PEN                      -1.879 0.452 34.5  -4.152  0.0026
#FARROWING - MOVING_FOSTERING                1.350 0.729 17.6   1.851  0.4612
#FARROWING - VISITING VETERINARIAN           1.764 0.729 17.6   2.418  0.2027
#GESTATION - SOW_GESTATION                  -2.645 0.862 36.2  -3.069  0.0432
#GESTATION - SWINE_PEN                      -2.211 0.597 26.9  -3.702  0.0112
#GESTATION - MOVING_FOSTERING                1.017 0.827 17.6   1.230  0.8167
#GESTATION - VISITING VETERINARIAN           1.431 0.827 17.6   1.731  0.5310
#SOW_GESTATION - SWINE_PEN                   0.434 0.802 42.6   0.541  0.9941
#SOW_GESTATION - MOVING_FOSTERING            3.662 0.985 31.7   3.717  0.0092
#SOW_GESTATION - VISITING VETERINARIAN       4.076 0.985 31.7   4.137  0.0030
#SWINE_PEN - MOVING_FOSTERING                3.228 0.765 23.0   4.222  0.0039
#SWINE_PEN - VISITING VETERINARIAN           3.642 0.765 23.0   4.764  0.0011
#MOVING_FOSTERING - VISITING VETERINARIAN    0.414 0.955 17.6   0.433  0.9977

VarCorr(model_logqpcrcensored)

qpcr_copies3 <- 
  final_metadata %>%
  filter(!TASKEXPOSURE %in% c("MockComm", "NegCtrl", "ENVIRONMENT")) %>%
  ggplot(aes(TASKEXPOSURE, y=log10(qpcr_16s_copies_ul), fill=TASKEXPOSURE))+
  geom_jitter(width = 0.25, shape=21, color= "black", size= 3) +
  geom_boxplot(alpha=0.3, shape=31, color= "black", outlier.colour = 'white', width=0.6, coef=5000, fatten=1, lwd=1)+ylim(0,NA)+
  scale_fill_tableau(palette = "Classic 10 Medium", type = "regular",direction = 1)+  #Tableau 10
  scale_color_tableau(palette = "Classic 10 Medium", type = "regular",direction = 1)+  #Tableau 10 +
  theme_bw()+
  theme(legend.position = "none")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  facet_wrap(~CollectionPhase,scales='free_x', nrow=1)
#topptx(filename = "rawreadbyCollectionPhase_plot2.pptx", width=5, height=4) 

plot(qpcr_copies3)

model_qpcr_copies3 <-subset (final_metadata, !TASKEXPOSURE %in% c("MockComm", "NegCtrl", "ENVIRONMENT"))
model_logqpcrcopies3<-
  lmer(DADAReadPairInput ~ TASKEXPOSURE +(1|Worker), data = model_qpcr_copies3, REML=F)
summary(model_logqpcrcopies3)
anova(model_logqpcrcopies3, type='III',test="F")

#Type III Analysis of Variance Table with Satterthwaite's method
              #Sum Sq      Mean Sq    NumDF DenDF F value Pr(>F)  
#TASKEXPOSURE 6407811377 2135937126     3    42  3.7906 0.0171 *

lsmeans(model_logqpcrcopies3, pairwise~TASKEXPOSURE, adjust="tukey")
#$contrasts
#contrast   estimate    SE df t.ratio p.value
#DIRECT - ENVIRONMENT     -50339 17566 42  -2.866  0.0316
#INDIRECT - ENVIRONMENT   -53667 18557 42  -2.892  0.0296

VarCorr(model_logqpcrcopies3)


###### ANALYSIS OF READ QUALITY AS A FUNCTION OF COLLECTION PHASE, OCCUPATIONAL TASK, AND TASK EXPOSURE (DIRCT VS. INDIRECT)  ######


#Unfiltered description of collection phase on mean read quality
quality_copies1 <- 
  ggplot(final_metadata, aes(CollectionPhase, y= MiSeqMeanQuality , fill=CollectionPhase))+
  #  stat_summary(fun= median, width=0.6, size=0.6)+
  #geom_dotplot(binaxis ="y", binwidth = 0.15, stackdir= "center")+
  geom_jitter(width = 0.25, shape=21, color= "black", size= 3) +
  geom_boxplot(alpha=0.3, shape=31, color= "black", outlier.colour = 'white', width=0.6, coef=5000, fatten=1, lwd=1)+ylim(30.5,NA)+
  #geom_point(aes(color = CollectionPhase), size=2, alpha=0.6, position = pd)+
  scale_fill_tableau(palette = "Classic 10 Medium", type = "regular",direction = 1)+  #Tableau 10
  scale_color_tableau(palette = "Classic 10 Medium", type = "regular",direction = 1)+  #Tableau 10
  theme_classic()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.line = element_line(size = 1))+
  theme(legend.position = "none")
plot(quality_copies1)
#topptx(filename = "rawreadbyCollectionPhase_plot.pptx", width=5, height=4) 
model_quality_reads1<-
  lmer(MiSeqMeanQuality ~ CollectionPhase +(1|Worker), data = final_metadata, REML=F)
summary(model_quality_reads1)
anova(model_quality_reads1)
#> anova(model_qpcr_copies1)
#Type III Analysis of Variance Table with Satterthwaite's method
#               Sum Sq Mean Sq NumDF  DenDF F value    Pr(>F)    
#CollectionPhase 0.00064617 0.0001077     6 26.423  2.3831 0.05678 
summary(anova)
lsmeans(model_quality_reads1, pairwise~CollectionPhase, adjust="tukey")
#no significant difference by CollectionPhase

#Filtered description of collection phase by mean read quality
quality_copies1_censored <-  
  final_metadata %>%
  filter(!CollectionPhase %in% c("MockComm", "NegCtrl", "ENVIRONMENT"))%>%
  ggplot(aes(CollectionPhase, y=MiSeqMeanQuality, fill=CollectionPhase))+
  geom_jitter(width = 0.25, shape=21, color= "black", size= 3) +
  geom_boxplot(alpha=0.3, shape=31, color= "black", outlier.colour = 'white', width=0.6, coef=5000, fatten=1, lwd=1)+ylim(30.5,NA)+
  #geom_point(aes(color = CollectionPhase), size=2, alpha=0.6, position = pd)+
  scale_fill_tableau(palette = "Classic 10 Medium", type = "regular",direction = 1)+  #Tableau 10
  scale_color_tableau(palette = "Classic 10 Medium", type = "regular",direction = 1)+  #Tableau 10
  theme_classic()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.line = element_line(size = 1))+
  theme(legend.position = "none")+ scale_x_discrete("Collection phase", labels= c("WORKDAY_START"= "Workday start", "WORKDAY_END"="Workday end", "POST_SHOWER"="Post-shower", "SWINE"="Swine"))
#geom_boxplot(alpha=0.45, outlier.colour = "white")+ylim(0,8)+
#geom_jitter(aes(color = CollectionPhase),size=2, alpha=0.6, position = pd)+
#scale_fill_tableau(palette = "Superfishel Stone", type = "regular",direction = 1)+  #Tableau 10
#scale_color_tableau(palette = "Superfishel Stone", type = "regular", direction = 1)+
#scale_color_manual(values=c("grey35", "red",  "#59A14F", "blue",  "purple"))+
#scale_fill_manual(values=c("grey35", "red",  "#59A14F", "blue",  "purple"))+
#theme_classic()+
#theme(axis.text.x = element_text(size = 12),
#      axis.text.y = element_text(size = 12),
#      text=element_text(family="Times New Roman",  size=12))+
#theme(legend.position = "none")+
#theme(axis.text.x = element_text(angle = 45, hjust = 1))
#facet_wrap(~Exposure,scales='free_x', nrow=1)

plot(quality_copies1_censored)

#topptx(filename = "qPCRbyCollectionPhase_plot.pptx", width=3.55, height=4) 

model_quality_reads1_censored <-subset (final_metadata, !CollectionPhase %in% c("MockComm", "NegCtrl","ENVIRONMENT"))
model_qualityreadscensored<-
  lmer(MiSeqMeanQuality ~ CollectionPhase +(1|Worker), data = model_quality_reads1_censored, REML=F)
summary(model_qualityreadscensored)

anova(model_qualityreadscensored, type='III',test="F")

#> anova(model_logqpcrcensored, type='III',test="F")
#Type III Analysis of Variance Table with Satterthwaite's method
#                Sum Sq Mean Sq NumDF  DenDF F value    Pr(>F)    
#CollectionPhase 0.5099 0.16997     3 23.605  4.0407 0.01874 *

lsmeans(model_qualityreadscensored, pairwise~CollectionPhase, adjust="tukey")
#$contrasts
#contrast                    estimate    SE   df t.ratio p.value
#POST_SHOWER - WORKDAY_END      0.300 0.0967 25.0   3.103  0.0228
VarCorr(model_logqpcrcensored)



quality_copies2 <-  
  final_metadata %>%
  filter(!OccupationalTask %in% c("MockComm", "NegCtrl", "SWINE", "ENV_FARROWING", "ENV_GESTATION"))%>%
  ggplot(aes(OccupationalTask, y=MiSeqMeanQuality, fill=OccupationalTask))+
  geom_jitter(width = 0.25, shape=21, color= "black", size= 3) +
  geom_boxplot(alpha=0.3, shape=31, color= "black", outlier.colour = 'white', width=0.6, coef=5000, fatten=1, lwd=1)+ylim(30.5,NA)+
  #geom_point(aes(color = CollectionPhase), size=2, alpha=0.6, position = pd)+
  scale_fill_tableau(palette = "Classic 10 Medium", type = "regular",direction = 1)+  #Tableau 10
  scale_color_tableau(palette = "Classic 10 Medium", type = "regular",direction = 1)+  #Tableau 10
  theme_classic()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.line = element_line(size = 1))+
  theme(legend.position = "none")+
  
  facet_wrap(~TASKEXPOSURE,scales='free_x', nrow=1)
#topptx(filename = "qpCRreadbyCollectionPhase_plot2.pptx", width=5, height=4) 
plot(quality_copies2)

model_quality_copies2_censored <-subset (final_metadata, !OccupationalTask %in% c("MockComm", "NegCtrl","ENV_FARROWING", "ENV_GESTATION"))
model_qualitycensored2<-
  lmer(MiSeqMeanQuality ~ OccupationalTask +(1|Worker), data = model_quality_copies2_censored, REML=F)
summary(model_qualitycensored2)
anova(model_qualitycensored2, type='III',test="F")
#> anova(model_logqpcrcensored, type='III',test="F")
#Type III Analysis of Variance Table with Satterthwaite's method
# Sum Sq Mean Sq NumDF  DenDF F value    Pr(>F)    
#OccupationalTask 0.17352 0.024788     7    42  0.4294 0.8782

quality_copies3 <- 
  final_metadata %>%
  filter(!TASKEXPOSURE %in% c("MockComm", "NegCtrl", "ENVIRONMENT")) %>%
  ggplot(aes(TASKEXPOSURE, y=MiSeqMeanQuality, fill=TASKEXPOSURE))+
  geom_jitter(width = 0.25, shape=21, color= "black", size= 3) +
  geom_boxplot(alpha=0.3, shape=31, color= "black", outlier.colour = 'white', width=0.6, coef=5000, fatten=1, lwd=1)+ylim(31,NA)+
  scale_fill_tableau(palette = "Classic 10 Medium", type = "regular",direction = 1)+  #Tableau 10
  scale_color_tableau(palette = "Classic 10 Medium", type = "regular",direction = 1)+  #Tableau 10 +
  theme_bw()+
  theme(legend.position = "none")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  facet_wrap(~CollectionPhase,scales='free_x', nrow=1)
#topptx(filename = "rawreadbyCollectionPhase_plot2.pptx", width=5, height=4) 

plot(quality_copies3)

model_quality_copies3 <-subset (final_metadata, !TASKEXPOSURE %in% c("MockComm", "NegCtrl", "ENVIRONMENT"))
model_qualitycopies3<-
  lmer(DADAReadPairInput ~ TASKEXPOSURE +(1|Worker), data = model_quality_copies3, REML=F)
summary(model_qualitycopies3)
anova(model_qualitycopies3, type='III',test="F")

#Type III Analysis of Variance Table with Satterthwaite's method
#Sum Sq      Mean Sq    NumDF DenDF F value Pr(>F)  
#TASKEXPOSURE 6407811377 2135937126     3    42  3.7906 0.0171 *

lsmeans(model_qualitycopies3, pairwise~TASKEXPOSURE, adjust="tukey")
#$contrasts
#contrast   estimate    SE df t.ratio p.value
#DIRECT - ENVIRONMENT     -50339 18468 43.9  -2.726  0.0438
#ENVIRONMENT - INDIRECT    53667 19509 40.0   2.751  0.0424

VarCorr(model_logqpcrcopies3)


###### REGRESSION ANALYSIS OF 16s copy by swine exposure ######

#visitingvet<- final_metadata%>% filter(!OccupationalTask %in% c("FARROWING", #"MOVING_FOSTERING", "GESTATION", "SWINE_PEN", "SOW_GESTATION", #"ENV_FARROWING", "MockComm", "NegCtrl", "SWINE", "ENVIRONMENT"))
swinecontactreg <-  
  final_metadata %>%
  filter(!CollectionPhase %in% c("MockComm", "NegCtrl", "SWINE", "ENVIRONMENT") & !TASKEXPOSURE %in% c("INDIRECT"))%>%
  mutate(CollectionPhase= fct_relevel(CollectionPhase,"WORKDAY_START", "WORKDAY_END", "POST_SHOWER", "SWINE"))%>% #performed to change x axis order 
  ggplot(aes(x= log10(DIRCOHRRATE), y= log10(qpcr_16s_copies_ul), color= CollectionPhase)) +
  coord_cartesian(xlim = c(0.5, 3), ylim = c(0, 6)) + geom_point(size=4, alpha=0.6) +
  geom_label(label="Veterinarian",
             x=750, y= 15000, label.size = 0, color= "black")+
  geom_smooth(formula = y~x, method = "glm", se=T, level=0.9)+

  #geom_label_repel(min.segment.length = 0, max.overlaps = Inf, label.size =0, label.padding = 0.1, size=3
  
  theme(legend.position = "none")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.line =element_line(size = 1))+
labs(x= "Estimated swine exposure per hour (log10)", y= "Skin 16S qPCR gene copies/ul (log10)", color=NULL,
    axis.tickks= element_blank(), axis.line =element_line(size = 0.5))+
  scale_fill_tableau(palette = "Classic 10 Medium", type = "regular",direction = 1)+  #Tableau 10
  scale_color_tableau(palette = "Classic 10 Medium", type = "regular",direction = 1)+  #Tableau 10
  theme_bw()+ theme(axis.text.x = element_text(size = 16),
                         axis.text.y = element_text(size = 16),
                         text=element_text(family="Helvetica Neue Medium",  size=18))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.line = element_line(size = 0.5))+
  theme(legend.position = "none")
swinecontactreg
# New facet label names for supp variable
supp.labs <- c("WORKDAY_START"= "Workday start", "WORKDAY_END"="Workday end", "POST_SHOWER"="Post-shower")
swinecontactreg
swinecontactreg + facet_wrap(~CollectionPhase, labeller = labeller(CollectionPhase = supp.labs)) + theme(legend.position="none")


 

no_vet<- filter(final_metadata, !CollectionPhase %in% c("MockComm", "NegCtrl", "SWINE", "ENVIRONMENT", "POST_SHOWER", "WORKDAY_END") & !TASKEXPOSURE %in% c("INDIRECT"))


fit<-glm(log10(qpcr_16s_copies_ul)~log10(DIRCOHRRATE), data = no_vet)
summary(fit) 
anova(fit, type='III',test="F")


#Call:
#  glm(formula = log10(qpcr_16s_copies_ul) ~ log10(DIRCOHRRATE), 
#      data = no_vet)

#Deviance Residuals: 
#  1        2        3        4        5        6        7  
#1.4789   1.4188   0.3258  -1.2575  -1.6624   0.7295  -1.0332  

#Coefficients:
#  Estimate Std. Error t value Pr(>|t|)  
#(Intercept)           5.517      2.086   2.645   0.0457 *
#  log10(DIRCOHRRATE)   -1.141      1.042  -1.095   0.3234  
#---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#(Dispersion parameter for gaussian family taken to be 2.050183)

#Null deviance: 12.709  on 6  degrees of freedom
#Residual deviance: 10.251  on 5  degrees of freedom
#AIC: 28.535

no_vet<- filter(final_metadata, !CollectionPhase %in% c("MockComm", "NegCtrl", "SWINE", "ENVIRONMENT", "POST_SHOWER", "WORKDAY_START") & !TASKEXPOSURE %in% c("INDIRECT"))

fit<-glm(log10(qpcr_16s_copies_ul)~log10(DIRCOHRRATE), data = no_vet)
summary(fit) 
anova(fit, type='III',test="F")


#Call:
#  glm(formula = log10(qpcr_16s_copies_ul) ~ log10(DIRCOHRRATE), 
      #data = no_vet)

#Deviance Residuals: 
#  1         2         3         4         5         6         7  
#-0.02026   0.91609  -0.38045  -0.44209  -0.39109   0.10627   0.21152  

#Coefficients:
#  Estimate Std. Error t value Pr(>|t|)    
#(Intercept)          6.3890     0.7679   8.321  0.00041 ***
#  log10(DIRCOHRRATE)  -0.9937     0.3835  -2.592  0.04875 *  
#  ---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#(Dispersion parameter for gaussian family taken to be 0.2777606)

#Null deviance: 3.2543  on 6  degrees of freedom
#Residual deviance: 1.3888  on 5  degrees of freedom
#AIC: 14.543



lsmeans(fit, pairwise~CollectionPhase, adjust="tukey")

shapiro.test(no_vet$qpcr_16s_copies_ul)
shapiro.test(no_vet$DIRCOSWINE)

spearman <- cor.test(no_vet$qpcr_16s_copies_ul, no_vet$DIRCOHRRATE, method = "spearman")

spearman

#Estimate Std. Error t value Pr(>|t|)    
#(Intercept)  4.0217324  0.3730692  10.780 1.58e-07 ***
#  DIRCOWKRATE -0.0011035  0.0003719  -2.967   0.0118 *  
#  ---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#Residual standard error: 0.9569 on 12 degrees of freedom
#Multiple R-squared:  0.4231,	Adjusted R-squared:  0.3751 
#F-statistic: 8.802 on 1 and 12 DF,  p-value: 0.01177




###### MICROBIOME ANALYSIS ######

## phyloseq object:
asv_mat<-read.csv(file.choose(),2)      ## 
tax_mat<- read.csv(file.choose(),3)
samples_df <- read.xlsx(file.choose(),4)



##need to have row.names
row.names(asv_mat) <- asv_mat$ASV
asv_mat <- asv_mat %>% select(-ASV) #remove the column asv since it is now used as a row name

row.names(tax_mat) <- tax_mat$ASV
tax_mat <- tax_mat %>% select (-ASV) 

row.names(samples_df) <- samples_df$SampleName
samples_df <- samples_df %>% select (-SampleName)

#Transform into matrixes otu and tax tables (sample table can be left as data frame)

asv_mat <- as.matrix(asv_mat)
tax_mat <- as.matrix(tax_mat)

##phyloseq:
ASV = otu_table(asv_mat, taxa_are_rows = F)
TAX = tax_table(tax_mat)
samples = sample_data(samples_df)
physeq <- phyloseq(ASV, TAX, samples)
saveRDS(physeq, "ilya_phyloseq_68samples_final_7.20")

physeq_microbiome <- readRDS("ilya_phyloseq_68samples_final_7.20")

###Exploring / accessing our ASV, TAXA, and SAMPLE dataframes saved as phyloseq object ###

sData <- as(sample_data(physeq_microbiome), 'data.frame')
head(sData)
tData <- as(tax_table(physeq_microbiome), 'matrix')
aData <- as(otu_table(physeq_microbiome), 'matrix')

#Since we will be using the final_metadata dataframe for the sample, we must prepare it for use as sample_data object

row.names(final_metadata) <- final_metadata$SampleName
final_updated_metadata <- final_metadata %>% dplyr::select(-SampleName)
samples = sample_data(final_updated_metadata)
head(final_updated_metadata)
ASV <- as.matrix(aData)
TAX <- as.matrix(tData)

### Adding the new metadata SAMPLE dataframe to be saved with ASV and TAXA table as a new phyloseq object ###
worker_microbiome.ps <- merge_phyloseq(otu_table(physeq_microbiome, taxa_are_rows = F),sample_data(final_updated_metadata),tax_table(physeq_microbiome))
saveRDS(worker_microbiome.ps, "ilya_worker_microbiome_phyloseq_7_21_22")

#We check to see if the ASV matrix row sample names are ordered in the same way as the sample data matrix
#rownames(final_updated_metadata)
#rownames(ASV)
#genomic_idx <- match(rownames(final_updated_metadata), rownames(ASV))
#ASV_ordered <- ASV[genomic_idx , ]#reordering the ASV table to match the metadata (sample) file
#View(ASV_ordered)
#physeq_final_microbiome<- phyloseq(ASV_ordered, TAX, samples) 
#physeq_final_microbiome

## remove USDA-blind samples:

physeq_edit <- subset_samples(worker_microbiome.ps, !Worker%in% c("USDA_Blind") ) # removed the USDA--does not have metadata
## total 9423 ASVs

physeq_bacteriome <- subset_taxa(physeq_edit,  !is.na(Kingdom) & !Kingdom %in% c("Eukaryota")) # removed the Eukaryota
## total 9401 after removing euk
## 1 asv didnot assigned to kingdom
## final 9400 asvs

#tax_glom(physeq, taxrank=rank_names(physeq)[1], NArm=TRUE, bad_empty=c(NA, "", " ", "\t"))

## relative abundance plots:

dna <- Biostrings::DNAStringSet(taxa_names(physeq_bacteriome))
names(dna) <- taxa_names(physeq_bacteriome)
physeq_bacteriome_edit <- merge_phyloseq(physeq_bacteriome, dna)
taxa_names(physeq_bacteriome_edit) <- paste0("ASV", seq(ntaxa(physeq_bacteriome_edit)))
physeq_bacteriome_edit
otu_table(physeq_bacteriome_edit)

## running decontaminant:
library(decontam)
df <- as.data.frame(sample_data(physeq_bacteriome_edit)) # Put sample_data into a ggplot-friendly data.frame
df$LibrarySize <- sample_sums(physeq_bacteriome_edit)
df <- df[order(df$LibrarySize),]
df$Index <- seq(nrow(df))
library(ggplot2)
ggplot(data=df, aes(x=Index, y=LibrarySize, color=Sample_or_Control)) + geom_point()

sample_data(physeq_bacteriome_edit)$qpcr_16s_copies_ul    <- as.numeric(sample_data(physeq_bacteriome_edit)$qpcr_16s_copies_ul)

contamdf.freq <- isContaminant(physeq_bacteriome_edit, method="frequency", conc="qpcr_16s_copies_ul")
head(contamdf.freq)

table(contamdf.freq$contaminant)
# FALSE  TRUE 
# 9358    42 

head(which(contamdf.freq$contaminant))
# head(which(contamdf.freq$contaminant))
# [1] 176 180 181 212 247 297

plot_frequency(physeq_bacteriome, taxa_names(physeq_bacteriome)[c(176,180,181)], conc="qpcr_16s_copies_ul") + 
  xlab("DNA Concentration (qPCR)")

ps.noncontam.worker <- prune_taxa(!contamdf.freq$contaminant, physeq_bacteriome_edit)
ps.noncontam.worker
taxonomytable_microbiome<-as.data.frame(tax_table(ps.noncontam.worker))
write.csv(taxonomytable_microbiome, "C:/Users/is233/OneDrive/Desktop/AG_WORKER_METAGENOMICS_STUDY/MANUSCRIPT/ANALYSIS/WORKER_MICROBIOME/final_taxonomy_table.csv", row.names = TRUE)
df$LibrarySize <- sample_sums(ps.noncontam.worker)
sum(sample_sums(ps.noncontam.worker))
#2391450 reads after removing contaminant
hist(contamdf.freq$p.freq, breaks = 30)

saveRDS(ps.noncontam.worker, "ilya_worker_microbiome_phyloseq_noncontam_7_21_22")
ps.noncontam.worker <- readRDS("ilya_worker_microbiome_phyloseq_noncontam_7_21_22")
otu_table(ps.noncontam.worker)
ps.noncontam.worker


########################################################
##                                                     #
###                                                    #
####Co-occurrence network analysis                     #
###                                                    # 
##                                                     #
########################################################

install.packages("https://cran.r-project.org/src/contrib/Archive/KMDA/KMDA_1.0.tar.gz", repos = NULL, type="source")
install_github("umerijaz/microbiomeSeq")
library(microbiomeSeq)
install.packages("impute")
library(adespatial)
BiocManager::install("Rhdf5lib")
library(Rhdf5lib)
install.packages("igraph", "nloptr", "Rcpp", "RCurl","Rhdf5lib","stringi","tibble", "KMDA", "adespatial")
BiocManager::install("DESeq2")
library(DESeq2)
BiocManager::install("GO.db")
library(GO.db)
BiocManager::install("impute")
library(impute)
BiocManager::install("preprocessCore")
library(preprocessCore)
BiocManager::install("graph")
library(graph)
library(phyloseq)
install_github("zdk123/SpiecEasi")
library(SpiecEasi)
install.packages("visNetwork")
library("igraph")

require(visNetwork)

#Trying above using SpiecEasi
#Using SpiecEasi ps object to explore package

#Select only samples that are WORKDAY_START in CollectionPhase
ps.noncontam.worker.WORKSTART<- subset_samples(ps.noncontam.worker, CollectionPhase == "WORKDAY_START")
ps.noncontam.worker.WORKSTART
#Running SpiecEasi on ps object using wrappers for phyloseq
head(tax_table(ps.noncontam.worker.WORKSTART))
head(otu_table(ps.noncontam.worker.WORKSTART))

#Reduce the number of ASVs to 100+ count
ps.noncontam.worker.WORKSTART.f<- prune_taxa(taxa_sums(ps.noncontam.worker.WORKSTART)>100, ps.noncontam.worker.WORKSTART)
ps.noncontam.worker.WORKSTART.f
saveRDS(ps.noncontam.worker.WORKSTART.f, "SE_ps.noncontam.worker.WORKSTART.f")
#> ps.noncontam.worker.WORKSTART.f
#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 385 taxa and 10 samples ]
#sample_data() Sample Data:       [ 10 samples by 90 sample variables ]
#tax_table()   Taxonomy Table:    [ 385 taxa by 7 taxonomic ranks ]
#refseq()      DNAStringSet:      [ 385 reference sequences ]

#Select only samples that are WORKDAY_END in CollectionPhase
ps.noncontam.worker.WORKEND<- subset_samples(ps.noncontam.worker, CollectionPhase == "WORKDAY_END")

#Reduce the number of ASVs to 100
ps.noncontam.worker.WORKEND.f<- prune_taxa(taxa_sums(ps.noncontam.worker.WORKEND)>100, ps.noncontam.worker.WORKEND)
saveRDS(ps.noncontam.worker.WORKEND.f, "SE_ps.noncontam.worker.WORKEND.f.rds")
ps.noncontam.worker.WORKEND.f
#> ps.noncontam.worker.WORKEND.f
#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 205 taxa and 10 samples ]
#sample_data() Sample Data:       [ 10 samples by 90 sample variables ]
#tax_table()   Taxonomy Table:    [ 205 taxa by 7 taxonomic ranks ]
#refseq()      DNAStringSet:      [ 205 reference sequences ]

#Select only samples that are POST_SHOWER in CollectionPhase
ps.noncontam.worker.POSTSHOWER<- subset_samples(ps.noncontam.worker, CollectionPhase == "POST_SHOWER")

#Reduce the number of ASVs to 100
ps.noncontam.worker.POSTSHOWER.f<- prune_taxa(taxa_sums(ps.noncontam.worker.POSTSHOWER)>100, ps.noncontam.worker.POSTSHOWER)
saveRDS(ps.noncontam.worker.POSTSHOWER.f, "SE_ps.noncontam.worker.POSTSHOWER.f.rds")
ps.noncontam.worker.POSTSHOWER.f
#> ps.noncontam.worker.POSTSHOWER.f
#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 254 taxa and 10 samples ]
#sample_data() Sample Data:       [ 10 samples by 90 sample variables ]
#tax_table()   Taxonomy Table:    [ 254 taxa by 7 taxonomic ranks ]
#refseq()      DNAStringSet:      [ 254 reference sequences ]

#Select only samples that are SWINE in CollectionPhase
ps.noncontam.worker.SWINE<- subset_samples(ps.noncontam.worker, CollectionPhase == "SWINE")

#Reduce the number ofASVs to 100
ps.noncontam.worker.SWINE.f<- prune_taxa(taxa_sums(ps.noncontam.worker.SWINE)>100, ps.noncontam.worker.SWINE)
saveRDS(ps.noncontam.worker.SWINE.f, "SE_ps.noncontam.worker.SWINE.f.rds")
ps.noncontam.worker.SWINE.f
#> ps.noncontam.worker.SWINE.f
#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 283 taxa and 10 samples ]
#sample_data() Sample Data:       [ 10 samples by 90 sample variables ]
#tax_table()   Taxonomy Table:    [ 283 taxa by 7 taxonomic ranks ]
#refseq()      DNAStringSet:      [ 283 reference sequences ]

#Select only samples that are ENVIRONMENT in CollectionPhase
ps.noncontam.worker.ENVIRONMENT<- subset_samples(ps.noncontam.worker, CollectionPhase == "ENVIRONMENT")

#Reduce the number ofASVs to 100
ps.noncontam.worker.ENVIRONMENT.f<- prune_taxa(taxa_sums(ps.noncontam.worker.ENVIRONMENT)>100, ps.noncontam.worker.ENVIRONMENT)
saveRDS(ps.noncontam.worker.ENVIRONMENT.f, "SE_ps.noncontam.worker.ENVIRONMENT.f.rds")
ps.noncontam.worker.ENVIRONMENT.f
#> ps.noncontam.worker.ENVIRONMENT.f
#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 272 taxa and 2 samples ]
#sample_data() Sample Data:       [ 2 samples by 90 sample variables ]
#tax_table()   Taxonomy Table:    [ 272 taxa by 7 taxonomic ranks ]
#refseq()      DNAStringSet:      [ 272 reference sequences ]

set.seed(1244)


## estimated richness from non-contaminant files
## Estimate the alpha diversity:

biome_phylum.ps <-  tax_glom(ps.noncontam.worker, "Phylum")
biome_class.ps <-   tax_glom(ps.noncontam.worker, "Class")
biome_family.ps <-  tax_glom(ps.noncontam.worker, "Family")
biome_genus.ps <-   tax_glom(ps.noncontam.worker, "Genus")
biome_species.ps <- tax_glom(ps.noncontam.worker, "Species")
biome_asv.ps <- ps.noncontam.worker

sample_df<- sample_data(ps.noncontam.worker)
head(sample_df)

####Rarefaction

mpse<- subsetted_biome_genus.ps %>% as.MPSE()
mpse
data(mpse)
mpse%<>%
  mp_cal_rarecurve(
    .abundance = Abundance,
    chunks = 100, action = "add", force=TRUE
  ) 


data(mpse)
set.seed(1234)
mpse %>% 
  mp_plot_rarecurve(
    .rare = AbundanceRarecurve,   
    .alpha = "Observe", 
    .group = CollectionPhase,
    plot.group = TRUE) + theme_bw() + coord_cartesian(xlim=c(0,100000), ylim = c(0,600))+  scale_color_manual(values=c("#789ED0", "#F79D40", "#72BF56", "#E4655A", "#AB8BCB")) +  scale_fill_manual(values=c("#789ED0", "#F79D40", "#72BF56", "#E4655A", "#AB8BCB"))  



subsetted_biome_genus.ps<- subset_samples(biome_asv.ps, CollectionPhase=="WORKDAY_START" | CollectionPhase=="WORKDAY_END" | CollectionPhase=="POST_SHOWER" |CollectionPhase=="SWINE" | CollectionPhase=="ENVIRONMENT")
rare.level <- min(sample_sums(subsetted_biome_genus.ps))
rareres <- get_rarecurve(obj=subsetted_biome_genus.ps, chunks=150)
prare1 <- ggrarecurve(obj=rareres, factorNames="CollectionPhase", indexNames=c("Observe"), se=T, shadow = TRUE) +
  scale_fill_manual(values=c("#789ED0", "#F79D40", "#72BF56", "#E4655A", "#AB8BCB"))+
  scale_color_manual(values=c("#789ED0", "#F79D40", "#72BF56", "#E4655A", "#AB8BCB"))+
  theme_bw()+ theme(axis.text=element_text(size=8), panel.grid=element_blank(),
        strip.background = element_rect(colour=NA,fill="grey"), strip.text = element_blank()) + xlab("Read depth") + ylab("Unique number of genera")+ theme(strip.text = element_blank())

prare1



p <- ggrare(subsetted_biome_genus.ps, sample=raremax, step = 200, color = "CollectionPhase", label = NULL, se = TRUE)
p <- p + facet_wrap(~CollectionPhase, nrow = 1) + theme_bw() + ylab("Genus richness")+  scale_fill_manual(values=c("#789ED0", "#F79D40", "#72BF56", "#E4655A", "#AB8BCB"))+ scale_color_manual(values=c("#789ED0", "#F79D40", "#72BF56", "#E4655A", "#AB8BCB")) 
plot(p)
p

subsetted_biome_genus.ps %>%
  otu_table() %>%
  t() %>%
  vegan::rarecurve()
source("https://raw.githubusercontent.com/mahendra-mariadassou/phyloseq-extended/master/load-extra-functions.R")




### Relative abudance-- phylum
biome_phylum_samples <- subset_samples(biome_phylum.ps, !Worker%in% c("NegCtrl", "MockComm") )

phy_relative <- transform_sample_counts(biome_phylum_samples, function(x) x / sum(x) )
phy_relative_long <- psmelt(phy_relative)
phy_relative_long <- phy_relative_long %>%
  group_by(Phylum) %>%
  mutate(mean_relative_abund = mean(Abundance))

phy_relative_long$Phylum <- as.character(phy_relative_long$Phylum)
phy_relative_long$mean_relative_abund <- as.numeric(phy_relative_long$mean_relative_abund)
phy_relative_long$Phylum[phy_relative_long$mean_relative_abund < 0.005] <- "Phyla (< 0.5%)"  ## mean_relative_abund < 0.005


phy_relative_long$CollectionPhase <- factor (phy_relative_long$CollectionPhase, levels= c("WORKDAY_START",
                                                                                          "WORKDAY_END",
                                                                                          "POST_SHOWER",
                                                                                          'SWINE',
                                                                                          "ENVIRONMENT"))


phy_relative_long$OccupationalTask <- factor (phy_relative_long$OccupationalTask, levels= c("FARROWING",
                                                                                            "GESTATION",
                                                                                            "SOW_GESTATION",
                                                                                            'SWINE_PEN',
                                                                                            "ENV_FARROWING",
                                                                                            "ENV_GESTATION",
                                                                                            "MOVING_FOSTERING",
                                                                                            "VISITING VETERINARIAN"))
phy_relative_long %>%
  filter(!CollectionPhase %in% c("NegCtrl"))%>%
  ggplot(aes(x = Sample, y = Abundance*100, fill = Phylum)) +
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
scale_fill_tableau(palette = "Color Blind", type = "regular",direction = 1) 

topptx(filename = "phy_relative_plot1.pptx", width=8, height=4)




collectphase_names <- list("WORKDAY_START"= "Workday start", "WORKDAY_END"="Workday end", "POST_SHOWER"="Post-shower", "SWINE"="Swine", "ENVIRONMENT"="Environment")

collectphase_labeller <- function(variable,value){
  return(collectphase_names[value])
}

subsetted_biome_asv.ps<- subset_samples(biome_asv.ps, CollectionPhase=="WORKDAY_START" | CollectionPhase=="WORKDAY_END" | CollectionPhase=="POST_SHOWER" |CollectionPhase=="SWINE" | CollectionPhase=="ENVIRONMENT")

phylumtaxa <- get_taxadf(obj=subsetted_biome_asv.ps, taxlevel=3)
phylumtaxa
# The 30 most abundant taxonomy will be visualized by default (parameter `topn=30`). 
pphylum <- ggbartax(obj=phylumtaxa, facetNames="CollectionPhase", topn=19) +
    scale_fill_tableau(palette = "Tableau 20", type = "regular",direction = 1) +
  guides(fill= guide_legend(keywidth = 0.5, keyheight = 0.5)) + xlab("Collection phase") + ylab("Relative abundance (%)")+ theme(strip.text = element_blank()) 
pphylum

Phylumoutcome <- ggbartax(obj=phylumtaxa, facetNames="CollectionPhase", plotgroup=TRUE, topn=10) + xlab(NULL) +
  ylab("relative abundance (%)") +
   scale_fill_tableau(palette = "Tableau 20", type = "regular",direction = 1) +
  guides(fill= guide_legend(keywidth = 0.5, keyheight = 0.5)) +xlab("Collection phase") + ylab("Relative abundance (%)")+ theme(strip.text = element_blank())+
  guides(fill= guide_legend(keywidth = 0.5, keyheight = 0.5, ncol=2)) + scale_x_discrete("Collection phase", labels= c("WORKDAY_START"= "Workday start", "WORKDAY_END"="Workday end", "POST_SHOWER"="Post-shower", "SWINE"="Swine", "ENVIRONMENT"="Environment"))

Phylumoutcome


## family:

biome_family_samples <- subset_samples(biome_family.ps, !Worker%in% c("NegCtrl") )

phy_relative <- transform_sample_counts(biome_family_samples, function(x) x / sum(x) )
phy_relative_long <- psmelt(phy_relative)
phy_relative_long <- phy_relative_long %>%
  group_by(Family) %>%
  mutate(mean_relative_abund = mean(Abundance))

phy_relative_long$Family <- as.character(phy_relative_long$Family)
phy_relative_long$mean_relative_abund <- as.numeric(phy_relative_long$mean_relative_abund)
phy_relative_long$Family[phy_relative_long$mean_relative_abund < 0.01] <- "Taxa (< 1%)"  ## mean_relative_abund < 0.005

phy_relative_long$CollectionPhase <- factor (phy_relative_long$CollectionPhase, levels= c("WORKDAY_START",
                                                                                          "WORKDAY_END",
                                                                                          "POST_SHOWER",
                                                                                          'SWINE',
                                                                                          "ENVIRONMENT"))

phy_relative_long %>%
  filter(!CollectionPhase %in% c("NegCtrl", "NA"))%>%
  ggplot(aes(x = Sample, y = Abundance*100, fill = Family)) +
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



## genus:

biome_genus_samples <- subset_samples(biome_genus.ps, !Worker%in% c("NegCtrl") )

phy_relative <- transform_sample_counts(biome_genus_samples, function(x) x / sum(x) )
phy_relative_long <- psmelt(phy_relative)
phy_relative_long <- phy_relative_long %>%
  group_by(Genus) %>%
  mutate(mean_relative_abund = mean(Abundance))

phy_relative_long$Genus <- as.character(phy_relative_long$Genus)
phy_relative_long$mean_relative_abund <- as.numeric(phy_relative_long$mean_relative_abund)
phy_relative_long$Genus[phy_relative_long$mean_relative_abund < 0.001] <- "Taxa (< 0.1%)"  ## mean_relative_abund < 0.005

phy_relative_long$CollectionPhase <- factor (phy_relative_long$CollectionPhase, levels= c("WORKDAY_START",
                                                                                          "WORKDAY_END",
                                                                                          "POST_SHOWER",
                                                                                          'SWINE',
                                                                                          "ENVIRONMENT"))

phy_relative_long %>%
  filter(!CollectionPhase %in% c("NegCtrl"))%>%
  ggplot(aes(x = Sample, y = Abundance*100, fill = Genus)) +
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


#diversity indices
div_phylum <- estimate_richness(biome_phylum.ps, measures = c("Observed", "Shannon", "InvSimpson"))
div_class <- estimate_richness(biome_class.ps, measures = c("Observed", "Shannon", "InvSimpson"))
div_family <- estimate_richness(biome_family.ps, measures = c("Observed", "Shannon", "InvSimpson"))
div_genus <- estimate_richness(biome_genus.ps, measures = c("Observed", "Shannon", "InvSimpson"))
div_species <- estimate_richness(biome_species.ps, measures = c("Observed", "Shannon", "InvSimpson"))
div_asv <- estimate_richness(biome_asv.ps, measures = c("Observed", "Shannon", "InvSimpson"))

install.packages("data.table")
library("data.table")
library(ggplot2)
install.packages("ggthemes")
library(ggthemes)
library(tibble)

#phylum

rownames(div_phylum)<- sub('X','',rownames(div_phylum))
div_phylum$SampleName <- row.names(div_phylum)
expanded_metadata <- left_join(updated_metadata, div_phylum, by = "SampleName")
row.names(expanded_metadata) <- expanded_metadata$SampleName
tail(expanded_metadata)

setnames(expanded_metadata, old = c('Observed','Shannon','InvSimpson'), new = c('Observed_phylum','shannon_phylum','InvSimpson_phylum'))
tail(expanded_metadata)

## class
rownames(div_class)<- sub('X','',rownames(div_class))
div_class$SampleName <- row.names(div_class)
expanded_metadata <- left_join(expanded_metadata, div_class, by = "SampleName")
row.names(expanded_metadata) <- expanded_metadata$SampleName
tail(expanded_metadata)

setnames(expanded_metadata, old = c('Observed','Shannon','InvSimpson'), new = c('Observed_class','shannon_class','InvSimpson_class'))
tail(expanded_metadata)


## family
rownames(div_family)<- sub('X','',rownames(div_family))
div_family$SampleName <- row.names(div_family)
expanded_metadata <- left_join(expanded_metadata, div_family, by = "SampleName")
row.names(expanded_metadata) <- expanded_metadata$SampleName
tail(expanded_metadata)

setnames(expanded_metadata, old = c('Observed','Shannon','InvSimpson'), new = c('Observed_family','shannon_family','InvSimpson_family'))
tail(expanded_metadata)


##genus
rownames(div_genus)<- sub('X','',rownames(div_genus))
div_genus$SampleName <- row.names(div_genus)
expanded_metadata <- left_join(expanded_metadata, div_genus, by = "SampleName")
row.names(expanded_metadata) <- expanded_metadata$SampleName
tail(expanded_metadata)
setnames(expanded_metadata, old = c('Observed','Shannon','InvSimpson'), new = c('Observed_genus','shannon_genus','InvSimpson_genus'))
tail(expanded_metadata)

##asv
rownames(div_asv)<- sub('X','',rownames(div_asv))
div_asv$SampleName <- row.names(div_asv)
expanded_metadata <- left_join(expanded_metadata, div_asv, by = "SampleName")
row.names(expanded_metadata) <- expanded_metadata$SampleName
tail(expanded_metadata)
setnames(expanded_metadata, old = c('Observed','Shannon','InvSimpson'), new = c('Observed_asv','shannon_asv','InvSimpson_asv'))
tail(expanded_metadata)

saveRDS(expanded_metadata, "expanded_metadata_decontam_ilya_16s.rds")
expanded_metadata<- readRDS("expanded_metadata_decontam_ilya_16s.rds")


## plot:

#Phylum richness by occupational task
obs_phy <-  
  expanded_metadata %>%
  filter(!OccupationalTask %in% c("MockComm", "NegCtrl", "ENV_FARROWING", "ENV_GESTATION"))%>%
 mutate(OccupationalTask= fct_relevel(OccupationalTask,"VISITING VETERINARIAN", "FARROWING", "GESTATION", "MOVING_FOSTERING", "SOW_GESTATION", "SWINE_PEN")) %>% #performed to change x axis order
   ggplot(aes(OccupationalTask, y=Observed_phylum, fill=OccupationalTask))+
  geom_jitter(shape=21, color= "black", size= 3, alpha=0.6, position = pd) +
  geom_boxplot(alpha=0.45, shape=31, color= "black", outlier.colour = 'white', width=0.6, coef=5000, fatten=1, lwd=1)+ylim(0,20)+
  scale_fill_tableau(palette = "Classic 10 Medium", type = "regular",direction = 1)+  #Tableau 10
  scale_color_tableau(palette = "Classic 10 Medium", type = "regular", direction = 1)+
  xlab("Occupational task")+ ylab("Phylum richness")+
  theme_classic()+
  theme(legend.position = "none")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.line = element_line(size = 1))+
  scale_x_discrete("Occupational description", labels= c("FARROWING"= "Farrowing", "GESTATION"="Gestation", "MOVING_FOSTERING"="Moving and fostering", "SOW_GESTATION"="Sow gestation", "SWINE_PEN"="Swine", "VISITING VETERINARIAN"="Visiting veterinarian"))

#facet_wrap(~Exposure,scales='free_x', nrow=1)
obs_phy

#topptx(filename = "obs_phybyOccupationalTask_plot.pptx", width=5, height=4) 
expanded_metadata_edit_occtask <- subset (expanded_metadata, !OccupationalTask%in% c(c("MockComm", "NegCtrl", "ENV_FARROWING", "ENV_GESTATION")))
str(expanded_metadata_edit_occtask)
model_obs_phy <-
  lmer(Observed_phylum ~ OccupationalTask +(1|Worker), data = expanded_metadata_edit_occtask, REML=T)
summary(model_obs_phy)
anova(model_obs_phy, type='III',test="F")

lsmeans(model_obs_phy, pairwise~OccupationalTask, adjust="tukey")
VarCorr(model_obs_phy)
#Type III Analysis of Variance Table with Satterthwaite's method
#                 Sum Sq Mean Sq NumDF  DenDF F value Pr(>F)
#OccupationalTask 22.752  4.5503     5 10.634  1.4692 0.2781
#No significant differences in observed phylum richness by occupational description
obs_phy + geom_text(data = tibble(x=1.75, y=3), size=3, aes(x=x, y=y, label= "Type III ANOVA P >0.05", fontface="italic"),inherit.aes=FALSE) 

#Phylum richness by Collection phase
obs_phy1 <-  
  expanded_metadata %>%
  filter(!CollectionPhase %in% c("MockComm", "NegCtrl"))%>%
  mutate(CollectionPhase= fct_relevel(CollectionPhase,"WORKDAY_START", "WORKDAY_END", "POST_SHOWER", "SWINE", "ENVIRONMENT"))%>% #performed to change x axis order 
  ggplot(aes(CollectionPhase, y=Observed_phylum, fill=CollectionPhase))+
  geom_boxplot(alpha=0.45, shape=31, color= "black", outlier.colour = 'white', width=0.6, coef=5000, fatten=1, lwd=1)+ylim(0,NA)+  xlab("Collection phase")+ ylab("Phylum richness")+
  geom_line(aes(group=Worker), size=0.05, alpha=0.7, color= "#bdbdbd")+
  geom_jitter(shape=21, color= "black", size= 3, alpha=0.6, position = pd)+
  scale_fill_tableau(palette = "Classic 10 Medium", type = "regular",direction = 1)+  #Tableau 10
  scale_color_tableau(palette = "Classic 10 Medium", type = "regular", direction = 1)+
  #scale_color_manual(values=c("grey35", "red",  "#59A14F", "blue",  "purple"))+
  #scale_fill_manual(values=c("grey35", "red",  "#59A14F", "blue",  "purple"))+
  theme_bw()+
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        text=element_text(family="Times New Roman",  size=12))+
  theme(legend.position = "none")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.line = element_line(size = 1))+
  scale_x_discrete("Collection phase", labels= c("WORKDAY_START"= "Workday start", "WORKDAY_END"="Workday end", "POST_SHOWER"="Post-shower", "SWINE"="Swine", "ENVIRONMENT"="Environment"))
#facet_wrap(~Exposure,scales='free_x', nrow=1)
obs_phy1
topptx(filename = "obs_phybyOccupationalTask_plot1.pptx", width=3, height=4) 

expanded_metadata_edit_occphase <- subset (expanded_metadata, !CollectionPhase%in% c(c("MockComm", "NegCtrl")))
str(expanded_metadata_edit_occphase)

model_obs_phy1 <-
  lmer(Observed_phylum ~ CollectionPhase +(1|Worker), data = expanded_metadata_edit_occphase, REML=T)
summary(model_obs_phy1)
anova(model_obs_phy1, type='III',test="F")
lsmeans(model_obs_phy1, pairwise~CollectionPhase, adjust="BH")
VarCorr(model_obs_phy1)
#Type III Analysis of Variance Table with Satterthwaite's method
#                Sum Sq Mean Sq NumDF  DenDF F value  Pr(>F)  
#CollectionPhase 44.862  11.215     4 26.051  3.4808 0.02095 *
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#$contrasts
#contrast                    estimate    SE   df t.ratio p.value
#ENVIRONMENT - POST_SHOWER        3.5 1.563 34.8   2.239  0.1895
#ENVIRONMENT - SWINE              5.2 1.563 34.8   3.327  0.0166
#ENVIRONMENT - WORKDAY_END        2.9 1.563 34.8   1.855  0.3596
#ENVIRONMENT - WORKDAY_START      3.3 1.563 34.8   2.111  0.2383
#POST_SHOWER - SWINE              1.7 0.902 34.8   1.884  0.3446
#POST_SHOWER - WORKDAY_END       -0.6 0.803 20.2  -0.747  0.9425
#POST_SHOWER - WORKDAY_START     -0.2 0.803 20.2  -0.249  0.9991
#SWINE - WORKDAY_END             -2.3 0.902 34.8  -2.549  0.1029
#SWINE - WORKDAY_START           -1.9 0.902 34.8  -2.105  0.2407
#WORKDAY_END - WORKDAY_START      0.4 0.803 20.2   0.498  0.9866
obs_phy1 + geom_text(data = tibble(x=1.75, y=5), size=4, aes(x=x, y=y, label= "Type III ANOVA P = 0.017", fontface="italic"),inherit.aes=FALSE) + geom_line(data=tibble(x=c(4,5), y=c(15.5,15.5)), aes(x=x, y=y), inherit.aes=FALSE) + geom_text(data = tibble(x=4.5, y=16), size=5, aes(x=x, y=y, label= "*"),inherit.aes=FALSE) + 
  theme(axis.text.x = element_text(size = 16), axis.text.y = element_text(size = 16), text=element_text(family="Helvetica Neue Medium",  size=18))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.line = element_line(size = 1))+
  theme(legend.position = "none") +  scale_x_discrete("Collection phase", labels= c("WORKDAY_START"= "Workday start", "WORKDAY_END"="Workday end", "POST_SHOWER"="Post-shower", "SWINE"="Swine", "ENVIRONMENT"="Environment"))


#Phylum Shannon diversity by Collection phase
shannon_phylum1 <-  
  expanded_metadata %>%
  filter(!CollectionPhase %in% c("MockComm", "NegCtrl"))%>%
  mutate(CollectionPhase= fct_relevel(CollectionPhase,"WORKDAY_START", "WORKDAY_END", "POST_SHOWER", "SWINE", "ENVIRONMENT"))%>% #performed to change x axis order 
  ggplot(aes(CollectionPhase, y=shannon_phylum, fill=CollectionPhase))+
  geom_boxplot(alpha=0.45, shape=31, color= "black", outlier.colour = 'white', width=0.6, coef=5000, fatten=1, lwd=1)+ylim(0,NA)+  xlab("Collection phase")+ ylab("Phylum Shannon diversity")+
  geom_line(aes(group=Worker), size=0.05, alpha=0.7, color= "#bdbdbd")+
  geom_jitter(shape=21, color= "black", size= 3, alpha=0.6, position = pd)+
  scale_fill_tableau(palette = "Classic 10 Medium", type = "regular",direction = 1)+  #Tableau 10
  scale_color_tableau(palette = "Classic 10 Medium", type = "regular", direction = 1)+
  #scale_color_manual(values=c("grey35", "red",  "#59A14F", "blue",  "purple"))+
  #scale_fill_manual(values=c("grey35", "red",  "#59A14F", "blue",  "purple"))+
  theme_bw()+
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        text=element_text(family="Times New Roman",  size=12))+
  theme(legend.position = "none")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.line = element_line(size = 1))+
  scale_x_discrete("Collection phase", labels= c("WORKDAY_START"= "Workday start", "WORKDAY_END"="Workday end", "POST_SHOWER"="Post-shower", "SWINE"="Swine", "ENVIRONMENT"="Environment"))
shannon_phylum1
topptx(filename = "Shannon_PhylumybyCollectionPhase_plot.pptx", width=3, height=4) 

model_shannon_phy1 <-
  lmer(shannon_phylum ~ CollectionPhase +(1|Worker), data = expanded_metadata_edit_occphase, REML=T)
summary(model_shannon_phy1)
shannon_phy1anova<- anova(model_shannon_phy1, type='III',test="F") #0.002048154

shannon_pairphy1anova<- lsmeans(model_shannon_phy1, pairwise~CollectionPhase, adjust="BH")
VarCorr(model_shannon_phy1)
#contrast                    estimate     SE   df t.ratio p.value
#ENVIRONMENT - POST_SHOWER     0.1040 0.1755 31.9   0.592  0.6972
#ENVIRONMENT - SWINE          -0.3861 0.1755 31.9  -2.199  0.0881
#ENVIRONMENT - WORKDAY_END    -0.0402 0.1755 31.9  -0.229  0.8203
#ENVIRONMENT - WORKDAY_START   0.0440 0.1755 31.9   0.250  0.8203
#POST_SHOWER - SWINE          -0.4901 0.1014 31.9  -4.835  0.0003
#POST_SHOWER - WORKDAY_END    -0.1442 0.0824 19.4  -1.751  0.1915
#POST_SHOWER - WORKDAY_START  -0.0600 0.0824 19.4  -0.729  0.6781
#SWINE - WORKDAY_END           0.3459 0.1014 31.9   3.412  0.0059
#SWINE - WORKDAY_START         0.4300 0.1014 31.9   4.243  0.0009
#WORKDAY_END - WORKDAY_START   0.0842 0.0824 19.4   1.022  0.5324
shannon_phylum1 + geom_text(data = tibble(x=1.75, y=0.1), size=3, aes(x=x, y=y, label= "Type III ANOVA P = 0.002", fontface="italic"),inherit.aes=FALSE) + geom_line(data=tibble(x=c(4,1), y=c(2.25,2.25)), aes(x=x, y=y), inherit.aes=FALSE)+ geom_line(data=tibble(x=c(4,2), y=c(2,2)), aes(x=x, y=y), inherit.aes=FALSE)+ geom_line(data=tibble(x=c(4,3), y=c(1.75,1.75)), aes(x=x, y=y), inherit.aes=FALSE) + geom_text(data = tibble(x=2.5, y=2.3), size=4, aes(x=x, y=y, label= "**"),inherit.aes=FALSE)+geom_text(data = tibble(x=3, y=2.05), size=4, aes(x=x, y=y, label= "**"),inherit.aes=FALSE)+geom_text(data = tibble(x=3.5, y=1.8), size=4, aes(x=x, y=y, label= "**"),inherit.aes=FALSE) + theme(axis.text.x = element_text(size = 16), axis.text.y = element_text(size = 16), text=element_text(family="Helvetica Neue Medium",  size=18))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.line = element_line(size = 0.5))+
  theme(legend.position = "none") +  scale_x_discrete("Collection phase", labels= c("WORKDAY_START"= "Workday start", "WORKDAY_END"="Workday end", "POST_SHOWER"="Post-shower", "SWINE"="Swine", "ENVIRONMENT"="Environment"))

obs_family1 <-  
  expanded_metadata %>%
  filter(!CollectionPhase %in% c("MockComm", "NegCtrl"))%>%
  ggplot(aes(CollectionPhase, y=Observed_family, fill=CollectionPhase))+
  geom_boxplot(alpha=0.45, outlier.colour = "white")+ylim(0,NA)+
  geom_jitter(aes(color = CollectionPhase),size=2, alpha=0.6, position = pd)+
  #scale_fill_tableau(palette = "Superfishel Stone", type = "regular",direction = 1)+  #Tableau 10
  #scale_color_tableau(palette = "Superfishel Stone", type = "regular", direction = 1)+
  scale_color_manual(values=c("grey35", "red",  "#59A14F", "blue",  "purple"))+
  scale_fill_manual(values=c("grey35", "red",  "#59A14F", "blue",  "purple"))+
  theme_classic()+
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        text=element_text(family="Times New Roman",  size=12))+
  theme(legend.position = "none")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
#facet_wrap(~Exposure,scales='free_x', nrow=1)
obs_family1
topptx(filename = "obs_FamilybyCollectionPhase_plot1.pptx", width=3, height=4) 

shannon_family1 <-  
  expanded_metadata %>%
  filter(!CollectionPhase %in% c("MockComm", "NegCtrl", "SWINE", "ENVIRONMENT"))%>%
  ggplot(aes(CollectionPhase, y=shannon_family, fill=CollectionPhase))+
  geom_boxplot(alpha=0.45, outlier.colour = "white")+ylim(0,4)+
  geom_jitter(aes(color = CollectionPhase),size=2, alpha=0.6, position = pd)+
  #scale_fill_tableau(palette = "Superfishel Stone", type = "regular",direction = 1)+  #Tableau 10
  #scale_color_tableau(palette = "Superfishel Stone", type = "regular", direction = 1)+
  scale_color_manual(values=c("grey35", "red",  "#59A14F", "blue",  "purple"))+
  scale_fill_manual(values=c("grey35", "red",  "#59A14F", "blue",  "purple"))+
  theme_classic()+
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        text=element_text(family="Times New Roman",  size=12))+
  theme(legend.position = "none")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
shannon_family1
topptx(filename = "Shannon_FamilybyCollectionPhase_plot1.pptx", width=3, height=4) 

obs_genus <-  
  expanded_metadata %>%
  #filter(!OccupationalTask %in% c("MockComm", "NegCtrl"))%>%
  ggplot(aes(OccupationalTask, y=Observed_genus, fill=OccupationalTask))+
  geom_boxplot(alpha=0.45, outlier.colour = "white")+ylim(0,NA)+
  geom_jitter(aes(color = OccupationalTask),size=2, alpha=0.6, position = pd)+
  scale_fill_tableau(palette = "Superfishel Stone", type = "regular",direction = 1)+  #Tableau 10
  scale_color_tableau(palette = "Superfishel Stone", type = "regular", direction = 1)+
  theme_bw()+
  theme(legend.position = "none")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
#facet_wrap(~Exposure,scales='free_x', nrow=1)
obs_genus
topptx(filename = "obs_GenusybyOccupationalTask_plot.pptx", width=5, height=4) 
obs_genus + theme(axis.text.x = element_text(size = 16), axis.text.y = element_text(size = 16), text=element_text(family="Helvetica Neue Medium",  size=18))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.line = element_line(size = 0.5))+
  theme(legend.position = "none") +  scale_x_discrete("Collection phase", labels= c("WORKDAY_START"= "Workday start", "WORKDAY_END"="Workday end", "POST_SHOWER"="Post-shower", "SWINE"="Swine", "ENVIRONMENT"="Environment"))

obs_genus1 <-  
  expanded_metadata %>%
  filter(!CollectionPhase %in% c("MockComm", "NegCtrl"))%>%
  mutate(CollectionPhase= fct_relevel(CollectionPhase,"WORKDAY_START", "WORKDAY_END", "POST_SHOWER", "SWINE", "ENVIRONMENT"))%>%
  ggplot(aes(CollectionPhase, y=Observed_genus, fill=CollectionPhase))+
  geom_boxplot(alpha=0.45, shape=31, color= "black", outlier.colour = 'white', width=0.6, coef=5000, fatten=1, lwd=1)+ylim(0,NA)+  xlab("Collection phase")+ ylab("Genus richness")+
  geom_line(aes(group=Worker), size=0.05, alpha=0.7, color= "#bdbdbd")+
  geom_jitter(shape=21, color= "black", size= 3, alpha=0.6, position = pd)+
  scale_fill_tableau(palette = "Classic 10 Medium", type = "regular",direction = 1)+  #Tableau 10
  scale_color_tableau(palette = "Classic 10 Medium", type = "regular", direction = 1)+
  theme_bw()+theme(axis.text.x = element_text(size = 16), axis.text.y = element_text(size = 16), text=element_text(family="Helvetica Neue Medium",  size=18))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.line = element_line(size = 0.5))+theme(legend.position = "none") +  scale_x_discrete("Collection phase", labels= c("WORKDAY_START"= "Workday start", "WORKDAY_END"="Workday end", "POST_SHOWER"="Post-shower", "SWINE"="Swine", "ENVIRONMENT"="Environment"))
#facet_wrap(~Exposure,scales='free_x', nrow=1)
obs_genus1

obs_genus1 + 


model_obs_genus1<-
  lmer(Observed_genus ~ CollectionPhase +(1|Worker), data = expanded_metadata_edit_occphase, REML=F)
summary(model_obs_genus1)
anova(model_obs_genus1)
#> anova(model_qpcr_copies1)
#Type III Analysis of Variance Table with Satterthwaite's method
#               Sum Sq Mean Sq NumDF  DenDF F value    Pr(>F)    
#CollectionPhase 0.00064617 0.0001077     6 26.423  2.3831 0.05678 
summary(anova)
lsmeans(model_quality_reads1, pairwise~CollectionPhase, adjust="tukey")
obs_genus1


shannon_genus <-  
  expanded_metadata %>%
  filter(!OccupationalTask %in% c("MockComm", "NegCtrl"))%>%
  ggplot(aes(OccupationalTask, y=shannon_genus, fill=OccupationalTask))+
  geom_boxplot(alpha=0.45, shape=31, color= "black", outlier.colour = 'white', width=0.6, coef=5000, fatten=1, lwd=1)+ylim(0,NA)+  xlab("Collection phase")+ ylab("Genus Shannon diversity")+
  geom_line(aes(group=Worker), size=0.05, alpha=0.7, color= "#bdbdbd")+
  geom_jitter(shape=21, color= "black", size= 3, alpha=0.6, position = pd)+
  scale_fill_tableau(palette = "Classic 10 Medium", type = "regular",direction = 1)+  #Tableau 10
  scale_color_tableau(palette = "Classic 10 Medium", type = "regular", direction = 1)+
  theme_bw()+
  theme(legend.position = "none")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
#facet_wrap(~Exposure,scales='free_x', nrow=1)
shannon_genus
topptx(filename = "Shannon_GenusybyOccupationalTask_plot.pptx", width=5, height=4) 

shannon_genus1 <-  
  expanded_metadata %>%
  filter(!CollectionPhase %in% c("MockComm", "NegCtrl"))%>%
  mutate(CollectionPhase= fct_relevel(CollectionPhase,"WORKDAY_START", "WORKDAY_END", "POST_SHOWER", "SWINE", "ENVIRONMENT"))%>% #performed to change x axis order 
  ggplot(aes(CollectionPhase, y=shannon_genus, fill=CollectionPhase))+
  geom_boxplot(alpha=0.45, shape=31, color= "black", outlier.colour = 'white', width=0.6, coef=5000, fatten=1, lwd=1)+ylim(0,NA)+  xlab("Collection phase")+ ylab("Genus Shannon diversity")+
  geom_line(aes(group=Worker), size=0.05, alpha=0.7, color= "#bdbdbd")+
  geom_jitter(shape=21, color= "black", size= 3, alpha=0.6, position = pd)+
  scale_fill_tableau(palette = "Classic 10 Medium", type = "regular",direction = 1)+  #Tableau 10
  scale_color_tableau(palette = "Classic 10 Medium", type = "regular", direction = 1)+
  theme_bw()+
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        text=element_text(family="Times New Roman",  size=12))+
  theme(legend.position = "none")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
shannon_genus1
#topptx(filename = "Shannon_GenusybyCollectionPhase_plot.pptx", width=4, height=5) 
topptx(filename = "Shannon_GenusybyCollectionPhase_plot1.pptx", width=3, height=4) 

shannon_genus1+ geom_text(data = tibble(x=1.75, y=0.1), size=3, aes(x=x, y=y, label= "Type III ANOVA P = 0.003", fontface="italic"),inherit.aes=FALSE) + geom_line(data=tibble(x=c(2,3), y=c(3.5,3.5)), aes(x=x, y=y), inherit.aes=FALSE)+ geom_line(data=tibble(x=c(4,3), y=c(3.8,3.8)), aes(x=x, y=y), inherit.aes=FALSE)+ geom_line(data=tibble(x=c(3,5), y=c(4.2,4.2)), aes(x=x, y=y), inherit.aes=FALSE) + geom_text(data = tibble(x=2.5, y=3.8), size=4, aes(x=x, y=y, label= "*"),inherit.aes=FALSE)+geom_text(data = tibble(x=3.5, y=4.0), size=4, aes(x=x, y=y, label= "**"),inherit.aes=FALSE)+geom_text(data = tibble(x=4, y=4.5), size=4, aes(x=x, y=y, label= "*"),inherit.aes=FALSE) + theme(axis.text.x = element_text(size = 16), axis.text.y = element_text(size = 16), text=element_text(family="Helvetica Neue Medium",  size=18))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.line = element_line(size = 0.5))+
  theme(legend.position = "none") +  scale_x_discrete("Collection phase", labels= c("WORKDAY_START"= "Workday start", "WORKDAY_END"="Workday end", "POST_SHOWER"="Post-shower", "SWINE"="Swine", "ENVIRONMENT"="Environment"))

## 

expanded_metadata$CollectionPhase


## model:
## lmer 

expanded_metadata_edit <- subset (expanded_metadata, !CollectionPhase%in% c(c("MockComm", "NegCtrl")))
str(expanded_metadata_edit)

##Observed_phylum 
##shannon_phylum 

##Observed_family 
##shannon_family

## Observed_genus 
## shannon_genus 


model_obs_phy <-
  lmer(Observed_phylum ~ CollectionPhase +(1|Worker), data = expanded_metadata_edit, REML=T)
summary(model_obs_phy)
anova(model_obs_phy, type='III',test="F")

lsmeans(model_obs_phy, pairwise~CollectionPhase, adjust="tukey")
VarCorr(model_obs_phy)


model_obs_family <-
  lmer(Observed_family ~ CollectionPhase +(1|Worker), data = expanded_metadata_edit, REML=T)
summary(model_obs_family)
anova(model_obs_family, type='III',test="F")

lsmeans(model_obs_family, pairwise~CollectionPhase, adjust="tukey")
VarCorr(model_obs_family)



model_obs_genus1 <-
  lmer(obs_genus1 ~ CollectionPhase +(1|Worker), data = expanded_metadata_edit_occphase, REML=T)
summary(model_obs_genus1)
anova(model_obs_genus1, type='III',test="F")

lsmeans(model_obs_genus, pairwise~CollectionPhase, adjust="tukey")
VarCorr(model_obs_genus)


model_obs_asv <-
  lmer(Observed_asv ~ CollectionPhase +(1|Worker), data = expanded_metadata_edit, REML=T)
summary(model_obs_asv)
anova(model_obs_asv, type='III',test="F")

lsmeans(model_obs_asv, pairwise~CollectionPhase, adjust="tukey")
VarCorr(model_obs_asv)


## shannon:
model_shannon_phy <-
  lmer(shannon_phylum ~ CollectionPhase +(1|Worker), data = expanded_metadata_edit, REML=T)
summary(model_shannon_phy)
anova(model_shannon_phy, type='III',test="F")

lsmeans(model_shannon_phy, pairwise~CollectionPhase, adjust="tukey")
VarCorr(model_shannon_phy)


model_shannon_family <-
  lmer(shannon_family ~ CollectionPhase +(1|Worker), data = expanded_metadata_edit, REML=T)
summary(model_shannon_family)
anova(model_shannon_family, type='III',test="F")

lsmeans(model_shannon_family, pairwise~CollectionPhase, adjust="tukey")
VarCorr(model_shannon_family)



model_shannon_genus <-
  lmer(shannon_genus ~ CollectionPhase +(1|Worker), data = expanded_metadata_edit, REML=T)
summary(model_shannon_genus)
anova(model_shannon_genus, type='III',test="F")

lsmeans(model_shannon_genus, pairwise~CollectionPhase, adjust="tukey")
VarCorr(model_shannon_genus)


model_shannon_asv <-
  lmer(shannon_asv ~ CollectionPhase +(1|Worker), data = expanded_metadata_edit, REML=T)
summary(model_shannon_asv)
anova(model_shannon_asv, type='III',test="F")

lsmeans(model_shannon_asv, pairwise~CollectionPhase, adjust="tukey")
VarCorr(model_shannon_asv)

ggplot(expanded_metadata_edit, aes(CollectionPhase, Observed_asv))+geom_boxplot()


## model after--removing swine and enviroment:

expanded_metadata_edit <- subset (expanded_metadata, !CollectionPhase%in% c("NegCtrl", "MockComm", "SWINE", "ENVIRONMENT"))
str(expanded_metadata_edit)

##Observed_phylum 
##shannon_phylum 

##Observed_family 
##shannon_family

## Observed_genus 
## shannon_genus 


model_obs_phy <-
  lmer(Observed_phylum ~ CollectionPhase +(1|Worker), data = expanded_metadata_edit, REML=T)
summary(model_obs_phy)
anova(model_obs_phy, type='III',test="F")

lsmeans(model_obs_phy, pairwise~CollectionPhase, adjust="tukey")
VarCorr(model_obs_phy)


model_obs_family <-
  lmer(Observed_family ~ CollectionPhase +(1|Worker), data = expanded_metadata_edit, REML=T)
summary(model_obs_family)
anova(model_obs_family, type='III',test="F")

lsmeans(model_obs_family, pairwise~CollectionPhase, adjust="tukey")
VarCorr(model_obs_family)



model_obs_genus <-
  lmer(Observed_genus ~ CollectionPhase +(1|Worker), data = expanded_metadata_edit, REML=T)
summary(model_obs_genus)
anova(model_obs_genus, type='III',test="F")

lsmeans(model_obs_genus, pairwise~CollectionPhase, adjust="tukey")
VarCorr(model_obs_genus)


model_obs_asv <-
  lmer(Observed_asv ~ CollectionPhase +(1|Worker), data = expanded_metadata_edit, REML=T)
summary(model_obs_asv)
anova(model_obs_asv, type='III',test="F")

lsmeans(model_obs_asv, pairwise~CollectionPhase, adjust="tukey")
VarCorr(model_obs_asv)


## shannon:
model_shannon_phy <-
  lmer(shannon_phylum ~ CollectionPhase +(1|Worker), data = expanded_metadata_edit, REML=T)
summary(model_shannon_phy)
anova(model_shannon_phy, type='III',test="F")

lsmeans(model_shannon_phy, pairwise~CollectionPhase, adjust="tukey")
VarCorr(model_shannon_phy)


model_shannon_family <-
  lmer(shannon_family ~ CollectionPhase +(1|Worker), data = expanded_metadata_edit, REML=T)
summary(model_shannon_family)
anova(model_shannon_family, type='III',test="F")

lsmeans(model_shannon_family, pairwise~CollectionPhase, adjust="tukey")
VarCorr(model_shannon_family)



model_shannon_genus <-
  lmer(shannon_genus ~ CollectionPhase +(1|Worker), data = expanded_metadata_edit_occphase, REML=T)
summary(model_shannon_genus)
anova(model_shannon_genus, type='III',test="F")

lsmeans(model_shannon_genus, pairwise~CollectionPhase, adjust="BH")
VarCorr(model_shannon_genus)


model_shannon_asv <-
  lmer(shannon_asv ~ CollectionPhase +(1|Worker), data = expanded_metadata_edit, REML=T)
summary(model_shannon_asv)
anova(model_shannon_asv, type='III',test="F")

lsmeans(model_shannon_asv, pairwise~CollectionPhase, adjust="tukey")
VarCorr(model_shannon_asv)



################################################################################
#                           BETA DIVERSITY                                     #      
################################################################################


ps.noncontam.worker <- readRDS("ilya_worker_microbiome_phyloseq_noncontam_7_21_22")
microbiome.data.file<-otu_table(ps.noncontam.worker)
microbiome.data.file.counts<-data.frame(microbiome.data.file)
microbiome.data.file.counts
write.csv(microbiome.data.file.counts, "Supplementary_datafile_2_microbiome.data.file.counts.csv")
microbiome_sample.data.file<-sample_data(ps.noncontam.worker)
microbiome.data.file.samples<-data.frame(microbiome_sample.data.file)
write.csv(microbiome.data.file.samples, "Supplementary_datafile_3_microbiome.data.file.counts.csv")
microbiome_taxa.data.file<-tax_table(ps.noncontam.worker)
microbiome.data.file.samples<-data.frame(microbiome_sample.data.file)
write.csv(microbiome.data.file.samples, "Supplementary_datafile_3_microbiome.data.file.counts.csv")

ps.noncontam.worker

#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 9358 taxa and 47 samples ]
#sample_data() Sample Data:       [ 47 samples by 90 sample variables ]
#tax_table()   Taxonomy Table:    [ 9358 taxa by 7 taxonomic ranks ]
#refseq()      DNAStringSet:      [ 9358 reference sequences ]

biome_phylum.ps <-  tax_glom(ps.noncontam.worker, "Phylum")
biome_class.ps <-   tax_glom(ps.noncontam.worker, "Class")
biome_family.ps <-  tax_glom(ps.noncontam.worker, "Family")
biome_genus.ps <-   tax_glom(ps.noncontam.worker, "Genus")
biome_species.ps <- tax_glom(ps.noncontam.worker, "Species")
biome_asv.ps <- ps.noncontam.worker

########################################
##########
#ANAYLYSIS OF GENUS LEVEL BETA DIVERSITY
########################################
###########

 
#Column(s) containing all zeros/unobserved values were found (check it out using zPatterns).
#So let us remove environmental samples
sum(taxa_sums(biome_genus.ps)==0) # taxa not present across all samples
#[1] 30
biome_genus.ps
biome_genus.ps.pruned <- prune_taxa(taxa_sums(biome_genus.ps) > 0.0000001, biome_genus.ps)#prunning low-abundance taxa of gene groups
sum(taxa_sums(biome_genus.ps.pruned)==0)
biome_genus.ps.pruned
ALLGENUSCounts <- as.data.frame(phyloseq::otu_table(biome_genus.ps.pruned))
if(phyloseq::taxa_are_rows(biome_genus.ps.pruned) == FALSE){ otus <- t(ALLASVCounts) }
ALLGENUSCounts

#zCompositions component of the anlaysis, using bayesian approach:


ALLGENUS_cmultRepl<-cmultRepl(ALLGENUSCounts,  label = 0, method = c("CZM"), output = c("p-counts"), frac = 0.2, threshold = 0.5, adjust = TRUE, t= NULL, s= NULL, z.warning = 0.8, suppress.print = FALSE)
ALLGENUS_cmultRepl

zPatterns(ALLARG_cmultRepl)

#Now we must reconvert this OTU table to the phyloseq object

biome_genus.ps.pruned <-phyloseq::otu_table(biome_genus.ps.pruned) <- phyloseq::otu_table(ALLGENUS_cmultRepl, taxa_are_rows = F)


SampData <- as(sample_data(biome_genus.ps), 'data.frame')
TAXData <- as(tax_table(biome_genus.ps), 'matrix')

#Create a new phyloseq object to include the imputed virus OTU table

IMPGENUS_MICROBIOME = otu_table(biome_genus.ps.pruned, taxa_are_rows = F)
TAX = tax_table(TAXData)
MICROBIOME_SAMPLES = sample_data(SampData)
IMPUTED_GENUS_MICROBIOME_physeq <- phyloseq(IMPGENUS_MICROBIOME, TAX, MICROBIOME_SAMPLES)
saveRDS(IMPUTED_GENUS_MICROBIOME_physeq, "pesky_code_plz_help.rds")
TAX

#Now let's retry the ordination w/ this new phyloseq object

sum(taxa_sums(IMPUTED_GENUS_MICROBIOME_physeq)==0) # 0 taxa not present across all samples
#IMPUTED_GENUS<- prune_taxa(taxa_sums(IMPUTED_GENUS_MICROBIOME_physeqSUB) > 0, IMPUTED_GENUS_MICROBIOME_physeqSUB)#prunning low-abundance taxa of gene groups

GENUS_physeq.css <- phyloseq_transform_css(IMPUTED_GENUS_MICROBIOME_physeqSUB)

GENUS_physeq.dist <- vegdist(decostand(t(otu_table(IMPUTED_GENUS_MICROBIOME_physeqSUB)),"rclr"), method = "euclidean")

set.seed(1999)
GENUS_physeq.dist
GENUS_physeq.ord<- ordinate(GENUS_physeq.css, method = "PCoA", distance = "euclidean")
GENUS_physeq.ord
GENUS_plot<-plot_ordination(IMPUTED_GENUS_MICROBIOME_physeqSUB,GENUS_physeq.ord, type = "samples", color="CollectionPhase")
str(GENUS_plot)

GENUScentroid<- GENUS_plot$data %>% group_by(CollectionPhase) %>% summarise(Axis.1=mean(Axis.1), Axis.2=mean(Axis.2), .groups = "drop")

#stressplot(plot_ord_all)
#str(plot_ord_all)

GENUS_ordplot <- ggplot(GENUS_plot$data,GENUS_plot$mapping) + stat_ellipse(geom = "polygon", type= "t", level = 0.85, lty =1, lwd = 1.3, alpha= 0.07, show.legend =FALSE) + theme_bw() +geom_point(aes(fill=CollectionPhase), size=6, alpha=0.3, stroke=1.0,) + geom_point(data=GENUScentroid, size=13, shape=16, fill="black",mapping= aes(colour=CollectionPhase), stroke=4, alpha= 0.8, show.legend = FALSE) + scale_color_manual(name="Collection phase", breaks = c("WORKDAY_START", "WORKDAY_END", "POST_SHOWER", "SWINE", "ENVIRONMENT"), values = c("#789ED0", "#F79D40", "#72BF56", "#E4655A", "#AB8BCB"), labels= c("WORKDAY_START"= "Workday start", "WORKDAY_END"="Workday end", "POST_SHOWER"="Post-shower", "SWINE"="Swine", "ENVIRONMENT"="Environment")) +
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
        axis.title.x = element_text(size = 18, vjust = -1.5)) + coord_fixed(xlim=c(-35,40), ylim = c(-35,40))  
GENUS_ordplot



ps.clr <- microbiome::transform(IMPUTED_GENUS_MICROBIOME_physeqSUB, "clr")
ps.clr

otu <- otu_table(ps.clr) %>% as.matrix()
meta <- sample_data(ps.clr) %>% data.frame()

ret <- vegan::anosim(otu, meta$CollectionPhase, distance = "euclidean", perm=1000)
ret
#ANOSIM statistic R: 0.4587 
#Significance: 0.000999 

d <- vegan::vegdist(otu, method = "euclidean")
str(d)
perm <- pairwise.adonis(d, meta$CollectionPhase, sim.method= 'euclidean', perm = 1000, p.adjust.m = "BH")

perm

#                           pairs Df SumsOfSqs  F.Model R2        p.value   p.adjusted sig
#1          ENVIRONMENT vs SWINE  1  5331.277 4.209700 0.2962554 0.015984016 0.017760018   .
#2  ENVIRONMENT vs WORKDAY_START  1  4378.171 2.452841 0.1969704 0.008991009 0.014985015   .
#3    ENVIRONMENT vs WORKDAY_END  1  4771.348 4.243250 0.2979130 0.015984016 0.017760018   .
#4    ENVIRONMENT vs POST_SHOWER  1  4558.650 3.227733 0.2440126 0.015984016 0.017760018   .
#5        SWINE vs WORKDAY_START  1  6554.868 4.635528 0.2047899 0.000999001 0.001998002   *
#6          SWINE vs WORKDAY_END  1  3477.533 3.321061 0.1557643 0.000999001 0.001998002   *
#7          SWINE vs POST_SHOWER  1  6030.596 4.996147 0.2172602 0.000999001 0.001998002   *
#8  WORKDAY_START vs WORKDAY_END  1  4407.020 3.300700 0.1549574 0.000999001 0.001998002   *
#9  WORKDAY_START vs POST_SHOWER  1  1614.516 1.079863 0.0565970 0.247752248 0.247752248    
#10   WORKDAY_END vs POST_SHOWER  1  3309.925 2.933871 0.1401495 0.000999001 0.001998002   * 



mod <- betadisper(d, meta$CollectionPhase,type = c("centroid"), bias.adjust=TRUE)
mod
anovamod <- anova(betadisper(d, meta$CollectionPhase,type = c("centroid"), bias.adjust=TRUE))
anovamod
#Analysis of Variance Table
#Response: Distances
#Df  Sum Sq Mean Sq F value    Pr(>F)    
#Groups     4  943.53 235.882  7.0516 0.0002522 ***
#  Residuals 37 1237.68  33.451    

boxplot(mod,medlwd=2,frame.plot = F,ylab="Distance to centroid", xlab="Collection phase",
        cex.axis = 0.5,
        col=(c("WORKDAY_START"="grey50", "WORKDAY_END" = "red", "POST_SHOWER" = "#59A14F",
               "SWINE" = "blue", "ENVIRONMENT" = "purple")))

anovamod.permdisp<- permutest(mod, permutations = 1000, pairwise=T, p.adjust.m="FDR")
anovamod.permdisp
#Swine samples relative to workday start samples had significant dispersion at 0.029.

box(bty="l")
axis(2)
text(2,2,family="A")
windowsFonts(A=windowsFont("Times New Roman"))
topptx(filename = "betadisper_genus.pptx",width=6, height=5)

#Procrustes analysis between MICROBIOME (Genus) and RESISTOME (Group)
proTest_microbiome_resistome<- protest(d,ALLARGS_physeq.dist.SUB , scores= "sites", symmetric = TRUE, scale=TRUE, type=segments(CollectionPhase), permutations = (nperm=999), choices=1:2)
proTest_microbiome_resistome
plot(proTest_microbiome_resistome, main="AMR group vs Microbiome Genus Procrustes")
text(0, 0.4, "M^2= 0.6866, p-value = 0.6",cex = 1)
#Procrustes Sum of Squares (m12 squared):        0.6866 
#Correlation in a symmetric Procrustes rotation: 0.5598 
#Significance:  0.6 
#Permutation: free
#Number of permutations: 999


#Procrustes analysis between MICROBIOME (Genus) and MED_IMP_RESISTOME (Group)
proTest_microbiome_MEDIMPresistome<- protest(d,IMPARGS_physeq.dist.SUB , scores= "sites", symmetric = TRUE, scale=TRUE, type=segments(CollectionPhase), permutations = (nperm=999), choices=1:2)
proTest_microbiome_MEDIMPresistome
plot(proTest_microbiome_MEDIMPresistome, main="Medically important ARG group vs Microbiome Genus Procrustes")
text(0, 0.4, "M^2= 0.6866, p-value = 0.6",cex = 1)
ade4::procruste(d, IMPARGS_physeq.dist.SUB, scale = TRUE) 
#Call:
#  protest(X = d, Y = IMPARGS_physeq.dist.SUB, scores = "sites",      symmetric = TRUE) 

#Procrustes Sum of Squares (m12 squared):        0.6866 
#Correlation in a symmetric Procrustes rotation: 0.5598 
#Significance:  0.607 

#Permutation: free
#Number of permutations: 999


#Procrustes analysis between MICROBIOME (Genus) and MOBILOME (Plasmid)
proTest_microbiome_resistome<- protest(d,PLASMID_physeq.ps.pruned.dist.SUB , scores= "sites", symmetric = TRUE, scale=TRUE, type=segments(CollectionPhase), permutations = (nperm=999), choices=1:2, display="rotated")
proTest_microbiome_resistome
plot(proTest_microbiome_resistome, main="AMR group vs Microbiome Genus Procrustes", 
     xlim=c(-0.17,0.17))
text(0, 0.4, "M^2= 0.6866, p-value = 0.6",cex = 1)
#Call:
#  protest(X = d, Y = PLASMID_physeq.ps.pruned.dist.SUB, scores = "sites",      symmetric = TRUE) 

#Procrustes Sum of Squares (m12 squared):        0.5856 
#Correlation in a symmetric Procrustes rotation: 0.6437 
#Significance:  0.251 

#Permutation: free
#Number of permutations: 999

#Procrustes analysis between MICROBIOME (Genus) and MOBILOME (Prophage)
proTest_microbiome_resistome<- protest(d,PROPHAGE_physeq.ps.pruned.dist.SUB , scores= "sites", symmetric = TRUE, type=segments(CollectionPhase), permutations = (nperm=999), choices=1:2)
proTest_microbiome_resistome
plot(proTest_microbiome_resistome, main="AMR group vs Microbiome Genus Procrustes", 
     xlim=c(-0.17,0.17))
text(0, 0.4, "M^2= 0.6866, p-value = 0.6",cex = 1)
 

#Call:
#  protest(X = d, Y = PROPHAGE_physeq.ps.pruned.dist.SUB, scores = "sites",      symmetric = TRUE) 
#Procrustes Sum of Squares (m12 squared):        0.7773 
#Correlation in a symmetric Procrustes rotation: 0.4719 
#Significance:  0.586 

#Procrustes analysis between RESISTOME (GROUP) and MOBILOME (ALLMOBILOME GENES)
proTest_microbiome_resistome<- protest(ALLARGS_physeq.dist.SUB,PLASMID_physeq.ps.pruned.dist.SUB , scores= "sites", symmetric = TRUE, scale=TRUE, type=segments(CollectionPhase), permutations = (nperm=999), choices=1:2)
proTest_microbiome_resistome
plot(proTest_microbiome_resistome, main="AMR group vs Microbiome Genus Procrustes", 
     xlim=c(-0.2,0.2), ylim=c(-0.28,0.25), kind=1)
text(0, 0.4, "M^2= 0.6866, p-value = 0.6",cex = 1)




############################################################################################
#DIFFERENTIAL ABUNDANCE TESTING#############################################################
############################################################################################


#Genus level analysis

sample_data(IMPUTED_GENUS_MICROBIOME_physeq)$CollectionPhase <- as.factor(sample_data(IMPUTED_GENUS_MICROBIOME_physeq)$CollectionPhase) # factorize for DESeq2
sample_data(IMPUTED_GENUS_MICROBIOME_physeq)$EVSMOK <- as.factor(sample_data(IMPUTED_GENUS_MICROBIOME_physeq)$EVSMOK) # factorize for DESeq2
sample_data(IMPUTED_GENUS_MICROBIOME_physeq)$DADAReadPairNonchim  <- as.numeric(sample_data(IMPUTED_GENUS_MICROBIOME_physeq)$DADAReadPairNonchim ) # factorize for DESeq2
sample_data(IMPUTED_GENUS_MICROBIOME_physeq)$DADAReadPairNonchim  <- scale(sample_data(IMPUTED_GENUS_MICROBIOME_physeq)$DADAReadPairNonchim , center = TRUE, scale=TRUE) # factorize for DESeq2
sample_data(IMPUTED_GENUS_MICROBIOME_physeq)$qpcr_16s_copies_ul    <- as.numeric(sample_data(IMPUTED_GENUS_MICROBIOME_physeq)$qpcr_16s_copies_ul) # factorize for DESeq2
sample_data(IMPUTED_GENUS_MICROBIOME_physeq)$qpcr_16s_copies_ul  <- scale(sample_data(IMPUTED_GENUS_MICROBIOME_physeq)$qpcr_16s_copies_ul , center = TRUE, scale=TRUE) # factorize for DESeq2
sample_data(IMPUTED_GENUS_MICROBIOME_physeq)$GEND   <- as.factor(sample_data(IMPUTED_GENUS_MICROBIOME_physeq)$GEND) # factorize for DESeq2
sample_data(IMPUTED_GENUS_MICROBIOME_physeq)$AGE   <- scale(sample_data(IMPUTED_GENUS_MICROBIOME_physeq)$AGE  , center = TRUE, scale=TRUE) # factorize for DESeq2
sample_data(IMPUTED_GENUS_MICROBIOME_physeq)$BMI   <- scale(sample_data(IMPUTED_GENUS_MICROBIOME_physeq)$BMI  , center = TRUE, scale=TRUE) # factorize for DESeq2

BiocManager::install('EnhancedVolcano')
library('EnhancedVolcano')
##############SUBSET FOR WORKSTART_WORK_END

IMPUTED_GENUS_MICROBIOME_physeqWSWE <- subset_samples(IMPUTED_GENUS_MICROBIOME_physeq, CollectionPhase=="WORKDAY_START" | CollectionPhase=="WORKDAY_END")

#Filter sparse features of IMPUTED_ALLRESISTOME_physeq object so that we only consider those features with less than 90% zeros
#sample_data(IMPUTED_ALLRESISTOME_physeq)$CollectionPhase<- as.factor(sample_data(IMPUTED_ALLRESISTOME_physeq))
IMPUTED_GENUS_MICROBIOME_physeqWSWE
#nosparse_GENUS_MICROBIOME<- prune_taxa(rowSums(otu_table(IMPUTED_GENUS_MICROBIOME_physeqWSWE)==0) < ncol(otu_table(IMPUTED_GENUS_MICROBIOME_physeqWSWE))*0.9, IMPUTED_GENUS_MICROBIOME_physeqWSWE)

DS_GENUSMICROBIOME= phyloseq_to_deseq2(IMPUTED_GENUS_MICROBIOME_physeqWSWE,design= ~CollectionPhase+ qpcr_16s_copies_ul + DADAReadPairNonchim + GEND + BMI + EVSMOK)
DS_ALLGENUS<- estimateSizeFactors(DS_GENUSMICROBIOME, type="poscounts")
DS_ALLGENUS$CollectionPhase <- relevel(DS_ALLGENUS$CollectionPhase, ref = "WORKDAY_START")
DS_ALLGENUS= DESeq(DS_ALLGENUS, test= "Wald", fitType= "parametric") 
resultsNames(DS_ALLGENUS)
alpha= 0.01
DS_ALLGENUS_RESULT= results(DS_ALLGENUS, alpha = alpha)
DS_ALLGENUS <- lfcShrink(DS_ALLGENUS, res= DS_ALLGENUS_RESULT, type="ashr")
DS_ALLGENUS_RESULT= DS_ALLGENUS_RESULT[order(DS_ALLGENUS_RESULT$padj, na.last = NA), ]
DS_ALLGENUS_RESULT
summary(DS_ALLGENUS_RESULT)
R.microbiome<-as.data.frame(DS_ALLGENUS_RESULT)
R_SIG_microbiome<- head(R.microbiome[, 1:6], 16)
R_SIG_microbiome<- R_SIG_microbiome[order(R_SIG_microbiome$baseMean, na.last = NA), ]
write.csv(R_SIG_microbiome, "C:/Users/is233/OneDrive/Desktop/AG_WORKER_METAGENOMICS_STUDY/MANUSCRIPT/ANALYSIS/WORKER_MICROBIOME/DESEQ_GENUS_WSWE.csv", row.names = TRUE)
DS_ALLGENUS <- estimateSizeFactors(DS_ALLGENUS)
DS_ALLGENUS <- estimateDispersions(DS_ALLGENUS)
plotDispEsts(DS_ALLGENUS)
par(mfrow=c(1,1))
# Make a basic volcano plot
with(R.microbiome, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot", xlim=c(-50,50)))
with(subset(R.microbiome, padj<.05 ), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
with(subset(R.microbiome, padj<.01 & abs(log2FoldChange)>2), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))
EnhancedVolcano(R.microbiome,
                lab = rownames(R.microbiome),
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

##############SUBSET FOR WORKSEND_POST_SHOWER

IMPUTED_GENUS_MICROBIOME_physeqWEPS <- subset_samples(IMPUTED_GENUS_MICROBIOME_physeq, CollectionPhase=="WORKDAY_END" | CollectionPhase=="POST_SHOWER")

#Filter sparse features of IMPUTED_ALLRESISTOME_physeq object so that we only consider those features with less than 90% zeros
#sample_data(IMPUTED_ALLRESISTOME_physeq)$CollectionPhase<- as.factor(sample_data(IMPUTED_ALLRESISTOME_physeq))
IMPUTED_GENUS_MICROBIOME_physeqWEPS
#nosparse_GENUS_MICROBIOME<- prune_taxa(rowSums(otu_table(IMPUTED_GENUS_MICROBIOME_physeqWSWE)==0) < ncol(otu_table(IMPUTED_GENUS_MICROBIOME_physeqWSWE))*0.9, IMPUTED_GENUS_MICROBIOME_physeqWSWE)

DS_GENUSMICROBIOME= phyloseq_to_deseq2(IMPUTED_GENUS_MICROBIOME_physeqWEPS,design= ~CollectionPhase+ qpcr_16s_copies_ul + DADAReadPairNonchim + GEND + BMI + EVSMOK)
DS_ALLGENUS<- estimateSizeFactors(DS_GENUSMICROBIOME, type="poscounts")
DS_ALLGENUS$CollectionPhase <- relevel(DS_ALLGENUS$CollectionPhase, ref = "WORKDAY_END")
DS_ALLGENUS= DESeq(DS_ALLGENUS, test= "Wald", fitType= "parametric") 
resultsNames(DS_ALLGENUS)
alpha= 0.01
DS_ALLGENUS_RESULT= results(DS_ALLGENUS, alpha = alpha)
DS_ALLGENUS <- lfcShrink(DS_ALLGENUS, res= DS_ALLGENUS_RESULT, type="ashr")
DS_ALLGENUS_RESULT= DS_ALLGENUS_RESULT[order(DS_ALLGENUS_RESULT$padj, na.last = NA), ]
DS_ALLGENUS_RESULT
summary(DS_ALLGENUS_RESULT)
R.microbiome<-as.data.frame(DS_ALLGENUS_RESULT)
R_SIG_microbiome<- head(R.microbiome[, 1:6], 15)
R_SIG_microbiome<- R_SIG_microbiome[order(R_SIG_microbiome$baseMean, na.last = NA), ]
write.csv(R_SIG_microbiome, "C:/Users/is233/OneDrive/Desktop/AG_WORKER_METAGENOMICS_STUDY/MANUSCRIPT/ANALYSIS/WORKER_MICROBIOME/DESEQ_GENUS_WEPS.csv", row.names = TRUE)
DS_ALLGENUS <- estimateSizeFactors(DS_ALLGENUS)
DS_ALLGENUS <- estimateDispersions(DS_ALLGENUS)
plotDispEsts(DS_ALLGENUS)
par(mfrow=c(1,1))
# Make a basic volcano plot
with(R.microbiome, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot", xlim=c(-50,50)))
with(subset(R.microbiome, padj<.05 ), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
with(subset(R.microbiome, padj<.01 & abs(log2FoldChange)>2), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))
EnhancedVolcano(R.microbiome,
                lab = rownames(R.microbiome),
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



##############SUBSET FOR POST_SHOWER_WORKSTART

IMPUTED_GENUS_MICROBIOME_physeqPSWS <- subset_samples(IMPUTED_GENUS_MICROBIOME_physeq, CollectionPhase=="WORKDAY_START" | CollectionPhase=="POST_SHOWER")

#Filter sparse features of IMPUTED_ALLRESISTOME_physeq object so that we only consider those features with less than 90% zeros
#sample_data(IMPUTED_ALLRESISTOME_physeq)$CollectionPhase<- as.factor(sample_data(IMPUTED_ALLRESISTOME_physeq))
IMPUTED_GENUS_MICROBIOME_physeqPSWS
#nosparse_GENUS_MICROBIOME<- prune_taxa(rowSums(otu_table(IMPUTED_GENUS_MICROBIOME_physeqWSWE)==0) < ncol(otu_table(IMPUTED_GENUS_MICROBIOME_physeqWSWE))*0.9, IMPUTED_GENUS_MICROBIOME_physeqWSWE)

DS_GENUSMICROBIOME= phyloseq_to_deseq2(IMPUTED_GENUS_MICROBIOME_physeqPSWS,design= ~CollectionPhase+ qpcr_16s_copies_ul + DADAReadPairNonchim + GEND + BMI + EVSMOK)
DS_ALLGENUS<- estimateSizeFactors(DS_GENUSMICROBIOME, type="poscounts")
DS_ALLGENUS$CollectionPhase <- relevel(DS_ALLGENUS$CollectionPhase, ref = "WORKDAY_START")
DS_ALLGENUS= DESeq(DS_ALLGENUS, test= "Wald", fitType= "parametric") 
resultsNames(DS_ALLGENUS)
alpha= 0.01
DS_ALLGENUS_RESULT= results(DS_ALLGENUS, alpha = alpha)
DS_ALLGENUS <- lfcShrink(DS_ALLGENUS, res= DS_ALLGENUS_RESULT, type="ashr")
DS_ALLGENUS_RESULT= DS_ALLGENUS_RESULT[order(DS_ALLGENUS_RESULT$padj, na.last = NA), ]
DS_ALLGENUS_RESULT
summary(DS_ALLGENUS_RESULT)
R.microbiome<-as.data.frame(DS_ALLGENUS_RESULT)
R_SIG_microbiome<- head(R.microbiome[, 1:6], 28)
R_SIG_microbiome<- R_SIG_microbiome[order(R_SIG_microbiome$baseMean, na.last = NA), ]
write.csv(R_SIG_microbiome, "C:/Users/is233/OneDrive/Desktop/AG_WORKER_METAGENOMICS_STUDY/MANUSCRIPT/ANALYSIS/WORKER_MICROBIOME/DESEQ_GENUS_PSWS.csv", row.names = TRUE)
DS_ALLGENUS <- estimateSizeFactors(DS_ALLGENUS)
DS_ALLGENUS <- estimateDispersions(DS_ALLGENUS)
plotDispEsts(DS_ALLGENUS)
par(mfrow=c(1,1))
# Make a basic volcano plot
with(R.microbiome, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot", xlim=c(-50,50)))
with(subset(R.microbiome, padj<.05 ), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
with(subset(R.microbiome, padj<.01 & abs(log2FoldChange)>2), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))
EnhancedVolcano(R.microbiome,
                lab = rownames(R.microbiome),
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

##############SUBSET FOR SWINE_WORKEND

IMPUTED_GENUS_MICROBIOME_physeqSWINE_WE <- subset_samples(IMPUTED_GENUS_MICROBIOME_physeq, CollectionPhase=="WORKDAY_END" | CollectionPhase=="SWINE")

#Filter sparse features of IMPUTED_ALLRESISTOME_physeq object so that we only consider those features with less than 90% zeros
#sample_data(IMPUTED_ALLRESISTOME_physeq)$CollectionPhase<- as.factor(sample_data(IMPUTED_ALLRESISTOME_physeq))
IMPUTED_GENUS_MICROBIOME_physeqSWINE_WE
#nosparse_GENUS_MICROBIOME<- prune_taxa(rowSums(otu_table(IMPUTED_GENUS_MICROBIOME_physeqWSWE)==0) < ncol(otu_table(IMPUTED_GENUS_MICROBIOME_physeqWSWE))*0.9, IMPUTED_GENUS_MICROBIOME_physeqWSWE)

DS_GENUSMICROBIOME= phyloseq_to_deseq2(IMPUTED_GENUS_MICROBIOME_physeqSWINE_WE,design= ~CollectionPhase+ qpcr_16s_copies_ul + DADAReadPairNonchim)
DS_ALLGENUS<- estimateSizeFactors(DS_GENUSMICROBIOME, type="poscounts")
DS_ALLGENUS$CollectionPhase <- relevel(DS_ALLGENUS$CollectionPhase, ref = "SWINE")
DS_ALLGENUS= DESeq(DS_ALLGENUS, test= "Wald", fitType= "parametric") 
resultsNames(DS_ALLGENUS)
alpha= 0.01
DS_ALLGENUS_RESULT= results(DS_ALLGENUS, alpha = alpha)
DS_ALLGENUS <- lfcShrink(DS_ALLGENUS, res= DS_ALLGENUS_RESULT, type="ashr")
DS_ALLGENUS_RESULT= DS_ALLGENUS_RESULT[order(DS_ALLGENUS_RESULT$padj, na.last = NA), ]
DS_ALLGENUS_RESULT
summary(DS_ALLGENUS_RESULT)
R.microbiome<-as.data.frame(DS_ALLGENUS_RESULT)
R_SIG_microbiome<- head(R.microbiome[, 1:6], 1)
R_SIG_microbiome<- R_SIG_microbiome[order(R_SIG_microbiome$baseMean, na.last = NA), ]
write.csv(R_SIG_microbiome, "C:/Users/is233/OneDrive/Desktop/AG_WORKER_METAGENOMICS_STUDY/MANUSCRIPT/ANALYSIS/WORKER_MICROBIOME/DESEQ_GENUS_SWINE_WE.csv", row.names = TRUE)
DS_ALLGENUS <- estimateSizeFactors(DS_ALLGENUS)
DS_ALLGENUS <- estimateDispersions(DS_ALLGENUS)
plotDispEsts(DS_ALLGENUS)
par(mfrow=c(1,1))
# Make a basic volcano plot
with(R.microbiome, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot", xlim=c(-50,50)))
with(subset(R.microbiome, padj<.05 ), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
with(subset(R.microbiome, padj<.01 & abs(log2FoldChange)>2), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))
EnhancedVolcano(R.microbiome,
                lab = rownames(R.microbiome),
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
