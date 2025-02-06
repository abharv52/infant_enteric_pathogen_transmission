#########################################
# Script to create figure 2
library(tidyverse)
library(gridExtra)
library(cowplot)
library(patchwork)

# Read in raw data
TAC <- read.csv('pathogen_detection_samples.csv')

# pathogen hits
TAC_hits <- TAC %>% group_by(type, pathogens) %>%
  summarize(N_pos = sum(PosNeg=="Positive"),
            N_neg = sum(PosNeg=="Negative"),
            n=n()) %>%
  mutate(frac_pos = N_pos/n)
#### Plot results

ggplot(TAC_hits,aes(x=type,y=pathogens,z=frac_pos)) +
  stat_summary_2d() +xlab('') + ylab('')+
  theme_bw() +
  scale_fill_gradient(low='#FFFFFF',
                      high='#004466',
                      limits=c(0,1)) +
  theme(axis.ticks.x=element_blank(),
        plot.margin = unit(c(0, 0, 0, 0), "cm"),
        legend.position='none')

#### Merge with household data to obtain child age data ####
cstool <- filter(TAC,type=="CSTOOL")

# Obtaining HH barcode numbers from the sample code
cstool$barcode_hh <- gsub("CSTOOL","HH",cstool$sample_barcode)
cstool$barcode_hh <- gsub("R2","R1",cstool$barcode_hh)
cstool$barcode_hh <- gsub("R3","R1",cstool$barcode_hh)
cstool$barcode_hh <- gsub("R4","R1",cstool$barcode_hh)
cstool$barcode_hh <- gsub("RE","MN",cstool$barcode_hh)

# merge with HH data
HH_survey <- read.csv("main_household_survey.csv")
HH_survey <- dplyr::select(HH_survey,barcode_hh,child_age)

cstoolHH <-merge(cstool,HH_survey,by="barcode_hh")

cstoolHH$agegroup <- factor(ifelse(cstoolHH$child_age<3,'0-2 months',
                                   ifelse(cstoolHH$child_age>=3&cstoolHH$child_age<6,'3-5 months',
                                          ifelse(cstoolHH$child_age>=6&cstoolHH$child_age<12,'6-11 months',
                                                 ifelse(cstoolHH$child_age>=12,'12-23 months',NA)))),
                            levels=c("0-2 months","3-5 months","6-11 months","12-23 months"))

# pathogen hits
cstool_hits <- cstoolHH %>% group_by(agegroup, pathogens) %>%
  summarize(N_pos = sum(PosNeg=="Positive"),
            N_neg = sum(PosNeg=="Negative"),
            n=n()) %>%
  mutate(frac_pos = N_pos/n)

typeorder_cstool <- c("0-2 months","3-5 months","6-11 months","12-23 months")

ggplot(cstool_hits,aes(x=agegroup,y=pathogens,z=frac_pos)) +
  stat_summary_2d() +xlab('') + ylab('')+
  theme_bw() +
  scale_x_discrete(limits=typeorder_cstool,labels=NULL) +
  scale_fill_gradient(low='#FFFFFF',
                      high='#004466',
                      limits=c(0,1)) +
  theme(legend.position='none',
        axis.ticks=element_blank(),
        plot.margin = unit(c(0, 0, 0, 0), "cm"),
        axis.text.y=element_blank())


############# Creating plots of number of pathogens ####################
# pulling out list of all pathogen targets
pathogens <- unique(TAC_hits$pathogens)
controls <- c("Bacterial 16S","Bacteroides fragilis","Crassphage","Total 18S")
pathogens2 <- pathogens[!pathogens%in%controls] #removing controls from list

TAC_pathogens <- TAC[order(TAC$Cardno,TAC$Sample,TAC$pathogens),]
TAC_pathogens <- filter(TAC_pathogens,pathogens%in%pathogens2)

TAC_pathogens <- TAC_pathogens %>%
  group_by(type,Sample) %>%
  summarize(npos=sum(PosNeg=="Positive")) %>%
  filter(type!="NTC"&type!="EX BLA")

TAC_pathogens$type <- factor(TAC_pathogens$type,levels=c("CSTOOL","CHHAND","MSTOOL","HHSOIL","CHFOOD","DRIWAT",
                                                         "CHCKNF","DOGSFS","CAMELF","CATTLE","GOATSF","SHEEPF"))

ggplot(TAC_pathogens,aes(x=type,y=npos,fill=type)) + geom_boxplot() +
  theme_bw() +
  xlab('') + ylab('Number of pathogens') + ylim(0,12.5) +
  theme(legend.position='none')

##### Compare number of pathogens by child age group ########
# remove non pathogens
cstoolHH3 <- filter(cstoolHH,!pathogens%in%controls)
cstoolHH2 <- cstoolHH3 %>%
  group_by(Sample) %>%
  summarize(npos=sum(PosNeg=="Positive"),
            type=first(type),
            agegroup=first(agegroup))

ggplot(cstoolHH2,aes(x=agegroup,y=npos)) + geom_boxplot(fill="#F8766D") +
  theme_bw() +
  xlab('') + ylab('Number of pathogens') + ylim(0,12.5) + ylab('')+
  geom_text(data=labelsfig,aes(x=agegroup,y=11,label=n)) +
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank())
