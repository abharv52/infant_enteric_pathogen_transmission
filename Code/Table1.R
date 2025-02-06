#### Script to examine risk factors for child pathogen carriage #####
#AHP 6/21/23

setwd('/Users/appaulo/Library/CloudStorage/Box-Box/Kenya_PathogenPathways')
library(tidyverse)
library(sandwich)
library(lmtest)


TAC <- read.csv('TAC/TAC results/full_tac_pathogen_nov7_2024.csv')


pathogens1 <- unique(TAC$pathogens)
pathogens1 <- pathogens1[-c(1,6,8,28)] #remove controls

TAC_pathogens <- filter(TAC,pathogens%in%pathogens1)

# filter to only child stool samples
cstool1 <- filter(TAC_pathogens,type=="CSTOOL")

# pre-process cstool to get desired outcome variables
# first, N pathogens carried
cstool <- cstool1 %>%
  group_by(Sample) %>%
  summarize(npath=sum(PosNeg=="Positive")) %>%
  mutate(inf_anypath=ifelse(npath>0,1,0))

# now, E coli pathogens carried
EC_pathogens <- pathogens1[c(2,5,14,17,19,25,26,33)]

cstool_EC <- cstool1 %>%
  filter(pathogens%in%EC_pathogens) %>%
  group_by(Sample) %>%
  summarize(npathEC=sum(PosNeg=="Positive")) %>%
  mutate(inf_anypathEC=ifelse(npathEC>0,1,0))

# Bacterial targets
# Pulling in full pathogen/assay names
test <- read.csv("TAC/SOPs/TAC card design.csv",row.names=NULL)
# adjusting names
test[46,c(2:3)] <- c("Total 18S","18s")
test$Gene[30] <- '5_UTR'
test$Gene[38] <- "G_18S"
test$Pathogen[25] <- "Total E. coli"
test$Gene[41] <- "TR_18S"
test$Gene[39] <- "EH_18S"
test$Gene[34] <- "CP_18S"
test$Gene[43] <- "AD_ITS2"
test$Gene[42] <- "NA_ITS2"
test[47,c(2:3)] <- c("Crassphage","orf00024")
test[48,c(2:3)] <- c("Bacterial 16S","B_16S")
# test$Gene[47] <- "B_16S"
test$Gene[36] <- "CP_LIB13"
test[45,c(2:3)] <- c("xeno","synethic construct")
test$Gene[31] <- "NGI_ORF1-ORF2"
test$Gene[32] <- "NGII_ORF1-ORF2"

#alphabetical order
test2 <- test[order(test$Gene),c(1:3)]

test3 <- test2 %>%
  group_by(Pathogen) %>%
  sample_n(.,1)
names(test3)[1] <- 'pathtype'

bacterial_targets <- test3 %>% filter(pathtype=="Bacteria") %>% select(Pathogen)
virus_targets <- test3 %>% filter(pathtype=="Viruses") %>% select(Pathogen)
ph_targets <- test3 %>% filter(pathtype=="Protozoa"|pathtype=="Helminths") %>% select(Pathogen)

cstool_bacteria <- cstool1 %>%
  filter(pathogens%in%bacterial_targets$Pathogen) %>%
  group_by(Sample) %>%
  summarize(npath_bact=sum(PosNeg=="Positive")) %>%
  mutate(inf_anypath_bact=ifelse(npath_bact>0,1,0))

cstool_virus <- cstool1 %>%
  filter(pathogens%in%virus_targets$Pathogen) %>%
  group_by(Sample) %>%
  summarize(npath_virus=sum(PosNeg=="Positive")) %>%
  mutate(inf_anypath_virus=ifelse(npath_virus>0,1,0))

cstool_ph <- cstool1 %>%
  filter(pathogens%in%ph_targets$Pathogen) %>%
  group_by(Sample) %>%
  summarize(npath_ph=sum(PosNeg=="Positive")) %>%
  mutate(inf_anypath_ph=ifelse(npath_ph>0,1,0))

common_path <- pathogens1[c(2,5,19,15,1,16,22)]

cstool_common <- cstool1 %>%
  filter(pathogens%in%common_path) %>%
  group_by(Sample,pathogens) %>%
  summarize(npath_common=sum(PosNeg=="Positive")) %>%
  mutate(inf_anypath_common=ifelse(npath_common>0,1,0),
         pathogen_var = paste0("c_",pathogens)) %>%
  select(Sample,pathogen_var,inf_anypath_common)  %>%
  pivot_wider(names_from="pathogen_var",values_from="inf_anypath_common")

# pull together all outcomes variables
cstool_risk <- cbind(cstool,cstool_EC[,c(2:3)],cstool_bacteria[,c(2:3)],
                     cstool_virus[,c(2:3)],cstool_ph[,c(2:3)],cstool_common[,c(2:8)])


# Load HH survey data
HH_survey <- read.csv("surveyctodata/clean_data/main_HH_cleaned.csv")

### PPI Score calculation #####
PPI <- read.csv('/Users/appaulo/Documents/GitHub/KenyaPathogenPathways/PPI_scorecard.csv')

HH_survey$PPI <- 0

HH_survey$resp_educ[HH_survey$resp_school==2] <- 1

HH_survey$PPI <- ifelse(HH_survey$county==1,HH_survey$PPI+2,HH_survey$PPI+3) # county
HH_survey$PPI <- ifelse(HH_survey$resp_educ==1,HH_survey$PPI+0,
                        ifelse(HH_survey$resp_educ==2|HH_survey$highest_educ==3,HH_survey$PPI+7,
                               ifelse(HH_survey$resp_educ==4,HH_survey$PPI+10,
                                      HH_survey$PPI+15))) # respondent education
HH_survey$PPI <- ifelse(HH_survey$highest_educ==1,HH_survey$PPI+0,
                        ifelse(HH_survey$highest_educ==2,HH_survey$PPI+3,
                               ifelse(HH_survey$highest_educ==3,HH_survey$PPI+1,
                                      HH_survey$PPI+5))) #any member education
HH_survey$PPI <- ifelse(HH_survey$acq_bread==1,HH_survey$PPI+10,HH_survey$PPI+0) # Acquired bread
HH_survey$PPI <- ifelse(HH_survey$acq_meat==1,HH_survey$PPI+12,HH_survey$PPI+0) # Acquired bread
HH_survey$PPI <- ifelse(HH_survey$acq_banana==1,HH_survey$PPI+9,HH_survey$PPI+0) # Acquired bread
HH_survey$PPI <- ifelse(HH_survey$towel==1,HH_survey$PPI+10,HH_survey$PPI+0) # Acquired bread
HH_survey$PPI <- ifelse(HH_survey$thermos==1,HH_survey$PPI+9,HH_survey$PPI+0) # Acquired bread
HH_survey$PPI <- ifelse(HH_survey$wall_material==1,HH_survey$PPI+8,
                        ifelse(HH_survey$wall_material==2,HH_survey$PPI+6,
                               HH_survey$PPI+0)) # Wall material
HH_survey$PPI <- ifelse(HH_survey$floor_material==1,HH_survey$PPI+0,
                        HH_survey$PPI+7)
HH_survey$PLprob <- PPI$National.Poverty.Line[HH_survey$PPI+1]

#### Food poverty ##########
HH_survey$FPL <- 0

HH_survey$resp_educ[HH_survey$resp_school==2] <- 1

HH_survey$FPL <- ifelse(HH_survey$county==1,HH_survey$FPL+0,HH_survey$FPL+2) # county
HH_survey$FPL <- ifelse(HH_survey$resp_educ==1,HH_survey$FPL+0,
                        ifelse(HH_survey$resp_educ==2|HH_survey$highest_educ==3,HH_survey$FPL+7,
                               ifelse(HH_survey$resp_educ==4,HH_survey$FPL+9,
                                      HH_survey$FPL+10))) # respondent education
HH_survey$FPL <- ifelse(HH_survey$highest_educ==1,HH_survey$FPL+0,
                        ifelse(HH_survey$highest_educ==2,HH_survey$FPL+5,
                               ifelse(HH_survey$highest_educ==3,HH_survey$FPL+5,
                                      HH_survey$FPL+7))) #any member education
HH_survey$FPL <- ifelse(HH_survey$acq_bread==1,HH_survey$FPL+11,HH_survey$FPL+0) # Acquired bread
HH_survey$FPL <- ifelse(HH_survey$acq_meat==1,HH_survey$FPL+11,HH_survey$FPL+0) # Acquired meat
HH_survey$FPL <- ifelse(HH_survey$acq_banana==1,HH_survey$FPL+11,HH_survey$FPL+0) # Acquired banana
HH_survey$FPL <- ifelse(HH_survey$towel==1,HH_survey$FPL+10,HH_survey$FPL+0) # Owns towel
HH_survey$FPL <- ifelse(HH_survey$thermos==1,HH_survey$FPL+9,HH_survey$FPL+0) # Owns thermos
HH_survey$FPL <- ifelse(HH_survey$wall_material==1,HH_survey$FPL+8,
                        ifelse(HH_survey$wall_material==2,HH_survey$FPL+8,
                               HH_survey$FPL+0)) # Wall material
HH_survey$FPL <- ifelse(HH_survey$floor_material==1,HH_survey$FPL+0,
                        HH_survey$FPL+7)
HH_survey$FPLprob <- PPI$National.Poverty.Line[HH_survey$FPL+1]

# Pre-process exploratory variables
HH_survey$pregnancy_term[HH_survey$pregnancy_term==999] <- NA

HH_survey <- HH_survey %>%
  mutate(BPL = ifelse(PLprob>0.5,1,0), #below poverty line
         BFPL = ifelse(FPLprob>0.5,1,0), #below food poverty line
         breastmilk=ifelse(food_types_1==1,1,0), # breastmilk consumption
         food=ifelse(food_types_3==1|food_types_4==1,1,0), # food consumption
         crawls=ifelse(crawl==1,1,0), # crawling status
         walks=ifelse(walk==1,1,0), # walking status
         resp_primary = ifelse(resp_educ>1,1,0), # Primary respondent education
         child_gender=ifelse(child_gender==1,1,0), # Gender
         pregnancy_ab=ifelse(pregnancy_ab==1,1,0), # Antibiotic consumption during pregnancy
         delivery_method=ifelse(delivery_method==1,1,0), # Birth delivery method
         latrine=ifelse(latrine==1,1,0), # Latrine ownership
         animal_feces = ifelse(animal_feces==1,3,ifelse(animal_feces==2|animal_feces==3,1,ifelse(animal_feces==4,2,NA))), # Animal feces disposal
         hh_own_livestock=ifelse(hh_own_livestock==1,1,0), # Livestock ownership
         own_cattle=ifelse(hh_cattle>0,1,0), # Cattle ownership
         own_sheep=ifelse(hh_sheep>0,1,0), # Sheep ownership
         own_goats=ifelse(hh_goats>0,1,0), # Goat ownership
         hh_poultry=rowSums(across(hh_chickens:hh_othpoultry),na.rm=T), # Poultry ownership
         own_poultry=ifelse(hh_poultry>0,1,0), # Poultry ownership
         own_dogs = ifelse(hh_dogs>0,1,0), # Dog ownership
         own_cats=ifelse(hh_cats>0,1,0), # Cat ownership
         enter_home=ifelse(enter_home==1,1,0), # animals enter the home
         members_travel=ifelse(members_travel==1,1,0), # Household members travel with animals
         food_insects=ifelse(food_insects==1,0,1), # Insects observed in food preparation area
         latrine_flies=ifelse(latrine_flies==1,1,0), # Flies observed in latrine area
         animal_feces_dispose=ifelse(animal_feces>1,1,0)) # Animal feces disposal in courtyard

HH_survey <- HH_survey %>%
  select(hhnum,date,village,county,child_age,PLprob,BPL,BFPL,breastmilk,food,crawls,walks,resp_primary,
         child_gender,total_members,compound_num,pregnancy_term,pregnancy_ab,delivery_method,
         walk_time,latrine,animal_feces,hh_own_livestock,total_animals_own,own_cattle,
         own_goats,own_sheep,own_poultry,own_dogs,own_cats,hh_cattle,hh_goats,hh_sheep,
         hh_poultry,hh_dogs,hh_cats,enter_home,members_travel,water_age,food_insects,
         latrine_flies,animal_feces_dispose) # selecting only variables needed for analysis



##### Merge cstool & HH survey data ######
cstool_risk$hhcode <- as.numeric(substr(cstool_risk$Sample,11,13))

cstool_HH <- merge(cstool_risk,HH_survey,by.x="hhcode",by.y="hhnum")
cstool_HH$agegroup <- factor(ifelse(cstool_HH$child_age<3,'0-2 months',
                                    ifelse(cstool_HH$child_age>=3&cstool_HH$child_age<6,'3-5 months',
                                           ifelse(cstool_HH$child_age>=6&cstool_HH$child_age<12,'6-11 months',
                                                  ifelse(cstool_HH$child_age>=12,'12-23 months',NA)))),
                             levels=c("0-2 months","3-5 months","6-11 months","12-23 months"))

outcomes <- names(cstool_HH)[3:19]
outcomes_bin <- outcomes[c(2,4,6,8,10:17)]
indep <- names(cstool_HH)[24:60]
indep_noadj <- indep

model_pairs <- expand.grid(outcomes_bin,indep_noadj)
model_pairs$Var1 <- as.character(model_pairs$Var1)
model_pairs$Var2 <- as.character(model_pairs$Var2)

model_results <- data.frame(matrix(ncol=8,nrow=0))

for (i in 1:nrow(model_pairs)){
  modeldf <- cstool_HH[,c(model_pairs[i,1],model_pairs[i,2],"hhcode")]
  names(modeldf)[1:2] <- c('var1','var2')
  test5 <- glm(var1 ~ var2, family=binomial,data = modeldf)
  test5_cl <- coeftest(test5, vcov = vcovCL, cluster = ~hhcode)
  
  modeldf2 <- cstool_HH[,c(model_pairs[i,1],model_pairs[i,2],"agegroup","hhcode")]
  names(modeldf2)[1:2] <- c('var1','var2')
  test6 <- glm(var1 ~ var2+agegroup, family = binomial, data = modeldf2)
  test6_cl <- coeftest(test6, vcov = vcovCL, cluster = ~hhcode)
  
  model_out <- c(model_pairs[i,1],model_pairs[i,2],
                 test5$coefficients[2],test5_cl[2,4],test5_cl[2,2]*1.96,
                 test6$coefficients[2],test6_cl[2,4],test6_cl[2,2]*1.96)
  model_results <- rbind(model_results,model_out)
}

names(model_results) <- c("var1","var2",
                          "estimate_unadj","p_unadj","ci_unadj",
                          "estimate_adj","p_adj","ci_adj")

# convert to odds ratios
model_results <- model_results %>%
  mutate(across(c(estimate_unadj,ci_unadj,estimate_adj,ci_adj),\(x) (as.numeric(x))))
model_results$ci_unadj[model_results$estimate_unadj<0] <- -1*model_results$ci_unadj[model_results$estimate_unadj<0]
model_results$ci_adj[model_results$estimate_adj<0] <- -1*model_results$ci_adj[model_results$estimate_adj<0]
model_results$lower_unadj <- model_results$estimate_unadj-model_results$ci_unadj
model_results$upper_unadj <- model_results$estimate_unadj+model_results$ci_unadj
model_results$lower_adj <- model_results$estimate_adj-model_results$ci_adj
model_results$upper_adj <- model_results$estimate_adj+model_results$ci_adj


model_results[,c(3,5,6,8:12)] <- exp(model_results[,c(3,5,6,8:12)])

modelsig <- model_results %>%
  filter(p_unadj<0.05) %>%
  arrange(var1,p_unadj)

# Make table to output results in for now

keyvars <- c("child_gender","breastmilk","PLprob","food","crawls","walks","pregnancy_ab",
             "resp_primary","latrine","food_insects","hh_own_livestock","own_cattle",
             "own_poultry","own_dogs","animal_feces_dispose","walk_time")
model_plot2 <- filter(model_plot,var2%in%keyvars)


#### Age group OR #####

model_results <- data.frame(matrix(ncol=5,nrow=0))

for (i in 1:length(outcomes_bin)){
  modeldf <- cstool_HH[,c(outcomes_bin[i],"agegroup","hhcode")]
  names(modeldf)[1:2] <- c('var1','var2')
  test5 <- glm(var1 ~ var2, family=binomial,data = modeldf)
  test5_cl <- coeftest(test5, vcov = vcovCL, cluster = ~hhcode)
 
  model_out <- c(outcomes_bin[i],"agegroup",
                 test5$coefficients[2],test5_cl[2,4],test5_cl[2,2]*1.96)
  model_out2 <- c(outcomes_bin[i],"agegroup",
                 test5$coefficients[3],test5_cl[3,4],test5_cl[3,2]*1.96)
  model_out3 <- c(outcomes_bin[i],"agegroup",
                 test5$coefficients[4],test5_cl[4,4],test5_cl[4,2]*1.96)
  model_results <- rbind(model_results,model_out,model_out2,model_out3)
}

names(model_results) <- c("var1","var2",
                          "estimate_unadj","p_unadj","ci_unadj")

# convert to odds ratios
model_results <- model_results %>%
  mutate(across(c(estimate_unadj,ci_unadj),\(x) (as.numeric(x))))
model_results$ci_unadj[model_results$estimate_unadj<0] <- -1*model_results$ci_unadj[model_results$estimate_unadj<0]
model_results$lower_unadj <- model_results$estimate_unadj-model_results$ci_unadj
model_results$upper_unadj <- model_results$estimate_unadj+model_results$ci_unadj

model_results[,c(3,5:7)] <- exp(model_results[,c(3,5:7)])

modelsig <- model_results %>%
  filter(p_unadj<0.05) %>%
  arrange(var1,p_unadj)

model_plot3 <- model_results %>%
  mutate(est = round(estimate_unadj,2),
         lower=round(lower_unadj,2),
         upper=round(upper_unadj,2),
         disp_text=paste0(est, " (",lower,", ",upper,")")) %>%
  select(var1,var2,disp_text)
ages <- c("3-5 months","6-11 months","12-23 months")
model_plot3$agegroup <- rep(ages,12)
model_plot4 <- model_plot3 %>%
  pivot_wider(names_from="var1",values_from="disp_text")
