# Script to make figure 3
library(tidyverse)
library(gridExtra)
library(cowplot)

TAC <- read.csv("pathogen_quantities_samples.csv")

TAC$Quant[TAC$Quant<=0] <- NA
TAC$log10_quant <- log10(TAC$Quant)

# Plot of quantities by pathogen
ggplot(TAC,aes(x=pathogens,y=log10_quant)) + geom_jitter(width=0.3,height=0,color="#F8766D",alpha=0.9) +
  geom_boxplot(fill=NA,outlier.shape=NA,width=0.75,linewidth=.75) + 
  geom_hline(yintercept=1) + theme_bw() +
  theme(legend.position = "none") + ylim(0,10)

# Plotting fraction of samples positive
positivity <- TAC %>% group_by(pathogens) %>%
  summarize(npos=sum(!is.na(log10_quant)),n=n(),
            frac_pos=npos/n*100)

ggplot(positivity,aes(x=pathogens,y=frac_pos)) + geom_col(color='black',fill='blue') +
  theme_bw() +ylim(0,100) +
  theme(legend.position = "none")



