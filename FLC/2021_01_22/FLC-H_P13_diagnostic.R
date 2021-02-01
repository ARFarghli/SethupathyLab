library(ggplot2)
library(tidyr)
library(dplyr)
library(ggpubr)
library(viridis)

install.packages("ggbeeswarm")
install.packages("viridis")

setwd("~/Documents/Sethupathylab/FLC/2021_01_22")

data_summary <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  data_sum <- rename(data_sum, c("mean" = varname))
  return(data_sum)
}

qPCR.df <- read.csv("FLC_H_qPCR_diagnostic_messy.csv", row.names = 1)

qPCR.tidy <- data.frame(qPCR.df) %>% 
  mutate(genotype = rownames(qPCR.df))

qPCR.tidy <- qPCR.tidy %>% 
  gather(key = sample, value = measurement, -genotype) %>% 
  separate(col = sample, into = c("experiment", "replicate"), sep = "_")

qPCR.summary <- data_summary(qPCR.tidy, varname = "measurement", 
                             groupnames = "genotype")

qPCR.tidy$genotype <- factor(qPCR.tidy$genotype, levels = c("Rps9", "DP_Fusion", "SLC16A", "CA12"))
qPCR.summary$genotype <- factor(qPCR.summary$genotype, levels = c("Rps9", "DP_Fusion", "SLC16A", "CA12"))

p <- ggplot(qPCR.summary, aes(x = genotype, y = measurement, fill = genotype)) + 
  geom_bar(stat="identity",
           color="black", 
           position=position_dodge(),
           width = .5) +
  scale_fill_brewer() +
  geom_jitter(data = qPCR.tidy, aes(x = genotype, y = measurement), width = 0.2) +
  geom_errorbar(aes(ymin = measurement - sd, ymax = measurement + sd),
                width = .4,
                position=position_dodge(.9)) +
  theme_classic() +
  theme(legend.position = "none") +
  labs(title = "FLC Marker Gene Expression - Passage 13", 
       x = "Target", 
       y = "Average Cq")
print(p)

png(units = 'in', width = 5, height = 5, res = 250)
p
dev.off()
