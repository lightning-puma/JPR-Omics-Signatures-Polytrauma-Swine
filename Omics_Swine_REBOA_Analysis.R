# Created By: Ian S. LaCroix


#### Libraries ####
library(tidyverse)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggsci)
library(plotrix)
library(ggpubr)
library(rstatix)
library(ggrepel)
library(stringr)
library(ggh4x)


#### Read in Data ####

dat.meta <- read.csv("20220324_DoD_Pig_Meta_Data.csv")
dat.clin <- read.csv("Swine Database_5.25.2022.csv")
dat.prot <- read.csv("20220720_DoD_Trauma_Pig_Plasma_ILC03_Protein_Report.csv")
dat.metab <- read.csv("20220922_DoD_Trauma_Pig_Metabolomics_Final_Report_RelativeQuant.csv")


# save metabolite names with correct symbol format

correct_metab_symbols <- dat.metab %>% 
  select(Col.1) %>% 
  slice(-1)


# fix time notation
Date <- dat.meta$Date
dat.meta$Date <- as.Date(dat.meta$Date, origin = "1899-12-30")


#### Proteomics Data Coercion ####
dat.prot <- dat.prot  %>% distinct(PG.Genes, .keep_all = TRUE)
# remove columns other than gene name and transpose
dat.prot.t <- t(dat.prot[,-c(1:3,5:10)])
# make column names gene names
colnames(dat.prot.t) <- dat.prot.t[1,]
# remove the gene name row
dat.prot.t <- dat.prot.t[-1,]
# create new data frame for tech mixes
dat.prot.tm <- dat.prot.t[c(334:343),]
# remove tech mix rows from prot data
dat.prot.t <- dat.prot.t[-c(334:343),]
# store the row names as separate data frame
temp <- data.frame("label" = rownames(dat.prot.t))
# separate label into columns
temp <- tidyr::separate(data = temp ,col = label, sep="_", 
                        into = c("col.1", "col.2", "col.3", "col.4","col.5"))
# keep the location and id columns
# location column stores control data labels, patient location: field, ED, OR
dat.prot.t<- data.frame("Vial" = temp[,1], dat.prot.t)
rownames(dat.prot.t) <- 1:nrow(dat.prot.t)
# remove rows that are not in the union of dat.prot and dat.meta
dat.prot.t <- dat.prot.t %>% subset(Vial %in% dat.meta$Vial)
dat.meta <- dat.meta %>% subset(Vial %in% dat.prot.t$Vial)
# merge meta data with proteomics data
dat <- cbind(dat.meta, dat.prot.t)
# check that merged data IDs are identical
# run identical function on the two columns labeled Vial
# output needs to be TRUE
identical(dat[,1], dat[,14])
# remove any columns that are not required for analysis
dat <- dat[,-c(10:14)]

swine_prot_no_imputation <- dat
write.csv(swine_prot_no_imputation,
          row.names= FALSE, 
          "Swine_Proteomics_No_Imputation.csv")


#### Imputation ####
# store meta data
ids <- dat[,c(1:9)]
# remove meta dat.prot from proteomics
dat <- dat[,-c(1:9)]
# replace missing values with 0
dat[is.na(dat)] <- 0
# convert dat.prot data frame to numeric
dat <-  as.data.frame(lapply(dat,as.numeric))

dat.no.imputation <- cbind(ids, dat)

for (i in 1:ncol(dat)) {
  dat[,i][dat[,i] == 0] <- min(dat[,i][which(dat[,i] > 0)])*0.20
}
# merge ids with data frame
dat_prot <- cbind(ids, dat)


#### Metabolomics Data Coercion ####

dat.metab.t <- t(dat.metab)
colnames(dat.metab.t) <- dat.metab.t[1, ]
dat.metab.t <- dat.metab.t[-1, ]
rownames(dat.metab.t) <- 1:nrow(dat.metab.t)
dat.metab.t <- as.data.frame(dat.metab.t)

swine_metab_no_imputation <- dat.metab.t
write.csv(swine_metab_no_imputation, 
          row.names = FALSE,
          "Swine_Metabolomics_No_Imputation.csv")

#### Imputation ####
# store meta data
ids <- dat.metab.t[,1]
# remove meta dat.prot from proteomics
dat.metab.t <- dat.metab.t[,-1]
# replace missing values with 0
dat.metab.t[is.na(dat.metab.t)] <- 0
# convert dat.prot data frame to numeric
dat.metab.t <-  as.data.frame(lapply(dat.metab.t,as.numeric))

dat.metab.t.no.imputation <- cbind(ids, dat.metab.t)

for (i in 1:ncol(dat.metab.t)) {
  dat.metab.t[,i][dat.metab.t[,i] == 0] <- min(dat.metab.t[,i][which(dat.metab.t[,i] > 0)])*0.20
}
# merge ids with data frame
dat_metab <- cbind(ids, dat.metab.t) %>% 
  dplyr::rename(Vial = ids)


#### Merge -Omics ####

dat_merged <- dat_prot %>% 
  inner_join(dat_metab, by = "Vial")

out <- dat_merged %>% 
  filter(Model == "REBOA Zone I + DCBI" | Model == "REBOA Zone III + DCBI"| Model == "DCBI")

write.csv(out, row.names = FALSE, "Swine_REBOA_Merged_Imputed.csv")


#### Figure 2 Experiment Design and Vital Measurements ####
##### Summary Bar Plot #####
dat_meta_sub <- dat.meta %>% 
  filter(Model == "REBOA Zone I + DCBI" | Model == "REBOA Zone III + DCBI"| Model == "DCBI") %>% 
  select(Pig.id, Model, Time.from.EOS.min)


dat_meta_sub$Model <- factor(dat_meta_sub$Model, levels = c("DCBI", "REBOA Zone I + DCBI", "REBOA Zone III + DCBI"))

bar_plot <- ggplot(dat_meta_sub, aes(x = Time.from.EOS.min, fill = Model))+
  geom_bar()+
  scale_x_continuous(limits = c(-140, 260), breaks = c(-120, 0, 30, 60, 120, 180, 240)) +
  scale_y_continuous(expand = c(0,0), limits = c(0,24), breaks = seq(0,24, by=4))+
  
  
  scale_fill_manual(values = c(
                                "#3B4992FF", #DCBI
                                "#BB0021FF", #REBOA Zone I + DCBI
                                "#008B45FF") #REBOA Zone III + DCBI
                    )+
  
  ggtitle("Sample Size")+
  
  xlab("Time from EOS (min)")+
  
  ylab("Count")+
  
  theme(legend.position = "none", 
        plot.margin = margin(t = 15, r = 10, b = 8, l = 8, unit = "pt"),
        plot.title = element_text(hjust = 0.5, size = 8), 
        axis.title = element_text(size = 6, color = 'black'), 
        axis.text.x = element_text(size = 6, vjust = 0.6, angle = 90, color = 'black'),
        axis.text.y = element_text(size = 6, color = 'black'),
        axis.ticks = element_line(size = 0.25),
        plot.background = element_rect(fill = "transparent",colour = NA),
        panel.background = element_rect(fill = NA),
        panel.border = element_rect(colour = "black", fill=NA, size = 0.25),
        #axis.line = element_line(color = 'black', size = 0.25, linetype = 'solid'),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
  

plot(bar_plot)

ggsave(bar_plot, 
       file = paste(Sys.Date(),"Pig_Count_Per_Model_Over_Time.svg", sep = "_"),
       width = 2, height = 2, units = "in", dpi = 300, bg = "transparent")


##### Sample Size Counts #####

nrow(filter(dat.meta, Model == "SHAM"))

dat.meta %>% 
  filter(Model == "SHAM") %>% 
  group_by(Time.from.EOS.min) %>% 
  dplyr::summarise(n = n())

dat.meta %>% 
  filter(Model == "DCBI") %>% 
  group_by(Time.from.EOS.min) %>% 
  dplyr::summarise(n = n())

dat.meta %>% 
  filter(Model == "REBOA Zone I + DCBI") %>% 
  group_by(Time.from.EOS.min) %>% 
  dplyr::summarise(n = n())

dat.meta %>% 
  filter(Model == "REBOA Zone III + DCBI") %>% 
  group_by(Time.from.EOS.min) %>% 
  dplyr::summarise(n = n())


##### Clin. Var. Line Plots #####
# Function for data summary table used for all longitudinal line plots
data_summary <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(median = median(x[[col]], na.rm=TRUE),
      sd = std.error(x[[col]], na.rm=TRUE))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  data_sum <- rename(data_sum, c("median" = varname))
  return(data_sum)
}


dat.line.plot <- dat.clin %>%
  filter(Model == "DCBI"| Model == "REBOA Zone I + DCBI" | Model == "REBOA Zone III + DCBI") %>% 
  dplyr::rename(Time.from.EOS.min = Time..Min..from.End.of.Shock, Pig.id = Pig.Number)


dat.line.plot$Model <- factor(dat.line.plot$Model, 
                              levels = c("DCBI", "REBOA Zone I + DCBI", "REBOA Zone III + DCBI"))


clin_vars <- c("SBP", "DBP", "MAP", "Pulse", "RR")

breaks = list(seq(0,150, by = 25), seq(0, 100, by = 20), seq(0, 100, by = 20), seq(40, 180, by = 20), seq(5, 25, by = 5))

limits = list(c(0, 150), c(0,100), c(0, 100), c(40, 180), c(10, 25))

# loop that exports temporal trend plot per clinical variable of interest
for(i in 1:length(clin_vars)){  
  
  # a pipe that does everything you need
  AOI <- dat.line.plot %>%
    select(c(clin_vars[i], 'Time.from.EOS.min', 'Model', 'Pig.id')) %>% 
    rename_(QUANT = names(.)[1]) %>% 
    group_by(Pig.id) %>% 
    dplyr::mutate(normQUANT = QUANT/QUANT[row_number()==1]) %>%
    dplyr::mutate(Pct.change = (normQUANT/normQUANT[row_number()==1]-1)*100) %>% 
    data_summary(varname = "QUANT",
                 groupnames = c("Model", "Time.from.EOS.min")) %>%  
    filter(!is.na(Model)== TRUE)
  
  p_1 <- ggplot(AOI, aes(x = Time.from.EOS.min, y = QUANT, group = Model, color = Model)) +
    #geom_errorbar(aes(ymin = QUANT - sd, ymax = QUANT + sd), width = 0.1) +
    
    geom_line(size = 0.5) + geom_point(size = 1) +
  
    scale_color_manual(values = c(
                                  "#3B4992FF", #DCBI
                                  "#BB0021FF", #REBOA Zone I + DCBI
                                  "#008B45FF") #REBOA Zone III + DCBI
                       )+
    
    labs(fill = "Method") + 
    
    scale_y_continuous(expand = c(0,0), breaks = breaks[[i]], limits = limits[[i]])+
    
    scale_x_continuous(breaks = c(-120, 0, 30, 60, 120, 180, 240), 
                       labels = c("-120", "0", "30","60","120","180","240")) +
    xlab("Time frrom EOS (min)") + ylab("") + ggtitle(clin_vars[i]) + #"abs(Median PG.Quantity)" 
    
    theme(legend.position = "none", 
          plot.margin = margin(t = 15, r = 10, b = 8, l = 8, unit = "pt"),
          plot.title = element_text(hjust = 0.5, size = 8), 
          axis.title = element_text(size = 6, color = 'black'), 
          axis.text.x = element_text(size = 6, vjust = 0.6, angle = 90, color = 'black'),
          axis.text.y = element_text(size = 6, color = 'black'),
          axis.ticks = element_line(size = 0.25),
          plot.background = element_rect(fill = "transparent",colour = NA),
          panel.background = element_rect(fill = NA),
          panel.border = element_rect(colour = "black", fill=NA, size = 0.25),
          #axis.line = element_line(color = 'black', size = 0.25, linetype = 'solid'),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
  
  
  plot(p_1)
  
  # Save the plot to working directory
  ggsave(p_1,
         file = paste(Sys.Date(),"DoD_Trauma_Pig_Line_Plot_without_std_error", clin_vars[i],".svg", sep = "_"),
         width = 2, height = 2, units = "in", dpi = 300, bg = "transparent")
  
}  


#### Figure 3 Global Heat Maps ####

out_1 <- dat_prot %>% 
  filter(Model == "DCBI"| Model == "REBOA Zone I + DCBI" | Model == "REBOA Zone III + DCBI") %>%
  select(-c(2,4:9))

write.csv(out_1,
          row.names = FALSE,
          paste(Sys.Date(), "MetaboAnalyst_Ready_Global_HeatMap_Proteomics.csv"))

temp <- out_1 %>% 
  select(Vial, Model)

out_2 <- dat_metab %>%
  right_join(temp, by = "Vial") %>% 
  relocate(Model, .after = Vial)

write.csv(out_2,
          row.names = FALSE,
          paste(Sys.Date(), "MetaboAnalyst_Ready_Global_HeatMap_Metabolomics.csv"))



#### Supp Fig 1-5 ####

# "#cc0000", TI + HS
# "#5a00b3", HS
# "#0066cd", TI
# "#00994d"  SHAM

##### MetaboAnalyst Export #####
# proteomics preparation for metabo analyst
out_1 <- dat_prot %>% 
  filter(Model == "DCBI"| Model == "REBOA Zone I + DCBI" | Model == "REBOA Zone III + DCBI") %>% 
  filter(!Time.from.EOS.min == -120) %>% 
  select(-c(2, 4:9))

# out_1 <- dat_prot %>% 
#   filter(Model == "DCBI" | Model == "REBOA Zone I + DCBI"| Model == "REBOA Zone III + DCBI") %>% 
#   filter(!Time.from.EOS.min == -120) %>% 
#   select(-c(2, 4:9))

# autoscale normalize this output file and use for PLSDA
write.csv(out_1, 
          row.names = FALSE,
          paste(Sys.Date(), "Proteomics_Comparison_of_Models_MetaboAnalyst_Ready_char.csv", sep = "_"))


out_2 <- dat_prot %>% 
  filter(Model == "DCBI"| Model == "REBOA Zone I + DCBI" | Model == "REBOA Zone III + DCBI") %>% 
  filter(Time.from.EOS.min == 240) %>% 
  select(-c(2, 4:9)) %>% 
  mutate(Model = case_when(
                           Model == "DCBI" ~ 0,
                           Model == "REBOA Zone I + DCBI" ~ 1,
                           Model == "REBOA Zone III + DCBI" ~ 2))


# autoscale normalize this output file and use for Heat Map
write.csv(out_2, 
          row.names = FALSE,
          paste(Sys.Date(), "Proteomics_Comparison_of_Models_MetaboAnalyst_Ready_EOS_min_240.csv", sep = "_"))

# metabolomics preparation for metabo analyst

meta <- select(dat.meta, Vial, Model, Time.from.EOS.min)

out_3 <- dat_metab %>% 
  inner_join(meta, by = "Vial") %>% 
  relocate(Model, Time.from.EOS.min, .after = Vial) %>% 
  filter(Model == "DCBI"| Model == "REBOA Zone I + DCBI" | Model == "REBOA Zone III + DCBI") %>% 
  filter(!Time.from.EOS.min == -120) %>% 
  select(-Time.from.EOS.min)

# autoscale normalize this output file and use for PLSDA
write.csv(out_3, 
          row.names = FALSE,
          paste(Sys.Date(), "Metabolomics_Comparison_of_Models_MetaboAnalyst_Ready_char.csv", sep = "_"))


out_4 <- dat_metab %>% 
  inner_join(meta, by = "Vial") %>% 
  relocate(Model, Time.from.EOS.min, .after = Vial) %>% 
  filter(Model == "DCBI"| Model == "REBOA Zone I + DCBI" | Model == "REBOA Zone III + DCBI") %>% 
  filter(Time.from.EOS.min == 240) %>% 
  select(-Time.from.EOS.min) %>% 
  mutate(Model = case_when(
                           Model == "DCBI" ~ 0,
                           Model == "REBOA Zone I + DCBI" ~ 1,
                           Model == "REBOA Zone III + DCBI" ~ 2))

# autoscale normalize this output file and use for Heat Map
write.csv(out_4, 
          row.names = FALSE,
          paste(Sys.Date(), "Metabolomics_Comparison_of_Models_MetaboAnalyst_Ready_EOS_min_240.csv", sep = "_"))






#### Figure 2 ####
# proteomics temporal trend line plot

# Function for data summary table used for all longitudinal line plots
data_summary <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(median = median(x[[col]], na.rm=TRUE),
      sd = std.error(x[[col]], na.rm=TRUE))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  data_sum <- rename(data_sum, c("median" = varname))
  return(data_sum)
}


dat.line.plot<- dat_prot %>% 
  filter(Model == "REBOA Zone I + DCBI" | Model == "REBOA Zone III + DCBI"| Model == "DCBI") 
  #filter(!Time.from.EOS.min == -120) 


dat.line.plot$Model <- factor(dat.line.plot$Model, 
                              levels = c("DCBI", "REBOA Zone I + DCBI", "REBOA Zone III + DCBI"))


analytes <- c("C4A",  "PRDX6", "CRYZ", "FAH", "GSTO1", "HPD", "ARG1", 
              "AKR1C3", "ENO3", "IDH1", "COL18A1", "COMT", "CFI", "FGG","FGA", "FGB", "CRP",
              "F2", "F5", "F10", "HP", "C4BPA", "C5", "C6", "C7", "C8A", "C8B", "C9", "PKM","ALDOB","ALDOC", "LDHA", "LDHB",
              "VWF", "SERPINA1", "SERPINE1", "H4", "H1.3")

analytes <- c("PKM")

# loop that exports temporal trend plot per analyte of interest
for(i in 1:length(analytes)){  
  
  # a pipe that does everything you need
  AOI <- dat.line.plot %>%
    select(c(analytes[i], 'Time.from.EOS.min', 'Model', 'Pig.id')) %>% 
    rename_(QUANT = names(.)[1]) %>% 
    group_by(Pig.id) %>% 
    dplyr::mutate(normQUANT = QUANT/QUANT[row_number()==1]) %>%
    dplyr::mutate(Pct.change = (normQUANT/normQUANT[row_number()==1]-1)*100) %>% 
    data_summary(varname = "Pct.change",
                 groupnames = c("Model", "Time.from.EOS.min")) %>%  
    filter(!is.na(Model)== TRUE)
  
  p_1 <- ggplot(AOI, aes(x = Time.from.EOS.min, y = Pct.change, group = Model, color = Model)) +
    
    geom_errorbar(aes(ymin = Pct.change - sd, ymax = Pct.change + sd), width = 0.1) +
    
    geom_line(size = 1) + geom_point(size = 1.5) +
    
    scale_color_manual(values = c(
      
                                  "#3B4992FF", #DCBI
                                  "#BB0021FF", #REBOA Zone I + DCBI
                                  "#008B45FF" #REBOA Zone III + DCBI
                                  
                                  ))+
  
  
    labs(fill = "Method") + 
    
    scale_x_continuous(breaks = c(-120, 0, 30, 60, 120, 180, 240), 
                       labels = c("-120", "0", "30","60","120","180","240")) +
  
    xlab("") + ylab("") + ggtitle(analytes[i]) + #"abs(Median PG.Quantity)" 

    theme(plot.title = element_text(hjust = 0.5, size = 8),
          axis.title.x = element_text(size = 8),
          axis.text.x = element_text(size = 6, angle = 90, hjust = 0.95, vjust = 0.5, color = "black"),
          axis.text.y = element_text(size = 6, color = "black"),
          axis.title.y = element_text(size = 8),
          panel.border = element_rect(colour = "black", fill=NA, size = 0.25),
          panel.background = element_rect(fill = "transparent",colour = NA),
          plot.background = element_rect(fill = "transparent",colour = NA),
          legend.position = "none"
    )
  
  plot(p_1)
  
  # Save the plot to working directory
  ggsave(p_1,
         file = paste(Sys.Date(),"Fig_2_Proteomics_Temp_Trends_REBOA_Models", analytes[i],".tiff", sep = "_"),
         width = 1.5, height = 1.5, units = "in", dpi = 300, bg = "transparent")
  
}  




#### Figure 3 ####
# metabolomics temporal trend line plot

# Function for data summary table used for all longitudinal line plots
data_summary <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(median = median(x[[col]], na.rm=TRUE),
      sd = std.error(x[[col]], na.rm=TRUE))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  data_sum <- rename(data_sum, c("median" = varname))
  return(data_sum)
}


meta <- select(dat.meta, Vial, Pig.id, Model, Time.from.EOS.min)

dat.line.plot <- dat_metab %>% 
  inner_join(meta, by = "Vial") %>% 
  relocate(Model, Time.from.EOS.min, .after = Vial) %>% 
  filter(Model == "REBOA Zone I + DCBI" | Model == "REBOA Zone III + DCBI"| Model == "DCBI") 
  #filter(!Time.from.EOS.min == -120) 


dat.line.plot$Model <- factor(dat.line.plot$Model, 
                              levels = c("DCBI", "REBOA Zone I + DCBI", "REBOA Zone III + DCBI"))


analytes <- c("Fumarate", "Malate", "Lactate", "Succinate", "L.lysine", "L.histidine", "L.proline", "L.alanine", "Bilirubin",
              "Biliverdin","Inosine", "kynurenine", "X5.Hydroxykynurenine", "X5.hydroxykynurenamine","S.Adenosyl.L.methionine",
              "Citrate", "Oxaloacetate", "Itaconate", "L.glutamine", "X5.hydroxykynurenamine")



# carnitines
analytes <- colnames(dat.line.plot[, c(94:113)])


analytes <- colnames(dat.line.plot[, c(4:113)])


# loop that exports temporal trend plot per analyte of interest
for(i in 1:length(analytes)){  
  
  # a pipe that does everything you need
  AOI <- dat.line.plot %>%
    select(c(analytes[i], 'Time.from.EOS.min', 'Model', 'Pig.id')) %>% 
    rename_(QUANT = names(.)[1]) %>% 
    group_by(Pig.id) %>% 
    dplyr::mutate(normQUANT = QUANT/QUANT[row_number()==1]) %>%
    dplyr::mutate(Pct.change = (normQUANT/normQUANT[row_number()==1]-1)*100) %>% 
    data_summary(varname = "Pct.change",
                 groupnames = c("Model", "Time.from.EOS.min")) %>%  
    filter(!is.na(Model)== TRUE)
  
  p_1 <- ggplot(AOI, aes(x = Time.from.EOS.min, y = Pct.change, group = Model, color = Model)) +
    
    geom_errorbar(aes(ymin = Pct.change - sd, ymax = Pct.change + sd), width = 0.1) +
    
    geom_line(size = 1) + geom_point(size = 1.5) +
    
    scale_color_manual(values = c(
      
                                    "#3B4992FF", #DCBI
                                    "#BB0021FF", #REBOA Zone I + DCBI
                                    "#008B45FF" #REBOA Zone III + DCBI
                                    
      ))+
    
    labs(fill = "Method") + 
    
    scale_x_continuous(breaks = c(-120, 0, 30, 60, 120, 180, 240), 
                       labels = c("-120", "0", "30","60","120","180","240")) +
    
    xlab("") + ylab("") + ggtitle(correct_metab_symbols[i,]) + #"abs(Median PG.Quantity)" 
    
    theme(plot.title = element_text(hjust = 0.5, size = 8),
          axis.title.x = element_text(size = 8),
          axis.text.x = element_text(size = 6, angle = 90, hjust = 0.95, vjust = 0.5, color = "black"),
          axis.text.y = element_text(size = 6, color = "black"),
          axis.title.y = element_text(size = 8),
          panel.border = element_rect(colour = "black", fill=NA, size = 0.25),
          panel.background = element_rect(fill = "transparent",colour = NA),
          plot.background = element_rect(fill = "transparent",colour = NA),
          legend.position = "none"
    )
  
  
  # Save the plot to working directory
  ggsave(p_1,
         file = paste(Sys.Date(),"Fig_3_Metabolomics_Temp_Trends_REBOA_Models", analytes[i],".tiff", sep = "_"),
         width = 1.5, height = 1.5, units = "in", dpi = 300, bg = "transparent")
  
}  





 #### Figure 4 ####

# Absolute Quant of Succinate and Lactate

# Function for data summary table used for all longitudinal line plots
data_summary <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(median = median(x[[col]], na.rm=TRUE),
      sd = std.error(x[[col]], na.rm=TRUE))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  data_sum <- rename(data_sum, c("median" = varname))
  return(data_sum)
}

# visualize temporal trends for absolute quant data

dat_abs_quant <- read.csv("20220922_DoD_Trauma_Pig_Metabolomics_Final_Report_AbsoluteQuant.csv") %>% 
  t() %>% 
  as.data.frame() %>% 
  tibble::rownames_to_column("Vial") %>% 
  janitor::row_to_names(row_number = 1) %>% 
  left_join(dat.meta, by = "Vial") %>% 
  relocate(Pig.id, Model, Time.from.EOS.min, .after = Vial) %>% 
  select(-c(7:15)) %>% 
  filter(Model == "REBOA Zone I + DCBI" | Model == "REBOA Zone III + DCBI"| Model == "DCBI") 


dat_abs_quant$Model <- factor(dat_abs_quant$Model, 
                              levels = c("DCBI", "REBOA Zone I + DCBI", "REBOA Zone III + DCBI"))


dat_abs_quant$Succinate <- as.numeric(dat_abs_quant$Succinate)
dat_abs_quant$Lactate <- as.numeric(dat_abs_quant$Lactate)

# a pipe that does everything you need
AOI <- dat_abs_quant %>%
  select(c(Succinate, 'Time.from.EOS.min', 'Model', 'Pig.id')) %>% 
  rename_(QUANT = names(.)[1]) %>% 
  group_by(Pig.id) %>% 
  dplyr::mutate(normQUANT = QUANT/QUANT[row_number()==1]) %>%
  dplyr::mutate(Pct.change = (normQUANT/normQUANT[row_number()==1]-1)*100) %>% 
  data_summary(varname = "QUANT",
               groupnames = c("Model", "Time.from.EOS.min")) %>%  
  filter(!is.na(Model)== TRUE)


# create the succinate plot
p_1 <- ggplot(AOI, aes(x = Time.from.EOS.min, y = QUANT, group = Model, color = Model)) +
  
  geom_errorbar(aes(ymin = QUANT - sd, ymax = QUANT + sd), width = 0.1) +
  
  geom_line(size = 1) + geom_point(size = 1.5) +
  
  scale_color_manual(values = c(
                                
                                "#3B4992FF", #DCBI
                                "#BB0021FF", #REBOA Zone I + DCBI
                                "#008B45FF" #REBOA Zone III + DCBI
                                
                                ))+
  
  labs(fill = "Method") + 
  
  scale_x_continuous(breaks = c(-120, 0, 30, 60, 120, 180, 240), 
                     labels = c("-120", "0", "30","60","120","180","240")) +
  
  xlab("") + ylab("") + ggtitle("Succinate") + 
  
  theme(plot.title = element_text(hjust = 0.5, size = 8),
        axis.title.x = element_text(size = 8),
        axis.text.x = element_text(size = 6, angle = 90, hjust = 0.95, vjust = 0.5, color = "black"),
        axis.text.y = element_text(size = 6, color = "black"),
        axis.title.y = element_text(size = 8),
        panel.border = element_rect(colour = "black", fill=NA, size = 0.25),
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA),
        legend.position = "none"
  )


plot(p_1)

# Save the plot to working directory
ggsave(p_1,
       file = paste(Sys.Date(),"Metabolomics_Temp_Trend_Absolute_Quant_Succinate.tiff", sep = "_"),
       width = 1.5, height = 1.5, units = "in", dpi = 300, bg = "transparent")




# a pipe that does everything you need
AOI <- dat_abs_quant %>%
  select(c(Lactate, 'Time.from.EOS.min', 'Model', 'Pig.id')) %>% 
  rename_(QUANT = names(.)[1]) %>% 
  group_by(Pig.id) %>% 
  data_summary(varname = "QUANT",
               groupnames = c("Model", "Time.from.EOS.min")) %>%  
  filter(!is.na(Model)== TRUE)

# create the succinate plot
p_2 <- ggplot(AOI, aes(x = Time.from.EOS.min, y = 10*QUANT, group = Model, color = Model)) +
  
  geom_errorbar(aes(ymin = 10*QUANT - 10*sd, ymax = 10*QUANT + 10*sd), width = 0.1) +
  
  geom_line(size = 1) + geom_point(size = 1.5) +
  
  scale_color_manual(values = c(
    
                                "#3B4992FF", #DCBI
                                "#BB0021FF", #REBOA Zone I + DCBI
                                "#008B45FF" #REBOA Zone III + DCBI
    
    ))+
  
  labs(fill = "Method") + 
  
  scale_x_continuous(breaks = c(-120, 0, 30, 60, 120, 180, 240), 
                     labels = c("-120", "0", "30","60","120","180","240")) +
  
  xlab("") + ylab("") + ggtitle("Lactate") + 
  
  theme(plot.title = element_text(hjust = 0.5, size = 8),
        axis.title.x = element_text(size = 8),
        axis.text.x = element_text(size = 6, angle = 90, hjust = 0.95, vjust = 0.5, color = "black"),
        axis.text.y = element_text(size = 6, color = "black"),
        axis.title.y = element_text(size = 8),
        panel.border = element_rect(colour = "black", fill=NA, size = 0.25),
        panel.background = element_rect(fill = "transparent",colour = NA),
        plot.background = element_rect(fill = "transparent",colour = NA),
        legend.position = "none"
  )


plot(p_2)

# Save the plot to working directory
ggsave(p_2,
       file = paste(Sys.Date(),"Metabolomics_Temp_Trend_Absolute_Quant_Lactate.tiff", sep = "_"),
       width = 1.5, height = 1.5, units = "in", dpi = 300, bg = "transparent")








#### Figure 5 ####

##### Proteomics #####


dat.box.plot <- dat_prot %>% 
  relocate(Model, Time.from.EOS.min, .after = Vial) %>% 
  filter(Model == "REBOA Zone I + DCBI" | Model == "REBOA Zone III + DCBI"| Model == "DCBI") %>% 
  filter(Time.from.EOS.min == 30) 

dat.box.plot$Model <- factor(dat.box.plot$Model, 
                             levels = c("DCBI", "REBOA Zone I + DCBI", "REBOA Zone III + DCBI"))


##### loop with significance filter #####

# select analytes
analytes <- colnames(dat.box.plot[, c(10:ncol(dat.box.plot))])


for (i in 1:length(analytes)) {
  
  a <- select(dat.box.plot, c("Model", analytes[i]))
  colnames(a) <- c("Model", "intensity")
  
  # if the variance is not zero, then continue with loop
  # zero variance means all rows have the same value 
  # indicates the protein values were imputed for this particular time point
  if(!var(a$intensity)==0){
  
    stat.test <- a %>%
      rstatix::dunn_test(intensity~Model, p.adjust.method = 'bonferroni') %>% 
      rstatix::add_significance(p.col = "p.adj") %>% 
      rstatix::add_xy_position(x = "Day", dodge = 0.8) %>%
      filter(p.adj < 0.05) 
    
    
    
    kt <- kruskal.test(intensity ~ Model, data = a)$p.value
    
    # if the kruskal test p-value is significant, then proceed with the loop
    if(kt < 0.05){
    
      print(analytes[i])  
      
      boxplot <- ggplot(a, aes(x = Model, y = intensity))+
        stat_boxplot(geom = "errorbar", width = 0.5) + 
        geom_boxplot(position="dodge", outlier.shape = 16, width = 0.5,aes(fill = Model) ) +
        
        scale_fill_manual(values = c(
                                     "#3B4992FF", #DCBI
                                     "#BB0021FF", #REBOA Zone I + DCBI
                                     "#008B45FF" #REBOA Zone III + DCBI
        )) +
        
        #scale_y_continuous(breaks = breaks[[i]], limits = limits[[i]])+
        
        ggtitle(analytes[i])+
        
        labs(x=" ", y = " ", family = "Arial")+
        
        ggprism::add_pvalue(data = stat.test, 
                            bracket.nudge.y = 1,
                            label.size = 8, 
                            tip.length = 0)+
        
        #stat_pvalue_manual(data = stat.test, size = 10, hide.ns = TRUE)+
        
        #stat_compare_means(data=a,label.x = 1.2, label.y = label.y[i], family = "Arial", size = 5)+
        
        theme(plot.title = element_text(hjust = 0.5, size = 8),
              axis.title.x = element_text(size = 8),
              axis.text.x = element_text(size = 6, angle = 90, hjust = 0.95, vjust = 0.5, color = "black"),
              axis.text.y = element_text(size = 6, color = "black"),
              axis.title.y = element_text(size = 8),
              panel.border = element_rect(colour = "black", fill=NA, size = 0.25),
              panel.background = element_rect(fill = "transparent",colour = NA),
              plot.background = element_rect(fill = "transparent",colour = NA),
              legend.position = "none"
        )
      
      # view the plot
      #plot(boxplot)
      
      # Save the plot as PNG to working directory
      ggsave(boxplot,
             file = paste(Sys.Date(),analytes[i],"Figure_5_proteomics.tiff", sep = "_"),
             width = 1.5, height = 1.5, units = "in", dpi = 300, bg = "transparent")
    }
  }
}




##### loop without significance filter #####

analytes <- c("H4", "VTN", "TPP1", "SDC1", "PRDX6", "PEBP1", "MYOC", "MMP9", "MDH1", "LTF", "GPX3", "FABP6", 
              "COL2A1", "CHGB", "FGA", "FGB", "FGG", "HP", "LDHA","GAPDH", "BLVRB", "HBE1")


for (i in 1:length(analytes)) {
  
  a <- select(dat.box.plot, c("Model", analytes[i]))
  colnames(a) <- c("Model", "intensity")
  
  # if the variance is not zero, then continue with loop
  # zero variance means all rows have the same value 
  # indicates the protein values were imputed for this particular time point
  if(!var(a$intensity)==0){
    
    stat.test <- a %>%
      rstatix::dunn_test(intensity~Model, p.adjust.method = "hochberg") %>% 
      rstatix::add_significance(p.col = "p.adj") %>% 
      rstatix::add_xy_position() %>%
      filter(p.adj < 0.05) 
    
    
    
    kt <- kruskal.test(intensity ~ Model, data = a)$p.value
    
    # if the kruskal test p-value is significant, then proceed with the loop
  
    
    print(analytes[i])  
    
    boxplot <- ggplot(a, aes(x = Model, y = intensity))+
      stat_boxplot(geom = "errorbar", width = 0.5, lwd = 0.25) + 
      geom_boxplot(position="dodge", outlier.shape = 16, outlier.size = 0.5, width = 0.5, lwd = 0.25, aes(fill = Model) ) +
      
      scale_fill_manual(values = c( 
                                    "#3B4992FF", #DCBI
                                    "#BB0021FF", #REBOA Zone I + DCBI
                                    "#008B45FF" #REBOA Zone III + DCBI
                                    
                                    )) +
      
      #scale_y_continuous(breaks = breaks[[i]], limits = limits[[i]])+
      
      ggtitle(analytes[i])+
      
      labs(x=" ", y = " ", family = "Arial")+
      
      ggprism::add_pvalue(data = stat.test, 
                          bracket.size = 0.25,
                          bracket.nudge.y = 0.5,
                          step.increase = 0.001,
                          label.size = 4, 
                          tip.length = 0)+
      
      
      theme(plot.title = element_text(hjust = 0.5, size = 8),
            axis.title.x = element_text(size = 8),
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            axis.text.y = element_text(size = 6, color = "black"),
            axis.title.y = element_text(size = 8),
            axis.ticks.y = element_line(size = 0.25),
            axis.line = element_line(color = "black", size = 0.25),
            panel.background = element_rect(fill = "transparent",colour = NA),
            plot.background = element_rect(fill = "transparent",colour = NA),
            #panel.border = element_rect(color = "black", fill=NA, size = 0.25),
            legend.position = "none"
      )
    
    # view the plot
    #plot(boxplot)
    
    # Save the plot as PNG to working directory
    ggsave(boxplot,
           file = paste(Sys.Date(),analytes[i],"Figure_5_proteomics.tiff", sep = "_"),
           width = 1.5, height = 1.9, units = "in", dpi = 300, bg = "transparent")
  
  }
}





##### Metabolomics #####

meta <- select(dat.meta, Vial, Pig.id, Model, Time.from.EOS.min)

dat.box.plot <- dat_metab %>% 
  inner_join(meta, by = "Vial") %>% 
  relocate(Model, Time.from.EOS.min, .after = Vial) %>% 
  filter(Model == "REBOA Zone I + DCBI" | Model == "REBOA Zone III + DCBI"| Model == "DCBI") %>% 
  filter(Time.from.EOS.min == 60) 

dat.box.plot$Model <- factor(dat.box.plot$Model, 
                             levels = c("DCBI", "REBOA Zone I + DCBI", "REBOA Zone III + DCBI"))

analytes <- colnames(dat.box.plot[, c(4:ncol(dat.box.plot))])

analytes <- "Ornithine"
 
for (i in 1:length(analytes)) {
  
  a <- select(dat.box.plot, c("Model", analytes[i]))
  colnames(a) <- c("Model", "intensity")
    
  stat.test <- a %>%
    rstatix::dunn_test(intensity~Model, p.adjust.method = "fdr") %>% 
    rstatix::add_significance(p.col = "p.adj") %>% 
    rstatix::add_y_position(step.increase = 0.4) %>%
    filter(p.adj < 0.05) 
  
  print(analytes[i])
  kruskal.test(intensity ~ Model, data = a)
  
  boxplot <- ggplot(a, aes(x = Model, y = intensity))+
    stat_boxplot(geom = "errorbar", width = 0.5, lwd = 0.25) + 
    geom_boxplot(position="dodge", outlier.shape = 16, outlier.size = 0.5, width = 0.5, lwd = 0.25, aes(fill = Model) ) +
    
    scale_fill_manual(values = c(
                                  "#3B4992FF", #DCBI
                                  "#BB0021FF", #REBOA Zone I + DCBI
                                  "#008B45FF" #REBOA Zone III + DCBI
                                )) +
    
    #scale_y_continuous(breaks = breaks[[i]], limits = limits[[i]])+
    
    ggtitle(correct_metab_symbols[i,])+
    ggtitle("Ornithine")+
    
    labs(x=" ", y = " ", family = "Arial")+
    
    ggprism::add_pvalue(data = stat.test, 
                        bracket.size = 0.25,
                        label.size = 4, 
                        tip.length = 0)+
    
    #stat_pvalue_manual(data = stat.test, size = 10, hide.ns = TRUE)+
    
    #stat_compare_means(data=a,label.x = 1.2, label.y = label.y[i], family = "Arial", size = 5)+
    
    theme(plot.title = element_text(hjust = 0.5, size = 8),
          axis.title.x = element_text(size = 8),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.y = element_text(size = 6, color = "black"),
          axis.title.y = element_text(size = 8),
          axis.ticks.y = element_line(size = 0.25),
          axis.line = element_line(color = "black", size = 0.25),
          panel.background = element_rect(fill = "transparent",colour = NA),
          plot.background = element_rect(fill = "transparent",colour = NA),
          #panel.border = element_rect(color = "black", fill=NA, size = 0.25),
          legend.position = "none"
    )
  # view the plot
  #plot(boxplot)
  
  # Save the plot as PNG to working directory
  ggsave(boxplot,
         file = paste(Sys.Date(),analytes[i],"Figure5_metabolomics.tiff", sep = "_"),
         width = 1.5, height = 1.9, units = "in", dpi = 300, bg = "transparent")
  
}



#### Figure 6 ####

##### LY30 Data Frame Prep ######

dat_teg_ly30 <- dat.clin %>% 
  select( Pig.Number, Time..Min..from.End.of.Shock, CN.LY30...) %>% 
  dplyr::rename(CN_LY30 = CN.LY30..., Pig.id = Pig.Number, Time.from.EOS.min = Time..Min..from.End.of.Shock) %>% 
  filter_at(vars(CN_LY30), all_vars(!is.na(.))) %>% 
  left_join(dat_prot, by = c("Pig.id", "Time.from.EOS.min")) %>% 
  select(-c(4, 6:10)) %>%  
  na.omit() %>% 
  tidyr::pivot_longer(5:858, 
                      names_to = "Protein", 
                      values_to = "Raw_Intensity") %>%
  mutate(Norm_Intensity = (Raw_Intensity - mean(Raw_Intensity))/sd(Raw_Intensity)) %>% 
  tidyr::pivot_wider(id_cols = 1:4, names_from = Protein, values_from = Norm_Intensity) %>% 
  filter(Model %in% c("DCBI", "REBOA Zone I + DCBI", "REBOA Zone III + DCBI")) %>% 
  mutate(Time_Group = case_when(Time.from.EOS.min %in% c(0, 30, 60) ~ "Early Resuscitation",
                                Time.from.EOS.min %in% c(120, 180, 240) ~ "Late Resuscitation", 
                                TRUE ~ "Baseline")) %>% 
  relocate(Time_Group, .after = Time.from.EOS.min)


##### Panel A: LY30 Temporal Trend #####


# Function for data summary table used for all longitudinal line plots
data_summary <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(median = median(x[[col]], na.rm=TRUE),
      sd = std.error(x[[col]], na.rm=TRUE))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  data_sum <- rename(data_sum, c("median" = varname))
  return(data_sum)
}



dat.line.plot <- dat_teg_ly30

dat.line.plot$Model <- factor(dat.line.plot$Model, 
                              levels = c("DCBI", "REBOA Zone I + DCBI", "REBOA Zone III + DCBI"))



# a pipe that does everything you need
AOI <- dat.line.plot %>%
  select(c(CN_LY30, 'Time.from.EOS.min', 'Model', 'Pig.id')) %>% 
  rename_(QUANT = names(.)[1]) %>% 
  group_by(Pig.id) %>% 
  data_summary(varname = "QUANT",
               groupnames = c("Model", "Time.from.EOS.min")) %>%  
  filter(!is.na(Model)== TRUE)

p_1 <- ggplot(AOI, aes(x = Time.from.EOS.min, y = QUANT, group = Model, color = Model)) +
  #geom_errorbar(aes(ymin = QUANT - sd, ymax = QUANT + sd), width = 0.1) +
  
  geom_line(size = 0.5) + geom_point(size = 1) +
  
  scale_color_manual(values = c(
                                "#3B4992FF", #DCBI
                                "#BB0021FF", #REBOA Zone I + DCBI
                                "#008B45FF"
                                ))+
  
  labs(fill = "Method") + 
  
  scale_y_continuous(expand = c(0,0), breaks = seq(1, 2.5, by=0.5), limits = c(1,2.5))+
  
  scale_x_continuous(breaks = c(-120, 0, 30, 60, 120, 180, 240), 
                     labels = c("-120", "0", "30","60","120","180","240")) +
  xlab("") + ylab("") + ggtitle("CN-LY30") + #"abs(Median PG.Quantity)" 
  
  theme(legend.position = "none", 
        plot.margin = margin(t = 4, r = 10, b = 8, l = 8, unit = "pt"),
        plot.title = element_text(hjust = 0.5, size = 8), 
        axis.title = element_text(size = 8, color = 'black'), 
        axis.text.x = element_text(size = 6, vjust = 0.6, angle = 90, color = 'black'),
        axis.text.y = element_text(size = 6, color = 'black'),
        axis.ticks = element_line(size = 0.25),
        plot.background = element_rect(fill = "transparent",colour = NA),
        panel.background = element_rect(fill = NA),
        panel.border = element_rect(colour = "black", fill=NA, size = 0.25),
        #axis.line = element_line(color = 'black', size = 0.25, linetype = 'solid'),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())


plot(p_1)

# Save the plot to working directory
ggsave(p_1,
       file = paste(Sys.Date(),"Temporal_Trend_CN_LY30.svg", sep = "_"),
       width = 3, height = 2, units = "in", dpi = 300, bg = "transparent")







##### Pearson Correlation Combined Models #####

time <- c("Early Resuscitation", "Late Resuscitation")

proteins <- colnames(dat_teg_ly30[,5:ncol(dat_teg_ly30)])

names.vect <- c()

cor.val.vect <- c()

cor_df_list <- list()


for(g in 1:length(time)){
  
  # remove meta data
  df.first <- dat_teg_ly30 %>% 
    filter(Time_Group == time[g]) %>% 
    select(-Pig.id, -Time.from.EOS.min, -Time_Group, -Model)
  
  # initialize empty vector
  p.val.vect <- c()
  
  # loop through dataframe, calculating p value corresponding to the calculation 
  for(i in 1:ncol(df.first)){
    
    
    s <- df.first[,c(1,i)]
    s <- na.omit(s)
    s <- as.matrix(s)
    
    if(nrow(s)>1){
      
      a <- colnames(s)[2]
      names.vect[i] <- a
      
      # run the test
      cor.result <- stats::cor.test(s[,1], s[,2], method = "spearman")
      # extract the estimate
      cor.val <- cor.result$estimate
      # store the estimate
      cor.val.vect[i] <- cor.val
      
      # extract the p value
      p.val <- cor.result$p.value
      
      # append the p value to the vector
      p.val.vect[i] <- p.val
    }
  }
  
  ## Data coercion ##
  
  df.p <- data.frame(Protein = names.vect, P_Value = p.val.vect, SCorrelation = cor.val.vect)
  
  df.final <- df.p %>% 
    na.omit() %>% 
    slice(-1) %>% 
    mutate(NegLogP = -log10(P_Value), 
           Time_Group = time[g],
           
           Significance = case_when(P_Value < 0.05 & SCorrelation > 0 ~ "positive significant",
                                    P_Value < 0.05 & SCorrelation < 0 ~ "negative significant",
                                    TRUE ~ "insignificant"),
           Point_Color = case_when(P_Value < 0.05 ~ Time_Group,
                                   P_Value > 0.05 ~ ""))
  
  
  cor_df_list[[g]] <- df.final
  
}


cor_df <- bind_rows(cor_df_list)

cor_df$Point_Color <- factor(cor_df$Point_Color, levels = c("", "Early Resuscitation", "Late Resuscitation"))


cor_df %>% 
  filter(Significance == "positive significant") %>% 
  filter(Point_Color == "Early Resuscitation") %>%
  nrow()

cor_df %>% 
  filter(Significance == "positive significant") %>% 
  filter(Point_Color == "Late Resuscitation") %>%
  nrow()

cor_df %>% 
  filter(Significance == "negative significant") %>% 
  filter(Point_Color == "Early Resuscitation") %>%
  nrow()

cor_df %>% 
  filter(Significance == "negative significant") %>% 
  filter(Point_Color == "Late Resuscitation") %>%
  nrow()



###### Combined Models U-Plot ######

combined_model_plot_df <- cor_df %>% 
  mutate(Label = case_when(#NegLogP > -log10(0.05) & SCorrelation > 0 & Time_Group == "Early Resuscitation" ~ Protein, 
                           #NegLogP > -log10(0.05) & SCorrelation > 0 & Time_Group == "Late Resuscitation" ~ Protein,
                           #NegLogP > -log10(0.05) & SCorrelation < 0 & Time_Group == "Early Resuscitation" ~ Protein, 
                           NegLogP > -log10(0.05) & SCorrelation < 0 & Time_Group == "Late Resuscitation" ~ Protein,
                           TRUE ~ ""))

## Create scatterplot using the ggplot function ##
combined_model_plot <- ggplot(combined_model_plot_df, aes(SCorrelation, NegLogP, col = Point_Color, label = Label)) + 
  
  geom_point(size = 0.5) + 
  
  scale_color_manual(values = c("grey", # Insignificant 
                                "#01bdf8", # Early Resuscitation
                                "#000080" # Late Resuscitation
  ))+
  
  #scale_x_continuous(breaks = seq(-0.8,0.8,0.4), limits = c(-0.8, 0.8))+
  scale_y_continuous(expand = c(0,0), limits = c(1,5.5), breaks = seq(1,5,1))+
  
  xlab("Pearson Correlation") + 
  ylab("-log10(Pvalue)") + 
  
  
  geom_abline(intercept = -log10(0.05), slope = 0, col = 'grey57', lwd = 0.25, lty = 2) +
  
  theme(legend.position = "none", 
        plot.margin = margin(t = 15, r = 10, b = 8, l = 8, unit = "pt"),
        plot.title = element_text(hjust = 0.5, size = 8), 
        axis.title = element_text(size = 8), 
        axis.text.x = element_text(size = 6, vjust = 0.6, angle = 90),
        axis.text.y = element_text(size = 6),
        axis.ticks = element_line(size = 0.25),
        plot.background = element_rect(fill = "transparent",colour = NA),
        panel.background = element_rect(fill = 'white'),
        panel.border = element_rect(colour = "black", fill=NA, size = 0.25),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  
  geom_text_repel(size = 5,
                  max.overlaps = getOption("ggrepel.max.overlaps", default = 100),
                  box.padding = 0.7)

plot(combined_model_plot)

ggsave(paste(Sys.Date(), "Pearson_Combined_Models_Corr_Plot_LY30.svg"), bg = "transparent",
       plot = combined_model_plot, width = 2, height = 2.5, dpi = 300)




##### Correlations Combined Models & Time #####

proteins <- colnames(dat_teg_ly30[,6:ncol(dat_teg_ly30)])

names.vect <- c()

cor.val.vect <- c()


# remove meta data
df.first <- dat_teg_ly30 %>% 
  select(-Pig.id, -Time.from.EOS.min, -Time_Group, -Model)


# initialize empty vector
p.val.vect <- c()

# loop through dataframe, calculating p value corresponding to the calculation 
for(i in 1:ncol(df.first)){
  
  
  s <- df.first[,c(1,i)]
  s <- na.omit(s)
  s <- as.matrix(s)
  
  if(nrow(s)>1){
    
    a <- colnames(s)[2]
    names.vect[i] <- a
    
    # run the test
    cor.result <- stats::cor.test(s[,1], s[,2], method = "pearson")
    # extract the estimate
    cor.val <- cor.result$estimate
    # store the estimate
    cor.val.vect[i] <- cor.val
    
    # extract the p value
    p.val <- cor.result$p.value
    
    # append the p value to the vector
    p.val.vect[i] <- p.val
  }
}

## Data coercion ##

df.p <- data.frame(Protein = names.vect, P_Value = p.val.vect, SCorrelation = cor.val.vect)

cor_df <- df.p %>% 
  na.omit() %>% 
  slice(-1) %>% 
  mutate(NegLogP = -log10(P_Value), 
         Significance = case_when(P_Value < 0.05 & SCorrelation > 0 ~ "positive significant",
                                  P_Value < 0.05 & SCorrelation < 0 ~ "negative significant",
                                  TRUE ~ "insignificant"))


# count significant correlates

cor_df %>% 
  filter(Significance == "positive significant") %>% 
  nrow()

cor_df %>% 
  filter(Significance == "negative significant") %>% 
  nrow()



out <- cor_df[order(cor_df$Significance), ] 

write.csv(out,
          row.names = FALSE,
          paste(Sys.Date(), "LY30_Pearson_Correlation_Proteomics_Results.csv"))



###### Combined Models U-Plot ######

combined_model_plot_df <- cor_df %>% 
  mutate(Label = case_when(#NegLogP > 3 & SCorrelation > 0  ~ Protein,
                           NegLogP > -log10(0.05) & SCorrelation < 0  ~ Protein,
                           TRUE ~ ""))

## Create scatterplot using the ggplot function ##
combined_model_plot <- ggplot(combined_model_plot_df, aes(SCorrelation, NegLogP, col = Significance, label = Label)) + 
  
  geom_point(size = 0.5) + 
  
  scale_color_manual(values = c("grey", # Insignificant 
                                "#4d004d", # Early Resuscitation
                                "#daa520" # Late Resuscitation
  ))+
  
  #scale_x_continuous(breaks = seq(-0.8, 0.8, 0.2), limits = c(-0.7, 0.7))+
  scale_y_continuous(expand = c(0,0), limits = c(1,20.5), breaks = c(1,5,10,15,20))+
  
  xlab("Pearson Correlation") + 
  ylab("-log10(p-value)") + 
  
  
  geom_abline(intercept = -log10(0.05), slope = 0, col = 'grey57', lwd = 0.25, lty = 2) +
  
  theme(legend.position = "none", 
        plot.margin = margin(t = 15, r = 10, b = 8, l = 8, unit = "pt"),
        plot.title = element_text(hjust = 0.5, size = 8), 
        axis.title = element_text(size = 8), 
        axis.text.x = element_text(size = 6, vjust = 0.6, angle = 90),
        axis.text.y = element_text(size = 6),
        axis.ticks = element_line(size = 0.25),
        plot.background = element_rect(fill = "transparent",colour = NA),
        panel.background = element_rect(fill = 'white'),
        panel.border = element_rect(colour = "black", fill=NA, size = 0.25),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) 

geom_text_repel(size = 5,
                max.overlaps = getOption("ggrepel.max.overlaps", default = 200),
                box.padding = 0.7)



plot(combined_model_plot)

ggsave(paste(Sys.Date(), "Combined_Models_Combined_Time_Corr_Plot_LY30.svg"), bg = "transparent",
       plot = combined_model_plot, width = 2, height = 2.5, dpi = 300)





##### Pearson Correlation Per Model #####

time <- c("Early Resuscitation", "Late Resuscitation")

model_vector <- c("DCBI", "REBOA Zone I + DCBI", "REBOA Zone III + DCBI")

proteins <- colnames(dat_teg_ly30[,5:ncol(dat_teg_ly30)])

names.vect <- c()

cor.val.vect <- c()

cor_df_list <- list()

cor_df_list_holder <- list()


for(f in 1:length(model_vector)) {
  
model_select <- filter(dat_teg_ly30, Model == model_vector[f])

  for(g in 1:length(time)){
    
    
    q_select <- filter(model_select, Time_Group == time[g])
    
    
    # remove meta data
    df.first <- q_select[ ,-c(1,2,3,5)]
    
    # run spearman correlation
    #data.cor <- as.matrix(cor(df.first, method = "pearson")[,1])
    
    # rename column
    #colnames(data.cor) <- c("SCorrelation")
    
    # initialize empty vector
    p.val.vect <- c()
    
    # loop through dataframe, calculating p value corresponding to the calculation 
    for(i in 1:ncol(df.first)){
      
      
      s <- df.first[,c(1,i)]
      s <- na.omit(s)
      s <- as.matrix(s)
      
      if(nrow(s)>1){
        
        a <- colnames(s)[2]
        names.vect[i] <- a
        
        # run the test
        cor.result <- stats::cor.test(s[,1], s[,2], method = "pearson")
        # extract the estimate
        cor.val <- cor.result$estimate
        # store the estimate
        cor.val.vect[i] <- cor.val
        
        # extract the p value
        p.val <- cor.result$p.value
        
        # append the p value to the vector
        p.val.vect[i] <- p.val
      }
    }
    
    ## Data coercion ##
    
    df.p <- data.frame(Protein = names.vect, P_Value = p.val.vect, SCorrelation = cor.val.vect)
    
    df.final <- df.p %>% 
      na.omit() %>% 
      slice(-1) %>% 
      mutate(NegLogP = -log10(P_Value), 
             Time_Group = time[g],
             Model = model_vector[f],
             Significance = case_when(P_Value < 0.05 & SCorrelation > 0 ~ "positive significant",
                                      P_Value < 0.05 & SCorrelation < 0 ~ "negative significant",
                                      TRUE ~ "insignificant"),
             Point_Color = case_when(P_Value < 0.05 ~ Time_Group,
                                     P_Value > 0.05 ~ ""))
    
    
    cor_df_list[[g]] <- df.final
    
    #assign(paste("corr_df_condition", model[h], sep = "_"), df.final)
    #write.csv(corr1.df, paste("COMBAT_PLT_unit_corr_Group_", h, ".csv"))
  }
cor_df_list_holder[[f]] <- cor_df_list
}

cor_df <- bind_rows(cor_df_list_holder)

cor_df$Point_Color <- factor(cor_df$Point_Color, levels = c("", "Early Resuscitation", "Late Resuscitation"))



###### Panel C: NO AO U-Plot ######

NO_AO_df_plot <- cor_df %>% 
  filter(Model == "DCBI") %>% 
  mutate(Label = case_when(NegLogP > 2.8 & SCorrelation > 0 & Time_Group == "Early Resuscitation" ~ Protein, 
                           NegLogP > 2.4 & SCorrelation > 0 & Time_Group == "Late Resuscitation" ~ Protein,
                           NegLogP > 2.5 & SCorrelation < 0 & Time_Group == "Early Resuscitation" ~ Protein, 
                           NegLogP > 2.5 & SCorrelation < 0 & Time_Group == "Late Resuscitation" ~ Protein,
                           TRUE ~ ""))

## Create scatterplot using the ggplot function ##
NO_AO_corr_plot <- ggplot(NO_AO_df_plot, aes(SCorrelation, NegLogP, col = Point_Color, label = Label)) + 
  
  geom_point(size = 0.5) + 
  
  scale_color_manual(values = c("grey", # Insignificant 
                                "#01bdf8", # Early Resuscitation
                                "#000080" # Late Resuscitation
                     ))+
   
  scale_x_continuous(breaks = seq(-0.8,0.8,0.4), limits = c(-0.8, 0.8))+
  scale_y_continuous(expand = c(0,0), limits = c(1,4.5), breaks = seq(1,4.5,0.5))+
  
  xlab("Pearson Correlation") + 
  ylab("-log10(Pvalue)") + 
  
  
  geom_abline(intercept = -log10(0.05), slope = 0, col = 'grey57', lwd = 0.25, lty = 2) +
  
  theme(legend.position = "none", 
        plot.margin = margin(t = 15, r = 10, b = 8, l = 8, unit = "pt"),
        plot.title = element_text(hjust = 0.5, size = 8), 
        axis.title = element_text(size = 8), 
        axis.text.x = element_text(size = 6, vjust = 0.6, angle = 90),
        axis.text.y = element_text(size = 6),
        axis.ticks = element_line(size = 0.25),
        plot.background = element_rect(fill = "transparent",colour = NA),
        panel.background = element_rect(fill = 'white'),
        panel.border = element_rect(colour = "black", fill=NA, size = 0.25),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  
  geom_text_repel(size = 2,
                  max.overlaps = getOption("ggrepel.max.overlaps", default = 100),
                  box.padding = 0.7)

plot(NO_AO_corr_plot)

ggsave(paste(Sys.Date(), "NO_AO_Corr_Plot_LY30.svg"), bg = "transparent",
       plot = NO_AO_corr_plot, width = 2, height = 2.5, dpi = 300)


###### Panel E: Zone I U-Plot ######

ZoneI_df_plot <- cor_df %>% 
  filter(Model == "REBOA Zone I + DCBI") %>% 
  mutate(Label = case_when(NegLogP > 3.5 & SCorrelation > 0 & Time_Group == "Early Resuscitation" ~ Protein, 
                           NegLogP > 2.8 & SCorrelation > 0 & Time_Group == "Late Resuscitation" ~ Protein,
                           NegLogP > 2.8 & SCorrelation < 0 & Time_Group == "Early Resuscitation" ~ Protein, 
                           NegLogP > 2.0 & SCorrelation < 0 & Time_Group == "Late Resuscitation" ~ Protein,
                           TRUE ~ ""))

## Create scatterplot using the ggplot function ##
ZoneI_corr_plot <- ggplot(ZoneI_df_plot, aes(SCorrelation, NegLogP, col = Point_Color, label = Label)) + 
  
  geom_point(size = 0.5) + 
  
  scale_color_manual(values = c("grey", # Insignificant 
                                "#01bdf8", # Early Resuscitation
                                "#000080" # Late Resuscitation
  ))+
  
  scale_x_continuous(breaks = seq(-1, 1, 0.5), limits = c(-1, 1))+
  scale_y_continuous(expand = c(0,0), limits = c(1,8), breaks = seq(1, 8, 1))+
  
  xlab("Pearson Correlation") + 
  ylab("-log10(Pvalue)") + 
  
  
  geom_abline(intercept = -log10(0.05), slope = 0, col = 'grey57', lwd = 0.25, lty = 2) +
  
  theme(legend.position = "none", 
        plot.margin = margin(t = 15, r = 10, b = 8, l = 8, unit = "pt"),
        plot.title = element_text(hjust = 0.5, size = 8), 
        axis.title = element_text(size = 8), 
        axis.text.x = element_text(size = 6, vjust = 0.6, angle = 90),
        axis.text.y = element_text(size = 6),
        axis.ticks = element_line(size = 0.25),
        plot.background = element_rect(fill = "transparent",colour = NA),
        panel.background = element_rect(fill = 'white'),
        panel.border = element_rect(colour = "black", fill=NA, size = 0.25),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  
  geom_text_repel(size = 2,
                  max.overlaps = getOption("ggrepel.max.overlaps", default = 100),
                  box.padding = 0.7)



plot(ZoneI_corr_plot)

ggsave(paste(Sys.Date(), "ZoneI_Corr_Plot_LY30.svg"), bg = "transparent",
       plot = ZoneI_corr_plot, width = 2, height = 2.5, dpi = 300)



###### Panel E: Zone III U-Plot ######

ZoneIII_df_plot <- cor_df %>% 
  filter(Model == "REBOA Zone III + DCBI") %>% 
  mutate(Label = case_when(NegLogP > 3 & SCorrelation > 0 & Time_Group == "Early Resuscitation" ~ Protein, 
                           NegLogP > 3 & SCorrelation > 0 & Time_Group == "Late Resuscitation" ~ Protein,
                           NegLogP > 3 & SCorrelation < 0 & Time_Group == "Early Resuscitation" ~ Protein, 
                           NegLogP > 3 & SCorrelation < 0 & Time_Group == "Late Resuscitation" ~ Protein,
                           TRUE ~ ""))

## Create scatterplot using the ggplot function ##
ZoneIII_corr_plot <- ggplot(ZoneIII_df_plot, aes(SCorrelation, NegLogP, col = Point_Color, label = Label)) + 
  geom_point(size = 0.5) + 
  
  scale_color_manual(values = c("grey", # Insignificant 
                                "#01bdf8", # Early Resuscitation
                                "#000080" # Late Resuscitation
  ))+
  
  scale_x_continuous(breaks = seq(-1, 1, 0.5), limits = c(-1, 1))+
  scale_y_continuous(expand = c(0,0), limits = c(1,5), breaks = seq(1, 5, 0.5))+
  
  xlab("Pearson Correlation") + 
  ylab("-log10(Pvalue)") + 
  
  
  geom_abline(intercept = -log10(0.05), slope = 0, col = 'grey57', lwd = 0.25, lty = 2) +
  
  theme(legend.position = "none", 
        plot.margin = margin(t = 15, r = 10, b = 8, l = 8, unit = "pt"),
        plot.title = element_text(hjust = 0.5, size = 8), 
        axis.title = element_text(size = 8), 
        axis.text.x = element_text(size = 6, vjust = 0.6, angle = 90),
        axis.text.y = element_text(size = 6),
        axis.ticks = element_line(size = 0.25),
        plot.background = element_rect(fill = "transparent",colour = NA),
        panel.background = element_rect(fill = 'white'),
        panel.border = element_rect(colour = "black", fill=NA, size = 0.25),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  
  geom_text_repel(size = 2,
                  max.overlaps = getOption("ggrepel.max.overlaps", default = 100),
                  box.padding = 0.7)


plot(ZoneIII_corr_plot)

ggsave(paste(Sys.Date(), "ZoneIII_Corr_Plot_LY30.svg"), bg = "transparent",
       plot = ZoneIII_corr_plot, width = 2, height = 2.5, dpi = 300)








#### LY30 Result Export for GO Analysis ####

##### NO AO (DCBI) ##### 

# DCBI Early Resuscitation
out <- cor_df %>% 
  filter(Model == "DCBI") %>% 
  filter(Time_Group == "Early Resuscitation") %>%
  filter(Significance == "positive significant") %>% 
  select(Protein) 

write.csv(out,
          row.names = FALSE,
          "GO_Analysis_DCBI_Early_Resus_Pos_Sig_Corr.csv")

out <- cor_df %>% 
  filter(Model == "DCBI") %>% 
  filter(Time_Group == "Early Resuscitation") %>%
  filter( Significance == "negative significant") %>% 
  select(Protein)

write.csv(out,
          row.names = FALSE,
          "GO_Analysis_DCBI_Early_Resus_Neg_Sig_Corr.csv")

# DCBI Late Resuscitation
out <- cor_df %>% 
  filter(Model == "DCBI") %>% 
  filter(Time_Group == "Late Resuscitation") %>%
  filter(Significance == "positive significant") %>% 
  select(Protein) 

write.csv(out,
          row.names = FALSE,
          "GO_Analysis_DCBI_Late_Resus_Pos_Sig_Corr.csv")

out <- cor_df %>% 
  filter(Model == "DCBI") %>% 
  filter(Time_Group == "Late Resuscitation") %>%
  filter( Significance == "negative significant") %>% 
  select(Protein)

write.csv(out,
          row.names = FALSE,
          "GO_Analysis_DCBI_Late_Resus_Neg_Sig_Corr.csv")



##### REBOA Zone I ##### 

# REBOA Zone I Early Resuscitation
out <- cor_df %>% 
  filter(Model == "REBOA Zone I + DCBI") %>% 
  filter(Time_Group == "Early Resuscitation") %>%
  filter(Significance == "positive significant") %>% 
  select(Protein) 

write.csv(out,
          row.names = FALSE,
          "GO_Analysis_Zone_I_Early_Resus_Pos_Sig_Corr.csv")

out <- cor_df %>% 
  filter(Model == "REBOA Zone I + DCBI") %>% 
  filter(Time_Group == "Early Resuscitation") %>%
  filter( Significance == "negative significant") %>% 
  select(Protein)

write.csv(out,
          row.names = FALSE,
          "GO_Analysis_Zone_I_Early_Resus_Neg_Sig_Corr.csv")

# REBOA Zone I Late Resuscitation
out <- cor_df %>% 
  filter(Model == "REBOA Zone I + DCBI") %>% 
  filter(Time_Group == "Late Resuscitation") %>%
  filter(Significance == "positive significant") %>% 
  select(Protein) 

write.csv(out,
          row.names = FALSE,
          "GO_Analysis_Zone_I_Late_Resus_Pos_Sig_Corr.csv")

out <- cor_df %>% 
  filter(Model == "REBOA Zone I + DCBI") %>% 
  filter(Time_Group == "Late Resuscitation") %>%
  filter( Significance == "negative significant") %>% 
  select(Protein)

write.csv(out,
          row.names = FALSE,
          "GO_Analysis_Zone_I_Late_Resus_Neg_Sig_Corr.csv")



##### REBOA Zone III ##### 

# REBOA Zone III Early Resuscitation
out <- cor_df %>% 
  filter(Model == "REBOA Zone III + DCBI") %>% 
  filter(Time_Group == "Early Resuscitation") %>%
  filter(Significance == "positive significant") %>% 
  select(Protein) 

write.csv(out,
          row.names = FALSE,
          "GO_Analysis_ZONE_III_Early_Resus_Pos_Sig_Corr.csv")

out <- cor_df %>% 
  filter(Model == "REBOA Zone III + DCBI") %>% 
  filter(Time_Group == "Early Resuscitation") %>%
  filter( Significance == "negative significant") %>% 
  select(Protein)

write.csv(out,
          row.names = FALSE,
          "GO_Analysis_ZONE_III_Early_Resus_Neg_Sig_Corr.csv")

# REBOA Zone III Late Resuscitation
out <- cor_df %>% 
  filter(Model == "REBOA Zone III + DCBI") %>% 
  filter(Time_Group == "Late Resuscitation") %>%
  filter(Significance == "positive significant") %>% 
  select(Protein) 

write.csv(out,
          row.names = FALSE,
          "GO_Analysis_ZONE_III_Late_Resus_Pos_Sig_Corr.csv")

out <- cor_df %>% 
  filter(Model == "REBOA Zone III + DCBI") %>% 
  filter(Time_Group == "Late Resuscitation") %>%
  filter( Significance == "negative significant") %>% 
  select(Protein)

write.csv(out,
          row.names = FALSE,
          "GO_Analysis_ZONE_III_Late_Resus_Neg_Sig_Corr.csv")


#### LY30 GO Bar Plots #####

##### Read In GO Result Data Frames ######
# created through gene ontology portal
# used the sus scrofa data base

## DCBI

# model DCBI Early Resuscitation , Negative Correlates 
dat_go_01 <- read.csv("GO_DCBI_Early_Resus_Neg_Sig.csv") %>% 
  mutate(Resus_Phase = "Early")

# model DCBI Early Resuscitation , Positive Correlates 
dat_go_02 <- read.csv("GO_DCBI_Early_Resus_Pos_Sig.csv") %>% 
  mutate(Resus_Phase = "Early")

# model DCBI Late Resuscitation , Negative Correlates 
dat_go_03 <- read.csv("GO_DCBI_Late_Resus_Pos_Sig.csv") %>% 
  mutate(Resus_Phase = "Late")


# model DCBI Late Resuscitation , Positive Correlates 
dat_go_04 <- read.csv("GO_DCBI_Late_Resus_Pos_Sig.csv") %>% 
  mutate(Resus_Phase = "Late")


# combine to use in bar plot
DCBI_go_neg <- rbind(dat_go_01, dat_go_03)
DCBI_go_pos <- rbind(dat_go_02, dat_go_04)



## REBOA Zone I

# model Zone I Early Resuscitation , Negative Correlates 
dat_go_05 <- read.csv("GO_ZoneI_Early_Neg_Sig.csv") %>% 
  mutate(Resus_Phase = "Early") %>% 
  mutate(Fold_Enrichment = -1*Fold_Enrichment)

# model DCBI Early Resuscitation , Positive Correlates 
dat_go_06 <- read.csv("GO_ZoneI_Early_Pos_Sig.csv") %>% 
  mutate(Resus_Phase = "Early")

# model Zone I Late Resuscitation , Negative Correlates 
dat_go_07 <- read.csv("GO_ZoneI_Late_Neg_Sig.csv") %>% 
  mutate(Resus_Phase = "Late") %>% 
  mutate(Fold_Enrichment = -1*Fold_Enrichment)


# model Zone I Late Resuscitation , Positive Correlates 
dat_go_08 <- read.csv("GO_ZoneI_Late_Pos_Sig.csv") %>% 
  mutate(Resus_Phase = "Late")


# combine to use in bar plot
ZoneI_go_neg <- rbind(dat_go_05, dat_go_07, dat_go_06, dat_go_08)
ZoneI_go_pos <- rbind(dat_go_06, dat_go_08)



## REBOA Zone III

# model Zone I Early Resuscitation , Negative Correlates 
dat_go_09 <- read.csv("GO_ZoneIII_Early_Neg_Sig.csv") %>% 
  mutate(Resus_Phase = "Early")

# model DCBI Early Resuscitation , Positive Correlates 
dat_go_10 <- read.csv("GO_ZoneIII_Early_Pos_Sig.csv") %>% 
  mutate(Resus_Phase = "Early")

# model Zone I Late Resuscitation , Negative Correlates 
dat_go_11 <- read.csv("GO_ZoneIII_Late_Neg_Sig.csv") %>% 
  mutate(Resus_Phase = "Late")


# model Zone I Late Resuscitation , Positive Correlates 
dat_go_12 <- read.csv("GO_ZoneIII_Late_Pos_Sig.csv") %>% 
  mutate(Resus_Phase = "Late")


# combine to use in bar plot
ZoneIII_go_neg <- rbind(dat_go_09, dat_go_11)
ZoneIII_go_pos <- rbind(dat_go_10, dat_go_12)



##### DCBI Negative Correlates #####

## Create Bar plot using the ggplot function ##
DCBI_neg_go_plot <- ggplot(DCBI_go_neg, aes(x = GO_BP, y = Fold_Enrichment, fill = Resus_Phase)) + 
  
  geom_bar(stat = "identity", width = 0.5) + 
  coord_flip() +

  #scale_y_continuous(expand = c(0,0), limits = c(0,100), breaks = seq(0,100,25))+
  
  xlab("") + 
  ylab("Fold Enrichment") + 
  
  scale_fill_manual(values = c(
                                "#01bdf8", # Early Resuscitation
                                "#000080" # Late Resuscitation
                                
                                ))+
  
  theme(legend.position = "none", 
        plot.margin = margin(t = 15, r = 10, b = 8, l = 8, unit = "pt"),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 16), 
        axis.title = element_text(size = 8, color = 'black'), 
        axis.text.x = element_text(size = 6, vjust = 0.6, angle = 90, color = 'black'),
        axis.text.y = element_text(size = 6, color = 'black'),
        axis.ticks = element_line(size = 0.25),
        plot.background = element_rect(fill = "transparent",colour = NA),
        panel.background = element_rect(fill = NA),
        #panel.border = element_rect(colour = "black", fill=NA),
        axis.line = element_line(color = 'black', size = 0.25, linetype = 'solid'),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
  

plot(DCBI_neg_go_plot)

ggsave(paste(Sys.Date(), "DCBI_Neg_Corr_GO_Plot.svg"), bg = "transparent",
       plot = DCBI_neg_go_plot, width = 3, height = 2.5, dpi = 300)


##### DCBI Positive Correlates #####

## Create Bar plot using the ggplot function ##
DCBI_pos_go_plot <- ggplot(DCBI_go_pos, aes(x = GO_BP, y = Fold_Enrichment, fill = Resus_Phase)) + 
  
  geom_bar(stat = "identity", width = 0.5) + 
  coord_flip() +
  
  scale_y_continuous(expand = c(0,0), limits = c(0,100), breaks = seq(0,100,25))+
  
  xlab("") + 
  ylab("Fold Enrichment") + 
  
  scale_fill_manual(values = c(
    "#01bdf8", # Early Resuscitation
    "#000080" # Late Resuscitation
    
  ))+
  
  theme(legend.position = "none", 
        plot.margin = margin(t = 15, r = 10, b = 8, l = 8, unit = "pt"),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 16), 
        axis.title = element_text(size = 8, color = 'black'), 
        axis.text.x = element_text(size = 6, vjust = 0.6, angle = 90, color = 'black'),
        axis.text.y = element_text(size = 6, color = 'black'),
        axis.ticks = element_line(size = 0.25),
        plot.background = element_rect(fill = "transparent",colour = NA),
        panel.background = element_rect(fill = NA),
        #panel.border = element_rect(colour = "black", fill=NA),
        axis.line = element_line(color = 'black', size = 0.25, linetype = 'solid'),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())


plot(DCBI_pos_go_plot)

ggsave(paste(Sys.Date(), "DCBI_Pos_Corr_GO_Plot.svg"), bg = "transparent",
       plot = DCBI_pos_go_plot, width = 3, height = 2.5, dpi = 300)







##### ZoneI Negative Correlates #####

## Create Bar plot using the ggplot function ##
ZoneI_neg_go_plot <- ggplot(ZoneI_go_neg, aes(x = GO_BP, y = Fold_Enrichment, fill = Resus_Phase)) + 
  
  geom_bar(stat = "identity", width = 0.5) + 
  coord_flip() +
  
  geom_hline(yintercept=0, color = "black", size = 0.25)+
  
  #scale_y_continuous(expand = c(0,0), limits = c(0,100), breaks = seq(0,100,25))+
  
  xlab("") + 
  ylab("Fold Enrichment") + 
  
  scale_fill_manual(values = c(
    "#01bdf8", # Early Resuscitation
    "#000080" # Late Resuscitation
    
  ))+
  
  theme(legend.position = "none", 
        plot.margin = margin(t = 15, r = 10, b = 8, l = 8, unit = "pt"),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 16), 
        axis.title = element_text(size = 8, color = 'black'), 
        axis.text.x = element_text(size = 6, vjust = 0.6, angle = 90, color = 'black'),
        axis.text.y = element_text(size = 6, color = 'black'),
        axis.ticks = element_line(size = 0.25),
        plot.background = element_rect(fill = "transparent",colour = NA),
        panel.background = element_rect(fill = NA),
        #panel.border = element_rect(colour = "black", fill=NA),
        axis.line = element_line(color = 'black', size = 0.25, linetype = 'solid'),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())


plot(ZoneI_neg_go_plot)

ggsave(paste(Sys.Date(), "ZoneI_Neg_Corr_GO_Plot.svg"), bg = "transparent",
       plot = ZoneI_neg_go_plot, width = 3, height = 2.5, dpi = 300)


##### ZoneI Positive Correlates #####

## Create Bar plot using the ggplot function ##
ZoneI_pos_go_plot <- ggplot(ZoneI_go_pos, aes(x = GO_BP, y = Fold_Enrichment, fill = Resus_Phase)) + 
  
  geom_bar(stat = "identity", width = 0.5) + 
  coord_flip() +
  
  scale_y_continuous(expand = c(0,0), limits = c(0,100), breaks = seq(0,100,25))+
  
  xlab("") + 
  ylab("Fold Enrichment") + 
  
  scale_fill_manual(values = c(
    "#01bdf8", # Early Resuscitation
    "#000080" # Late Resuscitation
    
  ))+
  
  theme(legend.position = "none", 
        plot.margin = margin(t = 15, r = 10, b = 8, l = 8, unit = "pt"),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 16), 
        axis.title = element_text(size = 8, color = 'black'), 
        axis.text.x = element_text(size = 6, vjust = 0.6, angle = 90, color = 'black'),
        axis.text.y = element_text(size = 6, color = 'black'),
        axis.ticks = element_line(size = 0.25),
        plot.background = element_rect(fill = "transparent",colour = NA),
        panel.background = element_rect(fill = NA),
        #panel.border = element_rect(colour = "black", fill=NA),
        axis.line = element_line(color = 'black', size = 0.25, linetype = 'solid'),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())


plot(ZoneI_pos_go_plot)

ggsave(paste(Sys.Date(), "ZoneI_Pos_Corr_GO_Plot.svg"), bg = "transparent",
       plot = ZoneI_pos_go_plot, width = 3, height = 2.5, dpi = 300)





##### ZoneIII Negative Correlates #####

## Create Bar plot using the ggplot function ##
ZoneIII_neg_go_plot <- ggplot(ZoneIII_go_neg, aes(x = GO_BP, y = Fold_Enrichment, fill = Resus_Phase)) + 
  
  geom_bar(stat = "identity", width = 0.5) + 
  coord_flip() +
  
  scale_y_continuous(expand = c(0,0), limits = c(0,100), breaks = seq(0,100,25))+
  
  xlab("") + 
  ylab("Fold Enrichment") + 
  
  scale_fill_manual(values = c(
    "#01bdf8", # Early Resuscitation
    "#000080" # Late Resuscitation
    
  ))+
  
  theme(legend.position = "none", 
        plot.margin = margin(t = 15, r = 10, b = 8, l = 8, unit = "pt"),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 16), 
        axis.title = element_text(size = 8, color = 'black'), 
        axis.text.x = element_text(size = 6, vjust = 0.6, angle = 90, color = 'black'),
        axis.text.y = element_text(size = 6, color = 'black'),
        axis.ticks = element_line(size = 0.25),
        plot.background = element_rect(fill = "transparent",colour = NA),
        panel.background = element_rect(fill = NA),
        #panel.border = element_rect(colour = "black", fill=NA),
        axis.line = element_line(color = 'black', size = 0.25, linetype = 'solid'),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())


plot(ZoneIII_neg_go_plot)

ggsave(paste(Sys.Date(), "ZoneIII_Neg_Corr_GO_Plot.svg"), bg = "transparent",
       plot = ZoneIII_neg_go_plot, width = 3, height = 2.5, dpi = 300)


##### ZoneIII Positive Correlates #####

## Create Bar plot using the ggplot function ##
ZoneIII_pos_go_plot <- ggplot(ZoneIII_go_pos, aes(x = GO_BP, y = Fold_Enrichment, fill = Resus_Phase)) + 
  
  geom_bar(stat = "identity", width = 0.5) + 
  coord_flip() +
  
  scale_y_continuous(expand = c(0,0), limits = c(0,100), breaks = seq(0,100,25))+
  
  xlab("") + 
  ylab("Fold Enrichment") + 
  
  scale_fill_manual(values = c(
    "#01bdf8", # Early Resuscitation
    "#000080" # Late Resuscitation
    
  ))+
  
  theme(legend.position = "none", 
        plot.margin = margin(t = 15, r = 10, b = 8, l = 8, unit = "pt"),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 16), 
        axis.title = element_text(size = 8, color = 'black'), 
        axis.text.x = element_text(size = 6, vjust = 0.6, angle = 90, color = 'black'),
        axis.text.y = element_text(size = 6, color = 'black'),
        axis.ticks = element_line(size = 0.25),
        plot.background = element_rect(fill = "transparent",colour = NA),
        panel.background = element_rect(fill = NA),
        #panel.border = element_rect(colour = "black", fill=NA),
        axis.line = element_line(color = 'black', size = 0.25, linetype = 'solid'),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())


plot(ZoneIII_pos_go_plot)

ggsave(paste(Sys.Date(), "ZoneIII_Pos_Corr_GO_Plot.svg"), bg = "transparent",
       plot = ZoneIII_pos_go_plot, width = 3, height = 2.5, dpi = 300)







##### LY30 Zone I vs Zone III Early, Pos Corr #####

dat_1 <- read.csv("ZoneI_vs_ZoneIII_Early_Pos.csv") 


dat_1$Model <- factor(dat_1$GO_BP)


## Create Bar plot using the ggplot function ##
pos_go_plot <- ggplot(dat_1, aes(x = Order, y = Log2_FC, fill = Resus_Phase)) + 
  
  geom_bar(stat = "identity", width = 0.5) + 
  coord_flip() +
  
  geom_hline(yintercept = 0, color = "black", size = 0.25)+
  
  scale_x_discrete(labels = dat_1$GO_BP)+
  scale_y_continuous(expand = c(0,0), limits = c(-15,15), breaks = seq(-15,15,5))+
  
  xlab("") + 
  ylab("Log2(Fold Change)") + 
  
  scale_fill_manual(values = c(
    "#01bdf8", # Early Resuscitation
    "#000080" # Late Resuscitation
    
  ))+
  
  theme(legend.position = "none", 
        plot.margin = margin(t = 15, r = 10, b = 8, l = 8, unit = "pt"),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 16), 
        axis.title = element_text(size = 8, color = 'black'), 
        axis.text.x = element_text(size = 6, vjust = 0.6, angle = 90, color = 'black'),
        axis.text.y = element_text(size = 6, color = 'black'),
        #axis.text.y = element_blank(),
        axis.ticks = element_line(size = 0.25),
        #axis.ticks.y = element_blank(),
        plot.background = element_rect(fill = "transparent",colour = NA),
        panel.background = element_rect(fill = NA),
        panel.border = element_rect(colour = "black", fill=NA, size = 0.25),
        axis.line = element_line(color = 'black', size = 0, linetype = 'solid'),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())


plot(pos_go_plot)

ggsave(paste(Sys.Date(), "Pos_Corr_GO_Plot.svg"), bg = "transparent",
       plot = pos_go_plot, width = 4, height = 1.5, dpi = 300)






##### LY30 Zone I vs Zone III Early, Neg Corr #####

dat_2 <- read.csv("ZoneI_vs_ZoneIII_neg_corr.csv") 


dat_2$Model <- factor(dat_2$GO_BP)


## Create Bar plot using the ggplot function ##
pos_go_plot <- ggplot(dat_2, aes(x = Order, y = Log2_FC, fill = Resus_Phase)) + 
  
  geom_bar(stat = "identity", width = 0.5) + 
  coord_flip() +
  
  geom_hline(yintercept = 0, color = "black", size = 0.25)+
  
  scale_x_discrete(labels = dat_2$GO_BP)+
  scale_y_continuous(expand = c(0,0), limits = c(-15,15), breaks = seq(-15,15,5))+
  
  xlab("") + 
  ylab("Log2(Fold Change)") + 
  
  scale_fill_manual(values = c(
    "#01bdf8", # Early Resuscitation
    "#000080" # Late Resuscitation
    
  ))+
  
  theme(legend.position = "none", 
        plot.margin = margin(t = 15, r = 10, b = 8, l = 32, unit = "pt"),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 16), 
        axis.title = element_text(size = 8, color = 'black'), 
        axis.text.x = element_text(size = 6, vjust = 0.6, angle = 90, color = 'black'),
        axis.text.y = element_text(size = 6, color = 'black'),
        #axis.text.y = element_blank(),
        axis.ticks = element_line(size = 0.25),
        #axis.ticks.y = element_blank(),
        plot.background = element_rect(fill = "transparent",colour = NA),
        panel.background = element_rect(fill = NA),
        panel.border = element_rect(colour = "black", fill=NA, size = 0.25),
        axis.line = element_line(color = 'black', size = 0, linetype = 'solid'),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())


plot(pos_go_plot)

ggsave(paste(Sys.Date(), "Neg_Corr_GO_Plot.svg"), bg = "transparent",
       plot = pos_go_plot, width = 4, height = 1.5, dpi = 300)





#### Proteomics Pearson Correlation MA #####

##### MA data frame prep #####

dat_teg_MA <- dat.clin %>% 
  select( Pig.Number, Time..Min..from.End.of.Shock, CN.MA.mm.) %>% 
  dplyr::rename(CN_MA = CN.MA.mm., Pig.id = Pig.Number, Time.from.EOS.min = Time..Min..from.End.of.Shock) %>% 
  filter_at(vars(CN_MA), all_vars(!is.na(.))) %>% 
  left_join(dat_prot, by = c("Pig.id", "Time.from.EOS.min")) %>% 
  select(-c(4, 6:10)) %>%  
  na.omit() %>% 
  tidyr::pivot_longer(5:858, 
                      names_to = "Protein", 
                      values_to = "Raw_Intensity") %>%
  mutate(Norm_Intensity = (Raw_Intensity - mean(Raw_Intensity))/sd(Raw_Intensity)) %>% 
  tidyr::pivot_wider(id_cols = 1:4, names_from = Protein, values_from = Norm_Intensity) %>% 
  filter(Model %in% c("DCBI", "REBOA Zone I + DCBI", "REBOA Zone III + DCBI")) %>% 
  mutate(Time_Group = case_when(Time.from.EOS.min %in% c(0, 30, 60) ~ "Early Resuscitation",
                                Time.from.EOS.min %in% c(120, 180, 240) ~ "Late Resuscitation", 
                                TRUE ~ "Baseline")) %>% 
  relocate(Time_Group, .after = Time.from.EOS.min)


##### MA Temporal Trend #####
# Function for data summary table used for all longitudinal line plots
data_summary <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(median = median(x[[col]], na.rm=TRUE),
      sd = std.error(x[[col]], na.rm=TRUE))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  data_sum <- rename(data_sum, c("median" = varname))
  return(data_sum)
}



dat.line.plot <- dat_teg_MA

dat.line.plot$Model <- factor(dat.line.plot$Model, 
                              levels = c("DCBI", "REBOA Zone I + DCBI", "REBOA Zone III + DCBI"))



# a pipe that does everything you need
AOI <- dat.line.plot %>%
  select(c(CN_MA, 'Time.from.EOS.min', 'Model', 'Pig.id')) %>% 
  rename_(QUANT = names(.)[1]) %>% 
  group_by(Pig.id) %>% 
  data_summary(varname = "QUANT",
               groupnames = c("Model", "Time.from.EOS.min")) %>%  
  filter(!is.na(Model)== TRUE)

p_1 <- ggplot(AOI, aes(x = Time.from.EOS.min, y = QUANT, group = Model, color = Model)) +
  #geom_errorbar(aes(ymin = QUANT - sd, ymax = QUANT + sd), width = 0.1) +
  
  geom_line(size = 0.5) + geom_point(size = 1) +
  scale_color_manual(values = c(
                                "#3B4992FF", #DCBI
                                "#BB0021FF", #REBOA Zone I + DCBI
                                "#008B45FF"))+
  labs(fill = "Method") + 
  
  scale_y_continuous(expand = c(0,0), breaks = seq(60,80, by=5), limits = c(60,80))+
  
  scale_x_continuous(breaks = c(-120, 0, 30, 60, 120, 180, 240), 
                     labels = c("-120", "0", "30","60","120","180","240")) +
  xlab("Time from EOS (min)") + ylab("MA(mm)") + ggtitle(" ") + #"abs(Median PG.Quantity)" 
  
  theme(legend.position = "none", 
        plot.margin = margin(t = 4, r = 10, b = 8, l = 8, unit = "pt"),
        plot.title = element_text(hjust = 0.5, size = 8), 
        axis.title = element_text(size = 6, color = 'black'), 
        axis.text.x = element_text(size = 6, vjust = 0.6, angle = 90, color = 'black'),
        axis.text.y = element_text(size = 6, color = 'black'),
        axis.ticks = element_line(size = 0.25),
        plot.background = element_rect(fill = "transparent",colour = NA),
        panel.background = element_rect(fill = NA),
        panel.border = element_rect(colour = "black", fill=NA, size = 0.25),
        #axis.line = element_line(color = 'black', size = 0.25, linetype = 'solid'),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())


plot(p_1)

# Save the plot to working directory
ggsave(p_1,
       file = paste(Sys.Date(),"Temporal_Trend_CN_MA.svg", sep = "_"),
       width = 2.2, height = 2, units = "in", dpi = 300, bg = "transparent")



##### Correlations Combined Models #####

time <- c("Early Resuscitation", "Late Resuscitation")

roteins <- colnames(dat_teg_MA[,6:ncol(dat_teg_MA)])

names.vect <- c()

cor.val.vect <- c()

cor_df_list <- list()

  
for(g in 1:length(time)){
  
  # remove meta data
  df.first <- dat_teg_MA %>% 
    filter(Time_Group == time[g]) %>% 
    select(-Pig.id, -Time.from.EOS.min, -Time_Group, -Model)
  
  
  # initialize empty vector
  p.val.vect <- c()
  
  # loop through dataframe, calculating p value corresponding to the calculation 
  for(i in 1:ncol(df.first)){
    
    
    s <- df.first[,c(1,i)]
    s <- na.omit(s)
    s <- as.matrix(s)
    
    if(nrow(s)>1){
      
      a <- colnames(s)[2]
      names.vect[i] <- a
      
      # run the test
      cor.result <- stats::cor.test(s[,1], s[,2], method = "pearson")
      # extract the estimate
      cor.val <- cor.result$estimate
      # store the estimate
      cor.val.vect[i] <- cor.val
      
      # extract the p value
      p.val <- cor.result$p.value
      
      # append the p value to the vector
      p.val.vect[i] <- p.val
    }
  }
  
  ## Data coercion ##
  
  df.p <- data.frame(Protein = names.vect, P_Value = p.val.vect, SCorrelation = cor.val.vect)
  
  df.final <- df.p %>% 
    na.omit() %>% 
    slice(-1) %>% 
    mutate(NegLogP = -log10(P_Value), 
           Time_Group = time[g],
           Significance = case_when(P_Value < 0.05 & SCorrelation > 0 ~ "positive significant",
                                    P_Value < 0.05 & SCorrelation < 0 ~ "negative significant",
                                    TRUE ~ "insignificant"),
           Point_Color = case_when(P_Value < 0.05 ~ Time_Group,
                                   P_Value > 0.05 ~ ""))
  
  
  cor_df_list[[g]] <- df.final
  
  #assign(paste("corr_df_condition", model[h], sep = "_"), df.final)
  #write.csv(corr1.df, paste("COMBAT_PLT_unit_corr_Group_", h, ".csv"))
}

cor_df <- bind_rows(cor_df_list)

cor_df$Point_Color <- factor(cor_df$Point_Color, levels = c("", "Early Resuscitation", "Late Resuscitation"))


cor_df %>% 
  filter(Significance == "positive significant") %>% 
  filter(Point_Color == "Early Resuscitation") %>%
  nrow()

cor_df %>% 
  filter(Significance == "positive significant") %>% 
  filter(Point_Color == "Late Resuscitation") %>%
  nrow()

cor_df %>% 
  filter(Significance == "negative significant") %>% 
  filter(Point_Color == "Early Resuscitation") %>%
  nrow()

cor_df %>% 
  filter(Significance == "negative significant") %>% 
  filter(Point_Color == "Late Resuscitation") %>%
  nrow()
  

###### Combined Models U-Plot ######

combined_model_plot_df <- cor_df %>% 
  mutate(Label = case_when(NegLogP > 1 & SCorrelation > 0 & Time_Group == "Early Resuscitation" ~ Protein, 
                           NegLogP > 1 & SCorrelation > 0 & Time_Group == "Late Resuscitation" ~ Protein,
                           NegLogP > 2.2 & SCorrelation < 0 & Time_Group == "Early Resuscitation" ~ Protein, 
                           NegLogP > 1 & SCorrelation < 0 & Time_Group == "Late Resuscitation" ~ Protein,
                           TRUE ~ ""))

## Create scatterplot using the ggplot function ##
combined_model_plot <- ggplot(combined_model_plot_df, aes(SCorrelation, NegLogP, col = Point_Color, label = Label)) + 
  
  geom_point(size = 0.5) + 
  
  scale_color_manual(values = c("grey", # Insignificant 
                                "#01bdf8", # Early Resuscitation
                                "#000080" # Late Resuscitation
  ))+
  
  #scale_x_continuous(breaks = seq(-1, 1, 0.5), limits = c(-1, 1))+
  #scale_y_continuous(expand = c(0,0), limits = c(1,6.5), breaks = seq(1, 6.5, 0.5))+
  
  xlab("Pearson Correlation") + 
  ylab("-log10(Pvalue)") + 
  
  
  geom_abline(intercept = -log10(0.05), slope = 0, col = 'grey57', lwd = 0.25, lty = 2) +
  
  theme(legend.position = "none", 
        plot.margin = margin(t = 15, r = 10, b = 8, l = 8, unit = "pt"),
        plot.title = element_text(hjust = 0.5, size = 8), 
        axis.title = element_text(size = 8), 
        axis.text.x = element_text(size = 6, vjust = 0.6, angle = 90),
        axis.text.y = element_text(size = 6),
        axis.ticks = element_line(size = 0.25),
        plot.background = element_rect(fill = "transparent",colour = NA),
        panel.background = element_rect(fill = 'white'),
        panel.border = element_rect(colour = "black", fill=NA, size = 0.25),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  
  geom_text_repel(size = 5,
                  max.overlaps = getOption("ggrepel.max.overlaps", default = 200),
                  box.padding = 0.7)



plot(combined_model_plot)

ggsave(paste(Sys.Date(), "Combined_Models_Corr_Plot_MA.svg"), bg = "transparent",
       plot = combined_model_plot, width = 2, height = 2.5, dpi = 300)








##### Correlations Combined Models & Time #####

proteins <- colnames(dat_teg_MA[,6:ncol(dat_teg_MA)])

names.vect <- c()

cor.val.vect <- c()


# remove meta data
df.first <- dat_teg_MA %>% 
  select(-Pig.id, -Time.from.EOS.min, -Time_Group, -Model)


# initialize empty vector
p.val.vect <- c()

# loop through dataframe, calculating p value corresponding to the calculation 
for(i in 1:ncol(df.first)){
  
  
  s <- df.first[,c(1,i)]
  s <- na.omit(s)
  s <- as.matrix(s)
  
  if(nrow(s)>1){
    
    a <- colnames(s)[2]
    names.vect[i] <- a
    
    # run the test
    cor.result <- stats::cor.test(s[,1], s[,2], method = "pearson")
    # extract the estimate
    cor.val <- cor.result$estimate
    # store the estimate
    cor.val.vect[i] <- cor.val
    
    # extract the p value
    p.val <- cor.result$p.value
    
    # append the p value to the vector
    p.val.vect[i] <- p.val
  }
}

## Data coercion ##

df.p <- data.frame(Protein = names.vect, P_Value = p.val.vect, SCorrelation = cor.val.vect)

cor_df <- df.p %>% 
  na.omit() %>% 
  slice(-1) %>% 
  mutate(NegLogP = -log10(P_Value), 
         Significance = case_when(P_Value < 0.05 & SCorrelation > 0 ~ "positive significant",
                                  P_Value < 0.05 & SCorrelation < 0 ~ "negative significant",
                                  TRUE ~ "insignificant"))

out <- cor_df[order(cor_df$Significance), ] 

write.csv(out,
          row.names = FALSE,
          paste(Sys.Date(), "MA_Pearson_Correlation_Results.csv"))


# count significant correlates

cor_df %>% 
  filter(Significance == "positive significant") %>% 
  nrow()

cor_df %>% 
  filter(Significance == "negative significant") %>% 
  nrow()




###### Combined Models U-Plot ######

combined_model_plot_df <- cor_df %>% 
  mutate(Label = case_when(NegLogP > 3 & SCorrelation > 0  ~ Protein,
                           NegLogP > 10 & SCorrelation < 0  ~ Protein,
                           TRUE ~ ""))

## Create scatterplot using the ggplot function ##
combined_model_plot <- ggplot(combined_model_plot_df, aes(SCorrelation, NegLogP, col = Significance, label = Label)) + 
  
  geom_point(size = 0.5) + 
  
  scale_color_manual(values = c("grey", # Insignificant 
                                "#daa520", # Early Resuscitation
                                "#daa520" # Late Resuscitation
  ))+
  
  scale_x_continuous(breaks = seq(-0.8, 0.8, 0.2), limits = c(-0.7, 0.7))+
  scale_y_continuous(expand = c(0,0), limits = c(1,20), breaks = c(1,5,10,15,20))+
  
  xlab("Pearson Correlation") + 
  ylab("-log10(p-value)") + 
  
  
  geom_abline(intercept = -log10(0.05), slope = 0, col = 'grey57', lwd = 0.25, lty = 2) +
  
  theme(legend.position = "none", 
        plot.margin = margin(t = 15, r = 10, b = 8, l = 8, unit = "pt"),
        plot.title = element_text(hjust = 0.5, size = 8), 
        axis.title = element_text(size = 8), 
        axis.text.x = element_text(size = 6, vjust = 0.6, angle = 90),
        axis.text.y = element_text(size = 6),
        axis.ticks = element_line(size = 0.25),
        plot.background = element_rect(fill = "transparent",colour = NA),
        panel.background = element_rect(fill = 'white'),
        panel.border = element_rect(colour = "black", fill=NA, size = 0.25),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  
  geom_text_repel(size = 5,
                  max.overlaps = getOption("ggrepel.max.overlaps", default = 200),
                  box.padding = 0.7)



plot(combined_model_plot)

ggsave(paste(Sys.Date(), "Combined_Models_Combined_Time_Corr_Plot_MA_OneColor.svg"), bg = "transparent",
       plot = combined_model_plot, width = 2, height = 2.5, dpi = 300)






##### Correlations Per Model #####

time <- c("Early Resuscitation", "Late Resuscitation")

model_vector <- c("DCBI", "REBOA Zone I + DCBI", "REBOA Zone III + DCBI")

proteins <- colnames(dat_teg_MA[,6:ncol(dat_teg_MA)])

names.vect <- c()

cor.val.vect <- c()

cor_df_list <- list()

cor_df_list_holder <- list()


for(f in 1:length(model_vector)) {
  
  model_select <- filter(dat_teg_MA, Model == model_vector[f])
  
  for(g in 1:length(time)){
    
    q_select <- filter(model_select, Time_Group == time[g])
    
    
    # remove meta data
    df.first <- q_select[ ,-c(1,2,3,5)]
    
    # run spearman correlation
    #data.cor <- as.matrix(cor(df.first, method = "pearson")[,1])
    
    # rename column
    #colnames(data.cor) <- c("SCorrelation")
    
    # initialize empty vector
    p.val.vect <- c()
    
    # loop through dataframe, calculating p value corresponding to the calculation 
    for(i in 1:ncol(df.first)){
      
      
      s <- df.first[,c(1,i)]
      s <- na.omit(s)
      s <- as.matrix(s)
      
      if(nrow(s)>1){
        
        a <- colnames(s)[2]
        names.vect[i] <- a
        
        # run the test
        cor.result <- stats::cor.test(s[,1], s[,2], method = "pearson")
        # extract the estimate
        cor.val <- cor.result$estimate
        # store the estimate
        cor.val.vect[i] <- cor.val
        
        # extract the p value
        p.val <- cor.result$p.value
        
        # append the p value to the vector
        p.val.vect[i] <- p.val
      }
    }
    
    ## Data coercion ##
    
    df.p <- data.frame(Protein = names.vect, P_Value = p.val.vect, SCorrelation = cor.val.vect)
    
    df.final <- df.p %>% 
      na.omit() %>% 
      slice(-1) %>% 
      mutate(NegLogP = -log10(P_Value), 
             Time_Group = time[g],
             Model = model_vector[f],
             Significance = case_when(P_Value < 0.05 & SCorrelation > 0 ~ "positive significant",
                                      P_Value < 0.05 & SCorrelation < 0 ~ "negative significant",
                                      TRUE ~ "insignificant"),
             Point_Color = case_when(P_Value < 0.05 ~ Time_Group,
                                     P_Value > 0.05 ~ ""))
    
    
    cor_df_list[[g]] <- df.final
    
    #assign(paste("corr_df_condition", model[h], sep = "_"), df.final)
    #write.csv(corr1.df, paste("COMBAT_PLT_unit_corr_Group_", h, ".csv"))
  }
  cor_df_list_holder[[f]] <- cor_df_list
}

cor_df <- bind_rows(cor_df_list_holder)

cor_df$Point_Color <- factor(cor_df$Point_Color, levels = c("", "Early Resuscitation", "Late Resuscitation"))






###### Panel E: Zone III U-Plot ######

ZoneIII_df_plot <- cor_df %>% 
  filter(Model == "REBOA Zone III + DCBI") %>% 
  mutate(Label = case_when(NegLogP > 2 & SCorrelation > 0 & Time_Group == "Early Resuscitation" ~ Protein, 
                           NegLogP > 1.5 & SCorrelation > 0 & Time_Group == "Late Resuscitation" ~ Protein,
                           NegLogP > 2 & SCorrelation < 0 & Time_Group == "Early Resuscitation" ~ Protein, 
                           NegLogP > 2 & SCorrelation < 0 & Time_Group == "Late Resuscitation" ~ Protein,
                           TRUE ~ ""))

## Create scatterplot using the ggplot function ##
ZoneIII_corr_plot <- ggplot(ZoneIII_df_plot, aes(SCorrelation, NegLogP, col = Point_Color, label = Label)) + 
  geom_point(size = 0.5) + 
  
  scale_color_manual(values = c("grey", # Insignificant 
                                "#01bdf8", # Early Resuscitation
                                "#000080" # Late Resuscitation
  ))+
  
  scale_x_continuous(breaks = seq(-1, 1, 0.5), limits = c(-1, 1))+
  scale_y_continuous(expand = c(0,0), limits = c(1,5), breaks = seq(1, 5, 0.5))+
  
  xlab("Pearson Correlation") + 
  ylab("-log10(Pvalue)") + 
  
  
  geom_abline(intercept = -log10(0.05), slope = 0, col = 'grey57', lwd = 0.25, lty = 2) +
  
  theme(legend.position = "none", 
        plot.margin = margin(t = 15, r = 10, b = 8, l = 8, unit = "pt"),
        plot.title = element_text(hjust = 0.5, size = 8), 
        axis.title = element_text(size = 8), 
        axis.text.x = element_text(size = 6, vjust = 0.6, angle = 90),
        axis.text.y = element_text(size = 6),
        axis.ticks = element_line(size = 0.25),
        plot.background = element_rect(fill = "transparent",colour = NA),
        panel.background = element_rect(fill = 'white'),
        panel.border = element_rect(colour = "black", fill=NA, size = 0.25),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  
  geom_text_repel(size = 2,
                  max.overlaps = getOption("ggrepel.max.overlaps", default = 100),
                  box.padding = 0.7)


plot(ZoneIII_corr_plot)

ggsave(paste(Sys.Date(), "ZoneIII_Corr_Plot_MA.svg"), bg = "transparent",
       plot = ZoneIII_corr_plot, width = 2, height = 2.5, dpi = 300)










#### LY30 Result Export for GO Analysis ####

##### NO AO (DCBI) ##### 

# DCBI Early Resuscitation
out <- cor_df %>% 
  filter(Model == "DCBI") %>% 
  filter(Time_Group == "Early Resuscitation") %>%
  filter(Significance == "positive significant") %>% 
  select(Protein) 

write.csv(out,
          row.names = FALSE,
          "GO_Analysis_MA_DCBI_Early_Resus_Pos_Sig_Corr.csv")

out <- cor_df %>% 
  filter(Model == "DCBI") %>% 
  filter(Time_Group == "Early Resuscitation") %>%
  filter( Significance == "negative significant") %>% 
  select(Protein)

write.csv(out,
          row.names = FALSE,
          "GO_Analysis_MA_DCBI_Early_Resus_Neg_Sig_Corr.csv")

# DCBI Late Resuscitation
out <- cor_df %>% 
  filter(Model == "DCBI") %>% 
  filter(Time_Group == "Late Resuscitation") %>%
  filter(Significance == "positive significant") %>% 
  select(Protein) 

write.csv(out,
          row.names = FALSE,
          "GO_Analysis_MA_DCBI_Late_Resus_Pos_Sig_Corr.csv")

out <- cor_df %>% 
  filter(Model == "DCBI") %>% 
  filter(Time_Group == "Late Resuscitation") %>%
  filter( Significance == "negative significant") %>% 
  select(Protein)

write.csv(out,
          row.names = FALSE,
          "GO_Analysis_MA_DCBI_Late_Resus_Neg_Sig_Corr.csv")

##### REBOA Zone I ##### 

# REBOA Zone I Early Resuscitation
out <- cor_df %>% 
  filter(Model == "REBOA Zone I + DCBI") %>% 
  filter(Time_Group == "Early Resuscitation") %>%
  filter(Significance == "positive significant") %>% 
  select(Protein) 

write.csv(out,
          row.names = FALSE,
          "GO_Analysis_MA_Zone_I_Early_Resus_Pos_Sig_Corr.csv")

out <- cor_df %>% 
  filter(Model == "REBOA Zone I + DCBI") %>% 
  filter(Time_Group == "Early Resuscitation") %>%
  filter( Significance == "negative significant") %>% 
  select(Protein)

write.csv(out,
          row.names = FALSE,
          "GO_Analysis_MA_Zone_I_Early_Resus_Neg_Sig_Corr.csv")

# REBOA Zone I Late Resuscitation
out <- cor_df %>% 
  filter(Model == "REBOA Zone I + DCBI") %>% 
  filter(Time_Group == "Late Resuscitation") %>%
  filter(Significance == "positive significant") %>% 
  select(Protein) 

write.csv(out,
          row.names = FALSE,
          "GO_Analysis_MA_Zone_I_Late_Resus_Pos_Sig_Corr.csv")

out <- cor_df %>% 
  filter(Model == "REBOA Zone I + DCBI") %>% 
  filter(Time_Group == "Late Resuscitation") %>%
  filter( Significance == "negative significant") %>% 
  select(Protein)

write.csv(out,
          row.names = FALSE,
          "GO_Analysis_MA_Zone_I_Late_Resus_Neg_Sig_Corr.csv")



##### REBOA Zone III ##### 

# REBOA Zone III Early Resuscitation
out <- cor_df %>% 
  filter(Model == "REBOA Zone III + DCBI") %>% 
  filter(Time_Group == "Early Resuscitation") %>%
  filter(Significance == "positive significant") %>% 
  select(Protein) 

write.csv(out,
          row.names = FALSE,
          "GO_Analysis_MA_ZONE_III_Early_Resus_Pos_Sig_Corr.csv")

out <- cor_df %>% 
  filter(Model == "REBOA Zone III + DCBI") %>% 
  filter(Time_Group == "Early Resuscitation") %>%
  filter( Significance == "negative significant") %>% 
  select(Protein)

write.csv(out,
          row.names = FALSE,
          "GO_Analysis_MA_ZONE_III_Early_Resus_Neg_Sig_Corr.csv")

# REBOA Zone III Late Resuscitation
out <- cor_df %>% 
  filter(Model == "REBOA Zone III + DCBI") %>% 
  filter(Time_Group == "Late Resuscitation") %>%
  filter(Significance == "positive significant") %>% 
  select(Protein) 

write.csv(out,
          row.names = FALSE,
          "GO_Analysis_MA_ZONE_III_Late_Resus_Pos_Sig_Corr.csv")

out <- cor_df %>% 
  filter(Model == "REBOA Zone III + DCBI") %>% 
  filter(Time_Group == "Late Resuscitation") %>%
  filter( Significance == "negative significant") %>% 
  select(Protein)

write.csv(out,
          row.names = FALSE,
          "GO_Analysis_MA_ZONE_III_Late_Resus_Neg_Sig_Corr.csv")



##### MA Zone I vs Zone III, Pos Corr #####

dat_1 <- read.csv("MA_ZoneI_vs_ZoneIII_pos_corr.csv") 


dat_1$Model <- factor(dat_1$GO_BP)


## Create Bar plot using the ggplot function ##
pos_go_plot <- ggplot(dat_1, aes(x = Order, y = Log2_FC, fill = Resus_Phase)) + 
  
  geom_bar(stat = "identity", width = 0.5) + 
  coord_flip() +
  
  geom_hline(yintercept = 0, color = "black", size = 0.25)+
  
  scale_x_discrete(labels = dat_1$GO_BP)+
  scale_y_continuous(expand = c(0,0), limits = c(-15,15), breaks = seq(-15,15,5))+
  
  xlab("") + 
  ylab("Log2(Fold Change)") + 
  
  scale_fill_manual(values = c(
    "#01bdf8", # Early Resuscitation
    "#000080" # Late Resuscitation
    
  ))+
  
  theme(legend.position = "none", 
        plot.margin = margin(t = 15, r = 10, b = 8, l = 24, unit = "pt"),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 16), 
        axis.title = element_text(size = 8, color = 'black'), 
        axis.text.x = element_text(size = 6, vjust = 0.6, angle = 90, color = 'black'),
        axis.text.y = element_text(size = 6, color = 'black'),
        #axis.text.y = element_blank(),
        axis.ticks = element_line(size = 0.25),
        #axis.ticks.y = element_blank(),
        plot.background = element_rect(fill = "transparent",colour = NA),
        panel.background = element_rect(fill = NA),
        panel.border = element_rect(colour = "black", fill=NA, size = 0.25),
        axis.line = element_line(color = 'black', size = 0, linetype = 'solid'),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())


plot(pos_go_plot)

ggsave(paste(Sys.Date(), "MA_Pos_Corr_GO_Plot.svg"), bg = "transparent",
       plot = pos_go_plot, width = 4, height = 1.5, dpi = 300)






##### MA Zone I vs Zone III, Neg Corr #####

dat_2 <- read.csv("MA_ZoneI_vs_ZoneIII_neg_corr.csv") 


dat_2$Model <- factor(dat_2$GO_BP)


## Create Bar plot using the ggplot function ##
pos_go_plot <- ggplot(dat_2, aes(x = Order, y = Log2_FC, fill = Resus_Phase)) + 
  
  geom_bar(stat = "identity", width = 0.5) + 
  coord_flip() +
  
  geom_hline(yintercept = 0, color = "black", size = 0.25)+
  
  scale_x_discrete(labels = dat_2$GO_BP)+
  scale_y_continuous(expand = c(0,0), limits = c(-15,15), breaks = seq(-15,15,5))+
  
  xlab("") + 
  ylab("Log2(Fold Change)") + 
  
  scale_fill_manual(values = c(
    "#01bdf8", # Early Resuscitation
    "#000080" # Late Resuscitation
    
  ))+
  
  theme(legend.position = "none", 
        plot.margin = margin(t = 15, r = 10, b = 8, l = 54, unit = "pt"),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 16), 
        axis.title = element_text(size = 8, color = 'black'), 
        axis.text.x = element_text(size = 6, vjust = 0.6, angle = 90, color = 'black'),
        axis.text.y = element_text(size = 6, color = 'black'),
        #axis.text.y = element_blank(),
        axis.ticks = element_line(size = 0.25),
        #axis.ticks.y = element_blank(),
        plot.background = element_rect(fill = "transparent",colour = NA),
        panel.background = element_rect(fill = NA),
        panel.border = element_rect(colour = "black", fill=NA, size = 0.25),
        axis.line = element_line(color = 'black', size = 0, linetype = 'solid'),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())


plot(pos_go_plot)

ggsave(paste(Sys.Date(), "MA_Neg_Corr_GO_Plot.svg"), bg = "transparent",
       plot = pos_go_plot, width = 4, height = 1.5, dpi = 300)





#### ******* END ****** ####

#### ***Unused Code*** ####
#### Waterfall Plot ####

##### Tissue Injury (TI) #####
dat_wat_ti <- dat_merged %>% 
  select(-Date, -Time, -SAS.time) %>% 
  filter(Model == "SHAM" | Model == "TI") %>% 
  filter(Time.from.EOS.min == 30) %>% 
  tidyr::pivot_longer(7:985, 
                      names_to = "Analyte", 
                      values_to = "Raw_Intensity") %>% 
  group_by(Analyte) %>% 
  dplyr::summarise(Ratio = log2(mean(Raw_Intensity[Model == "TI"])/mean(Raw_Intensity[Model == "SHAM"]))) %>% 
  arrange(Ratio) %>% 
  tibble::rowid_to_column("Analyte_Number") %>% 
  mutate(Label = case_when(Ratio > 6 ~ Analyte, 
                           Ratio < -5 ~ Analyte,
                           Analyte == "Succinate" | Analyte == "Lactate" ~ Analyte,
                           TRUE ~ ""))




FC_Plot <- ggplot(dat_wat_ti, aes(x = Analyte_Number, y = Ratio, label = Label))+
  
  #geom_hline(yintercept = 70, linetype = "dashed", color = '#999999')+
  #geom_hline(yintercept = -70, linetype = "dashed", color = '#999999')+
  geom_hline(yintercept = 0)+
  
  geom_point(size=1)+
  #geom_point(aes(size= Category_Numeric>0))+
  #scale_size_manual(values = c(1, 2))+
  
  #scale_color_manual(values=c('#5e3399','#b63d46', '#339991', 'black', '#2055da', '#daa520', '#c0c0c0'))+
  
  scale_y_continuous(expand = c(0,0),limits = c(-10,10), breaks = seq(-10,10, 5)) +
  #scale_x_continuous(expand = expansion(add = 20), limits = c(0,1400), breaks = seq(400,1400, 400))+
  
  theme(plot.margin = margin(0.5,0.5,0.5,0.5, "cm"),
        legend.position = "none", 
        axis.title.x = element_blank(), 
        axis.title.y = element_blank(),
        #axis.text.x = element_blank(),
        #axis.ticks.x = element_blank(),
        axis.text.x = element_text(size = 18, color = 'black', angle = 90, hjust=0.95, vjust=0.5),
        axis.text.y = element_text(size = 18, color = 'black'),
        plot.background = element_rect(fill = "transparent",colour = NA),
        #axis.line = element_line(colour = "black", size = 0.5, linetype = "solid"),
        panel.border = element_rect(colour = "black", size = 0.5, linetype = "solid", fill = NA),
        panel.background = element_rect(fill = 'white'),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  
  geom_text_repel(size = 4,
                  max.overlaps = getOption("ggrepel.max.overlaps", default = 500), 
                  box.padding = 0.7)

plot(FC_Plot)  


ggsave(FC_Plot, 
       file = paste0(Sys.Date(),"Log2FC_change_metabolites_and_proteins_TI.svg"),
       width = 3, height = 5, units = "in", dpi = 500, bg = "transparent")



##### Hemorrhagic Shock (HS) #####
dat_wat_hs <- dat_merged %>% 
  select(-Date, -Time, -SAS.time) %>% 
  filter(Model == "SHAM" | Model == "HS") %>% 
  filter(Time.from.EOS.min == 0) %>% 
  tidyr::pivot_longer(7:985, 
                      names_to = "Analyte", 
                      values_to = "Raw_Intensity") %>% 
  group_by(Analyte) %>% 
  dplyr::summarise(Ratio = log2(mean(Raw_Intensity[Model == "HS"])/mean(Raw_Intensity[Model == "SHAM"]))) %>% 
  arrange(Ratio) %>% 
  tibble::rowid_to_column("Analyte_Number") %>% 
  mutate(Label = case_when(Ratio > 7 ~ Analyte, 
                           Ratio < -5 ~ Analyte,
                           Analyte == "Succinate" | Analyte == "Lactate" ~ Analyte,
                           TRUE ~ ""))




FC_Plot <- ggplot(dat_wat_hs, aes(x = Analyte_Number, y = Ratio, label = Label))+
  
  #geom_hline(yintercept = 70, linetype = "dashed", color = '#999999')+
  #geom_hline(yintercept = -70, linetype = "dashed", color = '#999999')+
  geom_hline(yintercept = 0)+
  
  geom_point(size=1)+
  #geom_point(aes(size= Category_Numeric>0))+
  #scale_size_manual(values = c(1, 2))+
  
  #scale_color_manual(values=c('#5e3399','#b63d46', '#339991', 'black', '#2055da', '#daa520', '#c0c0c0'))+
  
  scale_y_continuous(expand = c(0,0),limits = c(-10,10), breaks = seq(-10,10, 5)) +
  #scale_x_continuous(expand = expansion(add = 20), limits = c(0,1400), breaks = seq(400,1400, 400))+
  
  theme(plot.margin = margin(0.5,0.5,0.5,0.5, "cm"),
        legend.position = "none", 
        axis.title.x = element_blank(), 
        axis.title.y = element_blank(),
        #axis.text.x = element_blank(),
        #axis.ticks.x = element_blank(),
        axis.text.x = element_text(size = 18, color = 'black', angle = 90, hjust=0.95, vjust=0.5),
        axis.text.y = element_text(size = 18, color = 'black'),
        plot.background = element_rect(fill = "transparent",colour = NA),
        #axis.line = element_line(colour = "black", size = 0.5, linetype = "solid"),
        panel.border = element_rect(colour = "black", size = 0.5, linetype = "solid", fill = NA),
        panel.background = element_rect(fill = 'white'),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  
  geom_text_repel(size = 4,
                  max.overlaps = getOption("ggrepel.max.overlaps", default = 200), 
                  box.padding = 0.7)

plot(FC_Plot)  


ggsave(FC_Plot, 
       file = paste0(Sys.Date(),"Log2FC_change_metabolites_and_proteins_HS.svg"),
       width = 3, height = 5, units = "in", dpi = 500, bg = "transparent")



##### Combined (TI/HS) #####
dat_wat_ti_hs <- dat_merged %>% 
  select(-Date, -Time, -SAS.time) %>% 
  filter(Model == "SHAM" | Model == "TI + HS") %>% 
  filter(Time.from.EOS.min == 0) %>% 
  tidyr::pivot_longer(7:985, 
                      names_to = "Analyte", 
                      values_to = "Raw_Intensity") %>% 
  group_by(Analyte) %>% 
  dplyr::summarise(Ratio = log2(mean(Raw_Intensity[Model == "TI + HS"])/mean(Raw_Intensity[Model == "SHAM"]))) %>% 
  arrange(Ratio) %>% 
  tibble::rowid_to_column("Analyte_Number") %>% 
  mutate(Label = case_when(Ratio > 6.5 ~ Analyte, 
                           Ratio < -5 ~ Analyte,
                           Analyte == "Succinate" | Analyte == "Lactate" ~ Analyte,
                           TRUE ~ ""))




FC_Plot <- ggplot(dat_wat_ti_hs, aes(x = Analyte_Number, y = Ratio, label = Label))+
  
  #geom_hline(yintercept = 70, linetype = "dashed", color = '#999999')+
  #geom_hline(yintercept = -70, linetype = "dashed", color = '#999999')+
  geom_hline(yintercept = 0)+
  
  geom_point(size=1)+
  #geom_point(aes(size= Category_Numeric>0))+
  #scale_size_manual(values = c(1, 2))+
  
  #scale_color_manual(values=c('#5e3399','#b63d46', '#339991', 'black', '#2055da', '#daa520', '#c0c0c0'))+
  
  #scale_y_continuous(expand = c(0,0),limits = c(-10,10), breaks = seq(-10,10, 5)) +
  #scale_x_continuous(expand = expansion(add = 20), limits = c(0,1400), breaks = seq(400,1400, 400))+
  
  theme(plot.margin = margin(0.5,0.5,0.5,0.5, "cm"),
        legend.position = "none", 
        axis.title.x = element_blank(), 
        axis.title.y = element_blank(),
        #axis.text.x = element_blank(),
        #axis.ticks.x = element_blank(),
        axis.text.x = element_text(size = 18, color = 'black', angle = 90, hjust=0.95, vjust=0.5),
        axis.text.y = element_text(size = 18, color = 'black'),
        plot.background = element_rect(fill = "transparent",colour = NA),
        #axis.line = element_line(colour = "black", size = 0.5, linetype = "solid"),
        panel.border = element_rect(colour = "black", size = 0.5, linetype = "solid", fill = NA),
        panel.background = element_rect(fill = 'white'),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  
  geom_text_repel(size = 4,
                  max.overlaps = getOption("ggrepel.max.overlaps", default = 200), 
                  box.padding = 0.7)

plot(FC_Plot)  


ggsave(FC_Plot, 
       file = paste0(Sys.Date(),"Log2FC_change_metabolites_and_proteins_TI_HS.svg"),
       width = 3, height = 5, units = "in", dpi = 500, bg = "transparent")












#### Metabolomics LY30 Correlation ####

##### Metabolomics Pearson Correlation LY30 ######

##### Grouped Time Points Corr #####

dat_teg_ly30 <- dat.clin %>% 
  select(Pig.Number, Time..Min..from.End.of.Shock, CN.LY30...) %>% 
  dplyr::rename(CN_LY30 = CN.LY30..., Pig.id = Pig.Number, Time.from.EOS.min = Time..Min..from.End.of.Shock) %>% 
  filter_at(vars(CN_LY30), all_vars(!is.na(.))) %>% 
  left_join(dat.meta, by = c("Pig.id", "Time.from.EOS.min")) %>% 
  select(-c(6:14)) %>%  
  na.omit() %>% 
  left_join(dat_metab, by = "Vial") %>% 
  tidyr::pivot_longer(6:130, 
                      names_to = "Metabolite", 
                      values_to = "Raw_Intensity") %>% 
  mutate(Norm_Intensity = (Raw_Intensity - mean(Raw_Intensity))/sd(Raw_Intensity)) %>% 
  tidyr::pivot_wider(id_cols = 1:5, names_from = Metabolite, values_from = Raw_Intensity) %>% 
  filter(Model %in% c("TI + HS", "HS", "TI", "SHAM")) %>% 
  mutate(Time_Group = case_when(Time.from.EOS.min %in% c(0, 30, 60) ~ "Early Resuscitation",
                                Time.from.EOS.min %in% c(120, 180, 240) ~ "Late Resuscitation", 
                                TRUE ~ "Baseline")) %>% 
  relocate(Time_Group, .after = Time.from.EOS.min)

time <- c("Early Resuscitation", "Late Resuscitation")

model_vector <- c("TI + HS", "HS", "TI", "SHAM")

proteins <- colnames(dat_teg_ly30[,5:ncol(dat_teg_ly30)])

names.vect <- c()

cor.val.vect <- c()

cor_df_list <- list()

cor_df_list_holder <- list()


for(f in 1:length(model_vector)) {
  
  model_select <- filter(dat_teg_ly30, Model == model_vector[f])
  
  for(g in 1:length(time)){
    
    q_select <- filter(model_select, Time_Group == time[g])
    
    
    # remove meta data
    df.first <- q_select[ ,-c(1,2,3,5,6)]
    
    
    # initialize empty vector
    p.val.vect <- c()
    
    # loop through dataframe, calculating p value corresponding to the calculation 
    for(i in 1:ncol(df.first)){
      
      
      s <- df.first[,c(1,i)]
      s <- na.omit(s)
      s <- as.matrix(s)
      
      if(nrow(s)>1){
        
        a <- colnames(s)[2]
        names.vect[i] <- a
        
        # run the test
        cor.result <- stats::cor.test(s[,1], s[,2], method = "spearman")
        # extract the estimate
        cor.val <- cor.result$estimate
        # store the estimate
        cor.val.vect[i] <- cor.val
        
        # extract the p value
        p.val <- cor.result$p.value
        
        # append the p value to the vector
        p.val.vect[i] <- p.val
      }
    }
    
    ## Data coercion ##
    
    df.p <- data.frame(Protein = names.vect, P_Value = p.val.vect, SCorrelation = cor.val.vect)
    
    df.final <- df.p %>% 
      na.omit() %>% 
      slice(-1) %>% 
      mutate(NegLogP = -log10(P_Value), 
             Time_Group = time[g],
             Model = model_vector[f],
             Significance = case_when(P_Value < 0.05 & SCorrelation > 0 ~ "positive significant",
                                      P_Value < 0.05 & SCorrelation < 0 ~ "negative significant",
                                      TRUE ~ "insignificant"),
             Point_Color = case_when(P_Value < 0.05 ~ Time_Group,
                                     P_Value > 0.05 ~ ""))
    
    
    cor_df_list[[g]] <- df.final
    
    #assign(paste("corr_df_condition", model[h], sep = "_"), df.final)
    #write.csv(corr1.df, paste("COMBAT_PLT_unit_corr_Group_", h, ".csv"))
  }
  cor_df_list_holder[[f]] <- cor_df_list
}

cor_df <- bind_rows(cor_df_list_holder)

cor_df$Point_Color <- factor(cor_df$Point_Color, levels = c("", "Early Resuscitation", "Late Resuscitation"))



###### SHAM Plot ######

SHAM_df_plot <- cor_df %>% 
  filter(Model == "SHAM") %>% 
  mutate(Label = case_when(NegLogP > 2.8 & SCorrelation > 0 & Time_Group == "Early Resuscitation" ~ Protein, 
                           NegLogP > 2.4 & SCorrelation > 0 & Time_Group == "Late Resuscitation" ~ Protein,
                           NegLogP > 2.5 & SCorrelation < 0 & Time_Group == "Early Resuscitation" ~ Protein, 
                           NegLogP > 2.5 & SCorrelation < 0 & Time_Group == "Late Resuscitation" ~ Protein,
                           TRUE ~ ""))

## Create scatterplot using the ggplot function ##
SHAM_corr_plot <- ggplot(SHAM_df_plot, aes(SCorrelation, NegLogP, col = Point_Color, label = Label)) + 
  scale_color_manual(values = c("grey", # Insignificant 
                                "#01bdf8", # Early Resuscitation
                                "#000080" # Late Resuscitation
  ))+
  
  scale_x_continuous(breaks = seq(-1,1,0.5), limits = c(-1, 1))+
  scale_y_continuous(expand = c(0,0), limits = c(0,4), breaks = seq(0,4,1))+
  
  xlab("Pearson Correlation") + 
  ylab("-log10(Pvalue)") + 
  
  geom_point() + 
  geom_abline(intercept = -log10(0.05), slope = 0, col = 'grey57', lwd = 0.5, lty = 2) +
  
  theme(legend.position = "none", 
        plot.margin = margin(t = 15, r = 10, b = 8, l = 8, unit = "pt"),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 16), 
        axis.title = element_text(face = "bold", size = 24), 
        axis.text.x = element_text(face = "bold", size = 20, vjust = 0.6, angle = 90),
        axis.text.y = element_text(face = "bold", size = 20),
        plot.background = element_rect(fill = "transparent",colour = NA),
        panel.background = element_rect(fill = 'white'),
        panel.border = element_rect(colour = "black", fill=NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  
  geom_text_repel(size = 5,
                  max.overlaps = getOption("ggrepel.max.overlaps", default = 100),
                  box.padding = 0.7)

plot(SHAM_corr_plot)

ggsave(paste(Sys.Date(), "SHAM_Corr_Plot_MA.svg"), bg = "transparent",
       plot = SHAM_corr_plot, width = 4.2, height = 6, dpi = 300)


###### TI Plot ######

TI_df_plot <- cor_df %>% 
  filter(Model == "TI") %>% 
  mutate(Label = case_when(NegLogP > 3.0 & SCorrelation > 0 & Time_Group == "Early Resuscitation" ~ Protein, 
                           NegLogP > 2.8 & SCorrelation > 0 & Time_Group == "Late Resuscitation" ~ Protein,
                           NegLogP > 3.0 & SCorrelation < 0 & Time_Group == "Early Resuscitation" ~ Protein, 
                           NegLogP > 3.0 & SCorrelation < 0 & Time_Group == "Late Resuscitation" ~ Protein,
                           TRUE ~ ""))

## Create scatterplot using the ggplot function ##
TI_corr_plot <- ggplot(TI_df_plot, aes(SCorrelation, NegLogP, col = Point_Color, label = Label)) + 
  scale_color_manual(values = c("grey", # Insignificant 
                                "#01bdf8", # Early Resuscitation
                                "#000080" # Late Resuscitation
  ))+
  
  scale_x_continuous(breaks = seq(-1,1,0.5), limits = c(-1, 1))+
  scale_y_continuous(expand = c(0,0), limits = c(0,8), breaks = seq(0,8,1))+
  
  xlab("Pearson Correlation") + 
  ylab("-log10(Pvalue)") + 
  
  geom_point(size = 1) + 
  geom_abline(intercept = -log10(0.05), slope = 0, col = 'grey57', lwd = 0.25, lty = 2) +
  
  theme(legend.position = "none", 
        plot.margin = margin(t = 15, r = 10, b = 8, l = 8, unit = "pt"),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 16), 
        axis.title = element_text(size = 8, color = 'black'), 
        axis.text.x = element_text(size = 6, vjust = 0.6, angle = 90, color = 'black'),
        axis.text.y = element_text(size = 6, color = 'black'),
        plot.background = element_rect(fill = "transparent",colour = NA),
        panel.background = element_rect(fill = NA),
        panel.border = element_rect(colour = "black", fill=NA, size = 0.25),
        #axis.line = element_line(color = 'black', size = 0.25, linetype = 'solid'),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  
  geom_text_repel(size = 2,
                  max.overlaps = getOption("ggrepel.max.overlaps", default = 100),
                  box.padding = 0.7,
                  segment.size = 0.15)


plot(TI_corr_plot)

ggsave(paste(Sys.Date(), "TI_Corr_Plot_MA.svg"), bg = "transparent",
       plot = TI_corr_plot, width = 2.5, height = 3.5, dpi = 300)



###### HS Plot ######

HS_df_plot <- cor_df %>% 
  filter(Model == "HS") %>% 
  mutate(Label = case_when(NegLogP > 2.5 & SCorrelation > 0 & Time_Group == "Early Resuscitation" ~ Protein, 
                           NegLogP > 2.0 & SCorrelation > 0 & Time_Group == "Late Resuscitation" ~ Protein,
                           NegLogP > 1.8 & SCorrelation < 0 & Time_Group == "Early Resuscitation" ~ Protein, 
                           NegLogP > 2.4 & SCorrelation < 0 & Time_Group == "Late Resuscitation" ~ Protein,
                           Protein == "Spermidine" ~ Protein,
                           TRUE ~ ""))

## Create scatterplot using the ggplot function ##
HS_corr_plot <- ggplot(HS_df_plot, aes(SCorrelation, NegLogP, col = Point_Color, label = Label)) + 
  scale_color_manual(values = c("grey", # Insignificant 
                                "#01bdf8", # Early Resuscitation
                                "#000080" # Late Resuscitation
  ))+
  
  scale_x_continuous(breaks = seq(-1, 1, 0.5), limits = c(-1, 1))+
  scale_y_continuous(expand = c(0,0), limits = c(0,4), breaks = seq(0,4,1))+
  
  xlab("Pearson Correlation") + 
  ylab("-log10(Pvalue)") + 
  
  geom_point(size = 1) + 
  geom_abline(intercept = -log10(0.05), slope = 0, col = 'grey57', lwd = 0.25, lty = 2) +
  
  theme(legend.position = "none", 
        plot.margin = margin(t = 15, r = 10, b = 8, l = 8, unit = "pt"),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 16), 
        axis.title = element_text(size = 8, color = 'black'), 
        axis.text.x = element_text(size = 6, vjust = 0.6, angle = 90, color = 'black'),
        axis.text.y = element_text(size = 6, color = 'black'),
        plot.background = element_rect(fill = "transparent",colour = NA),
        panel.background = element_rect(fill = NA),
        panel.border = element_rect(colour = "black", fill=NA, size = 0.25),
        #axis.line = element_line(color = 'black', size = 0.25, linetype = 'solid'),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  
  geom_text_repel(size = 2,
                  max.overlaps = getOption("ggrepel.max.overlaps", default = 100),
                  box.padding = 0.7,
                  segment.size = 0.15)


plot(HS_corr_plot)

ggsave(paste(Sys.Date(), "HS_Corr_Plot_MA.svg"), bg = "transparent",
       plot = HS_corr_plot, width = 2.5, height = 3.5, dpi = 300)



###### TI + HS Plot ######

TI_HS_df_plot <- cor_df %>% 
  filter(Model == "TI + HS") %>% 
  mutate(Label = case_when(NegLogP > 3.5 & SCorrelation > 0 & Time_Group == "Early Resuscitation" ~ Protein, 
                           NegLogP > 2.6 & SCorrelation > 0 & Time_Group == "Late Resuscitation" ~ Protein,
                           NegLogP > 2.8 & SCorrelation < 0 & Time_Group == "Early Resuscitation" ~ Protein, 
                           NegLogP > 2.6 & SCorrelation < 0 & Time_Group == "Late Resuscitation" ~ Protein,
                           TRUE ~ ""))

## Create scatter plot using the ggplot function ##
TI_HS_corr_plot <- ggplot(TI_HS_df_plot, aes(SCorrelation, NegLogP, col = Point_Color, label = Label)) + 
  scale_color_manual(values = c("grey", # Insignificant 
                                "#01bdf8", # Early Resuscitation
                                "#000080" # Late Resuscitation
  ))+
  
  scale_x_continuous(breaks = seq(-1,1,0.5), limits = c(-1, 1))+
  scale_y_continuous(expand = c(0,0), limits = c(0,7), breaks = seq(0,7,1))+
  xlab("Pearson Correlation") + 
  ylab("-log10(Pvalue)") + 
  
  geom_point(size = 1) + 
  geom_abline(intercept = -log10(0.05), slope = 0, col = 'grey57', lwd = 0.25, lty = 2) +
  
  theme(legend.position = "none", 
        plot.margin = margin(t = 15, r = 10, b = 8, l = 8, unit = "pt"),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 16), 
        axis.title = element_text(size = 8, color = 'black'), 
        axis.text.x = element_text(size = 6, vjust = 0.6, angle = 90, color = 'black'),
        axis.text.y = element_text(size = 6, color = 'black'),
        plot.background = element_rect(fill = "transparent",colour = NA),
        panel.background = element_rect(fill = NA),
        panel.border = element_rect(colour = "black", fill=NA, size = 0.25),
        #axis.line = element_line(color = 'black', size = 0.25, linetype = 'solid'),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  
  geom_text_repel(size = 2,
                  max.overlaps = getOption("ggrepel.max.overlaps", default = 100),
                  box.padding = 0.7,
                  segment.size = 0.15)


plot(TI_HS_corr_plot)

ggsave(paste(Sys.Date(), "TI + HS_Corr_Plot_MA.svg"), bg = "transparent",
       plot = TI_HS_corr_plot, width = 2.5, height = 3.5, dpi = 300)












#### Figure 9 ####

##### Proteomics Pearson Correlation MA ######

##### Grouped Time Points Corr #####

dat_teg_MA <- dat.clin %>% 
  select( Pig.Number, Time..Min..from.End.of.Shock, CN.MA.mm.) %>% 
  dplyr::rename(CN_MA = CN.MA.mm., Pig.id = Pig.Number, Time.from.EOS.min = Time..Min..from.End.of.Shock) %>% 
  filter_at(vars(CN_MA), all_vars(!is.na(.))) %>% 
  left_join(dat_prot, by = c("Pig.id", "Time.from.EOS.min")) %>% 
  select(-c(4, 6:10)) %>%  
  na.omit() %>% 
  tidyr::pivot_longer(5:858, 
                      names_to = "Protein", 
                      values_to = "Raw_Intensity") %>%
  mutate(Norm_Intensity = (Raw_Intensity - mean(Raw_Intensity))/sd(Raw_Intensity)) %>% 
  tidyr::pivot_wider(id_cols = 1:4, names_from = Protein, values_from = Norm_Intensity) %>% 
  filter(Model %in% c("DCBI", "REBOA Zone I + DCBI", "REBOA Zone III + DCBI")) %>% 
  mutate(Time_Group = case_when(Time.from.EOS.min %in% c(0, 30, 60) ~ "Early Resuscitation",
                                Time.from.EOS.min %in% c(120, 180, 240) ~ "Late Resuscitation", 
                                TRUE ~ "Baseline")) %>% 
  relocate(Time_Group, .after = Time.from.EOS.min)


##### MA Temporal Trend #####

#### Clin. Var. Line Plots ####
# Function for data summary table used for all longitudinal line plots
data_summary <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(median = median(x[[col]], na.rm=TRUE),
      sd = std.error(x[[col]], na.rm=TRUE))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  data_sum <- rename(data_sum, c("median" = varname))
  return(data_sum)
}



dat.line.plot <- dat_teg_MA

dat.line.plot$Model <- factor(dat.line.plot$Model, 
                              levels = c("TI + HS", "HS", "TI", "SHAM"))



# a pipe that does everything you need
AOI <- dat.line.plot %>%
  select(c(CN_MA, 'Time.from.EOS.min', 'Model', 'Pig.id')) %>% 
  rename_(QUANT = names(.)[1]) %>% 
  group_by(Pig.id) %>% 
  data_summary(varname = "QUANT",
               groupnames = c("Model", "Time.from.EOS.min")) %>%  
  filter(!is.na(Model)== TRUE)

p_1 <- ggplot(AOI, aes(x = Time.from.EOS.min, y = QUANT, group = Model, color = Model)) +
  #geom_errorbar(aes(ymin = QUANT - sd, ymax = QUANT + sd), width = 0.1) +
  
  geom_line(size = 0.5) + geom_point(size = 1) +
  scale_color_manual(values = c("#cc0000", "#5a00b3", "#0066cd","#00994d"))+
  labs(fill = "Method") + 
  
  scale_y_continuous(expand = c(0,0), breaks = seq(60,80, by=5), limits = c(60,80))+
  
  scale_x_continuous(breaks = c(-120, 0, 30, 60, 120, 180, 240), 
                     labels = c("-120", "0", "30","60","120","180","240")) +
  xlab("") + ylab("") + ggtitle("CN-MA") + #"abs(Median PG.Quantity)" 
  
  theme(legend.position = "none", 
        plot.margin = margin(t = 15, r = 10, b = 8, l = 8, unit = "pt"),
        plot.title = element_text(hjust = 0.5, size = 10), 
        axis.title = element_text(size = 8, color = 'black'), 
        axis.text.x = element_text(size = 6, vjust = 0.6, angle = 90, color = 'black'),
        axis.text.y = element_text(size = 6, color = 'black'),
        plot.background = element_rect(fill = "transparent",colour = NA),
        panel.background = element_rect(fill = NA),
        panel.border = element_rect(colour = "black", fill=NA, size = 0.25),
        #axis.line = element_line(color = 'black', size = 0.25, linetype = 'solid'),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())


plot(p_1)

# Save the plot to working directory
ggsave(p_1,
       file = paste(Sys.Date(),"DoD_Trauma_Pig_Line_Plot_without_std_error_CN_MA.svg", sep = "_"),
       width = 3.5, height = 2.5, units = "in", dpi = 300, bg = "transparent")


##### Correlations #####



time <- c("Early Resuscitation", "Late Resuscitation")

model_vector <- c("TI + HS", "HS", "TI", "SHAM")

proteins <- colnames(dat_teg_MA[,6:ncol(dat_teg_MA)])

names.vect <- c()

cor.val.vect <- c()

cor_df_list <- list()

cor_df_list_holder <- list()


for(f in 1:length(model_vector)) {
  
  model_select <- filter(dat_teg_MA, Model == model_vector[f])
  
  for(g in 1:length(time)){
    
    q_select <- filter(model_select, Time_Group == time[g])
    
    
    # remove meta data
    df.first <- q_select[ ,-c(1,2,3,5)]
    
    # run spearman correlation
    #data.cor <- as.matrix(cor(df.first, method = "pearson")[,1])
    
    # rename column
    #colnames(data.cor) <- c("SCorrelation")
    
    # initialize empty vector
    p.val.vect <- c()
    
    # loop through dataframe, calculating p value corresponding to the calculation 
    for(i in 1:ncol(df.first)){
      
      
      s <- df.first[,c(1,i)]
      s <- na.omit(s)
      s <- as.matrix(s)
      
      if(nrow(s)>1){
        
        a <- colnames(s)[2]
        names.vect[i] <- a
        
        # run the test
        cor.result <- stats::cor.test(s[,1], s[,2], method = "pearson")
        # extract the estimate
        cor.val <- cor.result$estimate
        # store the estimate
        cor.val.vect[i] <- cor.val
        
        # extract the p value
        p.val <- cor.result$p.value
        
        # append the p value to the vector
        p.val.vect[i] <- p.val
      }
    }
    
    ## Data coercion ##
    
    df.p <- data.frame(Protein = names.vect, P_Value = p.val.vect, SCorrelation = cor.val.vect)
    
    df.final <- df.p %>% 
      na.omit() %>% 
      slice(-1) %>% 
      mutate(NegLogP = -log10(P_Value), 
             Time_Group = time[g],
             Model = model_vector[f],
             Significance = case_when(P_Value < 0.05 & SCorrelation > 0 ~ "positive significant",
                                      P_Value < 0.05 & SCorrelation < 0 ~ "negative significant",
                                      TRUE ~ "insignificant"),
             Point_Color = case_when(P_Value < 0.05 ~ Time_Group,
                                     P_Value > 0.05 ~ ""))
    
    
    cor_df_list[[g]] <- df.final
    
    #assign(paste("corr_df_condition", model[h], sep = "_"), df.final)
    #write.csv(corr1.df, paste("COMBAT_PLT_unit_corr_Group_", h, ".csv"))
  }
  cor_df_list_holder[[f]] <- cor_df_list
}

cor_df <- bind_rows(cor_df_list_holder)

cor_df$Point_Color <- factor(cor_df$Point_Color, levels = c("", "Early Resuscitation", "Late Resuscitation"))



###### SHAM Plot ######

SHAM_df_plot <- cor_df %>% 
  filter(Model == "SHAM") %>% 
  mutate(Label = case_when(NegLogP > 2.8 & SCorrelation > 0 & Time_Group == "Early Resuscitation" ~ Protein, 
                           NegLogP > 2.4 & SCorrelation > 0 & Time_Group == "Late Resuscitation" ~ Protein,
                           NegLogP > 2.5 & SCorrelation < 0 & Time_Group == "Early Resuscitation" ~ Protein, 
                           NegLogP > 2.5 & SCorrelation < 0 & Time_Group == "Late Resuscitation" ~ Protein,
                           TRUE ~ ""))

## Create scatterplot using the ggplot function ##
SHAM_corr_plot <- ggplot(SHAM_df_plot, aes(SCorrelation, NegLogP, col = Point_Color, label = Label)) + 
  scale_color_manual(values = c("grey", # Insignificant 
                                "#01bdf8", # Early Resuscitation
                                "#000080" # Late Resuscitation
  ))+
  
  scale_x_continuous(breaks = seq(-1,1,0.5), limits = c(-1, 1))+
  scale_y_continuous(expand = c(0,0), limits = c(0,4), breaks = seq(0,4,1))+
  
  xlab("Pearson Correlation") + 
  ylab("-log10(Pvalue)") + 
  
  geom_point() + 
  geom_abline(intercept = -log10(0.05), slope = 0, col = 'grey57', lwd = 0.5, lty = 2) +
  
  theme(legend.position = "none", 
        plot.margin = margin(t = 15, r = 10, b = 8, l = 8, unit = "pt"),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 16), 
        axis.title = element_text(face = "bold", size = 24), 
        axis.text.x = element_text(face = "bold", size = 20, vjust = 0.6, angle = 90),
        axis.text.y = element_text(face = "bold", size = 20),
        plot.background = element_rect(fill = "transparent",colour = NA),
        panel.background = element_rect(fill = 'white'),
        panel.border = element_rect(colour = "black", fill=NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  
  geom_text_repel(size = 5,
                  max.overlaps = getOption("ggrepel.max.overlaps", default = 100),
                  box.padding = 0.7)

plot(SHAM_corr_plot)

ggsave(paste(Sys.Date(), "SHAM_Corr_Plot_MA.svg"), bg = "transparent",
       plot = SHAM_corr_plot, width = 4.2, height = 6, dpi = 300)


###### TI Plot ######

TI_df_plot <- cor_df %>% 
  filter(Model == "TI") %>% 
  mutate(Label = case_when(NegLogP > 3.0 & SCorrelation > 0 & Time_Group == "Early Resuscitation" ~ Protein, 
                           NegLogP > 2.8 & SCorrelation > 0 & Time_Group == "Late Resuscitation" ~ Protein,
                           NegLogP > 3.0 & SCorrelation < 0 & Time_Group == "Early Resuscitation" ~ Protein, 
                           NegLogP > 3.0 & SCorrelation < 0 & Time_Group == "Late Resuscitation" ~ Protein,
                           TRUE ~ ""))

## Create scatterplot using the ggplot function ##
TI_corr_plot <- ggplot(TI_df_plot, aes(SCorrelation, NegLogP, col = Point_Color, label = Label)) + 
  scale_color_manual(values = c("grey", # Insignificant 
                                "#01bdf8", # Early Resuscitation
                                "#000080" # Late Resuscitation
  ))+
  
  scale_x_continuous(breaks = seq(-1,1,0.5), limits = c(-1, 1))+
  scale_y_continuous(expand = c(0,0), limits = c(0,8), breaks = seq(0,8,1))+
  
  xlab("Pearson Correlation") + 
  ylab("-log10(Pvalue)") + 
  
  geom_point(size = 1) + 
  geom_abline(intercept = -log10(0.05), slope = 0, col = 'grey57', lwd = 0.25, lty = 2) +
  
  theme(legend.position = "none", 
        plot.margin = margin(t = 15, r = 10, b = 8, l = 8, unit = "pt"),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 16), 
        axis.title = element_text(size = 8, color = 'black'), 
        axis.text.x = element_text(size = 6, vjust = 0.6, angle = 90, color = 'black'),
        axis.text.y = element_text(size = 6, color = 'black'),
        plot.background = element_rect(fill = "transparent",colour = NA),
        panel.background = element_rect(fill = NA),
        panel.border = element_rect(colour = "black", fill=NA, size = 0.25),
        #axis.line = element_line(color = 'black', size = 0.25, linetype = 'solid'),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  
  geom_text_repel(size = 2,
                  max.overlaps = getOption("ggrepel.max.overlaps", default = 100),
                  box.padding = 0.7,
                  segment.size = 0.15)


plot(TI_corr_plot)

ggsave(paste(Sys.Date(), "TI_Corr_Plot_MA.svg"), bg = "transparent",
       plot = TI_corr_plot, width = 2.5, height = 3.5, dpi = 300)



###### HS Plot ######

HS_df_plot <- cor_df %>% 
  filter(Model == "HS") %>% 
  mutate(Label = case_when(NegLogP > 2.5 & SCorrelation > 0 & Time_Group == "Early Resuscitation" ~ Protein, 
                           NegLogP > 2.0 & SCorrelation > 0 & Time_Group == "Late Resuscitation" ~ Protein,
                           NegLogP > 1.8 & SCorrelation < 0 & Time_Group == "Early Resuscitation" ~ Protein, 
                           NegLogP > 2.4 & SCorrelation < 0 & Time_Group == "Late Resuscitation" ~ Protein,
                           TRUE ~ ""))

## Create scatterplot using the ggplot function ##
HS_corr_plot <- ggplot(HS_df_plot, aes(SCorrelation, NegLogP, col = Point_Color, label = Label)) + 
  scale_color_manual(values = c("grey", # Insignificant 
                                "#01bdf8", # Early Resuscitation
                                "#000080" # Late Resuscitation
  ))+
  
  scale_x_continuous(breaks = seq(-1, 1, 0.5), limits = c(-1, 1))+
  scale_y_continuous(expand = c(0,0), limits = c(0,4), breaks = seq(0,4,1))+
  
  xlab("Pearson Correlation") + 
  ylab("-log10(Pvalue)") + 
  
  geom_point(size = 1) + 
  geom_abline(intercept = -log10(0.05), slope = 0, col = 'grey57', lwd = 0.25, lty = 2) +
  
  theme(legend.position = "none", 
        plot.margin = margin(t = 15, r = 10, b = 8, l = 8, unit = "pt"),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 16), 
        axis.title = element_text(size = 8, color = 'black'), 
        axis.text.x = element_text(size = 6, vjust = 0.6, angle = 90, color = 'black'),
        axis.text.y = element_text(size = 6, color = 'black'),
        plot.background = element_rect(fill = "transparent",colour = NA),
        panel.background = element_rect(fill = NA),
        panel.border = element_rect(colour = "black", fill=NA, size = 0.25),
        #axis.line = element_line(color = 'black', size = 0.25, linetype = 'solid'),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  
  geom_text_repel(size = 2,
                  max.overlaps = getOption("ggrepel.max.overlaps", default = 100),
                  box.padding = 0.7,
                  segment.size = 0.15)


plot(HS_corr_plot)

ggsave(paste(Sys.Date(), "HS_Corr_Plot_MA.svg"), bg = "transparent",
       plot = HS_corr_plot, width = 2.5, height = 3.5, dpi = 300)



###### TI + HS Plot ######

TI_HS_df_plot <- cor_df %>% 
  filter(Model == "TI + HS") %>% 
  mutate(Label = case_when(NegLogP > 3.5 & SCorrelation > 0 & Time_Group == "Early Resuscitation" ~ Protein, 
                           NegLogP > 2.6 & SCorrelation > 0 & Time_Group == "Late Resuscitation" ~ Protein,
                           NegLogP > 2.8 & SCorrelation < 0 & Time_Group == "Early Resuscitation" ~ Protein, 
                           NegLogP > 2.6 & SCorrelation < 0 & Time_Group == "Late Resuscitation" ~ Protein,
                           TRUE ~ ""))

## Create scatter plot using the ggplot function ##
TI_HS_corr_plot <- ggplot(TI_HS_df_plot, aes(SCorrelation, NegLogP, col = Point_Color, label = Label)) + 
  scale_color_manual(values = c("grey", # Insignificant 
                                "#01bdf8", # Early Resuscitation
                                "#000080" # Late Resuscitation
  ))+
  
  scale_x_continuous(breaks = seq(-1,1,0.5), limits = c(-1, 1))+
  scale_y_continuous(expand = c(0,0), limits = c(0,7), breaks = seq(0,7,1))+
  xlab("Pearson Correlation") + 
  ylab("-log10(Pvalue)") + 
  
  geom_point(size = 1) + 
  geom_abline(intercept = -log10(0.05), slope = 0, col = 'grey57', lwd = 0.25, lty = 2) +
  
  theme(legend.position = "none", 
        plot.margin = margin(t = 15, r = 10, b = 8, l = 8, unit = "pt"),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 16), 
        axis.title = element_text(size = 8, color = 'black'), 
        axis.text.x = element_text(size = 6, vjust = 0.6, angle = 90, color = 'black'),
        axis.text.y = element_text(size = 6, color = 'black'),
        plot.background = element_rect(fill = "transparent",colour = NA),
        panel.background = element_rect(fill = NA),
        panel.border = element_rect(colour = "black", fill=NA, size = 0.25),
        #axis.line = element_line(color = 'black', size = 0.25, linetype = 'solid'),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  
  geom_text_repel(size = 2,
                  max.overlaps = getOption("ggrepel.max.overlaps", default = 100),
                  box.padding = 0.7,
                  segment.size = 0.15)


plot(TI_HS_corr_plot)

ggsave(paste(Sys.Date(), "TI + HS_Corr_Plot_MA.svg"), bg = "transparent",
       plot = TI_HS_corr_plot, width = 2.5, height = 3.5, dpi = 300)




#### Export for GO Analysis ####

##### Proteomics MA ##### 

# TI Early Resuscitation
out <- cor_df %>% 
  filter(Model == "TI") %>% 
  filter(Time_Group == "Early Resuscitation") %>%
  filter(Significance == "positive significant") %>% 
  select(Protein) 

write.csv(out,
          row.names = FALSE,
          "Protein_Names_for_GO_Analysis_TI_Early_Resus_Pos_Sig_Corr_to_MA.csv")

out <- cor_df %>% 
  filter(Model == "TI") %>% 
  filter(Time_Group == "Early Resuscitation") %>%
  filter( Significance == "negative significant") %>% 
  select(Protein)

write.csv(out,
          row.names = FALSE,
          "Protein_Names_for_GO_Analysis_TI_Early_Resus_Neg_Sig_Corr_to_MA.csv")

# TI Late Resuscitation
out <- cor_df %>% 
  filter(Model == "TI") %>% 
  filter(Time_Group == "Late Resuscitation") %>%
  filter(Significance == "positive significant") %>% 
  select(Protein) 

write.csv(out,
          row.names = FALSE,
          "Protein_Names_for_GO_Analysis_TI_Late_Resus_Pos_Sig_Corr_to_MA.csv")

out <- cor_df %>% 
  filter(Model == "TI") %>% 
  filter(Time_Group == "Late Resuscitation") %>%
  filter( Significance == "negative significant") %>% 
  select(Protein)

write.csv(out,
          row.names = FALSE,
          "Protein_Names_for_GO_Analysis_TI_Late_Resus_Neg_Sig_Corr_to_MA.csv")

# HS Early Resuscitation
out <- cor_df %>% 
  filter(Model == "HS") %>% 
  filter(Time_Group == "Early Resuscitation") %>%
  filter(Significance == "positive significant") %>% 
  select(Protein) 

write.csv(out,
          row.names = FALSE,
          "Protein_Names_for_GO_Analysis_HS_Early_Resus_Pos_Sig_Corr_to_MA.csv")

out <- cor_df %>% 
  filter(Model == "HS") %>% 
  filter(Time_Group == "Early Resuscitation") %>%
  filter( Significance == "negative significant") %>% 
  select(Protein)

write.csv(out,
          row.names = FALSE,
          "Protein_Names_for_GO_Analysis_HS_Early_Resus_Neg_Sig_Corr_to_MA.csv")

# HS Late Resuscitation
out <- cor_df %>% 
  filter(Model == "HS") %>% 
  filter(Time_Group == "Late Resuscitation") %>%
  filter(Significance == "positive significant") %>% 
  select(Protein) 

write.csv(out,
          row.names = FALSE,
          "Protein_Names_for_GO_Analysis_HS_Late_Resus_Pos_Sig_Corr_to_MA.csv")

out <- cor_df %>% 
  filter(Model == "HS") %>% 
  filter(Time_Group == "Late Resuscitation") %>%
  filter( Significance == "negative significant") %>% 
  select(Protein)

write.csv(out,
          row.names = FALSE,
          "Protein_Names_for_GO_Analysis_HS_Late_Resus_Neg_Sig_Corr_to_MA.csv")

# TI + HS Early Resuscitation
out <- cor_df %>% 
  filter(Model == "TI + HS") %>% 
  filter(Time_Group == "Early Resuscitation") %>%
  filter(Significance == "positive significant") %>% 
  select(Protein) 

write.csv(out,
          row.names = FALSE,
          "Protein_Names_for_GO_Analysis_TI + HS_Early_Resus_Pos_Sig_Corr_to_MA.csv")

out <- cor_df %>% 
  filter(Model == "TI + HS") %>% 
  filter(Time_Group == "Early Resuscitation") %>%
  filter( Significance == "negative significant") %>% 
  select(Protein)

write.csv(out,
          row.names = FALSE,
          "Protein_Names_for_GO_Analysis_TI + HS_Early_Resus_Neg_Sig_Corr_to_MA.csv")

# TI + HS Late Resuscitation
out <- cor_df %>% 
  filter(Model == "TI + HS") %>% 
  filter(Time_Group == "Late Resuscitation") %>%
  filter(Significance == "positive significant") %>% 
  select(Protein) 

write.csv(out,
          row.names = FALSE,
          "Protein_Names_for_GO_Analysis_TI + HS_Late_Resus_Pos_Sig_Corr_to_MA.csv")

out <- cor_df %>% 
  filter(Model == "TI + HS") %>% 
  filter(Time_Group == "Late Resuscitation") %>%
  filter( Significance == "negative significant") %>% 
  select(Protein)

write.csv(out,
          row.names = FALSE,
          "Protein_Names_for_GO_Analysis_TI + HS_Late_Resus_Neg_Sig_Corr_to_MA.csv")


#### GO Analysis Scatter Plots MA #####

##### Read In Data Frames ######
# created through gene ontology portal
# used the sus scrofa data base


# model HS + TI Early Resuscitation
dat_go_01 <- read.csv("Proteomics_HS+TI_Early_Resus_Pos_Sig_Corr_MA.csv") %>% 
  mutate(NegLogP = -log10(P_value),
         Condition = "HS + TI", 
         Time_Group = "Early_Resuscitation")

dat_go_02 <- read.csv("Proteomics_HS+TI_Early_Resus_Neg_Sig_Corr_MA.csv") %>% 
  mutate(Fold_Enrichment = -1*(Fold_Enrichment),
         NegLogP = -log10(P_value),
         Condition = "HS + TI", 
         Time_Group = "Early_Resuscitation") %>% 
  rbind(dat_go_01)


# model HS + TI Late Resuscitation
dat_go_03 <- read.csv("Proteomics_HS+TI_Late_Resus_Pos_Sig_Corr_MA.csv") %>% 
  mutate(NegLogP = -log10(P_value),
         Condition = "HS + TI", 
         Time_Group = "Late_Resuscitation")

dat_go_04 <- read.csv("Proteomics_HS+TI_Late_Resus_Neg_Sig_Corr_MA.csv") %>% 
  mutate(Fold_Enrichment = -1*(Fold_Enrichment),
         NegLogP = -log10(P_value),
         Condition = "HS + TI", 
         Time_Group = "Late_Resuscitation") %>% 
  rbind(dat_go_03, dat_go_02) %>% 
  tidyr::separate(col = GO_biological_process_complete, into = c('GO_biological_process_complete', 'GO_Accession'), sep = 'GO:' ) %>% 
  mutate(GO_biological_process_complete = str_replace(GO_biological_process_complete, "[(]", ""),
         GO_Accession = str_replace(GO_Accession, "[)]", ""))




# model HS Early Resuscitation
dat_go_05 <- read.csv("Proteomics_HS_Early_Resus_Pos_Sig_Corr_MA.csv") %>% 
  mutate(NegLogP = -log10(P_value),
         Condition = "HS", 
         Time_Group = "Early_Resuscitation")

dat_go_06 <- read.csv("Proteomics_HS_Early_Resus_Neg_Sig_Corr_MA.csv") %>% 
  mutate(Fold_Enrichment = -1*(Fold_Enrichment),
         NegLogP = -log10(P_value),
         Condition = "HS", 
         Time_Group = "Early_Resuscitation") %>% 
  rbind(dat_go_05)


# model HS Late Resuscitation
dat_go_07 <- read.csv("Proteomics_HS_Late_Resus_Pos_Sig_Corr_MA.csv") %>% 
  mutate(NegLogP = -log10(P_value),
         Condition = "HS", 
         Time_Group = "Late_Resuscitation")

dat_go_08 <- read.csv("Proteomics_HS_Late_Resus_Neg_Sig_Corr_MA.csv") %>% 
  mutate(Fold_Enrichment = -1*(Fold_Enrichment),
         NegLogP = -log10(P_value),
         Condition = "HS", 
         Time_Group = "Late_Resuscitation") %>% 
  rbind(dat_go_07, dat_go_06) %>% 
  tidyr::separate(col = GO_biological_process_complete, into = c('GO_biological_process_complete', 'GO_Accession'), sep = 'GO:' ) %>% 
  mutate(GO_biological_process_complete = str_replace(GO_biological_process_complete, "[(]", ""),
         GO_Accession = str_replace(GO_Accession, "[)]", ""))





# model TI Early Resuscitation
dat_go_09 <- read.csv("Proteomics_TI_Early_Resus_Pos_Sig_Corr_MA.csv") %>% 
  mutate(NegLogP = -log10(P_value),
         Condition = "TI", 
         Time_Group = "Early_Resuscitation")

dat_go_10 <- read.csv("Proteomics_TI_Early_Resus_Neg_Sig_Corr_MA.csv") %>% 
  mutate(Fold_Enrichment = -1*(Fold_Enrichment),
         NegLogP = -log10(P_value),
         Condition = "TI", 
         Time_Group = "Early_Resuscitation") %>% 
  rbind(dat_go_09)


# model TI Late Resuscitation
dat_go_11 <- read.csv("Proteomics_TI_Late_Resus_Pos_Sig_Corr_MA.csv") %>% 
  mutate(NegLogP = -log10(P_value),
         Condition = "TI", 
         Time_Group = "Late_Resuscitation")

dat_go_12 <- read.csv("Proteomics_TI_Late_Resus_Neg_Sig_Corr_MA.csv") %>% 
  mutate(Fold_Enrichment = -1*(Fold_Enrichment),
         NegLogP = -log10(P_value),
         Condition = "TI", 
         Time_Group = "Late_Resuscitation") %>% 
  rbind(dat_go_11, dat_go_10) %>% 
  tidyr::separate(col = GO_biological_process_complete, into = c('GO_biological_process_complete', 'GO_Accession'), sep = 'GO:' ) %>% 
  mutate(GO_biological_process_complete = str_replace(GO_biological_process_complete, "[(]", ""),
         GO_Accession = str_replace(GO_Accession, "[)]", ""))

dat_go_12$Fold_Enrichment <- as.numeric(dat_go_12$Fold_Enrichment)

##### TI + HS Scatter Plot #####


TI_HS_df_plot <- dat_go_04 %>% 
  filter(Time_Group == "Early_Resuscitation")

## Create scatter plot using the ggplot function ##
TI_HS_go_plot <- ggplot(TI_HS_df_plot, aes(Fold_Enrichment, NegLogP, col = Time_Group)) + 
  scale_color_manual(values = c("#01bdf8", # Early Resuscitation
                                "#000080" # Late Resuscitation
  ))+
  
  #scale_x_continuous(breaks = seq(-1,1,0.5), limits = c(-1, 1))+
  scale_y_continuous(expand = c(0,0), limits = c(-log10(0.05),18), breaks = seq(2,18,2))+
  
  xlab("Fold Enrichment") + 
  ylab("-log10(P-value)") + 
  
  geom_point(size = 1) + 
  
  
  theme(legend.position = "none", 
        plot.margin = margin(t = 15, r = 10, b = 8, l = 8, unit = "pt"),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 16), 
        axis.title = element_text(size = 8, color = 'black'), 
        axis.text.x = element_text(size = 6, vjust = 0.6, angle = 90, color = 'black'),
        axis.text.y = element_text(size = 6, color = 'black'),
        plot.background = element_rect(fill = "transparent",colour = NA),
        panel.background = element_rect(fill = NA),
        #panel.border = element_rect(colour = "black", fill=NA),
        axis.line = element_line(color = 'black', size = 0.25, linetype = 'solid'),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  
  geom_text_repel(label = TI_HS_df_plot$GO_biological_process_complete, size = 4,
                  max.overlaps = getOption("ggrepel.max.overlaps", default = 100),
                  box.padding = 0.7,
                  segment.size = 0.15)


plot(TI_HS_go_plot)

ggsave(paste(Sys.Date(), "TI + HS_GO_Plot_MA.svg"), bg = "transparent",
       plot = TI_HS_go_plot, width = 2.5, height = 3.5, dpi = 300)



##### HS Scatter Plot #####

HS_df_plot <- dat_go_08 %>% 
  filter(Time_Group == "Early_Resuscitation")

## Create scatter plot using the ggplot function ##
HS_go_plot <- ggplot(HS_df_plot, aes(Fold_Enrichment, NegLogP, col = Time_Group)) + 
  scale_color_manual(values = c("#01bdf8", # Early Resuscitation
                                "#000080" # Late Resuscitation
  ))+
  
  #scale_x_continuous(breaks = seq(-1,1,0.5), limits = c(-1, 1))+
  scale_y_continuous(expand = c(0,0), limits = c(-log10(0.05),12), breaks = seq(2,12,2))+
  
  xlab("Fold Enrichment") + 
  ylab("-log10(P-value)") + 
  
  geom_point(size = 1) + 
  
  
  theme(legend.position = "none", 
        plot.margin = margin(t = 15, r = 10, b = 8, l = 8, unit = "pt"),
        plot.title = element_text(hjust = 0.5, size = 16), 
        axis.title = element_text(size = 8, color = 'black'), 
        axis.text.x = element_text(size = 6, vjust = 0.6, angle = 90, color = 'black'),
        axis.text.y = element_text(size = 6, color = 'black'),
        plot.background = element_rect(fill = "transparent",colour = NA),
        panel.background = element_rect(fill = NA),
        #panel.border = element_rect(colour = "black", fill=NA),
        axis.line = element_line(color = 'black', size = 0.25, linetype = 'solid'),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  
  geom_text_repel(label = HS_df_plot$GO_biological_process_complete, size = 4,
                  max.overlaps = getOption("ggrepel.max.overlaps", default = 50),
                  box.padding = 0.7,
                  segment.size = 0.15)


plot(HS_go_plot)

ggsave(paste(Sys.Date(), "HS_GO_Plot_NO_LABELS.svg"), bg = "transparent",
       plot = HS_go_plot, width = 2.5, height = 3.5, dpi = 300)





##### TI Scatter Plot #####

TI_df_plot <- dat_go_12 %>% 
  filter(Time_Group == "Early_Resuscitation")

## Create scatter plot using the ggplot function ##
TI_go_plot <- ggplot(TI_df_plot, aes(Fold_Enrichment, NegLogP, col = Time_Group)) + 
  scale_color_manual(values = c("#01bdf8", # Early Resuscitation
                                "#000080" # Late Resuscitation
  ))+
  
  scale_x_continuous(breaks = seq(-100,100,50), limits = c(-100, 100))+
  scale_y_continuous(expand = c(0,0), limits = c(-log10(0.05),15), breaks = seq(2.5,15,2.5))+
  
  xlab("Fold Enrichment") + 
  ylab("-log10(P-value)") + 
  
  geom_point(size = 1) + 
  
  
  theme(legend.position = "none", 
        plot.margin = margin(t = 15, r = 10, b = 8, l = 8, unit = "pt"),
        plot.title = element_text(hjust = 0.5, size = 16), 
        axis.title = element_text(size = 8, color = 'black'), 
        axis.text.x = element_text(size = 6, vjust = 0.6, angle = 90, color = 'black'),
        axis.text.y = element_text(size = 6, color = 'black'),
        plot.background = element_rect(fill = "transparent",colour = NA),
        panel.background = element_rect(fill = NA),
        #panel.border = element_rect(colour = "black", fill=NA),
        axis.line = element_line(color = 'black', size = 0.25, linetype = 'solid'),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  
  geom_text_repel(label = TI_df_plot$GO_biological_process_complete, size = 4,
                  max.overlaps = getOption("ggrepel.max.overlaps", default = 80),
                  box.padding = 0.7,
                  segment.size = 0.15)


plot(TI_go_plot)

ggsave(paste(Sys.Date(), "TI_GO_Plot_NO_LABELS.svg"), bg = "transparent",
       plot = TI_go_plot, width = 2.5, height = 3.5, dpi = 300)




