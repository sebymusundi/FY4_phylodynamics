


# check if BiocManager is installed 
if(!require("BiocManager", quietly = T)) {install.packages("Biocmanager")}
if(!require("remotes", quietly = T)) {install.packages("remote")}


# List packages 


all_packages <- c("tidyverse", "cowplot", "ape", 
                  "phytools", "readxl", "reticulate", 
                  "colorspace", "diagram", "conflicted",
                  "geodata", "magick", "patchwork", "fastmap")


# Loop through each package
for (pkg in all_packages) {
  # Install only if the package is not installed
  if (!pkg %in% installed.packages()) {
    install.packages(pkg)
  }
}


# install packages from github
BiocManager::install("YuLab-SMU/treedataverse", ask = FALSE, force = T)
remotes::install_github("sdellicour/seraphim/unix_OS", force = T, ask=FALSE)


require(treedataverse)
require(tidyverse)
require(cowplot)
require(ape)
require(phytools)
require(readxl)
require(reticulate)
require(colorspace)
require(seraphim)
require(diagram)
require(conflicted)
require(geodata)
require(magick)
require(patchwork)


# Remove scientific notation
options(scipen = 1000000)
# Define colors to apply in subsequent images 

cols <- c( "gold4", "aquamarine4", "steelblue3", "black", "plum4", "salmon4",  
           "indianred", "sandybrown", "deeppink4", "olivedrab", "lightskyblue4", 
           "navajowhite2", "palevioletred", "cadetblue3")

###############################################################################
# Positivity plot for SARS-CoV-2 
################################################################################

# Load data 
positivity_data <- read_excel("~/Library/CloudStorage/OneDrive-KemriWellcomeTrust/draft_manuscript/FY4/Dataverse/FY4Rproject/Data/HF_weekly_SARs.xlsx")

# number of samples collected between January to July 
sum(positivity_data$tested)

# Number of positive cases
sum(positivity_data$SARS_CoV_2_bin)

# plot of SARS-CoV-2 positive samples
supplementary_plot_1 <- ggplot(positivity_data, aes(week, tested)) +
  geom_line()+
  geom_point() +
  geom_bar(aes(week,prop), stat = "identity", fill="black") +
  scale_y_continuous(breaks = seq(0,100,25),
                     limits = c(0,100),
                     sec.axis=sec_axis(~.*1,name="% Positive")) +
  theme_test()+
  ylab("Samples tested") +
  xlab("Week")+
  theme(axis.text.x = element_text(colour = "black"),
        axis.title = element_text(colour = "black", face = "bold"))+
  scale_x_continuous(breaks = seq(1,26,1))

supplementary_plot_1


# #Save output file 
#   ggsave("Result/global/supplementary_plot_1.tiff", 
#          height = 8, width = 10, dpi = 300)

##############################################################################
# Recombination aware graphs for graphs for XBB lineage 
#############################################################################


################################################################################
# Generate a segment of FY.4 
###############################################################################




segment1 <- read.tree("Data/filtered.segment1.fasta.treefile")


outliers <- c("hCoV-19/Kenya/C119499/2022", 
              "MN908947" , 
              "hCoV-19/Kenya/KNPHLGEN124/2022", 
              "hCoV-19/Kenya/KNPHLGEN084/2022", 
              "hCoV-19/Kenya/KNPHLGEN123/2022", 
              "hCoV-19/Kenya/SS11515/2022", 
              "hCoV-19/Kenya/SS11721/2022",
              "hCoV-19/Kenya/SS11724/2022"
)


segment1 <- ape::drop.tip(segment1, outliers)

segment1_metadata <-  read_tsv("Data/segment1_metadata.tsv")


# Assign color to specific lineages
cols <- c( "gold4", "aquamarine4", "steelblue3", "black", "plum4", "salmon4",  
           "indianred", "sandybrown", "deeppink4", "olivedrab", "lightskyblue4", 
           "navajowhite2", "palevioletred", "cadetblue3")



segment1_metadata <- segment1_metadata %>% 
  tidytree::filter(!strain %in% outliers) %>% 
  tidytree::select(strain, date, pangolin_lineage) %>%
  tidytree::mutate(lineage=case_when(str_starts(pangolin_lineage, "BA.1") ~ "BA.1*", 
                                     str_starts(pangolin_lineage, "FY.4") ~ "FY.4*", 
                                     str_starts(pangolin_lineage, "XBB.1.16") ~ "XBB.1.16*", 
                                     str_starts(pangolin_lineage, "BQ.1") ~ "BQ.1*", 
                                     str_starts(pangolin_lineage, "BA.4") ~ "BA.4*", 
                                     str_starts(pangolin_lineage, "XBB.1.5") ~ "XBB.1.5*", 
                                     str_starts(pangolin_lineage, "XBB.1.9|EG.5|FL|GW.5") ~ "XBB.1.9*",
                                     str_starts(pangolin_lineage, "GE.1|JE.1|XBB.2.3|GS.4") ~ "XBB.2.3*",
                                     str_starts(pangolin_lineage, "BA.5|BF") ~ "BA.5*", 
                                     str_starts(pangolin_lineage, "BA.2.75|BN.1|CH.1|BM.1") ~ "BA.2.75*", 
                                     str_starts(pangolin_lineage, "BA.2.86") ~ "BA.2.86*", 
                                     str_starts(pangolin_lineage, "BA.2") ~ "BA.2", 
                                     str_starts(pangolin_lineage, "JN.1|LE.1|BE.1") ~ "JN.1*", 
                                     str_starts(pangolin_lineage, "XBB.2|XBB.3|XBB.1.22|XBB.1.34|XBB.1.43|XBB|GA.6|FY.1|HH.1|GE.1|JE.1|XBB.2.3|GS.4") ~ "Other XBB*"),
                   lineage=as.factor(lineage)) %>%
  tidytree::mutate(date_class=format(date, "%b %Y"), 
                   date_class=as.factor(date_class)) %>%
  as.data.frame() 


levels(segment1_metadata$lineage)



segment1_tree <- ggtree(segment1, 
                        as.Date = T, 
                        color="grey80") %<+% segment1_metadata  +
  theme_tree() +
  geom_tippoint(aes(colour=lineage)) +
  
  geom_treescale(offset = -8)+
  #scale_x_date(date_breaks = "1 month", 
  #             date_labels = "%b-%Y" ) +
  labs(color="Lineage")  +
  theme(legend.title = element_text(colour = "black", 
                                    face = "bold")) +
  scale_color_manual(values = c("BA.2"= "gold4", 
                                "BA.2.75*" = "aquamarine4", 
                                "BA.4*" = "steelblue3" , 
                                "BA.5*" = "black", 
                                "BQ.1*" = "plum4",
                                "FY.4*" = "salmon4", 
                                "JN.1*" =  "indianred",
                                "Other XBB*" = "sandybrown", 
                                "XBB.1.16*" = "deeppink4", 
                                "XBB.1.5*"  = "olivedrab", 
                                "XBB.1.9*"  = "lightskyblue4", 
                                "XBB.2.3*"  = "navajowhite2" )) 


cols <- c( "gold4", "aquamarine4", "steelblue3", "black", "plum4", "salmon4",  
           "indianred", "sandybrown", "deeppink4", "olivedrab", "lightskyblue4", 
           "navajowhite2", "palevioletred", "cadetblue3")



segment3 <- read.tree("Data/filtered.segment3.fasta.treefile")

segment3 <- ape::drop.tip(segment3, outliers)

segment3_metadata <- read_tsv("Data/segment3_metadata.tsv")

segment3_metadata <- segment3_metadata %>% 
  tidytree::filter(!strain %in% outliers) %>% 
  tidytree::select(strain, date, pangolin_lineage) %>%
  tidytree::mutate(lineage=case_when(str_starts(pangolin_lineage, "BA.1") ~ "BA.1*", 
                                     str_starts(pangolin_lineage, "FY.4") ~ "FY.4*", 
                                     str_starts(pangolin_lineage, "XBB.1.16") ~ "XBB.1.16*", 
                                     str_starts(pangolin_lineage, "BQ.1") ~ "BQ.1*", 
                                     str_starts(pangolin_lineage, "BA.4") ~ "BA.4*", 
                                     str_starts(pangolin_lineage, "XBB.1.5") ~ "XBB.1.5*", 
                                     str_starts(pangolin_lineage, "XBB.1.9|EG.5|FL|GW.5") ~ "XBB.1.9*",
                                     str_starts(pangolin_lineage, "GE.1|JE.1|XBB.2.3|GS.4") ~ "XBB.2.3*",
                                     str_starts(pangolin_lineage, "BA.5|BF") ~ "BA.5*", 
                                     str_starts(pangolin_lineage, "BA.2.75|BN.1|CH.1|BM.1") ~ "BA.2.75*", 
                                     str_starts(pangolin_lineage, "BA.2.86") ~ "BA.2.86*", 
                                     str_starts(pangolin_lineage, "BA.2") ~ "BA.2", 
                                     str_starts(pangolin_lineage, "JN.1|LE.1|BE.1") ~ "JN.1*", 
                                     str_starts(pangolin_lineage, "XBB.2|XBB.3|XBB.1.22|XBB.1.34|XBB.1.43|XBB|GA.6|FY.1|HH.1|GE.1|JE.1|XBB.2.3|GS.4") ~ "Other XBB*"),
                   lineage=as.factor(lineage)) %>%
  tidytree::mutate(date_class=format(date, "%b %Y"), 
                   date_class=as.factor(date_class)) %>%
  as.data.frame() 



segment3_tree <- ggtree(segment3, 
                        color="grey80") %<+% segment3_metadata  +
  theme_tree() +
  geom_tippoint(aes(colour=lineage)) +
  #scale_color_manual(values = cols) +
  geom_treescale(offset = -8)+
  #scale_x_date(date_breaks = "1 month", 
  #             date_labels = "%b-%Y" ) +
  labs(color="Lineage")  +
  theme(legend.title = element_text(colour = "black", 
                                    face = "bold"), 
        legend.position = "none") +
  scale_color_manual(values = c("BA.2"= "gold4", 
                                "BA.2.75*" = "aquamarine4", 
                                "BA.4*" = "steelblue3" , 
                                "BA.5*" = "black", 
                                "BQ.1*" = "plum4",
                                "FY.4*" = "salmon4", 
                                "JN.1*" =  "indianred",
                                "Other XBB*" = "sandybrown", 
                                "XBB.1.16*" = "deeppink4", 
                                "XBB.1.5*"  = "olivedrab", 
                                "XBB.1.9*"  = "lightskyblue4", 
                                "XBB.2.3*"  = "navajowhite2" ))


#levels(segment3_metadata$lineage)


segments_analysis <- plot_grid(segment1_tree, segment3_tree, 
                               labels = c("A", "B"))

segments_analysis


###############################################################################

###############################################################################
# select samples for starting from May for phylogeography analysis 
################################################################################
start_date=as.Date("2023-05-01")

FY4_mayIDs <- metadata  %>%
  tidytree::filter(date>start_date) %>%
  tidytree::select(gisaid_epi_isl)

# wrtite out sequence IDS
write_tsv(FY4_mayIDs, 
          file = "Data/FY4_mayIDS.tsv")

# Load generated ML tree for may-jan2024 FY.4 sequences 

FY.4_may <- read.tree("Data/FY.4_may.treefile")


# write out metadata associated with the generated tree 
FY.4_may_metadata <-  metadata %>%
  tidytree::filter(date>start_date) %>%
  tidytree::select(strain, date, country) %>%
  mutate(location=case_when(country=="Kenya" ~ "Kenya", 
                            .default = "Others")) %>%
  tidytree::select(strain, date, location)



# Load tempest data 
may_fy.4_tempest <-  read_tsv("Data/FY.4_may_tempest_data.tsv")



# Run regression model with exported data 

model_subsampled_fy.4<- lm(distance~date, data = may_fy.4_tempest)

summary(model_subsampled_fy.4)




# Extract residuals
residuals_may <- resid(model_subsampled_fy.4)

# Subtract residuals from fitted values
adjusted_values_fy4 <- fitted(model_subsampled_fy.4) - residuals_may



# Refine temporal data 
may_fy.4_tempest<-  may_fy.4_tempest %>%
  mutate(date=as.Date(date_decimal(date)))

# add ajusted values to dataframe
adjusted_values_fy4 <- setNames(adjusted_values_fy4, "adjusted_values")
may_fy.4_tempest$adjusted_values <- adjusted_values_fy4



# Change names of the columns to match 

colnames(may_fy.4_tempest)[1:2] <- c("strain", "date")


FY.4_may_metadata <- FY.4_may_metadata %>%
  tidytree::select(strain, location)

finalized_temporal_may <- full_join(FY.4_may_metadata, 
                                    may_fy.4_tempest, 
                                    by="strain")

root_tip_plot_may <- ggplot(finalized_temporal_may, aes(date, distance, 
                                                        color=location)) +
  geom_point()+
  geom_point(aes(y=adjusted_values))+
  geom_smooth(method = "lm", linetype="dashed", color="black") +
  theme_classic()+
  xlab("Sampling Date") + 
  ylab("Root-to-tip divergence") +
  labs(color="Country") +
  theme(axis.text.x = element_text(colour = "black", angle = 90), 
        axis.title = element_text(colour = "black", face = "bold"), 
        legend.title = element_text(colour = "black", face = "bold")) +
  scale_y_continuous(labels = seq(0,0.007, 0.0001), 
                     breaks = seq(0,0.007, 0.0001), 
                     limits = c(0,0.0008)) +
  scale_x_date(date_breaks = "1 month", date_labels = "%b-%Y" ) +
  scale_color_manual(values = cols)+
  annotate(geom="text",
           x=as.Date("2023-06-30"), 
           y=0.00075, 
           label="Correlation coefficient = 0.63") +
  annotate(geom = "text",
           x=as.Date("2023-06-30"), 
           y=0.0007, 
           label="R-squared = 0.34") +
  theme(plot.title = element_text(colour = "black", hjust = 0.5), 
        legend.position = "bottom", 
        legend.box = "horizontal")


root_tip_plot_may

###############################################################################
# Supplementary 3 file showing local movement of FY.4 
################################################################################


may_FY.4 <- read.beast("Data/MCC_mayFY.4.tree")



may_preliminary_DTA <-  ggtree(may_FY.4,
                               aes(colour = location), 
                               mrsd = most_recent_sampling_date, 
                               as.Date = T) +
  theme_tree2()+
  scale_color_manual(values = new_colors) +
  geom_nodepoint(aes(node = !isTip, fill= location), size=3, shape=21) +
  scale_fill_manual(values = c("plum4", "salmon4"))+
  labs(fill="Ancestral node") +
  labs(color="Location") +
  theme(legend.title = element_text(colour = "black", face = "bold"), 
        axis.text.x = element_text(angle = 90, colour = "black")) +
  ylim(0,940) +
  scale_x_date(date_breaks = "1 month", 
               date_labels = "%b-%Y" )


may_preliminary_DTA
transitions_may_FY.4 <-  read.delim("Data/kenya_may_FY.4.out", 
                                    sep = "\t")


transitions_may_FY.4 <- transitions_may_FY.4[1:4, 1:4]


colnames(transitions_may_FY.4)[1:4] <- c("Transition", "mean", 
                                         "95% Low", "95% High")

may_transitions <-transitions_may_FY.4 %>%
  tidytree::select(Transition, mean) %>%
  tidytree::mutate(mean=round(mean,0)) %>%
  separate(Transition, sep = "=>", into = c("Origin", "Destination")) %>% 
  ggplot(., aes(Origin, Destination, fill = mean))+
  geom_tile(color="black") +
  geom_text(aes(label=mean)) +
  scale_fill_gradient(low = "gold", high = "aquamarine4") +
  labs(fill="Transitions")+
  theme(axis.text = element_text(colour = "black"), 
        axis.title = element_text(colour = "black", face = "bold"), 
        legend.title = element_text(colour = "black", face = "bold"))



may_transitions


supp3A <- plot_grid(root_tip_plot_may, may_transitions, 
                    labels = c("A", "B"), 
                    ncol = 1, 
                    rel_heights = c(1, 0.5), 
                    scale = 0.95)

supp3 <- plot_grid(supp3A, 
                   may_preliminary_DTA, 
                   labels = c(" ", "C"))


supp3

ggsave("Result/global/supp3.tiff", 
       height = 7, 
       width = 10)

