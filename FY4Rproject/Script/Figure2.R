# clear working environment

#rm(list = ls())


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

#############################################################################
# Tempest global analysis
############################################################################
tempest_data_subsampled <- read_tsv("Data/tempesT_export_data.tsv", 
                                    show_col_types = FALSE)


# Refine metadata 
tempest_metadata <- metadata %>%
  tidytree::filter (!strain %in% reference_outlier) %>%
  tidytree::mutate(country=as.factor(country), 
                   division=as.factor(division), 
                   region=as.factor(region), 
                   pangolin_lineage=as.factor(pangolin_lineage))


# Run regression model with exported data 

model_subsampled<- lm(distance~date, data = tempest_data_subsampled)

summary(model_subsampled)




# Extract residuals
residuals <- resid(model_subsampled)

# Subtract residuals from fitted values
adjusted_values <- fitted(model_subsampled) - residuals



# Refine temporal data 
tempest_data_subsampled<-  tempest_data_subsampled%>%
  mutate(date=as.Date(date_decimal(date)))

# add ajusted values to dataframe
adjusted_values <- setNames(adjusted_values, "adjusted_values")
tempest_data_subsampled$adjusted_values <- adjusted_values



# Change names of the columns to match 

colnames(tempest_data_subsampled)[1:2] <- c("strain", "date_tempest")

finalized_temporal <- full_join(tempest_metadata, 
                                tempest_data_subsampled, 
                                by="strain")

finalized_temporal <- finalized_temporal %>%
  mutate(county_class=case_when(country=="Kenya" ~ "Kenya", 
                                TRUE ~ "Global"))

root_tip_plot <- ggplot(finalized_temporal, aes(date, distance, 
                                                color=county_class)) +
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
           x=as.Date("2023-05-30"), 
           y=0.00075, 
           label="Correlation coefficient = 0.58") +
  annotate(geom = "text",
           x=as.Date("2023-05-30"), 
           y=0.0007, 
           label="R-squared = 0.34") +
  theme(plot.title = element_text(colour = "black", hjust = 0.5), 
        legend.position = "bottom", 
        legend.box = "horizontal")


root_tip_plot


######################################################################################
####################################################################################
#Fixed topology tree for beast analysis 
timetree_beast <- ape::write.tree(timetree, 
                                  file = "Data/timetree_all_beast.nwk")

# Dates file
timetree_beast_dates <- timetree_metadata %>%
  tidytree::select(strain,date)

write_tsv(timetree_beast_dates, 
          file = "Data/timetree_beast_dates.tsv")


# Traits file with kenya and other locations 

timetree_trait <- timetree_metadata %>%
  tidytree::select(strain,country) %>%
  tidytree::mutate(location=case_when(country=="Kenya" ~ "Kenyan", 
                                      .default = "Other")) %>%
  tidytree::select(strain,location)

write_tsv(timetree_trait, 
          file = "Data/timetree_trait.tsv")
#####################################################################################
# Identify kenyan clades from discrete trait analysis tree

clades_kenya <- read.beast("Data/DTA_mcc.tree")

# Convert to table
clades_kenya_table <- as_tibble(clades_kenya)

# assign colors to the tree 

new_colors <- c( "black", "grey90")

# plot tree 


preliminary_DTA <- ggtree(clades_kenya,
                          aes(colour = location), 
                          mrsd = most_recent_sampling_date, 
                          as.Date = T) +
  theme_tree2()+
  scale_color_manual(values = new_colors) +
  geom_nodepoint(aes(node = !isTip, fill= location), size=3, shape=21) +
  scale_fill_manual(values = c("plum4", "salmon"))+
  labs(fill="Ancestral node") +
  labs(color="Location") +
  theme(legend.title = element_text(colour = "black", face = "bold"), 
        axis.text.x = element_text(angle = 90, colour = "black")) +
  ylim(0,940) +
  scale_x_date(date_breaks = "1 month", 
               date_labels = "%b-%Y" )

# Save output file 
preliminary_DTA
# ggsave("Result/global/preliminary_DTA.tiff", 
#        height = 7, 
#        width = 10)


# number of transitions 
transition_counter <- read.delim("Data/transitions_DTA_tree.out", 
                                 sep = "\t") 

transition_counter <- transition_counter[1:4, 1:4]


# rename columns 

colnames(transition_counter)[1:4] <- c("Transition", "mean", "95% Low", "95% High")




prem_transitions <- transition_counter %>%
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

prem_transitions 

Figure_2AB_plot <- plot_grid(root_tip_plot, prem_transitions, 
                             labels = c("A", "B"), 
                             ncol = 1, 
                             rel_heights = c(1, 0.5), 
                             scale = 0.95)



Figure_2_plot <- plot_grid(Figure_2AB_plot, preliminary_DTA, 
                           labels = c("", "C"))



Figure_2_plot
# ggsave("Result/global/Figure_2_plot.tiff", 
#        height = 7, 
#        width = 10)