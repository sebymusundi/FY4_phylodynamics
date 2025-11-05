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

# Define colors to apply in subsequent images 

cols <- c( "gold4", "aquamarine4", "steelblue3", "black", "plum4", "salmon4",  
           "indianred", "sandybrown", "deeppink4", "olivedrab", "lightskyblue4", 
           "navajowhite2", "palevioletred", "cadetblue3")

###############################################################################
# Read kenyan_metadata file 
kenyan_metadata <- read_tsv("Data/kenyan_metadata.tsv")

# Make variable assigning dates 

# arrange dates  for distribution of  metadata from GISAID 

dates <- c( 
  "Mar 2023", "Apr 2023", "May 2023", "Jun 2023", "Jul 2023", "Aug 2023", 
  "Sep 2023", "Oct 2023", "Nov 2023", "Dec 2023", "Jan 2024")

# Plot frequency of the lineage from November 2022 to Jan 2024

kenyan_metadata <- kenyan_metadata %>%
  tidytree::filter(!date_class %in% c ("Sep 2022", "Oct 2022", "Nov 2022", 
                                       "Dec 2022", "Jan 2023", "Feb 2023" ))


# indicate the number of genomes that were deposited each month in Kenya 
genomes_numbers <- kenyan_metadata%>% 
  group_by(date_class) %>%
  summarise(count=n(), 
            .groups = "drop")


#str(genomes_numbers)

lineage_freq <- kenyan_metadata %>%
  tidytree::select(lineage, date_class) %>%
  group_by(lineage, date_class) %>%
  summarise(count=n()) %>%
  ggplot(., aes(date_class, count, fill=lineage))+
  geom_bar(stat = "identity", position = "fill") +
  geom_text(data = genomes_numbers, aes(x = date_class, y = 1, label = count), 
            inherit.aes = FALSE, vjust = -0.5, color = "black", size = 3.0) +
  xlab("Time") +
  ylab("Frequency") +
  labs(fill="Lineage") +
  theme_test() +
  theme(axis.text.x = element_text(colour = "black", angle = 90), 
        axis.text.y = element_text(colour = "black"), 
        axis.title = element_text(hjust = 0.5, color="black", face = "bold"), 
        plot.title = element_text(color = "black", hjust = 0.5, face = "bold"), 
        legend.title = element_text(color = "Black", face = "bold"))+
  scale_fill_manual(values = cols)+
  scale_x_discrete(limits=c(dates))



lineage_freq


#save the data
# ggsave("Result/global/lineage_freq.tiff", height = 6, width = 8, 
#        dpi = 300)



###############################################################################
# global distribution of FY4
###############################################################################

early_nov_seq <-  c("hCoV-19/Kenya/CDC-DLSP-Kisumu-COV223824/2022", 
                    "hCoV-19/Kenya/CDC-DLSP-Kisumu-COV223830/2022", 
                    "hCoV-19/Kenya/CDC-DLSP-Kisumu-COV223972/2022")

# Retrieve the global metadata from gisaid 

complete_metadata <- read_csv("Data/FY4_global_metadata.csv", 
                              col_types = cols(date=col_date(format = "%Y-%m-%d"))) %>%
  as.data.frame() %>%
  tidytree::filter(!strain %in% early_nov_seq)



# Arrange months in order 

month_order <- c( "Mar 2023",
                  "Apr 2023", "May 2023", "Jun 2023", "Jul 2023", 
                  "Aug 2023", "Sep 2023", "Oct 2023", "Nov 2023", "Dec 2023", 
                  "Jan 2024")

#  summarize global distribution of FY.4
FY4_data <-  complete_metadata %>%
  tidytree::select(strain, date, region,pangolin_lineage) %>% 
  tidytree::mutate(region=as.factor(region), 
                   pangolin_lineage=as.factor(pangolin_lineage), 
                   month_year=format(date, "%b %Y"), month_year=as.factor(month_year)) %>%
  na.omit() %>%
  group_by(region,pangolin_lineage, month_year) %>% 
  summarize(count=n()) %>%
  as.data.frame()


# Plot graph showing the distribution of FY4
FY4_distribution <- ggplot(FY4_data, aes(month_year, region, color=pangolin_lineage,
                                         size=count))+
  geom_point(alpha=0.8, shape=16) +
  theme_test() +
  xlab("Week") +
  ylab("Region") +
  labs(color="Pangolin lineage")+
  theme(axis.text.y=element_text(color = "black"), 
        axis.text.x=element_text(color = "black", angle = 90),
        axis.title = element_text(colour = "black", face = "bold"), 
        plot.title = element_text(color = "black", hjust = 0.5, face = "bold")) +
  scale_color_manual(values = cols) +
  scale_x_discrete(limits=month_order) +
  scale_size(range = c(2,15), name = "Number of genomes") 


FY4_distribution
# ggsave("Result/global/FY4_distribution.png", 
#        height = 6, width = 12,dpi = 300)


# combined local and global distribution of FY4 

FY.4_lineage_summary <- plot_grid(lineage_freq, 
                                  FY4_distribution,
                                  ncol = 1, nrow = 2, 
                                  labels = c ("A", "B"), 
                                  scale = 0.9)

################################################################################
# Spatial temporal analysis of FY.4 sequences 
################################################################################

# samples with incomplete collection date 
no_dates <- c("hCoV-19/England/CLIMB-CM7YG5X7/2023", 
              "hCoV-19/USA/OH-PLMI-HJCNT-20194/2023", 
              "hCoV-19/USA/un-NHIE001871N/2023", 
              "hCoV-19/England/CLIMB-CM7YEEZX/2023", 
              "hCoV-19/England/CLIMB-CM7Y8Y3I/2023", 
              "hCoV-19/USA/un-NHIE001884N/2023", 
              "hCoV-19/England/CLIMB-CM7Y8XQ9/2023", 
              "hCoV-19/USA/OH-PLMI-HJC5G-19976/2023", 
              "hCoV-19/USA/un-NHIE002457N/2023", 
              "hCoV-19/England/CLIMB-CM7YF9QX/2023", 
              "hCoV-19/SaudiArabia/KFSHRC_28AF/2023")

less_coverage <- "hCoV-19/Kenya/C122561/2023"


not_FY.4 <- c("hCoV-19/USA/TX-HMH-M-132390/2023", 
              "hCoV-19/USA/VA-CDC-LC1045259/2023", 
              "hCoV-19/USA/CA-HLX-STM-GNCG28RMC/2023", 
              "hCoV-19/Kenya/SS11907/2023", 
              "hCoV-19/Kenya/CDC-DLSP-Kisumu-COV233835/2023")

raw_tree_outliers <- c("hCoV-19/USA/CA-HLX-STM-NRKZZT4Q9/2023", 
                       "hCoV-19/USA/TX-HMH-M-133133/2023", 
                       "hCoV-19/Canada/SK-RRPL-657955/2023", 
                       "hCoV-19/USA/MS-MSPHL-0319/2023", 
                       "hCoV-19/USA/IL-S23WGS2025/2023", 
                       "hCoV-19/USA/KY-KSPHL-23822/2023", 
                       "hCoV-19/USA/VA-CAV_VAS3N_00019310_01/2023", 
                       "hCoV-19/Spain/IB-HUSE-09211/2024", 
                       "hCoV-19/USA/CA-HLX-STM-9ZUNUW3RB/2023", 
                       "hCoV-19/USA/CA-HLX-STM-GNCG28RMC/2023", 
                       "hCoV-19/Kenya/C122615/2023", 
                       "hCoV-19/Kenya/SS11821/2023",
                       "hCoV-19/USA/TX-HHD-2308094526/2023", 
                       "hCoV-19/USA/TX-HHD-2308230312/2023", 
                       "hCoV-19/USA/MN-HLX-VSX-A04503/2023", 
                       "hCoV-19/Japan/PG-561600/2023")



# all outlier sequences IDs
all_outliers <- c(no_dates, less_coverage,
                  not_FY.4, raw_tree_outliers) 

# include reference outliers

reference_outlier <- c("XBB", "hCoV-19/Kenya/C123075/2023", early_nov_seq, 
                       all_outliers)

# load Original tree containing all sequences 
results_all <- read.tree("Data/cleaned_FY.4.treefile")


# remove all outliers before proceeding to generate a time-resolved phylogeny
results_all <-  ape::drop.tip(results_all, reference_outlier)

results_all

ape::write.tree(results_all, "Data/cleaned_FY.4ML.nwk" )


# load metadata
metadata<- complete_metadata %>%
  tidytree::filter(!strain %in% reference_outlier)


write_tsv(metadata, 
          file = "Data/metadata_timetree_FY.4.tsv")

# Generate output files for timetree



# Load timetree outliers  
#timetree_outliers <- read_tsv("Data/FY4_outliers.tsv") %>%
#  tidytree::select(1)%>%
#  tidytree::pull()


# append all outlier sequences 
#timetree_outliers <- c("XBB", timetree_outliers)


# Load time-resolved timetree 
timetree <- read.nexus("Data/FY.4_timetree.nexus")

timetree


# remove outliers 
#timetree <- ape::drop.tip(timetree, timetree_outliers)



# timetree metadata 

#timetree_metadata <-  metadata %>%
#  tidytree::filter(!strain %in% c(timetree_outliers))

# retrive the most recent sampling date 

most_recent_sampling_date <- metadata %>%
  tidytree::select(date) %>%
  tidytree::pull(date) %>%
  max()

# Classify locations as Kenya and global 

timetree_metadata <-  metadata %>%
  tidytree::mutate(global_locations=case_when(country=="Kenya" ~ "Kenya", 
                                              .default = "Global"), 
                   global_locations=as.factor(global_locations))



# plot time-resolved tree 


timetree

timetree_locations <- ggtree(timetree, color="grey80",
                             mrsd = most_recent_sampling_date, 
                             as.Date = T) %<+% timetree_metadata + 
  theme_tree2() +
  geom_tippoint(aes(color=global_locations, shape = pangolin_lineage), size=2) +
  scale_color_manual(values = cols) +
  # scale_color_brewer(palette = "Paired")+
  scale_x_date(date_breaks = "1 month", 
               date_labels = "%b-%Y" ) +
  scale_shape_manual(values = c(15,16,18,24, 23)) +
  labs(color="Region", shape= "Pangolin lineage") +
  theme(axis.text.x = element_text(angle = 90, colour = "black"), 
        legend.title = element_text(colour = "black", face = "bold")) +
  ylim(0,935) 


timetree_locations

#ggsave("Result/global/timetree_locations.tiff", 
#       height = 6, 
#     width = 8)


Fig1_summary <- plot_grid(FY.4_lineage_summary,
                          timetree_locations, labels = c ("", "C"))

Fig1_summary

#ggsave("Result/global/Fig1_summary.pdf", 
#       height = 12, width = 16)


