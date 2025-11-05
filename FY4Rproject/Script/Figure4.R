# clear working environment


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
# DISCRETE TRAIT ANALYSIS 
###########################################################################
# Number of transitions
markov_jumps <- read_csv("Data/markov_mcc_complete.csv")

#str(markov_jumps)

# Plot the number of markov jumps
markov_jumps_data <- markov_jumps %>%
  tidytree::mutate(startLocation=as.factor(startLocation), 
                   endLocation=as.factor(endLocation), 
                   period_time=date_decimal(time), 
                   period_time=as.Date(period_time,  format="%Y-%m-%d"),
                   year_month=format(period_time, "%b %Y"), 
                   week_number=week(period_time)) %>%
  tidytree::select(startLocation, endLocation, year_month, period_time,
                   week_number) %>%
  tidytree::filter(startLocation=="Kenya") %>%
  group_by(endLocation, year_month) %>%
  summarize(count=n()) %>%
  as.data.frame()


# assign dates to be included in x-axis
dates_scaled <- c("Feb 2023", "Mar 2023", "Apr 2023", "May 2023", "Jun 2023", "Jul 2023", "Aug 2023")

# plot number of markov jumps
markov_jump_plot <-  ggplot(markov_jumps_data, aes(year_month, count, fill=endLocation)) +
  geom_bar(stat = "identity") +
  theme_test() +
  xlab("Time") +
  ylab("Number of exports") +
  scale_fill_manual(values = c("gold4", "aquamarine4","black"))+
  scale_y_continuous(limits = c(0,20), 
                     breaks = seq(0,20,5), 
                     labels = seq(0,20,5)) +
  scale_x_discrete(limits=dates_scaled) +
  labs(fill="Region") +
  theme(axis.text.x = element_text(colour = "black", angle = 90), 
        axis.title = element_text(colour = "black", face = "bold"), 
        legend.title = element_text(colour = "black", face = "bold")) 


markov_jump_plot 


##############################################################################
# Change in effective population size 

population_size <- read_tsv("Data/population_size_uncorrelated.tsv", 
                            show_col_types = TRUE, skip = 1) %>%
  as.data.frame() %>%
  mutate(mean=as.numeric(mean), 
         median=as.numeric(median), 
         upper=as.numeric(upper), 
         lower=as.numeric(lower))



# Set the start date based on the tree’s time range (adjust if needed)
start_date <- as.Date(date_decimal(2022.928))
end_date <-  as.Date(date_decimal(2023.967))


# Plot the effective population size plot

population_size_plot <- ggplot(population_size, aes(date)) +
  geom_ribbon(aes(ymin=lower, ymax = upper), fill="steelblue") +
  geom_line(aes(y=median), linewidth=1) +
  ylab("Effective population size (Ne)") +
  xlab("") +
  scale_y_log10() +
  # geom_vline(xintercept = 2023.1171, linetype="dashed") +
  # geom_vline(xintercept = 2023.032)+
  theme_test() +
  theme(axis.text.x = element_blank(),
        axis.title = element_text(color = "black", face = "bold"))+
  scale_x_date(date_breaks = "1 months", 
               date_labels = "%b-%Y", 
               limits = c(start_date, end_date)) 

population_size_plot



###############################################################################
# DISCRETE TRAIT ANALYSIS 
DTA_metadata <-  timetree_metadata %>%
  tidytree::filter(!country =="Uganda" & !region=="South America") %>% 
  tidytree::mutate(region_new=case_when(country=="Kenya" ~ "Kenya", 
                                        region %in% c ("Asia", "Oceania") ~ "Asia-Pacific",
                                        region == "Europe" ~ "Europe", 
                                        region=="North America" ~ "North America"), 
                   region_new <- as.factor(region_new))




# Discrete trait analysis file global

DTA_mcc <- read.beast("Data/subsampled_DTA_mcc.tree")


DTA_MCC_TREE_PLOT  <- ggtree(DTA_mcc,
                             as.Date = T, 
                             mrsd = most_recent_sampling_date,
                             aes(color=location)) +
  theme_tree2() +
  scale_color_manual(values = cols) +
  geom_tippoint() +
  ylim(0,520)+
  xlab("Time") +
  ylab("") +
  #theme_test() +
  scale_x_date(date_breaks = "1 months", 
               date_labels = "%b-%Y", 
               limits = c(start_date, end_date))+
  labs(color="Region") +
  theme(axis.text.x = element_text(colour = "black", angle = 90),
        plot.title = element_text(colour = "black",  hjust = 0.5, face="bold"), 
        axis.title.x = element_text(color = "black", face = "bold"), 
        legend.title = element_text(color = "black", face = "bold"), 
        legend.position = "bottom") 


DTA_MCC_TREE_PLOT 



grid_2 <- plot_grid( population_size_plot, DTA_MCC_TREE_PLOT,
                     labels = c("A", "B"), 
                     ncol = 1, 
                     rel_heights = c(1,2), 
                     scale = 0.98, 
                     label_size = 10)

markov_jump_plot_edit <- plot_grid(markov_jump_plot, 
                                   labels = "C", 
                                   label_size = 10, 
                                   hjust = 0.5, 
                                   vjust = 0.5, 
                                   scale = 0.7, 
                                   label_x = 0.1, 
                                   label_y = 0.9)



markov_jump_plot_edit

Figure_4_plot <- plot_grid(grid_2, markov_jump_plot_edit, 
                           rel_heights = c(1,1.5))



# plot overall figure 
Figure_4_plot





# Examined loaded tree 
DTA_mcc_table <-  as_tibble(DTA_mcc)

 x <- DTA_mcc_table %>%
  filter(label %in% str_detect(label, "hCoV")) %>%
  select(parent, node, location.set, location.set.prob) %>%
  unnest_wider(location.set, names_sep = ",") %>%
  mutate(location.set.prob=lapply(location.set.prob, as.numeric)) %>%
   unnest_wider(location.set.prob, names_sep = ",", 
                names_repair = "universal") %>%
   mutate(dominant_locations=case_when(
       location.set.prob.1 > 0.5 ~ location.set.1,
       location.set.prob.2 > 0.5 ~ location.set.2,
       location.set.prob.3 > 0.5 ~ location.set.3,
       location.set.prob.4 > 0.5 ~ location.set.4
   )) %>%
   select(parent,node, dominant_locations)
 
 
 transitions <-  x %>%
   select(child_node=node, 
          child_location=dominant_locations) %>%
   left_join(x%>%
               )


