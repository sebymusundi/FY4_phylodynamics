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

################################################################################
# preliminary phylgeography analysis using the main ancestral kenyan node

# extract kenyan samples only and kenyan timetree 
kenyan_samples_dates <- timetree_metadata %>%
  tidytree::filter(country=="Kenya") %>%
  tidytree::select(strain,date) 


write_tsv(kenyan_samples_dates, 
          file = "Data/kenya_samples_dates.tsv")


coordinates_kenya <- timetree_metadata %>%
  tidytree::filter(country=="Kenya") %>%
  tidytree::mutate( coords = case_when(division == "Embu" ~ list(c(0.53112, 37.45061)), 
                                       division == "Kakamega" ~ list(c(0.28422, 34.75229)), 
                                       division == "Kiambu" ~ list(c(-1.16667, 36.83333)),
                                       division == "Kilifi" ~ list(c(-3.63045, 39.84992)),
                                       division == "Kisumu" ~ list(c(-0.10221, 34.76171)),
                                       division == "Lamu" ~ list(c(-2.2666656, 40.916663)), 
                                       division == "Machakos" ~ list(c(-1.51667, 37.26667)), 
                                       division == "Makueni" ~ list(c(-1.80409, 37.62034)), 
                                       division == "Migori" ~ list(c(-1.06343, 34.47313)),
                                       division == "Mombasa" ~ list(c(-4.05466, 39.66359)),
                                       division == "Nairobi" ~ list(c(-1.28333, 36.81667)), 
                                       division == "Nakuru" ~ list(c(-0.2833322, 36.0666664)),
                                       division == "Narok" ~ list(c(-1.07829, 35.86012)),
                                       division == "Nyeri" ~ list(c(-0.42013, 36.94759)),
                                       division == "Trans Nzoia" ~ list(c(1.045, 34.979)), 
                                       division =="Busia" ~ list(c(0.451831526, 34.12166618)), 
                                       division== "Laikipia" ~ list(c(0.3606, 36.7820)), 
                                       division== "Siaya Surveillance" ~ list(c(0.06116, 34.28823))
                                       
  ), 
  latitude = map_dbl(coords, 1),  # Extract the first element (latitude)
  longitude = map_dbl(coords, 2)  # Extract the second element (longitude)
  ) %>%
  tidytree::select(-coords) %>%
  tidytree::select(strain, latitude, longitude)

write_tsv(coordinates_kenya, 
          file = "Data/coordinates_kenya.tsv")


# retrieve newick tree using kenyan samples

kenya_only_sequences <- timetree_metadata %>%
  tidytree::filter(!country=="Kenya") %>%
  tidytree::select(strain) %>%
  pull()

timetree_kenya <- ape::drop.tip(timetree, kenya_only_sequences)

ape::write.tree(timetree_kenya, 
                file = "Data/timetree_kenya.nwk")
################################################################################

kenyan_samples_dates %>%
  tidytree::select(date) %>%
  pull(date) %>%
  max() %>%
  decimal_date()

#extracting spatio-temporal information embedded in posterior trees


localTreesDirectory = "Extracted_trees"
allTrees = scan(file="Data/continous_sampled.trees",
                what="", sep="\n", quiet=TRUE)
burnIn = 0
randomSampling = FALSE
nberOfTreesToSample = 1000
mostRecentSamplingDatum = 2023.693
coordinateAttributeName = "coordinates"



treeExtractions(localTreesDirectory, allTrees, burnIn, randomSampling,
                nberOfTreesToSample, mostRecentSamplingDatum, coordinateAttributeName)


#extracting spatio-temporal information embedded in the MCC tree
mcc_tre <- readAnnotatedNexus("Data/continuous_mcc.tree")


source("Script/mccExtractions.r")

mcc_tab = mccTreeExtraction(mcc_tre, mostRecentSamplingDatum)



#estimating the HPD region for each time slice

nberOfExtractionFiles = nberOfTreesToSample
prob = 0.95
precision = 0.025
startDatum = min(mcc_tab[,"startYear"])

polygons = suppressWarnings(spreadGraphic2(localTreesDirectory,
                                           nberOfExtractionFiles, prob, startDatum, precision))
polygons
#defining the different colour scales to use

colour_scale = colorRampPalette(brewer.pal(11,"RdYlGn"))(141)[21:121]
minYear = min(mcc_tab[,"startYear"]); maxYear = max(mcc_tab[,"endYear"])
endYears_indices = (((mcc_tab[,"endYear"]-minYear)/(maxYear-minYear))*100)+1
endYears_colours = colour_scale[endYears_indices]
polygons_colours = rep(NA, length(polygons))
for (i in 1:length(polygons)) {
  date = as.numeric(names(polygons[[i]]))
  polygon_index = round((((date-minYear)/(maxYear-minYear))*100)+1)
  polygons_colours[i] = paste0(colour_scale[polygon_index],"40") }

#co-plotting the HPD regions and MCC tree

template_raster = raster("Data/map_kenya/DIF.asc", 
                         xmn=33.90959, xmx=41.92622, ymn=-4.720417, ymx=5)

template_raster

# Get GADM data for Kenya using geodata
# borders = crop(getData("GADM", country="BRA", level=1), extent(template_raster))
border=gadm(country = "KEN", level = 1, path=tempdir())

#plot(border)

# Crop the borders to the extent of the template raster
cropped_borders = crop(border, extent(template_raster))




#plot(cropped_borders, add = TRUE)


dev.new(width=6.0, height=6.3)
par(mar=c(0,0,0,0), oma=c(1.2,3.5,1,0), mgp=c(0,0.4,0), lwd=0.2, bty="o")
plot(template_raster, col="white", box=F, axes=F, colNA="grey90", legend=F)


for (i in 1:length(polygons)) {
  plot(polygons[[i]], axes=F, col=polygons_colours[i], add=T, cropped_borders=NA) }

plot(cropped_borders, add=T, lwd=0.1, border="gray10")


for (i in dim(mcc_tab)[1]:1) {
  curvedarrow(cbind(mcc_tab[i,"startLon"],mcc_tab[i,"startLat"]),
              cbind(mcc_tab[i,"endLon"],mcc_tab[i,"endLat"]), arr.length=0,
              arr.width=0, lwd=0.2, lty=1, lcol="gray10", arr.col=NA,
              arr.pos=F, curve=0.1, dr=NA, endhead=F)
  
}  


for (i in dim(mcc_tab)[1]:1) {
  if (i == 1) {
    points(mcc_tab[i,"startLon"], mcc_tab[i,"startLat"], pch=16,
           col=colour_scale[1], cex=0.8)
    points(mcc_tab[i,"startLon"], mcc_tab[i,"startLat"], pch=1,
           col="gray10", cex=0.8)
  }
  points(mcc_tab[i,"endLon"], mcc_tab[i,"endLat"], pch=16,
         col=endYears_colours[1], cex=0.8)
  points(mcc_tab[i,"endLon"], mcc_tab[i,"endLat"], pch=1,
         col="gray10", cex=0.8)
}


rect(xmin(template_raster), ymin(template_raster), xmax(template_raster), 
     ymax(template_raster), xpd=T, lwd=0.2)

axis(1, c(ceiling(xmin(template_raster)), floor(xmax(template_raster))),
     pos=ymin(template_raster), mgp=c(0,0.2,0), cex.axis=0.5, lwd=0, lwd.tick=0.2,
     padj=-0.8, tck=-0.01, col.axis="gray30")

axis(2, c(ceiling(ymin(template_raster)), floor(ymax(template_raster))),
     pos=xmin(template_raster), mgp=c(0,0.5,0), cex.axis=0.5, lwd=0, lwd.tick=0.2,
     padj=1, tck=-0.01, col.axis="gray30")

rast = raster(matrix(nrow=1, ncol=2))
rast[1] = min(mcc_tab[,"startYear"])
rast[2] = max(mcc_tab[,"endYear"])
plot(rast, legend.only=T, add=T, col=colour_scale, legend.width=0.5, legend.shrink=0.3,
     smallplot=c(0.10,0.40,0.14,0.155), legend.args=list(text="", cex=0.7, line=0.3,
                                                         col="gray30"), horizontal=T, 
     axis.args=list(cex.axis=0.6, lwd=0, lwd.tick=0.2,
                    tck=-0.5, col.axis="gray30", line=0, mgp=c(0,-0.02,0), at=seq(2022.6,2023.8,0.2)))




