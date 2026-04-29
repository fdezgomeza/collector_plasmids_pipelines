#!/usr/bin/env Rscript

# ==============================================================================
# Collector Plasmid Trajectory Visualizer (Chord Diagram)
#
# Description:
#   This script generates a chord diagram to visualize the  
#   trajectories (taxonomic transitions) of a Collector plasmid. It takes the 
#   LCA taxonomic assignments, filters them to the genus level, and determines 
#   the chronological transitions between consecutive spacers within each array. 
#   Transitions within the same genus are excluded.
#
# Author: Andrea Fernández
# Date: April 2026
# ==============================================================================

library(tidyverse)
library(circlize)
library(viridis)

# --- 1. Load and process data ---
df <- read_tsv("file.tsv", 
               col_names = c("QueryID", "LCA_TaxID", "Scientific_Name", "Rank", "Lineage"),
               skip = 1) 

processed_df <- df %>%
  separate(QueryID, into = c("ArrayID", "SpacerID"), sep = "_", convert = TRUE) %>%
  mutate(
    Taxon_Name = str_extract(Lineage, "[^;]+$"),
    Genus = case_when(
      str_detect(Lineage, "Klebsiella/Raoultella group") ~ "Klebsiella/Raoultella",
      str_detect(Taxon_Name, "Klebsiella pneumoniae complex") ~ "Klebsiella",
      str_detect(Taxon_Name, "unclassified ") ~ str_trim(str_remove(Taxon_Name, "unclassified ")),
      Rank == "genus" ~ Taxon_Name,
      Rank == "species" ~ word(Taxon_Name, 1), 
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(Genus))

# --- 2. Determine trajectories ---
trajectories <- processed_df %>%
  arrange(ArrayID, desc(SpacerID)) %>%
  group_by(ArrayID) %>%
  mutate(Next_Genus = lead(Genus)) %>%
  ungroup() %>%
  # Filter out terminal spacers and auto-transitions (same genus)
  filter(!is.na(Next_Genus), Genus != Next_Genus)

flow_summary <- trajectories %>%
  count(Source = Genus, Target = Next_Genus, name = "Frequency")

# Optional: View the summary table in the RStudio viewer
# view(flow_summary)

# --- 3. Colors ---
sectors <- unique(c(flow_summary$Source, flow_summary$Target))
set.seed(9) 
# Generate colors using the turbo palette
palette_colors <- sample(turbo(length(sectors), begin = 0.1, end = 0.95))
names(palette_colors) <- sectors

# --- 4. Dynamic parameters ---
total_flow <- sum(flow_summary$Frequency) * 2 
width_threshold <- total_flow * 0.04 

# --- 5. DRAW PLOT ---
png("ChordDiagram_Plasmid_Flows.png", width = 14, height = 14, units = "in", res = 300)

circos.clear()
circos.par(gap.after = 4, canvas.xlim = c(-1.4, 1.4), canvas.ylim = c(-1.4, 1.4))

chordDiagram(
  x = flow_summary,
  grid.col = palette_colors,
  directional = 1,
  direction.type = "diffHeight", 
  diffHeight = -0.05,            
  transparency = 0.35,
  annotationTrack = "grid", 
  preAllocateTracks = list(track.height = 0.25) 
)

# TRACK 1: Dynamic names (Shifted outward and larger)
circos.track(track.index = 1, panel.fun = function(x, y) {
  sector_width <- CELL_META$xlim[2] - CELL_META$xlim[1]
  
  if(sector_width > width_threshold) {
    # Names parallel to the ring for large sectors
    circos.text(x = CELL_META$xcenter, 
                y = CELL_META$ylim[1] + mm_y(6), 
                labels = CELL_META$sector.index, 
                facing = "bending.inside", 
                niceFacing = TRUE, 
                cex = 1.15) 
  } else {
    # Names pointing outward for small sectors to avoid overlapping
    circos.text(x = CELL_META$xcenter, 
                y = CELL_META$ylim[1] + mm_y(6), 
                labels = CELL_META$sector.index, 
                facing = "clockwise", 
                niceFacing = TRUE, 
                adj = c(0, 0.5), 
                cex = 1.0) 
  }
}, bg.border = NA)

# TRACK 2: Axes and numbers (Smart scaling based on quantity)
circos.track(track.index = 2, panel.fun = function(x, y) {
  max_val <- round(CELL_META$xlim[2]) # Extract the total flow for this genus
  
  if(max_val < 20) {
    # If there are fewer than 20, explicitly force ticks for 0 and the max value
    circos.axis(h = "top", 
                labels.cex = 0.6, 
                labels.niceFacing = TRUE,
                major.at = c(0, max_val),  # Where to put the ticks
                labels = c(0, max_val),    # What numbers to write
                major.tick.length = mm_y(1.5), 
                sector.index = CELL_META$sector.index, 
                track.index = 2)
  } else {
    # For larger sectors, calculate automatically (e.g., 20, 40, 60...)
    circos.axis(h = "top", 
                labels.cex = 0.6, 
                labels.niceFacing = TRUE,
                major.tick.length = mm_y(1.5), 
                sector.index = CELL_META$sector.index, 
                track.index = 2)
  }
}, bg.border = NA)

dev.off()
circos.clear()