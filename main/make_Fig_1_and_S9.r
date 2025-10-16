# Fig 1 - 
# A - similarity between MPXV and VACV antigens
# B - virion diagram showing location
# C - Correlation heatmap between MPXV and VACV orthologs

# Fig S9 -
# Correlation plot between MPXV and VACV orthologs

# Greninger Lab - Jon Reed
# 2025-Oct-15


source("MPXV_paper_functions.r")

package_list <- c("tidyverse",
                  "ggplot2",
                  "scales",
                  "openxlsx2",
                  "ggtext",
                  "ggnewscale",
                  "patchwork")

          
invisible(lapply(package_list, check_and_load))

######################### Fig 1A Similarity plots ##############################

# Load Biostrings for this section, and unload at the end
check_and_load("Biostrings")

calc_similarity <- function(window_matrix) {
  
  #for each column in the window first tally frequency of each amino acid using table
  #get count of the most frequent a.a. with max
  #calculate %similarity by count most freq. a.a./length of column
  #finally take the mean of these value to get %similarity of the window
  
  colnames(window_matrix) <- paste0("V", 1:ncol(window_matrix))
  
  #gaps are counted as a character
  x <- as_tibble(window_matrix) %>%
    summarise(
      across(everything(), ~if_else(length(unique(.x)) == 1, 1, 0))
    ) %>%
    rowMeans()
  
}

## This function will shrink window near the begining and end
set_window <- function(position, length, window_size) {
  
  if(window_size %% 2 == 0) {
    window_l = position - window_size/2
    window_r = position + ((window_size/2) - 1)
  } else {
    window_l = position - floor(window_size/2)
    window_r = position + floor(window_size/2)
  }
  
  if(window_l < 1) {
    window_l <- 1
  }
  
  if(window_r > length) {
    window_r <- length
  }
  
  return(list(left = window_l, right = window_r))
  
}

## This function returns windows size based on fixed window length
## This will return window positions that are negative or greater
## than the length of the alignment, which need to be filtered out
set_window_fixed <- function(position, length, window_size) {
  
  if(window_size %% 2 == 0) {
    window_l = position - window_size/2
    window_r = position + ((window_size/2) - 1)
  } else {
    window_l = position - floor(window_size/2)
    window_r = position + floor(window_size/2)
  }
  
  return(list(left = window_l, right = window_r))
  
}

# Load alignments
alignment <- 
  tibble(fasta = list.files("similarity files/MSD antigen pairwise", full.names = TRUE)) %>%
  mutate(
    gene = str_extract(fasta, "VACV.+[^.fasta]"),
    alignment = map(fasta, ~readAAStringSet(.x)),
    alignment_matrix = map(alignment, ~as.matrix(.x)),
    window_size = 15,
    alignment_length = map_int(alignment, ~width(.x)[1]),
    num_seqs = map_int(alignment, ~length(.x)),
    similarity = map(alignment_length, ~tibble(position = 1:.x)),
  ) %>%
  unnest(similarity) %>%
  group_by(fasta) %>%
  mutate(
    window_pos = pmap_dfr(list(position, alignment_length, window_size), function(x,y,z) {set_window_fixed(x, y, z)})
  ) %>% 
  unnest_wider(window_pos) %>%
  group_by(gene) %>%
  filter(left > 0 & right <= alignment_length) %>%
  mutate(
    window_width = (right - left) + 1,
    window = pmap(list(alignment, left, right), function(alignment, left, right) {subseq(alignment, start = left, end = right)}),
    window_matrix = map(window, ~as.matrix(.x)),
    similarity = map_dbl(window_matrix, ~calc_similarity(.x))
  )

dat <- read_xlsx("similarity files/MSD_antigen_domain.xlsx") %>%
  mutate(
    level = case_when(
      str_detect(antigen, "VACV") ~ 0,
      str_detect(antigen, "MPXV") ~ 0.4
    )
  )

make_topo_plot <- function(ortho_pair) {
  
  #Luminal domain – Light Blue #A6CEE3
  #Transmembrane domain – Dark Gray #636363
  #Cytosolic domain – Light Green #B2DF8A
  
  labels <- tribble(
    ~virus, ~antigen, ~level,
    "VACV", str_extract(ortho_pair, "VACV_[^_]+") %>% str_replace("_"," "), 0 + 0.15,
    "MPXV", str_extract(ortho_pair, "MPXV_[^_]+") %>% str_replace("_"," "), 0.4 + 0.15
  )
  
  domain_color <- c("CT" = "#B2DF8A", "TM" = "lightpink", "ED" = "#A6CEE3")
  
  plot <- ggplot(dat %>% filter(ortholog_pair_name == ortho_pair), aes(xmin = start, xmax = end, ymin = level, ymax = level + 0.3, fill = domain_ab)) +
    geom_rect(
      show.legend = FALSE
    ) +
    geom_text(
      aes(x = (end-start)/2 + start, y = level + 0.15, label = domain_ab),
      color = "black",
      size = 4/.pt
    ) +
    geom_text(
      data = labels,
      inherit.aes = FALSE,
      aes(x = 0, y = level, label = antigen, color = virus),
      hjust = 1,
      size = 5/.pt,
      show.legend = FALSE
    ) +
    scale_color_manual(
      values = c("#1F78B4","#B15928")
    ) +
    scale_fill_manual(
      values = domain_color
    ) +
    coord_cartesian(
      clip = "off",
      xlim = c(1, max(alignment$alignment_length))
    ) +
    scale_x_continuous(
      breaks = c(1, seq(40, max(alignment$alignment_length), 40)),
      limits = c(-20, max(alignment$alignment_length)),
      expand = c(0,0)
    ) +
    labs(
      x = "Alignment Position",
      y = "Amino Acid Pairwise Similarity (%)"
    ) +
    theme(
      plot.margin = margin(t = 0, r = 0, b = 0, l = 30, unit = "pt"),
      panel.background = element_blank(),
      panel.grid = element_blank(),
      axis.line.y = element_blank(),
      axis.line = element_line(color = "black"),
      axis.ticks.y = element_blank(),
      axis.text.y = element_blank(),
      axis.text.x = element_text(size = 5),
      axis.title.x = element_text(size = 8),
      axis.title.y = element_text(size = 8, margin = margin(r = 20, unit = "pt"))
      
    )
  
  plot
}

topo_plots <- tibble(
  pair = dat$ortholog_pair_name %>% unique(),
  graph = map(pair, make_topo_plot)
)

make_sim_plot <- function(ortho_pair, overall_sim) {
  
  ggplot(alignment %>% filter(gene == ortho_pair), aes(x = position, y = similarity)) + 
    scale_y_continuous(
      limits = c(0.65,1),
      breaks = seq(0.70,1, 0.1),
      label = percent
    ) +
    scale_x_continuous(
      breaks = c(1, seq(40, max(alignment$alignment_length), 40)),
      limits = c(1, max(alignment$alignment_length)),
      expand = c(0,0)
    ) +
    geom_line(
      linewidth = 0.25
    ) +
    labs(
      x = "Alignment Position",
      y = "Amino Acid Pairwise Similarity (%)"
    ) +
    annotate(
      geom = "text",
      x = 5,
      y = 0.7,
      hjust = 0,
      size = 5/.pt,
      label = paste0("Overall Similarity: ", overall_sim)
    ) +
    theme(
      panel.grid = element_blank(),
      panel.background = element_blank(),
      axis.line = element_line(color = "black"),
      axis.text = element_text(size = 5),
      axis.title.x = element_text(size = 8),
      axis.title.y = element_text(size = 8, margin = margin(r = 20, unit = "pt"))
    )
}

sim_plots <- tibble(
  pair = alignment$gene %>% unique()
)

overall_sim <- tribble(~pair, ~over_sim,
                       "VACV_A27L_MPXV_A29L", "94.5%",
                       "VACV_A33R_MPXV_A35R", "95.6%",
                       "VACV_B5R_MPXV_B6R", "96.5%",
                       "VACV_D8L_MPXV_E8L", "94.4%",
                       "VACV_L1R_MPXV_M1R", "98.8%"
)

sim_plots <- sim_plots %>% full_join(overall_sim, by = "pair")

sim_plots <- sim_plots %>%
  mutate(
    graph = map2(pair, over_sim, make_sim_plot)
  )


# remove these libraries due to conflicts with tidyverse
detach("package:Biostrings", unload = TRUE)
detach("package:GenomeInfoDb", unload = TRUE)
detach("package:XVector", unload = TRUE)
detach("package:IRanges", unload = TRUE)
unloadNamespace("UCSC.utils")
detach("package:S4Vectors", unload = TRUE)


############################# Fig 1C #########################################

dat_MSD <- read_csv("2025-Oct-15_compiled_data.csv")

dat_MSD_wide <- dat_MSD %>%
  select(Sample_ID, Assay, MSD_mean, group) %>%
  mutate(
    antigen_group = case_when(
      Assay %in% c("MPXV A29L", "VACV A27L") ~ "MPXV A29L/VACV A27L",
      Assay %in% c("MPXV A35R", "VACV A33R") ~ "MPXV A35R/VACV A33R",
      Assay %in% c("MPXV B6R", "VACV B5R") ~ "MPXV B6R/VACV B5R",
      Assay %in% c("MPXV E8L", "VACV D8L") ~ "MPXV E8L/VACV D8L",
      Assay %in% c("MPXV M1R", "VACV L1R") ~ "MPXV M1R/VACV L1R",
      TRUE ~ NA_character_
    ),
    axis = if_else(str_detect(Assay, "MPXV"), "MPXV", "VACV")
  ) %>%
  group_by(antigen_group) %>%
  pivot_wider(
    id_cols = c(Sample_ID, group, antigen_group),
    names_from = axis,
    values_from = MSD_mean
  ) %>%
  nest(data = c(Sample_ID, MPXV, VACV)) %>%
  mutate(
    correlation = map(data, ~cor.test(log2(.x$MPXV),
                                      log2(.x$VACV),
                                      conf.level = 0.95,
                                      method = "pearson",
    )),
    estimate = map_dbl(correlation, ~.x$estimate),
    p_value = map_dbl(correlation, ~.x$p.value),
    lower_est = map_dbl(correlation, ~.x$conf.int[1]),
    upper_est = map_dbl(correlation, ~.x$conf.int[2]),
    log10_p_value = -log10(p_value),
    plot_p_value = case_when(
      log10_p_value >= 20 ~ 20,
      log10_p_value < 1 ~ 0,
      TRUE ~ log10_p_value
    ),
    x = str_extract(antigen_group, "MPXV.+(?=\\/)"),
    y = str_extract(antigen_group, "VACV.+")
  )

#create a tibble for just the facet labels.  
x_labels = c(
  "<span style = 'color:#1F78B4;'>MPXV A29L</span><span style = 'color:#000000;'><br></span><span style = 'color:#B15928;'>VACV A27L</span>",
  "<span style = 'color:#1F78B4;'>MPXV A35R</span><span style = 'color:#000000;'><br></span><span style = 'color:#B15928;'>VACV A33R</span>",
  "<span style = 'color:#1F78B4;'>MPXV B6R</span><span style = 'color:#000000;'><br></span><span style = 'color:#B15928;'>VACV B5R</span>",
  "<span style = 'color:#1F78B4;'>MPXV E8L</span><span style = 'color:#000000;'><br></span><span style = 'color:#B15928;'>VACV D8L</span>",
  "<span style = 'color:#1F78B4;'>MPXV M1R</span><span style = 'color:#000000;'><br></span><span style = 'color:#B15928;'>VACV L1R</span>"
)


colors <- brewer_pal(palette = "YlOrRd")(9)

dat_MSD_wide$group_plot <- factor(dat_MSD_wide$group,
                                  levels =  c("MPXV - 10 month", "MPXV - 2-6 month", "Vaccine B - 9 month",
                                              "Vaccine A - 2 month", "Vaccine A - 1 month", "Vaccine A - 0 month",
                                              "rubella IgG+", "pediatric"),
                                  labels =  c("MPXV - 10 mo", "MPXV - 2-6 mo", "vaccine B - 9 mo",
                                              "vaccine A - 2 mo", "vaccine A - 1 mo", "vaccine A - 0 mo",
                                              "rubella", "pediatric"))


#cohort labels
group_labels <- c("MPXV - 10 mo", "MPXV - 2-6 mo", "Vaccine B - 12 mo",
                  "Vaccine A - 2 mo", "Vaccine A - 1 mo",  "Vaccine A - 0 mo", 
                  expression(Rubella~IgG^"+"),  "Pediatric")


Fig_1C <- ggplot(dat_MSD_wide, aes(x = antigen_group, y = group_plot, fill = plot_p_value)) + 
  geom_tile(
    data = dat_MSD_wide %>% filter(p_value <= 0.1),
    color = "black",
  ) +
  scale_fill_gradient(
    low = colors[1],
    high = colors[7],
    name = "p-value",
    breaks = c(1, 5, 10, 15, 20),
    labels = c(expression(10^-1), expression(10^-5), expression(10^-10), expression(10^-15), expression(paste("<",10^-20))),
    guide = guide_colorbar(order = 1)
  ) +
  new_scale_fill() +
  geom_tile(
    data = dat_MSD_wide %>% filter(p_value > 0.05),
    aes(fill = factor(plot_p_value)),
    color = "black"
  ) +
  scale_fill_manual(
    values = "grey",
    breaks = 0,
    name = NULL,
    label = ">0.05",
    guide = guide_legend(order = 2)
  ) +
  geom_text(
    aes(label = format(estimate,digits = 1)),
    fontface = "bold",
    color = "black",
    size = 6/.pt
  ) +
  scale_x_discrete(
    labels = x_labels,
    name = "Antigen Group",
    expand = c(0,0)
  ) +
  scale_y_discrete(
    name = "Cohort",
    labels = group_labels,
    expand = c(0,0)
  ) +
  theme(
    axis.text.y = element_text(size = 6, vjust = 0.3),
    axis.title = element_text(size = 8),
    axis.text.x = element_markdown(size = 6, angle = 90, vjust = 0.5, hjust = 1),
    legend.text = element_text(size = 6),
    legend.key.size = unit(10, "pt"),
    legend.title = element_text(size = 6),
    legend.spacing.y = unit(3, units = "pt")
  )

list <- ls()
rm(list = list[!list %in% c("gcv", "gmean", "gsd", "check_and_load", "sim_plots", "topo_plots", "Fig_1C", "dat_MSD_wide")])

# Put together Fig 1A, 1B (diagram of virion), 1C

# Load your SVG file
img <- png::readPNG("similarity files/virion image/virion_antigen_location.png")
grob_img <- grid::rasterGrob(img, interpolate = TRUE)
Fig_1B <- wrap_elements(full = grob_img)

layout <- c(
  area(t = 1, l = 1, b = 2, r = 3),
  area(t = 3, l = 1, b = 6, r = 3),
  
  area(t = 7, l = 1, b = 8, r = 3),
  area(t = 9, l = 1, b = 12, r = 3),
  
  area(t = 13, l = 1, b = 14, r = 3),
  area(t = 15, l = 1, b = 18, r = 3),
  
  area(t = 19, l = 1, b = 20, r = 3),
  area(t = 21, l = 1, b = 24, r = 3),
  
  area(t = 25, l = 1, b = 26, r = 3),
  area(t = 27, l = 1, b = 30, r = 3)
)

p1 <- topo_plots[topo_plots$pair == "VACV_A27L_MPXV_A29L",]$graph[[1]]
p2 <- sim_plots[sim_plots$pair == "VACV_A27L_MPXV_A29L",]$graph[[1]]
p3 <- topo_plots[topo_plots$pair == "VACV_A33R_MPXV_A35R",]$graph[[1]]
p4 <- sim_plots[sim_plots$pair == "VACV_A33R_MPXV_A35R",]$graph[[1]]
p5 <- topo_plots[topo_plots$pair == "VACV_B5R_MPXV_B6R",]$graph[[1]]
p6 <- sim_plots[sim_plots$pair == "VACV_B5R_MPXV_B6R",]$graph[[1]]
p7 <- topo_plots[topo_plots$pair == "VACV_D8L_MPXV_E8L",]$graph[[1]]
p8 <- sim_plots[sim_plots$pair == "VACV_D8L_MPXV_E8L",]$graph[[1]]
p9 <- topo_plots[topo_plots$pair == "VACV_L1R_MPXV_M1R",]$graph[[1]]
p10 <- sim_plots[sim_plots$pair == "VACV_L1R_MPXV_M1R",]$graph[[1]]

sub_plot <- p1 + p2 + p3 + p4 + p5 + p6 + p7 + p8 + p9 + p10 + 
  plot_layout(design = layout, axis_titles = "collect", axes = "collect") &
  theme(
    plot.margin = margin(t = 0, b = 0, l = 0, r = 0, unit = "pt")
  )

Fig_1A <- wrap_elements(sub_plot)

layout <- c(
  area(t = 1, l = 1, b = 5, r = 4),
  area(t = 1, l = 5, b = 2, r = 6),
  area(t = 3, l = 5, b = 5, r = 6)
)

Fig_1 <- free(Fig_1A, type = "panel") + Fig_1B + Fig_1C +
  plot_layout(design = layout) +
  plot_annotation(tag_levels = "A")

ggsave("fig_output/Fig_1.pdf", width = 6.5, height = 4)

################################ Fig_S9 ###############################

list <- ls()
rm(list = list[!list %in% c("gcv", "gmean", "gsd", "check_and_load", "dat_MSD_wide")])

group_labels <- c(
  "pediatric" = "Pediatric",
  "rubella IgG+" = "Rubella~IgG^'+'",
  "Vaccine A - 0 month" = "Vaccine~A~'-'~0~mo",
  "Vaccine A - 1 month" = "Vaccine~A~'-'~1~mo",
  "Vaccine A - 2 month" = "Vaccine~A~'-'~2~mo",
  "Vaccine B - 9 month" = "Vaccine~B~'-'~9~mo",
  "MPXV - 2-6 month" = "MPXV~'-'~2-6~mo",
  "MPXV - 10 month" = "MPXV~'-'~10~mo")


dat_MSD_wide$group_plot <- factor(dat_MSD_wide$group,
                                  levels = c("pediatric",  "rubella IgG+",
                                             "Vaccine A - 0 month", "Vaccine A - 1 month", "Vaccine A - 2 month",
                                             "Vaccine B - 9 month", "MPXV - 2-6 month","MPXV - 10 month"))



#create a tible for just the facet labels.  
facet_labels <- tibble(
  antigen_group =  c("MPXV A29L/VACV A27L",
                     "MPXV A35R/VACV A33R",
                     "MPXV B6R/VACV B5R",
                     "MPXV E8L/VACV D8L",
                     "MPXV M1R/VACV L1R"),
  label = c(
    "<span style = 'color:#1F78B4;'>MPXV A29L</span><span style = 'color:#000000;'><br></span><span style = 'color:#B15928;'>VACV A27L</span>",
    "<span style = 'color:#1F78B4;'>MPXV A35R</span><span style = 'color:#000000;'><br></span><span style = 'color:#B15928;'>VACV A33R</span>",
    "<span style = 'color:#1F78B4;'>MPXV B6R</span><span style = 'color:#000000;'><br></span><span style = 'color:#B15928;'>VACV B5R</span>",
    "<span style = 'color:#1F78B4;'>MPXV E8L</span><span style = 'color:#000000;'><br></span><span style = 'color:#B15928;'>VACV D8L</span>",
    "<span style = 'color:#1F78B4;'>MPXV M1R</span><span style = 'color:#000000;'><br></span><span style = 'color:#B15928;'>VACV L1R</span>"
  )
)

label_map <- setNames(facet_labels$label, facet_labels$antigen_group)

#Removing one MSD observation that was undetected 
plot_dat <- dat_MSD_wide %>% unnest(data) %>% filter(!(Sample_ID == "P61" & antigen_group == "MPXV B6R/VACV B5R"))
message("P61, MPXV B6R/VACV B5R, was removed due to undetected")

Fig_S9 <- ggplot(plot_dat, aes(x = log2(MPXV), y = log2(VACV))) +
  facet_grid(
    group_plot ~ antigen_group,
    axes = "all_x",
    switch = "y",
    labeller = labeller(
      antigen_group = label_map,
      group_plot = as_labeller(group_labels, label_parsed)
    )
  ) +
  geom_point(
    shape = 21,
    size = 1
  ) +
  scale_y_continuous(
    limits = c(0, 18),
    name = expression(paste("Antibody Level to VACV Antigen (",log[2]," AU/mL)"))
  ) +
  scale_x_continuous(
    limits = c(0,18),
    name = expression(paste("Antibody Level to MPXV Antigen (", log[2], " AU/mL)"))
  ) +
  geom_abline(
    intercept = 0,
    slope = 1,
    linetype = "dashed",
    linewidth = 0.5,
    color = "grey"
  ) +
  geom_label(
    data = dat_MSD_wide %>% select(antigen_group, group_plot, estimate) %>% unique() 
    %>% mutate(estimate = sprintf("%0.3f", estimate)),
    aes(label = paste0("rho==\"",estimate,"\"")),
    y = 0,
    x = Inf,
    hjust = "right",
    vjust = "bottom",
    parse = TRUE,
    size = 6/.pt,
    color = "blue",
    fill = "grey",
    label.padding = unit(0, "pt"),
    alpha = 0.75
  ) +
  theme(
    panel.background = element_blank(),
    axis.line = element_line(color = "black"),
    aspect.ratio = 1,
    strip.text.x = element_markdown(size = 6),
    strip.text.y = element_text(size = 6)
  )

ggsave("fig_output/Fig_S9.pdf", width = 7, height = 8, unit = "in")

