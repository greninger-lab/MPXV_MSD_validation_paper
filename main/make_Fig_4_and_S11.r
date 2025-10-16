#Fig 4A neut data - categorical scatter
#Fig 4B - heat map of correlations - Spearman correlation between antigens and neutralization

#Figure S11 - plots of data with correlations

# Greninger Lab - Jon Reed
# 2025-Oct-15

source("MPXV_paper_functions.r")

package_list <- c("tidyverse",
                  "ggplot2",
                  "rstatix",
                  "scales",
                  "ggbeeswarm",
                  "ggnewscale",
                  "ggtext",
                  "patchwork",
                  "cowplot",
                  "ggh4x")

invisible(lapply(package_list, check_and_load))

# load data
dat_MSD <- read_csv("2025-Oct-15_compiled_data.csv")

#Tibble for plotting groups names and p-values above the plot
groups <- tibble(
  group = c(
    "pediatric", "rubella IgG+",
    "Vaccine A - 0 month", "Vaccine A - 1 month", "Vaccine A - 2 month",
    "Vaccine B - 9 month", "MPXV - 2-6 month", "MPXV - 10 month"),
  p_labels = c(
    "pediatric" = "Pediatric",
    "rubella IgG+" = "Rubella~IgG^'+'",
    "Vaccine A - 0 month" = "Vaccine~A~'-'~0~mo",
    "Vaccine A - 1 month" = "Vaccine~A~'-'~1~mo",
    "Vaccine A - 2 month" = "Vaccine~A~'-'~2~mo",
    "Vaccine B - 9 month" = "Vaccine~B~'-'~9~mo",
    "MPXV - 2-6 month" = "MPXV~'-'~2-6~mo",
    "MPXV - 10 month" = "MPXV~'-'~10~mo"),
  y = 2^seq(from = log2(1400), by = 0.4, length.out = 8),
  group_num = seq(1,8,1),
)

#cohort labels for x-axis, same order as groups
cohort_labels = c("Pediatric", expression(Rubella~IgG^"+"), "0 mo", "1 mo", "2 mo", "9 mo", "2-6 mo", "10 mo")

# prep neut results
dat_neut <- dat_MSD %>%
  select(Sample_ID, group, ND50_mean, vaccine_history) %>%
  distinct() %>%
  #Add group_num for plotting on continuous axis
  left_join(groups %>% select(group, group_num), by = "group")


#### Create all unique pairwise combination, do t-test ####

# combos
combn_mat <- t(combn(groups$group, 2))
colnames(combn_mat) <- c("group1", "group2")

# t-test
dat_t_test <- as_tibble(combn_mat) %>%
  # for each combo, pull relavent data from dat_neut
  # Assign an index for each pairwise comparison, for latter grouping
  mutate(
    data = map2(group1, group2, ~dat_neut %>% filter(group %in% c(.x,.y))),
    index = seq(1, length(group1),1)
  ) %>%
  unnest(data) %>%
  # Log transform ND50 before analysis
  mutate(
    ND50_mean_log = log(ND50_mean)
  ) %>%
  # By index to group each pairwaise comparison
  group_by(index) %>%
  # Do t-test
  t_test(ND50_mean_log ~ group, var.equal = FALSE, alternative = "two.sided", detailed = TRUE) %>%
  # Adjust for multiple comparisons 
  adjust_pvalue(method = "bonferroni") %>%
  # Create joining column group1_group2
  mutate(
    join_col = map2_chr(group1, group2, ~paste0(sort(c(.x, .y)), collapse = "_"))
  ) %>%
  # delog ND50 mean estimates
  mutate(
    estimate1_delog = exp(estimate1),
    estimate2_delog = exp(estimate2),
    .after = estimate2
  ) 

#Create a data frame of pairwise comparisons for plotting purposes

pairwise_df <- as_tibble(combn_mat) %>%
  # create join columns
  mutate(
    join_col = map2_chr(group1, group2, ~paste0(sort(c(.x, .y)), collapse = "_"))
  ) %>%
  # join with dat_t_test to get p.adj
  left_join(dat_t_test %>% select(join_col, p.adj), by = "join_col") %>%
  # cap lowest p-value to 1e-7 for plotting
  mutate(
    p.adj_plot = if_else(p.adj < 1e-7, 1e-7, p.adj)
  ) %>%
  # add group number for plotting along x-axis and y-position
  left_join(groups %>% select(group, group_num), by = c("group1" = "group")) %>%
  left_join(groups %>% select(group, y), by = c("group2" = "group"))


pairwise_df$group1 <- factor(pairwise_df$group1,
                             levels = groups$group)

pairwise_df$group2 <- factor(pairwise_df$group2,
                             levels = groups$group)


#=====================================================

# order the prev. vaccine history levels
dat_neut$vaccine_history <- factor(dat_neut$vaccine_history, levels = c("yes", "no"))

# tibble for plotting x-axis annotations
yseg <- -0.9  #set y-position for grouping segment
ytext <- -1.0 #set y-position for x-axis labels

x_axis_annotation <- tribble(
  ~group, ~text, ~x_text, ~y_text, ~x, ~y, ~xend, ~yend,
  "negative", "Negative Cohort", 1.5, ytext, 0.5, yseg, 2.4, yseg,
  "vac_cohort_A", "Vaccine Cohort A", 4, ytext, 2.6, yseg, 5.4, yseg,
  "vac_cohort_B", "Vaccine Cohort B", 6,  ytext, 5.6, yseg, 6.4, yseg,
  "MPXV", "MPXV-infected Cohorts", 7.5, ytext, 6.6, yseg, 8.4, yseg
)

gmean_plot <- dat_neut %>%
  group_by(group) %>%
  summarize(
    ND50_gmean = log2(exp(mean(log(ND50_mean))))
  ) %>%
  left_join(groups %>% select(group, group_num), by = "group") %>%
  mutate(
    x = group_num - 0.4,
    xend = group_num + 0.4,
    y = ND50_gmean,
    yend = ND50_gmean
  )

#picking and viewing colors
my_pal <- pal_brewer(palette = "YlOrRd")(9)

# Set colors for prev. vaccination key
prior_vac_colors = c("yes" = "red", "no" = "black")
dat_neut$vaccine_history <- factor(dat_neut$vaccine_history, levels = c("yes", "no"))

# tibble for plotting self-pairwise diagonal in black
diag <- tibble(
  x = seq(1,8,1),
  y = groups$y
)


Fig_4A <- ggplot(dat_neut, aes(x = group_num, y = log2(ND50_mean), color = vaccine_history, group = group_num)) + 
  geom_violin(
    aes(y = log2(ND50_mean)),
    color = "grey42",
    linewidth = 0.1,
    trim = TRUE
  ) +
  geom_quasirandom(
    shape = 21,
    varwidth = TRUE,
    size = 1,
    width = 0.2,
    show.legend = FALSE
  ) +
  coord_cartesian(
    clip = "off",
    ylim = c(0,10),
    xlim = c(0.5,8.5)
  ) +
  scale_color_manual(
    values = prior_vac_colors,
    name = "Smallpox vaccination history: ",
    labels = ~str_wrap(.x,10),
    guide = guide_legend(order = 3)
  ) +
  scale_y_continuous(
    labels = ~comma(2^.x),
    limits = c(-20,15),
    breaks = seq(0,10,2),
    name = "Neutralizing Antibody Level (ND50)",
    expand = expansion(mult = c(0.025,0))
  ) +
  scale_x_continuous(
    labels = cohort_labels,
    breaks = seq(1,8,1),
    limits = c(-4,9)
  ) +
  #x-axis annotation for negative cohort
  geom_segment(
    inherit.aes = FALSE,
    data = x_axis_annotation %>% filter(group == "negative"),
    aes(x = x, y = y, xend = xend, yend = yend),
    linewidth = 0.25
  ) +
  geom_text(
    inherit.aes = FALSE,
    data = x_axis_annotation %>% filter(group == "negative"),
    aes(x = x_text, y = y_text, label = text),
    vjust = 1,
    size = 6/.pt
  ) +
  # #x-axis annotation for vaccine cohort A
  geom_segment(
    inherit.aes = FALSE,
    data = x_axis_annotation %>% filter(group == "vac_cohort_A"),
    aes(x = x, y = y, xend = xend, yend = yend),
    linewidth = 0.25
  ) +
  geom_text(
    inherit.aes = FALSE,
    data = x_axis_annotation %>% filter(group == "vac_cohort_A"),
    aes(x = x_text, y = y_text, label = text),
    vjust = 1,
    size = 6/.pt
  ) +
  #x-axis annotation for vaccine cohort B
  geom_segment(
    inherit.aes = FALSE,
    data = x_axis_annotation %>% filter(group == "vac_cohort_B"),
    aes(x = x, y = y, xend = xend, yend = yend),
    linewidth = 0.25
  ) +
  geom_text(
    inherit.aes = FALSE,
    data = x_axis_annotation %>% filter(group == "vac_cohort_B"),
    aes(x = x_text, y = y_text, label = str_wrap(text,10)),
    vjust = 1,
    size = 6/.pt
  ) +
  # #x-axis annotation for MPXV infected
  geom_segment(
    inherit.aes = FALSE,
    data = x_axis_annotation %>% filter(group == "MPXV"),
    aes(x = x, y = y, xend = xend, yend = yend),
    linewidth = 0.25
  ) +
  geom_text(
    inherit.aes = FALSE,
    data = x_axis_annotation %>% filter(group == "MPXV"),
    aes(x = x_text, y = y_text, label = str_wrap(text,10)),
    vjust = 1,
    size = 6/.pt
  ) +
  #plotting median lines
  geom_segment(
    data = gmean_plot,
    aes(x = x, xend = xend, y = y, yend = yend),
    inherit.aes = FALSE,
    color = "grey20",
    linewidth = 0.3
  ) +
  #plot heat map
  geom_tile(
    data = pairwise_df %>% filter(p.adj_plot <= 0.05),
    inherit.aes = FALSE,
    aes(x = group_num, y = log2(y), fill = -log10(p.adj_plot)),
    color = "black",
    show.legend = FALSE
  ) +
  #add black tiles representing pairwise comparisons to self
  geom_tile(
    data = diag,
    inherit.aes = FALSE,
    aes(x = x, y = log2(y)),
    fill = "black",
    color = "black"
  ) +
  geom_text(
    inherit.aes = FALSE,
    data = groups,
    aes(label = p_labels, y = log2(y), x = 0.4),
    parse = TRUE,
    hjust = 1,
    size = 6/.pt
  ) +
  scale_fill_gradient(
    low = my_pal[1],
    high = my_pal[7],
    na.value = "green",   #There should be no NA values, highlight with green
    name = "p-value",
    breaks = -log10(c(1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7)),
    labels = c(expression(10^-2), expression(10^-3), expression(10^-4), expression(10^-5), expression(10^-6), expression(paste("<",10^-7))),
    guide = guide_colorbar(order = 1)
  ) +
  new_scale_fill() +
  scale_fill_manual(
    values = "grey",
    breaks = "0",
    name = NULL,
    labels = ">0.05",
    guide = guide_legend(order = 2)
  ) +
  geom_tile(
    data = pairwise_df %>% filter(p.adj_plot > 0.05),
    inherit.aes = FALSE,
    aes(x = group_num, y = log2(y), fill = factor(0)),
    color = "black",
    show.legend = FALSE
  ) +
  theme(
    panel.background = element_blank(),
    panel.grid = element_blank(),
    axis.line = element_line(),
    axis.text.x = element_text(size = 6, vjust = 0.5, margin = margin(t = -1)),
    axis.text.y = element_text(size = 6),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 8),
    plot.margin = margin(t = 77, l = 10, r = -20, b = 25, unit = "pt"),
  )


##### Figure 4B ####

list <- ls()
rm(list = list[!list %in% c("gcv", "gmean", "gsd", "check_and_load", "Fig_4A")])

dat_neut_MSD <- read_csv("2025-Oct-15_compiled_data.csv")

dat_neut_MSD_corr <- dat_neut_MSD %>%
  group_by(group, Assay) %>%
  mutate(
    ND50_rank = rank(ND50_mean),
    MSD_rank = rank(MSD_mean)
  ) %>%
  ungroup() %>%
  nest(data = !c(group, Assay)) %>%
  mutate(
    #Spearman correlation, pull estiamte and p-value
    correlation_spe = map(data, ~cor.test(log(.x$ND50_mean), log(.x$MSD_mean), method = "spearman")),
    estimate_spe = map_dbl(correlation_spe, ~.x$estimate),
    p_value_spe = map_dbl(correlation_spe, ~.x$p.value),
    
    #For plotting p-values
    log10_p_value_spe = -log10(p_value_spe),
    plot_p_value_spe = case_when(
      log10_p_value_spe >= 7 ~ 7, 
      log10_p_value_spe <= -log10(0.05) ~ 0,
      TRUE ~ log10_p_value_spe
    )
  )

message("Cannot compute exact p-values with ties, ok")


##### Setup plot for Fig 4B

#Colors p-value plotting
colors <- brewer_pal(palette = "YlOrRd")(9)
#show_col(colors)

#Order prev. vaccine history levels
dat_neut_MSD$vaccine_history <- factor(dat_neut_MSD$vaccine_history, levels = c("yes", "no")) 

#Setup expression labels for groups
group_labels <- c(
  "MPXV - 10 month" = "MPXV~'-'~10~mo",
  "MPXV - 2-6 month" = "MPXV~'-'~2-6~mo",
  "Vaccine B - 9 month" = "Vaccine~B~'-'~9~mo",
  "Vaccine A - 2 month" = "Vaccine~A~'-'~2~mo",
  "Vaccine A - 1 month" = "Vaccine~A~'-'~1~mo",
  "Vaccine A - 0 month" = "Vaccine~A~'-'~0~mo",
  "rubella IgG+" = "Rubella~IgG^'+'",
  "pediatric" = "Pediatric")

#Set the order of group levels
group_levels <- c("MPXV - 10 month", 
                  "MPXV - 2-6 month", 
                  "Vaccine B - 9 month",
                  "Vaccine A - 2 month", 
                  "Vaccine A - 1 month", 
                  "Vaccine A - 0 month",
                  "rubella IgG+", 
                  "pediatric")

#Assign group expressions and group levels to group_plot
dat_neut_MSD_corr$group_plot <- factor(dat_neut_MSD_corr$group, 
                                       levels = group_levels
)

#Setup markdown labels for x axis
x_labels = c(
  "<span style = 'color:#1F78B4;'>MPXV A29L</span>",
  "<span style = 'color:#B15928;'>VACV A27L</span>",
  "<span style = 'color:#1F78B4;'>MPXV A35R</span>",
  "<span style = 'color:#B15928;'>VACV A33R</span>",
  "<span style = 'color:#1F78B4;'>MPXV B6R</span>",
  "<span style = 'color:#B15928;'>VACV B5R</span>",
  "<span style = 'color:#1F78B4;'>MPXV E8L</span>",
  "<span style = 'color:#B15928;'>VACV D8L</span>",
  "<span style = 'color:#1F78B4;'>MPXV M1R</span>",
  "<span style = 'color:#B15928;'>VACV L1R</span>"
)

#Define order of x labels
x_levels = c("MPXV A29L","VACV A27L",
             "MPXV A35R", "VACV A33R",
             "MPXV B6R", "VACV B5R",
             "MPXV E8L", "VACV D8L",
             "MPXV M1R", "VACV L1R")


#Assign x levels and labels to Assay_plot
dat_neut_MSD_corr$Assay_plot <- factor(dat_neut_MSD_corr$Assay, levels = x_levels,
                                       labels = x_labels)


##### Plot 4B

Fig_4B <- ggplot(dat_neut_MSD_corr, aes(x = Assay_plot, y = group_plot, fill = plot_p_value_spe)) + 
  geom_tile(
    color = NA,
    fill = NA
  ) +
  geom_tile(
    data = dat_neut_MSD_corr %>% filter(p_value_spe <= 0.05),
    color = "black",
    show.legend = FALSE
  ) +
  scale_fill_gradient(
    low = colors[1],
    high = colors[7],
    name = "p-value",
    breaks = c(2, 3, 4, 5, 6, 7),
    labels = c(expression(10^-2), expression(10^-3), expression(10^-4), expression(10^-5), expression(10^-6), expression(paste("<", 10^-7))),
    limits = c(-log10(0.05),7),
    guide = guide_colorbar(order = 1)
  ) +
  new_scale_fill() +
  scale_fill_manual(
    values = "grey",
    breaks = "0",
    name = NULL,
    labels = ">0.05",
    guide = guide_legend(order = 2)
  ) +
  geom_tile(
    data = dat_neut_MSD_corr %>% filter(p_value_spe > 0.05),
    aes(fill = factor(0)),
    color = "black",
    show.legend = FALSE
  ) +
  geom_text(
    aes(label = sprintf("%0.2f", estimate_spe)),
    fontface = "bold",
    color = "black",
    size = 5/.pt,
    hjust = 0.5
  ) +
  scale_x_discrete(
    labels = x_labels,
    name = "Assay",
    expand = c(0,0)
  ) +
  scale_y_discrete(
    name = "Cohort",
    expand = c(0,0),
    labels = label_parse()(group_labels),
  ) +
  theme(
    panel.background = element_blank(),
    panel.grid = element_blank(),
    axis.text.y = element_text(size = 6, vjust = 0.4),
    axis.text.x = element_markdown(size = 6, angle = 90, vjust = 0.5, hjust = 1, margin = margin(l = 1, t = 2, unit = "pt")),
    axis.title = element_text(size = 8),
    plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"),
    legend.text = element_text(size = 6),
    legend.key.size = unit(10, "pt"),
    legend.title = element_text(size = 6)
  )

############# Build Figure 4 #####################


########## Create Legends ####

test <- tibble(
  x = 1,
  y = 1,
  vaccine_history = rep(c("yes","no"),4),
  p_value = -log10(c(0.1, 0.05, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7))
)

# Create the plot
legend_plot_1 <- ggplot(test, aes(x = x, y = y, fill = p_value)) +
  geom_tile() +  # triggers the legend
  scale_fill_gradient(
    low = colors[1],
    high = colors[7],
    name = "p-value",
    breaks = c(2, 3, 4, 5, 6, 7),
    labels = c(expression(10^-2), expression(10^-3), expression(10^-4), expression(10^-5), expression(10^-6), expression(paste("<", 10^-7))),
    limits = c(-log10(0.05),7),
    guide = guide_colorbar(order = 1),
  ) +
  new_scale_fill() +
  scale_fill_manual(
    values = "grey",
    breaks = "0",
    name = NULL,
    labels = ">0.05",
    guide = guide_legend(order = 2)
  ) +
  geom_tile(aes(fill = factor(0))) +
  theme_void() +
  theme(
    legend.title = element_text(size = 7),
    legend.text = element_text(size = 6),
    legend.key.size = unit(12, "pt"),
    legend.position = "right",
    legend.box.margin = margin(t = 0, r = 20, b = 10, l = 0),
    legend.margin = margin(t = 0, r = 0, b = 0, l = 0)
  ) 


#Assign color for prev. vaccine history
prior_vac_colors <- c("yes" = "red", "no" = "black")
test$vaccine_history <- factor(test$vaccine_history, levels = c("yes", "no"))


legend_plot_2 <- ggplot(test, aes(x = x, y = y, color = vaccine_history)) +
  geom_point(
    size = 1.5,
    shape = 21,
    inherit.aes = FALSE,
    aes(x = x, y = y, color = vaccine_history),
    show.legend = TRUE
  ) +
  scale_color_manual(
    values = prior_vac_colors,
    name = "Smallpox\nvaccination\nhistory:",
    drop = FALSE
  ) +
  theme_void() +
  theme(
    legend.title = element_text(size = 7),
    legend.text = element_text(size = 7),
    legend.key.size = unit(12, "pt"),
    legend.position = "right",
    legend.box.margin = margin(t = 0, r = 50, b = 45, l = 0),
    legend.margin = margin(t = 0, r = 0, b = 0, l = 0)
  ) 



# Extract just the legend using cowplot
legend_only_1 <- cowplot::get_legend(legend_plot_1)
legend_only_2 <- cowplot::get_legend(legend_plot_2)

#####

layout <- c(
  area(t = 1, l = 1, b = 5, r = 7),
  area(t = 1, l = 8, b = 3, r = 10),
  area(t = 4, l = 8, b = 5, r = 10)
)



Fig_4 <- wrap_elements(Fig_4A) + Fig_4B + (wrap_elements(legend_only_2) + wrap_elements(legend_only_1)) +
  plot_layout(design = layout) +
  plot_annotation(tag_levels = list(c("A","B","","")))


ggsave("fig_output/Fig_4.pdf", width = 8, height = 4.5, units = "in", dpi = 600, device = cairo_pdf)


################################ Fig S11 ##################################

##### Setup plot for Fig S11

#reverse the levels and labels
group_labels <- rev(group_labels)
group_levels <- rev(group_levels)

#reassign group_plot factor reversed levels
dat_neut_MSD_corr$group_plot <- factor(dat_neut_MSD_corr$group,
                                       levels = group_levels)


#create a tibble for just the facet labels.  
facet_labels <- tibble(
  Assay = c("MPXV A29L","VACV A27L",
            "MPXV A35R", "VACV A33R",
            "MPXV B6R", "VACV B5R",
            "MPXV E8L", "VACV D8L",
            "MPXV M1R", "VACV L1R"),
  label = c(
    "<span style = 'color:#1F78B4;'>MPXV A29L</span>",
    "<span style = 'color:#B15928;'>VACV A27L</span>",
    "<span style = 'color:#1F78B4;'>MPXV A35R</span>",
    "<span style = 'color:#B15928;'>VACV A33R</span>",
    "<span style = 'color:#1F78B4;'>MPXV B6R</span>",
    "<span style = 'color:#B15928;'>VACV B5R</span>",
    "<span style = 'color:#1F78B4;'>MPXV E8L</span>",
    "<span style = 'color:#B15928;'>VACV D8L</span>",
    "<span style = 'color:#1F78B4;'>MPXV M1R</span>",
    "<span style = 'color:#B15928;'>VACV L1R</span>"
  )
)

label_map <- setNames(facet_labels$label, facet_labels$Assay)


dat_neut_MSD_corr_unnested <- dat_neut_MSD_corr %>%
  unnest(data)

dat_neut_MSD_corr_unnested$vaccine_history <- factor(dat_neut_MSD_corr_unnested$vaccine_history, levels = c("yes", "no"))

Fig_S11 <- ggplot(dat_neut_MSD_corr_unnested, aes(x = ND50_rank, y = MSD_rank, color = vaccine_history)) +
  facet_grid2(
    group_plot ~ Assay_plot,
    axes = "all",
    switch = "y",
    scales = "free",
    independent = "x",
    labeller = labeller(
      group_plot = as_labeller(group_labels, default = label_parsed), 
      Assay = label_map
    )  
  ) +
  geom_point(
    size = 1
  ) +
  scale_color_manual(
    values = prior_vac_colors,
    name = "Smallpox\nvaccination\nhistory",
    label = ~str_wrap(.x,10)
  ) +
  scale_y_continuous(
    name = "rank(MSD result)",
    breaks = pretty_breaks(n = 4),
  ) +
  scale_x_continuous(
    name = "rank(ND50)",
    breaks = pretty_breaks(n = 4)
  ) +
  geom_label(
    data = dat_neut_MSD_corr %>% select(group_plot, Assay_plot, estimate_spe) %>% unique() 
    %>% mutate(estimate_spe = sprintf("%0.2f", estimate_spe)),
    aes(label = paste0("rho==\"",estimate_spe,"\"")),
    y = Inf,
    x = -Inf,
    hjust = "left",
    vjust = 1,
    parse = TRUE,
    size = 10/.pt,
    color = "blue",
    fill = "grey",
    label.padding = unit(0, "pt"),
    alpha = 0.75
  ) +
  geom_abline(
    intercept = 0,
    slope = 1,
    linetype = "dashed",
    linewidth = 0.5,
    color = "grey"
  ) +
  theme(
    panel.background = element_blank(),
    axis.line = element_line(color = "black"),
    aspect.ratio = 1,
    strip.text.x = element_markdown(size = 8),
    strip.text.y = element_text(size = 8),
    strip.placement = "outside"
  )

ggsave("fig_output/Fig_S11.pdf", width = 15, height = 15, unit = "in", dpi = 600, device = cairo_pdf)

