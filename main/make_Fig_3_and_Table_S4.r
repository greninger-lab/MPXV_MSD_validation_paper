# Fig 3 ratio plots between MPXV and VACV orthologs
# Table S4 summary

# Greninger Lab - Jon Reed
# 2025-Oct-15


source("MPXV_paper_functions.r")

package_list <- c("tidyverse",
                  "ggplot2",
                  "rstatix",
                  "scales",
                  "openxlsx2",
                  "ggbeeswarm",
                  "ggpubr",
                  "ggtext")

invisible(lapply(package_list, check_and_load))

dat_MSD <- read_csv("2025-Oct-15_compiled_data.csv")

dat_MSD_wide <- dat_MSD %>%
  select(Sample_ID, Assay, MSD_mean, group, vaccine_history) %>%
  mutate(
    antigen_group = case_when(
      Assay %in% c("MPXV A29L", "VACV A27L") ~ "MPXV A29L/VACV A27L",
      Assay %in% c("MPXV A35R", "VACV A33R") ~ "MPXV A35R/VACV A33R",
      Assay %in% c("MPXV B6R", "VACV B5R") ~ "MPXV B6R/VACV B5R",
      Assay %in% c("MPXV E8L", "VACV D8L") ~ "MPXV E8L/VACV D8L",
      Assay %in% c("MPXV M1R", "VACV L1R") ~ "MPXV M1R/VACV L1R",
      TRUE ~ NA_character_
    ),
    virus = if_else(str_detect(Assay, "MPXV"), "MPXV", "VACV")
  ) %>%
  group_by(antigen_group) %>%
  pivot_wider(
    id_cols = c(Sample_ID, group, antigen_group, vaccine_history),
    names_from = virus,
    values_from = MSD_mean
  ) %>%
  mutate(
    ratio = MPXV/VACV
  )

cohorts <- factor(c("MPXV - 10 month", "MPXV - 2-6 month", "Vaccine B - 9 month",
                    "Vaccine A - 2 month", "Vaccine A - 1 month", "Vaccine A - 0 month",
                    "rubella IgG+", "pediatric"), ordered = TRUE)


combinations <- tibble(status_1 = cohorts) %>%
  mutate(
    x = map(status_1, ~tibble(status_2 = cohorts))
  ) %>%
  unnest(x) %>%
  mutate(
    combo = paste0(pmax(status_1, status_2), "_", pmin(status_1, status_2))
  ) %>%
  distinct(combo, .keep_all = TRUE) %>%
  select(status_1, status_2)


data <- dat_MSD_wide %>%
  mutate(
    log2_ratio = log2(ratio)
  ) %>%
  nest(data = !c(antigen_group, group))



#Check if data_max is needed
data_max <- tibble(antigen_group = c("MPXV A29L/VACV A27L","MPXV A35R/VACV A33R",
                                     "MPXV B6R/VACV B5R", "MPXV E8L/VACV D8L",
                                     "MPXV M1R/VACV L1R"),
                   combinations = list(combinations)) %>%
  unnest(combinations) %>%
  left_join(data, by = c("antigen_group", "status_1" = "group")) %>%
  left_join(data, by = c("antigen_group", "status_2" = "group"), suffix = c("_1", "_2")) %>%
  mutate(
    max_1 = map_dbl(data_1, ~max(.x$log2_ratio)),
    max_2 = map_dbl(data_2, ~max(.x$log2_ratio)),
    max = pmax(max_1, max_2)
  ) %>%
  rename(
    group1 = status_1,
    group2 = status_2
  ) %>%
  mutate(
    groupID = paste0(pmax(group1,group2), "_", pmin(group1,group2))
  ) %>%
  select(antigen_group, groupID, max)

p_value_set <- tibble(antigen_group = c("MPXV A29L/VACV A27L","MPXV A35R/VACV A33R",
                                        "MPXV B6R/VACV B5R", "MPXV E8L/VACV D8L",
                                        "MPXV M1R/VACV L1R"),
                      combinations = list(combinations)) %>%
  unnest(combinations) %>%
  filter(
    #remove self comparisons
    status_1 != status_2,
    #status_1 %in% c("MPXV - 2-6 month", "MPXV - 10 month") & status_2 %in% c("Vaccine A - 1 month", "Vaccine A - 2 month", "Vaccine B - 9 month")
  ) %>%
  mutate(
    index = seq(1,length(antigen_group),1)
  ) %>%
  pivot_longer(
    cols = c(status_1, status_2),
    names_to = "status",
  ) %>%
  left_join(data, by = c("antigen_group", "value" = "group")) %>%
  unnest(data) %>%
  group_by(antigen_group, index) %>%
  t_test(log2_ratio ~ value, alternative = "two.sided", var.equal = FALSE, detailed = TRUE) %>%
  adjust_pvalue(method = "bonferroni") %>%
  mutate(
    groupID = paste0(pmax(group1,group2), "_", pmin(group1,group2))
  ) %>%
  left_join(data_max, by = c("antigen_group", "groupID")) %>%
  mutate(
    in_figure = case_when(
      (group1 %in% c("MPXV - 2-6 month", "MPXV - 10 month") & group2 %in% c("Vaccine A - 1 month", "Vaccine A - 2 month", "Vaccine B - 9 month")) & p.adj < 0.07 ~ TRUE,
      TRUE ~ FALSE
    ),
    in_figure_label = if_else(in_figure, label_pvalue()(p.adj), "")
  )

####### write supp table with t-test results ####

cohort_key <- tribble(
  ~table_value, ~display_value,
  "MPXV - 10 month", "MPXV-infected - 10 month",
  "Vaccine B - 9 month", "Vaccine cohort B - 9 month",
  "MPXV - 2-6 month", "MPXV-infected - 2-6 month",
  "Vaccine A - 0 month", "Vaccine cohort A - 0 month",
  "Vaccine A - 1 month", "Vaccine cohort A - 1 month",
  "Vaccine A - 2 month", "Vaccine cohort A - 2 month",
  "pediatric", "pediatric",
  "rubella IgG+", "rubella IgG+"
)     

###### colors for gradient fill
colors <- pal_brewer(palette = "YlOrRd")(9)

# Configure table for writing to excel
test_pairs_wb <- p_value_set %>% select(
  antigen_group, 
  group1, 
  group2, 
  n1, n2,
  estimate1, estimate2,
  df, statistic,
  p, p.adj,
  in_figure,
  in_figure_label
) %>%
  mutate(
    group1 = map_chr(group1, ~cohort_key[cohort_key$table_value == .x,]$display_value),
    group2 = map_chr(group2, ~cohort_key[cohort_key$table_value == .x,]$display_value)
  ) %>%
  rename(
    `Antigen Group` = antigen_group,
    `Group 1 - cohort` = group1,
    `Group 2 - cohort` = group2,
    `Group 1 - log2 MPXV/VACV ratio` = estimate1,
    `Group 1 - n` = n1,
    `Group 2 - log2 MPXV/VACV ratio` = estimate2,
    `Group 2 - n` = n2,
    `degrees of freedem` = df,
    `t-statistic` = statistic,
    `P-value` = p,
    `P-value adj. MC` = p.adj,
    `printed in figure` = in_figure,
    `display value in figure` = in_figure_label
  ) 

# Create workbook, add data
wb <- wb_workbook() %>% 
  wb_add_worksheet() %>%
  wb_add_data(x = test_pairs_wb)

wb_add_dxfs_style(
  wb,
  name = "grey_fill",
  bg_fill = wb_color(hex = "#BEBEBE")
)

wb <- wb_add_conditional_formatting(
  wb,
  sheet = "Sheet 1",
  dims = wb_dims(
    rows = 2:(nrow(test_pairs_wb)+1),
    cols = 11
  ),
  style = "grey_fill",
  rule = "K2 > 0.05",
  type = "expression"
)

wb <- wb_add_conditional_formatting(
  wb,
  sheet = "Sheet 1",
  dims = wb_dims(
    rows = 2:(nrow(test_pairs_wb)+1),
    cols = 11
  ),
  style = c(colors[9], colors[1]),
  rule = c(1e-7, 0.049999),
  type = "colorScale"
)

# Reduce the number of digits after decimal to 2, for certain columns
wb <- wb_add_numfmt(
  wb,
  sheet = "Sheet 1",
  dims = wb_dims(
    rows = 2:(nrow(test_pairs_wb)+1),
    cols = c(6:9)
  ),
  numfmt = "0.00"
)

# Configure p-values to scientific
wb <- wb_add_numfmt(
  wb,
  sheet = "Sheet 1",
  dims = wb_dims(
    rows = 2:(nrow(test_pairs_wb)+1),
    cols = 10:11
  ),
  numfmt = "0.00E+00"
)

#Center columns
wb <- wb_add_cell_style(
  wb,
  sheet = "Sheet 1",
  dims = wb_dims(
    rows = 1:(nrow(test_pairs_wb)+1),
    cols = c(4:13)
  ),
  vertical = "center",
  horizontal = "center"
)

#Auto adjust column widths
wb_set_col_widths(
  wb,
  sheet = "Sheet 1",
  cols = c(1:13),
  widths = "auto"
)

#Freeze the top row
wb_freeze_pane(wb, sheet = "Sheet 1", first_active_row = 2, first_row = TRUE)

wb_save(wb, "fig_output/Table_S4.xlsx")

########


#cohort labels
group_labels <- c(
  "pediatric", expression(rubella~IgG^"+"),
  "Vaccine A - 0 mo", "Vaccine A - 1 mo", "Vaccine A - 2 mo",
  "Vaccine B - 9 mo", "MPXV - 2-6 mo", "MPXV - 10 mo")


p_value_set$group1 <- factor(p_value_set$group1,
                             levels = c("pediatric", "rubella IgG+",
                                        "Vaccine A - 0 month","Vaccine A - 1 month","Vaccine A - 2 month",
                                        "Vaccine B - 9 month","MPXV - 2-6 month","MPXV - 10 month"),
                             labels = group_labels
)


p_value_set$group2 <- factor(p_value_set$group2,
                             levels = c("pediatric", "rubella IgG+",
                                        "Vaccine A - 0 month","Vaccine A - 1 month","Vaccine A - 2 month",
                                        "Vaccine B - 9 month","MPXV - 2-6 month","MPXV - 10 month"),
                             labels = group_labels
)

dat_MSD_wide$group_plot <- factor(dat_MSD_wide$group, 
                                  levels = c("pediatric", "rubella IgG+",
                                             "Vaccine A - 0 month","Vaccine A - 1 month","Vaccine A - 2 month",
                                             "Vaccine B - 9 month","MPXV - 2-6 month","MPXV - 10 month"),
                                  labels = group_labels
)


means_median_labels <- dat_MSD_wide %>%
  group_by(antigen_group, group_plot) %>%
  summarize(
    mean_ratio = mean(log2(ratio)),
    median_ratio = median(log2(ratio))
  ) %>%
  mutate(
    label_mean = format(round(2^mean_ratio,2), nsmall = 2),
    label_median = format(round(2^median_ratio,2), nsmall = 2),
    direction = if_else(mean_ratio > 0, "pos", "neg")
  )

color_labels <- function(x) {
  
  x <- tibble(breaks = x)
  
  x <- x %>%
    mutate(
      breaks_color = case_when(
        breaks == 0 ~ paste0("<span style = 'color:#000000;'>", breaks, "</span>"),
        breaks < 0 ~ paste0("<span style = 'color:#B15928;'>", breaks, "</span>"),
        breaks > 0 ~ paste0("<span style = 'color:#1F78B4;'>", breaks, "</span>"),
      )
    )
  return(x$breaks_color)
  
}


#create a tible for just the facet labels.  
facet_labels <- tibble(
  antigen_group =  c("MPXV A29L/VACV A27L",
                     "MPXV A35R/VACV A33R",
                     "MPXV B6R/VACV B5R",
                     "MPXV E8L/VACV D8L",
                     "MPXV M1R/VACV L1R"),
  label = c(
    "<span style = 'color:#1F78B4;'>MPXV A29L</span><span style = 'color:#000000;'> / </span><span style = 'color:#B15928;'>VACV A27L</span>",
    "<span style = 'color:#1F78B4;'>MPXV A35R</span><span style = 'color:#000000;'> / </span><span style = 'color:#B15928;'>VACV A33R</span>",
    "<span style = 'color:#1F78B4;'>MPXV B6R</span><span style = 'color:#000000;'> / </span><span style = 'color:#B15928;'>VACV B5R</span>",
    "<span style = 'color:#1F78B4;'>MPXV E8L</span><span style = 'color:#000000;'> / </span><span style = 'color:#B15928;'>VACV D8L</span>",
    "<span style = 'color:#1F78B4;'>MPXV M1R</span><span style = 'color:#000000;'> / </span><span style = 'color:#B15928;'>VACV L1R</span>"
  ),
  x = 0, xend = 9
)

label_map <- setNames(facet_labels$label, facet_labels$antigen_group)

## For adjusting the position of the y-axis label
adj <- 2.5


## filter p_value_set for only values to be plotted, set y.position for plotting
p_value_set_fil <- p_value_set %>%
  filter(in_figure) %>%
  group_by(antigen_group) %>%
  mutate(
    y = seq(from = 1, length.out = length(antigen_group), by = 1.025),
    overall_max = max(max),
    y.position = y + overall_max
  ) %>%
  ungroup() 

Fig_3 <- ggplot(dat_MSD_wide, aes(x = group_plot, y = log2(ratio))) +
  facet_wrap(
    vars(antigen_group),
    ncol = 5,
    labeller = as_labeller(label_map)
  ) +
  geom_segment(
    data = facet_labels,
    aes(x = x, xend = xend, y = 13.9, yend = 13.9),
    inherit.aes = FALSE,
    linewidth = 0.5
  ) +
  geom_quasirandom(
    shape = 21,
    size = 0.75,
    color = "black"
  ) +
  geom_boxplot(
    alpha = 0.5,
    linewidth = 0.25,
    fill = "lightblue",
    color = "blue",
    outliers = FALSE
  ) +
  scale_y_continuous(
    breaks = seq(-8,13,1),
    labels = color_labels(seq(-8,13)),
    limits = c(-10,14),
  ) +
  scale_x_discrete(
    name = "Cohort",
    labels = group_labels
  ) +
  coord_cartesian(
    clip = "off",
    xlim = c(0,8),
    ylim = c(-10,13)
  ) +
  #### Add mean labels
  geom_text(
    data = tibble(antigen_group = "MPXV A29L/VACV A27L", label = "Mean Ratio:"),
    inherit.aes = FALSE,
    aes(label = label),
    x = -3,
    y = -9.5,
    size = 7/.pt
  ) +
  ###### print mean ratio values on plot and color them according to direction
  geom_text(
    data = means_median_labels,
    inherit.aes = FALSE,
    aes(x = group_plot, y = -9.5, label = label_mean, color = direction),
    size = 7/.pt,
    angle = 90,
    show.legend = FALSE
  ) +
  scale_color_manual(
    values = c("pos" = "#1F78B4", "neg" = "#B15928")
  ) +
  
  ###### plot P-values
  stat_pvalue_manual(
    data = p_value_set_fil,
    label = "{pvalue(p.adj)}",
    inherit.aes = FALSE,
    tip.length = 0,
    size = 6/.pt,
  ) +
  
  ####### Create custom y-axis title
  geom_segment(
    data = tibble(antigen_group = "MPXV A29L/VACV A27L"),
    x = -3.4,
    xend = -3.4,
    y = 4.2 + adj,
    yend = -4.2 + adj,
    linewidth = 0.5
  ) +
  geom_text(
    data = tibble(antigen_group = "MPXV A29L/VACV A27L"),
    x = -3.9,
    y = 0 + adj,
    label = "MPXV Ortholog (AU/mL)",
    color = "#1F78B4",
    angle = 90,
    size = 8/.pt
  ) +
  geom_text(
    data = tibble(antigen_group = "MPXV A29L/VACV A27L"),
    x = -2.9,
    y = 0 + adj,
    label = "VACV Ortholog (AU/mL)",
    color = "#B15928",
    angle = 90,
    size = 8/.pt
  ) +
  geom_text(
    data = tibble(antigen_group = "MPXV A29L/VACV A27L"),
    x = -3.6,
    y = -5.2 + adj,
    label = "(",
    angle = 90,
    size = 15/.pt
  ) +
  geom_text(
    data = tibble(antigen_group = "MPXV A29L/VACV A27L"),
    x = -3.6,
    y = 5.2 + adj,
    label = ")",
    angle = 90,
    size = 15/.pt
  ) +
  geom_text(
    data = tibble(antigen_group = "MPXV A29L/VACV A27L"),
    x = -3.6,
    y = -6.3 + adj,
    label = "log[2]",
    angle = 90,
    size = 9/.pt,
    parse = TRUE
  ) +
  theme(
    #remove panel background and grid
    panel.background = element_blank(),
    panel.grid = element_blank(),
    
    #to allow using markdown for facet lables
    strip.text = element_markdown(size = 7, margin = margin(l = 8, unit = "pt")),
    strip.background = element_blank(),
    
    axis.text.x = element_text(angle = 90, size = 8, hjust = 1, vjust = 0.5),
    axis.line = element_line(color = "black"),
    legend.title = element_text(size = 8),
    axis.title.x = element_text(size = 10),
    axis.title.y = element_blank(),
    axis.text.y = element_markdown(size = 8),
    # Increase plot margin for custom y-axis label
    plot.margin = margin(l = 30, t = 20, r = 10, unit = "pt")
  )

ggsave("fig_output/Fig_3.pdf", width = 7.5, height = 5, unit = "in", device = cairo_pdf)


