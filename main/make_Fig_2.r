# Fig 2 - MSD data - categorical scatter plot

# Greninger Lab - Jon Reed
# 2025-Oct-15

source("MPXV_paper_functions.r")

package_list <- c("tidyverse",
                  "ggplot2",
                  "scales",
                  "ggbeeswarm",
                  "ggtext")

invisible(lapply(package_list, check_and_load))

dat_MSD <- read_csv("2025-Oct-15_compiled_data.csv")

#make status a factor variable and order the categories
dat_MSD$group <- factor(dat_MSD$group,
                        levels = c("pediatric",  "rubella IgG+",
                                   "Vaccine A - 0 month", "Vaccine A - 1 month", "Vaccine A - 2 month",
                                   "Vaccine B - 9 month", "MPXV - 2-6 month","MPXV - 10 month"),
                        labels = c("pediatric", "rubella",
                                   "Vaccine A - 0 mo", "Vaccine A - 1 mo", "Vaccine A - 2 mo",
                                   "Vaccine B - 9 mo", "MPXV 2-6 mo", "MPXV 10 mo"))



dat_MSD$Assay <- factor(dat_MSD$Assay, levels = c("MPXV A29L", "VACV A27L",
                                                  "MPXV A35R", "VACV A33R",
                                                  "MPXV B6R", "VACV B5R",
                                                  "MPXV E8L", "VACV D8L",
                                                  "MPXV M1R", "VACV L1R"))


#picking and viewing colors
my_pal <- pal_brewer(palette = "Paired")(12)

dat_MSD <- dat_MSD %>%
  mutate(
    cat = case_when(
      Assay %in% c("MPXV A29L", "VACV A27L") ~ "A29L/A27L",
      Assay %in% c("MPXV A35R", "VACV A33R") ~ "A35R/A33R",
      Assay %in% c("MPXV B6R", "VACV B5R") ~ "B6R/B5R",
      Assay %in% c("MPXV E8L", "VACV D8L") ~ "E8L/D8L",
      Assay %in% c("MPXV M1R", "VACV L1R") ~ "M1R/L1R",
      TRUE ~ NA
    ),
    #markdown for facet labels
    cat_md = case_when(
      Assay %in% c("MPXV A29L", "VACV A27L") ~ "<span style = 'color:#1F78B4;'>MPXV A29L</span><span style = 'color:#000000;'> </span><span style = 'color:#B15928;'>VACV A27L</span>",
      Assay %in% c("MPXV A35R", "VACV A33R") ~ "<span style = 'color:#1F78B4;'>MPXV A35R</span><span style = 'color:#000000;'> </span><span style = 'color:#B15928;'>VACV A33R</span>",
      Assay %in% c("MPXV B6R", "VACV B5R") ~ "<span style = 'color:#1F78B4;'>MPXV B6R</span><span style = 'color:#000000;'> </span><span style = 'color:#B15928;'>VACV B5R</span>",
      Assay %in% c("MPXV E8L", "VACV D8L") ~ "<span style = 'color:#1F78B4;'>MPXV E8L</span><span style = 'color:#000000;'> </span><span style = 'color:#B15928;'>VACV D8L</span>",
      Assay %in% c("MPXV M1R", "VACV L1R") ~ "<span style = 'color:#1F78B4;'>MPXV M1R</span><span style = 'color:#000000;'> </span><span style = 'color:#B15928;'>VACV L1R</span>",
      TRUE ~ NA
    ),
    ortholog = if_else(str_detect(Assay, "MPXV"), "MPXV", "VACV"),
    ortho_group_plot = paste0(ortholog, "_", group),
    ortho_group_num = case_when(
      ortho_group_plot == "MPXV_pediatric" ~ 1,
      ortho_group_plot == "VACV_pediatric" ~ 2,
      ortho_group_plot == "MPXV_rubella" ~ 3,
      ortho_group_plot == "VACV_rubella" ~ 4,
      ortho_group_plot == "MPXV_Vaccine A - 0 mo" ~ 6,
      ortho_group_plot == "VACV_Vaccine A - 0 mo" ~ 7,
      ortho_group_plot == "MPXV_Vaccine A - 1 mo" ~ 8,
      ortho_group_plot == "VACV_Vaccine A - 1 mo" ~ 9,
      ortho_group_plot == "MPXV_Vaccine A - 2 mo" ~ 10,
      ortho_group_plot == "VACV_Vaccine A - 2 mo" ~ 11,
      ortho_group_plot == "MPXV_Vaccine B - 9 mo" ~ 13,
      ortho_group_plot == "VACV_Vaccine B - 9 mo" ~ 14,
      ortho_group_plot == "MPXV_MPXV 2-6 mo" ~ 16,
      ortho_group_plot == "VACV_MPXV 2-6 mo" ~ 17,
      ortho_group_plot == "MPXV_MPXV 10 mo" ~ 18,
      ortho_group_plot == "VACV_MPXV 10 mo" ~ 19,
      TRUE ~ 25
    )
  ) %>%
  nest(data = !c(Assay)) %>%
  mutate(
    baseline = map_dbl(data, function(x) {exp(mean(log(x %>% filter(group == "rubella") %>% pull(MSD_mean))))}),
  ) %>%
  unnest(data) %>%
  mutate(
    fold_baseline = MSD_mean/baseline,
    .after = MSD_mean
  )


#order the display of the cohorts
dat_MSD$ortho_group_plot <- factor(dat_MSD$ortho_group_plot, c("MPXV_pediatric","VACV_pediatric",
                                                               "MPXV_rubella","VACV_rubella",
                                                               "MPXV_Vaccine A - 0 mo","VACV_Vaccine A - 0 mo",
                                                               "MPXV_Vaccine A - 1 mo","VACV_Vaccine A - 1 mo",
                                                               "MPXV_Vaccine A - 2 mo","VACV_Vaccine A - 2 mo",
                                                               "MPXV_Vaccine B - 9 mo","VACV_Vaccine B - 9 mo",
                                                               "MPXV_MPXV 2-6 mo","VACV_MPXV 2-6 mo",
                                                               "MPXV_MPXV 10 mo","VACV_MPXV 10 mo"))

#create a tible for just the facet labels.  Could probably just create the markdown here instead of above.
facet_labels <- tibble(label = dat_MSD %>% pull(cat_md) %>% unique(), cat = dat_MSD %>% pull(cat) %>% unique())

custom_lab <- function(x) {
  
  #Helper function to format number labels:
  # <1 - format for two figures after the decimal (e.g. 0.25)
  # >= 1 - format with comma, no figures after the decimal 
  y <- c()
  
  for(i in 1:length(x)) {
    if(x[i] < 1) {
      y[i] <- round(x[i],2)
    } else {
      y[i] <- comma(x[i])
    }
  }
  return(y)
}

#define fill colors for h/o vaccination
prior_vac_col = c("yes" = "red", "no" = "white")
dat_MSD$vaccine_history <- factor(dat_MSD$vaccine_history, levels = c("yes", "no"))


# data frame for defining annotations along the bottom of the X-axis.
x_axis_annotation <- tribble(
  ~group, ~cat, ~text, ~x_text, ~y_text, ~x, ~y, ~xend, ~yend,
  "negative", "M1R/L1R","Negative Cohorts",2.5,-9.5,0.5,-8,4.5,-8,
  "vac_cohort_A", "M1R/L1R", "Vaccine Cohort A",8.5,-9.5,5.5,-8,11.5,-8,
  "vac_cohort_B","M1R/L1R","Vaccine Cohort B",13.5,-9.5,12,-8,15,-8,
  "MPXV", "M1R/L1R", "MPXV-infected Cohorts",17.5,-9.5,15.5,-8,19.5,-8
)

#cohort labels
cohort_labels <- c(
  "pediatric", expression(rubella~IgG^"+"), 
  "0 month", "1 month", "2 month", 
  "9 month", "2-6 month", "10 month")


median_plot <- dat_MSD %>%
  group_by(cat, ortho_group_num, ortho_group_plot) %>%
  summarize(
    median = log2(exp(mean(log(MSD_mean))))
  ) %>%
  mutate(
    x = ortho_group_num - 0.25,
    xend = ortho_group_num + 0.25,
    y = median,
    yend = median
  )


Fig2 <- ggplot(dat_MSD %>% filter(!(Sample_ID == "P61" & MSD_dr_final == "not detected")), aes(x = ortho_group_num, y = log2(MSD_mean), group = ortho_group_plot)) +
  #Note: one undetected observation in pediatric set is not plotted
  facet_wrap(
    vars(cat),
    ncol = 1
  ) +
  scale_color_manual(
    values = c(my_pal[2], my_pal[12])
  ) +
  scale_fill_manual(
    values = prior_vac_col
  ) +
  geom_violin(
    trim = FALSE,
    color = "grey42",
    linewidth = 0.1
  ) +
  geom_quasirandom(
    aes(color = ortholog, fill = vaccine_history),
    method = "tukeyDense",
    shape = 21,
    size = 0.5,
    stroke = 0.25,
    varwidth = TRUE,
  ) +
  #plotting median lines
  geom_segment(
    data = median_plot,
    aes(x = x, xend = xend, y = y, yend = yend),
    color = "grey20",
    linewidth = 0.3
  ) +
  scale_y_continuous(
    label = ~custom_lab(2^.x),
    breaks = seq(-2,19,2),
    limits = c(-10,20),
    name = "Antibody Level (AU/mL)"
  ) +
  scale_x_continuous(
    breaks = c(1.5, 3.5, 6.5, 8.5, 10.5, 13.5, 16.5, 18.5),
    expand = c(0,0),
    labels = cohort_labels
  ) +
  #turn of clipping to add additional x-axis annotation outside of plot
  coord_cartesian(
    clip = "off",
    ylim = c(-3, 20)
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
    size = 7/.pt
  ) +
  #x-axis annotation for vaccine cohort A
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
    size = 7/.pt
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
    aes(x = x_text, y = y_text, label = text),
    size = 7/.pt
  ) +
  #x-axis annotation for MPXV infected
  geom_segment(
    inherit.aes = FALSE,
    data = x_axis_annotation %>% filter(group == "MPXV"),
    aes(x = x, y = y, xend = xend, yend = yend),
    linewidth = 0.25
  ) +
  geom_text(
    inherit.aes = FALSE,
    data = x_axis_annotation %>% filter(group == "MPXV"),
    aes(x = x_text, y = y_text, label = text),
    size = 7/.pt
  ) +
  guides(
    color = "none",
    #to increase size of points in legend
    fill = guide_legend(override.aes = list(size=1.5), title = "Smallpox vaccination history:")
  ) +
  geom_richtext(
    data = facet_labels,
    inherit.aes = FALSE,
    label.color = NA,
    aes(y = log2(2^19.5), x = 0.1, label = label, hjust = 0, vjust = 0.5),
    size = 7/.pt
  ) +
  theme(
    #remove these elements
    panel.background = element_blank(),
    panel.grid = element_blank(),
    strip.background = element_blank(),
    strip.text.x = element_blank(),
    #plot margin on bottom increased for extra x-axis annotation
    plot.margin = margin(t = 5, r = 5, b = 15, l = 5, unit = "pt"),
    
    #axis adjustments
    axis.line = element_line(),
    axis.text.y = element_text(size = 7, hjust = 1),
    axis.title.y = element_text(size = 8),
    axis.text.x = element_text(vjust = 0, size = 7, margin = margin(t = -2)),
    axis.title.x = element_blank(),
    
    #legends adjustments
    legend.title = element_text(size = 8, margin = margin(t = -1, r = 0, b = 0, l = 0)),
    legend.text = element_text(size = 8, margin = margin(t = -1, l = -3, unit = "pt")),
    legend.position = "top",
    legend.justification = "left",
    legend.margin = margin(0,0,0,0),
    legend.box.margin = margin(0,0,0,0),
    legend.box.spacing = unit(0, "pt"),
  )

ggsave("fig_output/Fig_2.pdf", height = 7, width = 5.5, unit = "in", device = cairo_pdf)


