# Fig S10
# Table S3

# Greninger Lab - Jon Reed
# 2025-Oct-15

source("MPXV_paper_functions.r")

package_list <- c("tidyverse",
                  "ggplot2", 
                  "scales",
                  "rstatix",
                  "openxlsx2",
                  "ggtext",
                  "ggnewscale")

invisible(lapply(package_list, check_and_load))

dat_MSD <- read_csv("2025-Oct-15_compiled_data.csv")

#define a few lists for building pairwise tables
cohorts <- c("pediatric",  "rubella IgG+", 
             "Vaccine A - 0 month", "Vaccine A - 1 month", "Vaccine A - 2 month",
             "Vaccine B - 9 month", "MPXV - 2-6 month","MPXV - 10 month")

virus <- c("MPXV", "VACV")
antigen_list = c("MPXV A29L", "VACV A27L", "MPXV A35R", "VACV A33R", "MPXV B6R", "VACV B5R", "MPXV E8L", "VACV D8L", "MPXV M1R", "VACV L1R")
antigen_group_list = c("MPXV A29L, VACV A27L", "MPXV A35R, VACV A33R", "MPXV B6R, VACV B5R", "MPXV E8L, VACV D8L", "MPXV M1R, VACV L1R")

#generate list of virus and cohort combos
virus_cohorts <- expand_grid(cohorts, virus) %>%
  mutate(
    temp = paste0(virus,"_",cohorts)
  ) %>%
  pull(temp)

#get results
results <- dat_MSD %>%
  mutate(
    mean_log = log(MSD_mean)
  ) %>%
  rename(
    cohort = group,
    antigen = Assay
  ) %>%
  select(cohort, antigen, mean_log)

test_pairs <- tibble(
  antigen_group = antigen_group_list,
  data = map(1:length(antigen_group_list), ~tibble(
    factor_1_virus_cohort = map(0:15, ~ virus_cohorts[seq_len(length(virus_cohorts) - .x)]) %>% unlist(),
    factor_2_virus_cohort = map2(16:1, 16:1, ~rep(virus_cohorts[.x], .y)) %>% unlist()
  ))
) %>%
  unnest(data) %>%
  #create columns for joining with results
  mutate(
    factor_1_cohort = str_extract(factor_1_virus_cohort, "[^_]+$"),
    factor_1_virus = str_extract(factor_1_virus_cohort, "^[^_]+"),
    factor_1_antigen = if_else(factor_1_virus == "MPXV", str_extract(antigen_group, "^[^,]+"), str_extract(antigen_group, "(?<=, )[^,]+")),
    factor_2_cohort = str_extract(factor_2_virus_cohort, "[^_]+$"),
    factor_2_virus = str_extract(factor_2_virus_cohort, "^[^_]+"),
    factor_2_antigen = if_else(factor_2_virus == "MPXV", str_extract(antigen_group, "^[^,]+"), str_extract(antigen_group, "(?<=, )[^,]+"))
  ) %>% 
  ## get results
  mutate(
    factor_1_results = map2(factor_1_cohort, factor_1_antigen, function(x,y) {results %>% filter(cohort == x, antigen == y)}),
    factor_2_results = map2(factor_2_cohort, factor_2_antigen, function(x,y) {results %>% filter(cohort == x, antigen == y)})
  ) %>%
  pivot_longer(
    cols = c(factor_1_results, factor_2_results),
    names_to = "factor"
  ) %>%
  unnest(value) %>%
  group_by(antigen_group, factor_1_virus_cohort, factor_2_virus_cohort, factor_1_cohort, factor_2_cohort, factor_1_virus, factor_2_virus, factor_1_antigen, factor_2_antigen) %>%
  t_test(mean_log ~ factor, alternative = "two.sided", var.equal = FALSE, detailed = TRUE) %>%
  adjust_pvalue(method = "bonferroni") %>%
  ## limit upper end of p-value for filling
  mutate(
    p_value_plot = case_when(
      p.adj < 1e-7 ~ 10^-7,
      p.adj <= 0.05 & p.adj >= 1e-7 ~ p.adj,
      TRUE ~ NA_real_
    )
  )


###### colors for gradient fill
colors <- pal_brewer(palette = "YlOrRd")(9)

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

####### write supp table with t-test results

# Configure table for writing to excel
test_pairs_wb <- test_pairs %>% select(
  antigen_group, 
  factor_1_cohort, factor_1_antigen, 
  factor_2_cohort, factor_2_antigen, 
  n1, n2,
  estimate1, estimate2,
  df, statistic,
  p, p.adj
) %>%
  mutate(
    factor_1_cohort = map_chr(factor_1_cohort, ~cohort_key[cohort_key$table_value == .x,]$display_value),
    factor_2_cohort = map_chr(factor_2_cohort, ~cohort_key[cohort_key$table_value == .x,]$display_value)
  ) %>%
  rename(
    `Antigen Group` = antigen_group,
    `Group 1 - cohort` = factor_1_cohort,
    `Group 1 - antigen` = factor_1_antigen,
    `Group 2 - cohort` = factor_2_cohort,
    `Group 2 - antigen` = factor_2_antigen,
    `Group 1 - geometric mean (AU/mL)` = estimate1,
    `Group 1 - n` = n1,
    `Group 2 - geometric mean (AU/mL)` = estimate2,
    `Group 2 - n` = n2,
    `degrees of freedem` = df,
    `t-statistic` = statistic,
    `P-value` = p,
    `P-value adj. MC` = p.adj
  ) %>%
  mutate(
    `log10(P-value adj. MC)` = log10(`P-value adj. MC`),
    `Group 1 - geometric mean (AU/mL)` = exp(`Group 1 - geometric mean (AU/mL)`),
    `Group 2 - geometric mean (AU/mL)` = exp(`Group 2 - geometric mean (AU/mL)`),
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
    cols = 14
  ),
  style = "grey_fill",
  rule = "N2 > log10(0.05)",
  type = "expression"
)

wb <- wb_add_conditional_formatting(
  wb,
  sheet = "Sheet 1",
  dims = wb_dims(
    rows = 2:(nrow(test_pairs_wb)+1),
    cols = 14
  ),
  style = c(colors[9], colors[1]),
  rule = c(log10(1e-7), log10(0.049999)),
  type = "colorScale"
)

# Reduce the number of digits after decimal to 2, for certain columns
wb <- wb_add_numfmt(
  wb,
  sheet = "Sheet 1",
  dims = wb_dims(
    rows = 2:(nrow(test_pairs_wb)+1),
    cols = c(8:11,14)
  ),
  numfmt = "0.00"
)

# Configure p-values to scientific
wb <- wb_add_numfmt(
  wb,
  sheet = "Sheet 1",
  dims = wb_dims(
    rows = 2:(nrow(test_pairs_wb)+1),
    cols = 12:13
  ),
  numfmt = "0.00E+00"
)

#Center columns
wb <- wb_add_cell_style(
  wb,
  sheet = "Sheet 1",
  dims = wb_dims(
    rows = 1:(nrow(test_pairs_wb)+1),
    cols = c(3,5,6:14)
  ),
  vertical = "center",
  horizontal = "center"
)

#Auto adjust column widths
wb_set_col_widths(
  wb,
  sheet = "Sheet 1",
  cols = c(1:14),
  widths = "auto"
)

#Freeze the top row
wb_freeze_pane(wb, sheet = "Sheet 1", first_active_row = 2, first_row = TRUE)

wb_save(wb, "fig_output/Table_S3.xlsx")

########

test_pairs$factor_1_virus_cohort <- 
  factor(
    test_pairs$factor_1_virus_cohort, 
    virus_cohorts
  )

test_pairs$factor_2_virus_cohort <- 
  factor(
    test_pairs$factor_2_virus_cohort, 
    rev(virus_cohorts)
  )

labels <- c("<span style = 'color:#1F78B4;'>pediatric</span>", 
            "<span style = 'color:#B15928;'>pediatric</span>", 
            "<span style = 'color:#1F78B4;'>rubella IgG<sup>+</sup></span>",
            "<span style = 'color:#B15928;'>rubella IgG<sup>+</sup></span>", 
            "<span style = 'color:#1F78B4;'>Vac A - 0</span>",
            "<span style = 'color:#B15928;'>Vac A - 0</span>", 
            "<span style = 'color:#1F78B4;'>Vac A - 1 mo</span>",
            "<span style = 'color:#B15928;'>Vac A - 1 mo</span>", 
            "<span style = 'color:#1F78B4;'>Vac A - 2 mo</span>",
            "<span style = 'color:#B15928;'>Vac A - 2 mo</span>", 
            "<span style = 'color:#1F78B4;'>Vac B - 9 mo</span>",
            "<span style = 'color:#B15928;'>Vac B - 9 mo</span>", 
            "<span style = 'color:#1F78B4;'>MPXV - 2 - 6 mo</span>",
            "<span style = 'color:#B15928;'>MPXV - 2 - 6 mo</span>", 
            "<span style = 'color:#1F78B4;'>MPXV - 10 mo</span>",
            "<span style = 'color:#B15928;'>MPXV - 10 mo</span>") 

test_pairs$antigen_group <- factor(test_pairs$antigen_group, c("MPXV A29L, VACV A27L", "MPXV A35R, VACV A33R", "MPXV B6R, VACV B5R", "MPXV E8L, VACV D8L", "MPXV M1R, VACV L1R"))

facet_labels <- tribble(~antigen_group, ~label,
                        "MPXV A29L, VACV A27L", "<span style = 'color:#1F78B4;'>MPXV A29L</span><span style = 'color:#000000;'> </span><span style = 'color:#B15928;'>VACV A27L</span>",
                        "MPXV A35R, VACV A33R", "<span style = 'color:#1F78B4;'>MPXV A35R</span><span style = 'color:#000000;'> </span><span style = 'color:#B15928;'>VACV A33R</span>",
                        "MPXV B6R, VACV B5R", "<span style = 'color:#1F78B4;'>MPXV B6R</span><span style = 'color:#000000;'> </span><span style = 'color:#B15928;'>VACV B5R</span>",
                        "MPXV E8L, VACV D8L", "<span style = 'color:#1F78B4;'>MPXV E8L</span><span style = 'color:#000000;'> </span><span style = 'color:#B15928;'>VACV D8L</span>",
                        "MPXV M1R, VACV L1R", "<span style = 'color:#1F78B4;'>MPXV M1R</span><span style = 'color:#000000;'> </span><span style = 'color:#B15928;'>VACV L1R</span>")

Fig_S10 <- ggplot(test_pairs, aes(x = factor_1_virus_cohort, y = factor_2_virus_cohort)) + 
  facet_wrap(
    vars(antigen_group),
    axes = "all",
    axis.labels = "margins",
    ncol = 1
  ) +
  geom_tile(
    color = "black",
    fill = NA
  ) +
  geom_tile(
    data = test_pairs %>% filter(p.adj <= 0.05),
    aes(fill = -log10(p_value_plot)),
    color = "black"
  ) +
  scale_x_discrete(
    labels = labels,
    expand = expansion(mult = c(0.02,0.06)),
    name = "Cohort"
  ) +
  scale_y_discrete(
    labels = rev(labels),
    expand = expansion(mult = c(0,0.1)),
    name = "Cohort"
  ) +
  geom_richtext(
    data = facet_labels,
    inherit.aes = FALSE,
    label.color = NA,
    fill = NA,
    aes(y = 17.2, x = 0, label = label, hjust = 0, vjust = 0.5),
    size = 6/.pt
  ) +
  scale_fill_gradient(
    low = colors[1],
    high = colors[9],
    na.value = "green",   #There should be no NA values, highlight with green
    name = "p-value",
    breaks = -log10(c(1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7)),
    labels = c(expression(10^-2), expression(10^-3), expression(10^-4), expression(10^-5), expression(10^-6), expression(paste("< ",10^-7))),
    guide = guide_colorbar(order = 1)
  ) +
  coord_cartesian(
    clip = "off",
    ylim = c(0,16)
  ) +
  new_scale_fill() +
  scale_fill_manual(
    values = "grey",
    breaks = 0,
    name = NULL,
    label = "> 0.05",
    guide = guide_legend(order = 2)
  ) +
  geom_tile(
    data = test_pairs %>% filter(p.adj > 0.05),
    aes(fill = factor(0)),
    color = "black"
  ) +
  theme(
    axis.text.x = element_markdown(size=5, angle=90, vjust=0.5, hjust = 1),
    axis.text.y = element_markdown(size=5),
    axis.title = element_text(size = 8),
    strip.background = element_blank(),
    strip.text.x = element_blank(),
    panel.spacing = unit(10, unit = "pt"),
    panel.border = element_blank(),
    aspect.ratio = 1,
    legend.text = element_text(size = 7),
    legend.title = element_text(size = 7),
    legend.key.size = unit(10, unit = "pt"),
    panel.background = element_blank(),
    panel.grid = element_blank(),
    plot.margin = margin(t = 10, unit = "pt"),
    axis.line = element_line(color = "black"),
  ) 

ggsave("fig_output/Fig_S10.pdf", height = 7, width = 3, unit = "in", device = cairo_pdf)

