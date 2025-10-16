# Fig S12 ROC analysis of MPXV-infected vs. unexposed and vaccinated vs. unexposed
# Tables S5-S6

# Greninger Lab - Jon Reed
# 2025-Oct-15

source("MPXV_paper_functions.r")

package_list <- c("tidyverse",
                  "ggplot2",
                  "gt",
                  "pROC",
                  "scales",
                  "patchwork",
                  "openxlsx2")

invisible(lapply(package_list, check_and_load))

dat_neut_MSD <- read_csv("2025-Oct-15_compiled_data.csv")


#set for CI bootstrapping, needs to be 10,000 for final figure, but can
#use smaller numbers for testing
boot_num = 10000


dat_neut_MSD <- dat_neut_MSD %>%
  mutate(
    vaccine_case = case_when(
      group %in% c("Vaccine A - 1 month", "Vaccine A - 2 month", "Vaccine B - 9 month") ~ TRUE,
      group == "Vaccine A - 0 month" & vaccine_history == "yes" ~ TRUE,
      TRUE ~ FALSE
    ),
    MPXV_case = case_when(
      group %in% c("MPXV - 10 month", "MPXV - 2-6 month") ~ TRUE,
      TRUE ~ FALSE
    ),
    negative_case = case_when(
      group %in% c("pediatric", "rubella IgG+") ~ TRUE,
      group == "Vaccine A - 0 month" & vaccine_history == "no" ~ TRUE,
      TRUE ~ FALSE
    )
  )




#A - MPXV-infected vs. pediatric/rubella
#===============================================================================

dat_ROC <- dat_neut_MSD %>%
  ungroup() %>%
  filter(
    #!group %in% c("Vaccine A - 0 month", "Vaccine A - 1 month", "Vaccine A - 2 month", "Vaccine B - 9 month")
    MPXV_case | negative_case
  ) %>%
  mutate(
    cases = if_else(MPXV_case == TRUE, TRUE, FALSE),
    count = 1
  ) 

count_A <- dat_ROC %>%
  select(group, cases, count)%>%
  pivot_wider(
    id_cols = group,
    names_from = cases,
    values_from = count,
    values_fn = sum
  ) %>%
  mutate(
    `TRUE` = if_else(is.na(`TRUE`), 0, `TRUE`/10),
    `FALSE` = if_else(is.na(`FALSE`), 0, `FALSE`/10),
  ) %>%
  rowwise() %>%
  mutate(
    total = sum(`TRUE`, `FALSE`)
  ) %>%
  gt() %>%
  grand_summary_rows(fns = list(TOTAL = "sum"), 
                     columns = c(`TRUE`, `FALSE`, total)) 


####### best.policy = "random" was set for determining the CI, due multiple best
####### threshold identified for some of the resampled ROC curves.
####### I think the random selection will avoid bias, since it will not favor
####### one choice over another - plus curves are resampled 10,000 times.

ROC <- dat_ROC %>%
  nest(data = !Assay) %>%
  mutate(
    roc_curve = map(data, ~roc(.x$cases, .x$MSD_mean, levels = c(TRUE, FALSE), direction = ">")),
    values = map(roc_curve, ~tibble(sensitivity = .x$sensitivities, specificity = .x$specificities)),
    AUC = map_dbl(roc_curve, auc),
    best_threshold = map_dbl(roc_curve, ~coords(.x, x = "best", best.method = "youden")[["threshold"]]),
    best_specificity = map_dbl(roc_curve, ~coords(.x, x = "best", best.method = "youden")[["specificity"]]),
    best_sensitivity = map_dbl(roc_curve, ~coords(.x, x = "best", best.method = "youden")[["sensitivity"]]),
    CIs = map(roc_curve, ~ci.coords(.x, x = "best", best.policy = "random", ret = c("specificity", "sensitivity"), boot.n = boot_num)),
    best_specificity_range = map(CIs, ~.x$specificity[c(1,3)]),
    best_sensitivity_range = map(CIs, ~.x$sensitivity[c(1,3)]),
    AUC_range = pmap(list(roc_curve, AUC), ~ci.auc(.x, boot.n = boot_num)[c(1,3)])
  ) 

#format as ratio
format_CI <- function(x) {
  paste0("[",format(round(x[[1]],3), nsmall = 3), " - ", format(round(x[[2]],3), nsmall = 3),"]")
}


ROC$Assay <- factor(
  ROC$Assay, 
  levels = c(
    "MPXV A35R", "VACV A33R",
    "MPXV B6R", "VACV B5R",
    "MPXV E8L", "VACV D8L",
    "MPXV M1R", "VACV L1R",
    "MPXV A29L", "VACV A27L"
  )
)


A_Table <- ROC %>%
  select(Assay, AUC, best_threshold, best_specificity, best_specificity_range, best_sensitivity, best_sensitivity_range, AUC_range) %>%
  rowwise() %>%
  mutate(
    best_threshold = comma(best_threshold, accuracy = 0.1),
    best_specificity_range = format_CI(best_specificity_range),
    best_sensitivity_range = format_CI(best_sensitivity_range),
    AUC_range = format_CI(AUC_range),
    best_specificity = format(round(best_specificity, 3), nsmall = 3),
    best_sensitivity = format(round(best_sensitivity, 3), nsmall = 3),
    AUC = format(round(AUC, 3), nsmall = 3)
  ) %>%
  rename(
    `AUC 95% CI` = AUC_range,
    `specificity 95% CI` = best_specificity_range,
    `sensitivity 95% CI` = best_sensitivity_range,
    `threshold (AU/mL)` = best_threshold, 
    specificity = best_specificity,
    sensitivity = best_sensitivity
  ) %>%
  select(Assay, AUC, `AUC 95% CI`, `threshold (AU/mL)`, specificity, `specificity 95% CI`, sensitivity, `sensitivity 95% CI`) %>%
  arrange(Assay)



ROC_plot <- ROC %>%
  select(Assay, values) %>%
  unnest(values)


cols <- brewer_pal(palette = "Paired")(10)

colors <- c(
  "MPXV A35R" = cols[1], "VACV A33R" = cols[2], #blue pair
  "MPXV B6R" = cols[3], "VACV B5R" = cols[4], #green pair
  "MPXV E8L" = cols[9], "VACV D8L" = cols[10], #purple pair
  "MPXV M1R" = cols[7], "VACV L1R" = cols[8], #orange pair
  "MPXV A29L" = "hotpink", "VACV A27L" = cols[6] #red pair
)

A <- ggplot(ROC_plot, aes(x = 1-specificity, y = sensitivity, color = Assay)) +
  geom_path(
    linewidth = 0.3
  ) +
  scale_color_manual(
    values = colors,
    name = "Antigen"
  ) +
  geom_abline(
    slope = 1,
    linetype = "dashed",
    color = "grey"
  ) +
  labs(
    title = "MPXV-infected vs. Negative",
    subtitle = "Individual Antigen",
  ) +
  xlab("1-Specificity") +
  ylab("Sensitivity") +
  theme(
    aspect.ratio = 1,
    panel.background = element_blank(),
    panel.grid = element_blank(),
    axis.line = element_line(color = "black"),
    plot.title = element_text(size = 8, face = "bold"),
    plot.subtitle = element_text(size = 8),
    axis.text = element_text(size = 8),
    axis.title = element_text(size = 9),
    plot.margin = unit(c(t = 0, r = 60, b = 0, l = 0), "pt")
  )

#B - vaccinated vs. pediatric/rubella 
#===============================================================================

dat_ROC <- dat_neut_MSD %>%
  ungroup() %>%
  filter(
    !group %in% c("MPXV - 2-6 month", "MPXV - 10 month")
  ) %>%
  mutate(
    cases = if_else(vaccine_case == TRUE, TRUE, FALSE),
    count = 1
  ) 

count_B <- dat_ROC %>%
  select(group, cases, count)%>%
  pivot_wider(
    id_cols = group,
    names_from = cases,
    values_from = count,
    values_fn = sum
  ) %>%
  mutate(
    `TRUE` = if_else(is.na(`TRUE`), 0, `TRUE`/10),
    `FALSE` = if_else(is.na(`FALSE`), 0, `FALSE`/10),
  ) %>%
  rowwise() %>%
  mutate(
    total = sum(`TRUE`, `FALSE`)
  ) %>%
  gt() %>%
  grand_summary_rows(fns = list(TOTAL = "sum"), 
                     columns = c(`TRUE`, `FALSE`, total)) 


####### best.policy = "random" was set for determining the CI, due multiple best
####### threshold identified for some of the resampled ROC curves.
####### I think the random selection will avoid bias, since it will not favor
####### one choice over another - plus curves are resampled 10,000 times.

ROC <- dat_ROC %>%
  nest(data = !Assay) %>%
  mutate(
    roc_curve = map(data, ~roc(.x$cases, .x$MSD_mean, levels = c(TRUE, FALSE), direction = ">")),
    values = map(roc_curve, ~tibble(sensitivity = .x$sensitivities, specificity = .x$specificities)),
    AUC = map_dbl(roc_curve, auc),
    best_threshold = map_dbl(roc_curve, ~coords(.x, x = "best", best.method = "youden")[["threshold"]]),
    best_specificity = map_dbl(roc_curve, ~coords(.x, x = "best", best.method = "youden")[["specificity"]]),
    best_sensitivity = map_dbl(roc_curve, ~coords(.x, x = "best", best.method = "youden")[["sensitivity"]]),
    CIs = map(roc_curve, ~ci.coords(.x, x = "best", best.policy = "random", ret = c("specificity", "sensitivity"), boot.n = boot_num)),
    best_specificity_range = map(CIs, ~.x$specificity[c(1,3)]),
    best_sensitivity_range = map(CIs, ~.x$sensitivity[c(1,3)]),
    AUC_range = pmap(list(roc_curve, AUC), ~ci.auc(.x, boot.n = boot_num)[c(1,3)])
  ) 

#format as ratio
format_CI <- function(x) {
  paste0("[",format(round(x[[1]],3), nsmall = 3), " - ", format(round(x[[2]],3), nsmall = 3),"]")
}

ROC$Assay <- factor(
  ROC$Assay, 
  levels = c(
    "MPXV A35R", "VACV A33R",
    "MPXV B6R", "VACV B5R",
    "MPXV E8L", "VACV D8L",
    "MPXV M1R", "VACV L1R",
    "MPXV A29L", "VACV A27L"
  )
)


B_Table <- ROC %>%
  select(Assay, AUC, best_threshold, best_specificity, best_specificity_range, best_sensitivity, best_sensitivity_range, AUC_range) %>%
  rowwise() %>%
  mutate(
    best_threshold = comma(best_threshold, accuracy = 0.1),
    best_specificity_range = format_CI(best_specificity_range),
    best_sensitivity_range = format_CI(best_sensitivity_range),
    AUC_range = format_CI(AUC_range),
    best_specificity = format(round(best_specificity, 3), nsmall = 3),
    best_sensitivity = format(round(best_sensitivity, 3), nsmall = 3),
    AUC = format(round(AUC, 3), nsmall = 3)
  ) %>%
  rename(
    `AUC 95% CI` = AUC_range,
    `specificity 95% CI` = best_specificity_range,
    `sensitivity 95% CI` = best_sensitivity_range,
    `threshold (AU/mL)` = best_threshold, 
    specificity = best_specificity,
    sensitivity = best_sensitivity
  ) %>%
  select(Assay, AUC, `AUC 95% CI`, `threshold (AU/mL)`, specificity, `specificity 95% CI`, sensitivity, `sensitivity 95% CI`) %>%
  arrange(Assay)



ROC_plot <- ROC %>%
  select(Assay, values) %>%
  unnest(values)


cols <- brewer_pal(palette = "Paired")(10)

colors <- c(
  "MPXV A35R" = cols[1], "VACV A33R" = cols[2], #blue pair
  "MPXV B6R" = cols[3], "VACV B5R" = cols[4], #green pair
  "MPXV E8L" = cols[9], "VACV D8L" = cols[10], #purple pair
  "MPXV M1R" = cols[7], "VACV L1R" = cols[8], #orange pair
  "MPXV A29L" = "hotpink", "VACV A27L" = cols[6] #red pair
)

B <- ggplot(ROC_plot, aes(x = 1-specificity, y = sensitivity, color = Assay)) +
  geom_path(
    linewidth = 0.3
  ) +
  scale_color_manual(
    values = colors,
    name = "Antigen"
  ) +
  geom_abline(
    slope = 1,
    linetype = "dashed",
    color = "grey"
  ) +
  labs(
    title = "Vaccinated vs. Negative",
    subtitle = "Individual Antigen",
  ) +
  xlab("1-Specificity") +
  ylab("Sensitivity") +
  theme(
    aspect.ratio = 1,
    panel.background = element_blank(),
    panel.grid = element_blank(),
    axis.line = element_line(color = "black"),
    plot.title = element_text(size = 8, face = "bold"),
    plot.subtitle = element_text(size = 8),
    axis.text = element_text(size = 8),
    axis.title = element_text(size = 9),
    plot.margin = unit(c(t = 0, r = 60, b = 0, l = 0), "pt")
  )


#===============================================================================
#Figure
#===============================================================================

layout = c(
  area(t = 1, l = 1, b = 1, r = 4),
  area(t = 1, l = 5, b = 1, r = 8))


Fig_S12 <- A + B +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = "A") &
  theme(
    legend.title = element_text(size = 8, face = "bold"),
    legend.text = element_text(size = 6),
    legend.key.size = unit(10, "pt"),
  )


ggsave("fig_output/Fig_S12.pdf", width = 7, height = 5.5, unit = "in")

table_list <- list("Table S5" = A_Table, "Table S6" = B_Table)
table_list <- map(table_list, as_tibble)
writexl::write_xlsx(table_list, "fig_output/Fig_S12_Tables_S5_S6.xlsx")

