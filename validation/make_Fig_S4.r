# Fig S4 - linearity and parallelism


# Greninger Lab - Jon Reed
# 2025-Oct-15

source("../main/MPXV_paper_functions.r")

package_list <- c("tidyverse",
                  "ggplot2",
                  "rcartocolor",
                  "scales",
                  "patchwork")

invisible(lapply(package_list, check_and_load))

dat <- read_csv("2025-Oct-15_MPXV_paper_validation_data.csv")
dat_assay_para <- read_csv("assay_limits/assay_para_validated.csv")

######################## linear regression of calibrators ######################

#pull out calibrator data, assign days to the validation runs
dat_CAL <- dat %>%
  filter(str_detect(Sample, "CAL-0[1-8]")) %>%
  filter(!is.na(day))


#Add LoD, LLoQ, and ULoQ values to the table
dat_CAL <- full_join(dat_CAL, dat_assay_para %>% select(Assay, LoD, LLoQ_final, ULoQ_final), by ="Assay") 

dat_CAL <- dat_CAL %>%
  mutate(
    LLOQ_flag = if_else(conc_calc < LLoQ_final, TRUE, FALSE),
    ULOQ_flag = if_else(conc_calc > ULoQ_final, TRUE, FALSE),
    within_AMR = if_else(LLOQ_flag == FALSE & ULOQ_flag == FALSE, TRUE, FALSE)
  )

#perform linear regression
models_CAL <- dat_CAL %>%
  filter(within_AMR == TRUE) %>%
  group_by(Assay) %>%
  nest() %>%
  mutate(
    model = map(data, ~lm(log10(conc_calc) ~ log10(conc_exp), data = .x)),
    model_tidy = map(model, ~broom::tidy(.x)),
    model_data = map(model, ~broom::augment(.x)),
    intercept = map_dbl(model_tidy, ~.x[[1,2]]),
    intercept_se = map_dbl(model_tidy, ~.x[[1,3]]),
    intercept_pvalue = map_dbl(model_tidy, ~.x[[1,5]]),
    slope = map_dbl(model_tidy, ~.x[[2,2]]),
    slope_se = map_dbl(model_tidy, ~.x[[2,3]]),
    R_sq = map_dbl(model, ~summary(.x)$r.squared)
  )

models_CAL_simple <- models_CAL %>%
  select(Assay, intercept, slope, R_sq) %>%
  mutate(
    intercept = if_else(abs(intercept) < 0.1, "< +/- 0.1 AU/mL", paste0(round(intercept,2), "AU/mL")),
    slope = formatC(slope, format = "f", digits = 2),
    R_sq = formatC(R_sq, format = "f", digits = 2),
    label = paste0("y-int.: ", intercept, "\nslope: ", slope, "\nR^2: ", R_sq),
  )

models_CAL_annot <- models_CAL_simple %>%
  select(Assay,label)

#data frame for plotting limit lines
lines <- dat_assay_para %>%
  select(Assay, LoD, LLoQ_final, ULoQ_final) %>%
  pivot_longer(
    cols = c("LoD", "LLoQ_final", "ULoQ_final"),
    names_to = "Limit"
  ) %>%
  mutate(
    Limit = factor(
      Limit,
      levels = c("ULoQ_final", "LLoQ_final", "LoD"),
      labels = c("ULoQ", "LLoQ", "LoD")
    )
  )

dat_CAL$day <- factor(dat_CAL$day, levels = 1:10)

Fig_4A <- ggplot(dat_CAL %>% filter(!str_detect(Sample, "CAL-08")), aes(x = log2(conc_exp), y = log2(conc_calc), fill = day)) +
  facet_wrap(
    vars(Assay),
    nrow = 2,
    axes = "all_x"
  ) +
  geom_smooth(
    data = dat_CAL %>% filter(within_AMR),
    inherit.aes = FALSE,
    aes(x = log2(conc_exp), y = log2(conc_calc), linetype = "Calibrator\nregression"),
    method = "lm",
    se = FALSE,
    color = "black",
    linewidth = 0.4,
    fullrange = TRUE
  ) +
  scale_linetype_manual(
    values = "dotted",
    name = NULL
  ) +
  geom_point(
    shape = 21,
    size = 2,
    stroke = 0.15,
    alpha = 0.6
  ) +
  ylab(expression(observed~AU/mL)) +
  xlab(expression(expected~AU/mL)) +
  geom_hline(
    data = lines,
    aes(yintercept = log2(value), color = Limit),
    linetype = "dashed",
    linewidth = 0.4
  ) +
  scale_color_manual(
    values = c("LoD" = "red", "LLoQ" = "firebrick4", "ULoQ" = "black")
  ) +
  scale_fill_manual(
    name = "run day",
    values = carto_pal(10, "Safe"),
    drop = FALSE
  ) +
  coord_cartesian(
    clip = "off",
    ylim = c(-12, 8)
  ) +
  scale_y_continuous(
    guide = guide_axis(minor.ticks = TRUE),
    minor_breaks = c(-11,-10, -9, -7, -6, -5, -3, -2, -1, 1, 2, 3, 5, 6, 7),
    limits = c(-30,30),
    breaks = seq(-12,8,4),
    labels = label_math(expr = 2^.x),
  ) +
  scale_x_continuous(
    guide = guide_axis(minor.ticks = TRUE),
    minor_breaks = c(-11,-10, -9, -7, -6, -5, -3, -2, -1, 1, 2, 3, 5, 6, 7),
    limits = c(-12,8),
    breaks = seq(-12,8,4),
    labels = label_math(expr = 2^.x)
  ) +
  geom_text(
    inherit.aes = FALSE,
    data = models_CAL_annot,
    mapping = aes(x = -12, y = 18, label = label),
    size = 8/.pt,
    hjust = 0,
    color = "blue"
  ) +
  theme(
    aspect.ratio = 1,
    
    panel.grid = element_blank(),
    panel.background = element_blank(),
    panel.spacing.y = unit(40, unit = "pt"),
    
    axis.text.x = element_text(angle = 90, vjust = 0.5),
    axis.line = element_line(color = "black"),
    axis.text = element_text(size = 6),
    axis.title = element_text(size = 8),
    
    plot.title = element_text(size = 8),
    plot.margin = margin(b = 0, t = 20, unit = "pt"),
    
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 8),
    
  )

#Need this for plotting CAL line over other lines next
dat_CAL_sum <- dat_CAL %>%
  group_by(Sample, Assay, within_AMR) %>%
  summarize(
    exp = gmean(conc_exp),
    obs = gmean(conc_calc)
  )


######################## linear regression of VIG ######################

rm(list=setdiff(ls(), c("gmean", "gsd", "gcv", "gmed", "read_data", "dat_CAL_sum", "Fig_4A", "lines")))

dat <- read_csv("2025-Oct-15_MPXV_paper_validation_data.csv")
dat_assay_para <- read_csv("assay_limits/assay_para_validated.csv")

#filter only VIG linearity data
dat_VIG <- dat %>%
  filter(str_detect(Sample,"VIG-0[1-8]")) %>%
  mutate(
    day = case_when(
      str_detect(run_ID, "20240802") ~ "run 1",
      str_detect(run_ID, "20240806") ~ "run 2",
      TRUE ~ NA_character_
    ),
    .after = run_ID
  )


#Get VIG expected values, based on previous testing
VIG_exp <- read_csv("expected_value_tables/VIG_exp_vales.csv", n_max = 10) 

VIG_exp <- VIG_exp %>%
  mutate(
    `VIG-01` = exp_AU_mL/500,
    `VIG-02` = `VIG-01`/4,
    `VIG-03` = `VIG-02`/4,
    `VIG-04` = `VIG-03`/4,
    `VIG-05` = `VIG-04`/4,
    `VIG-06` = `VIG-05`/4,
    `VIG-07` = `VIG-06`/4,
    `VIG-08` = `VIG-07`/4
  ) %>%
  select(Assay, `VIG-01`, `VIG-02`, `VIG-03`, `VIG-04`, `VIG-05`, `VIG-06`, `VIG-07`, `VIG-08`) %>%
  pivot_longer(
    cols = starts_with("VIG"),
    names_to = "Sample",
    values_to = "exp"
  )

#Average replicates from each day
dat_VIG_sum <- dat_VIG %>%
  group_by(Sample, Assay, Dilution, day) %>%
  summarize(
    obs = gmean(conc_calc),
    obs_sd = gsd(conc_calc)
  ) 

#Join average data with expected values, calc fold diff between obs and exp
dat_VIG_sum <- inner_join(dat_VIG_sum, VIG_exp, by = c("Sample", "Assay")) %>%
  mutate(
    ratio = obs/exp
  )

dat_VIG_sum <- full_join(dat_VIG_sum, dat_assay_para %>% select(Assay, LoD, LLoQ_final, ULoQ_final), by ="Assay")

dat_VIG_sum <- dat_VIG_sum %>%
  mutate(
    LLOQ_flag = if_else(obs/Dilution < LLoQ_final, TRUE, FALSE),
    ULOQ_flag = if_else(obs/Dilution > ULoQ_final, TRUE, FALSE),
    within_AMR = if_else(LLOQ_flag == FALSE & ULOQ_flag == FALSE, TRUE, FALSE)
  )

models_VIG <- dat_VIG_sum %>%
  filter(within_AMR == TRUE) %>%
  filter(! Sample %in% c("VIG-08", "VIG-07")) %>%
  group_by(Assay) %>%
  nest() %>%
  mutate(
    model = map(data, ~lm(log10(obs) ~ log10(exp), data = .x)),
    model_tidy = map(model, ~broom::tidy(.x)),
    model_data = map(model, ~broom::augment(.x)),
    intercept = map_dbl(model_tidy, ~.x[[1,2]]),
    intercept_se = map_dbl(model_tidy, ~.x[[1,3]]),
    intercept_pvalue = map_dbl(model_tidy, ~.x[[1,5]]),
    slope = map_dbl(model_tidy, ~.x[[2,2]]),
    slope_se = map_dbl(model_tidy, ~.x[[2,3]]),
    R_sq = map_dbl(model, ~summary(.x)$r.squared)
  )

models_VIG_simple <- models_VIG %>%
  select(Assay, intercept, slope, R_sq) %>%
  mutate(
    intercept = if_else(abs(intercept) < 0.1, "< +/- 0.1 AU/mL", paste0(round(intercept,2), "AU/mL")),
    slope = formatC(slope, format = "f", digits = 2),
    R_sq = formatC(R_sq, format = "f", digits = 2),
    label = paste0("y-int.: ", intercept, "\nslope: ", slope, "\nR^2: ", R_sq),
  )

models_VIG_annot <- models_VIG_simple %>%
  select(Assay,label)

#Not used in supplementary material
a <- ggplot(dat_VIG_sum, aes(x = log2(exp), y = log2(obs), fill = day)) +
  facet_wrap(
    vars(Assay),
    nrow = 2,
    axes = "all_x"
  ) +
  geom_smooth(
    data = dat_CAL_sum %>% filter(within_AMR == TRUE),
    inherit.aes = FALSE,
    aes(x = log2(exp), y = log2(obs), linetype = "Calibrator\nregression"),
    method = lm,
    se = FALSE,
    #linetype = "dotted",
    linewidth = 0.4,
    color = "black",
    fullrange= TRUE
  ) +
  scale_linetype_manual(
    values = "dotted",
    name = NULL
  ) +
  scale_fill_manual(
    values = c("#009E73","#D55E00")
  ) +
  geom_point(
    shape = 21,
    alpha = 0.5
  ) +
  ylab("observed AU/mL") +
  xlab("expected AU/mL") +
  geom_hline(
    data = lines,
    aes(yintercept = log2(value), color = Limit),
    linetype = "dashed",
    linewidth = 0.4
  ) +
  scale_color_manual(
    values = c("LoD" = "red", "LLoQ" = "firebrick4", "ULoQ" = "black")
  ) +
  coord_cartesian(
    clip = "off",
    ylim = c(-10,10)
  ) +
  scale_y_continuous(
    guide = guide_axis(minor.ticks = TRUE),
    minor_breaks = c(-9, -8, -7, -5, -4, -3, -1, 0, 1, 3, 4, 5, 7, 8, 9),
    limits = c(-30,30),
    breaks = seq(-10,10,4),
    labels = label_math(expr = 2^.x)
  ) +
  scale_x_continuous(
    guide = guide_axis(minor.ticks = TRUE),
    minor_breaks = c(-9, -8, -7, -5, -4, -3, -1, 0, 1, 3, 4, 5, 7, 8, 9),
    limits = c(-10,10),
    breaks = seq(-10,10,4),
    labels = label_math(expr = 2^.x)
  ) +
  geom_text(
    inherit.aes = FALSE,
    data = models_VIG_annot,
    mapping = aes(x = -10, y = 18, label = label),
    size = 8/.pt,
    hjust = 0,
    color = "blue"
  ) +
  theme(
    aspect.ratio = 1,
    axis.text.x = element_text(angle = 90, vjust = 0.5),
    panel.grid = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(color = "black"),
    axis.text = element_text(size = 6),
    plot.title = element_text(size = 8),
    axis.title = element_text(size = 8),
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 8, hjust = 0),
    panel.spacing.y = unit(40, unit = "pt"),
    plot.margin = margin(b = 0, t = 20, unit = "pt")
  ) 


######################## linear regression of M32 (specimen) ######################

rm(list=setdiff(ls(), c("gmean", "gsd", "gcv", "gmed", "read_data", "dat_CAL_sum", "Fig_4A", "lines")))

dat <- read_csv("2025-Oct-15_MPXV_paper_validation_data.csv")
dat_assay_para <- read_csv("assay_limits/assay_para_validated.csv")

dat_M32 <- dat %>%
  filter(str_detect(Sample,"M32-0[1-8]")) %>%
  mutate(
    day = case_when(
      str_detect(run_ID, "20240802") ~ 1,
      str_detect(run_ID, "20240806") ~ 2,
      TRUE ~ NA_integer_
    ),
    .after = run_ID
  ) %>%
  filter(!is.na(conc_calc))


M32_exp <- read_csv("expected_value_tables/M32_exp_vales.csv") 

M32_exp <- M32_exp %>%
  mutate(
    `M32-01` = exp_AU_mL/50,
    `M32-02` = `M32-01`/4,
    `M32-03` = `M32-02`/4,
    `M32-04` = `M32-03`/4,
    `M32-05` = `M32-04`/4,
    `M32-06` = `M32-05`/4,
    `M32-07` = `M32-06`/4,
    `M32-08` = `M32-07`/4
  ) %>%
  select(Assay, `M32-01`, `M32-02`, `M32-03`, `M32-04`, `M32-05`, `M32-06`, `M32-07`, `M32-08`) %>%
  pivot_longer(
    cols = starts_with("M32"),
    names_to = "Sample",
    values_to = "exp"
  )


dat_M32_sum <- dat_M32 %>%
  group_by(Sample, Assay, day, Dilution) %>%
  summarize(
    obs = gmean(conc_calc),
    obs_sd = gsd(conc_calc),
    obs_signal = mean(Signal)
  ) 

dat_M32_sum <- inner_join(dat_M32_sum, M32_exp, by = c("Sample", "Assay")) %>%
  mutate(
    ratio = obs/exp
  )

rm(M32_exp)

dat_M32_sum <- full_join(dat_M32_sum, dat_assay_para %>% select(Assay, LoD, LLoQ_final, ULoQ_final), by ="Assay")

dat_M32_sum <- dat_M32_sum %>%
  mutate(
    LLOQ_flag = if_else(obs/Dilution < LLoQ_final, TRUE, FALSE),
    ULOQ_flag = if_else(obs/Dilution > ULoQ_final, TRUE, FALSE),
    within_AMR = if_else(LLOQ_flag == FALSE & ULOQ_flag == FALSE, TRUE, FALSE)
  )

models_M32 <- dat_M32_sum %>%
  filter(within_AMR == TRUE) %>%
  group_by(Assay) %>%
  nest() %>%
  mutate(
    model = map(data, ~lm(log10(obs) ~ log10(exp), data = .x)),
    model_tidy = map(model, ~broom::tidy(.x)),
    model_data = map(model, ~broom::augment(.x)),
    intercept = map_dbl(model_tidy, ~.x[[1,2]]),
    intercept_se = map_dbl(model_tidy, ~.x[[1,3]]),
    intercept_pvalue = map_dbl(model_tidy, ~.x[[1,5]]),
    slope = map_dbl(model_tidy, ~.x[[2,2]]),
    slope_se = map_dbl(model_tidy, ~.x[[2,3]]),
    R_sq = map_dbl(model, ~summary(.x)$r.squared)
  )

models_M32_simple <- models_M32 %>%
  select(Assay, intercept, slope, R_sq) %>%
  mutate(
    intercept = if_else(abs(intercept) < 0.1, "< +/- 0.1 AU/mL", paste0(round(intercept,2), "AU/mL")),
    slope = formatC(slope, format = "f", digits = 2),
    R_sq = formatC(R_sq, format = "f", digits = 2),
    label = paste0("y-int.: ", intercept, "\nslope: ", slope, "\nR^2: ", R_sq),
  )

models_M32_annot <- models_M32_simple %>%
  select(Assay,label)

dat_M32_sum$day <- factor(dat_M32_sum$day, levels = 1:10)

Fig_4B <- ggplot(dat_M32_sum, aes(x = log2(exp), y = log2(obs), fill = day)) +
  facet_wrap(
    vars(Assay),
    nrow = 2,
    axes = "all_x"
  ) +
  geom_point(
    shape = 21,
    size = 2,
    stroke = 0.15,
    alpha = 0.6,
    show.legend = FALSE #Going to use legend from Fig 4A - could not get patchwork to merge legneds
  ) +
  scale_fill_manual(
    name = "run day",
    values = carto_pal(10, "Safe"),
    drop = FALSE
  ) +
  ylab("observed AU/mL") +
  xlab("expected AU/mL") +
  geom_hline(
    data = lines,
    aes(yintercept = log2(value), color = Limit),
    linetype = "dashed",
    linewidth = 0.4
  ) +
  scale_color_manual(
    values = c("LoD" = "red", "LLoQ" = "firebrick4", "ULoQ" = "black")
  ) +
  geom_smooth(
    data = dat_CAL_sum %>% filter(within_AMR == TRUE),
    inherit.aes = FALSE,
    aes(x = log2(exp), y = log2(obs), linetype = "Calibrator\nregression"),
    method = lm,
    se = FALSE,
    linewidth = 0.4,
    color = "black",
    fullrange = TRUE
  ) +
  scale_linetype_manual(
    values = "dotted",
    name = NULL
  ) +
  coord_cartesian(
    clip = "off",
    ylim = c(-13,10)
  ) +
  scale_y_continuous(
    guide = guide_axis(minor.ticks = TRUE),
    minor_breaks = c(-11, -10, -9, -7, -6, -5, -3, -2, -1, 1, 2, 3, 5, 6, 7, 9, 10),
    limits = c(-30,30),
    breaks = seq(-12,10,4),
    labels = label_math(expr = 2^.x)
  ) +
  scale_x_continuous(
    guide = guide_axis(minor.ticks = TRUE),
    minor_breaks = c(-11, -10, -9, -7, -6, -5, -3, -2, -1, 1, 2, 3, 5, 6, 7, 9, 10),
    limits = c(-12,10),
    breaks = seq(-12,10,4),
    labels = label_math(expr = 2^.x)
  ) +
  geom_text(
    inherit.aes = FALSE,
    data = models_M32_annot,
    mapping = aes(x = -12, y = 21, label = label),
    size = 8/.pt,
    hjust = 0,
    color = "blue"
  ) +
  theme(
    aspect.ratio = 1,
    
    panel.grid = element_blank(),
    panel.background = element_blank(),
    panel.spacing.y = unit(40, unit = "pt"),
    
    axis.text.x = element_text(angle = 90, vjust = 0.5),
    axis.line = element_line(color = "black"),
    axis.text = element_text(size = 6),
    
    axis.title = element_text(size = 8),
    
    plot.title = element_text(size = 8),
    plot.margin = margin(b = 0, t = 20, unit = "pt"),
    
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 8),
  )

Fig_S4 <- Fig_4A / Fig_4B +
  plot_annotation(tag_levels = "A") +
  plot_layout(guides = "collect")

ggsave("fig_output/Fig_S4.pdf", width = 8, height = 9, units = "in", device = cairo_pdf)



