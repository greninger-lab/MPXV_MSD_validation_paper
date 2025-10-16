# requires running LoD/LLoQ/ULoQ section to get assay limits to make 
# assay_para_validated.csv file

# Fig S5A - Results - samples
# Fig S5B - Results - Control 1,2,3
# Fig S6A - Precision - samples
# Fig S6B - Precision - Control 1,2,3
# Precision - Control 1, 2, 3

# Greninger Lab - Jon Reed
# 2025-Oct-15

source("../main/MPXV_paper_functions.r")

package_list <- c("tidyverse",
                  "ggplot2",
                  "ggbeeswarm",
                  "rcartocolor",
                  "scales",
                  "ggrepel",
                  "VCA",
                  "patchwork")

invisible(lapply(package_list, check_and_load))

dat <- read_csv("2025-Oct-15_MPXV_paper_validation_data.csv")
dat_assay_para <- read_csv("assay_limits/assay_para_validated.csv")

precision_Sample <- tribble(
  ~Sample, ~status,
  "M4", "Vaccine",
  "M31", "MPXV",
  "M37", "Vaccine",
  "M54", "Vaccine",
  "M60", "MPXV"
)

dat_p <- dat %>%
  filter(
    Sample %in% paste0(precision_Sample$Sample, "_rep_1") |
      Sample %in% paste0(precision_Sample$Sample, "_rep_2")
  ) %>%
  filter(!is.na(day)) %>%
  mutate(
    rep = case_when(
      str_detect(Sample,"rep_1") ~ 1,
      str_detect(Sample,"rep_2") ~ 2,
      TRUE ~ NA_integer_
    ),
    Sample_alt = str_remove(Sample, "_rep_1|_rep_2"),
    .after = Sample
  ) %>%
  mutate(
    conc_calc_dil = conc_calc/Dilution
  ) %>%
  mutate(
    Sample_Assay = paste0(Sample_alt, " - ", Assay),
    .after = Assay
  ) %>%
  mutate(
    day = factor(
      day, 
      levels = c(1,2,3,4,5,6,7,8,9,10),
      labels = c("day 1","day 2","day 3","day 4","day 5",
                 "day 6", "day 7","day 8","day 9","day 10")
    )
  ) %>%
  # Get Vaccine or MPXV inf for each specimen
  full_join(precision_Sample, by = c("Sample_alt" = "Sample")) %>%
  # Set Sample name for plotting
  mutate(
    Sample_alt_2 = case_when(
      status == "Vaccine" ~ str_replace(Sample_alt,"M","V"),
      status == "MPXV" ~ Sample_alt,
      TRUE ~ NA_character_
    ),
    .after = Sample_alt
  ) %>%
  #Join with dat_assay_para to get LoD, LLoQ, and ULoQ
  full_join(dat_assay_para, by = "Assay") %>%
  mutate(
    LoD_flag = if_else(Signal < LoD, "below LoD", "above LoD"),
    LLoQ_flag = if_else(conc_calc_dil < LLoQ_final, "below LLoQ", "above LLoQ"),
    ULoQ_flag = if_else(conc_calc_dil > ULoQ_final, "above ULoQ", "below ULLoQ")
  ) %>%
  mutate(
    day_rep = factor(
      paste0(day, ", rep: ", rep), 
      levels = c("day 1, rep: 1", "day 1, rep: 2",
                 "day 2, rep: 1", "day 2, rep: 2",
                 "day 4, rep: 1", "day 4, rep: 2",
                 "day 8, rep: 1", "day 8, rep: 2",
                 "day 9, rep: 1", "day 9, rep: 2",
                 "day 10, rep: 1", "day 10, rep: 2")
    )
  )



dat_p_sum <- dat_p %>%
  group_by(Assay, Sample_alt_2, day, rep, day_rep, LoD, LLoQ_final, ULoQ_final) %>%
  summarize(
    gmean_conc = gmean(conc_calc),
    gmean_signal = gmean(Signal),
    gcv_conc = gcv(conc_calc),
    max = max(conc_calc),
    min = min(conc_calc),
    fold = max/min
  ) 

# filter out results where GCV between wells was > 0.37

high_interwell_GCV_filter <- dat_p_sum %>%
  ungroup() %>%
  filter(gcv_conc > 0.37)

dat_p_filtered <- dat_p %>%
  anti_join(high_interwell_GCV_filter, by = c("Assay", "Sample_alt_2", "day", "rep"))

dat_p_sum_filtered <- dat_p_sum %>%
  anti_join(high_interwell_GCV_filter, by = c("Assay", "Sample_alt_2", "day", "rep"))

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
    ),
    #For precision the specimens were diluted 1/500, so adjust limits accordingly
    value_dil_adj = value*500
  )

# average of two wells
Fig_S5A <- ggplot(dat_p_sum_filtered, aes(y = log2(gmean_conc), x = Sample_alt_2, fill = day_rep)) +
  facet_wrap(
    vars(Assay),
    nrow = 2,
    axes = "all_x"
  ) +
  geom_quasirandom(
    shape = 21,
    size = 1.75,
    varwidth = TRUE,
    alpha = 0.8,
    stroke = 0
  ) +
  geom_hline(
    data = lines,
    aes(yintercept = log2(value_dil_adj), color = Limit),
    linetype = "dashed",
    linewidth = 0.4
  ) +
  scale_color_manual(
    values = c("LoD" = "red", "LLoQ" = "firebrick4", "ULoQ" = "black")
  ) +
  scale_fill_manual(
    name = "Test (day and replicate)",
    values = carto_pal(12, "Safe")
  ) +
  scale_y_continuous(
    limits = c(-2,16),
    breaks = seq(-2,20,2),
    labels = label_math(expr = 2^.x),
    expand = c(0,0),
    guide = guide_axis(minor.ticks = TRUE),
    minor_breaks = seq(-1,15,2)
  ) +
  ylab("AU/mL") +
  xlab("Specimen - MPXV-infected (M) or smallpox vaccinated (V)") +
  guides(
    fill = guide_legend(order = 1),
    color = guide_legend(order = 2)
  ) +
  theme(
    panel.grid = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(color = "black"),
    axis.text = element_text(size = 7),
    axis.title = element_text(size = 8),
    
    legend.title = element_text(hjust = 0, size = 8),
    legend.text = element_text(size = 7),
    legend.key.height = unit(10, "pt")
  )

get_anova_var <- function(dat_sub, group_key, var) {
  
  dat_sub <- dat_sub %>%
    mutate(
      value = log2({{ var }})
    )
  
  fit <- anovaVCA(value ~ day/rep,
                  Data = as.data.frame(dat_sub),
                  NegVC = TRUE
  )
  
  tibble(
    mean = 2^fit$Mean,
    sd_inter_well = 2^fit$aov.tab[4,6],
    sd_intra = 2^fit$aov.tab[3,6],
    sd_inter = 2^fit$aov.tab[2,6],
    sd_lab = 2^fit$aov.tab[1,6],
    VCA = list(fit)
  )
}

dat_p_anova <- dat_p_filtered %>%
  group_by(Sample_alt_2, Assay, Sample_Assay) %>%
  group_modify(~get_anova_var(.x,.y, conc_calc)) %>%
  mutate(
    cv_inter_well = sqrt(exp(log(sd_inter_well)^2) - 1),
    cv_intra = sqrt(exp(log(sd_intra)^2) - 1),
    cv_inter = sqrt(exp(log(sd_inter)^2) - 1),
    cv_lab = sqrt(exp(log(sd_lab)^2) - 1)
  )

Fig_S6A <- ggplot(dat_p_anova, aes(x = cv_intra, y = cv_lab, label = Sample_alt_2, fill = cv_inter)) +
  geom_point(
    shape = 21,
  ) +
  facet_wrap(
    vars(Assay),
    nrow = 2,
    axes = "all_x"
  ) +
  geom_label_repel(
    box.padding = 0.1,
    label.padding = 0.1,
    size = 2,
    max.overlaps = Inf,
    force = 10
  ) +
  scale_fill_continuous(
    type = "viridis",
    alpha = 0.5,
    label = percent,
    limits = c(0,0.4),
    name = str_wrap("inter-assay geometric CV (%)", width = 20)
  ) +
  geom_hline(
    yintercept = 0.37,
    linetype = "dashed"
  ) +
  scale_y_continuous(
    breaks = seq(0,0.6,0.1),
    limits = c(0,0.6),
    expand = c(0,0),
    label = percent
  ) +
  scale_x_continuous(
    breaks = seq(0,0.4,0.1),
    limits = c(0,0.4),
    expand = c(0.05,0.05),
    label = percent
  ) +
  ylab("within-lab geometric CV (%)") +
  xlab("intraassay, geometric CV (%)") +
  theme(
    panel.grid = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(color = "black"),
    axis.text = element_text(size = 7),
    axis.title = element_text(size = 8),
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 8)
  )


####################### Precision of Controls ###############################

rm(list=setdiff(ls(), c("gmean", "gsd", "gcv", "gmed", "read_data", "lines", "Fig_S5A", "Fig_S6A")))

dat <- read_csv("2025-Oct-15_MPXV_paper_validation_data.csv")
dat_assay_para <- read_csv("assay_limits/assay_para_validated.csv")

precision_Sample <- c("Control 1", "Control 2", "Control 3")

dat_p <- dat %>%
  filter(Sample %in% precision_Sample) %>%
  filter(!is.na(day)) %>%
  mutate(
    rep = 1,
    .after = Sample
  ) %>%
  mutate(
    conc_calc_dil = conc_calc/Dilution
  ) %>%
  mutate(
    Sample_Assay = paste0(Sample, " - ", Assay),
    .after = Assay
  ) %>%
  mutate(
    day = factor(
      day,
      levels = c(1,2,3,4,5,6,7,8,9,10),
      labels = c("day 1","day 2","day 3","day 4","day 5","day 6",
                 "day 7","day 8","day 9","day 10")
    )
  )


#Join with dat_assay_para to get LoD and LLoQ

dat_p <- full_join(dat_p, dat_assay_para, by = "Assay")

dat_p <- dat_p %>%
  mutate(
    LoD_flag = if_else(Signal < LoD, "below LoD", "above LoD"),
    LLoQ_flag = if_else(conc_calc_dil < LLoQ_final, "below LLoQ", "above LLoQ"),
    ULoQ_flag = if_else(conc_calc_dil > ULoQ_final, "above ULoQ", "below ULLoQ")
  )

dat_p_sum <- dat_p %>%
  group_by(Assay, Sample, day, LoD, LLoQ_final, ULoQ_final, conc_exp) %>%
  summarize(
    gmean_conc = gmean(conc_calc),
    gmean_signal = gmean(Signal),
    count = n()
  ) 


Fig_S5B <- ggplot(dat_p_sum, aes(y = log2(gmean_conc), x = Sample, fill = day)) +
  facet_wrap(
    vars(Assay),
    nrow = 2,
    axes = "all_x"
  ) +
  geom_quasirandom(
    shape = 21,
    size = 1.75,
    varwidth = TRUE,
    alpha = 0.8,
    stroke = 0
  ) +
  geom_hline(
    data = lines,
    aes(yintercept = log2(value), color = Limit),
    linetype = "dashed",
    linewidth = 0.4,
    show.legend = FALSE
  ) +
  scale_color_manual(
    values = c("LoD" = "red", "LLoQ" = "firebrick4", "ULoQ" = "black")
  ) +
  scale_fill_manual(
    name = "Test (day only)",
    values = carto_pal(12, "Safe")
  ) +
  scale_y_continuous(
    limits = c(-12,8),
    breaks = seq(-12,8,2),
    labels = label_math(expr = 2^.x),
    expand = c(0,0),
    guide = guide_axis(minor.ticks = TRUE),
    minor_breaks = seq(-11,7,2)
  ) +
  ylab("AU/mL") +
  xlab("Control specimen") +
  theme(
    panel.grid = element_blank(),
    panel.background = element_blank(),
    
    axis.line = element_line(color = "black"),
    axis.text.x = element_text(size = 6),
    axis.title = element_text(size = 8),
    axis.text.y = element_text(size = 7),
    
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 7),
    legend.key.height = unit(10, "pt")
  )

get_anova_var <- function(dat_sub, group_key, var) {
  
  dat_sub <- dat_sub %>%
    mutate(
      value = log2({{ var }})
    )
  
  fit <- anovaVCA(value ~ day,
                  Data = as.data.frame(dat_sub),
                  NegVC = TRUE)
  
  tibble(
    mean = 2^fit$Mean,
    sd_intra = 2^fit$aov.tab[3,6],
    sd_inter = 2^fit$aov.tab[2,6],
    sd_lab = 2^fit$aov.tab[1,6],
    VCA = list(fit)
  )
}

dat_p_anova <- dat_p %>%
  group_by(Sample, Assay, Sample_Assay) %>%
  group_modify(~get_anova_var(.x,.y, conc_calc)) %>%
  mutate(
    cv_intra = sqrt(exp(log(sd_intra)^2) - 1),
    cv_inter = sqrt(exp(log(sd_inter)^2) - 1),
    cv_lab = sqrt(exp(log(sd_lab)^2) - 1)
  )

Fig_S6B <- ggplot(dat_p_anova, aes(x = cv_intra, y = cv_lab, label = Sample, fill = cv_inter)) +
  geom_point(
    shape = 21,
  ) +
  facet_wrap(
    vars(Assay),
    nrow = 2,
    axes = "all_x"
  ) +
  geom_label_repel(
    box.padding = 0.1,
    label.padding = 0.1,
    size = 2,
    max.time = 10
  ) +
  scale_fill_continuous(
    type = "viridis",
    alpha = 0.5,
    name = str_wrap("inter-assay geometric CV (%)", width = 20),
    limits = c(0,0.4),
    label = percent
  ) +
  geom_hline(
    yintercept = 0.37,
    linetype = "dashed"
  ) +
  scale_y_continuous(
    breaks = seq(0,0.4,0.1),
    limits = c(0,0.4),
    expand = c(0,0),
    label = percent
  ) +
  scale_x_continuous(
    breaks = seq(0,0.3,0.1),
    limits = c(0,0.3),
    expand = c(0.05,0.05),
    label = percent
  ) +
  ylab("within-lab geometric CV (%)") +
  xlab("intra-assay, between-well, geometric CV (%)") +
  theme(
    panel.grid = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(color = "black"),
    axis.text = element_text(size = 7),
    axis.title = element_text(size = 8),
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 8)
  )

Fig_S5 <- Fig_S5A / Fig_S5B +
  plot_annotation(tag_levels = "A")

ggsave("fig_output/Fig_S5.pdf", width = 8, height = 7, units = "in", device = cairo_pdf)


Fig_S6 <- Fig_S6A / Fig_S6B +
  plot_annotation(tag_levels = "A") +
  plot_layout(guides = "collect")

ggsave("fig_output/Fig_S6.pdf", width = 8, height = 7, units = "in", device = cairo_pdf)

