# Fig S1 (A and B) - CAL Signal and AU/mL, LoD
# Fig S3 (A and B) - Precision
# Table S1 - assay parameters
# Table S2 - 4PL fits

# LOD - based on AU/mL
# Verification of calibration model
# ULoQ and LLoQ

# Greninger Lab - Jon Reed
# 2025-Oct-15

source("../main/MPXV_paper_functions.r")

package_list <- c("tidyverse",
                  "ggplot2",
                  "scales",
                  "ggbeeswarm",
                  "patchwork",
                  "VCA",
                  "openxlsx2")

invisible(lapply(package_list, check_and_load))


# create table of manf. LLoQ and ULoQ values and expected AU/mL for CAL-01
dat_assay_para <- tribble(~Assay, ~manf_A0080357, ~manf_LLoQ, ~manf_ULoQ,
                          "VACV A27L", 25.8, 0.028, 25.8,
                          "VACV A33R", 89.8, 0.026, 89.8,
                          "VACV B5R", 18.9, 0.025, 18.9,
                          "VACV D8L", 17.6, 0.028, 17.6,
                          "VACV L1R", 9.83, 0.039, 9.83,
                          "MPXV M1R", 8.56, 0.022, 8.56,
                          "MPXV E8L", 11.8, 0.025, 11.8,
                          "MPXV B6R", 14.2, 0.058, 14.2,
                          "MPXV A35R", 60.1, 0.021, 60.1,
                          "MPXV A29L", 20.4, 0.033, 20.4)


dat <- read_csv("2025-Oct-15_MPXV_paper_validation_data.csv")

filter_lod <- tribble(~Assay, ~Sample,
                      "MPXV A29L", "CAL-08",
                      "MPXV A35R", "CAL-08",
                      "MPXV B6R", "CAL-08",
                      "MPXV E8L", "CAL-08",
                      "MPXV M1R", "CAL-08",
                      "VACV A27L", "CAL-08",
                      "VACV A33R", "CAL-08",
                      "VACV B5R", "CAL-08",
                      "VACV D8L", "CAL-08",
                      "VACV L1R", "CAL-08",
                      "MPXV A29L", "CAL-07",
                      "MPXV A35R", "CAL-07",
                      "MPXV B6R", "CAL-07",
                      "MPXV E8L", "CAL-07",
                      "MPXV M1R", "CAL-07",
                      "VACV A27L", "CAL-07",
                      "VACV A33R", "CAL-07",
                      "VACV B5R", "CAL-07",
                      "VACV D8L", "CAL-07",
                      "VACV L1R", "CAL-07")

#To calculate AU/mL
#y_pred is assay Signal
#x is AU/mL
# b1 <- upper
# b2 <- lower
# b3 <- midpoint
# b4 <- slope

#y_pred = b1 +  ((b2-b1)/(1+(x/b3)^b4))
#x = b3*(((b2-b1)/(y-b1)) - 1)^(1/b4)

dat_LOD <- dat %>%
  filter(str_detect(Sample, "CAL-0[1-8]")) %>% #Pull out calibrator data
  mutate(
    #for Signal below the cal curve fit, calc. minimum AU/mL
    Signal = if_else( `Detection Range` == "Below Fit Curve Range", b2+1, Signal),
    conc_calc_alt = if_else(
      `Detection Range` == "Below Fit Curve Range",
      b3*(((b2-b1)/(Signal-b1)) - 1)^(1/b4),
      conc_calc
    ),
    .after = conc_calc
  )

#Save all the calibrator data
dat_CAL <- dat_LOD

#calculate geometric mean of LoD as determined by manf.
dat_LOD_manf <- dat_LOD %>%
  group_by(Assay) %>%
  summarize(
    manf_LoD = gmean(`Detection Limits: Calc. Low`),
  )

#join LoD estimates with assay parameter table
dat_assay_para <- full_join(dat_assay_para, dat_LOD_manf, by = "Assay") %>%
  select(Assay, manf_LoD, manf_LLoQ, manf_ULoQ)

CAL_08_stats <- dat_LOD %>%
  filter(Sample == "CAL-08") %>%
  group_by(Assay) %>%
  summarize(
    days = n()/2,
    signal_mean = round(mean(Signal),2),
    signal_sd = round(sd(Signal),2),
    CI_95_upper = round(signal_mean + 1.96*signal_sd,2)
  )

#Number of assays for display on graph
number_of_assays <- nrow(dat_CAL %>% select(run_ID) %>% unique())

Fig_S1A <- ggplot(dat_CAL, aes(x = "CAL-1 thru CAL-08", y = log2(Signal), fill = Sample)) +
  facet_wrap(
    vars(Assay),
    axes = "all_x",
    nrow = 2
  ) +
  geom_beeswarm(
    shape = 21,
    size = 1,
    cex = 4
  ) +
  scale_y_continuous(
    label = label_math(expr = 2^.x),
    guide = guide_axis(minor.ticks = TRUE),
    limits = c(4,21),
    breaks = seq(4,21,2),
    expand = c(0,0)
  ) +
  ylab(expression(Signal)) +
  ggtitle(paste0("number of assays analyzed: ",number_of_assays)) +
  theme(
    panel.grid = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(color = "black"),
    axis.title.x = element_blank(),
    axis.text = element_text(size = 6),
    plot.title = element_text(size = 8),
    axis.title.y = element_text(size = 8),
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 8, hjust = 0)
  )

Fig_S1B <- ggplot(dat_CAL, aes(x = "CAL-1 thru CAL-08", y = log2(conc_calc_alt), fill = Sample)) +
  facet_wrap(
    vars(Assay),
    axes = "all_x",
    nrow = 2
  ) +
  geom_beeswarm(
    shape = 21,
    size = 1,
    cex = 5
  ) +
  geom_hline(
    data = dat_LOD_manf,
    aes(yintercept = log2(manf_LoD), color = "LoD"),
    linetype = "dashed",
  ) +
  scale_color_manual(
    values = "red",
    name = NULL
  ) +
  scale_y_continuous(
    label = label_math(expr = 2^.x),
    guide = guide_axis(minor.ticks = TRUE),
    limits = c(-18,8),
    breaks = seq(-18,8,2)
  ) +
  ylab("Antibody Level (AU/mL)") +
  ggtitle(paste0("number of assays analyzed: ",number_of_assays)) +
  theme(
    panel.grid = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(color = "black"),
    axis.title.x = element_blank(),
    axis.text = element_text(size = 6),
    plot.title = element_text(size = 8),
    axis.title.y = element_text(size = 8),
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 8, hjust = 0)
    
  )

AB <- (Fig_S1A / Fig_S1B) + 
  plot_annotation(tag_levels = "A") +
  plot_layout(guides = "collect")

ggsave("fig_output/Fig_S1.pdf", width = 7, height = 8, dpi = 600, units = "in", device = cairo_pdf)


# Assign final LoD in assay parameter table
dat_assay_para <- dat_assay_para %>%
  mutate(
    LoD = 0.004
  )

# cleanup 
rm(dat_CAL, dat_LOD, filter_lod, number_of_assays, dat_LOD_manf, CAL_08_stats)

# setup for next section

cal_model_process_dat <- function(x) {
  
  x <- x %>%
    filter(str_detect(Sample, "CAL-0[1-7]")) %>%
    filter(!is.na(day)) %>%
    mutate(
      rep = NA_integer_
    ) %>%
    nest(.by = c(Sample, Assay, day)) %>%
    mutate(
      data = map(data, function(x) {
        x$rep <- c(1,2)
        x
      }),
    ) %>%
    unnest(data)
  
  return(x)
}

dat_p <- cal_model_process_dat(dat)

############################ % Relative Error #################################

# Obtain mean value from each day (two technical replicates)
# Calc %RE for each CAL from each day of the assay

dat_p_sum <- dat_p %>%
  group_by(Sample, Assay, day, conc_exp) %>%
  summarize(
    mean = gmean(conc_calc),
    count = n()
  ) %>%
  mutate(
    RE = ((mean - conc_exp)/conc_exp)
  )

#Obtain mean of %RE across 10 days of testing

dat_p_sum_avg <- dat_p_sum %>%
  group_by(Assay, Sample) %>%
  summarize(
    mean = mean(RE)
  )

# plot %RE for each day of testing

# make day a factor and order, for plotting
dat_p$day <- factor(dat_p$day, levels = c(1,2,3,4,5,6,7,8,9,10),
                    labels = c("day 1","day 2","day 3","day 4","day 5",
                               "day 6","day 7","day 8","day 9","day 10")
)

# This figure was not included in supplemental
a <- ggplot(dat_p_sum, aes(x = Sample, y = RE, fill = Assay)) +
  geom_beeswarm(
    shape = 21,
    size = 1.5,
    alpha = 0.6,
    cex = 1.5
  ) +
  facet_wrap(
    vars(day),
    nrow = 2,
    axes = "all_x",
    labeller = labeller(day = function(x) paste("day: ", x))
  ) +
  geom_hline(
    yintercept = -.30,
    linetype = "dashed"
  ) +
  geom_hline(
    yintercept = .30,
    linetype = "dashed"
  ) +
  scale_y_continuous(
    limits = c(-.60,.60),
    label = percent
  ) +
  ylab("relative error (%)")
  
  

#plot %RE mean across ten days of testing

number_of_assays <- dat_p %>% select(run_ID) %>% unique() %>% nrow()

Fig_S3A <- ggplot(dat_p_sum_avg, aes(x = Assay, y = mean, fill = Sample)) +
  geom_beeswarm(
    shape = 21,
    size = 2,
    alpha = 1,
    cex = 1.5
  ) +
  geom_hline(
    yintercept = -.15,
    linetype = "dashed"
  ) +
  geom_hline(
    yintercept = .15,
    linetype = "dashed"
  ) +
  scale_y_continuous(
    limits = c(-.20,.20),
    breaks = seq(-0.2,0.2, 0.05),
    label = percent
  ) +
  ylab("Mean Relative Error (%)") +
  ggtitle(paste0("number of assays: ", number_of_assays)) +
  theme(
    panel.grid = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(color = "black"),
    axis.text.x = element_text(angle = 90, size = 6, vjust = 0.5),
    axis.text.y = element_text(size = 8),
    axis.title.y = element_text(size = 8),
    plot.title = element_text(size = 8),
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 8)
  )


################################## Precision ###################################

get_anova_var <- function(dat_sub, group_key, var) {
  
  dat_sub <- dat_sub %>%
    mutate(
      value = log2({{var}})
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

dat_anova <- dat_p %>%
  group_by(Sample, Assay) %>%
  group_modify(~get_anova_var(.x,.y, conc_calc))


# calculate gcv from sd estimate fron ANOVA
dat_anova <- dat_anova %>%
  mutate(
    gcv_lab = sqrt(exp(log(sd_lab)^2)-1)
  )

# get mean Signal value for each CAL for each antigen
temp <- dat_p %>%
  group_by(Sample, Assay) %>%
  summarize(
    mean_signal = mean(Signal)
  )

# combine with anova stats
dat_anova <- full_join(dat_anova, temp, by = c("Sample", "Assay"))

rm(temp)

number_of_assays <- nrow(dat_p %>% select(run_ID) %>% unique())

# plot of standard deviation, but not using this for figure
a <- ggplot(dat_anova, aes(x = log2(mean), y = sd_lab, fill = Sample)) +
  geom_point(
    shape = 21,
    size = 2.5
  ) +
  facet_wrap(
    vars(Assay),
    nrow = 2,
    axes = "all_x"
  ) +
  geom_hline(
    yintercept = 1.14,
    linetype = "dashed"
  ) +
  scale_x_continuous(
    limits = c(-15,8),
    labels = label_math(expr = 2^.x)
  ) +
  scale_y_continuous(
    limits = c(1,1.7),
    breaks = seq(1,1.6,0.2),
  ) +
  ylab("geometric standard deviation") +
  theme(
    panel.grid = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(color = "black")
  ) +
  ggtitle("Within-lab geometric standard deviation")

# plot of within-lab GCV
Fig_S3B <- ggplot(dat_anova, aes(x = log2(mean), y = gcv_lab, fill = Sample)) +
  geom_point(
    shape = 21,
    size = 2
  ) +
  facet_wrap(
    vars(Assay),
    nrow = 2,
    axes = "all_x"
  ) +
  geom_hline(
    yintercept = 0.15,
    linetype = "dashed"
  ) +
  scale_x_continuous(
    limits = c(-15,8),
    labels = label_math(expr = 2^.x)
  ) +
  scale_y_continuous(
    label = percent
  ) +
  ylab("Within-lab Geometric Coefficient of Variation") +
  xlab(expression(Mean~Antibody~Level~(log[2]~AU/mL))) +
  ggtitle(paste0("number of assays analyzed: ",number_of_assays)) +
  theme(
    panel.grid = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(color = "black"),
    axis.text = element_text(size = 6),
    plot.title = element_text(size = 8),
    axis.title = element_text(size = 8),
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 8)
  )

Fig_S3 <- Fig_S3A / Fig_S3B + 
  plot_annotation(tag_levels = "A") +
  plot_layout(guides = "collect")

ggsave("fig_output/Fig_S3.pdf", width = 7, height = 8, dpi = 600, units = "in", device = cairo_pdf)



################################## LLoQ ###################################

#get LLoQ values

dat_LLoQ <- full_join(dat_anova, dat_assay_para, by = "Assay") %>%
  filter(gcv_lab < 0.15) %>%
  filter(mean_signal > LoD) %>%
  nest(data = !Assay) %>%
  mutate(
    LLoQ_cal = map_chr(data, ~ .x[.x$mean == min(.x$mean),]$Sample),
    LLoQ = map_dbl(data, ~ .x[.x$mean == min(.x$mean),]$mean),
    ULoQ_cal = map_chr(data, ~ .x[.x$mean == max(.x$mean),]$Sample),
    ULoQ = map_dbl(data, ~ .x[.x$mean == max(.x$mean),]$mean),
  ) %>%
  select(!data)


#join empiric LLoQ with dat_cal table

dat_assay_para <- inner_join(dat_assay_para, dat_LLoQ, by = "Assay")

#LLoQ is the highest of our empiric estimate or LLoQ provided by manf.

dat_assay_para <- dat_assay_para %>%
  mutate(
    LLoQ_final = if_else(manf_LLoQ > LLoQ, manf_LLoQ, LLoQ),
    ULoQ_final = ULoQ,
    assays = number_of_assays
  )

# Create workbook, add data
wb <- wb_workbook() %>% 
  wb_add_worksheet(sheet = "Table S1") %>%
  wb_add_data(x = dat_assay_para) %>%
  # format manf LoD and final LoD
  wb_add_numfmt(
    sheet = "Table S1",
    dims = wb_dims(
      rows = 2:(nrow(dat_assay_para)+1),
      cols = c(2,5)
    ),
    numfmt = "0.0000"
  ) %>%
  # format manf ULoQ, empiric ULoQ, final ULoQ
  wb_add_numfmt(
    sheet = "Table S1",
    dims = wb_dims(
      rows = 2:(nrow(dat_assay_para)+1),
      cols = c(4,9,11)
    ),
    numfmt = "0.0"
  ) %>%
  # format manf LLoQ, empiric LLoQ, LLoQ final
  wb_add_numfmt(
    sheet = "Table S1",
    dims = wb_dims(
      rows = 2:(nrow(dat_assay_para)+1),
      cols = c(3,7,10)
    ),
    numfmt = "0.000"
  ) %>%
  #Auto adjust column widths
  wb_set_col_widths(
    sheet = "Table S1",
    cols = c(1:12),
    widths = "auto"
  ) %>%
  #Center columns
  wb_add_cell_style(
    sheet = "Table S1",
    dims = wb_dims(
      rows = 1:(nrow(dat_p_sum)+1),
      cols = c(1:12)
    ),
    vertical = "center",
    horizontal = "center"
  ) 

write_csv(dat_assay_para, "assay_limits/assay_para_validated.csv")

rm(dat_p_sum, dat_p_sum_avg, number_of_assays, dat, dat_anova, dat_LLoQ, dat_assay_para)



################################### 4PL stats ##################################


# b1 <- upper
# b2 <- lower
# b3 <- midpoint
# b4 <- slope

#y_pred = b1 +  ((b2-b1)/(1+(x/b3)^b4))
#x = b3*(((b2-b1)/(y-b1)) - 1)^(1/b4)

dat_p_sum <- dat_p %>%
  group_by(Assay, day) %>%
  summarize(
    min = unique(b2),
    max = unique(b1),
    hill = unique(b4),
    R = unique(`Fit Statistic: RSquared`),
    count = n()
  ) 

dat_p_sum <- dat_p %>%
  group_by(Assay) %>%
  summarize(
    min_min = round(min(b2),1),
    max_min = round(max(b2),1),
    min_max = round(min(b1),1),
    max_max = round(max(b1),1),
    min_hill = round(min(b4),2),
    max_hill = round(max(b4),2),
    min_R = round(min(`Fit Statistic: RSquared`),4),
    max_R = round(max(`Fit Statistic: RSquared`),4),
    days = length(unique(day))
  ) 

# Create workbook, add data
wb <- wb %>% 
  wb_add_worksheet(sheet = "Table S2") %>%
  wb_add_data(x = dat_p_sum) %>%
  # format low asy values
  wb_add_numfmt(
    sheet = "Table S2",
    dims = wb_dims(
      rows = 2:(nrow(dat_p_sum)+1),
      cols = 2:3
    ),
    numfmt = "0.0"
  ) %>%
  # format upper asy values
  wb_add_numfmt(
    sheet = "Table S2",
    dims = wb_dims(
      rows = 2:(nrow(dat_p_sum)+1),
      cols = 4:5
    ),
    numfmt = "0.00E+00"
  ) %>%
  #format hill slopes
  wb_add_numfmt(
    sheet = "Table S2",
    dims = wb_dims(
      rows = 2:(nrow(dat_p_sum)+1),
      cols = 6:7
    ),
    numfmt = "0.00"
  ) %>%
  #format R^2
  wb_add_numfmt(
    sheet = "Table S2",
    dims = wb_dims(
      rows = 2:(nrow(dat_p_sum)+1),
      cols = 8:9
    ),
    numfmt = "0.000"
  ) %>%
  #Auto adjust column widths
  wb_set_col_widths(
    sheet = "Table S2",
    cols = c(1:10),
    widths = "auto"
  ) %>%
  #Center columns
  wb_add_cell_style(
    sheet = "Table S2",
    dims = wb_dims(
      rows = 1:(nrow(dat_p_sum)+1),
      cols = c(1:10)
    ),
    vertical = "center",
    horizontal = "center"
  ) %>%
  wb_save("fig_output/Table_S1_S2.xlsx")










