# Fig S2 - Review of CAL-08 high Signal

# Greninger Lab - Jon Reed
# 2025-Oct-15


source("../main/MPXV_paper_functions.r")


package_list <- c("tidyverse",
                  "ggplot2",
                  "scales")

invisible(lapply(package_list, check_and_load))

dat <- read_csv("2025_Aug_26_all_CAL_Control_data.csv") %>%
  filter(!str_detect(files, "20250827")) #exclude this run outside of testing range

#Check on any runs that failed QC

QC <- read_csv("2025_Aug_26_all_QC.csv")

#get all runs with a failure.
QC <- QC %>%
  filter(`final QC` == "FAIL")

################################################################################

#These failed for reasons outside of high CAL-08, so will be excluded.
exclude <- c("20241015_kadenmc_V-PLEX_Orthopox_Plate2", #Assay failed CAL GSD
             "20241101_kadenmc_V-PLEX_Orthopox_Plate1" #high blank signal in E8L and D8L, which is unusual 
)

# remove excluded runs
QC <- QC %>%
  filter(!run_ID %in% exclude)

# Manually confirmed remaining runs only failed due to high CAL-08 signal
# Generate list for markup in plot
high_blank_list <- QC %>%
  select(run_ID, Assay) %>%
  mutate(
    plate_num = case_when(
      str_detect(run_ID, "Plate1") ~ 1,
      str_detect(run_ID, "Plate2") ~ 2,
      TRUE ~ 1
    ),
    run_date = ymd(str_extract(run_ID, "^[^_]+")),
    operator = if_else(str_detect(run_ID, "cfdowns"), "cfdowns", "kadenmc"),
    blank_low = "no"
  )

#List of validation runs to get validation set
validation_runs <- tribble(
  ~run_num, ~run_ID,
  1, "20240730_cfdowns_V-PLEX_Orthopox",
  2, "20240731_cfdowns_V-PLEX_Orthopox",
  3, "20240801_cfdowns_V-PLEX_Orthopox",
  4, "20240802_cfdowns_V-PLEX_Orthopox",
  5, "20240806_cfdowns_V-PLEX_Orthopox",
  6, "20240814_cfdowns_V-PLEX_Orthopox",
  7, "20240820_cfdowns_V-PLEX_Orthopox",
  8, "20240827_cfdowns_V-PLEX_Orthopox",
  9, "20240829_cfdowns_V-PLEX_Orthopox",
  10, "20240830_cfdowns_V-PLEX_Orthopox"
)

cal01_08 <- dat %>%
  filter(Sample %in% c("CAL-01","CAL-08")) %>%
  mutate(
    run_ID = str_extract(files, "[^/]+(?=\\.csv$)"),
    plate_num = case_when(
      str_detect(run_ID, "Plate1") ~ 1,
      str_detect(run_ID, "Plate2") ~ 2,
      TRUE ~ 1
    ),
    run_date = ymd(str_extract(run_ID, "^[^_]+")),
    operator = if_else(str_detect(run_ID, "cfdowns"), "cfdowns", "kadenmc")
  ) %>%
  filter(!run_ID %in% exclude)

cal01_08_sum <- cal01_08 %>%
  group_by(run_ID, Assay, Sample, run_date, plate_num, operator) %>%
  summarize(
    mean_signal = gmean(Signal)
  ) %>%
  full_join(
    high_blank_list %>% select(run_date, Assay, plate_num, operator, blank_low)
  )


#get ggplot default colors for 8 groups (8 CAL)
pal <- hue_pal()(8)

#tibble for indicating position for max CAL-08 signal
annotate_300 <- tribble(
  ~Assay, ~label, ~x, ~y,
  "MPXV A29L", 300, ymd(20240810), log2(500),
  "VACV A27L", 300, ymd(20240810), log2(500)
)

#tibble for marking validation plots
annotate_val <- tibble(
  Assay = c("MPXV A29L", "MPXV A35R", "MPXV B6R", "MPXV E8L", "MPXV M1R",
            "VACV A27L", "VACV A33R", "VACV B5R", "VACV D8L", "VACV L1R"),
  xmin = ymd(20240725), xmax = ymd(20240905),
  ymin = log2(16), ymax = log2(24)
)

Fig_S2 <- ggplot(cal01_08_sum, aes(x = run_date, y = log2(mean_signal), color = Sample)) +
  facet_wrap(
    vars(Assay),
    nrow = 2
  ) +
  geom_point(
    size = 1,
    shape = 21
  ) +
  geom_rect(
    data = annotate_val,
    inherit.aes = FALSE,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    fill = "cornflowerblue",
    color = NA
  ) +
  scale_y_continuous(
    labels = ~comma(2^.x),
    breaks = seq(4,log2(2000000),2),
    name = "Signal",
    limits = c(-10,25)
  ) +
  coord_cartesian(
    clip = "off",
    xlim = ymd(c("20240730", "20250228")),
    ylim = c(log2(16),log2(2000000))
  ) +
  geom_text(
    data = annotate_300,
    inherit.aes = FALSE,
    aes(x = x, y = y, label = label),
    color = "red",
    size = 6/.pt
  ) +
  scale_color_manual(
    values = c(pal[1], pal[8])
  ) +
  geom_point(
    data = cal01_08_sum %>% filter(blank_low == "no"),
    aes(x = run_date, y = log2(mean_signal), fill = blank_low),
    inherit.aes = FALSE,
    shape = 21
  ) +
  scale_fill_manual(
    values = "red",
    name = "CAL-01/CAL-08\nSignal ratio > 200"
  ) +
  scale_x_date(
    date_labels = "%b-%Y",
    date_breaks = "1 month",
    limits = ymd(c("20240620", "20250131")),
    name = "run date"
  ) +
  geom_hline(
    yintercept = log2(300),
    linetype = "dashed"
  ) +
  theme(
    panel.background = element_blank(),
    panel.grid = element_blank(),
    axis.line = element_line(color = "black"),
    axis.text.x = element_text(angle = 90, vjust = 0.5),
    axis.text = element_text(size = 6),
    plot.title = element_text(size = 8),
    axis.title = element_text(size = 8),
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 8, hjust = 0),
  )

ggsave("fig_output/Fig_S2.pdf", width = 6.5, height = 4, units = "in", device = cairo_pdf)



