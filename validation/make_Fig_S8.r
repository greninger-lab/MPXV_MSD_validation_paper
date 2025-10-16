# Freeze-Thaw
# Fig S8

# Greninger Lab - Jon Reed
# 2025-Oct-15

source("../main/MPXV_paper_functions.r")

package_list <- c("tidyverse",
                  "ggplot2")

invisible(lapply(package_list, check_and_load))

dat <- read_csv("2025-Oct-15_MPXV_paper_validation_data.csv")
dat_assay_para <- read_csv("assay_limits/assay_para_validated.csv")

ft_Sample <- tribble(
  ~Sample, ~status,
  "M14", "Vaccine",
  "M42", "MPXV",
  "M48", "MPXV",
  "M57", "Vaccine"
)

dat_ft <- dat %>%
  filter(str_detect(Sample,"_FT"))

dat_ft_sum <- dat_ft %>%
  group_by(Sample, Assay) %>%
  summarise(
    geo_mean = gmean(conc_calc)
  ) %>%
  mutate(
    freeze_thaw = str_extract(Sample,"(?<=FT_)."),
    Sample = str_extract(Sample,"M.."),
    .after = Sample
  ) %>%
  full_join(ft_Sample, by = "Sample")


dat_ft_sum <- dat_ft_sum %>%
  mutate(
    Sample_alt = case_when(
      status == "Vaccine" ~ str_replace(Sample,"M","V"),
      status == "MPXV" ~ Sample,
      TRUE ~ NA_character_
    ),
    .after = Sample
  )

#Not included in supplement
a <- ggplot(dat_ft_sum, aes(x = freeze_thaw, y = log2(geo_mean), group = Sample_alt, fill = Sample_alt)) +
  facet_wrap(
    vars(Assay),
    nrow = 2,
    axes = "all_x"
  ) +
  geom_point(
    shape = 21
  ) +
  xlab("freeze-thaw cycles") +
  ylab("titer geometric mean (AU/mL)") +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5),
    panel.grid = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(color = "black")
  )


dat_ft_sum_ratio <- dat_ft_sum %>%
  pivot_wider(
    names_from = freeze_thaw,
    names_prefix = "FT_",
    values_from = geo_mean
  ) %>%
  mutate(
    `ratio FT 0/1` = FT_0/FT_1,
    `ratio FT 0/2` = FT_0/FT_2,
    `ratio FT 0/3` = FT_0/FT_3
  ) %>%
  select(Assay,Sample_alt,`ratio FT 0/1`, `ratio FT 0/2`, `ratio FT 0/3`) %>%
  pivot_longer(
    cols = c("ratio FT 0/1" , "ratio FT 0/2", "ratio FT 0/3"),
    names_to = "freeze_thaw",
    values_to = "ratio"
  )


Fig_S8 <- ggplot(dat_ft_sum_ratio, aes(x = freeze_thaw, y = ratio, fill = Sample_alt)) +
  facet_wrap(
    vars(Assay),
    axes = "all_x",
    nrow = 2
  ) +
  geom_point(
    shape = 21
  ) +
  labs(
    fill = "Sample"
  ) +
  scale_y_continuous(
    limits = c(0,2)
  ) +
  ylab("ratio of freeze-thaw 0 AU/mL over subsequent freeze-thaws 1,2, and 3 AU/mL") +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, size = 8),
    panel.grid = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(color = "black"),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 8)
  )

ggsave("fig_output/Fig_S8.pdf", width = 7.5, height = 4.5, units = "in", device = cairo_pdf)


