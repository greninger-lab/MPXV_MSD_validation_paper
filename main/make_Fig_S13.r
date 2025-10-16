# Fig S13

# Greninger Lab - Jon Reed
# 2025-Oct-15

source("MPXV_paper_functions.r")

package_list <- c("tidyverse",
                  "ggplot2",
                  "scales",
                  "ggrepel",
                  "ggtext")

invisible(lapply(package_list, check_and_load))


dat_neut_MSD <- read_csv("2025-Oct-15_compiled_data.csv")

dat_neut_MSD$Assay <- factor(dat_neut_MSD$Assay,
                             levels = c("MPXV A29L", "VACV A27L",
                                        "MPXV A35R", "VACV A33R",
                                        "MPXV B6R", "VACV B5R",
                                        "MPXV E8L", "VACV D8L",
                                        "MPXV M1R", "VACV L1R"))


x_labels <- tibble(
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



yoden_threshold <- c(
  "MPXV A35R" = 84.3,
  "VACV A33R" = 90.7,
  "MPXV B6R" = 79.9,
  "VACV B5R" = 108.8,
  "MPXV E8L" = 216.9,
  "VACV D8L" =	233.1,
  "MPXV M1R" = 45.1,
  "VACV L1R" = 40.5,
  "MPXV A29L" = 67.2,
  "VACV A27L" = 53.9)


dat_yoden_threshold <- enframe(yoden_threshold, name = "Assay", value = "threshold")

temp <- hcl.colors(10, "Greens")
green_pal <- c(temp[1], temp[2], temp[3], temp[4], temp[5])

temp <- hcl.colors(10, "Oranges")
orange_pal <- c(temp[1], temp[2], temp[3], temp[4], temp[5])

test <- c(green_pal, rev(orange_pal))


dat_neut_MSD2 <- dat_neut_MSD %>%
  filter(group %in% c("rubella IgG+", "pediatric")) %>%
  mutate(
    exceed_threshold = if_else(MSD_mean > yoden_threshold[Assay], TRUE, FALSE),
    .after = MSD_mean
  ) %>%
  nest(data = !Sample_ID) %>%
  mutate(
    exceed_threshold_any = map_lgl(data, ~any(.x$exceed_threshold)),
    exceed_threshold_count = map_int(data, ~sum(.x$exceed_threshold))
  ) %>%
  unnest(data)

list_ID_gt_5 <- dat_neut_MSD2 %>%
  nest(data = !c(Sample_ID, exceed_threshold_count)) %>%
  filter(exceed_threshold_count >= 5) %>%
  unnest(data) %>%
  filter(Assay == "VACV L1R") %>%
  select(Sample_ID, group, exceed_threshold_count, Assay,MSD_mean, ND50_mean) %>%
  mutate(
    label = paste0(Sample_ID, ", ND50 = ", sprintf("%0.1f", ND50_mean))
  )

dat_neut_MSD2_gt_1 <- dat_neut_MSD2 %>%  
  filter(exceed_threshold_any)

number_ped <- dat_neut_MSD2_gt_1 %>% filter(group == "pediatric") %>% pull(Sample_ID) %>% unique() %>% length()
number_rub <- dat_neut_MSD2_gt_1 %>% filter(group == "rubella IgG+") %>% pull(Sample_ID) %>% unique() %>% length()
total_ped <- dat_neut_MSD %>% filter(group == "pediatric") %>% pull(Sample_ID) %>% unique() %>% length()
total_rub <- dat_neut_MSD %>% filter(group == "rubella IgG+") %>% pull(Sample_ID) %>% unique() %>% length()

facet_label_parsed <- tribble(
  ~group, ~label,
  "pediatric", paste0(number_ped, "/", total_ped, '~"specimens in pediatric group exceed the Youden\'s threshold for at least one antigen"'),
  "rubella IgG+", paste0(number_rub, "/", total_rub, '~"specimens in " * rubella~IgG^"+" * " group exceed the Youden\'s threshold for at least one antigen"')
)



facet_map_parsed <- set_names(facet_label_parsed$label, facet_label_parsed$group)



Fig_S13 <- ggplot(dat_neut_MSD2_gt_1, aes(x = Assay, y = log2(MSD_mean), group = Sample_ID, color = factor(exceed_threshold_count))) +
  facet_wrap(
    vars(group),
    nrow = 2,
    labeller = labeller(group = as_labeller(facet_map_parsed, default = label_parsed))
  ) +
  geom_point(
    size = 0.5
  ) +
  geom_line(
    linewidth = 0.5,
    alpha = 0.5
  ) +
  scale_y_continuous(
    label = ~comma(2^.x),
    breaks = seq(1,13,1),
    limits = c(0,13),
    name = "Antibody Level (AU/mL)",
    expand = c(0,0)
  ) +
  scale_color_manual(
    values = test,
    name = str_wrap("Number of antigens exceeding threshold",20)
  ) +
  scale_x_discrete(
    labels = x_labels$label
  ) +
  geom_col(
    data = dat_yoden_threshold,
    aes(x = Assay, y = log2(threshold), fill = "Youden's threshold"),
    inherit.aes = FALSE,
    alpha = 0.25
  ) +
  scale_fill_manual(
    values = "gray42",
    name = NULL
  ) +
  geom_label_repel(
    data = list_ID_gt_5,
    aes(x = Assay, y = log2(MSD_mean), label = label),
    direction = "y",
    xlim = c(11,11.8),
    show.legend = FALSE,
    size = 5/.pt,
    max.overlaps = 20
  ) +
  coord_cartesian(
    clip = "off"
  ) +
  theme(
    axis.text.x = element_markdown(angle = 90, vjust = 0.5),
    panel.background = element_blank(),
    panel.grid = element_blank(),
    axis.line = element_line(color = "black"),
    axis.text = element_text(size = 7),
    axis.title = element_text(size = 8),
    legend.key.height = unit(10, "pt"),
    legend.title = element_text(size = 7),
    legend.text = element_text(size = 7),
    strip.text = element_text(hjust = 0, size = 8),
    strip.background = element_blank(),
    legend.position = "bottom",
    plot.margin = margin(t = 0, b = 0, l = 0, r = 50, unit = "pt")
  )

ggsave("fig_output/Fig_S13.pdf", width = 6, height = 6, units = "in", device = cairo_pdf)


tally <- dat_neut_MSD2 %>%
  nest(data = !c(Sample_ID, exceed_threshold_count, group)) %>%
  select(Sample_ID, exceed_threshold_count, group)


table <- tally %>%
  group_by(group) %>%
  count(exceed_threshold_count) %>%
  mutate(percent = n / sum(n) * 100)

