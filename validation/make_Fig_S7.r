# Accuracy control 1,2,3
# Figure S7A - RE of control 1,2
# Figure S7B - mean results of control 3, set control 3 limit

# Greninger Lab - Jon Reed
# 2025-Oct-15

source("../main/MPXV_paper_functions.r")

package_list <- c("tidyverse",
                  "ggplot2",
                  "ggbeeswarm",
                  "scales",
                  "patchwork")

invisible(lapply(package_list, check_and_load))

dat <- read_csv("2025-Oct-15_MPXV_paper_validation_data.csv")
dat_assay_para <- read_csv("assay_limits/assay_para_validated.csv")

controls <- c("Control 1", "Control 2", "Control 3")

dat_acc <- dat %>%
  filter(Sample %in% controls) %>%
  #Join with dat_assay_para to get LoD, LLoQ, and ULoQ
  full_join(dat_assay_para, by = "Assay")

dat_acc_sum <- dat_acc %>%
  group_by(Sample, Assay, day, LoD, LLoQ_final, manf_ULoQ, conc_exp) %>%
  summarize(
    gmean_conc = gmean(conc_calc),
    gmean_signal = gmean(Signal),
    count = n()
  ) %>%
  mutate(
    RE = ((gmean_conc - conc_exp)/conc_exp)
  )


dat_acc_sum_cont12 <- dat_acc_sum %>%
  filter(Sample %in% c("Control 1", "Control 2"))


Fig_S7A <- ggplot(dat_acc_sum_cont12, aes(x = Sample, y = RE, fill = Assay)) +
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
  ylab("relative error (%)") +
  theme(
    panel.grid = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(color = "black"),
    axis.text.x = element_text(angle = 90, size = 6, vjust = 0.5),
    axis.text.y = element_text(size = 8),
    legend.text = element_text(size = 6),
    legend.title = element_text(size = 8),
    axis.title.y = element_text(size = 8),
    axis.title.x = element_blank()
  ) 

dat_acc_control3 <- dat_acc %>%
  filter(Sample == "Control 3") %>%
  group_by(Assay) %>%
  summarise(
    control3_geo_mean = gmean(conc_calc),
    control3_geo_sd = gsd(conc_calc),
    control3_limit = exp(log(control3_geo_mean)+(log(control3_geo_sd)*3))
  ) %>%
  select(Assay, control3_geo_mean, control3_geo_sd, control3_limit)

number_of_assays <- dat_acc %>% filter(Sample == "Control 3") %>% select(day) %>% unique() %>% nrow()

Fig_S7B <- ggplot(dat_acc %>% filter(Sample == "Control 3"), aes(y = log2(conc_calc), x = "Control 3")) +
  geom_beeswarm(
    shape = 21,
    cex = 3
  ) +
  facet_wrap(
    vars(Assay),
    nrow = 2,
    axes = "all_x"
  ) +
  geom_hline(
    data = dat_acc_control3,
    aes(yintercept = log2(control3_limit)),
    linetype = "dashed"
  ) +
  scale_y_continuous(
    #labels = ~ formatC(round(2^.x,2), format ="f", digits = 2),
    limits = c(-6, -3),
    labels = label_math(expr = 2^.x)
  ) +
  geom_text(
    data = dat_acc_control3,
    aes(x = 1, y = log2(control3_limit) + 0.25, label = round(control3_limit,3)),
    color = "red"
  ) +
  ggtitle(paste0("Number of assays: ", number_of_assays)) +
  ylab("AU/mL") +
  xlab("Control specimen") +
  theme(
    panel.grid = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(color = "black"),
    plot.title = element_text(size = 8),
    axis.title = element_text(size = 8),
    axis.text = element_text(size = 7)
  )

Fig_S7 <- Fig_S7A / Fig_S7B +
  plot_annotation(tag_levels = "A")

ggsave("fig_output/Fig_S7.pdf", width = 7.5, height = 7, units = "in", device = cairo_pdf)

