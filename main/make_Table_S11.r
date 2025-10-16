# Table_S11

# Greninger Lab - Jon Reed
# 2025-Oct-15

source("MPXV_paper_functions.r")

package_list <- c("tidyverse",
                  "openxlsx2")

invisible(lapply(package_list, check_and_load))


dat_neut_MSD <- read_csv("2025-Oct-15_compiled_data.csv")


dat_neut_MSD <- dat_neut_MSD %>%
  mutate(
    Sample_ID = if_else(group == "Vaccine B - 9 month", str_replace(Sample_ID, "M", "V"), Sample_ID)
  )

wb <- wb_workbook()
wb$add_worksheet("Table S11")
wb$add_data("Table S11", dat_neut_MSD, na.strings = NULL)


# Reduce the number of digits after decimal to 2, for certain columns
wb <- wb_add_numfmt(
  wb,
  sheet = "Table S11",
  dims = wb_dims(
    rows = 2:(nrow(dat_neut_MSD)+1),
    cols = c(10:14,18:20)
  ),
  numfmt = "0.00"
)

#Center columns
wb <- wb_add_cell_style(
  wb,
  sheet = "Table S11",
  dims = wb_dims(
    rows = 1:(nrow(dat_neut_MSD)+1),
    cols = c(4:20)
  ),
  vertical = "center",
  horizontal = "center"
)

#Auto adjust column widths
wb_set_col_widths(
  wb,
  sheet = "Table S11",
  cols = c(1:22),
  widths = "auto"
)

#Freeze the top row
wb_freeze_pane(wb, sheet = "Table S11", first_active_row = 2, first_row = TRUE)

wb$save(file.path("fig_output", "Table_S11.xlsx"))