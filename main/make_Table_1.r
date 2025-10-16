# Table 1 - Cohort characteristics and demographics

# Greninger Lab - Jon Reed
# 2025-Oct-15

# load functions
source("MPXV_paper_functions.r")

# install packaages
package_list <- c("tidyverse",
                  "openxlsx2",
                  "scales")
                
invisible(lapply(package_list, check_and_load))

# load data
dat <- read_csv("2025-Oct-15_compiled_data.csv")

# formatting helper
format_meta <- function(dat) {
  
  col_name <- dat$group[[1]]
  
  dat %>%
    summarize(
      count = length(unique(PTID)),
      male = paste0(sum(Sex == "M"), " (", percent_format(accuracy = 1)(sum(Sex == "M")/count),")"),
      female = paste0(sum(Sex == "F"), " (", percent_format(accuracy = 1)(sum(Sex == "F")/count),")"),
      mean = paste0(round(mean(Age),0), " (", min(Age), " - ", max(Age), ")"),
      HIV_pos = if_else(is.na(sum(HIV == "yes")),
                        NA_character_,
                        paste0(sum(HIV == "yes"), " (", percent_format(accuracy = 1)(sum(HIV == "yes")/count), ")")),
      born_lt_1972 = paste0(sum(YOB_est < 1972), " (", percent_format(accuracy = 1)(sum(YOB_est < 1972)/count), ")"),
      history_vacc = paste0(sum(vaccine_history == "yes"), " (", percent_format(accuracy = 1)(sum(vaccine_history == "yes")/count), ")"),
      months_post_MPXV = if_else(is.na(mean(months_post_MPXV)),
                                 NA_character_,
                                 paste0(round(mean(months_post_MPXV), 1), " (", round(min(months_post_MPXV), 1), " - ", round(max(months_post_MPXV), 1), ")")), 
      months_post_vaccine = if_else(is.na(mean(months_post_vaccine)),
                                    NA_character_,
                                    paste0(round(mean(months_post_vaccine),1), " (", round(min(months_post_vaccine),1), " - ", round(max(months_post_vaccine),1), ")"))
    ) %>%
    mutate(
      across(everything(), as.character)
    ) %>%
    pivot_longer(
      cols = everything(), 
      names_to = "Metric", 
      values_to = col_name
    )
}

# columsn of metadata to select
cols_to_select <- c("group", "PTID", "Sex", "Age", "HIV", "YOB_est", "vaccine_history", "months_post_MPXV", "months_post_vaccine")

# read in metadata - pediatric
temp1 <- dat %>% filter(group == "pediatric") %>% select(all_of(cols_to_select)) %>% distinct()

# read in metadata - rubella IgG+
temp2 <- dat %>% filter(group == "rubella IgG+") %>% select(all_of(cols_to_select)) %>% distinct()

# read in metadata - VRC, vaccine set
temp3 <- dat %>% 
  # filter all the Vaccine A sampling, don't calculate the months_post_vaccine since it
  # varied, then get distinct rows (should boil down to individuals)
  filter(group %in% c("Vaccine A - 0 month", "Vaccine A - 1 month", "Vaccine A - 2 month")) %>% 
  select(all_of(cols_to_select)) %>% 
  mutate(
    group = "Vaccine A",
    months_post_vaccine = NA_integer_
  ) %>%
  distinct()

#read in metadata - golden, vaccine set
temp4 <- dat %>% filter(group == "Vaccine B - 9 month") %>% select(all_of(cols_to_select)) %>% distinct()

#read in metadata - VRC, MPXV set
temp5 <- dat %>% filter(group == "MPXV - 2-6 month") %>% select(all_of(cols_to_select)) %>% distinct()         

#read in metadata - golden, MPXV set
temp6 <- dat %>% filter(group == "MPXV - 10 month") %>% select(all_of(cols_to_select)) %>% distinct()   

#format table
all_meta_formated <- bind_cols(
  format_meta(temp1),
  format_meta(temp2) %>% select(!Metric),
  # Manually modified months_post_vaccine to "multiple collections" since collections were at 0, 1, 2 months.
  # other cohorts had single collection - except for on patid in golden set
  format_meta(temp3) %>% mutate(`Vaccine A` = if_else(Metric == "months_post_vaccine", "multiple collections", `Vaccine A`)) %>% select(!Metric),
  format_meta(temp4) %>% select(!Metric),
  format_meta(temp5) %>% select(!Metric),
  # Manually modified male row since total samples is 18, but total individuals is 17
  format_meta(temp6) %>% mutate(`MPXV - 10 month` = if_else(Metric == "male", "17 (100%)", `MPXV - 10 month`)) %>% select(!Metric)
) 

#update row names 
x <- all_meta_formated[all_meta_formated$Metric == "count",][-1]
y <- c("Pediatric", "Rubella IgG+", "Vaccine cohort A", "Vaccine cohort B", "MPXV 2-6 mo", "MPXV 10 mo")
new_col_name <- c("Metric", paste0(y, " (n = ", x, ")"))

all_meta_formated <- all_meta_formated %>%
  rename_with(~ new_col_name) %>%
  filter(Metric != "count")

#write to .xlsx file - finish formatting in Word
wb <- wb_workbook()
wb$add_worksheet("Table 1")
wb$add_data("Table 1", all_meta_formated, na.strings = "NA")

#Auto adjust column widths
wb_set_col_widths(
  wb,
  sheet = "Table 1",
  cols = c(1:7),
  widths = "auto"
)

wb$save(file.path("fig_output", "Table_1.xlsx"))


