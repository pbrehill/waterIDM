source("generalhelpers.R")
library(dplyr)


clean_date <- function(date) {
  str_date <- as.character(date)
  if (nchar(str_date) < 6) {
    str_date <- paste0("0", str_date)
  }
  proper_date <- as.Date.character(str_date, format="%d%m%y")
  return(proper_date)
}


import.clean <- function(rds)
{
  # Import data
  fiji_scored <- readRDS(file=rds)
  
   fiji_jst_scores <- fiji_scored # %>%
  #    dplyr::select(sex, age.categories, area, EA, sector, ethnicity, disability3, relationship,
  #                  score1, score2, score3,
  #                  score4, score5, score6, score7, score8, score9, score12, score13, score14, score15)
  
  for (i in 1:ncol(fiji_jst_scores)){
    if (grepl("score", colnames(fiji_jst_scores)[i]))
    {
      fiji_jst_scores[i] <- fiji_jst_scores[i] %>%
      round() %>%
      sapply(paste0) %>%
      factor(ordered = TRUE, exclude = c(NA, "NA"))
    }
  }
  
  fiji_jst_scores <- fiji_scored %>%
    select(HHID_unique, age) %>%
    bind_cols(fiji_jst_scores)
  
  # Set household head
  fiji_jst_scores <- fiji_scored %>%
    filter(relationship == 1) %>%
    select(sex, HHID_unique) %>%
    dplyr::rename(hhh_sex = sex) %>%
    right_join(fiji_jst_scores, by = "HHID_unique")
  
  # Set ethnicity dummy variables
  fiji_jst_scores <- fiji_scored %>%
    mutate(itaukei = ethnicity == 1, indian = ethnicity == 2, other.ethnicity = ethnicity > 2) %>%
    select(itaukei, indian, other.ethnicity) %>%
    bind_cols(fiji_jst_scores)
  
  # Set sector dummy variables
  fiji_jst_scores <- fiji_scored %>%
    mutate(rural = sector == 0, urban = sector == 1, informal = sector == 2) %>%
    select(rural, urban, informal) %>%
    bind_cols(fiji_jst_scores)
  
  # Clean dates
  newdatevec <- c()
  for(i in 1:length(fiji_jst_scores$date.interview)){
    newdate <- clean_date(fiji_jst_scores$date.interview[i])
    newdatevec[i] <- newdate
  }
  fiji_jst_scores$date.interview <- as.Date(fiji_jst_scores$date.interview, origin = "1970-01-01")
  fiji_jst_scores$date.interview <- newdatevec
  fiji_jst_scores$date.interview <- as.Date(fiji_jst_scores$date.interview, origin="1970-01-01")
  
  # Convert factors
  fiji_jst_scores$hhh_sex <- as.factor(fiji_jst_scores$hhh_sex)
  fiji_jst_scores$HHID_unique <- as.factor(fiji_jst_scores$HHID_unique)
  
  # Add weather stations
  fiji_jst_scores <- nearest_weather(fiji_jst_scores)
  
  return(fiji_jst_scores)
}
