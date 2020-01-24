library(tidyverse)
library(ordinal)
library(UpSetR)
library(gridExtra)
library(lubridate)
library(stargazer)
library(VGAM)
source('importclean.R')
source('generalhelpers.R')
source('regressions.R')

### Note due to a bug, this script sometimes throws a "subscript out of bounds" error when
### calculating Type III ANOVA. Restarting the R session should resolve this issue.

# Import data
rdsname <- "fiji_newscored_cleaned_Oct19.rds"
fiji_cleaned <- import.clean(rdsname)

## Descriptive statistics
freq_table <- apply(fiji_cleaned[c("score2", "score4", "score7", "score14")], 2, table)
# stargazer(freq_table, type = 'html', title = 'Frequency table for dimension scores', out = 'freqs.html',
#           summary = FALSE)
# base::summary(fiji_cleaned)

graph_dimensions <- function (data, title, ncol = 2) {
  require("gridExtra")
  # Make graphs
  graph_list <- list()
  for (i in 1:ncol(data)) {
    graph_list[[i]] <- ggplot(fiji_cleaned, aes_string(colnames(data)[i])) + geom_bar() + 
      theme_bw()
  }
  # Lay out graphs
  references <- paste0("graph_list[[", as.character(1:length(graph_list)), "]]")
  collapse_references <- paste(references, collapse = ", ")
  grid_expression <- paste0("grid.arrange(", collapse_references, ", ncol = ", ncol, ", top = title)")
  graph_grid <- eval(parse(text = grid_expression))
  
  return(graph_grid)
}

desc_graphs <- graph_dimensions(fiji_cleaned[c("score2", "score4", "score7", "score14")], 
                                'Distribution of responses in the dimensions studied', 2)

## Spearman correlation
corr_table <-fiji_cleaned[c(paste0('score', 1:9), paste0('score', 12:15), 
                                   'rainfallperc', 'rainfallobs', 'urban', 'rural', 'informal', 'sex')] %>%
  sapply( as.numeric ) %>%
  cor(use = 'pairwise.complete.obs', method = 'spearman') %>%
  round(2)

write.csv(corr_table, 'spearmans.csv')
  
# # Get deprivations df
# deprivations_corr <- fiji_cleaned %>%
#   ungroup() %>%
#   select(score1, score2, score3, score4, score5, score6, score7, score8, score9,
#          score12, score13, score14, score15, sex, urban, rural, age) %>%
#   sapply(as.numeric) %>%
#   cor(use = 'pairwise.complete.obs', method = 'spearman')
# 
# # Write output
# write.csv(deprivations_corr, 'spearman_output.csv')

deprivations <- fiji_cleaned %>%
  ungroup() %>%
  select(score1, score2, score3, score4, score5, score6, score7, score8, score9,
         score12, score13, score14, score15) %>%
  sapply(as.numeric)

deprivations <- deprivations < 2
for (i in length(colnames(deprivations)))
  deprivations[i] <- as.numeric(deprivations[i])

deprivations <- as.data.frame(deprivations)
deprivations["sector"] <- fiji_cleaned$sector

# Make UpSet

sets1 <- c("score2", "score4", "score7", "score14")
rural <- filter(deprivations, sector == 0)
urban <- filter(deprivations, sector == 1)
informal <- filter(deprivations, sector == 2)

upplot_rural <- upset(urban, nintersects = NA, intersections = get_combs(sets1, key_set = "score2"),
                order.by = "freq")

upplot_urban <- upset(rural, nintersects = NA, intersections = get_combs(sets1, key_set = "score2"),
                      order.by = "freq")

upplot_informal <- upset(informal, nintersects = NA, intersections = get_combs(sets1, key_set = "score2"),
                         order.by = "freq")

# upplot_rural
# upplot_urban
# upplot_informal

## Regression

# Get rural data
rural <- filter(fiji_cleaned, sector == 1)

# Water
clm_poa <- clm(score2 ~ (sex * rainfallperc) + rural + urban + age, data = fiji_cleaned)
waterclm <- regress_po(fiji_cleaned, "score2", interaction = TRUE)[["ppo_model"]]
clm_poa <- clm(score2 ~ sex + rainfallperc + rural + urban + age, data = fiji_cleaned)
waterclm_noint <- regress_po(fiji_cleaned, "score2", interaction = FALSE)[["ppo_model"]]

# # Health
clm_poa <- clm(score4 ~ (sex * rainfallperc) + rural + urban + age, data = fiji_cleaned)
healthclm <- regress_po(fiji_cleaned, "score4", interaction = TRUE)[["ppo_model"]]
clm_poa <- clm(score4 ~ sex + rainfallperc + rural + urban + age, data = fiji_cleaned)
healthclm_noint <- regress_po(fiji_cleaned, "score4", interaction = FALSE)[["ppo_model"]]

# # Time use
clm_poa <- clm(score14 ~ (sex * rainfallperc) + rural + urban + age, data = fiji_cleaned)
timeclm <- regress_po(fiji_cleaned, "score14", interaction = TRUE)[["ppo_model"]]
clm_poa <- clm(score14 ~ sex + rainfallperc + rural + urban + age, data = fiji_cleaned)
timeclm_noint <- regress_po(fiji_cleaned, "score14", interaction = FALSE)[["ppo_model"]]

# # Sanitation
clm_poa <- clm(score7 ~ (sex * rainfallperc) + rural + urban + age, data = fiji_cleaned)
sanitationclm <- regress_po(fiji_cleaned, "score7", interaction = TRUE)[["ppo_model"]]
clm_poa <- clm(score7 ~ sex + rainfallperc + rural + urban + age, data = fiji_cleaned)
sanitationclm_noint <- regress_po(fiji_cleaned, "score7", interaction = FALSE)[["ppo_model"]]

# 
# # Write coefficients to csvs
# write_coefficients(waterclm, 'water.csv')
# write_coefficients(waterclm_noint, 'water_noint.csv')
# 
# write_coefficients(healthclm, 'health.csv')
# write_coefficients(healthclm_noint, 'health_noint.csv')
# 
# write_coefficients(timeclm, 'time.csv')
# write_coefficients(timeclm_noint, 'time_noint.csv')
# 
# write_coefficients(sanitationclm, 'sanitation.csv')
# write_coefficients(sanitationclm_noint, 'sanitation_noint.csv')
# 
# # HTML table
# 
## Define function
# tables_from_keyword <- function(dimension_keyword) {
#   print("Reminder: clm models named {{keyword}}clm and {{keyword}}clm_noint must be available.")
#   eval(parse(text = paste0(dimension_keyword, "_models <- list()")))
#   # eval(parse(text = paste0(dimension_keyword, "_models[[\'", dimension_keyword, "clm_nocont\']] <- ", dimension_keyword, "clm_nocont")))
#   eval(parse(text = paste0(dimension_keyword, "_models[[\'", dimension_keyword, "clm_noint\']] <-", dimension_keyword, "clm_noint")))
#   eval(parse(text = paste0(dimension_keyword, "_models[[\'", dimension_keyword, "clm\']] <- ", dimension_keyword, "clm")))
#   eval(parse(text = paste0("stargazer_clm(", dimension_keyword, "_models, '", dimension_keyword, "_out.html')")))
# }
# # 
# tables_from_keyword('health')
# tables_from_keyword('water')
# tables_from_keyword('time')
# tables_from_keyword('sanitation')

combn_fewer <- function(x) {
  len <- length(x)
  vec <- vector('character', 0)
  for (i in 1:len) {
    vec <- c(vec, combn(x, i))
  }
  return(vec)
}

# Export tables
stargazer(vglm_gaze(waterclm), vglm_gaze(waterclm_noint), type = 'html', out ='water.html', summary = FALSE)
stargazer(vglm_gaze(healthclm), vglm_gaze(healthclm_noint), type = 'html', out ='health.html', summary = FALSE)
stargazer(vglm_gaze(sanitationclm), vglm_gaze(sanitationclm_noint), type = 'html', out ='sanitation.html', summary = FALSE)
stargazer(vglm_gaze(timeclm), vglm_gaze(timeclm_noint), type = 'html', out ='time.html', summary = FALSE)


# print("Score 4")
# clm_poa <- clm(score4 ~ (sex * rainfallperc) + rural + urban + age, data = fiji_cleaned)
# result4 <- regress_po(fiji_cleaned, "score4")
# 
# print("Score 7")
# clm_poa <- clm(score7 ~ (sex * rainfallperc) + rural + urban + age, data = fiji_cleaned)
# result7 <- regress_po(fiji_cleaned, "score7")
# 
# print("Score 14")
# clm_poa <- clm(score14 ~ (sex * rainfallperc) + rural + urban + age, data = fiji_cleaned)
# result14 <- regress_po(fiji_cleaned, "score14")