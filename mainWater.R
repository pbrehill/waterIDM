library(tidyverse)
library(ordinal)
library(UpSetR)
library(gmodels)
library(gridExtra)
library(lubridate)
library(stargazer)
library(VGAM)
source('importclean.R')
source('generalhelpers.R')
source('regressions.R')

# Import data
rdsname <- "fiji_newscored_cleaned_Oct19.rds"
fiji_cleaned <- import.clean(rdsname)
fiji_cleaned <- filter(fiji_cleaned, sector == 0)

# Descriptive statistics
freq_table <- apply(fiji_cleaned[c("score2", "score4", "score7", "score14")], 2, table)

# Age brackets
fiji_cleaned[["brackets"]] <- cut(fiji_cleaned$age, breaks=c(17, 36, 51, 66, 200), right = FALSE)
CrossTable(fiji_cleaned$brackets, fiji_cleaned$sex, prop.r = FALSE, prop.t = FALSE, prop.chisq = FALSE)

# Ethnicities
CrossTable(fiji_cleaned$ethnicity, fiji_cleaned$sex, prop.r = FALSE, prop.t = FALSE, prop.chisq = FALSE)

# Relationship to respondent
CrossTable(fiji_cleaned$relationship, fiji_cleaned$sex, prop.r = FALSE, prop.t = FALSE, prop.chisq = FALSE)

# Disability
CrossTable(fiji_cleaned$disability3, fiji_cleaned$sex, prop.r = FALSE, prop.t = FALSE, prop.chisq = FALSE)


# Graph dimensions
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

desc_graphs <- graph_dimensions(fiji_cleaned[c("score1", "score2", "score4", "score7", "score14")], 
                                'Distribution of responses in the dimensions studied', 2)

## Spearman correlation
corr_table <-fiji_cleaned[c(paste0('score', 1:9), paste0('score', 12:15), 
                            'rainfallperc', 'rainfallobs', 'urban', 'rural', 'informal', 'sex', 'score15.1.1')] %>%
  sapply( as.numeric ) %>%
  cor(use = 'pairwise.complete.obs', method = 'spearman') %>%
  round(2)


write.csv(corr_table, 'spearmans.csv')

# Water
clm_poa <- clm(score2 ~ (sex * rainfallperc) + age + disability3, data = fiji_cleaned)
waterclm <- regress_po(fiji_cleaned, "score2", interaction = TRUE)[["ppo_model"]]
clm_poa <- clm(score2 ~ sex + rainfallperc + age+ disability3, data = fiji_cleaned)
waterclm_noint <- regress_po(fiji_cleaned, "score2", interaction = FALSE)[["ppo_model"]]

# # Health
clm_poa <- clm(score4 ~ (sex * rainfallperc) + age+ disability3, data = fiji_cleaned)
healthclm <- regress_po(fiji_cleaned, "score4", interaction = TRUE)[["ppo_model"]]
clm_poa <- clm(score4 ~ sex + rainfallperc + age+ disability3, data = fiji_cleaned)
healthclm_noint <- regress_po(fiji_cleaned, "score4", interaction = FALSE)[["ppo_model"]]

# # Time use
clm_poa <- clm(score14 ~ (sex * rainfallperc) + age + disability3, data = fiji_cleaned)
timeclm <- regress_po(fiji_cleaned, "score14", interaction = TRUE)[["ppo_model"]]
clm_poa <- clm(score14 ~ sex + rainfallperc + age + disability3, data = fiji_cleaned)
timeclm_noint <- regress_po(fiji_cleaned, "score14", interaction = FALSE)[["ppo_model"]]

# # Sanitation
clm_poa <- clm(score7 ~ (sex * rainfallperc) + age + disability3, data = fiji_cleaned)
sanitationclm <- regress_po(fiji_cleaned, "score7", interaction = TRUE)[["ppo_model"]]
clm_poa <- clm(score7 ~ sex + rainfallperc + age + disability3, data = fiji_cleaned)
sanitationclm_noint <- regress_po(fiji_cleaned, "score7", interaction = FALSE)[["ppo_model"]]

# # Food
clm_poa <- clm(score1 ~ (sex * rainfallperc) + age + disability3, data = fiji_cleaned)
foodclm <- regress_po(fiji_cleaned, "score1", interaction = TRUE)[["ppo_model"]]
clm_poa <- clm(score1 ~ sex + rainfallperc + age + disability3, data = fiji_cleaned)
foodclm_noint <- regress_po(fiji_cleaned, "score1", interaction = FALSE)[["ppo_model"]]

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
stargazer(vglm_gaze(foodclm), vglm_gaze(foodclm_noint), type = 'html', out ='food.html', summary = FALSE)

# Export interaction plots
int_plot_list <- list()
int_plot_list[[1]] <- interaction.ggplot(fiji_cleaned, "score1", foodclm)
int_plot_list[[2]] <- interaction.ggplot(fiji_cleaned, "score2", waterclm)
int_plot_list[[3]] <- interaction.ggplot(fiji_cleaned, "score4", healthclm)
int_plot_list[[4]] <- interaction.ggplot(fiji_cleaned, "score7", sanitationclm)
int_plot_list[[5]] <- interaction.ggplot(fiji_cleaned, "score14", timeclm)

for (i in 1:length(int_plot_list)) {
  ggsave(paste0("plot", i, ".png") ,int_plot_list[[i]])
}
