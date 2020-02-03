library(dplyr)
library(geosphere)
library(stargazer)
library(gtools)

dim_names <- c("Food", "Water", "Shelter", "Health", "Education", "Energy", "Sanitation",
               "Relationships", "Clothing", "Violence", "Family Planning", "Environment",
               "Voice", "Time-use", "Work")

get_dim_estimates <- function(dimensions_list=dimensions_list, statistic='point', model, variable){
  
  # Set row index (which statistic we're getting)
  if (statistic == 'point'){
    selector <- 1
  } else if (statistic == 'lower') {
    selector <- 2
  } else if (statistic == 'upper') {
    selector <- 3
  } else if (is.numeric(statistic) && statistic >= 1 && statistic <= 3) {
    selector <- statistic
  } else {
    stop("Invalid statistic selector")
  }
  
  # Get estimates
  best_estimates <- c()
  for (i in 1:length(dimensions_list)){
    new_best <- try(dimensions_list[[i]][[model]][selector ,variable])
    if (length(new_best == 1) && is.numeric(new_best)) {
      best_estimates[i] <- new_best
    } else {
      best_estimates[i] <- NA
    }
  }
  # names(best_estimates) <- dim_names
  
  return(best_estimates)
  
}

# Paste problem
paste2 <- function(string) {
  if(is.na(string)) {
    return(NA)
  }
  else{
    return(paste0(string))
  }
}

anova_estimates <- function(anovaslist=anovas_list){
  df <- data.frame(matrix(nrow = length(anovaslist), ncol = 4))
  
  for (i in 1:length(anovaslist)){
    df[i, 1] <- as.numeric(try(anovaslist[[i]][1, "logLik"]))
    df[i, 2] <- as.numeric(try(anovaslist[[i]][2, "logLik"]))
    df[i, 3] <- as.numeric(try(anovaslist[[i]][1, "Pr(>Chisq)"]))
    df[i, 4] <- as.numeric(try(anovaslist[[i]][2, "Pr(>Chisq)"]))
  }
  
  colnames(df) <- c("sex_model logLik", "hhh_model logLik", "sex_model Pr(>Chisq)",
                    "hhh_model Pr(>Chisq)")
  return(df)
}

extract_coefficient <- function(modelno, variable, dimlist = dimensions_list){
  hhh_ests_point <- get_dim_estimates(dimlist, "point", modelno, variable)
  hhh_ests_lower <- get_dim_estimates(dimlist, "lower", modelno, variable)
  hhh_ests_upper <- get_dim_estimates(dimlist, "upper", modelno, variable)
  hhh_ests_df <- data.frame(point = hhh_ests_point, lower = hhh_ests_lower, upper = hhh_ests_upper)
  row.names(hhh_ests_df) <- names(dimlist)
  
  return(hhh_ests_df)
}

compile_all_coef <- function(scores, coefficients_list, dimlist, modelno){
  sex_all_coef_df <- data.frame(matrix(nrow = 15, ncol = length(coefficients_list)))
  row.names(sex_all_coef_df) <- all_themes
  colnames(sex_all_coef_df) <- colnames(dimensions_list[[1]][[modelno]])
  for (i in 1:length(coefficients_list)){
    coef_name <- colnames(dimensions_list[[1]][[modelno]])[i]
    sex_all_coef_df[coef_name] <- paste(round(coefficients_list[[coef_name]]$point, 2),
                                        ' (',
                                        round(coefficients_list[[coef_name]]$lower, 2),
                                        ' / ',
                                        round(coefficients_list[[coef_name]]$upper, 2), 
                                        ')')
  }
  
  return(sex_all_coef_df)
}

# coefficients_list_ind <- list()
# for(i in 1:length(colnames(dimensions_list[[1]][[1]]))){
#   coef_name <- colnames(dimensions_list[[1]][[1]])[i]
#   try(coefficients_list_ind[[coef_name]] <- extract_coefficient(dimensions_list, modelno = 1, variable = coef_name))
# }

get_combs <- function(vector, key_set = NA) {
  ls <- list()
  index <- 1
  
  n = length(vector)
  for (i in 1:n) {
    if (i == 1) {
      for (j in 1:n) {
        ls[[index]] <- (vector[j])
        index <- index + 1
      }
    }
    else {
      for (j in 1:(n - i + 1)) {
        for (k in (j + i - 1):n) {
          newval <- c(vector[j:(j + i - 2)], vector[k])
          ls[[index]] <- (newval)
          index <- index + 1
        }
      }
    }
  }
  
  new_ls <- list()
  new_indexer <- 1
  
  if (!is.na(key_set)){
    for (i in 1:length(ls)){
      if (key_set %in% ls[[i]]) {
        new_ls[[new_indexer]] <- ls[[i]]
        new_indexer <- new_indexer + 1
      }
    }
    return(new_ls)
  }
  return(ls)
}

# Nearest weather station
nearest_weather <- function(fiji_data){
  # Read in list of points
  weather_coords <- read.csv('weathercoords.csv')
  weather_stations <- filter(weather_coords, type == "WS")
  tikinas <- filter(weather_coords, type == "TK")
  
  # Find closest weather station to each tikina
  for (i in 1:nrow(tikinas)){
    
    # Check valid values
    if (is.na(tikinas$longitude[i]) || is.na(tikinas$longitude[i])){
      print(paste("Skipped tikina", i, "due to missing values"))
      next
    }
    
    # Set variable to trash values
    closest_ws <- Inf
    names(closest_ws) <- "None"
    
    #Iterate through weather stations for each Tikina
    for (j in 1:nrow(weather_stations)) {
      temp_ws <- distHaversine(c(tikinas$longitude[i], tikinas$latitude[i]), 
                    c(weather_stations$longitude[j], weather_stations$latitude[j]))
      if (temp_ws < closest_ws) {
        closest_ws <- temp_ws
        names(closest_ws) <- weather_stations$name[j]
      }
    }
    
    # Set closest weather station for Tikina
    tikinas[i, "closest_WS"] <- names(closest_ws)[1]
    print(paste0("Closest WS to tikina ", i, " is ", names(closest_ws)[1]))
  }
  
  # Join tikina coords to fiji_data df
  colnames(tikinas)[1] <- "tikina"
  fiji_data$tikina <- as.factor(fiji_data$tikina)
  tikinas <- select(tikinas, tikina, closest_WS)
  fiji_data <- left_join(fiji_data, tikinas, by = "tikina")
  
  # Join rainfall
  rainfall <- read_csv('rainfallsests.csv')
  fiji_data <- fiji_data %>%
    left_join(rainfall, by = "closest_WS")
  
  for (i in 1:nrow(fiji_data)){
    if (is.na(month(fiji_data$date.interview[i]))){
      next
    } else if (month(fiji_data$date.interview[i]) == 7){
      fiji_data[i, "rainfallperc"] <- fiji_data[i, "percjul"] * 100
      fiji_data[i, "rainfallobs"] <- fiji_data[i, "obsjul"]
    } else if (month(fiji_data$date.interview[i]) == 8){
      fiji_data[i, "rainfallperc"] <- fiji_data[i, "percaug"] * 100
      fiji_data[i, "rainfallobs"] <- fiji_data[i, "obsaug"]
    }
  }
  
  # Return dataframe
  return (fiji_data)
}

write_coefficients <- function(model, filename, anova_type = 1) {
  options(contrasts=c("contr.sum", "contr.poly"))
  model_coefs <- coef(summary(model))
  # data_df <- as.data.frame(model_coefs, row.names = rownames(model_coefs))
  rownames(model_coefs) <- sapply(rownames(model_coefs), str_replace, "TRUE", "")
  anova_results <- as.data.frame(anova(model, type = anova_type))
  test <- merge(model_coefs, anova_results, by="row.names", all.y = TRUE)
  write.csv(test, filename)
}

stargazer_clm <- function(model_list, output_file, ...) {
  options(contrasts=c("contr.sum", "contr.poly"))
  pvals <- list()
  
  # Get p vals
  for (i in 1:length(model_list)) {
    mo_name <- names(model_list)[i]
    pvals[[mo_name]] <- eval(parse(text = paste0("anova(model_list[[\"", mo_name ,"\"]], type = 3)$`Pr(>Chi)`")))
  }
  
  # Format stargazer
  collapsed_mo_names <- paste(names(model_list), collapse = " , ")
  star_text <- paste0("stargazer(", collapsed_mo_names, ", type = 'html', out = output_file, 
                      p = pvals, report=('vc*ps'), ci = T)")
  eval(parse(text = star_text))
}

vglm_gaze <- function (model, output_file) {
  sumtable <- summary(model)
  sumtable <- as.data.frame(sumtable@coef3)
  restable <- sumtable[c("Estimate","Pr(>|z|)")]
  restable["Sig"] <- restable["Pr(>|z|)"] <= 0.05
  
  for (i in 1:nrow(restable)) {
    if (restable[i, "Sig"]) restable[i, "Sig"] <- "*"
    else restable[i, "Sig"] <- " "
  }
  
  colnames(restable) <- c("Estimate", "Pr(>|z|)", "Sig")
  return(restable)
}

interaction.ggplot <- function (test_data, score, clm) {
  eval(parse(text=paste0("fiji_cleaned %>% select(sex, rainfallperc, age, disability3,", score , ")"))) %>%
    na.omit() -> test_data
  
  test_data$sex <- as.character(test_data$sex) %>%
    str_replace("1", "Male") %>%
    str_replace("2", "Female")
  
  preds <- as.data.frame(predictvglm(clm))
  test_data <- bind_cols(preds, as.data.frame(test_data)) %>%
    gather(key = "Threshold",
           value = "Logit",
           `logitlink(P[Y>=2])`,
           `logitlink(P[Y>=3])`,
           `logitlink(P[Y>=4])`)
  
  int_plot <- ggplot(data = test_data,
         aes(x = rainfallperc, y = Logit, linetype = sex, group=sex)) +
    stat_summary(fun.y=mean, geom="point") +
    stat_summary(fun.y=mean, geom="line") +
    facet_wrap(~Threshold,  ncol=1) +
    ggtitle(paste0(score, " interaction plots")) +
    theme_bw()
  
  return(int_plot)
}
