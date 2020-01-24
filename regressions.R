library(tidyverse)

regression <- function(score, model, model_name, ivs, data=fiji_jst_scores){
  reg_ind_vars <- ivs
  if (model == "clmm") {
    reg_ind_vars <- append(reg_ind_vars, "(1|HHID_unique)")
  }
  fiji_reg_data <- data
  # # Remove incomplete cases to allow comparison of models
  # complete_cases_text <- paste0(deparse(substitute(data)), "[complete.cases(fiji_jst_scores[c('sex', 'hhh_sex')]), ]")
  # fiji_reg_data <- eval(parse(text = complete_cases_text))
  
  #Do an adjusted cumulative link mixed model (ordinal logistic regression mixed model). This is evaluating a weird string.
  expression_text <- paste0(model, "(", score, " ~", paste(reg_ind_vars, collapse = "+"), ", data=fiji_reg_data, Hess=T, nAGQ=10)")
  model <- eval(parse(text=expression_text))
  #get confidence intervals
  ci <- exp(confint(model))
  #exponentiate to get in terms of OR
  dfraw2 <- model[["coefficients"]][-1:-3] %>%
    exp() %>%
    as.data.frame()
  if (i == 15){
    ci
  }
  dfraw2$lower <- ci[rownames(dfraw2),1]
  dfraw2$upper <- ci[rownames(dfraw2),2]
  colnames(dfraw2) <- c(paste0(score, "mixed_sex_point_est"), paste0(score, "mixed_sex_2.5_est"), paste0(score, "mixed_sex_97.5_est"))
  coef.df <- as.data.frame(t(dfraw2))
  
  # Check random effects are normally distributed
  if (model == "clmm"){
    try(qqnorm(ranef(model)))
  }
  
  return(list(coef.df, model))
}

check_regression <- function(df){
  if (class(df) == "try-error"){
    return(as.data.frame(matrix(ncol=7, nrow=3)))
  } else {
    return(df)
  }
}
# 
# control_vars <- c("age", "disability3", "sector", "itaukei", "indian")
# test <- try(regression(score="score15", model='clm', model_name='clmm_sex', ivs=append(control_vars, "sex")))
# 

regress_po <- function(data=fiji_cleaned, score , interaction = TRUE) {
  # Get models
  if (interaction) {
    po_model <- vglm(eval(parse(text = score)) ~ (sex * rainfallperc) + rural + urban + age, data = data, 
                     family=cumulative(parallel=TRUE, reverse = TRUE))
    gen_model <- vglm(eval(parse(text = score)) ~ (sex * rainfallperc) + rural + urban + age, data = data, 
                      family=cumulative(parallel=FALSE, reverse = TRUE))
  } else {
    po_model <- vglm(eval(parse(text = score)) ~ sex + rainfallperc + rural + urban + age, data = data, 
                     family=cumulative(parallel=TRUE, reverse = TRUE))
    gen_model <- vglm(eval(parse(text = score)) ~ sex + rainfallperc + rural + urban + age, data = data, 
                      family=cumulative(parallel=FALSE, reverse = TRUE))
  }
  
  # Compare models
  results <- lrtest(gen_model, po_model)
  result <- results@Body[2, "Pr(>Chisq)"]
  
  # Throw warning for non-poa validate
  if (result <= 0.05) warning(paste0("Regression on ", score, " does not meet the proportional odds assumption with a 
                                     chisq result of ", result))
  
  # Get a clm and test clm
  # clm_poa <- clm(eval(parse(text=score)) ~ (sex * rainfallperc) + rural + urban + age, data = data)
  nomtest <- nominal_test(clm_poa)
  not_validated <- na.omit(rownames(nomtest)[nomtest$`Pr(>Chi)` <= 0.05])
  if ("rural" %in% not_validated) c(not_validated, "urban")
  if ("sex:rainfallperc" %in% not_validated) str_replace(not_validated, "sex:rainfallperc", "sex*rainfallperc")
  
  if (interaction) {
    ppo_expression <- paste0("vglm(eval(parse(text = score)) ~ (sex * rainfallperc) + rural + urban + age, data = data, 
                             family=cumulative(parallel=FALSE ~ ", paste(not_validated, collapse = " + "), ",reverse = TRUE))")
  } else {
    ppo_expression <- paste0("vglm(eval(parse(text = score)) ~ sex + rainfallperc + rural + urban + age, data = data, 
                             family=cumulative(parallel=FALSE ~ ", paste(not_validated, collapse = " + "), ",reverse = TRUE))")
  }
  
  
  ppo_model <- eval(parse(text=ppo_expression))
  resultsppo <- lrtest(gen_model, ppo_model)
  resultppo <- resultsppo@Body[2, "Pr(>Chisq)"]
  if (resultppo <= 0.05) warning(paste0("PPO regression on ", score, " does not meet the proportional odds assumption with a 
                   chisq result of ", resultppo))
  
  ls <- list()
  ls[["po_model"]] <- po_model
  ls[["gen_model"]] <- gen_model
  ls[["ppo_model"]] <- ppo_model
  ls[["poa_test"]] <- result
  
  return(ls)
}
