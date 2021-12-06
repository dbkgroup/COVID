library(gtools)

# This function splits the data into train / test at proportion 'prop' and balanced on a label 'balance'
get_train_rows <- function(data, prop = 0.8, balance =as.character(NA)) {
  
  train_rows <- NA  
  N = nrow(data)
  if(!is.na(balance)) {
    data2 <- data %>% rownames_to_column("RowNb") %>% filter(pull(data, !!balance) == F)
    n <- nrow(data2)
    train_rows_temp <- sample(n,size=n*prop)
    train_rows_F <- data2[train_rows_temp,]$RowNb
    
    data2 <- data %>% rownames_to_column("RowNb") %>% filter(pull(data, !!balance) == T)
    n <- nrow(data2)
    train_rows_temp <- sample(n,size=n*prop)
    train_rows_T <- data2[train_rows_temp,]$RowNb
    
    train_rows <- sort(as.numeric(c(train_rows_F, train_rows_T)))
  }
  else {
    print("No balancing")
    train_rows<-sample(N,size=N*prop)
  }
  
  return(train_rows)
}

# This function splits the data into train / test at proportion 'prop' and balanced on a 2 labels simultaneously 'balance1' and 'balance2'
get_train_rows2 <- function(data, prop = 0.8, balance1 =as.character(NA), balance2 =as.character(NA)) {
  
  train_rows <- NA  
  N = nrow(data)
  if(!is.na(balance1) & !is.na(balance2)) {
    data2 <- data %>% rownames_to_column("RowNb") %>% filter(pull(data, !!balance1) == F & pull(data, !!balance2) == F)
    n <- nrow(data2)
    train_rows_temp <- sample(n,size=n*prop)
    train_rows_F <- data2[train_rows_temp,]$RowNb
    
    data2 <- data %>% rownames_to_column("RowNb") %>% filter(pull(data, !!balance1) == T & pull(data, !!balance2) == F)
    n <- nrow(data2)
    train_rows_temp <- sample(n,size=n*prop)
    train_rows_TF <- data2[train_rows_temp,]$RowNb
    
    data2 <- data %>% rownames_to_column("RowNb") %>% filter(pull(data, !!balance1) == T & pull(data, !!balance2) == T)
    n <- nrow(data2)
    train_rows_temp <- sample(n,size=n*prop)
    train_rows_TT <- data2[train_rows_temp,]$RowNb
    
    train_rows <- sort(as.numeric(c(train_rows_F, train_rows_TF, train_rows_TT)))
  }
  else {
    print("No balancing")
    train_rows<-sample(N,size=N*prop)
  }
  
  return(train_rows)
}

# This function is used in the Monte-Carlo accuracy for rstan model
build_metrics_stan2<-function(data, model) {
  res <- posterior_predict(model, newdata = (data %>% select(-Label)))
  res_mean <- apply(res, MARGIN = 2, mean)
  prior <- 0.5
  levels <- c(T,F)
  data.frame(pred = res_mean,
             truth = factor(data$Label, levels=levels),
             pred_class = factor(res_mean>prior, levels=levels)) 
}

# Get Monte-Carlo accuracy for rstan model
get_stan_test_accuracy <- function(data, weight = NA) { 
  
  tt <- get_train_rows(data, 0.8, "Label")
  data_train <- data[tt,]
  data_test <- data[-tt,]
  
  options(mc.cores = parallel::detectCores())
  if(!is.na(weight)) {
    sglm_r <- stan_glm(Label ~ ., family = binomial(), data = data_train, prior = normal(scale = 1), weights = weight)
    
  }
  else {
    sglm_r <- stan_glm(Label ~ ., family = binomial(), data = data_train, prior = normal(scale = 1))
  }
  
  test_metrics <- build_metrics_stan2(data_test, sglm_r)
  sum <- test_metrics %>% conf_mat(truth,pred_class) %>% summary
  accuracy <- sum %>% filter(.metric == "accuracy") %>% pull(.estimate)
  bal_accuracy <- sum %>% filter(.metric == "bal_accuracy") %>% pull(.estimate)
  auc <- test_metrics %>% roc_auc(truth, pred) %>% pull(.estimate)
  
  results <- data.frame(Metric = c("accuracy", "bal_accuracy", "AUC"), Value = c(accuracy, bal_accuracy, auc))
  
  return(results)
}

# Get Monte-Carlo accuracy for rstan model
get_mc_stan_accuracy <- function(nb_iter, data, weight = NA) {
  
  pb <- txtProgressBar(min = 0, max = nb_iter, style = 3, width = 50, char = "=") 
  set.seed(256)
  res <- data.frame()
  for(i in 1:nb_iter) {
    setTxtProgressBar(pb, i)
    res <- res %>% rbind(get_stan_test_accuracy(data, weight))
  }
  return(res %>% group_by(Metric) %>% dplyr::summarise(mean = mean(Value), var = var(Value), sd = sd(Value)))
}



# This function return a results for logistic regression adjusted by all columns in 'corrected_for'
get_compound_logistic_reg <- function(eset, condition, compound, corrected_for = c()) {
  
  set.seed(4569)
  eset <- eset[eset@featureData$Compound == compound,]
  
  p_retain <- append(corrected_for, c("Sample", {{condition}}))
  df <- as.data.frame(t(exprs(eset))) %>%
    dplyr::rename(Area := {{compound}}) %>%
    rownames_to_column(var = "Sample") %>%
    left_join(pData(eset) %>% select(p_retain)) %>%
    dplyr::rename(Condition := {{condition}}) %>% 
    drop_na()
  
  df_z <- df %>% mutate_if(is.numeric, scale) %>%
    mutate_if(is.character, factor)
  
  if(length(corrected_for) > 0)
    fmla_t <- "Condition ~ 0 + Area"
  else
    fmla_t <- "Condition~Area"
  for(c in corrected_for) {
    fmla_t <- paste(fmla_t, "+", c, sep = " ")
  }
  fmla <- as.formula(fmla_t)
  
  logit_model<-glm(fmla,data=df_z,family=binomial(link='logit'))
  if(logit_model$converged==T ) {
    
    if(length(corrected_for) > 0) {
      or_res <- exp(coef(logit_model))[[1]]
      ci_res <- exp(confint(logit_model))
      ci_2_5 <- ci_res[1, '2.5 %']
      ci_97_5 <- ci_res[1, '97.5 %']
      
      coefs <- coef(summary(logit_model))
      p_value <- coefs[1,'Pr(>|z|)']
      
    }
    else{
      or_res <- exp(coef(logit_model))[[2]]
      ci_res <- exp(confint(logit_model))
      ci_2_5 <- ci_res[2, '2.5 %']
      ci_97_5 <- ci_res[2, '97.5 %']
      
      coefs <- coef(summary(logit_model))
      
      if(nrow(coefs)==2)
        p_value <- coefs[2,'Pr(>|z|)']
    }
    return (data.frame("Compound" = compound, "OR_97.5CI" = paste(format(or_res, digits = 2), " (", format(ci_2_5, digits = 2), "-", format(ci_97_5, digits = 2), ")", sep = ""), "P_value" =  stars.pval(p_value)))
  }
  return (data.frame("Compound" = compound, "OR(97.5%CI)" = 0, "P_value" =  0, "P_value val" = 0))
}

get_compound_stan_logistic_reg <- function(eset, condition, compound, corrected_for = c()) {
  
  set.seed(4569)
  eset <- eset[eset@featureData$Compound == compound,]
  
  p_retain <- append(corrected_for, c("Sample", {{condition}}))
  df <- as.data.frame(t(exprs(eset))) %>%
    dplyr::rename(Area := {{compound}}) %>%
    rownames_to_column(var = "Sample") %>%
    left_join(pData(eset) %>% select(p_retain)) %>%
    dplyr::rename(Condition := {{condition}}) %>% 
    drop_na()
  
  df_z <- df %>% mutate_if(is.numeric, scale) %>%
    mutate_if(is.character, factor)
  
  if(length(corrected_for) > 0)
    fmla_t <- "Condition ~ 0 + Area"
  else
    fmla_t <- "Condition~Area"
  for(c in corrected_for) {
    fmla_t <- paste(fmla_t, "+", c, sep = " ")
  }
  fmla <- as.formula(fmla_t)
  
  logit_model<-stan_glm(fmla,data=df_z,family=binomial(link='logit'))
  
    if(length(corrected_for) > 0) {
      or_res <- exp(coef(logit_model))[[1]]
      ci_res <- exp(posterior_interval(logit_model, prob = 0.95))
      ci_2_5 <- ci_res[1, '2.5%']
      ci_97_5 <- ci_res[1, '97.5%']
      
    }
    else{
      or_res <- exp(coef(logit_model))[[2]]
      ci_res <- exp(posterior_interval(logit_model, prob = 0.95))
      ci_2_5 <- ci_res[2, '2.5%']
      ci_97_5 <- ci_res[2, '97.5%']
      
    }
    return (data.frame("Compound" = compound, "OR_97.5CI" = paste(format(or_res, digits = 2), " (", format(ci_2_5, digits = 2), "-", format(ci_97_5, digits = 2), ")", sep = "")))
  
}
