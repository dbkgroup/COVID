
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

