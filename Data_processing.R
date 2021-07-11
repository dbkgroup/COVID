library(tidyverse)
library(fANCOVA)


# This function plots the time drift normalization for a compound  
# use_median: defines the use of mean or median as normalization objective
# deg_freedom: feed into loess.as degree of freedom

normalize_time_drift_plot_corrections <- function(df, qc_names, b, use_median = T, deg_freedom = 2) {
  ensure_positive_area <- function(x) if_else(x<0, 0, x)
  
  res <- df %>% dplyr::select(c(Sample, Run)) %>% as_tibble()
  
  comp <- colnames(df %>% select(-c(Sample, Run, Group)))
  bb_qc <- df %>% filter(Sample %in% qc_names)
  qc_run <- bb_qc$Run
  
  for(c in comp) {
    qc_areas <- pull(bb_qc, !!c)
    df_med <- df %>% filter(Group %in% c("Sample", "QC"))
    gmed <- median(pull(df_med, !!c))
    if(!use_median)
      gmed <- mean(pull(df_med, !!c))
    
    fit <- loess.as(qc_run, qc_areas, criterion="aicc",degree = deg_freedom,plot = T,control=loess.control(surface="direct"))
    
    xpred <- predict(fit,df$Run)
    norm <- ensure_positive_area(pull(df, !!c) + gmed  - xpred)
    
    # visuals 
    print(fit)
    
    temp_df2 <- data.frame(normalized = norm, original = pull(df, !!c), run = df$Run, pred = xpred, median = gmed) %>% 
      mutate(QC = run %in% qc_run)
    plot_corr <- ggplot(temp_df2, aes(x= run)) + 
      geom_point(aes(y = normalized, color = "Normalized", shape = QC)) +
      geom_point(aes(y = original, color = "Original", shape = QC)) +
      geom_line(aes(y = pred, color = "Correction")) +
      geom_line(aes(y = median, color = "Median")) +
      labs(title = paste("Normalization for", b, c), x = "Run order", y = "Normalization correciton") +
      scale_colour_manual(breaks = c("Normalized", "Original", "Correction", "Median", "QC"), values = c("#F8766D", "#00BFC4", "#FF6C91", "black", "green"))+ 
      scale_shape_manual(breaks = c("TRUE", "FALSE"), values=c(19, 1)) +
      theme_pubr(base=9, base_family = "sans", border = F) +
      theme(plot.title = element_text(hjust = 0.5), axis.text = element_text(color = "gray20", size = 8)) + 
      theme(legend.direction = 'vertical', legend.position = 'right')
    print(plot_corr)
    
    #end of visuals
    
    temp_bb_norm <- tibble({{c}} := norm)
    res <- res %>% cbind(temp_bb_norm)
  }
  res
}

# This function plots the time drift normalization for a compound per batch 
# use_median: defines the use of mean or median as normalization objective
# deg_freedom: feed into loess.as degree of freedom
covid_normalize_by_batch_plot_correction <- function(eset, use_median = T, deg_freedom = 2) {
  
  #batch loop
  data_norm2 <- tibble()
  batches <- unique(eset@phenoData$Batch)
  #batches <- "B2"
  for(b in batches) {
    eset2 <- eset[,eset@phenoData$Batch==b]
    #create a tibble with areas, run # and sample names
    
    dataT <- as.data.frame(t(exprs(eset2))) %>%
      rownames_to_column(var = "Sample") %>%
      left_join(pData(eset2) %>% select(c(Run, Sample, Group)))
    
    
    qc_names <- (dataT %>% filter(str_detect(Sample, "QC")))$Sample
    batch_norm <- normalize_time_drift_plot_corrections(dataT, qc_names, b, use_median, deg_freedom)
    
    data_norm2 <- data_norm2 %>% rbind(batch_norm)
  }
  
  #rownames(data_norm2) <- data_norm2$Sample
  xData <- data_norm2 %>% select(-Run) %>% column_to_rownames("Sample") %>% as.matrix()
  xData2 <- t(xData)
  
  pD <- pData(eset) %>%
    arrange(factor(Sample, levels = colnames(xData2)))
  rownames(pD) <- pD$Sample
  
  #debug code
  #print("In xData but not in pData")
  #print(setdiff(colnames(xData2), pD$Sample))
  #print("In pData but not in xData")
  #print(setdiff(pD$Sample, colnames(xData2)))
  #end of debug code
  
  eset2 <- ExpressionSet(assayData = xData2, phenoData = AnnotatedDataFrame(pD), featureData = AnnotatedDataFrame(fData(eset)))
  
  eset2
}

# This funciton take a data.frame or tibble containing
# areas, 
# Run:  run_order
# and vector with QC sample names (sample to normalize on)
normalize_time_drift <- function(df, qc_names) {
  ensure_positive_area <- function(x) if_else(x<0, 0, x)
  
  res <- df %>% select(c(Sample, Run)) %>% as_tibble()
  
  comp <- colnames(df %>% select(-c(Sample, Run, Group)))
  bb_qc <- df %>% filter(Sample %in% qc_names)
  qc_run <- bb_qc$Run
  
  #print("1.1")
  #print(length(comp))
  
  mc = length(comp)
  pb <- txtProgressBar(min = 0, max = mc, style = 3, width = 50, char = "=")
  i <- 0
  
  for(c in comp) {
    i <- i + 1
    setTxtProgressBar(pb, i)
    
    qc_areas <- pull(bb_qc, !!c)
    df_med <- df %>% filter(Group %in% c("Sample", "QC"))
    gmed <- median(pull(df_med, !!c))
    fit <- loess.as(qc_run, qc_areas, criterion="aicc",degree = 2,plot = F,control=loess.control(surface="direct"))
    
    xpred <- predict(fit,df$Run)
    
    norm <- ensure_positive_area(pull(df, !!c) + gmed  - xpred)
    temp_bb_norm <- tibble({{c}} := norm)
    res <- res %>% cbind(temp_bb_norm)
  }
  res
}

# This fucntion take eSet of covid formatted data and normalizes for QC based time drift per batch
# Expects a Batch column in phenoData
# Expects that the sample list contains QC 
covid_normalize_by_batch <- function(eset, use_median = T, ged_freedom = 2) {
  
  #batch loop
  data_norm2 <- tibble()
  batches <- unique(eset@phenoData$Batch)
  #batches <- "B2"
  
  mc = length(batches)
  pb <- txtProgressBar(min = 0, max = mc, style = 3, width = 50, char = "=")
  i <- 0
  
  for(b in batches) {
    i <- i + 1
    setTxtProgressBar(pb, i)
    
    eset2 <- eset[,eset@phenoData$Batch==b]
    #create a tibble with areas, run # and sample names
    
    dataT <- as.data.frame(t(exprs(eset2))) %>%
      rownames_to_column(var = "Sample") %>%
      left_join(pData(eset2) %>% select(c(Run, Sample, Group)))
    
    qc_names <- pData(eset2) %>% filter(Group == "QC") %>% pull(Sample)
    qc_names2 <- (dataT %>% filter(str_detect(Sample, "QC")))$Sample
    
    batch_norm <- normalize_time_drift(dataT, qc_names)
    data_norm2 <- data_norm2 %>% rbind(batch_norm)
    
  }
  
  #rownames(data_norm2) <- data_norm2$Sample
  xData <- data_norm2 %>% select(-Run) %>% column_to_rownames("Sample") %>% as.matrix()
  xData2 <- t(xData)
  
  
  pD <- pData(eset) %>%
    arrange(factor(Sample, levels = colnames(xData2)))
  rownames(pD) <- pD$Sample
  
  eset2 <- ExpressionSet(assayData = xData2, phenoData = AnnotatedDataFrame(pD), featureData = AnnotatedDataFrame(fData(eset)))
  
  eset2
}

# Constant normalization between batches
# df : data.frame or tibble with areas, Samples names and batches names (it's expected that samples and batch names are not numeric)
normalize_between_batches <- function(df, qc_names) {
  ensure_positive_area <- function(x) if_else(x<0, 0, x)
  
  df1 <- df %>% drop_na()
  
  #Get mean by compound across batches for IQA samples
  iqa <- df %>% filter(Sample_ori %in% qc_names) 
  iqa_mean <- iqa %>% 
    mutate(across(where(is.numeric), ~mean(.x))) %>%
    select_if(is.numeric) %>%
    head(1)
  
  batches <- unique(df$Batch)
  data_norm_res2 <- data.frame()
  
  #print(batches)
  #print(iqa_mean)
  
  for(b in batches) {
    #print(b)
    bb <- df %>% filter(Batch == b) 
    #print(bb)
    #IQA mean in the batch
    bb_iqa_mean <- bb %>% filter(Sample_ori %in% qc_names) %>%
      mutate(across(where(is.numeric), ~mean(.x))) %>%
      select_if(is.numeric) %>%
      head(1)
    
    #correction to be applied in this batch (verctorized for all compounds)
    
    bb_corr <- iqa_mean - bb_iqa_mean
    
    
    bbc <- bb %>% 
      mutate(ID = paste(Batch, Sample, sep = "_"))
    
    rownames(bbc) <- bbc$ID
    
    
    bbc %<>% select_if(is.numeric)
    
    
    
    bb_corr2 <- data.frame()
    
    #create correction data.frame of the same size as batch samples to allow for vectorized operation
    for(i in 1:nrow(bbc)){
      bb_corr2 <- bb_corr2 %>% rbind(bb_corr)
    }
    #apply correction to all samples
    data_norm_res2 <- data_norm_res2 %>% rbind(bbc + bb_corr2)
  }
  
  data_norm_res22 <- data_norm_res2 %>% 
    rownames_to_column(var = "Sample") %>%
    mutate(across(where(is.numeric), ensure_positive_area))
  
  data_norm_res22
}


# This function normalized between batches
covid_normalize_between_batches <- function(eset) {
  
  dataT <- as.data.frame(t(exprs(eset))) %>%
    rownames_to_column(var = "Sample_ori") %>%
    separate(Sample_ori, c("Batch", "Sample"), sep = "_", extra = "merge", remove = F)
  
  qc_names <- pData(eset) %>% filter(Group == "IQA") %>% pull(Sample)
  #qc_names2 <- (dataT %>% filter(startsWith(Sample, "IQA")) %>% filter(!Sample %in% c("IQA50", "IQA75", "IQA125", "IQA150")))$Sample
  
  res <- normalize_between_batches(dataT, qc_names)
  
  xData <- res %>% column_to_rownames("Sample") %>% as.matrix()
  xData2 <- t(xData)
  eset2 <- ExpressionSet(assayData = xData2, phenoData = AnnotatedDataFrame(pData(eset)), featureData = AnnotatedDataFrame(fData(eset)))
  
  eset2
}





