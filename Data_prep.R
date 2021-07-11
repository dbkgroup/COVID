library(tidyverse)
library(magrittr)
library(Biobase)


# This function extracts only the compound information of interest and standardizes the column names. 
# get_normalized_areas: if True the data subset will include Thermo's Compound Discoverer normalized areas, if False raw areas will be retained 
get_areas_and_info <- function(df, get_normalized_areas = T) {
  
  area_by <- "Area:"
  if(get_normalized_areas) area_by <- "Norm."
  
  temp <- df %>% select(starts_with("Norm."))
  
  
  temp2 <- df %>% 
    mutate(RT = df$`RT [min]`) %>%
    mutate(ChemSpider = df$`# ChemSpider Results`) %>%
    mutate(mzCloud = df$`# mzCloud Results`) %>%
    mutate(mzcloudBest = df$`mzCloud Best Match`) %>%
    mutate(Area_mean = df$`Mean Area`) %>%
    mutate(AreaMax = df$`Area (Max.)`)
  
  CD3.2 <- F
  if("Calc. MW" %in% colnames(df))
    CD3.2 <- T
  
  if(CD3.2) {
    temp2 %<>% mutate(MW = df$`Calc. MW`)
    temp2 %<>% mutate(Predicted_Composition = df$`Annot. Source: Predicted Compositions`)
  }
  else {
    temp2 %<>% mutate(MW = df$`Molecular Weight`)
    temp2 %<>% mutate(Predicted_Composition = df$`Annotation Source: Predicted Compositions`)
  }
  
  temp2 %<>% 
    mutate(Has_norm = rowSums(is.na(temp))==0) %>%
    mutate(AreaCV = df$`Area CV [%]`) %>%
    mutate(RSD_QC_corr = `RSD Corr. QC Areas [%]`)
  
  subDF <- temp2 %>% dplyr::select(Name, 
                                   Formula, 
                                   MW, 
                                   RT, 
                                   ChemSpider, 
                                   mzCloud, 
                                   mzcloudBest, 
                                   MS2, 
                                   Background, 
                                   AreaCV, 
                                   Has_norm,
                                   Area_mean,
                                   RSD_QC_corr,
                                   Predicted_Composition,
                                   AreaMax,
                                   starts_with(area_by))
  
  return(subDF)
} 

# This function contains extracts the basic info about a run from Thermo's Compound Discoverer 'Compounds table' export
# df: dataFrame containing Thermo's Compound Discoverer 'Compounds table' export
get_basic_stats <- function(df, title = NA) {
  
  if(is.na(title))
    title <- "Unknown"
  
  Criteria <- c("Data file", "Total # compounds", "Background compounds", "Other excluded compounds", "Retained compounds", 
                "Has MS2", "Has Predicted composition match", "Has ChemSpider ID", "Has MassList full match", "Predicted Composition + mzCloud ID >=70", "Predicted Composition + mzValut ID >=75")
  Value <- c(title, 
             nrow(df),
             nrow(filter(df, Background == T)), 
             nrow(filter(df, (Background == F & Norm == F))), 
             nrow(filter(df, (Background == F & Norm == T))), 
             nrow(filter(df, (Background == F & Norm == T & MS2 != "No MS2"))), 
             nrow(filter(df, (Background == F & Norm == T & Predicted_Composition == "Full match"))),
             nrow(filter(df, (Background == F & Norm == T & ChemSpider > 0))), 
             nrow(filter(df, (Background == F & Norm == T & MassList == "Full match"))),
             nrow(filter(df, (Background == F & Norm == T & mzcloudBest >= 70 & Predicted_Composition == "Full match"))),
             nrow(filter(df, (Background == F & Norm == T & mzVaultBest >= 70 & Predicted_Composition == "Full match"))))
  
  
  data.frame(Criteria, Value)
}

# This function selects a subset of basic information from Thermo's Compound Discoverer 'Compounds table' export and renames the column in R friendly names
# df: dataFrame containing Thermo's Compound Discoverer 'Compounds table' export
get_compound_info_short_list <- function(df) {
  
  temp <- df %>% dplyr::select(starts_with('Norm.'))
  
  temp2 <- df %>% 
    mutate(RT = df$`RT [min]`) %>%
    mutate(ChemSpider = df$`# ChemSpider Results`) %>%
    mutate(mzCloud = df$`# mzCloud Results`) %>%
    mutate(mzcloudBest = df$`mzCloud Best Match`) %>%
    mutate(mzVaultBest = df$`mzVault Best Match`) %>%
    mutate(Area = df$`Mean Area`)
  
  
  CD3.2 <- F
  if("Calc. MW" %in% colnames(df)) ### this will now be Calc. MW rather than Acquired MW
    CD3.2 <- T
  
  if(CD3.2) {
    temp2 %<>% mutate(MW = df$`Calc. MW`) ### this will now be Calc. MW rather than Acquired MW
    temp2 %<>% mutate(Predicted_Composition = df$`Annot. Source: Predicted Compositions`)
    temp2 %<>% mutate(MassList = df$`Annot. Source: MassList Search`)
  }
  else {
    temp2 %<>% mutate(MW = df$`Molecular Weight`)
    temp2 %<>% mutate(Predicted_Composition = df$`Annotation Source: Predicted Compositions`)
    temp2 %<>% mutate(MassList = df$`Annotation Source: MassList Search`)
    
  }
  
  temp2 %<>% 
    mutate(Norm = rowSums(is.na(temp))==0) %>%
    mutate(AreaCV = df$`Area CV [%]`) %>%
    mutate(AreaSampleMean = rowMeans(dplyr::select(df, starts_with("Area: SerumMS1")), na.rm = TRUE)) %>%
    mutate(NormAreaSampleMean = rowMeans(dplyr::select(df, starts_with("Norm. Area: SerumMS1")), na.rm = TRUE))
  
  subDF <- temp2 %>% dplyr::select(Name, 
                                   Formula, 
                                   MW, 
                                   RT, 
                                   ChemSpider, 
                                   MassList,
                                   mzCloud, 
                                   mzcloudBest, 
                                   mzVaultBest,
                                   MS2, 
                                   Background, 
                                   AreaCV, 
                                   Norm,
                                   Area,
                                   AreaSampleMean,
                                   NormAreaSampleMean,
                                   Predicted_Composition)
  
  return(subDF)
}

# This function takes compounds areas, compounds information as provided by Compound Discoverer and samples metadata to 
# create and ExpressionSet object (see Biobase documentation)

# compounds_data: extract from Thermo's Compound Discoverer compound table
# sample_info: samples metadata including run order
# normalized_areas: if True the CD normalized areas will be extracted, if False the raw areas will be extracted (when we want to perform custom normalization)
# background: if True only background labeled compounds will be extracted, if False only non-backgroud compounds will be extracted
# normalized: if True only compounds that have fit the normalization criterias of CD workflow will be extracted 
# (QC RSD < 30%, QC decetion > 80%, etc..), if False all compounds will be exctracted and the custom normalizaiton should include those filter (not implemented in this code base) 

covid_transform_to_eset <- function(compounds_data, sample_info, normalized_areas = T, background = F, normalized = T) {
  
  sample_name <- function(name) sub(".*?CV1_S1_(.*?)\\..*", "\\1", name, T)  
  
  area_by <- "Area:"
  if(normalized_areas) area_by <- "Norm."
  
  dff <- get_areas_and_info(compounds_data, normalized_areas) %>% 
    filter(Background == background) %>% # & Has_norm == normalized) %>%
    mutate(Compound = paste(MW, RT, sep = "_")) %>%
    arrange(Compound)
  
  if(normalized)
    dff %<>% filter(Has_norm == T)
  
  
  fData <- dff %>% select(Compound, !starts_with(area_by)) %>% as.data.frame()
  xData <- dff %>% 
    select(starts_with(area_by)) %>%
    rename_with(sample_name) %>%
    as.data.frame()
  
  
  rownames(xData) <- dff$Compound
  rownames(fData) <- dff$Compound
  
  
  pData <- sample_info %>%
    mutate(Group = ifelse(startsWith(Sample, 'QC'), "QC", 
                          ifelse(startsWith(Sample, 'B') | startsWith(Sample, 'GB'), "Blank", 
                                 ifelse(startsWith(Sample, "IQA"), "IQA", 
                                        ifelse(startsWith(Sample, "SQA"), "SQA",
                                               "Sample"))))) %>%
    mutate(Sample_name = Sample) %>%
    mutate(Sample = paste(Batch, Sample, sep = "_")) %>%
    filter(!str_detect(Sample, "DD") & !str_detect(Sample, "dd") & !str_detect(Sample, "CQC")) %>%
    arrange(factor(Sample, levels = colnames(xData))) %>%
    filter(Sample %in% colnames(xData)) %>%
    as.data.frame() # we need this to be a dataframe so we can set rownames, the rownames settings is deprecated on tibbles
  
  rownames(pData) <- pData$Sample 
  
  # ExpressionSet object requires phenoData to match xData columns names (in the same order) and featureData to match xData row names (in the same order)
  # Mismatch names or order in those is the most common error in creating this object 
  eset <- ExpressionSet(assayData = as.matrix(xData), phenoData = AnnotatedDataFrame(pData), featureData = AnnotatedDataFrame(fData))
  return(eset)
}


# This function remove compounds that reach 'tresh' percentage missing in the subgoup of samples defined by group_filter (usually QC samples)
# group_filter: subgoup(s) of samples to be used for the cut_off (usually QC)
# tresh: threshold of tolerance for missing detection percentage (usually 0.2 -> 80% present or 20% missing) 

remove_compounds <- function(eset, group_filter, tresh) {
  
  eset_filter <- eset[,eset@phenoData$Group %in% group_filter]
  
  df <- data.frame(exprs(eset_filter)) %>% rownames_to_column(var = "Compound") %>%
    na_if(0) 
  df %<>% mutate(Absent = apply(is.na(df), 1, sum)) %>%
    select(c(Compound, Absent))
  
  count_tresh <- (length(colnames(exprs(eset_filter)))-2)*tresh
  df_retain <- df %>% filter(Absent < count_tresh)
  qc_retain_list <- df_retain$Compound
  
  eset2 <- eset[eset@featureData$Compound %in% df_retain$Compound,]
  return(eset2)
}


