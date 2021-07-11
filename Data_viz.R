#To install ggbiplot, uncomment line 3 and 4 and run only once, comment them back once done

#library(devtools)
#install_github("vqv/ggbiplot")
library(ggbiplot)

library(ggpubr)
library(gridExtra)
library(ggrepel)
library(RColorBrewer)


# Basic PCA plot function based on ggbiplot
i_plotPCA <- function (pca, title, groups = NULL, colorsName = NULL, colors = NULL, ellipse = F, circle = F, labels = NULL, alpha = 1, title_center = T) {
  
  plot <- ggbiplot(pca, obs.scale = 1, var.scale = 1, alpha = alpha,
                   groups = as.factor(groups), ellipse = ellipse, circle = circle, var.axes = F, labels = labels) + 
    ggtitle(title) +
    theme_pubr(base=9, base_family = "sans", border = F) +
    theme(legend.direction = 'horizontal', legend.position = 'none')
  
  if(! is.null(colors))
    plot <- plot + scale_color_manual(name = colorsName, breaks = colors$breaks, values = colors$values) 
  
  if(title_center) 
    plot <- plot + theme(plot.title = element_text(hjust = 0.5), axis.text = element_text(color = "gray20", size = 8))
  
  return(plot)
}


# Overload of i_plotPCA with focus on LC-MS batches
get_pca_by_batch <- function(eset, group_filter, labels = F, title = "By batch", title_center = T) {
  
  eset2 <- eset[,eset@phenoData$Group %in% group_filter]
  eset2 <- eset2[,!(eset2@phenoData$Sample_name %in% c("IQA50", "IQA75", "IQA125", "IQA150", "SQA50", "SQA75", "SQA125", "SQA150"))]
  
  pca_pos_a <- prcomp(t(exprs(eset2)), scale. = TRUE)
  
  #PCAcolors <- data.frame(breaks = c("B1", "B2", "B3", "B4", "B5"), 
  #                         values = c("midnightblue", "lightcoral", "gray50", "mediumaquamarine", "dodgerblue"))
  
  pD <- pData(eset2)
  col <- pull(pD, Batch)
  PCAcolors <- data.frame(breaks = unique(col), 
                          values = brewer.pal(length(unique(col)), "Dark2")[1:length(unique(col))])
  
  if(labels) 
    p1 <- i_plotPCA(pca_pos_a, title, groups = eset2@phenoData$Batch, colorsName = "Sample type", colors = PCAcolors, alpha = 0.7, labels = eset2@phenoData$Sample_name, title_center = title_center)
  
  else 
    p1 <- i_plotPCA(pca_pos_a, title, groups = eset2@phenoData$Batch, colorsName = "Sample type", colors = PCAcolors, alpha = 0.7, title_center = title_center)
  
  p1 <- p1 + theme(legend.direction = 'horizontal', legend.position = 'top', legend.title = element_blank())
  
  p1
}

# Overload of i_plotPCA with focus on sample groups (Sample, QC, IQA, SQA, Blanks)
get_pca_by_group <- function(eset, group_filter, labels = F, title = "By group", title_center = T) {
  
  
  
  eset2 <- eset[,eset@phenoData$Group %in% group_filter]
  eset2 <- eset2[,!(eset2@phenoData$Sample_name %in% c("IQA50", "IQA75", "IQA125", "IQA150", "SQA50", "SQA75", "SQA125", "SQA150"))]
  
  pca_pos_a <- prcomp(t(exprs(eset2)), scale. = TRUE)
  
  #PCAcolors <- data.frame(breaks = c("QC", "Sample", "Blank", "IQA", "SQA"), 
  #                         values = c("mediumaquamarine", "lightcoral", "gray50", "midnightblue", "dodgerblue"))
  
  pD <- pData(eset2)
  col <- pull(pD, Group)
  uni <- unique(col)
  
  #return()
  PCAcolors <- data.frame(breaks = uni, 
                          values = brewer.pal(length(uni), "Dark2")[1:length(uni)])
  
  if(labels) 
    p1 <- i_plotPCA(pca_pos_a, title, groups = eset2@phenoData$Group, colorsName = "Sample type", colors = PCAcolors, alpha = 0.7, labels = eset2@phenoData$Sample_name, title_center = title_center)
  else
    p1 <- i_plotPCA(pca_pos_a, title, groups = eset2@phenoData$Group, colorsName = "Sample type", colors = PCAcolors, alpha = 0.7, title_center = title_center)
  
  p1 <- p1 + theme(legend.direction = 'horizontal', legend.position = 'top', legend.title = element_blank())
  
  p1
}


# Overload of i_plotPCA with focus on whatever desired column in the Expression_set@phenoData
get_pca_by_custom <- function(eset, group_filter, by, labels = F, title = as.character(NA), title_center = T, breaks_list = c()) {
  
  if(is.na(title))
    title <- paste("By", by, sep = " ")
  
  eset2 <- eset[,eset@phenoData$Group %in% group_filter]
  eset2 <- eset2[,!(eset2@phenoData$Sample_name %in% c("IQA50", "IQA75", "IQA125", "IQA150", "SQA50", "SQA75", "SQA125", "SQA150"))]
  
  pca_pos_a <- prcomp(t(exprs(eset2)), scale. = TRUE)
  
  pD <- pData(eset2)
  col <- pull(pD, {{by}})
  brs <- unique(col)
  if(length(breaks_list) > 0)
    brs <- breaks_list
  
  PCAcolors <- data.frame(breaks = brs, 
                          values = brewer.pal(length(unique(col)), "Dark2")[1:length(brs)])
  
  
  if(labels) 
    p1 <- i_plotPCA(pca_pos_a, title = title, groups = col, colorsName = "Sample type", colors = PCAcolors, alpha = 0.8, labels = eset2@phenoData$Sample_name, title_center = title_center)
  else
    p1 <- i_plotPCA(pca_pos_a, title = title, groups = col, colorsName = "Sample type", colors = PCAcolors, alpha = 0.8, title_center = title_center)
  
  p1 <- p1 + theme(legend.direction = 'horizontal', legend.position = 'top', legend.title = element_blank())
  
  p1
}

# Volcano plot, expects p-values and fold change info
# use_q_value: if True q-values will be used, if False p-values will be used
get_volcano_plot <- function(data_fold_pvalue, x_max = NA, add_lables = T, foldChange_lable_lim = 0.5, p_value_lable_lim = 0.05, use_q_value = T) {
  
  y_label <- expression(-log[10]("P_value"))
  data_fold_pvalue_clean <- data_fold_pvalue
  if(use_q_value) {
    data_fold_pvalue_clean %<>% mutate(P_value = Q_value)
    y_label <- expression(-log[10]("Q_value"))
  }
  
  data_fold_pvalue_clean %<>% filter(is.finite(FoldChangeLog2) & !is.na(FoldChangeLog2) & is.finite(P_value) & !is.na(P_value))
  
  data_fold_pvalue_clean$Up_or_Down = case_when(
    (data_fold_pvalue_clean$P_value <=0.05 & data_fold_pvalue_clean$FoldChangeLog2 >= 0.5) ~ "U",
    (data_fold_pvalue_clean$P_value <= 0.05 & data_fold_pvalue_clean$FoldChangeLog2  <= -0.5) ~ "D",
    (data_fold_pvalue_clean$P_value >= 0.05 | data_fold_pvalue_clean$FoldChangeLog2  <= 0.5) ~ "NC",
    TRUE ~ as.character(NA))
  
  p <- ggplot(data_fold_pvalue_clean, aes(x = FoldChangeLog2, y = -log10(P_value), color = Up_or_Down)) +
    geom_point() +
    xlab(expression(log[2]("Fold Change"))) + ylab(y_label) +
    
    scale_colour_manual(values=c("blue", "grey","red"),name=" ",breaks=c("D", "NC", "U"),labels=c("Down", "No Change", "Up")) +
    geom_vline(
      xintercept = c(-0.5,0.5),
      col = "red",
      linetype = "dotted",
      size = 1) +
    geom_hline(
      yintercept = c(-log10(0.01),-log10(0.05)),
      col = "red",
      linetype = "dotted",
      size = 1)+
    theme_bw() +
    theme(legend.position = "none")+
    theme_pubr()+
    theme(panel.grid.minor = element_line(colour = "black", linetype = "dotted")) +
    theme(panel.grid.major = element_line(colour = "black", linetype = "dotted")) +
    theme(panel.background = element_rect(colour = "black")) +
    theme(axis.title.x = element_text( face="bold",size=16)) +
    theme(axis.title.y = element_text(face="bold",size=16)) +
    theme(axis.text.x=element_text(face="bold",colour='black', size=14)) +
    theme(axis.text.y=element_text(face="bold",colour='black', size=14))
  
  if(add_lables) {
    p <- p +
      geom_text_repel(data=subset(data_fold_pvalue_clean,abs(FoldChangeLog2) >= foldChange_lable_lim & P_value < p_value_lable_lim),
                      aes(FoldChangeLog2, -log10(P_value), label = Compound),size = 3, fontface = "bold", color="steelblue")
  }
  
  if(!is.na(x_max)) {
    p <- p +
      scale_x_continuous(limits = c(-x_max, x_max))
  }
  p
}

# Visulaize a compound in boxplot fashion
# eset: Expressionset
# compound: compound ID m/z_RT example "174.11151_0.679" for Arginine
# by: boxplot x axis , must be a column name present in eset@phenoData
# filter_na: if True remove NA in 'by' column if False NAs will apear as a box in the plot
# y_log: if True - log transform areas
getBoxPlotForCompound <- function(eset, compound, by, filter_na = T, subtitle = "", y_log = F) {
  
  eset2 <- eset[,eset@phenoData$Group == "Sample"]
  eset2 <- eset2[eset2@featureData$Compound == compound,]
  
  df <- as.data.frame(t(exprs(eset2))) %>%
    dplyr::rename(Area := {{compound}}) %>%
    rownames_to_column(var = "Sample") %>%
    left_join(pData(eset2) %>% select(c(Sample, {{by}}))) %>%
    dplyr::rename(Label := {{by}}) 
  
  if(filter_na) 
    df %<>% filter(!is.na(Label))
  
  
  gg <- ggplot(data = df, aes(x = Label, y = Area)) +
    geom_boxplot(outlier.colour="black", outlier.shape=16, outlier.size=2, notch=FALSE) +
    ggtitle(paste(compound, eset2@featureData$Name), subtitle = subtitle) +
    labs(x = by) + 
    theme_pubr(base=9, base_family = "sans", border = F) +
    theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5),  axis.text = element_text(color = "gray20", size = 8))
  
  if(y_log) {
    gg <- gg + 
      scale_y_log10() +
      labs(y = "Area (log10 scaled)")
  }
  gg
}




