rm(list=ls())
# server.R

# load useful packages/libraries ----------------------

if(!require("d3heatmap")){
  install.packages("d3heatmap")
  library(d3heatmap)
}

if(!require("DT")){
  install.packages("DT")
  library(DT)
}

if(!require("plotly")){
  install.packages("plotly")
  library(plotly)
}

if(!require("viridis")){
  install.packages("viridis")
  library(viridis)
}

if(!require("dplyr")){
  install.packages("dplyr")
  library(dplyr)
}

# load useful functions ----------------------------
epi_names <- function(de_table){
  if (nrow(de_table) == 0){
    epi_names <- "-"
  } else if (nrow(de_table) > 0){
    epi_names <- de_table$Gene_name }
  return(epi_names)}

no_peaks_in_selected_cluster <- function(cluster_table, inputCluster){
  print(inputCluster)
  if (inputCluster != "none"){
    print("inputCluster != none")
    s <- c(inputCluster, "row.names")
    selected_cluster <- cluster_table[,match(s, colnames(cluster_table))]
    print(head(selected_cluster))
    selected_cluster <- selected_cluster[selected_cluster[,c(1:ncol(selected_cluster)-1)] > 0,]
    print(head(selected_cluster))
    nr_selected_clusters <- nrow(selected_cluster)
    print(nr_selected_clusters)}
  else{
    nr_selected_clusters <- 0 }
  print(paste("Number finally delivered", nr_selected_clusters))
  return(nr_selected_clusters)}

motif_heatmap <- function(ucol_Sums_Total_df,colSums.no_clusters,raw_motif_table_per_peak, name){
  pval.Tot.df <- as.data.frame(ucol_Sums_Total_df)
  max_row=nrow(ucol_Sums_Total_df)
  max_col=ncol(ucol_Sums_Total_df)-1
  for (i in (1:max_row)){
    for (j in (1:max_col)){
      q=ucol_Sums_Total_df[i,j]
      m=colSums.no_clusters[colnames(ucol_Sums_Total_df)[j]]
      k=ucol_Sums_Total_df[i,ncol(ucol_Sums_Total_df)]
      n=nrow(raw_motif_table_per_peak)-m
      x <- phyper(q,m,n,k, lower.tail=FALSE, log.p=FALSE)
      pval.Tot.df[i,j] <- x
    }
  }
  pval.Tot.df$Nrows <- NULL
  write.table(as.data.frame(pval.Tot.df), file=paste0(name, ".txt"), sep="\t", quote=F)
  write.table(as.data.frame(-log10(pval.Tot.df +1)), file=paste0(name, "-log10pval1.txt"), sep="\t", quote=F)
  par(mar=c(1,1,1,1), mfrow=c(1,1))
  #result <- list("pval.Tot.df" = pval.Tot.df, "pval.adj.Tot.df" = pval.adj.Tot.df)
  #print(result)
  return(pval.Tot.df)
}

# load necessary data ----------------------
load('RData/genes2names_tfbs.RData')
load("RData/FoxD3_Shiny_Test.RData")
load("RData/oriTFBS.RData")
load("RData/tfbs_enh.RData")
load("RData/annot.RData")
load("RData/clusters.RData")
load("RData/enh2gene_peaks.RData")
load("RData/cluster_lists.RData")

# Define server logic ------------------------

shinyServer(function(input, output, session){

  
#  RNA-seq Analysis Part -------------------------------------------------------------
#  updateSelectizeInput(session, 'Gene_Name', choices = annot$Associated.Gene.Name, server = TRUE)
  output$Gene_Name_UI <- renderUI({
    selectInput(
      'Gene_Name',
      label = NULL,
      choices = annot$Associated.Gene.Name,
      selected = "foxd5",
      multiple = FALSE
    )
  })
  
  output$SelectedGene <- renderText({ 
    paste("You have selected the gene: ", input$Gene_Name )
  })
  
  EnsId <- reactive({
    ensembl <- annot[ annot$Associated.Gene.Name == input$Gene_Name,"Ensembl.Gene.ID"]
    return(ensembl) 
  })
  
  output$Ensembl_Geneid <- renderText({ 
    paste("The Ensembl Gene Id for that is: ", EnsId() )
  })
  
  GeneDescription <- reactive({
    description <- annot[ annot$Associated.Gene.Name == input$Gene_Name,"Description"]
    return(description) 
  })
  
  output$Full_Gene_Name <- renderText({ 
    paste("Full name of this gene is: ", GeneDescription() )
  })
  
  FPKMs <- reactive({
    fpkms.gene <- subset(fpkm, fpkm$geneid == EnsId())
    return(fpkms.gene)
  })
  
  output$plot <- renderPlotly({
    #plot_ly(FPKMs(), x = colnames(FPKMs()), y = FPKMs()[1,])
    df <-  subset(melted_fpkms, melted_fpkms$geneid==EnsId())
    #print(df)
    
    plot_ly(data=df, x =df$value, y=df$variable, 
            type="bar", 
            marker = list(color=c("green","green","firebrick","firebrick","grey", "grey",
                                  "green","green","firebrick","firebrick","grey", "grey",
                                  "green","green","firebrick","firebrick",
                                  "green","green","firebrick","firebrick","grey", "grey",
                                  "purple"))) %>% layout(margin = list(l=100))
  })
  
  output$event <- renderPrint({
    d <- event_data("plotly_hover")
  })    
  
  
  # Subset log2FC table (C vs negs) ----------------------------
  table_negs_gene <- reactive({
    # If missing input, return to avoid error later in function
    if(is.null(input$Gene_Name))
      return(NULL)
    # subset df
    table_negs_gene <- subset(all_vs_negs.de, all_vs_negs.de$geneid == EnsId())
  })
  
  # Citrine vs negs table ----------------------------
  
  table_negs.f <- reactive({
    # If missing input, return to avoid error later in function
    if(is.null(input$Gene_Name))
      return(NULL)
    # subset df
    table_negs1 <- table_negs_gene()[,c("Epi_vs_negs_log2FoldChange", "Epi_vs_negs_padj", "Significant_at_Epiboly")]
    table_negs2 <- table_negs_gene()[,c("NCC5.6ss_vs_negs_log2FoldChange", "NCC5.6ss_vs_negs_padj", "Significant_at_5.6ss")]      
    table_negs3 <- table_negs_gene()[,c("NCC14ss_vs_negs_log2FoldChange", "NCC14ss_vs_negs_padj", "Significant_at_14ss")]
    colnames(table_negs1) <- c("Log2FoldChange", "Adjusted P-value", "Significant?")
    colnames(table_negs2)  <- c("Log2FoldChange", "Adjusted P-value", "Significant?")
    colnames(table_negs3) <- c("Log2FoldChange", "Adjusted P-value", "Significant?")
    table_negs.f <- rbind(table_negs1,table_negs2,table_negs3)
    table_negs.f$Stage <- c("Epiboly", "5-6ss", "14-16ss")
    table_negs.f <- table_negs.f[,c(4,1,2,3)]
    #print(table_negs.f)
  })
  output$table_negs.f <- renderTable({table_negs.f()})
  
  # Subset log2FC table (CC vs C) ----------------------------
  
  table_mut_gene <- reactive({
    # If missing input, return to avoid error later in function
    if(is.null(input$Gene_Name))
      return(NULL)
    # subset df
    table_mut_gene <- subset(all_cc_vs_c.de, all_vs_negs.de$geneid == EnsId())
  })
  
  
  # CC vs C table ----------------------------
  
  table_mut.f <- reactive({
    # If missing input, return to avoid error later in function
    if(is.null(input$Gene_Name))
      return(NULL)
    # subset df
    table_mut1 <- table_mut_gene()[,c("Epi_CC_vs_C_log2FoldChange", "Epi_CC_vs_C_padj", "Significant_at_Epiboly")]
    table_mut2 <- table_mut_gene()[,c("NCC5.6ss_CC_vs_C_log2FoldChange", "NCC5.6ss_CC_vs_C_padj", "Significant_at_5.6ss")]
    table_mut3 <- table_mut_gene()[,c("NCC12ss_CC_vs_C_log2FoldChange", "NCC12ss_CC_vs_C_padj", "Significant_at_12ss")]
    table_mut4 <- table_mut_gene()[,c("NCC14ss_CC_vs_C_log2FoldChange", "NCC14ss_CC_vs_C_padj", "Significant_at_14ss")]
    colnames(table_mut1) <- c("Log2FoldChange", "Adjusted P-value", "Significant?")
    colnames(table_mut2)  <- c("Log2FoldChange", "Adjusted P-value", "Significant?")
    colnames(table_mut3) <- c("Log2FoldChange", "Adjusted P-value", "Significant?")
    colnames(table_mut4) <- c("Log2FoldChange", "Adjusted P-value", "Significant?")
    table_mut.f <- rbind(table_mut1,table_mut2,table_mut3, table_mut4)
    table_mut.f$Stage <- c("Epiboly", "5-6ss", "12-14ss", "14-16ss")
    table_mut.f <- table_mut.f[,c(4,1,2,3)]
    #print(table_mut.f)
  })
  
  output$table_mut.f <- renderTable({table_mut.f()})
  
  # Subset Promoter peaks --------------------------------------
  
  table_promoter_peaks <- reactive({
    # If missing input, return to avoid error later in function
    if(is.null(input$Gene_Name))
      return(NULL)
    # subset df
    table_promoter_peaks <- subset(Prom_peaks, Prom_peaks$Nearest.PromoterID == EnsId())
    class(table_promoter_peaks)
    table_promoter_peaks <- table_promoter_peaks[,c(1:3)]
  })
  output$table_promoter_peaks <- renderTable({table_promoter_peaks()})
  
  # Subset Enhancer peaks --------------------------------------
  
  table_enhancer_peaks <- reactive({
    
    # If missing input, return to avoid error later in function
    if(is.null(input$Gene_Name))
      return(NULL)

    ## Process genes association with only one gene per peaks 
    
    table_enhancer_peaks1 <- subset(Enh_peaks.singles, Enh_peaks.singles$EnsemblGeneID == EnsId())
    table_enhancer_peaks1 <- table_enhancer_peaks1[,c(2:4)]
    table_enhancer_peaks1$Peak_name <- paste(table_enhancer_peaks1[,1],table_enhancer_peaks1[,2],table_enhancer_peaks1[,3],sep="_")

    clust_corresp1 <- subset(cluster_table, cluster_table$row.names %in% table_enhancer_peaks1$Peak_name)
    
    #nrow(clust_corresp1)
    
    rownames(clust_corresp1) <- clust_corresp1$row.names
    clust_corresp1$row.names <- NULL
    
    if (nrow(clust_corresp1) > 0){
      #print("yep")
      t_clust_corresp1 <- as.data.frame(t(clust_corresp1))
      head(t_clust_corresp1)
      table_enhancer_peaks1$Peak_name <- NULL
      table_enhancer_peaks1$Clusters <- ""
      head(table_enhancer_peaks1)
      for (i in (1:ncol(t_clust_corresp1))){
        which_clust_row_in <- rownames(t_clust_corresp1[t_clust_corresp1[i,] > 0 ,])
        which_clust_row_in_str <- toString(which_clust_row_in)
        table_enhancer_peaks1[i,4] <- which_clust_row_in_str
        }
      }
    
    ## Process genes association with multiple genes per peaks  (set1)
    
    table_enhancer_peaks2 <- subset(Enh_peaks.mults1, Enh_peaks.mults1$EnsemblGeneID == EnsId())
    #print(paste("Mults1 size:", nrow(table_enhancer_peaks2)))
    
    table_enhancer_peaks2 <- table_enhancer_peaks2[,c(2:4)]
    table_enhancer_peaks2$Peak_name <- paste(table_enhancer_peaks2[,1],table_enhancer_peaks2[,2],table_enhancer_peaks2[,3],sep="_")
    
    clust_corresp2 <- subset(cluster_table, cluster_table$row.names %in% table_enhancer_peaks2$Peak_name)
    table_enhancer_peaks_fin <- table_enhancer_peaks1
    
    #nrow(clust_corresp2)
    
    rownames(clust_corresp2) <- clust_corresp2$row.names
    clust_corresp2$row.names <- NULL
    
    if (nrow(clust_corresp2) > 0){
      #print("yep")
      t_clust_corresp2 <- as.data.frame(t(clust_corresp2))
      head(t_clust_corresp2)
      table_enhancer_peaks2$Peak_name <- NULL
      table_enhancer_peaks2$Clusters <- ""
      head(table_enhancer_peaks2)
      for (i in (1:ncol(t_clust_corresp2))){
        which_clust_row_in2 <- rownames(t_clust_corresp2[t_clust_corresp2[i,] > 0 ,])
        which_clust_row_in_str2 <- toString(which_clust_row_in2)
        table_enhancer_peaks2[i,4] <- which_clust_row_in_str2
      }
      table_enhancer_peaks_fin <- cbind(table_enhancer_peaks_fin, table_enhancer_peaks2)
      }
    ## Process genes association with multiple genes per peaks  (set2)
    table_enhancer_peaks3 <- subset(Enh_peaks.mults2, Enh_peaks.mults2$EnsemblGeneID == EnsId())
    table_enhancer_peaks3 <- table_enhancer_peaks3[,c(2:4)]
    table_enhancer_peaks3$Peak_name <- paste(table_enhancer_peaks3[,1],table_enhancer_peaks3[,2],table_enhancer_peaks3[,3],sep="_")
    clust_corresp3 <- subset(cluster_table, cluster_table$row.names %in% table_enhancer_peaks3$Peak_name)
    rownames(clust_corresp3) <- clust_corresp3$row.names
    clust_corresp3$row.names <- NULL
    if (nrow(clust_corresp3) > 0){
      t_clust_corresp3 <- as.data.frame(t(clust_corresp3))
      head(t_clust_corresp3)
      table_enhancer_peaks3$Peak_name <- NULL
      table_enhancer_peaks3$Clusters <- ""
      head(table_enhancer_peaks3)
      for (i in (1:ncol(t_clust_corresp3))){
        which_clust_row_in3 <- rownames(t_clust_corresp2[t_clust_corresp3[i,] > 0 ,])
        which_clust_row_in_str3 <- toString(which_clust_row_in3)
        table_enhancer_peaks3[i,4] <- which_clust_row_in_str3
        }
      table_enhancer_peaks_fin <- unique(cbind(table_enhancer_peaks_fin, table_enhancer_peaks3))
      }
    return(table_enhancer_peaks_fin)
  })
  
  output$table_enhancer_peaks <- renderTable({table_enhancer_peaks()})

  # Subset TFBS in peaks -------------------------------------
  table_tfbs_gene <- reactive({
     list_enh <- table_enhancer_peaks()
     list_enh <- paste(table_enhancer_peaks()$Chrom,table_enhancer_peaks()$Start, table_enhancer_peaks()$End, sep="_")  
     tfbs2 <- subset(tfbs, tfbs$geneid %in% list_enh)
     rownames(tfbs2) <- tfbs2$geneid
     tfbs2$geneid <- NULL
     t_tfbs <- as.data.frame(t(tfbs2))
     t_tfbs$Sums <- rowSums(t_tfbs)
     t_tfbs <- t_tfbs[t_tfbs$Sums > 0,]
     t_tfbs$Sums <- NULL
     table_tfbs_gene <- as.data.frame(t(t_tfbs))
   })
       output$table_tfbs_gene <- renderTable({table_tfbs_gene()})
  
  # # TFBS Heatmap --------------------------------------------
   output$HM <- renderD3heatmap({
     if(nrow(table_tfbs_gene()) < 2)
       return(NULL)
     myPalette <- colorRampPalette(c("white", "red", "black"))(100)
     d3heatmap(table_tfbs_gene(),
               show_grid = FALSE,
               anim_duration = 0,
               colors = myPalette,
               xaxis_font_size = "8px",
               yaxis_font_size = "8px",
               xaxis_height = 150,
               yaxis_width = 150
     )
   })
   

# TFBS Analysis Part -----------------------------------------------------------
  # Print selection ----------------------------- 
  output$SelectedTFBS <- renderText({
    paste("You have selected ", input$tf_select)
  })
  
  # Get tf_id
  tf_id <- reactive({
    tf_id <- corresp[corresp$tf_name == input$tf_select,"tf_id"]
    tf_id <- get(as.character(tf_id))
    print("Got tf_id")
    print(tf_id)
    return(tf_id)
  })
  
  output$Number_Paralogs <- renderText({
    paste("There are ",length(tf_id()) ," paralogs corresponding to this TF family in danRer10. They are:" )
  })
  
  table_paralog_list <- reactive({
    if(is.null(input$tf_select))
      return(NULL)
    table_paralog_list <- as.data.frame(tf_id())
    table_paralog_list$Gene_name <- rownames(table_paralog_list)
    rownames(table_paralog_list) <- NULL
    table_paralog_list <- table_paralog_list[,c(2,1)]
    colnames(table_paralog_list) <- c("Gene_name", "EnsemblID")
    #print("Paralog list:")
    #print(table_paralog_list)
    return(table_paralog_list)
  })
  
  output$table_paralog_list <- renderTable({table_paralog_list()})
  
  # Preparse diff exp. 
  
  epi_up <- reactive({
    epi_up <- all_cc_vs_c.de[all_cc_vs_c.de$Epi_CC_vs_C_log2FoldChange > 0 & all_cc_vs_c.de$Epi_CC_vs_C_padj < 0.05 & !is.na(all_cc_vs_c.de$Epi_CC_vs_C_padj),]
    #print(dim(epi_up))
    return(epi_up)
  })
  epi_dn <- reactive({
    epi_dn <- all_cc_vs_c.de[all_cc_vs_c.de$Epi_CC_vs_C_log2FoldChange < 0 & all_cc_vs_c.de$Epi_CC_vs_C_padj < 0.05 & !is.na(all_cc_vs_c.de$Epi_CC_vs_C_padj),]
    #print(dim(epi_dn))
    return(epi_dn)
  })
  
  ncc5.6ss_up <- reactive({
    ncc5.6ss_up <- all_cc_vs_c.de[all_cc_vs_c.de$NCC5.6ss_CC_vs_C_log2FoldChange > 0 & all_cc_vs_c.de$NCC5.6ss_CC_vs_C_padj < 0.05 & !is.na(all_cc_vs_c.de$NCC5.6ss_CC_vs_C_padj),]
    return(ncc5.6ss_up)
  })
  ncc5.6ss_dn <- reactive({
    ncc5.6ss_dn <- all_cc_vs_c.de[all_cc_vs_c.de$NCC5.6ss_CC_vs_C_log2FoldChange < 0 & all_cc_vs_c.de$NCC5.6ss_CC_vs_C_padj < 0.05 & !is.na(all_cc_vs_c.de$NCC5.6ss_CC_vs_C_padj),]
    return(ncc5.6ss_dn)
  })
  
  ncc12ss_up <- reactive({
    ncc12ss_up <- all_cc_vs_c.de[all_cc_vs_c.de$NCC12ss_CC_vs_C_log2FoldChange > 0 & all_cc_vs_c.de$NCC12ss_CC_vs_C_padj < 0.05 & !is.na(all_cc_vs_c.de$NCC12ss_CC_vs_C_padj),]
    return(ncc12ss_up)
  })
  ncc12ss_dn <- reactive({
    ncc12ss_dn <- all_cc_vs_c.de[all_cc_vs_c.de$NCC12ss_CC_vs_C_log2FoldChange < 0 & all_cc_vs_c.de$NCC12ss_CC_vs_C_padj < 0.05 & !is.na(all_cc_vs_c.de$NCC12ss_CC_vs_C_padj),]
    return(ncc12ss_dn)
  })
  
  ncc14ss_up <- reactive({
    ncc14ss_up <- all_cc_vs_c.de[all_cc_vs_c.de$NCC14ss_CC_vs_C_log2FoldChange > 0 & all_cc_vs_c.de$NCC14ss_CC_vs_C_padj < 0.05 & !is.na(all_cc_vs_c.de$NCC14ss_CC_vs_C_padj),]
    return(ncc14ss_up)
  })
  ncc14ss_dn <- reactive({
    ncc14ss_dn <- all_cc_vs_c.de[all_cc_vs_c.de$NCC14ss_CC_vs_C_log2FoldChange < 0 & all_cc_vs_c.de$NCC14ss_CC_vs_C_padj < 0.05 & !is.na(all_cc_vs_c.de$NCC14ss_CC_vs_C_padj),]
    return(ncc14ss_dn)
  })
  
  # Subset paralogs in each diff exp list -------------------------------------
  
  epi_up.de <- reactive({
    epi_up.de <- subset(table_paralog_list(), table_paralog_list()$EnsemblID %in% epi_up()$geneid)
  })
  epi_dn.de <- reactive({  
    epi_dn.de <- subset(table_paralog_list(), table_paralog_list()$EnsemblID %in% epi_dn()$geneid)
  })
  
  ncc5.6s_up.de <- reactive({ 
    ncc5.6s_up.de <- subset(table_paralog_list(), table_paralog_list()$EnsemblID %in% ncc5.6ss_up()$geneid)
  })
  ncc5.6s_dn.de <- reactive({ 
    ncc5.6s_dn.de <- subset(table_paralog_list(), table_paralog_list()$EnsemblID %in% ncc5.6ss_dn()$geneid)
  })
  
  ncc12s_up.de  <- reactive({ 
    ncc12s_up.de <- subset(table_paralog_list(), table_paralog_list()$EnsemblID %in% ncc12ss_up()$geneid)
  })
  ncc12s_dn.de <- reactive({ 
    ncc12s_dn.de <- subset(table_paralog_list(), table_paralog_list()$EnsemblID %in% ncc12ss_dn()$geneid)
  })
  
  ncc14s_up.de <- reactive({ 
    ncc14s_up.de <- subset(table_paralog_list(), table_paralog_list()$EnsemblID %in% ncc14ss_up()$geneid)
  })
  ncc14s_dn.de <- reactive({ 
    ncc14s_dn.de <- subset(table_paralog_list(), table_paralog_list()$EnsemblID %in% ncc14ss_dn()$geneid)
  })
  
  # Get names & numbersof Diff exp paralogs  
  
  diff_exp_names <- reactive({
    diff_exp_names <- list(toString(epi_names(epi_up.de())), toString(epi_names(epi_dn.de())), toString(epi_names(ncc5.6s_up.de())), toString(epi_names(ncc5.6s_dn.de())), toString(epi_names(ncc12s_up.de())), toString(epi_names(ncc12s_dn.de())), toString(epi_names(ncc14s_up.de())), toString(epi_names(ncc14s_dn.de())))
    #print("Names")
    #print(diff_exp_names)
    return(diff_exp_names)
  })
  
  diff_exp_nos <- reactive({
    diff_exp_nos <- list(nrow(epi_up.de()), nrow(epi_dn.de()), nrow(ncc5.6s_up.de()), nrow(ncc5.6s_dn.de()),nrow(ncc12s_up.de()), nrow(ncc12s_dn.de()),nrow(ncc14s_up.de()), nrow(ncc14s_dn.de()))
    #print("Numbers")
    #print(diff_exp_nos)
    return(diff_exp_nos)
  })
  
  
  output$epi_upreg <- renderText({
    n1 <- toString(diff_exp_nos()[[1]])
    t1 <- diff_exp_names()[[1]]
    paste(n1, "x upregulated at epiboly: ", t1 , collapse=" ")
  })
  
  output$epi_dnreg <- renderText({
    n2 <- toString(diff_exp_nos()[[2]])
    t2 <- diff_exp_names()[[2]]
    paste(n2, "x downregulated at epiboly: ",t2)
  })
  
  output$ncc5.6ss_upreg <- renderText({
    n3 <-toString(diff_exp_nos()[[3]])
    t3 <- diff_exp_names()[[3]]
    paste(n3, "x upregulated at 5-6ss:",t3 )
  })
  
  output$ncc5.6ss_dnreg <- renderText({
    n4 <- toString(diff_exp_nos()[[4]])
    t4 <- diff_exp_names()[[4]]
    paste(n4, "x downregulated at 5-6ss:", t4)
  })
  
  output$ncc12ss_upreg <- renderText({
    n5 <- toString(diff_exp_nos()[[5]])
    t5 <- diff_exp_names()[[5]]
    paste(n5, "x upregulated at 12-14ss:", t5)
  })
  
  output$ncc12ss_dnreg <- renderText({
    n6 <- toString(diff_exp_nos()[[6]])
    t6 <- diff_exp_names()[[6]]
    paste(n6, "x downregulated at 12-14ss:", t6)
  })
  
  output$ncc14ss_upreg <- renderText({
    n7 <- toString(diff_exp_nos()[[7]])
    t7 <- diff_exp_names()[[7]]
    paste(n7, "x upregulated at 14-16ss:", t7 )
  })
  
  output$ncc14ss_dnreg <- renderText({
    n8 <- toString(diff_exp_nos()[[8]])
    t8 <- diff_exp_names()[[8]]
    paste(n8, "x downregulated at 14-16ss:",t8 )
  })
  
  diffExp_table <-  reactive({
    diffExp_table <- rbind(
      c("Epiboly", toString(diff_exp_nos()[[1]]), toString(diff_exp_nos()[[2]]), toString(diff_exp_names()[[1]]), toString(diff_exp_names()[[2]])),
      c("5-6ss",   toString(diff_exp_nos()[[3]]), toString(diff_exp_nos()[[4]]), toString(diff_exp_names()[[3]]), toString(diff_exp_names()[[4]])),
      c("12-14ss", toString(diff_exp_nos()[[5]]), toString(diff_exp_nos()[[6]]), toString(diff_exp_names()[[5]]), toString(diff_exp_names()[[6]])),
      c("14-16ss", toString(diff_exp_nos()[[7]]), toString(diff_exp_nos()[[8]]), toString(diff_exp_names()[[7]]), toString(diff_exp_names()[[8]])))
    colnames(diffExp_table) <- c("Stage","# Upregulated"," #Downregulated", "Names Upregulated", "Names Downregulated")
    return(diffExp_table)
  })
  
  output$diffExp_table  <- renderTable({diffExp_table()})
  
  log_tf_fpkms <-  reactive({
    tf_fpkms <- fpkm[rownames(fpkm) %in% tf_id(),]
    tf_fpkms <- tf_fpkms[,c(1:4,7:10,13:20,23)]
    log_tf_fpkms <- log(tf_fpkms[,c(1:ncol(tf_fpkms)-1)]+ 0.1)
    log_tf_fpkms$genes <- annot[match(rownames(log_tf_fpkms), annot$Ensembl.Gene.ID),]$Associated.Gene.Name
    return(log_tf_fpkms)
  })
  
  # Heatmap all stages --------------------------
  
   output$heatmap1 <- renderD3heatmap({
     d3heatmap(as.matrix(log_tf_fpkms()[,c(1:ncol(log_tf_fpkms())-1)]),
               show_grid = FALSE,
               anim_duration = 0,
               colors = plasma(100),
               xaxis_font_size = "8px",
               yaxis_font_size = "8px",
               xaxis_height = 150,
               yaxis_width = 150,
               labRow=log_tf_fpkms()[,ncol(log_tf_fpkms())]
     )
   })
   
  # Table with original PWMs to cluster 
   
   ori_table <-  reactive({
     ori_table <- get(paste0(input$tf_select,".ori"))
     ori_table <- as.data.frame(ori_table)
     #print(ori_table)
   })
   
   output$ori_table <- DT::renderDataTable({
     DT::datatable(ori_table(),escape=FALSE)
   })      
   
   tfbs_enh <- reactive({
     #print(tfbs)
     sel_colname <- paste(input$tf_select, unique(ori_table()$Clustered.Motif.ID), sep="_")
     
     tfbs_cols <- colnames(tfbs) %in% c(sel_colname, "geneid")
     print ("tfbs_cols (logical):")
     #print(tfbs_cols)
     
     enh_with_tf <- tfbs[,tfbs_cols]
     enh_with_tf <- enh_with_tf[enh_with_tf > 0,]
     print(enh_with_tf)
     return(enh_with_tf$geneid)
   })
  
 # output$tfbs_enh  <- DT::renderDataTable({
#    if(ncol(tfbs_enh()) < 1)
#      return(NULL)
#    DT::datatable(tfbs_enh())
#  })      
  

  # Cluster parsing part ------------------------------------------------
  
  selected_cluster_list <- reactive({
    selected_cluster_list <- input$cluster_select
    selected_cluster_list <- selected_cluster_list[selected_cluster_list != "none"]
    print(selected_cluster_list)
    return(selected_cluster_list) 
  })
  
  output$cluster_text <- renderText({
    if (length(selected_cluster_list()) > 0 ){
      paste0("You have selected the following clusters: ",toString(clust_nbs[clust_nbs$real_name %in% selected_cluster_list(),]$name))
    }else{
      paste0("You haven't selected any clusters ...")
    }
  })

  nb_clusters_sel <- reactive({
    if (length(selected_cluster_list()) == 0){
      nb_clusters_sel <- rep(0, 6)
      }else if(length(selected_cluster_list()) > 0){
        nb_clusters_sel <- as.list(c(nrow(subset(clust_nbs[1:13,], clust_nbs[1:13,]$real_name %in% selected_cluster_list())),
                             nrow(subset(clust_nbs[14:29,], clust_nbs[14:29,]$real_name %in% selected_cluster_list())),
                             nrow(subset(clust_nbs[30:40,], clust_nbs[30:40,]$real_name %in% selected_cluster_list())),
                             nrow(subset(clust_nbs[41:60,], clust_nbs[41:60,]$real_name %in% selected_cluster_list())),
                             nrow(subset(clust_nbs[61:80,], clust_nbs[61:80,]$real_name %in% selected_cluster_list())),
                             nrow(subset(clust_nbs[81:100,], clust_nbs[81:100,]$real_name %in% selected_cluster_list()))))
        }
    print("Number clusters selected=")
    print(nb_clusters_sel)
  })
  
  cluster_names <- reactive({
    if (length(selected_cluster_list()) == 0){
      cluster_names <- rep("none", 6)
      }else if(length(selected_cluster_list()) > 0){
        cluster_names <- as.list(c(toString(subset(clust_nbs[1:13,], clust_nbs[1:13,]$real_name %in% selected_cluster_list())$name),
                           toString(subset(clust_nbs[14:29,], clust_nbs[14:29,]$real_name %in% selected_cluster_list())$name),
                           toString(subset(clust_nbs[30:40,], clust_nbs[30:40,]$real_name %in% selected_cluster_list())$name),
                           toString(subset(clust_nbs[41:60,], clust_nbs[41:60,]$real_name %in% selected_cluster_list())$name),
                           toString(subset(clust_nbs[61:80,], clust_nbs[61:80,]$real_name %in% selected_cluster_list())$name),
                           toString(subset(clust_nbs[81:100,], clust_nbs[81:100,]$real_name %in% selected_cluster_list())$name)))
      }
    print("Names of selected clusters =")
    print(cluster_names)
  })
  
    nb_positions <- reactive({
      if (length(selected_cluster_list()) == 0){
        nb_positions <- rep(0, 6)
        }else if(length(selected_cluster_list()) > 0){
          nb_positions <- as.list(c(sum(subset(clust_nbs[1:13,], clust_nbs[1:13,]$real_name %in% selected_cluster_list())$nb_peaks_in_cluster),
                            sum(subset(clust_nbs[14:29,], clust_nbs[14:29,]$real_name %in% selected_cluster_list())$nb_peaks_in_cluster),
                            sum(subset(clust_nbs[30:40,], clust_nbs[30:40,]$real_name %in% selected_cluster_list())$nb_peaks_in_cluster),
                            sum(subset(clust_nbs[41:60,], clust_nbs[41:60,]$real_name %in% selected_cluster_list())$nb_peaks_in_cluster),
                            sum(subset(clust_nbs[61:80,], clust_nbs[61:80,]$real_name %in% selected_cluster_list())$nb_peaks_in_cluster),
                            sum(subset(clust_nbs[81:100,], clust_nbs[81:100,]$real_name %in% selected_cluster_list())$nb_peaks_in_cluster)))
        }
      print("Number enhancers selected=")
      print(nb_positions)
    })
    
    positions_subset_cluster <- reactive({
      cluster_table2 <- cluster_table
      rownames(cluster_table2) <- cluster_table2$row.names
      cluster_table2$row.names <- NULL
      subset_cluster_table <- cluster_table2[,colnames(cluster_table2) %in% selected_cluster_list()] # subset only the columns from table that were selected
      if (length(selected_cluster_list())==1){
        print("Only 1 cluster selected")
        positions_subset_cluster <- cluster_table$row.names[subset_cluster_table > 0]
        print(dim(positions_subset_cluster))
      }else if(length(selected_cluster_list()) > 1){
        print("Multiple clusters selected")
        subset_cluster_table$sum <- rowSums(subset_cluster_table)
        print(head(subset_cluster_table))
        positions_subset_cluster <- rownames(subset_cluster_table[subset_cluster_table$sum == (ncol(subset_cluster_table)-1),][,1:(ncol(subset_cluster_table)-1)])
        print(dim(positions_subset_cluster))
      }
      return(positions_subset_cluster)
      })
    
    last_line <- reactive({
      if (length(selected_cluster_list())==0){
        last_line <- c("Total overlap in selection", 0, "NA", 0)
      }else{
      last_line <- c("Total overlap in selection", as.character(length(selected_cluster_list())), "NA", as.character(length (positions_subset_cluster())))
      print(paste0("last_line ", last_line))
      }
      return(last_line)
    })

  table_tss_cluster <- reactive({
    Cluster_group <- as.list(c("NCC Peak Clustering", "Sox10 Peak Clustering", "H3K27Ac Peak Clustering", "Epiboly Nucleosome Clustering", "HB Nucleosome Clustering", "July Nucleosome Clustering"))
    df <- cbind(Cluster_group, nb_clusters_sel(), cluster_names(), nb_positions())
    df <- as.data.frame(df)
    colnames(df) <- c("Group", "Number Selected Cluster(s)", "Name(s) of selected clusters", "Number cis-reg. elements in cluster")
    df[7,] <- last_line()
    return(df)
  })
   
  table_positions_subset_clusters <- reactive({
    # If missing input, return to avoid error later in function
    if(length(selected_cluster_list()) == 0)
      return(NULL)
    table_positions_subset_clusters <- Enh_peaks[Enh_peaks$name %in% positions_subset_cluster(),]
    annot_table_positions_subset_clusters <- merge(table_positions_subset_clusters, annot,by.x="EnsemblGeneID",by.y="Ensembl.Gene.ID")
    annot_table_positions_subset_clusters <- annot_table_positions_subset_clusters[,c(2,3,4,1,6)]
    return(annot_table_positions_subset_clusters)
  })

  output$table_clusters_selected <- renderTable({table_tss_cluster()})
  
  output$table_cluster_positions <- DT::renderDataTable({
    DT::datatable(table_positions_subset_clusters(),escape=FALSE)
  })      
  
  # Boxplot of expression levels of genes associated to selected clusters --------------------------------
  output$boxplot <- renderPlotly({
    df <-  subset(melted_fpkms, melted_fpkms$geneid %in% unique(table_positions_subset_clusters()[,4]))
    df$col_cond <- ifelse(grepl("_CC", df$variable, ignore.case = T), "CC", 
                          ifelse(grepl("_C[[:digit:]]", df$variable, ignore.case = T), "wt",
                                 ifelse(grepl("Sox", df$variable,ignore.case = T), "Sox","Other")))
    p <- plot_ly(data=df, 
                 x =df$value, 
                 y=df$variable,
                 type = "box", 
                 marker=list( size=2 , opacity=0.7), 
                 color = ~col_cond,
                 colors=c("firebrick","grey","purple","darkgreen"))  %>% layout(margin = list(l=110))
    
  })
  
  
  output$event <- renderPrint({
    d <- event_data("plotly_hover")
  })  
  
})
