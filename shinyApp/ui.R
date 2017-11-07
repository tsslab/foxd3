#ui.R

# load useful packages/libraries ----------------------

library(shiny)

if(!require("shinythemes")){
  install.packages("shinythemes")
  library(shinythemes)
}

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

if(!require("knitr")){
  install.packages("knitr")
  library(knitr)
}

# load necessary data ----------------------
load("RData/annot.RData")

# User Interface  ----------------------

shinyUI(fluidPage(theme = shinytheme("cosmo"),
  titlePanel("FoxD3"),
  sidebarLayout(
    # Sidebar -----------------------------------------
    sidebarPanel(position = "left",
                 radioButtons("chto_delat", "Select one:",
                              choices=c("Search for genes of interest"=1,
                                "Search putative cis-regulatory clusters"=2,
                                "Parse TFBS"=3))
                 ),
    # Main  -----------------------------------------
                 mainPanel(
                   tabsetPanel(
                     tabPanel("Data", value=1, 
                              # RNA-seq  -----------------------------------------
                              
                              conditionalPanel(condition="input.chto_delat==1", 
                                               verbatimTextOutput("RNA-seq"),
                                               helpText("Select a gene that you are  
               interested in"),
                                               selectInput('Gene_Name',
                                                              label=NULL,
                                                              choices=annot$Associated.Gene.Name,
                                                              selected = "foxd5",
                                                              multiple = FALSE),
                                               # uiOutput("Gene_Name_UI"),
                                               
                                               h3("Gene Information"),
                                               helpText("Some information about the gene you have selected."),   
                                               textOutput("SelectedGene"),
                                               textOutput("Ensembl_Geneid"),
                                               textOutput("Full_Gene_Name"),
                                               br(),
                                               h3("Expression levels"),
                                               helpText("Expression levels measured in FPKMs (Fragments Per Kilobase of 
                                                        transcript per Million mapped reads)"),
                                               plotlyOutput("plot", width = "80%", height = "400px"),
                                               verbatimTextOutput("event"),
                                               br(),
                                               h3("Differential Expression"),
                                               helpText("Differential expression performed using R package DESeq2. Genes were 
                                                        deemed significantly differentially expressed if adjusted p-value < 0.05"),
                                               h4("FoxD3+ cells vs FoxD3 - cells"),
                                               tableOutput('table_negs.f'),
                                               br(),
                                               h4("FoxD3+ cells expressing mutant version of FoxD3 vs FoxD3 + wt cells"),
                                               tableOutput('table_mut.f'),
                                               #plotlyOutput("HM_CC_vs_C", width = "850px", height = "800px")),
                                               br(),
                                               h3("Associated proximal ATAC-seq peaks"),
                                               helpText("ATAC-seq peaks were deemed associated to a promoter using annotatePeaks.pl from Homer (v.4.8)"),
                                               tableOutput('table_promoter_peaks'),
                                               br(),
                                               h3("Associated distal cis-regulatory elements"),
                                               helpText("Enhancer associations were determined using closestBed function from 
                                                        bedtools (v.2.15.0) to genes expressed above 1 FPKM at 5-6ss in Citrine samples "),
                                               br(),
                                               tableOutput('table_enhancer_peaks'), 
                                               br(),
                                               h3("TFBS in associated distal cis-regulatory elements"),
                                               #tableOutput('table_tfbs_gene')
                                               d3heatmapOutput("HM",  width = "100%", height = "1200px")     
                                               ),
                              # Cluster parsing -----------------------------------------
                              conditionalPanel(condition="input.chto_delat==2", 
                                               verbatimTextOutput("Clusters"),
                                               br(),
                                               h3("Cluster Selection"),
                                               br(),
                                               helpText("Select a pre-parsed cluster to examine:"),
                                               
                                               selectizeInput(
                                                 'cluster_select', NULL, 
                                                 choices = c("none" = "none",
                                                             "NCC cluster1"= "july_clust_cl1_6777",
                                                             "NCC cluster2" = "july_clust_cl2_10751_control",
                                                             "NCC cluster3" = "july_clust_cl3_17390",
                                                             "NCC cluster4" = "july_clust_cl4_4493",
                                                             "NCC cluster5" = "july_clust_cl5_8182_inverted",
                                                             "NCC clusters567"="july_clust_cl5_6_7_inverted",
                                                             "NCC cluster67" = "july_clust_cl6_7_inverted",
                                                             "NCC cluster8" = "july_clust_cl8_8245",
                                                             "NCC clusters148" = "july_clust_cl1_4_8",
                                                             "NCC cluster3_13468" = "july_clust_subcl3_Cl1_3_4_6_8",
                                                             "NCC cluster3_1368" = "july_clust_subcl3_Cl1_3_6_8",
                                                             "NCC cluster3_257" = "july_clust_subcl3_Cl2_5_7_control",
                                                             "NCC cluster3_4" = "july_clust_subcl3_Cl4",
                                                             "Sox10 clusters125" = "sox10_clust_Cl1_2_5",
                                                             "Sox10 cluster1"="sox10_clust_Cl1_8572",
                                                             "Sox10 cluster2" = "sox10_clust_Cl2_6311",
                                                             "Sox10 cluster3" = "sox10_clust_Cl3_13541_ctrl", 
                                                             "Sox10 cluster4" = "sox10_clust_Cl4_3341_inverted",
                                                             "Sox10 cluster5" = "sox10_clust_Cl5_11184",
                                                             "Sox10 cluster6" = "sox10_clust_Cl6_14524",
                                                             "Sox10 cluster7" = "sox10_clust_Cl7_3156",
                                                             "Sox10 cluster6_12345" = "cluster6_12345",
                                                             "Sox10 cluster6_1235" = "cluster6_1235",
                                                             "Sox10 cluster6_235" = "sox10_clust_subcl6_Cl2_3_5",
                                                             "Sox10 cluster6_1" = "sox10_clust_subcl6_Cl1_4045",
                                                             "Sox10 cluster6_2" = "sox10_clust_subcl6_Cl2_1080",
                                                             "Sox10 cluster6_4" = "sox10_clust_subcl6_Cl4_2622",
                                                             "Sox10 cluster6_5" = "sox10_clust_subcl6_Cl5_1617",
                                                             "Sox10 clusters6_678"="sox10_clust_subcl6_Cl6_7_8_control",
                                                             "H3K27Ac clusters12" = "h3k27ac_clust_1_2_pooled_control",
                                                             "H3K27Ac clusters36" = "h3k27ac_clust_3_6",
                                                             "H3K27Ac cluster4" = "h3k27ac_clust_4_6726_inverted",
                                                             "H3K27Ac cluster5" = "h3k27ac_clust_5_10538",
                                                             "H3K27Ac cluster8" = "h3k27ac_clust_8_2667_inverted",
                                                             "H3K27Ac cluster9" = "h3k27ac_clust_9_2287",
                                                             "H3K27Ac cluster10" = "h3k27ac_clust_10_2071",
                                                             "H3K27Ac cluster5_1" = "h3k27ac_subcl5_cl1_2395",
                                                             "H3K27Ac cluster5_2" = "h3k27ac_subcl5_cl2_3151",
                                                             "H3K27Ac cluster5_3" = "h3k27ac_subcl5_cl3_846",
                                                             "H3K27Ac cluster5_4" = "h3k27ac_subcl5_cl4_4146",
                                                             "Epiboly nucleosomes cluster1" = "cluster1_EpiNucK20",
                                                             "Epiboly nucleosomes cluster2" = "cluster2_EpiNucK20",
                                                             "Epiboly nucleosomes cluster3" = "cluster3_EpiNucK20", 
                                                             "Epiboly nucleosomes cluster4" = "cluster4_EpiNucK20",
                                                             "Epiboly nucleosomes cluster5" = "cluster5_EpiNucK20",
                                                             "Epiboly nucleosomes cluster6" = "cluster6_EpiNucK20",
                                                             "Epiboly nucleosomes cluster7" = "cluster7_EpiNucK20",
                                                             "Epiboly nucleosomes cluster8" = "cluster8_EpiNucK20",
                                                             "Epiboly nucleosomes cluster9" = "cluster9_EpiNucK20",
                                                             "Epiboly nucleosomes cluster10" = "cluster10_EpiNucK20",
                                                             "Epiboly nucleosomes cluster11" = "cluster11_EpiNucK20",
                                                             "Epiboly nucleosomes cluster12" = "cluster12_EpiNucK20",
                                                             "Epiboly nucleosomes cluster13" = "cluster13_EpiNucK20",
                                                             "Epiboly nucleosomes cluster14" = "cluster14_EpiNucK20",
                                                             "Epiboly nucleosomes cluster15" = "cluster15_EpiNucK20",
                                                             "Epiboly nucleosomes cluster16" = "cluster16_EpiNucK20",
                                                             "Epiboly nucleosomes cluster17" = "cluster17_EpiNucK20",
                                                             "Epiboly nucleosomes cluster18" = "cluster18_EpiNucK20",
                                                             "Epiboly nucleosomes cluster19" = "cluster19_EpiNucK20",
                                                             "Epiboly nucleosomes cluster20" = "cluster20_EpiNucK20",
                                                             "HB nucleosomes cluster1" = "cluster_1_HBNucK20",
                                                             "HB nucleosomes cluster2" = "cluster_2_HBNucK20",
                                                             "HB nucleosomes cluster3" = "cluster_3_HBNucK20", 
                                                             "HB nucleosomes cluster4" = "cluster_4_HBNucK20",
                                                             "HB nucleosomes cluster5" = "cluster_5_HBNucK20",
                                                             "HB nucleosomes cluster6" = "cluster_6_HBNucK20",
                                                             "HB nucleosomes cluster7" = "cluster_7_HBNucK20",
                                                             "HB nucleosomes cluster8" = "cluster_8_HBNucK20",
                                                             "HB nucleosomes cluster9" = "cluster_9_HBNucK20",
                                                             "HB nucleosomes cluster10" = "cluster_10_HBNucK20",
                                                             "HB nucleosomes cluster11" = "cluster_11_HBNucK20",
                                                             "HB nucleosomes cluster12" = "cluster_12_HBNucK20",
                                                             "HB nucleosomes cluster13" = "cluster_13_HBNucK20",
                                                             "HB nucleosomes cluster14" = "cluster_14_HBNucK20",
                                                             "HB nucleosomes cluster15" = "cluster_15_HBNucK20",
                                                             "HB nucleosomes cluster16" = "cluster_16_HBNucK20",
                                                             "HB nucleosomes cluster17" = "cluster_17_HBNucK20",
                                                             "HB nucleosomes cluster18" = "cluster_18_HBNucK20",
                                                             "HB nucleosomes cluster19" = "cluster_19_HBNucK20",
                                                             "HB nucleosomes cluster20" = "cluster_20_HBNucK20",
                                                             "NCC nucleosomes cluster1" = "cluster_1_JulyNucK20",
                                                             "NCC nucleosomes cluster2" = "cluster_2_JulyNucK20",
                                                             "NCC nucleosomes cluster3" = "cluster_3_JulyNucK20", 
                                                             "NCC nucleosomes cluster4" = "cluster_4_JulyNucK20",
                                                             "NCC nucleosomes cluster5" = "cluster_5_JulyNucK20",
                                                             "NCC nucleosomes cluster6" = "cluster_6_JulyNucK20",
                                                             "NCC nucleosomes cluster7" = "cluster_7_JulyNucK20",
                                                             "NCC nucleosomes cluster8" = "cluster_8_JulyNucK20",
                                                             "NCC nucleosomes cluster9" = "cluster_9_JulyNucK20",
                                                             "NCC nucleosomes cluster10" = "cluster_10_JulyNucK20",
                                                             "NCC nucleosomes cluster11" = "cluster_11_JulyNucK20",
                                                             "NCC nucleosomes cluster12" = "cluster_12_JulyNucK20",
                                                             "NCC nucleosomes cluster13" = "cluster_13_JulyNucK20",
                                                             "NCC nucleosomes cluster14" = "cluster_14_JulyNucK20",
                                                             "NCC nucleosomes cluster15" = "cluster_15_JulyNucK20",
                                                             "NCC nucleosomes cluster16" = "cluster_16_JulyNucK20",
                                                             "NCC nucleosomes cluster17" = "cluster_17_JulyNucK20",
                                                             "NCC nucleosomes cluster18" = "cluster_18_JulyNucK20",
                                                             "NCC nucleosomes cluster19" = "cluster_19_JulyNucK20",
                                                             "NCC nucleosomes cluster20" = "cluster_20_JulyNucK20"),
                                                 multiple = TRUE),
                                               br(),
                                               #actionButton("hideshow", "Hide/show plot of last selected cluster"),
                                               textOutput("cluster_text"),
                                               br(),
                                               #htmlOutput('cluster_image'),
                                               br(),
                                               h3("Numbers of putative cis-regulatory elements"),
                                               br(),
                                               tableOutput('table_clusters_selected'),
                                               br(),
                                               h3("Gene associations of clustered putative cis-regulatory elements"),
                                               br(),
                                               DT::dataTableOutput('table_cluster_positions', width="80%"), 
                                               h3("Expression levels of gene associated of clustered putative cis-regulatory elements "),
                                               br(),
                                               plotlyOutput("boxplot", width = "80%", height = "400px"),
                                               verbatimTextOutput("event2"),
                                               br(),
                                               h3("TFBS enrichment in selected clusters")
                              ),
                     
                     # TFBS parsing -----------------------------------------
                     
                              conditionalPanel(condition="input.chto_delat==3", 
                                               verbatimTextOutput("TFBS"),
                                               helpText("Select a transcription factor family to explore binding site"),
                                               selectInput("tf_select", "Choose TF family :", 
                                                           choices = c("alx" = "alx", 
                                                                       "cdx" = "cdx",
                                                                       "cux" = "cux",
                                                                       "cxxc" = "cxxc",
                                                                       "ctcf" = "ctcf",
                                                                       "dlx" = "dlx",
                                                                       "ebf" = "ebf",
                                                                       "egr" = "egr",
                                                                       "emx" = "emx",
                                                                       "ets" = "ets",
                                                                       "fli" = "fli",
                                                                       "fox" = "fox",
                                                                       "gata" = "gata",
                                                                       "gli" = "gli", 
                                                                       "hand" = "hand", 
                                                                       "hif" = "hif", 
                                                                       "irx" = "irx", 
                                                                       "jun" = "jun", 
                                                                       "klf" = "klf", 
                                                                       "lmx" = "lmx", 
                                                                       "mef" = "mef", 
                                                                       "msx" = "msx", 
                                                                       "nfat" = "nfat", 
                                                                       "nkx" = "nkx", 
                                                                       "nr2f" = "nr2f", 
                                                                       "otx" = "otx", 
                                                                       "pax" = "pax", 
                                                                       "POU" = "pou",
                                                                       "prdm" = "prdm", 
                                                                       "six" = "six", 
                                                                       "snail" = "snail", 
                                                                       "sox" = "sox", 
                                                                       "stat" = "stat", 
                                                                       "tal" = "tal", 
                                                                       "tbx" = "t", 
                                                                       "tcf" = "tcf", 
                                                                       "tead" = "tead", 
                                                                       "tfap2" = "tfap2", 
                                                                       "tfe" = "tfe",
                                                                       "twist" = "twist", 
                                                                       "ventx" = "ventx", 
                                                                       "zic" = "zic",
                                                                       "znf" = "znf")),
                                h3("Gene Family information"),
                                helpText("Some information regarding the genes whose protein products may be binding to this motif"),   
                                textOutput("text1"),
                                textOutput("text2"),
                                br(),
                                tableOutput('table_paralog_list'),
                                br(),
                                h3("Differentially expressed?"),
                                tableOutput('diffExp_table'),
                                br(),
                                h3("Expression levels"),
                                d3heatmapOutput("heatmap1",  width = "80%", height = "800px"),
                                h3("Binding Sites Parsing"),
                                helpText("Position Weighted Matrices were downloaded from CIS-BP (Weirauch et al. 2014). Binding sites were clustered using gimme suite’s cluster option"),
                                h5("Original PWMs"),
                                DT::dataTableOutput('ori_table', width="70%")#,
                                #h3("Number of binding sites found"),
                                #DT::dataTableOutput('tfbs_enh', width="80%")
                                )
                     
                      ),
                     # About us -----------------------------------------
                     tabPanel("About us", value=2,
                              includeMarkdown(knitr::knit("FoxD3.Rmd"))),
                     id = "tabselected"))
                   )
                 )
    )

