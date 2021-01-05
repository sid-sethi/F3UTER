library(shiny)
library(shinyjs)
library(shinythemes)
library(shinycssloaders)
library(tidyverse)
library(ggpubr)
library(DBI)
library(pool)
library(glue)
library(grid)
library(ggpubr)
library(GenomicRanges)
library(GenomicFeatures)
library(waiter)

options(spinner.color="#808080", spinner.type = 5)

waiting_screen <- tagList(
  h2("Loading  F3UTER"), br(),
  spin_balance()
  #spin_flower()
  #spin_dots()
) 


########### loading datasets ########################

#feuterdb <- DBI::dbConnect(RSQLite::SQLite(), "data/f3uterdb.sqlite")

feuterdb <- dbPool(
  drv = RSQLite::SQLite(),
  dbname = "data/f3uterdb.sqlite"
)

onStop(function() {
  poolClose(feuterdb)
})

# Load TxDb object
txdb_v92 <- loadDb(file = "data/gtf_TxDb_v92.Rsqlite")



############## functions ############
source("app_functions.R")
####### functions END ###############


####### tissue name formatting #########
tissue_GTEx_choices <- c("Adipose - subcutaneous" =	"adipose_subcutaneous",
                         "Adipose - visceral" =	"adipose_visceral_omentum",
                         "Adrenal gland" =	"adrenal_gland",
                         "Aorta" =	"artery_aorta",
                         "Artery - coronary" =	"artery_coronary",
                         "Artery - tibial" =	"artery_tibial",
                         "Amygdala" =	"brain_amygdala",
                         "Anterior cingulate cortex" =	"brain_anterior_cingulate_cortex_ba24",
                         "Caudate" =	"brain_caudate_basal_ganglia",
                         "Cerebellar hemisphere" =	"brain_cerebellar_hemisphere",
                         "Frontal Cortex" =	"brain_frontal_cortex_ba9",
                         "Hippocampus" =	"brain_hippocampus",
                         "Hypothalamus"	= "brain_hypothalamus",
                         "Nucleus accumbens" =	"brain_nucleus_accumbens_basal_ganglia",
                         "Putamen" =	"brain_putamen_basal_ganglia",
                         "Spinal cord" =	"brain_spinal_cord_cervical_c_1",
                         "Substantia nigra" =	"brain_substantia_nigra",
                         "Sigmoid" =	"colon_sigmoid",
                         "Transverse" =	"colon_transverse",
                         "Gastroesophageal junction" =	"esophagus_gastroesophageal_junction",
                         "Mucosa" =	"esophagus_mucosa",
                         "Muscularis" =	"esophagus_muscularis",
                         "Atrial appendage" =	"heart_atrial_appendage",
                         "Left ventricle" =	"heart_left_ventricle",
                         "Kidney" =	"kidney_cortex",
                         "Liver" =	"liver",
                         "Lung" =	"lung",
                         "Minor salivary gland" =	"minor_salivary_gland",
                         "Skeletal muscle" =	"muscle_skeletal",
                         "Nerve - tibial" =	"nerve_tibial",
                         "Pancreas" =	"pancreas",
                         "Pituitary" =	"pituitary",
                         "Skin (suprapubic)" =	"skin_not_sun_exposed_suprapubic",
                         "Skin (lower leg)"	= "skin_sun_exposed_lower_leg",
                         "Small Intestine" =	"small_intestine_terminal_ileum",
                         "Spleen" =	"spleen",
                         "Stomach" =	"stomach",
                         "Thyroid" =	"thyroid",
                         "Whole blood" =	"whole_blood",
                         "All" = "All")

tissue_GTEx_choices_alphabetical <- tissue_GTEx_choices[names(tissue_GTEx_choices) %>% order()]

#########################################################



ui <- fluidPage( 
  
  use_waiter(spinners = 5),
  waiter_show_on_load(html = waiting_screen, color = "black"),
  use_waitress(percent_color = "#006600", color = "#E3EBEE"),
  
  
  theme = shinytheme("flatly"),
                 
  tags$head(tags$style(".shiny-notification {position: fixed; bottom: 3% ; right: 2%; background-color:#E0E0E0; color:black}")),
  
  
  navbarPage(title = "Finding 3' Un-translated Expressed Regions", windowTitle = "F3UTER",
        
            
             # main tab 1 
             tabPanel(title = "Home",
                      fluidRow(column(width = 8, offset = 2, h1(strong("F3UTER"), align="center"), hr(), h1(strong("Finding 3' Un-translated Expressed Regions"), align="center"),
                               wellPanel(
                                 h5("This platform can be used to identify unannotated 3'UTRs associated with a gene predicted by F3UTER."), br(),
                                 h5("F3UTER is a machine learning-based framework which leverages both genomic and tissue-specific transcriptomic features to predict previously unannotated 3’UTRs in the human genome. F3UTER was applied to transcriptomic data covering 39 human tissues studied within GTEx, enabling the identification of tissue-specific unannotated 3’UTRs for 1,513 genes."), br(),
                                 h5("This resource takes a 'Gene identifier' as input and provides extensive information about unannotated 3’UTR predictions associated with that gene. For an unannotated 3'UTR of interest, the user can visualise it in a genome-browser based view, and inspect the omic features driving F3UTER’s prediction."), br(),
                                 div(img(src="shiny_f3uter_homePage.png", width="85%"), style="text-align: center;"),
                                 br(),
                                 h5(strong("For more deails on F3UTER, please see: "), "Leveraging omic features with F3UTER enables identification of unannotated 3’UTRs for synaptic genes"),
                                 h5(strong("For more details on unannotated ERs, please see: "), a(href="https://advances.sciencemag.org/content/6/24/eaay8299.full", "Incomplete annotation of disease-associated genes is limiting our understanding of Mendelian and complex neurogenetic disorders", target="_blank"))
                              
                              ), # well panel
                              hr(),
                              h4(strong("Developed in partnership with", a(href="https://astx.com", img(src="astexlogo.gif", width="15%"), target="_blank"), " under an Astex Sustaining Innovation Postdoctoral Fellowship"), align="center"), hr()
                      ))
             ),
             
             
             # main tab 2
             tabPanel(
               
               title = "F3UTER",
  
                  titlePanel("F3UTER"),
                  useShinyjs(),
  
                  sidebarLayout(
                    sidebarPanel(
                      width = 3, 
                      textInput(inputId = "gene_input", label = "Gene identifier", placeholder = "Enter gene symbol or ENSG ID"),
                      selectInput(inputId = "tissue_input", label = "Tissue", selected = "All", choices = tissue_GTEx_choices_alphabetical),
                      actionButton("submit_gene", label = "Submit gene", icon("arrow-alt-circle-right"), 
                                   style="color: #fff; background-color: #337ab7; border-color: #2e6da4"),
                      verbatimTextOutput("gene_output"),
                      hr(),
                      shinyjs::hidden(helpText(id = "select_row_text", "Select a row from the Gene-table and submit ER to visualise it")),
                      verbatimTextOutput("er_output"),
                      shinyjs::hidden(actionButton("submit_er", label = "Submit ER", icon("arrow-alt-circle-right"), 
                                                     style="color: #fff; background-color: #337ab7; border-color: #2e6da4")),
                      hr(),
                      shinyjs::hidden(sliderInput(inputId = "zoom_factor_input", label = "'Zoom-in' value for genomic view:",min = 1, max = 40, value = 30, step = 1 )),
                      hr(),
                      shinyjs::hidden(downloadButton("downloadGenomeView", label = "GenomeView plot",
                                                     style="color: #fff; background-color: #000000; border-color: #000000")),
                      shinyjs::hidden(downloadButton("downloadFeature", label = "Feature plot",
                                                     style="color: #fff; background-color: #000000; border-color: #000000"))
                    ),
    
                    
                    mainPanel(
                      tabsetPanel(id = "tabs",
                                  type = "tabs", selected = "Gene-table",
                                  tabPanel("Gene-table", DT::dataTableOutput("genetable") %>% withSpinner()),
                                  #tabPanel("ER genomeView plot", plotOutput("genomeView", width = "90%", height = "500px") %>% withSpinner()),
                                  tabPanel("ER genomeView plot", plotOutput("genomeView", width = "90%", height = "500px")),
                                  #tabPanel("ER features plot", plotOutput("featureResult", width = "90%", height = "700px") %>% withSpinner()),
                                  tabPanel("ER features plot", plotOutput("featureResult", width = "90%", height = "700px") ),
                                  tabPanel("Download",
                                    fluidRow( column(width=7,
                                      br(),
                                      h4("Gene data download"),
                                      wellPanel(
                                        h5("Download Gene-table (.csv)"),
                                        uiOutput("download_gene_table")
                                      ),
                                      br(),
                                      h4("ER data download"),
                                      wellPanel(
                                        h5("Download machine learning features of ERs in Gene-table (.csv)"),
                                        uiOutput("download_er_feature_table"),
                                        br(),
                                        h5("Download data frame of ER feature plot (.csv)"),
                                        uiOutput("download_er_feature_plot_df"),
                                        br()
                                      )
                                    ))
                                           
                                  )
                                  
                      ) # tabsetPanel end
                    ) # mainPanel end
                  ) # sidebar layout
            ), # F3UTER tab panel
            
            
            
            # main tab 3 
            tabPanel(title = "Documentation",
                     fluidRow(column(width = 8, offset = 2, 
                      h3(strong("Using F3UTER to identify potentially unannotated 3'UTRs")), br(),
                      shiny::tags$ul(
                        shiny::tags$li(h4(strong("Accessing F3UTER"))),
                        p("F3UTER can be accessed by clicking on the tab named 'F3UTER' on the navigation bar at the top of the page."),
                        p(img(src="f3uter_accessing_tab.png", width="65%"), align = "center"), hr(),
                        
                        shiny::tags$li(h4(strong("Input parameters"))),
                        p("F3UTER requires a gene idenitifier and tissue as input to retrieve the ERs associated with a gene."),
                        p(img(src="f3uter_input.png", width="65%"), align = "center"),
                        p(strong("1. Gene identifier: "), "You may input either a gene symbol or Ensembl ID (based on Ensembl version 94). You can only query one gene at a time. Only a valid gene identifier is accepted, otherwise an error is shown."),
                        p(strong("2. Tissue: "), "This parameter selects the tissue of interest. Unannotated 3'UTR predictions are available in 39 GTEx tissues and can be accessed the drop down menu. You may select a single tissue of choice or query 'All' tissues. The default is to retrieve 3'UTR predictions from all tissues."),
                        p(strong("3. Submit gene: "), "Click this button to submit your input choices to the tool."), hr(),
                        
                        shiny::tags$li(h4(strong("Gene-table output"))),
                        p(img(src="f3uter_gene_table.png", width="100%"), align = "center"),
                        p(strong("1. Gene information: "), "This provides the gene symbol and the Ensembl ID of the submitted gene for confirmation."),
                        p(strong("2. Gene-table: "), "Gene-table output provides all the unannotated 3'UTR predictions associated with the submitted gene in the selected tissue. An empty table will be shown if no data is found for the submitted gene. The Gene-table can be downloaded in csv format from the 'Download' tab. The output columns of Gene-table are explained below:"),
                        shiny::tags$ul(
                          shiny::tags$li(strong("seqnames"), ": Chromosome expressed region (ER) is found on (1-22, X or Y)."),
                          shiny::tags$li(strong("start"), ": Start position of expressed region (hg38)."),
                          shiny::tags$li(strong("end"), ": End position of expressed region (hg38)."),
                          shiny::tags$li(strong("width"), ": Length of expressed region in base pairs."),
                          shiny::tags$li(strong("tissue"), ": The GTEx tissue in which the expressed region was detected."),
                          shiny::tags$li(strong("predicted.response"), ": Prediction outcome of expressed region from F3UTER. Possible values: '3'UTR' or 'non-3'UTR'."),
                          shiny::tags$li(strong("predicted.probability"), ": Prediction probability associated with the F3UTER classification result. Predictions with prob>=0.6 are classified as 3'UTRs, while the remaining are classified as non-3'UTRs."),
                          shiny::tags$li(strong("ER_tissue_specificity_category"), ": Tissue-specificity of the expressed region across the 39 GTEx tissues. Possible values: 'Absolute tissue-specific', 'Highly brain-specific', 'Shared', 'Ambiguous'. Please see the methods section of the F3UTER paper for more details."),
                          shiny::tags$li(strong("distance_from_geneEnd"), ": Distance (in bp) between the expressed region and the end of the associated gene."),
                          shiny::tags$li(strong("annotationType_split_read_annot"), ": This annotation shows whether the expressed region has a connecting junction-read (split read) or not. Possible values: 'partially annotated split read' - the expressed region has a split read connected to an annotated gene; 'uannotated split read' - the expressed region has a split read connected to an unannotated region; 'none' - the expressed region has no split read."),
                          shiny::tags$li(strong("uniq_genes_split_read_annot"), ": If applicable, the gene connected to the expressed region via junction read."),
                          shiny::tags$li(strong("associated_gene"), ": The gene associated with the expressed region. The expressed regions are associated to a protein-coding gene within a 10 kb window, either via a junction read, or to the nearest protein-coding gene when no junction read exists. See the methods section of the F3UTER paper for more details."),
                          shiny::tags$li(strong("id"), ": A unique identifier of the expressed region (ER) formed by combining tissue of ER, chromosome of ER, start position of ER, end position of ER and strand of ER.")
                        
                        ), hr(),
                        
                        
                        shiny::tags$li(h4(strong("Selecting an expressed region (ER) to visualise"))),
                        p("You can further inspect and visulaise interesting ERs from the Gene-table output."),
                        p(img(src="f3uter_selecting_er.png", width="100%"), align = "center"),
                        p(strong("1. Select ER: "), "You can select an ER by clicking on an ER row in the Gene-table. Only one ER can be selected at a time. The unique 'id' of the selected ER is shown in the side panel. For instance, in the above picture, we selected the ER in the first row, as it was predicted as a 3'UTR by F3UTER. Furthermore, this ER has a high prediction probability (0.89) and has a partially annotated split read, which adds more confidence to the prediction."),
                        p(strong("2. Submit ER: "), "Click this button to submit your selected ER for visualisation. This automatically switches the ouput tab to 'ER genomeView plot'."), hr(),
                      
                        shiny::tags$li(h4(strong("ER genome-view plot"))),
                        p("The ER genome-view plot provides a genome browser based schematic visualisation of the ER with respect to its associated gene. Please note that this schematic is specific to the submitted gene and does not show other genes in the genomic window of the plot. The ER genome-view plot is explained in the picture below."),
                        p(img(src="f3uter_genome_view.png", width="100%"), align = "center"),
                        hr(),
                        
                        shiny::tags$li(h4(strong("ER features plot"))),
                        p("The ER features plot provides a visualisation of the features associated with the selected ER, alongside the features of 'known 3'UTRs' and 'known non-3'UTRs'. This can help you to understand the variables driving the prediction outcome by F3UTER. The features are displayed in decreasing order of their importance in 3'UTR classification. The feature importance was calculated from the training data. For more information on the individual features, please refer to the methods section of the F3UTER paper. The underlying data structure of this plot, and the ER features can be downloaded in a csv format under the 'Download' tab. The ER features plot is explained in the picture below."), 
                        p(img(src="f3uter_feature_plot.png", width="100%"), align = "center"),
                        hr()
                          
                      )

                                           
                )) 
            ),
            
            
            
            
            # main tab 4 
            tabPanel(title = "Data",
                     fluidRow(column(width = 8, offset = 2, br(), h3(strong("Downloads")), hr(),
                                     h4("Download all 3' intergenic ER predictions and its associated meta-data (all_ER_predictions.csv - 39.2 MB)"),
                                     downloadButton(outputId = "download_er_data", label = "ER predictions", style="color: #fff; background-color: #000000; border-color: #000000"), br(), br(), br(),
                                     h4("Download feature matrix of all 3' intergenic ERs (ER_feature_matrix.csv - 44.7 MB)"),
                                     downloadButton(outputId = "download_er_feature_matrix", label = "ER feature matrix", style="color: #fff; background-color: #000000; border-color: #000000"), br(), br(), br(),
                                     h4("Download training feature matrix of F3UTER (training_feature_matrix.csv - 127 MB)"),
                                     downloadButton(outputId = "download_training_feature_matrix", label = "Training feature matrix", style="color: #fff; background-color: #000000; border-color: #000000"), br(), br(), br(), hr(),
                                     p("It may take 5-10 seconds to start the download")
                        
                     
                     ))
            ),
            
            
            
            # main tab 5 
            tabPanel(title = "Terms of Use",
                     fluidRow(column(width = 8, offset = 2, br(), h3(strong("F3UTER - Terms of Use")), hr(),
                                     h4(strong("User Submissions")),
                                     p("By submitting data to or using the F3UTER App, you agree to comply with and be bound by these Terms of Use."), br(),
                                     h4(strong("Use of Data by Astex")),
                                     p("Astex has no visibility of any data that you submit to the F3UTER App and can only view App usage metrics provided by shinyapps.io such as network usage, CPU usage, memory usage and number of hours of App usage."), br(),
                                     h4(strong("Use of Data by User")),
                                     p("A user has the rights granted by the licence below, along with the right to download the data file generated by the F3UTER App when used by the user."), br(),
                                     h4(strong("Trademarks")),
                                     p("The trademarks, trade names, logos and service marks (collectively the 'Trademarks') displayed on the F3UTER App are registered and unregistered Trademarks of Astex and others. Nothing contained on the website should be construed as granting, by implication, estoppel or otherwise, any license or right to use any Trademark displayed on the F3UTER App without the written permission of Astex or such third party that may own the Trademarks displayed on the F3UTER App."), br(),
                                     h4(strong("Modification of these Terms")),
                                     p("Astex may change these above Terms of Use from time to time and at any time, without notice to the user. Unless otherwise noted by Astex, changes to these Terms will become effective immediately after they are posted on this page. Your continued use of the F3UTER App will mean that you accept such changes. If the user does not agree to the revised Terms, please stop using the F3UTER App."), br(),
                                     
                                     h4(strong("Licence")),
                                     h4("Copyright 2020 Astex Therapeutics Ltd."),
                                     p("This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version."),
                                     p("This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details (", a(href="https://www.gnu.org/licenses/gpl-3.0.en.html", "https://www.gnu.org/licenses/gpl-3.0.en.html", target="_blank"), ").")
                                     
                                     
                     ))
            ),
            
            
            
            # main tab 6 
            tabPanel(title = "Contact",
                     fluidRow(
                       column(width = 10, offset = 1, br(), h4("F3UTER is developed by Astex Pharmaceuticals in collaboration with  Ryten lab and Botia lab, under an Astex Sustaining Innovation Postdoctoral Fellowship"),hr()
                       ),
                       column(width = 2, offset = 1,
                                     h3(strong("Astex Pharmaceuticals")),
                                     p("436 Cambridge Science Park"),
                                     p("Cambridge, CB4 0QA"),
                                     a(href="https://astx.com", "https://astx.com", target="_blank")
                       ),
                       column(width = 3, offset = 1, 
                                     h3(strong("Ryten Lab")),
                                     p("UCL Great Ormond Street Institute of Child Health"),
                                     p("30 Guilford Street, London WC1N 1EH"),
                                     a(href="http://www.rytenlab.com", "http://www.rytenlab.com", target="_blank")
                       ),
                       column(width = 3, offset = 1, 
                              h3(strong("Juan A. Botía")),
                              p("Intelligent Systems and Telematics research group"),
                              p("Department of Information and Telecommunications Engineering"),
                              p("University of Murcia, Spain"),
                              a(href="https://www.um.es/gsit/?lg=en", "https://www.um.es/gsit/?lg=en", target="_blank")
                       ),
                       column(width = 10, offset = 1, 
                                    br(), br(), hr(), br(),
                                     h4("For questions related to this resource, please contact ", strong("Siddharth Sethi"), "-", a(href="mailto:siddharth.sethi@astx.com","siddharth.sethi@astx.com")), br(), 
                                     h4("Source code: ", a(href="https://github.com/sid-sethi/F3UTER/F3UTER_app", "https://github.com/sid-sethi/F3UTER/F3UTER_app", target="_blank")), br(), br(), br(), br(), br(),
                              p(img(src="astexlogo.gif", width="15%"), img(src="ucl_logo.png", width="15%"), img(src="umu_logo.jpg", width="14%"), align="center")
                       )
                     
                  )
            )
            
            
      ) # navbar page
) # fluid page





server <- function(input, output, session) {
  
  waiter_hide()
  
  
  ############ SUBMIT GENE ############################
  # extracting gene-er data
  gene_er_data <- eventReactive(input$submit_gene, {
    
    
    validate(need(input$gene_input != "", "Please input gene symbol or ENSG ID"))
    
    if(str_detect(input$gene_input, "ENSG")){
      gene_id <- input$gene_input
      geneName <- glue_sql("SELECT gene_name FROM geneData WHERE gene_id = {input$gene_input}", .con = feuterdb) %>% 
        dbGetQuery(feuterdb, .) %>% 
        .$gene_name
    }else{
      geneName <- input$gene_input
      gene_id <- glue_sql("SELECT gene_id FROM geneData WHERE gene_name = {geneName} LIMIT 1", .con = feuterdb) %>% 
        dbGetQuery(feuterdb, .) %>% 
        .$gene_id
    }
    
    
    # if(length(gene_id) == 0){
    #   stop(str_c(input$gene_input, " not found in the ensembl database (v94)"))
    # }
    validate(need(gene_id != 0, str_c("ERROR: ", input$gene_input, " not found in the ensembl database (v94). Please input a valid gene identifier")))
    
    
    # outout gene identifier info
    output$gene_output <- renderText({
      paste("Gene name: ", geneName, "\nGene ID: ", gene_id, sep="")
    })
    
    
    # show help text for next steps
    shinyjs::show("select_row_text")  
    shinyjs::show("submit_er")
    
    
    extract_gene_er_data(gene_id = gene_id, tissue_input = input$tissue_input, datadb = feuterdb)
    
    
  }) # gene_er_data end
  ########### submit gene end ####################  
  
  
  
  ################# GENE-TABLE OUTPUT ###############
  output$genetable <- DT::renderDataTable(
    gene_er_data(), width = "30%", class = 'cell-border stripe compact hover', 
    selection='single'
  )
  ################## END ###########################
  
  
  ######################## SUBMIT ER ###############
  er_data <- eventReactive(input$submit_er, {
    
    validate(
      need(input$genetable_rows_selected != "", "Please select a row from the Gene-table and submit ER")
    )

    sel = input$genetable_rows_selected
    er_id = gene_er_data() %>% dplyr::slice(sel) %>% .$id

    extract_er_data(er_id = er_id, datadb = feuterdb)

  }) # er_data end

  

  # outout ER selection info
  output$er_output <- renderText({
    sel = input$genetable_rows_selected
    if(length(sel)) {
      er_id = gene_er_data() %>% dplyr::slice(sel) %>% .$id
      str_c("ER_id = ", er_id)
    }
  })
  
  
  
  
  ################# ER FEATURE PLOT OUTPUT ###############
  output$featureResult <- renderPlot({
    
    input$submit_er
    
    isolate({
      
    sel = input$genetable_rows_selected
    if(length(sel) & input$submit_er != 0) {

      er_id = gene_er_data() %>% dplyr::slice(sel) %>% .$id
      
      waitress <- Waitress$new("#featureResult", hide_on_render = TRUE, theme = "overlay-percent") # call the waitress
      waitress$auto(value = 1, ms = 15)
      
      #gc()
    

    }else{
      validate(
        need(input$submit_er != 0 & input$genetable_rows_selected != "", "Please select a row from the Gene-table and submit ER")
      )
    }

    #featurePlot <- plot_features(df = er_data())

    # call plot()
    #suppressWarnings(plot(featurePlot))
    plot_features(df = er_data()) %>% plot() %>% suppressWarnings()

    # prepare download
    shinyjs::show("downloadFeature")
    output$downloadFeature <- downloadHandler(
      filename = function() {
        str_c(input$gene_input, "_", er_id, "features.png")
      },
      content = function(file) {
        png(file, bg="transparent",units="in",width = 10.25, height= 7.75 ,res=600)
        suppressWarnings(plot(
          #featurePlot +
          plot_features(df = er_data()) +
            theme(strip.text = element_text(size=7, face="bold"),
                  axis.text.y = element_text(size = 7),
                  legend.text = element_text(size = 10)
            )
        ))
        dev.off()
      }
    )
    })

  })
  ################## ER FEATURE PLOT END ############################
  
  
  
  ####################### DOWNLOAD TABLES ########################
  #######  gene-table #######
  # activate gene table download
  output$download_gene_table <- renderUI({
    
    validate(need(input$gene_input != "" & input$submit_gene != 0, "Please input gene identifier to activate this download"))
    
    waitress <- Waitress$new("#download_gene_table", hide_on_render = TRUE, theme = "overlay-radius") 
    waitress$auto(value = 1, ms = 2)
    
    if(nrow(gene_er_data()) > 0) {
      downloadButton("Gene_table_csv", label = "Gene-table", style="color: #fff; background-color: #000000; border-color: #000000")
    }else{
      h6("Gene-table output is empty - no data generated for download")
    }
    
  })
  
  # gene-table csv download data
  output$Gene_table_csv <- downloadHandler(
    filename = function() { str_c(input$gene_input, "_erTable.csv")},
    content = function(file) { write.csv(gene_er_data(), file, row.names = FALSE)}
  )
  
  
  ####### ER feature table #######
  # activate er feature table download
  output$download_er_feature_table <- renderUI({
    
    validate(need(input$gene_input != "" & input$submit_gene != 0, "Please input gene identifier to activate this download"))
    
    waitress <- Waitress$new("#download_er_feature_table", hide_on_render = TRUE, theme = "overlay-radius") # call the waitress
    waitress$auto(value = 1, ms = 2)
    
    if(nrow(gene_er_data()) > 0) {
      downloadButton("er_feature_table_csv", label = "ER feature table", style="color: #fff; background-color: #000000; border-color: #000000")
    }else{
      h6("Gene-table output is empty - no data generated for download")
    }
    
  })
  
  # er feature table csv download data
  output$er_feature_table_csv <- downloadHandler(
    filename = function() { str_c(input$gene_input, "_erFeatureTable.csv")},
    content = function(file) { 
      
      #erFeatureTable <- dplyr::semi_join(erData, gene_er_data(), by = "id")
      erFeatureTable <- glue_sql("SELECT * FROM erData WHERE id IN ({gene_er_data()$id*})", .con = feuterdb) %>% 
        dbGetQuery(feuterdb, .)
      
      write.csv(erFeatureTable, file, row.names = FALSE)
      rm(erFeatureTable)
      }
  )
  
  
  #######  ER feature plot data frame #######
  # activate er_feature_plot df download
  output$download_er_feature_plot_df <- renderUI({
    
    validate(need(input$gene_input != "" & input$submit_gene != 0, "Please input gene identifier to activate this download"))
    
    validate(need(input$submit_er != 0 & input$genetable_rows_selected != "", "Please submit an ER to activate this download"))
    
    waitress <- Waitress$new("#download_er_feature_plot_df", hide_on_render = TRUE, theme = "overlay-radius") 
    waitress$auto(value = 1, ms = 1.5)
    
    if(nrow(er_data()) > 0) {
      downloadButton("er_feature_plot_df_csv", label = "ER feature plot df", style="color: #fff; background-color: #000000; border-color: #000000")
    }else{
      h6("ER feature plot data frame output is empty - no data generated for download")
    }
    
  })
  
  # er feature plot df csv download
  output$er_feature_plot_df_csv <- downloadHandler(
    filename = function() { 
      er_id = gene_er_data() %>% dplyr::slice(input$genetable_rows_selected) %>% .$id
      str_c(input$gene_input, "_", er_id, "FeaturePlotDf.csv")
      },
    content = function(file) { write.csv(er_data() %>% dplyr::select(-dummy), file, row.names = FALSE)}
  )
  
  ##########################################################
  
  
  

  ########################## OBSERVE EVENTS ########################
  
  observeEvent(input$submit_er, {
    updateTabsetPanel(session = session, inputId = "tabs", selected = "ER genomeView plot")
    updateSliderInput(session = session, inputId = "zoom_factor_input", value = 30)
  }, ignoreNULL = FALSE)
  
  
  observeEvent(input$submit_gene, {
    updateTabsetPanel(session = session, inputId = "tabs", selected = "Gene-table")
    shinyjs::hide("zoom_factor_input")
    shinyjs::hide("downloadFeature")
    shinyjs::hide("downloadGenomeView")
    gc()
  }, ignoreNULL = FALSE)
  
  
  
  
  ################################## Data Downloads #################################
  
  output$download_er_data <- downloadHandler(
    filename = "all_ER_predictions.csv",
    content = function(file) {
      #write.csv(dbGetQuery(feuterdb, 'SELECT * FROM er_raw_all'), file, row.names = FALSE)
      file.copy('data/all_ER_predictions.csv', file)
    }
  )
  
  output$download_er_feature_matrix <- downloadHandler(
    filename = "ER_feature_matrix.csv",
    content = function(file) {
      #write.csv(dbGetQuery(feuterdb, 'SELECT * FROM erData'), file, row.names = FALSE)
      file.copy('data/ER_feature_matrix.csv', file)
    }
  )
  
  output$download_training_feature_matrix <- downloadHandler(
    filename = "training_feature_matrix.csv",
    content = function(file) {
      file.copy('data/training_feature_matrix.csv', file)
    }
  )
  #######################################################
  
  
  
  #################### ER GENOME VIEW PLOT ################
  
  output$genomeView <- renderPlot({
    
    input$submit_er
    input$zoom_factor_input
    
    isolate({
      
      sel = input$genetable_rows_selected
      if(length(sel) & input$submit_er != 0) {
        
        er_id = gene_er_data() %>% dplyr::slice(sel) %>% .$id
        
        waitress <- Waitress$new("#genomeView", hide_on_render = TRUE, theme = "overlay-percent") 
        waitress$start()
        
        gc()
        
        # progress <- shiny::Progress$new(min = 0, max=10)
        # on.exit(progress$close())
        # progress$set(value = 0, message = paste("Making genome view plot... ", er_id, sep=""))
        
      }else{
        validate(
          need(input$submit_er != 0 & input$genetable_rows_selected != "", "Please select a row from the Gene-table and submit ER")
        )
      }
      
      
      if(str_detect(input$gene_input, "ENSG")){
        gene_id <- input$gene_input
        geneName <- glue_sql("SELECT gene_name FROM geneData WHERE gene_id = {input$gene_input}", .con = feuterdb) %>% 
          dbGetQuery(feuterdb, .) %>% 
          .$gene_name
      }else{
        geneName <- input$gene_input
        gene_id <- glue_sql("SELECT gene_id FROM geneData WHERE gene_name = {geneName} LIMIT 1", .con = feuterdb) %>% 
          dbGetQuery(feuterdb, .) %>% 
          .$gene_id
      }
      
      
      # create named list of gene/tx id
      gene_tx_id_list <- list(gene_id)
      names(gene_tx_id_list) <- "gene_id"
      # get gene cordinates
      gene_tx_to_plot <- GenomicFeatures::genes(txdb_v92, filter = gene_tx_id_list)
      
      
      # get plot coordinates
      coords_to_plot <- .get_coords_to_plot(gene_tx_to_plot, er_id, geneName, zoom_factor = input$zoom_factor_input)
      coords_to_plot_gr <- GRanges(seqnames = coords_to_plot[["chr"]], ranges = IRanges(start = coords_to_plot[["x_min"]], end = coords_to_plot[["x_max"]]))
      
      #progress$set(value = 1)
      waitress$inc(10)
      
      # get gene exons to plot
      gene_exons_to_plot <- .get_gene_exons_to_plot(gene_tx_id_list, coords_to_plot_gr, ref = txdb_v92)
      rm(coords_to_plot_gr)
      
      #progress$set(value = 2)
      waitress$inc(10)
      
      # get ERs to plot
      ers_to_plot <- dbGetQuery(feuterdb, 'SELECT seqnames, start, end, width, strand,  [predicted.response] FROM er_raw_all WHERE tissue = ? AND seqnames = ? AND start >= ? AND end <= ?', params = c(coords_to_plot[["er_tissue"]], coords_to_plot[["chr"]], coords_to_plot[["x_min"]], coords_to_plot[["x_max"]]))
      
      #progress$set(value = 3)
      waitress$inc(10)
      
      # get PAS to plot
      pas_to_plot <- dbGetQuery(feuterdb, 'SELECT * FROM pa WHERE chr = ? AND start >= ? AND end <= ?', params = c(coords_to_plot[["chr"]], coords_to_plot[["x_min"]], coords_to_plot[["x_max"]]))
      
      #progress$set(value = 4)
      waitress$inc(10)
      
      # gene track
      genome_view_plot <- .plot_gene_track(coords_to_plot, gene_exons_to_plot)
      rm(gene_exons_to_plot)
      
      #progress$set(value = 5)
      waitress$inc(10)
      
      # add ers to plot
      genome_view_plot <- .plot_add_ers(coords_to_plot, genome_view_plot, ers_to_plot)
      rm(ers_to_plot)
      
      #progress$set(value = 6)
      waitress$inc(10)
      
      # add PAS to plot
      genome_view_plot <- .plot_add_pas(coords_to_plot, genome_view_plot, pas_to_plot)
      rm(pas_to_plot)
      #gc()
      
      #progress$set(value = 7)
      waitress$inc(10)
      
      # calculate junction points to plot
      junctions_points <- get_junction_points(er_id, feuterdb, coords_to_plot)
      
      
      #progress$set(value = 8)
      waitress$inc(10)
      
      # add junctions to plot
      genome_view_plot <- .plot_add_junctions(coords_to_plot, genome_view_plot, junctions_points)
      rm(junctions_points)
      
      #progress$set(value = 9)
      waitress$inc(10)
      
      # add annotation to plot
      genome_view_plot <- .plot_annotation(genome_view_plot, coords_to_plot)
      rm(coords_to_plot)
      
      # call plot()
      suppressWarnings(plot(genome_view_plot))
      #progress$set(value = 10)
      waitress$inc(10)

      # slider input
      shinyjs::show("zoom_factor_input")
      
      
      # prepare download
      shinyjs::show("downloadGenomeView")
      output$downloadGenomeView <- downloadHandler(
        filename = function() {
          str_c(input$gene_input, "_", er_id, "genomeView.png")
        },
        content = function(file) {
          png(file, bg="transparent",units="in",width = 9.25, height= 4.75 ,res=600)
          suppressWarnings(plot(
            genome_view_plot
          ))
          dev.off()
        }
      )
      
      
    })
  

  })

  
  

  
} # end server


# Run the application 
shinyApp(ui = ui, server = server)

