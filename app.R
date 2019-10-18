library(shiny)
library(tidyverse)
library(biophysr)
library(shinythemes)
library(gtools)
library(gridExtra)
library(dygraphs)
library(shinycssloaders)


# Define UI for app  ----
ui <- fluidPage(theme = shinytheme("slate"),

  navbarPage("Muscle Biophysics Lab: Data Analyzer v2.0",
             
    #Unregulated Motility Panel
    tabPanel("Motility Data",
       sidebarLayout(
        sidebarPanel(
           fileInput(inputId = "file",
                  label = "Upload Motility Files (.xls)",
                  multiple = TRUE,
                  accept = ".xls",
                  buttonLabel = "Browse...",
                  placeholder = "No files selected"),

           radioButtons(inputId = "myo_type",
                   label = "Myosin Type:",
                   choiceNames = list("Skeletal (4um cutoff)", "Cardiac (1um cutoff)"),
                   choiceValues = list("skeletal", "cardiac"),
                   inline = TRUE),

           actionButton(inputId = "action_button",
                   class="btn btn-info", 
                   label = "Analyze",
                   width = 175),

           downloadButton("download_by_video",
                   class="btn btn-success", 
                   "Download 'By Video'",
                   width = 175),

           downloadButton("download_summary",
                   class="btn btn-success", 
                   "Download Summary",
                   width = 175)
                   ), #sidebar panel


    
        mainPanel(

            tabsetPanel(type = "tabs",
                  tabPanel("By Video", tableOutput("video_table_output")),
                  tabPanel("Summary", tableOutput("summary_table_output"))

                 )
             ) #main panel close
         ) #side bar layout close
     ), #motility tab Panel close
 
 
 
 
 
     #unregulated Motility plot Tab
     tabPanel("Motility Plot",
              tabPanel("Motility Plot",
                       sidebarLayout(
                         sidebarPanel(
                           fileInput(inputId = "file",
                                     label = "Upload Summary Motility File (.xls)",
                                     multiple = TRUE,
                                     accept = ".xls",
                                     buttonLabel = "Browse...",
                                     placeholder = "No files selected")
                         ), #sidebar panel
                           mainPanel(
                             "Plot", plotOutput("motility_plot")

                             
                           ) #main panel close
                         ) #sidebar panel close
                       ) #tab panel layout close
              ), #motility tab Panel close
           
 
 
 
 #regulated Motility plot Tab
 tabPanel("Regulated Motility Plot",
          tabPanel("Regulated Motility Plot",
                   sidebarLayout(
                     sidebarPanel(
                       fileInput(inputId = "file",
                                 label = "Upload Summary Motility File (.xls)",
                                 multiple = TRUE,
                                 accept = ".xls",
                                 buttonLabel = "Browse...",
                                 placeholder = "No files selected")
                     ), #sidebar panel
                     mainPanel(
                       "Plot", plotOutput("regulated_motility_plot")
                       
                       
                     ) #main panel close
                   ) #sidebar panel close
          ) #tab panel layout close
 ), #motility tab Panel close
 
 
 
 
 
  #Mini_ensmble trap Tab
  tabPanel("Mini-Ensemble Trap",
          
           sidebarLayout(
             sidebarPanel(
             fileInput(inputId = "mini_file",
                     label = "Upload Mini-Ensemble Data (.txt)",
                     multiple = FALSE,
                     accept = ".txt",
                     buttonLabel = "Browse...",
                     placeholder = "No files selected"),
             
           #  actionButton(inputId = "view_mini_raw",
            #              label = "View Raw Trace"),
             
             checkboxInput(inputId = "detrend_mini",
                           label = "Detrend",
                           value = FALSE),
              actionButton(inputId = "mini_action_button",
                       label = "Analyze")
           
          
             ),
                  
             
             mainPanel(
               
               tabsetPanel(type = "tabs",
                           tabPanel("Raw Plot", dygraphOutput("explore_mini_raw") %>% withSpinner()),
                           tabPanel("Measured Events", tableOutput("mini_table_output") %>% withSpinner()),
                           tabPanel("Plot", plotOutput("mini_plot") %>% withSpinner())
                           
               )
             ) #main panel close
    
             )
  ), #mini-ensemble tab close
  
  tabPanel("Single-Molecule Trap",
           print("This is for single molecule analysis")
           ) #single molecule tab panel close
        
)#navbar page close
)#fluid page close



# Define server logic required to draw a histogram ----
server <- function(input, output) {

  options(shiny.maxRequestSize=30*1024^2)


  data_by_video <- eventReactive(input$action_button,{


    filenames <- input$file$name
    data_upload <- map(input$file$datapath, read.delim)
    for (i in seq_along(data_upload)){
      data_upload[[i]]<-cbind(data_upload[[i]], filenames[i])}

    combined <- bind_rows(data_upload) %>%
      separate('filenames[i]', c("Condition", "Rest"), sep = "_00") %>%
      separate("Rest", c("Video", "Extenstion"), sep = "_")

    if(input$myo_type == "skeletal"){

    data_analyzed <- combined %>%
      group_by(Condition, Video)%>%
      analyze_motility() %>%
      rename("Velocity" = average_velocity,
             "Percent_Moving" =  percent_moving,
             "Moving_Filaments" = moving_filaments,
             "Total_Filaments" = total_filaments)
    } else if(input$myo_type == "cardiac"){


      data_analyzed <- combined %>%
        group_by(Condition, Video)%>%
        analyze_cardiac_motility() %>%
        rename("Velocity" = average_velocity,
               "Percent_Moving" =  percent_moving,
               "Moving_Filaments" = moving_filaments,
               "Total_Filaments" = total_filaments)

    }


  })

  output$video_table_output <- renderTable({
    if(is.null(input$file)) return(tibble("No data currently selected" =  "Please use 'Browse...' to upload"))
    else return(data_by_video())
  })



  data_summary <- eventReactive(input$action_button, {

    filenames <- input$file$name
    data_upload <- map(input$file$datapath, read.delim)
    for (i in seq_along(data_upload)){
      data_upload[[i]]<-cbind(data_upload[[i]], filenames[i])}

    combined <- bind_rows(data_upload) %>%
      separate('filenames[i]', c("Condition", "Video"), sep = "_00")

    if(input$myo_type == "skeletal"){

    data_analyzed <- combined %>%
      group_by(Condition, Video)%>%
      analyze_motility() %>%
      rename("Velocity" = average_velocity,
             "Percent_Moving" =  percent_moving,
             "Moving_Filaments" = moving_filaments,
             "Total_Filaments" = total_filaments)

    } else if(input$myo_type == "cardiac"){

      data_analyzed <- combined %>%
        group_by(Condition, Video)%>%
        analyze_cardiac_motility() %>%
        rename("Velocity" = average_velocity,
               "Percent_Moving" =  percent_moving,
               "Moving_Filaments" = moving_filaments,
               "Total_Filaments" = total_filaments)
    }

   summary_analysis <- data_analyzed %>%
     group_by(Condition) %>%
     summarize("Average Velocity" = mean(Velocity),
               "Average Percent Moving" = mean(Percent_Moving))


  })

  output$summary_table_output <- renderTable({
    if(is.null(input$file)) return(NULL)
          else return(data_summary())
  })


  output$download_by_video <- downloadHandler(
    filename = function() {
      paste("motility_analyzed_videos.csv", sep = "")
    },
    content = function(file1) {
      write.csv(data_by_video(), file1, row.names = FALSE)
    })

  output$download_summary <- downloadHandler(
    filename = function() {
      paste("motility_summary.csv", sep = "")
    },
    content = function(file1) {
      write.csv(data_summary(), file1, row.names = FALSE)
    })
  
  ###################### MINI ANALYSIS ######################
  
  
 # mini_raw_data <- reactive({
 ## mini <- scan(input$mini_file$datapath)
  #})
  
  explore_mini <- reactive({ 
    
    mini <- scan(input$mini_file$datapath)
    
    if(input$detrend_mini == FALSE){
    
    raw_mini_trace <- data.frame(x = 1:length(mini),
                                 y = mini)
    

    } else if(input$detrend_mini == TRUE){
      
      mini_detrend <- detrend_trap(mini, 25000)
    
      raw_mini_trace <- data.frame(x = 1:length(mini_detrend),
                                 y = mini_detrend)
      
    }
    
    
    mini_raw_plot <- dygraph(raw_mini_trace) %>% 
      dyRangeSelector() %>% 
      dyOptions(colors = "black") 
    
   return(mini_raw_plot)
    
    
  })
  


  
  
  output$explore_mini_raw <- renderDygraph({
    if(is.null(input$mini_file)){
      
      return(NULL)
      
    } else {
      
      return(explore_mini())
    }
  })
  
  mini_data <- eventReactive(input$mini_action_button,{
  
    
   #mini <- scan(input$mini_file$datapath)
    
    mini <- scan(input$mini_file$datapath)
    
    if(input$detrend_mini == FALSE){
      
      analyzed_trap <- analyze_mini_ensemble(mini)
      
      
    } else if(input$detrend_mini == TRUE){
      
      mini_detrend <- detrend_trap(mini, 25000)
      
      analyzed_trap <- analyze_mini_ensemble(mini_detrend)
      
    }

    return(analyzed_trap)
    
    
  })


  plot_mini <- eventReactive(input$mini_action_button,{ 
    
    mini <- scan(input$mini_file$datapath)
    
    
    if(input$detrend_mini == FALSE){
      
      plots <- plot_analyze_mini_ensemble(mini)
      
      
      
    } else if(input$detrend_mini == TRUE){
      
      mini_detrend <- detrend_trap(mini, 25000)
      plots <- plot_analyze_mini_ensemble(mini_detrend)
      
      
    }
    

    plot(plots)
  
 
        
    
    
  })
  


  output$mini_table_output <- renderTable({
    if(is.null(input$mini_file)){
      
    return (tibble("No data currently selected" =  "Please use 'Browse...' to upload"))
      
    } else {
     return(mini_data())
    }
  })
  
 output$mini_plot <- renderPlot({
    if(is.null(input$mini_file)){
      
     return(NULL)
      
    } else {
      
    return(plot_mini())
    }
 })

  
  
  ###################### MINI END ##################
  
  
  
} # server close

shinyApp(ui = ui, server = server)


