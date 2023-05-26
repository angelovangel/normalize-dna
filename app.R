library(rhandsontable)
library(shiny)
library(shinydashboard)
library(shinyjs)
#library(shinyalert)
library(plater)
library(tibble)
library(stringr)
library(reactable)
library(dplyr)
library(rmarkdown)
library(curl)
library(waiter)

wells_colwise <- lapply(1:12, function(x) {str_c(LETTERS[1:8], x)}) %>% unlist()

# creates base empty dataframe to view and fill later
make_dest <- function() {
    dest <- tibble(source_well = wells_colwise, dest_well = wells_colwise,
                   conc = NA, dna_size = NA, vol1 = NA, vol2 = NA,
                   ng = NA, fmoles = NA)
    dest
}

example_table <- make_dest()
example_table$dna_size = c(10000, 20000, 20000, rep(NA, 93))
example_table$conc = c(120, 80, 35, rep(NA, 93))
example_table$vol1 = c(12, 8, 35, rep(NA, 93))

tab1 <-  fluidRow(
  box(width = 12, height = 2800, status = "info", solidHeader = FALSE, 
      title = "Noramlize DNA Opentrons protocol", collapsible = F,
      fluidRow(
        column(12, tags$p('This app generates an OT-2 Python script 
                          to normalize DNA (by target conc or molarity) or to transfer 
                          liquid from source to destination plate.')
               ),
        
        column(2, selectizeInput('protocol_type', 'Select protocol', 
                                 choices = list('Normalize nanograms' = 'normalize_ng', 
                                                'Normalize fmoles' = 'normalize_fmoles', 
                                                'Simple transfer' = 'transfer'), 
                                 selected = 'normalize_ng')),
        column(2, uiOutput('target_amount')),
        column(2, numericInput('target_vol', 'Final volume', min = 1, max = 200, value = 10)),
        column(1),
        column(2, selectizeInput('left_m', 'Left pipette mount', 
                                 choices = c('p20_single_gen2', 'p20_multi_gen2'), 
                                 selected = 'p20_single_gen2')
        ),
        column(2, selectizeInput('right_m', 'Right pipette mount', 
                                 choices = c('p20_single_gen2', 'p20_multi_gen2'), 
                                 selected = 'p20_multi_gen2')
        )
      ),
      fluidRow(
        column(2, actionButton('deck', 'Show deck layout', width = '100%', style = 'margin-top:20px')),
        column(3, downloadButton('download_script', 'Opentrons script', width = '100%', style = 'margin-top:20px'),
                  downloadButton('download_samples', 'Sample sheet', width = '100%', style = 'margin-top:20px')
               )
      ),
      tags$hr(),
      column(5, 
             tags$p(HTML("Source plate - <b>vol1</b> is DNA (source plate), <b>vol2</b> is water up to final volume")),
             rHandsontableOutput('hot')),
      column(7, 
             tags$p("Reaction plate preview"),
             reactableOutput('plate'), 
             tags$hr()
             )
  )
)

tab2 <- fluidRow(
  box(width = 12, status = "info", solidHeader = FALSE, title = "Opentrons protocol preview", collapsible = F,
      verbatimTextOutput('protocol_preview')
      )
)

ui <- dashboardPage(skin = 'yellow',
  #useShinyalert(),
  
  header = dashboardHeader(title = 'Normalize DNA by concentration or molarity', titleWidth = 800),
  sidebar = dashboardSidebar(disable = T),
  body = dashboardBody(
   useShinyjs(),
   tabsetPanel(
     tabPanel(title = "Enter samples data", icon = icon("list"), tab1),
     tabPanel(title = "Opentrons script preview", icon = icon('code'), tab2)
   )
  )
)
  

# server #
server = function(input, output, session) {
  
  ### read template
  protocol_url <- "https://raw.githubusercontent.com/angelovangel/opentrons/main/protocols/06-normalize-dna.py"
  
  if (curl::has_internet()) {
    waiter_show(html = spin_wave())
    con <- url(protocol_url)
    protocol_template <- readLines(con, warn = F)
    close(con)
    waiter_hide()
  } else {
    protocol_template <- readLines('06-normalize-dna.py', warn = F)
  }
  
  
  ### REACTIVES
    protocol <- reactiveValues(rxn_vol = 0, sample_vol = 0, total_fmoles = 0, total_ng = 0)
    
    # target ng
    hot <- reactive({
      if(!is.null(input$hot)) {
        as_tibble(hot_to_r(input$hot)) %>%
          mutate(
            vol1 = case_when(
              input$protocol_type == 'normalize_ng' ~ input$amount/conc,
              input$protocol_type == 'normalize_fmoles' ~ input$amount/(conc/((dna_size*617.96) + 36.04) * 1000000),
              input$protocol_type == 'transfer' ~ vol1,
              .default = NA
            )
          ) %>%
          mutate(
            vol1 = case_when(
              vol1 < 0 ~ 0,
              (vol1 < 0.5 & vol1 > 0) ~ 0.5,
              (input$protocol_type == 'transfer' & vol1 <= 200) ~ vol1,
              vol1 > 200 ~ 200,
              vol1 > input$target_vol ~ input$target_vol,
              TRUE ~ vol1
              ),
            vol2 = case_when(
              vol2 < 0 ~ 0,
              (input$protocol_type == 'transfer' & vol2+vol1 <= 200) ~ vol2,
              (input$protocol_type == 'transfer' & vol2+vol1 > 200) ~ 200 - vol1,
              TRUE ~ input$target_vol - vol1
              ),
            fmoles = (conc/((dna_size*617.96) + 36.04) * 1000000) * vol1,
            ng = conc*vol1
          )
          
      } else {
         example_table
        }
    })
  
    
    plate <- reactive({
      if(!is.null(input$hot)) {
        df <- hot() %>% 
          mutate(label = str_c(source_well, "<br>", round(vol1 + vol2, 0), ' ul'))
        
        plater::view_plate(
          #hot_to_r(input$hot) , 
          df,
          well_ids_column = 'dest_well', columns_to_display = c('label')
        )
      } else {
        plater::view_plate(example_table %>% 
                             mutate(label = str_c(source_well, "<br>", conc)), 
                           well_ids_column = 'dest_well', columns_to_display = c('label'))
      }
    })
    
  myvalues <- reactive({
    
    source_wells <- hot()$source_well %>% str_replace_na(replacement = ' ')
    volume1 <- str_replace_na(hot()$vol1, '0') # replace NA with 0, gDNA
    # water wells is always A1
    volume2 <- str_replace_na(hot()$vol2, '0') # water
      c(
        str_flatten(source_wells, collapse = "','"),  
        str_flatten(volume1, collapse = ", "),
        str_flatten(volume2, collapse = ", ")
        ) 
  })
      
  myprotocol <- reactive({
    
    str_replace(string = protocol_template, 
                pattern = 'sourcewells=.*', 
                replacement = paste0("sourcewells=['", myvalues()[1], "']")) %>%
      str_replace(pattern = 'volume1=.*', replacement = paste0('volume1=[', myvalues()[2], ']')) %>%
      str_replace(pattern = 'volume2=.*', replacement = paste0('volume2=[', myvalues()[3], ']')) %>%
      str_replace(pattern = 'left_mount =.*', replacement = paste0('left_mount = ', "'", input$left_m, "'")) %>%
      str_replace(pattern = 'right_mount =.*', replacement = paste0('right_mount = ', "'", input$right_m, "'"))
    })
  
      
  ### OBSERVERS
  
  observe({
    if(input$protocol_type == 'transfer') {
      shinyjs::hide('target_vol')
      shinyjs::hide('target_amount')
    } else {
      shinyjs::show('target_vol')
      shinyjs::show('target_amount')
    }
  })
  
  observeEvent(input$deck, {
      showModal(
        modalDialog(title = 'Opentrons deck preview',
                    HTML('<img src="deck.png">'),
                    size = 'l', easyClose = T, 
        )
      )
  })
    
  ## OUTPUTS
    
  output$target_amount <- renderUI({
    switch(input$protocol_type,
           'normalize_ng' = numericInput('amount', 'Target ng per well', min = 1, max = 5000, value = 200),
           'normalize_fmoles' = numericInput('amount', 'Target fmoles per well', min = 1, max = 5000, value = 20)
           #'transfer' = numericInput('amount', 'Transfer volume', min = 0.5, max = 200, value = 10)
              )
  })
    
  renderer <- function() {
    
    "function(instance, td, row, col, prop, value, cellProperties) {
    Handsontable.renderers.NumericRenderer.apply(this, arguments);
    tbl = this.HTMLWidgets.widgets[0]
    
    if (value >= 200 || value <= 0.5) {
    td.style.color = 'red'
      }
    }
    "}
    # renders first column well in grey for better plate overview
    rendergrey <- function() {
      "
    function(instance, td, row, col, prop, value, cellProperties) {
      Handsontable.renderers.TextRenderer.apply(this, arguments);
      
      tbl = this.HTMLWidgets.widgets[0]
      //hrows = tbl.params.greylines
      hrows = [0 ,8, 16, 24, 32, 40, 48, 56, 64, 72, 80, 88]
      hrows = hrows instanceof Array ? hrows : [hrows] 

      
      if (hrows.includes(row)) {
        td.style.background = '#D6EAF8';
      }
      
      return td;
  }"
    }
    
    output$hot <- renderRHandsontable({
      df <- hot()
      if(!is.null(df)) {
      rhandsontable(df,
                    stretchH  = 'all',  
                    #svol = 9,
                    height = 2800,
                    rowHeaders = NULL) %>%
        #hot_cols(renderer = rendergrey() ) %>%
        hot_col('ng', readOnly = T, type = 'numeric', format = '0.0') %>%
        hot_col('fmoles', readOnly = T, type= 'numeric', format = '0.0') %>%
        # depends on protocol
        hot_col('vol1', readOnly = if_else(input$protocol_type == 'transfer', F, T), 
                type = 'numeric', allowInvalid = F, renderer = renderer() ) %>% # highlight volumes > max
        hot_col('vol2', readOnly = if_else(input$protocol_type == 'transfer', F, T), 
                type = 'numeric', allowInvalid = F, renderer = renderer() ) %>% # highlight volumes > max
        hot_col('dna_size', format = '0') %>%
        #hot_cell(1, 3, 'test') %>%
        hot_validate_numeric('conc', min = 1, max = 5000, allowInvalid = T) %>%
        hot_col('dest_well', readOnly = T, renderer = rendergrey()) %>%
        hot_col('source_well', readOnly = F, type = 'dropdown', source = wells_colwise, allowInvalid = F, strict = T) 
        }
    })
    
    
    output$plate <- renderReactable({
      reactable(plate()$label, 
                highlight = T, wrap = F, 
                bordered = T, compact = T, fullWidth = T, sortable = F,
                defaultColDef = colDef(minWidth = 50,
                  html = TRUE, 
                  headerStyle = list(background = "#f7f7f8", fontSize = '80%')
                  # color barcodes if duplicate
                  # style = function(value) {
                  #   myvalue <- str_extract(value, 'barcode.*')
                  #   mydf <- hot()
                  #   if(!is.na(myvalue)) {
                  #     textcolor <-mydf$mycolor[match(myvalue, mydf$barcode)]
                  #   } else {
                  #     textcolor <- 'black'
                  #   }
                  #   #print(myvalue)
                  #   list(color = textcolor, fontSize = '80%')
                  # }
                  )
                )
    })
    
    output$protocol_preview <- renderPrint({
      write(myprotocol(), file = "")
    })
    
  ### DOWNLOADS
  output$download_script <- downloadHandler(
      filename = function() {
        paste0(format(Sys.time(), "%Y%m%d-%H%M%S"), '-normalize-dna.py')
        },
      content = function(con) {
        # at download time, replace name so that it appears on the Opentrons app
        replacement <- paste0(format(Sys.time(), "%Y%m%d-%H%M%S"), '-normalize-dna.py')
        write(myprotocol() %>%
                str_replace(pattern = "06-normalize-dna.py", 
                            replacement = replacement), 
              con)
      }
    )
    
  output$download_samples <- downloadHandler(
      filename = function() {
        paste0(format(Sys.time(), "%Y%m%d-%H%M%S"), '-samplesheet.csv')
      },
      content = function(con) {
        write.csv(hot(), con, row.names = F)
      }
      )
}
  
  
shinyApp(ui, server)