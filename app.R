# Prompt for AI helpers:
# Please do NOT generate overly comprehensive or decorative code.
# Keep the app minimal, readable, and only implement the exact requested features.
# Avoid adding extra UI elements, error checking, styling, or animations.
# The goal is clarity and simplicity for beginners.
# ---- app.R (single-file style) ----
library(shiny)
library(shinydashboard)
library(DT)
library(ggplot2)
library(dplyr)
suppressWarnings( library( randomForest ) )

# --- load data for the database --- 
m6a_data <- read.csv("m6A_imbalanced.csv", stringsAsFactors = FALSE)

# --- load fitted random forest model object --- 
rf_fit <- readRDS("rf_fit.rds")

# --- load examplar data CSV --- 
example_data <- read.csv("m6A_input_example.csv")

# --- define necessary functions for ML prediction --- 
dna_encoding <- function(dna_strings){
nn <- nchar( dna_strings[1] )
seq_m <- matrix( unlist( strsplit(dna_strings, "") ), ncol = nn, byrow = TRUE)
colnames(seq_m) <- paste0("nt_pos", 1:nn)
seq_df <- as.data.frame(seq_m)
seq_df[] <- lapply(seq_df, factor, levels = c("A", "T", "C", "G"))
return(seq_df)
}

prediction_multiple <- function(ml_fit, feature_df, positive_threshold = 0.5){
stopifnot(all(c("gc_content", "RNA_type", "RNA_region", "exon_length", "distance_to_junction", "evolutionary_conservation", "DNA_5mer") %in% colnames(feature_df))) #Check errors if incorrect column names of input data.frame
# 1. Set the factor levels exactly as used during model training
feature_df$RNA_type <- factor(feature_df$RNA_type, levels = c("mRNA", "lincRNA", "lncRNA", "pseudogene"))
feature_df$RNA_region <- factor(feature_df$RNA_region, levels = c("CDS", "intron", "3'UTR", "5'UTR"))

# 2. Build the model input: numeric features + categorical features + encoded nucleotide positions
model_input <- cbind(
  feature_df[, c("gc_content", "RNA_type", "RNA_region", "exon_length",
                 "distance_to_junction", "evolutionary_conservation")],
  dna_encoding(feature_df$DNA_5mer)
)

# 3. Predict the positive-class probability and apply the decision threshold
prob <- predict(ml_fit, newdata = model_input, type = "prob")[, "Positive"]

feature_df$predicted_m6A_prob <- prob
feature_df$predicted_m6A_status <- ifelse(prob > positive_threshold, "Positive", "Negative")

return(feature_df) #return a data.frame with supplied columns of predicted m6A prob and predicted m6A status
}

prediction_single <- function(ml_fit, gc_content, RNA_type, RNA_region, exon_length, distance_to_junction, evolutionary_conservation, DNA_5mer, positive_threshold = 0.5){
# 1. Build a one-row data frame with the same factor levels as in training
single_df <- data.frame(
  gc_content = gc_content,
  RNA_type = factor(RNA_type, levels = c("mRNA", "lincRNA", "lncRNA", "pseudogene")),
  RNA_region = factor(RNA_region, levels = c("CDS", "intron", "3'UTR", "5'UTR")),
  exon_length = exon_length,
  distance_to_junction = distance_to_junction,
  evolutionary_conservation = evolutionary_conservation,
  DNA_5mer = DNA_5mer,
  stringsAsFactors = FALSE
)

# 2. Reuse the multiple-row prediction function on this single row
result <- prediction_multiple(ml_fit, single_df, positive_threshold = positive_threshold)

# 3. Keep only the two predicted values
returned_vector <- c(predicted_m6A_prob = round(result$predicted_m6A_prob[1], 3),
                     predicted_m6A_status = result$predicted_m6A_status[1])

return(returned_vector) #return a named vector with values for predicted m6A prob and predicted m6A status
}

# --- create combined shiny code bellow ---
# ---- numbers shown on the home page ----
n_sites <- nrow(m6a_data)
n_positive <- sum(m6a_data$m6A_status == "Positive")
n_rna_types <- length(unique(m6a_data$RNA_type))

# ---- choices for the visualization filters ----
status_choices <- c("No filter", sort(unique(m6a_data$m6A_status)))
region_choices <- c("No filter", sort(unique(m6a_data$RNA_region)))
type_choices   <- c("No filter", sort(unique(m6a_data$RNA_type)))

ui <- dashboardPage(
  dashboardHeader(title = "m6A Database & Prediction Server", titleWidth = 360),
  dashboardSidebar(
    sidebarMenu(
      menuItem("Home", tabName = "home", icon = icon("home")),
      menuItem("Database", icon = icon("database"), startExpanded = TRUE,
        menuSubItem("Data Browser", tabName = "data", icon = icon("table")),
        menuSubItem("Visualization", tabName = "viz", icon = icon("bar-chart")),
        menuSubItem("Download", tabName = "download", icon = icon("download"))
      ),
      menuItem("Prediction Server", icon = icon("flask"), startExpanded = TRUE,
        menuSubItem("m6A Prediction Multiple", tabName = "multiple", icon = icon("upload")),
        menuSubItem("m6A Prediction Single", tabName = "single", icon = icon("edit"))
      )
    )
  ),
  dashboardBody(
    tags$head(
      tags$style(HTML("
        .content-wrapper p, .content-wrapper li { font-size: 14px; line-height: 1.65; }
        .content-wrapper h2 { font-size: 22px; font-weight: 600; margin-top: 4px; }
        .content-wrapper h4 { font-size: 15px; font-weight: 600; margin-top: 4px; margin-bottom: 8px; }
        .box-title { font-weight: 600; }
        .content-wrapper ul, .content-wrapper ol { padding-left: 20px; margin-bottom: 8px; }
        .content-wrapper li { margin-bottom: 6px; }
      "))
    ),
    tabItems(
      # ================= home page =================
      tabItem(
        tabName = "home",
        fluidRow(
          box(
            width = 12,
            imageOutput("home_image", width = "100%", height = "240px"),
            p("m6A (N6-methyladenosine) methylation of an RNA transcript, catalysed by the methyltransferase complex.",
              style = "color:#666; font-size:12px; margin-top:6px; margin-bottom:0;")
          )
        ),
        fluidRow(
          valueBox(format(n_sites, big.mark = ","), "m6A sites",
                   icon = icon("database"), width = 4, color = "blue"),
          valueBox(format(n_positive, big.mark = ","), "positive sites",
                   icon = icon("bullseye"), width = 4, color = "green"),
          valueBox(n_rna_types, "RNA types",
                   icon = icon("list"), width = 4, color = "purple")
        ),
        fluidRow(
          box(
            title = "a) Database Features", status = "primary", solidHeader = TRUE, width = 6,
            p("The database collects m6A sites and their sequence, annotation and conservation information in one table."),
            tags$ul(
              tags$li(strong("48,055 m6A sites"), " described by 16 variables."),
              tags$li(strong("11,248 positive"), " and ", strong("36,807 negative"),
                      " sites, marked by ", strong("m6A_status"), "."),
              tags$li(strong("4 RNA types:"), " mRNA (41,365), lncRNA (5,439), lincRNA (869), pseudogene (382)."),
              tags$li(strong("4 RNA regions:"), " 3'UTR (19,307), CDS (15,968), intron (9,971), 5'UTR (2,809).")
            )
          ),
          box(
            title = "b) Data Features", status = "info", solidHeader = TRUE, width = 6,
            tags$ul(
              tags$li(strong("Genomic annotation (5):"), " chromosome, start, end, width, strand."),
              tags$li(strong("m6A measurements (4):"), " m6A_status, m6A_coverage, total_coverage, m6A_level."),
              tags$li(strong("Sequence features (2):"), " gc_content, DNA_5mer."),
              tags$li(strong("Transcript annotation (4):"), " RNA_type, RNA_region, exon_length, distance_to_junction."),
              tags$li(strong("Evolutionary conservation (1):"), " evolutionary_conservation.")
            ),
            p(strong("m6A_status"), " is the label of a site (Positive or Negative) and ",
              strong("m6A_level"), " is its methylation level, a value between 0 and 1.")
          )
        ),
        fluidRow(
          box(
            title = "c) Prediction Server", status = "success", solidHeader = TRUE, width = 6,
            p("A random forest model trained on all 48,055 sites predicts whether a new site is m6A methylated."),
            tags$ul(
              tags$li(strong("7 input features:"), " gc_content, RNA_type, RNA_region, exon_length, distance_to_junction, evolutionary_conservation and DNA_5mer."),
              tags$li(strong("Multiple mode:"), " upload a CSV of many sites, choose the probability threshold, predict every row and download the result."),
              tags$li(strong("Single mode:"), " set the features of one site and read its probability and status instantly.")
            )
          ),
          box(
            title = "d) How to quick start?", status = "warning", solidHeader = TRUE, width = 6,
            tags$ol(
              tags$li(strong("Data Browser:"), " search or filter the 16-column table."),
              tags$li(strong("Visualization:"), " read the summary counts, then filter and compare the distributions."),
              tags$li(strong("Download:"), " export the complete dataset or a single RNA type subset."),
              tags$li(strong("Prediction Multiple:"), " click Example data, then Predict, then download the predictions."),
              tags$li(strong("Prediction Single:"), " adjust the seven features and watch the prediction update.")
            )
          )
        )
      ),
      # ================= data page =================
      tabItem(
        tabName = "data",
        h2("m6A Data Browser"),
        p("Use the Search box or the column filters of the table to look for specific sites. Drag the horizontal scroll bar to see all 16 columns."),
        box(
          title = "m6A sites (16 variables)", status = "primary",
          solidHeader = TRUE, width = 12,
          DTOutput("data_table")
        )
      ),
      # ================= visualization page =================
      tabItem(
        tabName = "viz",
        h2("Visual exploration"),
        fluidRow(
          box(
            title = "Summary Statistics", status = "info", solidHeader = TRUE, width = 6,
            verbatimTextOutput("count_status")
          ),
          box(
            title = "Summary Statistics", status = "info", solidHeader = TRUE, width = 6,
            verbatimTextOutput("count_type")
          )
        ),
        fluidRow(
          box(
            title = "Pie chart of a categorical column", status = "primary",
            solidHeader = TRUE, width = 6,
            selectInput("pie_var", "Categorical column:",
                        choices = c("m6A_status", "RNA_region", "RNA_type")),
            plotOutput("pie_plot")
          ),
          box(
            title = "Histogram of m6A level after filtering", status = "primary",
            solidHeader = TRUE, width = 6,
            selectInput("f_status", "m6A_status:", choices = status_choices),
            selectInput("f_region", "RNA_region:", choices = region_choices),
            selectInput("f_type", "RNA_type:", choices = type_choices),
            plotOutput("level_hist")
          )
        )
      ),
      # ================= download page =================
      tabItem(
        tabName = "download",
        h2("Download m6A datasets"),
        p("All files are comma-separated values (CSV) and can be opened in Excel or R. The first file contains the whole database; the other four files contain only the sites located on one RNA type."),
        fluidRow(
          box(
            title = "Complete database", status = "primary", solidHeader = TRUE, width = 12,
            p("All ", strong("48,055"), " m6A sites with all 16 columns and no filtering. Use this file when you need the complete dataset, for example to train a machine learning model."),
            downloadButton("dl_total", "Download complete dataset")
          )
        ),
        fluidRow(
          box(
            title = "mRNA subset", status = "success", solidHeader = TRUE, width = 3,
            p(strong("41,365"), " sites with ", tags$code("RNA_type = mRNA"), "."),
            downloadButton("dl_mrna", "Download mRNA")
          ),
          box(
            title = "lncRNA subset", status = "success", solidHeader = TRUE, width = 3,
            p(strong("5,439"), " sites with ", tags$code("RNA_type = lncRNA"), "."),
            downloadButton("dl_lnc", "Download lncRNA")
          ),
          box(
            title = "lincRNA subset", status = "success", solidHeader = TRUE, width = 3,
            p(strong("869"), " sites with ", tags$code("RNA_type = lincRNA"), "."),
            downloadButton("dl_linc", "Download lincRNA")
          ),
          box(
            title = "pseudogene subset", status = "success", solidHeader = TRUE, width = 3,
            p(strong("382"), " sites with ", tags$code("RNA_type = pseudogene"), "."),
            downloadButton("dl_pseudo", "Download pseudogene")
          )
        )
      ),
      # ================= multiple-row prediction page =================
      tabItem(
        tabName = "multiple",
        h2("m6A Prediction Multiple"),
        p("Upload a CSV file that contains the seven input feature columns, then run the model on every row."),
        box(
          title = "Step 1. Upload an input feature CSV file", status = "primary",
          solidHeader = TRUE, width = 12,
          fileInput("file", "Choose a CSV file:", accept = ".csv"),
          actionButton("example", "Example data")
        ),
        box(
          title = "Step 2. Select a decision threshold for the positive class probability",
          status = "primary", solidHeader = TRUE, width = 12,
          sliderInput("threshold", "Positive probability threshold:",
                      min = 0, max = 1, value = 0.5, step = 0.01)
        ),
        box(
          title = "Step 3. Click Predict to run inference", status = "primary",
          solidHeader = TRUE, width = 12,
          actionButton("predict", "Predict"),
          downloadButton("download", "Download predictions")
        ),
        box(
          title = "Input data and prediction results", status = "info",
          solidHeader = TRUE, width = 12,
          DTOutput("table")
        )
      ),
      # ================= single-observation prediction page =================
      tabItem(
        tabName = "single",
        h2("m6A Prediction Single"),
        p("Set the feature values below. The model will automatically predict whether this site is m6A ('Positive' or 'Negative'). The probability score reflects model confidence."),
        fluidRow(
          box(
            title = "Feature values", status = "primary", solidHeader = TRUE, width = 6,
            selectInput("s_rna_type", "RNA_type:",
                        choices = c("mRNA", "lincRNA", "lncRNA", "pseudogene"), selected = "mRNA"),
            selectInput("s_rna_region", "RNA_region:",
                        choices = c("CDS", "intron", "3'UTR", "5'UTR"), selected = "CDS"),
            sliderInput("s_gc", "gc_content:", min = 0, max = 1, value = 0.5, step = 0.01),
            sliderInput("s_exon", "exon_length:", min = 0, max = 15, value = 10),
            sliderInput("s_dist", "distance_to_junction:", min = 0, max = 15, value = 8),
            sliderInput("s_cons", "evolutionary_conservation:", min = 0, max = 1, value = 0.5, step = 0.01),
            textInput("s_5mer", "DNA_5mer (five nucleotides of A/T/C/G):", value = "GGACA")
          ),
          box(
            title = "Prediction result", status = "success", solidHeader = TRUE, width = 6,
            p("The result updates automatically whenever a feature value changes."),
            valueBoxOutput("prob_box", width = 12),
            valueBoxOutput("status_box", width = 12)
          )
        )
      )
    )
  )
)

server <- function(input, output, session) {
  # ---- home page illustration ----
  output$home_image <- renderImage({
    list(src = "m6A_home.jpg",
         contentType = "image/jpeg",
         width = "100%",
         height = 240,
         alt = "Illustration of m6A RNA methylation")
  }, deleteFile = FALSE)

  # ---- data table ----
  output$data_table <- renderDT({
    datatable(m6a_data, options = list(scrollX = TRUE, pageLength = 15))
  })

  # ---- visualization: pie chart ----
  output$pie_plot <- renderPlot({
    var <- input$pie_var
    counts <- as.data.frame(table(m6a_data[[var]]), stringsAsFactors = FALSE)
    colnames(counts) <- c("category", "count")
    counts$pct <- counts$count / sum(counts$count)

    ggplot(counts, aes(x = "", y = count, fill = category)) +
      geom_col(width = 1, colour = "white") +
      coord_polar(theta = "y") +
      geom_text(aes(label = paste0(category, "\n", round(pct * 100, 1), "%")),
                position = position_stack(vjust = 0.5), size = 4.5) +
      labs(title = paste("Distribution of", var), fill = var) +
      theme_void() +
      theme(plot.title = element_text(hjust = 0.5, face = "bold"))
  })

  # ---- visualization: reactive filtering ----
  filtered_data <- reactive({
    d <- m6a_data
    if (input$f_status != "No filter") d <- d[d$m6A_status == input$f_status, ]
    if (input$f_region != "No filter") d <- d[d$RNA_region == input$f_region, ]
    if (input$f_type   != "No filter") d <- d[d$RNA_type   == input$f_type, ]
    d
  })

  # ---- visualization: count of sites per m6A status ----
  output$count_status <- renderPrint({
    d <- filtered_data()
    cat("Counts:\n")
    for (s in sort(unique(m6a_data$m6A_status))) {
      cat(sprintf("%s : %s\n", s, sum(d$m6A_status == s)))
    }
  })

  # ---- visualization: count of sites per RNA type ----
  output$count_type <- renderPrint({
    d <- filtered_data()
    cat("Counts:\n")
    for (t in sort(unique(m6a_data$RNA_type))) {
      cat(sprintf("%s : %s\n", t, sum(d$RNA_type == t)))
    }
  })

  # ---- visualization: histogram of m6A level ----
  output$level_hist <- renderPlot({
    d <- filtered_data()
    ggplot(d, aes(x = m6A_level)) +
      geom_histogram(aes(y = after_stat(density)), bins = 30,
                     fill = "steelblue", colour = "white", alpha = 0.8) +
      geom_density(colour = "darkred", linewidth = 1) +
      labs(title = "m6A level distribution",
           subtitle = paste(nrow(d), "sites after filtering"),
           x = "m6A level", y = "Density") +
      theme_bw() +
      theme(plot.title = element_text(hjust = 0.5, face = "bold"),
            plot.subtitle = element_text(hjust = 0.5))
  })

  # ---- database downloads ----
  output$dl_total <- downloadHandler(
    filename = "m6A_total.csv",
    content = function(file) write.csv(m6a_data, file, row.names = FALSE)
  )

  output$dl_mrna <- downloadHandler(
    filename = "m6A_mRNA.csv",
    content = function(file) {
      write.csv(m6a_data[m6a_data$RNA_type == "mRNA", ], file, row.names = FALSE)
    }
  )

  output$dl_linc <- downloadHandler(
    filename = "m6A_lincRNA.csv",
    content = function(file) {
      write.csv(m6a_data[m6a_data$RNA_type == "lincRNA", ], file, row.names = FALSE)
    }
  )

  output$dl_lnc <- downloadHandler(
    filename = "m6A_lncRNA.csv",
    content = function(file) {
      write.csv(m6a_data[m6a_data$RNA_type == "lncRNA", ], file, row.names = FALSE)
    }
  )

  output$dl_pseudo <- downloadHandler(
    filename = "m6A_pseudogene.csv",
    content = function(file) {
      write.csv(m6a_data[m6a_data$RNA_type == "pseudogene", ], file, row.names = FALSE)
    }
  )

  # ---- multiple-row prediction ----
  input_data <- reactiveVal(NULL)
  prediction_results <- reactiveVal(NULL)

  observeEvent(input$file, {
    req(input$file)
    input_data(read.csv(input$file$datapath, stringsAsFactors = FALSE))
    prediction_results(NULL)
  })

  observeEvent(input$example, {
    input_data(example_data)
    prediction_results(NULL)
  })

  observeEvent(input$predict, {
    req(input_data())
    prediction_results(prediction_multiple(rf_fit, input_data(), positive_threshold = input$threshold))
  })

  observeEvent(input$threshold, {
    req(input_data(), prediction_results())
    prediction_results(prediction_multiple(rf_fit, input_data(), positive_threshold = input$threshold))
  }, ignoreInit = TRUE)

  output$table <- renderDT({
    displayed <- if (is.null(prediction_results())) input_data() else prediction_results()
    req(displayed)
    datatable(displayed, options = list(scrollX = TRUE, pageLength = 15))
  })

  output$download <- downloadHandler(
    filename = "predicted.csv",
    content = function(file) {
      req(input_data())
      result <- if (is.null(prediction_results())) input_data() else prediction_results()
      write.csv(result, file, row.names = FALSE)
    }
  )

  # ---- single-observation prediction ----
  single_result <- reactive({
    validate(need(nchar(input$s_5mer) == 5 && grepl("^[ATCG]{5}$", toupper(input$s_5mer)),
                  "DNA_5mer must be 5 characters, each one of A, T, C or G."))
    prediction_single(rf_fit,
                      gc_content = input$s_gc,
                      RNA_type = input$s_rna_type,
                      RNA_region = input$s_rna_region,
                      exon_length = input$s_exon,
                      distance_to_junction = input$s_dist,
                      evolutionary_conservation = input$s_cons,
                      DNA_5mer = toupper(input$s_5mer),
                      positive_threshold = 0.5)
  })

  output$prob_box <- renderValueBox({
    valueBox(round(as.numeric(single_result()["predicted_m6A_prob"]), 3),
             "Predicted m6A probability", icon = icon("percent"), color = "blue")
  })

  output$status_box <- renderValueBox({
    status <- single_result()["predicted_m6A_status"]
    valueBox(status, "Predicted m6A status",
             icon = icon(if (status == "Positive") "check" else "times"),
             color = if (status == "Positive") "green" else "red")
  })
}

shinyApp(ui, server)
