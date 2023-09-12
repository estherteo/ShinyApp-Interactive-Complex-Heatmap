# PROBLEM DEFINITION
# in the project details.docx - Generate Heat map of all markers (exclude DNA) for all patient ids, highlight indication(treated/untreated)

# APPROACH TO SOLVE IT
# Calculate the average expression value for each marker, for each patient.
# Then normalize the data across all patients per marker
# Generate three heatmaps (Original, Standard Normalization and Z-Score Normalisation) in Shiny App
# Add column annotation for patientID and indication (NACT or NAIVE)
# Let the user select the type of available heatmap and turn on or off the col annotation (patient ID and indication) through dropdown

install.packages("tidyr")

install.packages("viridis")  # Install
library("viridis")           # Load

install.packages("viridisLite")
library(viridisLite)

install.packages("colorRamp2")
library(colorRamp2)

install.packages("randomcoloR")
library(randomcoloR)

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("ComplexHeatmap")
library(ComplexHeatmap)
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("InteractiveComplexHeatmap")
library(InteractiveComplexHeatmap)

library("data.table")
library("tidyr")
library(tidyverse)
library(shiny)

setwd("~/Documents/R")
set.seed(2020)

my_data <- fread('cells.csv')
View(my_data)

# take out the last 7 cols from the dataset
my_data <- my_data[, 1:(ncol(my_data) - 7)]

View(my_data)
my_data$Image <- gsub("\\.tiff$", "", my_data$Image)

# avg out the values for each marker for each patientID 
my_data <- my_data[, lapply(.SD, mean, na.rm=TRUE), by=Image ]
my_data[, "Object" := NULL]
View(my_data)
colnames(my_data)[colnames(my_data) == "Image"] <- "PatientID"

# Data Manipulation
# convert to long data format 
# Three cols: PatientID, Marker & Value(avg)
long <- melt(setDT(my_data), id.vars = c("PatientID"), variable.name = "Marker", value.name = "Gene Expression")    

View(long)    

# try to do dcast to convert it to wide data format for the matrix
# By converting it to wide, transpose the rows and cols such that Marker is now the rows and PatientID is now the cols
data_for_matrix <- dcast(long, Marker ~ PatientID, value.var = "Gene Expression")

View(data_for_matrix)
summary(data_for_matrix) # displays 57 columns
ncol(data_for_matrix) # 57
data_for_matrix[,57]

# up until dcast, still have all the cols intact

# aft dcast
# make the first column as the rows of that matrix and the remaining columns as the columns of the matrix
rnames <- data_for_matrix[,1]              
mat_data <- data.matrix(data_for_matrix[,2:57])  
rownames(mat_data) <- unlist( rnames )               

View(mat_data)

is.matrix(mat_data)

ncol(mat_data)  #56 cols 

# Get the first 6 characters from column names and create unique groups (14 of them)
grouped_patient_ids_heatmap <- unique(substr(colnames(mat_data), 1, 6))

# Group grouped_patient_ids into two groups 
patient_groups_heatmap <- ifelse(grouped_patient_ids_heatmap %in% c("TAR006", "TAR018", "TAR020", "TAR027", "TAR028", "TAR050"),
                         "RESPONDER",
                         "NON-RESPONDER")

# Assign orange to RESPONDER and blue to NON-RESPONDER
group_colors_indication <- c("RESPONDER" = "orange", "NON-RESPONDER"="lightblue")

# Create a vector of group labels for each column
group_labels_indication_heatmap <- rep(patient_groups_heatmap, each = 4) # Repeat each group label 4 times

# Create a vector of group labels for each column
group_labels_heatmap <- sapply(colnames(mat_data), function(col_name) {
  substr(col_name, 1, 6)
})

group_labels_heatmap

n_groups_heatmap <- length(grouped_patient_ids_heatmap)

# Generate the color palette for patient_id annotation
color_palette_heatmap <- randomColor(n_groups_heatmap,luminosity="dark")

# Create a named vector of colors for each unique group
group_colors_heatmap <- setNames(color_palette_heatmap, grouped_patient_ids_heatmap)

# Map group_colors to each group_label  
top_annotation_colors_heatmap <- group_colors_heatmap[group_labels_heatmap]

# Added in this part --
# Create a vector of group labels for each column
patientid_annotation_heatmap <- data.frame(patient_id = group_labels_heatmap)

indication_annotation_heatmap <- data.frame(indication = group_labels_indication_heatmap)

colors_patientid_annotation <- list(patient_id = top_annotation_colors_heatmap)
colors_indication_annotation <- list(indication = group_colors_indication)

patientid_ha_heatmap <- HeatmapAnnotation(
  df = patientid_annotation_heatmap,
  which = 'col',
  col = colors_patientid_annotation, 
  show_legend = TRUE,
  annotation_width = unit(c(1, 4), 'cm'),
  gap = unit(1, 'mm')
)

indication_ha_heatmap <- HeatmapAnnotation(
  df = indication_annotation_heatmap,
  which = 'col',
  col = colors_indication_annotation, 
  show_legend = TRUE,
  annotation_width = unit(c(1, 4), 'cm'),
  gap = unit(1, 'mm')
)

# Added in above part --

combi_annotation_heatmap <- data.frame(patient_id = group_labels_heatmap,indication = group_labels_indication_heatmap)
colors_heatmap <- list(patient_id = top_annotation_colors_heatmap, indication = group_colors_indication)

column_ha_heatmap <- HeatmapAnnotation(
  df = combi_annotation_heatmap,
  which = 'col',
  col = colors_heatmap, 
  show_legend = TRUE,
  annotation_width = unit(c(1, 4), 'cm'),
  gap = unit(1, 'mm')
)

# FEATURE SCALING
# 1. Z-score Normalisation -subtract the mean from the indiv values and divide by sd (z-miu / sigma)
# mean is 0 and sd is 1
transposed_data <- t(mat_data)  # Transpose the data

standard_normalized_data <- scale(transposed_data)

View(standard_normalized_data)

# need to transpose it again so the format matches 
transposed_standard_normalised <- t(standard_normalized_data)

is.matrix(transposed_standard_normalised)

View(transposed_standard_normalised)

# 2. Standard Normalisation 
# rescales a dataset so each value falls in btw 0 and 1 

normalize <- function(x, na.rm = TRUE) {
  return((x- min(x)) /(max(x)-min(x)))
}

# Apply min-max scaling to each column of the matrix
scaled_data <- apply(mat_data, 2, normalize)

View(scaled_data)

is.matrix(scaled_data)

# SHINY APP
# Let user select to turn on or off the col annotation (patient ID and indication)
ui = fluidPage(
  # App title ----
  titlePanel("Shiny App for Heatmap Visualisation of Spatial Omics Data"),
  sidebarLayout(
    
    # Sidebar panel for inputs ----
    sidebarPanel(
      fileInput('file', 'Select a file',
                accept=c('text/csv','text/comma-separated-values,text/plain','.csv')),
      selectInput("heatmap_type", "Select type of Heatmap", 
                  choices = c("Original", "Z-Score Normalisation", "MinMax Normalisation")),
      selectInput("col_annotation", "Select the column annotation", choices=c("None", "Patient ID", "Indication", "Patient ID and Indication")),
    ),
    mainPanel(
      plotOutput("heatmap") 
    )  
  )
)

server <- function(input,output, session){
  options(shiny.maxRequestSize = 100*1024^2)
  sessioninfo::session_info()
  writeLines(capture.output(sessionInfo()), "sessionInfo.txt")
  observe({
    file = input$file
    if (is.null(file)){
      return(NULL)
    }
    data = read.csv(file$datapath)
    
  })
    
  datasetInput <- reactive({
    switch(input$heatmap_type,
           "Original" = mat_data,
           "Z-Score Normalisation" = transposed_standard_normalised,
           "MinMax Normalisation" = scaled_data)
  })
  
  colAnnotation <- reactive({
    switch(input$col_annotation,
           "None" = NULL,
           "Patient ID" = patientid_ha_heatmap,
           "Indication" = indication_ha_heatmap,
           "Patient ID and Indication" = column_ha_heatmap
          )
  })
  
  datasetTitle <- reactive({
    switch(input$heatmap_type,
           "Original" = "Single−cell marker expression per image",
           "Z-Score Normalisation" = "Single−cell marker expression per image - Z-Score Normalisation",
           "MinMax Normalisation" = "Single−cell marker expression per image - MinMax Normalisation")
  })
  
  output$heatmap <- renderPlot({
    heatmap <- Heatmap(datasetInput(), name = "Gene Expression Values",column_title = datasetTitle(),
                                 row_names_gp = gpar(fontsize = 8),
                                 column_names_gp = gpar(fontsize = 8),
                                 col = viridis(25),
                                 top_annotation=colAnnotation())
    heatmap <- draw(heatmap)
    
  }, res = 96)
}

shiny::shinyApp(ui, server)


