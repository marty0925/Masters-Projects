#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    https://shiny.posit.co/
#

library(shiny)
library(ggplot2)
library(plotly)
library(dplyr)
library(DT)
library(openxlsx)


#--------Setup section--------------
#Loading RNA Data

input_data <- read.xlsx("C:/Users/invate/OneDrive - University of Southampton/Desktop/KCL Applied Bio 25-26/7BBG1002 Applied Bio and Cloud Computing/Group project/Coursework Replication/Georgia Work/DESeq2_control_vs_evolved_padj0.01_genes.xlsx",
                        sep = ",", 
                        na.strings = c("NA", "NULL", "null", ""))
clean_df <- na.omit(input_data) # Removes empty rows
clean_df$Gene_Name <- ifelse(is.na(clean_df$Gene_Name) | clean_df$Gene_Name == "", clean_df$Yeast_id, clean_df$Gene_Name)
head(input_data)

#Loading Proteomics data
prot_raw <- read.csv("C:/Users/invate/OneDrive - University of Southampton/Desktop/KCL Applied Bio 25-26/7BBG1002 Applied Bio and Cloud Computing/Group project/Coursework Replication/Proteomics Process/Original Attempt 2 Results.csv", header = TRUE, stringsAsFactors = FALSE)
prot_df <- prot_raw[-c(1,2), ] # Remove first two no needed rows

prot_df$Diff <- as.numeric(prot_df$Student.s.T.test.Difference.Control_Reconstruction) # Making columns numeric
prot_df$NegLogP <- as.numeric(prot_df$X.Log.Student.s.T.test.p.value.Control_Reconstruction)
prot_df$Protein_ID <- prot_df$Majority.protein.IDs

prot_df <- na.omit(prot_df[, c("Protein_ID", "Diff", "NegLogP")]) # Removing the rows that didn't parse correctly

#-----UI side-----
# Define UI for application that draws a histogram
ui <- tagList(
  tags$head(
    tags$style(HTML("
      /* OUtline for tab buttoms */
      .navbar-nav > li > a {
          border: 2px solid #1465AC; /* Blue Outline */
          border-radius: 5px;                   /* Rounded Corners */
          margin-right: 10px;                   /* Space between tabs */
          font-weight: bold;                    /* Bold Text */
          background-color: white;   /* White background */
          color: black;              /* Black Text */
      }

      /* 2. style of activly selected tab */
      .navbar-nav > .active > a {
          background-color: #1465AC; /* Blue background when selected */
          color: white;              /* White text when selected */
          border: 2px solid #0f4c81; /* Darker blue border */
      }
    "))
  ),

  navbarPage("Multi-omics Dashbaord",
 tabPanel("RNA-Seq Analysis",
  sidebarLayout(
    sidebarPanel(
      h4("RNA Parameters"),
      sliderInput("pval_cutoff", "Adjusted P-Value Cutoff",
                  min = 0, max = 0.1, value = 0.01, step = 0.01),
      
      sliderInput("lfc_cutoff", "Log2 Fold Change Cutoff:", 
                  min = 0, max = 3, value = 1, step = 0.5)
    ),
    
    mainPanel(
      h3("Gene Expression Exploration"), #Header for clarity
      
      tabsetPanel(
        # First sub=Tab
        tabPanel("Volcano Explorer",
                 br(),
                 plotlyOutput("volcanoPlot"),
                 br(), hr(),
                 DT::DTOutput("resultsTable")
        ),
        
        # Second tab
        tabPanel("Highlighted Genes by Expression", 
                 br(),
                 selectInput("bar_direction", "Show Genes:", 
                             choices = c("Most Upregulated", "Most Downregulated")),
                 
                 sliderInput("bar_count", "Number of Genes:", 
                             min = 5, max = 20, value = 10),
                 
                 plotlyOutput("barPlot")
        )
      )
    )
  )
),
tabPanel("Proteomics Analysis",
         sidebarLayout(
           sidebarPanel(
             h4("Proteomics Parameters"),
             p("Note: Protein Id is in the SGD systematic name."),
             
             # Independent sliders for Proteomics
             sliderInput("prot_pval", "-Log10 P-value Cutoff:", min = 0, max = 6, value = 0.05, step = 0.1),
             sliderInput("prot_fc", "Log2 Difference Cutoff:", min = 0, max = 5, value = 1, step = 0.1),
             hr(),
             helpText("Positive Difference = Enriched in Control (Downregulated).")
           ),
           
           mainPanel(
             tabsetPanel(
               # Proteomics Sub-Tab 1: Volcano
               tabPanel("Proteins Volcano Plot",
                        br(),
                        plotlyOutput("protVolcano"),
                        br(), hr(),
                        h4("Detailed Proteomics Results"),
                        DT::DTOutput("protTable")
            )
           )
          )
         )
        )
)
)


#-------Server side-----------
# Defining the server for the interactive figure
server <- function(input, output) {
  
  
#---------------------------------RNA Side  
  # Making the sliders reactive
  filtered_data <- reactive({
    data <- clean_df  # Uses the global dataframe loaded at the top
    
    data %>%
      mutate(
        # Create a NEW 'Significance' column based on the SLIDER inputs
        Significance = case_when(
          padj < input$pval_cutoff & abs(log2FoldChange) > input$lfc_cutoff ~ "Significant",
          TRUE ~ "Not Significant"
        ),
        # Distinguish Up vs Down for coloring
        ColorGroup = case_when(
          padj < input$pval_cutoff & log2FoldChange > input$lfc_cutoff ~ "Upregulated",
          padj < input$pval_cutoff & log2FoldChange < -input$lfc_cutoff ~ "Downregulated",
          TRUE ~ "NS"
        ),
        # Recalculate -log10(p) for the plot y-axis
        negLogP = -log10(padj)
      )
    })
  
  #adding in the tables for better clarity
  output$resultsTable <- DT::renderDT({
    
    data <- filtered_data() # for reactive data tables
    
    clean_table <- data %>% # For selecting specific tables
      filter(Significance == "Significant") %>%
      select(Gene_Name, log2FoldChange, padj, Significance) %>%
      arrange(padj)
    
    # creates the interactive table
    DT::datatable(clean_table,
                  rownames = FALSE,
                  options = list(pageLength = 10, scrollX = TRUE),
                  colnames = c("Gene Name", "Log2 Fold Change", "Adj. P-value", "Status")) %>%
      # numbers are easier to read
      DT::formatRound(columns = c("log2FoldChange"), digits = 2) %>%
      DT::formatSignif(columns = c("padj"), digits = 3)
  })
  
  output$barPlot <- renderPlotly({
    
#Getring the data
    data <- filtered_data() %>%
      filter(Significance == "Significant")
    
  #Sorting the data
    if (input$bar_direction == "Most Upregulated") {
      top_data <- data %>%
        arrange(desc(log2FoldChange)) %>% # Highest positive values
        head(input$bar_count)
    } else {
      top_data <- data %>%
        arrange(log2FoldChange) %>%       # Lowest negative values
        head(input$bar_count)
    }

    
#creating the bar chart
    p <- ggplot(top_data, aes(x = reorder(Gene_Name, log2FoldChange), 
                              y = log2FoldChange, 
                              fill = log2FoldChange)) +
      geom_bar(stat = "identity") +
      coord_flip() + #FLips to make it easier to read
      
      scale_fill_gradient2(low = "#1465AC", mid = "pink", high = "#B31B21", midpoint = 0) +
      
      labs(title = paste("Top", input$bar_count, input$bar_direction),
           x = "Gene Name", 
           y = "Log2 Fold Change") +
      theme_minimal(base_size = 14) +
      theme(legend.position = "none")
    
    ggplotly(p, tooltip = c("x", "y"))
})

  
#creating the plot
  output$volcanoPlot <- renderPlotly({
    
    plot_data <- filtered_data()
    
    # Createing the base ggplot
    p <- ggplot(plot_data, aes(x = log2FoldChange, 
                               y = negLogP, 
                               color = ColorGroup,
                               # 'text' aesthetic defines what you see when hovering
                               text = paste("Gene:", Gene_Name, 
                                            "<br>Log2FC:", round(log2FoldChange, 2), 
                                            "<br>Padj:", formatC(padj, format = "e", digits = 2)))) +
      
      geom_point(alpha = 0.8, size = 1.5) +
      
      geom_vline(xintercept = c(-input$lfc_cutoff, input$lfc_cutoff), 
                 linetype = "dashed", color = "black") +
      
      geom_hline(yintercept = -log10(input$pval_cutoff), 
                 linetype = "dashed", color = "black") +
      
      scale_color_manual(values = c("Upregulated" = "#F4B142", 
                                    "Downregulated" = "#89183C", 
                                    "NS" = "gray70")) +
      
      theme_minimal(base_size = 14) +
      
  
      labs(title = "Volcano Plot: Evolved vs Control", 
           x = "Log2 Fold Change", 
           y = "-Log10 Adjusted P-value") +
      
      theme(legend.position = "none")
    
    #Converts ggplot to interactive Plotly
    ggplotly(p, tooltip = "text")
  })
    
    
#-----------------------------------------------DNA Side

    prot_filtered <- reactive({
      prot_df %>%
        mutate(
          # Determine status based on proteomics sliders
          ColorGroup = case_when(
            NegLogP > input$prot_pval & Diff > input$prot_fc ~ "Enriched in Control",
            NegLogP > input$prot_pval & Diff < -input$prot_fc ~ "Enriched in Evolved",
            TRUE ~ "NS"
          )
        )
    })
    
    # Proteomics Volcano plot 
    output$protVolcano <- renderPlotly({
      data <- prot_filtered()
      p <- ggplot(data, aes(x = Diff, y = NegLogP, color = ColorGroup,
                            text = paste("Protein:", Protein_ID, "<br>Diff:", round(Diff, 2)))) +
        geom_point(alpha = 0.7, size = 1.5) +
        scale_color_manual(values = c("Enriched in Control"="#1465AC", "Enriched in Evolved"="#B31B21", "NS"="gray80")) +
        geom_vline(xintercept = c(-input$prot_fc, input$prot_fc), linetype="dashed") +
        geom_hline(yintercept = input$prot_pval, linetype="dashed") +
        theme_minimal() + theme(legend.position = "none") +
        labs(title = "Volcano Plot: Proteomics", x = "Log2 Fold Change (Control - Evolved)", y = "-Log10 P-value")
      ggplotly(p, tooltip = "text")
    })
    
    # Proteomics Table
    output$protTable <- renderDT({
      prot_filtered() %>% filter(ColorGroup != "NS") %>%
        select(Protein_ID, Diff, NegLogP, ColorGroup) %>%
        datatable(rownames=F) %>% formatRound(c("Diff", "NegLogP"), 2)
    
  })
}
# Run the application 
shinyApp(ui = ui, server = server)