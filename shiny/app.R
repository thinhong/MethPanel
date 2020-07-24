library(shiny)
library(shinyauthr)
library(shinyjs)
library(plotly)
library(gplots)
library(ggplot2)
library(reshape2)
library(processx)
library(withr)

multi.ttest <- function(mat, ...) {
  mat <- as.matrix(mat)
  n <- ncol(mat)
  p.mat<- matrix(NA, n, n)
  diag(p.mat) <- 1
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      test <- t.test(mat[, i], mat[, j], ...)
      p.mat[i, j] <- p.mat[j, i] <- test$p.value
    }
  }
  colnames(p.mat) <- rownames(p.mat) <- colnames(mat)
  signif(p.mat,4)
}


# Default path
p = "/home/thinh/bioinfo/proj_methpanel/test"

# Dataframe that holds usernames, passwords and other user data
user_base <- data.frame(
  user = c("user1", "user2", "user3"),
  password = c("pass1", "pass2", "pass3"), 
  permissions = c("admin", "standard", "standard"),
  name = c("User One", "User Two", "User Three"),
  stringsAsFactors = FALSE,
  row.names = NULL
)

ui <- navbarPage(title = "MethPanel",
  tabPanel(title = "Login",
    # must turn shinyjs on
    shinyjs::useShinyjs(),
    # add logout button UI 
    div(class = "pull-right", shinyauthr::logoutUI(id = "logout")),
    # add login panel UI function
    shinyauthr::loginUI(id = "login"),
    # setup table output to show user info after login
    strong(textOutput(outputId = "welcome")),
    uiOutput(outputId = "project")
  ),
  ##### FASTQ QC #####
  tabPanel(title = "1. Trimmed FASTQ QC",
    sidebarLayout(
      sidebarPanel(
        uiOutput(outputId = "qc_sample"),
        uiOutput(outputId = "download_button"),
        p(),
        uiOutput(outputId = "qc_dl_all_guide"),
        p(),
        uiOutput(outputId = "qc_dl_all_button")
      ),
      mainPanel(
        dataTableOutput(outputId = "qc_table") 
      )
    )
  ),
  ##### Alignment #####
  tabPanel(title = "2. Alignment",
    tabsetPanel(
      tabPanel(title = "Overview",
        sidebarLayout(
          sidebarPanel(
            uiOutput(outputId = "align_var"),
            uiOutput(outputId = "align_guide"),
            p(),
            uiOutput(outputId = "button_metrics")
          ),
          mainPanel(
            plotlyOutput(outputId = "alignment_violin")
          )
        )
      ),
      tabPanel(title = "Metrics for group",
        sidebarLayout(
          sidebarPanel(
            uiOutput(outputId = "align_guide2"),
            uiOutput(outputId = "align_group_button"),
            p(),
            uiOutput(outputId = "align_upload"),
            uiOutput(outputId = "align_var_group")
          ),
          mainPanel(
            plotlyOutput(outputId = "align_group_plot")
          )
        )
      )
    )
  ),
  ##### DNA methylation #####
  tabPanel(title = "3. DNA methylation",
    tabsetPanel(
      tabPanel(title = "Overview",
        sidebarLayout(
          sidebarPanel(
            uiOutput(outputId = "overview_select"),
            conditionalPanel(
              condition = "input.overview_select == 1| input.overview_select == 2",
              radioButtons(inputId = "heatmap_table",
                           label = "Select an option to display",
                           choices = c("Heatmap" = 1, "Table" = 2, "MDS plot" = 3))
            ),
            conditionalPanel(
              condition = "input.heatmap_table == 1| input.heatmap_table == 3",
              downloadButton(outputId = "downloadDNAMetPlot",
                             label = "Download plot as PDF file")
            )
          ),
          mainPanel(
            conditionalPanel(
              condition = "input.heatmap_table == 1",
              plotlyOutput(outputId = "overview_heatmap")
            ),
            conditionalPanel(
              condition = "input.heatmap_table == 2",
              dataTableOutput(outputId = "align_table")
            ),
            conditionalPanel(
              condition = "input.heatmap_table == 3",
            plotlyOutput(outputId = "mds_overview")
            )
          )
        )
      ),
      tabPanel(title = "By sample",
        sidebarLayout(
          sidebarPanel(
            uiOutput(outputId = "meth_sample_select")
          ),
          mainPanel(
            plotlyOutput(outputId = "meth_sample_plot")
          )
        )
      ),
      tabPanel(title = "By amplicon",
       sidebarLayout(
         sidebarPanel(
           uiOutput(outputId = "bigtable_guide1"),
           p(),
           uiOutput(outputId = "bigtable_dlbutton"),
           uiOutput(outputId = "single_all"),
           uiOutput(outputId = "Amplicon_methylation")
         ),
         mainPanel(
           plotlyOutput("Methylation_violin")
         )
       )
      ),
      tabPanel(title = "By group",
        sidebarLayout(
          sidebarPanel(
            uiOutput(outputId = "meth_guide"),
            uiOutput(outputId = "download_template"),
            p(),
            uiOutput(outputId = "file1"),
            uiOutput(outputId = "meth_group_mode"),
            conditionalPanel(
              condition = "input.meth_group_mode == 2",
              downloadButton(outputId = "pdf_mds_group",
                             label = "Download plot as PDF")
            ),
            conditionalPanel(
              condition = "input.meth_group_mode == 1",
              radioButtons(inputId = "meth_group_box",
                           label = "Select an option",
                           choices = c("By amplicon" = 1, 
                                       "DNA methylation in sample groups" = 2,
                                       "All amplicons" = 3))
            ),
            conditionalPanel(
              condition = "input.meth_group_mode == 1 && input.meth_group_box == 1",
              uiOutput(outputId = "meth_amplicon")
            )
          ),
          mainPanel(
            conditionalPanel(
              condition = "input.meth_group_mode == 1 && input.meth_group_box == 1",
              plotlyOutput(outputId = "Amplicon_group_violin"),
              verbatimTextOutput(outputId = "meth_group_ttest")
            ),
            conditionalPanel(
              condition = "input.meth_group_mode == 1 && input.meth_group_box == 2",
              plotlyOutput(outputId = "group_violin"),
              verbatimTextOutput(outputId = "meth_group_ttest2")
            ),
            conditionalPanel(
              condition = "input.meth_group_mode == 1 && input.meth_group_box == 3",
              plotlyOutput(outputId = "Amplicon_all")
            ),
            conditionalPanel(
              condition = "input.meth_group_mode == 2",
              plotlyOutput(outputId = "mds_group")
            )
          )
        )
      )
    )
  ),
  ##### Pattern #####
  tabPanel(title = "4. Pattern",
    sidebarLayout(
      sidebarPanel(
        uiOutput(outputId = "selectproj"),
        uiOutput(outputId = "result1"),
        downloadButton(outputId = "downloadPatternHeatmap",
                       label = "Download plot as PDF")
      ),
      mainPanel(
        plotlyOutput(outputId = "Heatmap")
      )
    )
  ),
  ##### Polymorphism #####
  tabPanel(title = "5. Polymorphism",
    sidebarLayout(
      sidebarPanel(
        uiOutput(outputId = "poly_single_all"),
        uiOutput(outputId = "polymorphism_amp"),
        downloadButton(outputId = "pdf_polymorphism",
                       label = "Download plot as PDF")
      ),
      mainPanel(
        plotlyOutput("Polymorphism_scatter_plot"),
        DT::dataTableOutput(outputId = "Polymorphism_table")
      )
    )
  ),
  ##### PCR bias correction #####
  tabPanel(title = "6. PCR bias correction",
    tabsetPanel(
      tabPanel(title = "1. For amplicons",
        sidebarLayout(
          sidebarPanel(
            uiOutput(outputId = "pcr_guide"),
            uiOutput(outputId = "download_pcr"),
            p(),
            uiOutput(outputId = "upload_expected"),
            p(),
            uiOutput(outputId = "pcr_guide_3"),
            p(),
            uiOutput(outputId = "pcr_amplicon"),
            p(),
            uiOutput(outputId = "download_corrected_amps_guide"),
            p(),
            uiOutput(outputId = "download_corrected_amps")
          ),
          mainPanel(
            plotlyOutput(outputId = "plot_b_uncorrected"),
            p(),
            fluidRow(
              column(6, plotlyOutput(outputId = "spikesin_before")),
              column(6, plotlyOutput(outputId = "spikesin_after"))
            ),
            p(),
            plotlyOutput(outputId = "beforeafter_plot")
          )
        )
      ),
      tabPanel(title = "2. For CpG sites",
        sidebarLayout(
          sidebarPanel(
            uiOutput(outputId = "pcr_guide_2"),
            p(),
            uiOutput(outputId = "pcr_amplicon_CpGs"),
            p(),
            uiOutput(outputId = "download_corrected_CpGs_guide"),
            p(),
            uiOutput(outputId = "download_corrected_CpGs")
          ),
          mainPanel(
            plotlyOutput(outputId = "plot_b_uncorrected_CpGs"),
            p(),
            fluidRow(
              column(6, plotlyOutput(outputId = "spikesin_before_CpGs")),
              column(6, plotlyOutput(outputId = "spikesin_after_CpGs"))
            ),
            p(),
            plotlyOutput(outputId = "beforeafter_plot_CpGs")
          )
        )
      )
    )
  ),
  ##### DNA methylation after correction #####
  tabPanel(title = "7. DNA methylation after correction",
    tabsetPanel(
      tabPanel(title = "By sample",
        sidebarLayout(
          sidebarPanel(
            uiOutput(outputId = "after_sample_select")
          ),
          mainPanel(
            plotlyOutput(outputId = "after_sample_plot")
          )
        )
      ),
      tabPanel(title = "By CpG site",
        sidebarLayout(
          sidebarPanel(
            uiOutput(outputId = "after_single_all"),
            uiOutput(outputId = "after_CpG_select")
          ),
          mainPanel(
            plotlyOutput(outputId = "after_CpG_plot")
          )
        )
      ),
      tabPanel(title = "By group",
        sidebarLayout(
          sidebarPanel(
            uiOutput(outputId = "after_group_mode"),
            conditionalPanel(
              condition = "input.after_group_mode == 1",
              radioButtons(inputId = "after_group_box",
                           label = "Select an option",
                           choices = c("By amplicon" = 1, 
                                       "DNA methylation in sample groups" = 2,
                                       "All amplicons" = 3))
            ),
            conditionalPanel(
              condition = "input.after_group_mode == 1 && input.after_group_box == 1",
              uiOutput(outputId = "after_group_select")
            )
          ),
          mainPanel(
            conditionalPanel(
              condition = "input.after_group_mode == 1 && input.after_group_box == 1",
              plotlyOutput(outputId = "after_group_single"),
              verbatimTextOutput(outputId = "after_group_ttest")
            ),
            conditionalPanel(
              condition = "input.after_group_mode == 1 && input.after_group_box == 2",
              plotlyOutput(outputId = "after_group_sum"),
              verbatimTextOutput(outputId = "after_group_ttest2")
            ),
            conditionalPanel(
              condition = "input.after_group_mode == 1 && input.after_group_box == 3",
              plotlyOutput(outputId = "after_group_all")
            ),
            conditionalPanel(
              condition = "input.after_group_mode == 2",
              plotlyOutput(outputId = "mds_group_after")
            )
          )
        )
      )
    )
  )
)

server <- function(input, output, session) {
  ##### Login #####
  logout_init <- callModule(shinyauthr::logout, 
                            id = "logout", 
                            active = reactive(credentials()$user_auth))
  credentials <- callModule(shinyauthr::login, 
                            id = "login", 
                            data = user_base,
                            user_col = user,
                            pwd_col = password,
                            log_out = reactive(logout_init()))
  output$welcome <- renderText({
    req(credentials()$user_auth)
    paste("Welcome", credentials()$info$name, sep = " ")
  })
  output$project <- renderUI({
    req(credentials()$user_auth)
    selectInput(inputId = "project",
                label = "Please select a project",
                choices = list.dirs(paste(p, credentials()$info$user, sep = "/"),
                                    full.names = F, recursive = F))
  })
  ##### Load data #####
  bigTable_filtered <- reactive({
    req(credentials()$user_auth)
    read.table(paste(p, credentials()$info$user, input$project, "bigTable/bigTable_filtered.tsv", sep = "/"), 
               header = T, stringsAsFactors = F, check.names = F, na.strings = ".")
  })
  ##### FASTQ QC #####
  output$qc_sample <- renderUI({
    req(credentials()$user_auth)
    selectInput(inputId = "qc_sample",
                label = "Select a sample to download the trimmed FASTQ QC report",
                choices = list.dirs(paste(p, credentials()$info$user, input$project, "raw_trimmed", sep = "/"),
                                    full.names = F, recursive = F))
  })
  output$download_button <- renderUI({
    req(credentials()$user_auth)
    downloadButton(outputId = "qc_download",
                   label = "Download")
  })
  output$qc_download <- downloadHandler(
    filename = function() {
      paste0(input$qc_sample, "_QC.html")
    },
    content = function(con) {
      file.copy(paste(p, credentials()$info$user, input$project, "raw_trimmed", input$qc_sample, 
                              paste0(input$qc_sample, "_R1_trimmed_fastqc.html"), sep = "/"), con)
    }
  )
  output$qc_dl_all_guide <- renderUI({
    req(credentials()$user_auth)
    strong("Click here to download FASTQ QC reports of all samples")
  })
  output$qc_dl_all_button <- renderUI({
    req(credentials()$user_auth)
    downloadButton(outputId = "qc_dl_all",
                   label = "Download all FASTQ QC reports")
  })
  output$qc_dl_all <- downloadHandler(
    filename = "fastqc_html.zip",
    content = function(con) {
      file.copy(paste(p, credentials()$info$user, input$project, "bigTable/metrics/fastqc_html.zip", sep = "/"), con)
    }
  )
  df_qc <- reactive({
    req(credentials()$user_auth)
    read.table(paste(p, credentials()$info$user, input$project, "bigTable/metrics/trimmed_matrix.tsv", sep = "/"),
               header = T, stringsAsFactors = F, check.names = F, sep="\t", comment.char = "\\")
  })
  output$qc_table <- renderDataTable({
    df_qc()
  })
  ##### Template data generate #####
  # Generate template data
  group_template <- reactive({
    df_temp <- as.data.frame(df_qc()$Sample)
    colnames(df_temp)[1] <- "Sample"
    df_temp$Group <- ""
    df_temp
  })
  ##### Alignment #####
  ##### Overview #####
  df_align_download <- reactive({
    req(credentials()$user_auth)
    read.table(file.path(p, credentials()$info$user, input$project, "bigTable/metrics/metrics_table.tsv"), 
               header = T, stringsAsFactors = F, check.names = F)
  })
  df_alignment <- reactive({
    req(df_align_download())
    df_alignment <- df_align_download()
    df_alignment$`Number of CpG dropout` <- as.numeric(sapply(df_alignment$`No. CpG dropout`, function(x) strsplit(x, "/")[[1]][[1]]))
    df_alignment$`Number of Amplicon dropout` <- as.numeric(sapply(df_alignment$`No. Amplicon dropout`, function(x) strsplit(x, "/")[[1]][[1]]))
    df_alignment$`Number of CpG dropout with full read` <- as.numeric(sapply(df_alignment$`No. CpG dropout with full read`, function(x) strsplit(x, "/")[[1]][[1]]))
    df_alignment$`Number of Amplicon dropout with full read` <- as.numeric(sapply(df_alignment$`No. Amplicon dropout with full read`, function(x) strsplit(x, "/")[[1]][[1]]))
    df_alignment
  })
  output$align_var <- renderUI({
    selectInput(inputId = "align_var",
                label = "Select a metric to plot",
                choices = colnames(Filter(is.numeric, df_alignment())))
  })
  output$align_guide <- renderUI({
    req(credentials()$user_auth)
    strong("Click here to download full metrics table")
  })
  output$button_metrics <- renderUI({
    req(credentials()$user_auth)
    downloadButton(outputId = "download_metrics",
                   label = "Download metrics table")
  })
  output$download_metrics <- downloadHandler(
    filename = "metrics_table.tsv",
    content = function(con) {
      write.table(df_align_download(), sep ="\t", row.names = F, con)
    }
  )
  output$alignment_violin <- renderPlotly({
    req(input$align_var)
    df <- reshape(df_alignment(), idvar = "Sample", varying = colnames(Filter(is.numeric, df_alignment())),
                  v.names = "Value", times = colnames(Filter(is.numeric, df_alignment())), direction = "long")
    if (input$align_var %in% c("Number of CpG dropout", "Number of Amplicon dropout", 
                               "Number of CpG dropout with full read", "Number of Amplicon dropout with full read")) {
      total <- strsplit(df[1, grep(paste0(gsub("Number of ", "No. ", input$align_var), "$"), colnames(df))], "/")[[1]][2]
      plot_ly(df[df$time == input$align_var,], y = ~Value, 
              type = "box", boxmean = T, hoverinfo = "y")%>%  
        layout(title = paste0(input$align_var, "\n", "Total number of ", 
                              gsub("Number of | dropout", "", input$align_var), ": ", total),
               xaxis = list(title = "", showticklabels = F),
               yaxis = list(title = "", zeroline = F))
    } else {
      plot_ly(df[df$time == input$align_var,], y = ~Value, 
              type = "box", boxmean = T, hoverinfo = "y")%>%
        layout(title = input$align_var,
               xaxis = list(title = "", showticklabels = F),
               yaxis = list(title = "", zeroline = F))
    }
  })
  output$align_table <- renderDataTable({
    req(df_align_download())
    df_align_download()[,c(1, 10, 11, 12, 13)]
  })
  ##### Metrics for group #####
  output$align_guide2 <- renderUI({
    req(credentials()$user_auth)
    p("1. You will need to annotate the groups of samples with the template in the download file")
  })
  output$align_group_button <- renderUI({
    req(credentials()$user_auth)
    downloadButton(outputId = "align_template",
                   label = "Download template")
  })
  output$align_template <- downloadHandler(
    filename = "align_template.tsv",
    content = function(con) {
      write.table(group_template(), sep ="\t", row.names = F, con)
    }
  )
  output$align_upload <- renderUI({
    req(credentials()$user_auth)
    fileInput(inputId = "align_upload", 
              label = "2. Upload the annotated file here",
              multiple = F)
  })
  # Process uploaded file
  align_group <- reactive({
    req(input$align_upload)
    t <- read.table(file = input$align_upload$datapath, 
                    header = T, stringsAsFactors = F, check.names = F, fill = T)
    t <- t[t$Group!="",]
    t
  })
  df_align_merge <- reactive({
    t <- merge(df_alignment(), align_group(), by = "Sample")
    t
  })
  output$align_var_group <- renderUI({
    selectInput(inputId = "align_var_group",
                label = "Select a metric to plot",
                choices = colnames(Filter(is.numeric, df_align_merge())))
  })
  output$align_group_plot <- renderPlotly({
    req(input$align_var_group)
    df <- reshape(df_align_merge(), idvar = "Sample", varying = colnames(Filter(is.numeric, df_align_merge())),
                  v.names = "Value", times = colnames(Filter(is.numeric, df_align_merge())), direction = "long")
    if (input$align_var_group %in% c("Number of CpG dropout", "Number of Amplicon dropout", 
                                     "Number of CpG dropout with full read", "Number of Amplicon dropout with full read")) {
      total <- strsplit(df[1, grep(paste0(gsub("Number of ", "No. ", input$align_var_group), "$"), colnames(df))], "/")[[1]][2]
      plot_ly(df[df$time == input$align_var_group,], x = ~Group, y = ~Value, 
              type = "box", boxmean = T, color = ~Group, hoverinfo = "y")%>%  
        layout(title = paste0(input$align_var_group, "\n", "Total number of ", 
                              gsub("Number of | dropout", "", input$align_var_group), ": ", total),
               xaxis = list(title = ""),
               yaxis = list(title = "", zeroline = F),
               showlegend = F)
    } else {
      plot_ly(df[df$time == input$align_var_group,], x = ~Group, y = ~Value, 
              type = "box", boxmean = T, color = ~Group, hoverinfo = "y")%>%
        layout(title = input$align_var_group,
               xaxis = list(title = ""),
               yaxis = list(title = "", zeroline = F),
               showlegend = F)
    }
  })
  ##### DNA methylation #####
  ##### Overview #####
  output$overview_select <- renderUI({
    req(credentials()$user_auth)
    selectInput(inputId = "overview_select",
                label = "Select a metric",
                choices = c("Coverage" = 1, "Methylation" = 2))
  })
  overview_heatmapDNAMet <- reactive({
    req(input$overview_select)
    if (input$overview_select == 1) {
      df_cover <- read.table(paste(p, credentials()$info$user, input$project, "bigTable/bigTable_filtered.tsv.coverage.tsv", sep = "/"), 
                             header = T, stringsAsFactors = F, check.names = F)
      rownames(df_cover) <- df_cover$Group.1
      df_cover <- df_cover[,-1]
      df_cover <- as.matrix(df_cover)
      colnames(df_cover) <- gsub(".cov", "", colnames(df_cover))
      p <- plot_ly(x=rownames(df_cover), y=colnames(df_cover), z = t(df_cover), 
                   type = "heatmap", 
                   hovertemplate = paste0("<b>Sample:</b> %{y}<br>",
                                          "<b>Amplicon:</b> %{x}<br>",
                                          "<b>Coverage:</b> %{z}<br>",
                                          "<extra></extra>"),
                   height = 800)%>%
        layout(title = "Coverage of all samples")
    } else {
      req(credentials()$user_auth)
      df_meth <- read.table(paste(p, credentials()$info$user, input$project, "bigTable/bigTable_fullRead.tsv.methylation.csv", sep = "/"), 
                            header = T, stringsAsFactors = F, check.names = F)
      rownames(df_meth) <- df_meth$Group.1
      df_meth <- df_meth[,-1]
      df_meth <- as.matrix(df_meth)
      colnames(df_meth) <- gsub(".ratio", "", colnames(df_meth))
      p <- plot_ly(x=rownames(df_meth), y=colnames(df_meth), z = t(df_meth), 
                   type = "heatmap", 
                   hovertemplate = paste0("<b>Sample:</b> %{y}<br>",
                                          "<b>Amplicon:</b> %{x}<br>",
                                          "<b>Methylation:</b> %{z}<br>",
                                          "<extra></extra>"),
                   height = 800)%>%
        layout(title = "DNA methylation of all samples")
    }
  })
  output$overview_heatmap <- renderPlotly({
    overview_heatmapDNAMet()
  })
  ##### MDS overview #####
  df_mds <- reactive({
    df <- bigTable_filtered()
    df1 <- df[, grep("ratio", colnames(df))]
    dm <- 1 - abs(cor(df1, use = "complete.obs"))
    fit <- cmdscale(dm, eig = TRUE, k = 2)
    df_plot <- as.data.frame(fit$points)
    df_plot$Sample <- row.names(df_plot)
    df_plot$Sample <- gsub(".ratio", "", df_plot$Sample)
    df_plot
  })
  mds_overviewDNAMet <- reactive({
    plot_ly(df_mds(), x = ~V1, y = ~V2, text = ~Sample, color = "Set2",
            type = "scatter", mode = "markers", 
            marker = list(size = 12), 
            hovertemplate = paste0("<b>Sample:</b> %{text}",
                                   "<extra></extra>"),
            height = 600)%>%
      add_text(textposition = "top right")%>%
      layout(title = "MDS plot",
             xaxis = list(title = "PC1",
                          zeroline = F, range = c(min(df_mds()$V1) - 0.2,
                                                  max(df_mds()$V1) + 0.2)),
             yaxis = list(title = "PC2", zeroline = F),
             showlegend = F)
  })
  output$mds_overview <- renderPlotly({
    mds_overviewDNAMet()
  })
  output$downloadDNAMetPlot <- downloadHandler(
    filename = "DNAMetplot.pdf",
    content = function(file){
      if (input$heatmap_table == 1) {
        orca(overview_heatmapDNAMet(), file = file.path(p, credentials()$info$user, input$project, "tempPlot.pdf"))
        file.copy(file.path(p, credentials()$info$user, input$project, "tempPlot.pdf"), file)
      } else if (input$heatmap_table == 3) {
        orca(mds_overviewDNAMet(), file = file.path(p, credentials()$info$user, input$project, "tempPlot.pdf"))
        file.copy(file.path(p, credentials()$info$user, input$project, "tempPlot.pdf"), file)
      }
    }
  )
  ##### By group #####
  # Guides and download button
  output$meth_guide <- renderUI({
    req(credentials()$user_auth)
    p("1. You will need to annotate the groups of samples with the template in the download file")
  })
  output$download_template <- renderUI({
    req(credentials()$user_auth)
    downloadButton(outputId = "meth_template",
                   label = "Download template")
  })
  # Process to download template data
  output$meth_template <- downloadHandler(
    filename = "meth_template.tsv",
    content = function(con) {
      write.table(group_template(), sep ="\t", row.names = F, con)
    }
  )
  # Upload annotation file
  output$file1 <- renderUI({
    req(credentials()$user_auth)
    fileInput(inputId = "file1", 
              label = "2. Upload the annotated file here",
              multiple = F)
  })
  output$meth_group_mode <- renderUI({
    req(input$file1)
    radioButtons(inputId = "meth_group_mode",
                 label = "Select a data visualization",
                 choices = c("MDS plot" = 2, "Box plot" = 1),
                 selected = 2)
  })
  # Process uploaded file
  dataset_group <- reactive({
    req(input$file1)
    t <- read.table(file = input$file1$datapath, 
                    header = T, stringsAsFactors = F, check.names = F, fill = T)
    t <- t[t$Group!="",]
    t
  })
  df_meth <- reactive({
    cov <- grep(".ratio", colnames(bigTable_filtered()), value = T)
    t <- bigTable_filtered()[bigTable_filtered()$InAmp==1, cov]
    t[t=="."] <- NA
    t <- as.data.frame(sapply(t, as.integer))
    cov_ave <- aggregate(t, list(bigTable_filtered()[bigTable_filtered()$InAmp==1,"Amplicon"]), function(x){round(mean(x),0)})
    colname <- gsub("\\.ratio", "", grep("\\.ratio", colnames(cov_ave), value=T))
    rowname <- cov_ave$Group.1
    a <- melt(cov_ave)
    a$Amp <- sapply(a$Group.1, function(x){ strsplit(x, "::")[[1]][1] })
    a <- na.omit(a)
    a$variable <-  gsub(".ratio$", "", a$variable)
    names(a)[names(a) == "Group.1"] <- "Amplicon"
    names(a)[names(a) == "variable"] <- "Sample"
    names(a)[names(a) == "value"] <- "Methylation"
    a
  })
  # Select an amplicon for plotting
  output$meth_amplicon <- renderUI({
    req(input$file1)
    selectInput(inputId = "Amplicon_in_each_group",
                label = "Select an amplicon to plot DNA methylation of the groups of samples", 
                choices = unique(df_meth()$Amplicon), 
                multiple = F)
  })
  # Plotting DNA methylation for selected amplicon by group
  output$Amplicon_group_violin <- renderPlotly({
    req(dataset_group())
    a <- df_meth()
    a$Group <- dataset_group()$Group[match(a$Sample, dataset_group()$Sample)]
    a2 <- a
    plot_ly(a2[a2$Amplicon==input$Amplicon_in_each_group,],x= ~Group ,y = ~Methylation, 
            type = "box", boxmean = T, color = ~Group,
            hoverinfo = "y")%>%
      layout(title = paste0("Amplicon: ", input$Amplicon_in_each_group),
             xaxis = list(title = ""),
             showlegend = F)
  })
  ##### Meth t test 1 #####
  output$meth_group_ttest <- renderPrint({
    req(input$Amplicon_in_each_group)
    a <- df_meth()[df_meth()$Amplicon == input$Amplicon_in_each_group,]
    a$Group <- dataset_group()$Group[match(a$Sample, dataset_group()$Sample)]
    a <- na.omit(a)
    a$idnew <- paste(a$Amplicon, a$Sample, sep = "_")
    if (length(unique(a$Group)) == 2) {
      t_test <- t.test(a$Methylation ~ a$Group, var.equal = F, alternative = "two.sided")
      cat(paste("Two-sided t test:",
                paste("p-value =", t_test$p.value),
                sep = "\n"))
    } else {
      subset <- a[,c("idnew", "Methylation", "Group")]
      df_wide <- reshape(subset, idvar = "idnew", timevar = "Group", direction = "wide")
      colnames(df_wide) <- gsub("Methylation.", "", colnames(df_wide))
      p.mat <- multi.ttest(df_wide[,-1], var.equal = F, alternative = "two.sided")
      upper <- p.mat
      upper[upper.tri(p.mat)] <- ""
      upper <- as.data.frame(upper)
      cat("Two-sided t test:", "\n")
      print(upper)
    }
  })
  # Plotting DNA methylation by group
  output$group_violin <- renderPlotly({
    req(dataset_group())
    a <- df_meth()
    a$Group <- dataset_group()$Group[match(a$Sample, dataset_group()$Sample)]
    plot_ly(a,x= ~Group ,y = ~Methylation, 
            type = "box", boxmean = T, color = ~Group,
            hoverinfo = "y")%>%
      layout(title = "All amplicons: DNA methylation in sample groups",
             xaxis = list(title = ""),
             showlegend = F)
  })
  ##### Meth t test 2 #####
  output$meth_group_ttest2 <- renderPrint({
    req(dataset_group())
    a <- df_meth()
    a$Group <- dataset_group()$Group[match(a$Sample, dataset_group()$Sample)]
    a <- na.omit(a)
    a$idnew <- paste(a$Amplicon, a$Sample, sep = "_")
    if (length(unique(a$Group)) == 2) {
      t_test <- t.test(a$Methylation ~ a$Group, var.equal = F, alternative = "two.sided")
      cat(paste("Two-sided t test:",
                paste("p-value =", t_test$p.value),
                sep = "\n"))
    } else {
      subset <- a[,c("idnew", "Methylation", "Group")]
      df_wide <- reshape(subset, idvar = "idnew", timevar = "Group", direction = "wide")
      colnames(df_wide) <- gsub("Methylation.", "", colnames(df_wide))
      p.mat <- multi.ttest(df_wide[,-1], var.equal = F, alternative = "two.sided")
      upper <- p.mat
      upper[upper.tri(p.mat)] <- ""
      upper <- as.data.frame(upper)
      cat("Two-sided t test:", "\n")
      print(upper)
    }
  })
  # Plotting DNA methylation for all amplicon by group
  output$Amplicon_all <- renderPlotly({
    req(dataset_group())
    a <- df_meth()
    a$Group <- dataset_group()$Group[match(a$Sample, dataset_group()$Sample)]
    a4 <- na.omit(a)
    plot_ly(a4, x = ~Amp, y = ~Methylation, color = ~Group, 
            type = "box", boxmean = T) %>%
      layout(boxmode = "group",
             title = "Box Plots for all amplicons",
             xaxis = list(title = ""),
             showlegend = F)
  })
  ##### MDS plot by group #####
  df_mds_merge <- reactive({
    merge(df_mds(), dataset_group(), by = "Sample")
  })
  plot_mds_group <- reactive({
    plot_ly(df_mds_merge(), x = ~V1, y = ~V2, text = ~Sample, color = ~Group, 
            type = "scatter", mode = "markers", 
            marker = list(size = 12), 
            hovertemplate = paste0("<b>Sample:</b> %{text}<br>",
                                   "<b>Group:</b> %{fullData.name}<br>",
                                   "<extra></extra>"),
            height = 600)%>%
      add_text(textposition = "top right", showlegend = F)%>%
      layout(title = "MDS plot",
             xaxis = list(title = "PC1",
                          zeroline = F, range = c(min(df_mds_merge()$V1) - 0.2,
                                                  max(df_mds_merge()$V1) + 0.2)),
             yaxis = list(title = "PC2", zeroline = F),
             showlegend = T)
  })
  output$mds_group <- renderPlotly({
    plot_mds_group()
  })
  output$pdf_mds_group <- downloadHandler(
    filename = "MDS_group.pdf",
    content = function(file){
      orca(plot_mds_group(), file = "tempPlot.pdf")
      file.copy("tempPlot.pdf", file)
    }
  )
  ##### By amplicon #####
  output$bigtable_guide1 <- renderUI({
    req(credentials()$user_auth)
    p("Click here to download the big table")
  })
  output$bigtable_dlbutton <- renderUI({
    req(credentials()$user_auth)
    downloadButton(outputId = "bigtable_download",
                   label = "Download big table")
  })
  output$bigtable_download <- downloadHandler(
    filename = "bigTable_filtered.tsv",
    content = function(con) {
      file.copy(paste(p, credentials()$info$user, input$project, "bigTable/bigTable_filtered.tsv", sep = "/"), con)
    }
  )
  output$single_all <- renderUI({
    req(credentials()$user_auth)
    radioButtons(inputId = "single_all", 
                 label = "", 
                 choices = c("All amplicons" = 2, "Single amplicon" = 1))
  })
  output$Amplicon_methylation <- renderUI({
    req(df_meth())
    req(input$single_all == 1)
    selectInput(inputId = "Amplicon_methylation", 
                label = "Select an amplicon", 
                choices = unique(df_meth()$Amplicon))
  })
  dataset_methylation <- reactive({
    df_meth()[df_meth()$Amplicon == input$Amplicon_methylation,]
  })
  #Plot for methylation in 1 amplicon
  output$Methylation_violin <- renderPlotly({
    req(input$single_all)
    if(input$single_all == 1) {
      plot_ly(dataset_methylation(), y = ~Methylation, 
              type = "box", boxmean = T, hoverinfo = "y")%>%
        layout(title = paste0("Amplicon: ", input$Amplicon_methylation),
               xaxis = list(title = "", showticklabels = F))
    } else {
      plot_ly(df_meth(), x = ~Amp , y = ~Methylation, 
              type = "box", boxmean = T, color = ~Amp)%>%
        layout(title = "DNA methylation of all amplicons",
               xaxis = list(title = ""))
    }
  })
  ##### By sample #####
  output$meth_sample_select <- renderUI({
    req(credentials()$user_auth)
    selectInput(inputId = "meth_sample_select",
                label = "Select a sample to plot",
                choices = unique(df_meth()$Sample))
  })
  output$meth_sample_plot <- renderPlotly({
    plot_ly(df_meth()[df_meth()$Sample == input$meth_sample_select,], x = ~Amplicon, y = ~Methylation, 
            type = "scatter", mode = "lines+markers", 
            hovertemplate = paste0("<b>Amplicon:</b> %{x}<br>",
                                   "<b>Methylation:</b> %{y}<br>",
                                   "<extra></extra>"),
            height = 700)%>%
      layout(title = paste0("Sample: ", input$meth_sample_select),
             xaxis = list(title = ""),
             yaxis = list(range = c(0, 105)))
  })
  ##### Pattern #####
  output$selectproj <- renderUI({
    req(credentials()$user_auth)
    selectInput(inputId = "pattern_sample",
                label = "Select a sample to plot DNA methylation pattern",
                choices = list.dirs(paste(p, credentials()$info$user, input$project, "patterned", sep = "/"),
                                    full.names = F, recursive = FALSE))
  })
  output$result1 <- renderUI({
    req(input$pattern_sample)
    selectInput(inputId = "pattern_file", 
                label = "Select an amplicon to plot DNA methylation pattern", 
                choices = gsub(".tsv", "", 
                                  list.files(paste(p, credentials()$info$user, input$project, "patterned", input$pattern_sample, "pattern", sep = "/"), 
                                          full.names = F, pattern = "*.tsv")))
  })
  # Read pattern data
  data_pattern_1 <- reactive({
    req(input$pattern_file)
    data_pattern <- read.table(paste(p, credentials()$info$user, input$project, "patterned", input$pattern_sample, "pattern", 
                                     paste0(input$pattern_file, ".tsv"), sep = "/"), 
                               header = T, stringsAsFactors = F, check.names = F)
  })
  # Pattern plot
  heatmapPattern <- reactive({
    req(data_pattern_1())
    df <- data_pattern_1()
    if(dim(df)[2] == 3) {
      # Pattern part
      df_plot <- as.data.frame(df[-dim(df)[1], 1])
      colnames(df_plot) <- "base"
      df_plot$position <- colnames(df)[1]
      df_plot$mat <- ifelse(df_plot$base == "T", 0, ifelse(df_plot$base == "C", 100, 50))
      df_plot$id <- seq_along(df_plot[,1])
      
      pattern <- plot_ly(df_plot, x = ~position, y = ~id, z = ~mat, 
                         type = "heatmap", showscale = F, hoverinfo = "none", colors = "Set3",
                         height = 300)%>%
        add_annotations(text = ~base, font = list(color = "black"),
                        showarrow = F)%>%
        layout(xaxis = list(title = ""),
               yaxis = list(title = "", autorange = "reversed", 
                            tickcolor = toRGB("white"), showticklabels = F, zeroline = F))
      # Frequency part
      df2 <- df[-dim(df)[1], grep("Frequency|Percentage", colnames(df))]
      df_freq <- reshape(df2, varying = colnames(df2),
                         v.names = "value", times = colnames(df2),
                         direction = "long")
      colnames(df_freq)[colnames(df_freq) == "time"] <- "metric"
      freq <- plot_ly(df_freq, x = ~metric, y = ~id, z = ~value, 
                      type = "heatmap", showscale = F, hoverinfo = "none",
                      height = 300)%>%
        add_annotations(text = ~value, font = list(color = "black"),
                        showarrow = F)%>%
        layout(xaxis = list(title = ""),
               yaxis = list(title = "", autorange = "reversed", showticklabels = F, zeroline = F))
      
      # Combine 2 parts
      subplot(pattern, freq, widths = c(0.7, 0.3), nrows = 1, shareX = F, shareY = T)
    } else {
      # Pattern part
      df1 <- df[-dim(df)[1], grep("Frequency|Percentage", colnames(df), invert = T)]
      df_plot <- reshape(df1, varying = colnames(df1),
                         v.names = "base", times = colnames(df1),
                         direction = "long")
      colnames(df_plot)[colnames(df_plot) == "time"] <- "position"
      df_plot$position <- sprintf("%03d", as.numeric(df_plot$position))
      df_plot$mat <- ifelse(df_plot$base == "T", 0, ifelse(df_plot$base == "C", 100, 50))
      
      pattern <- plot_ly(df_plot, x = ~position, y = ~id, z = ~mat, 
                         type = "heatmap", showscale = F, hoverinfo = "none", colors = "Set3",
                         height = 600)%>%
        add_annotations(text = ~base, font = list(color = "black"),
                        showarrow = F)%>%
        layout(xaxis = list(title = ""),
               yaxis = list(title = "", autorange = "reversed", 
                            tickcolor = toRGB("white"), showticklabels = F, zeroline = F))
      
      # Frequency part
      df2 <- df[-dim(df)[1], grep("Frequency|Percentage", colnames(df))]
      df_freq <- reshape(df2, varying = colnames(df2),
                         v.names = "value", times = colnames(df2),
                         direction = "long")
      colnames(df_freq)[colnames(df_freq) == "time"] <- "metric"
      
      freq <- plot_ly(df_freq, x = ~metric, y = ~id, z = ~value, 
                      type = "heatmap", showscale = F, hoverinfo = "none", colors = "Set3",
                      height = 600)%>%
        add_annotations(text = ~value, font = list(color = "black"),
                        showarrow = F)%>%
        layout(xaxis = list(title = ""),
               yaxis = list(title = "", autorange = "reversed", showticklabels = F, zeroline = F))
      
      # Combine 2 parts
      subplot(pattern, freq, widths = c(0.7, 0.3), nrows = 1, shareX = F, shareY = T)
    }
  })
  output$Heatmap <- renderPlotly({
    heatmapPattern()
  })
  output$downloadPatternHeatmap <- downloadHandler(
    filename = "heatmapPattern.pdf",
    content = function(file){
      orca(heatmapPattern(), file = "tempPlot.pdf")
      file.copy("tempPlot.pdf", file)
    }
  )
  ##### Polymorphism #####
  data_polymorphism <- reactive({
    req(credentials()$user_auth)
    read.table(paste(p, credentials()$info$user, input$project, "bigTable/bigTable_fullRead_polymorphism.tsv", sep = "/"), 
               header=TRUE, sep="\t", comment.char = "\\", stringsAsFactors=F, check.names=FALSE)
  })
  output$poly_single_all <- renderUI({
    req(credentials()$user_auth)
    radioButtons(inputId = "poly_single_all", 
                 label = "", 
                 choices = c("All amplicons" = 2, "Single amplicon" = 1))
  })
  output$polymorphism_amp <- renderUI({
    req(input$poly_single_all == 1)
    selectInput(inputId = "poly_amplicon", 
                label = "Please select an amplicon", 
                choices = unique(data_polymorphism()$Amplicon))
  })
  # Data for Polymorphism_all
  tb <- reactive({
    req(data_polymorphism())
    data_polymorphism <- data_polymorphism()
    dt <- data_polymorphism[data_polymorphism$Type %in% c("mean_ratio", "polymorphism"),]
    l <- lapply(unique(dt$Amplicon), function(x) {
      dt1 <- dt[dt$Amplicon==x,-c(1,2)]
      dt1 <- as.data.frame(t(dt1))
      colnames(dt1) <- dt[dt$Amplicon==x,2]
      dt1$Sample <- row.names(dt1)
      dt1$Amplicon <- x
      return(dt1)
    })
    tb <- Reduce(function(x, y) rbind(x,y), l)
    colnames(tb)[colnames(tb) == "mean_ratio"] <- "Methylation"
    colnames(tb)[colnames(tb) == "polymorphism"] <- "Polymorphism"
    tb
  })
  # Plot for polymorphism
  plot_polymorphism <- reactive({
    req(input$poly_single_all)
    if (input$poly_single_all == 1) {
      df <- tb()[tb()$Amplicon == input$poly_amplicon,]
      plot_ly(df, x = ~Methylation, y = ~Polymorphism, 
              color = ~Sample, opacity = 0.7,
              type = "scatter", mode = "markers", marker = list(size = 12),
              hovertemplate = paste0("<b>Sample:</b> %{fullData.name}<br>",
                                     "<b>Polymorphism:</b> %{y}<br>",
                                     "<b>Methylation:</b> %{x}<br>",
                                     "<extra></extra>")) %>%
        layout(title = paste0("Amplicon: ", input$poly_amplicon),
               xaxis = list(title = "Methylation",
                            range = c(-3, 103)),
               yaxis = list(title = "Polymorphism score",
                            range = c(-3, 103))
        )
    } else {
      plot_ly(tb(), x = ~Methylation, y = ~Polymorphism, text = ~Amplicon,
              color = ~Sample, opacity = 0.7,
              type = "scatter", mode = "markers", marker = list(size = 12),
              hovertemplate = paste0("<b>Sample:</b> %{fullData.name}<br>",
                                     "<b>Amplicon:</b> %{text}<br>",
                                     "<b>Polymorphism:</b> %{y}<br>",
                                     "<b>Methylation:</b> %{x}<br>",
                                     "<extra></extra>")) %>%
        layout(title = "Polymorphism score of all amplicons",
               xaxis = list(title = "Methylation",
                            range = c(-3, 103)),
               yaxis = list(title = "Polymorphism score",
                            range = c(-3, 103))
        )
    }
  })
  output$Polymorphism_scatter_plot <- renderPlotly({
    plot_polymorphism()
  })
  output$Polymorphism_table <- DT::renderDataTable({
    req(input$poly_single_all == 1)
    df <- tb()[tb()$Amplicon == input$poly_amplicon,]
    row.names(df) <- df$Sample
    df[, grep("Amplicon|Sample", colnames(df), invert = T)]
  })
  output$pdf_polymorphism <- downloadHandler(
    filename = "polymorphism.pdf",
    content = function(file){
      withr::with_dir("/srv/shiny-server/sample-apps/temp/", orca(plot_polymorphism(), file = "tempPlot.pdf"))
      file.copy("tempPlot.pdf", file)
    }
  )
  ##### PCR bias correction #####
  ##### For amplicons #####
  output$pcr_guide <- renderUI({
    req(credentials()$user_auth)
    p("1. Please download this template and fill in the DNA methylation values (Expectation) 
      of the control spike-in samples")
  })
  output$download_pcr <- renderUI({
    req(credentials()$user_auth)
    downloadButton(outputId = "pcr_control",
                   label = "Download template")
  })
  # Generate template data
  pcr_template <- reactive({
    pct_temp <- as.data.frame(df_qc()$Sample)
    colnames(pct_temp)[1] <- "Sample"
    pct_temp$Expectation <- ""
    pct_temp
  })
  # Process to download template data
  output$pcr_control <- downloadHandler(
    filename = "pcr_template.tsv",
    content = function(con) {
      write.table(pcr_template(), sep ="\t", row.names = F, con)
    }
  )
  output$upload_expected <- renderUI({
    req(credentials()$user_auth)
    fileInput(inputId = "expected_file", 
              label = "2. Upload your filled template here to perform PCR bias correction",
              multiple = F)
  })
  expected <- reactive({
    req(input$expected_file)
    t <- read.table(file = input$expected_file$datapath, 
                    header = T, stringsAsFactors = F, check.names = F, fill = T)
    t <- t[!is.na(t$Expectation),]
    t
  })
  output$pcr_amplicon <- renderUI({
    req(input$expected_file)
    selectInput(inputId = "pcr_amplicon",
                label = "Select an amplicon",
                choices = unique(bigTable_filtered()$Amplicon))
  })
  observed <- reactive({
    table <- bigTable_filtered()[bigTable_filtered()$Amplicon == input$pcr_amplicon, 
                                 grep("ratio", colnames(bigTable_filtered()))]
    table <- as.data.frame(apply(table, 2, function(x) mean(as.numeric(x), na.rm = T)))
    colnames(table) <- "Observation"
    table$Sample <- row.names(table)
    table$Sample <- gsub(".ratio", "", table$Sample)
    table
  })
  merge_oe <- reactive({
    df_merge <- merge(expected(), observed(), by = "Sample")
    df_merge <- df_merge[order(df_merge$Expectation),]
    df_merge
  })
  b_uncorrected <- reactive({
    x <- merge_oe()$Expectation
    y <- merge_oe()$Observation
#    y0 <- min(y)
#    y1 <- max(y)
    y0 <- 0
    y1 <- 100
    nls(y ~ ((b*y1-y0)*x+100*y0)/(b*x-x+100), start = list(b=1),control=nls.control(maxiter = 500, warnOnly=TRUE),na.action=na.exclude)
  })
  # Display all b_uncorrected
  output$plot_b_uncorrected <- renderPlotly({
    As <- unique(bigTable_filtered()$Amplicon)
    names(As) <- unique(bigTable_filtered()$Amplicon)
    bs <- sapply(As, function(A) {
      t <- bigTable_filtered()[bigTable_filtered()$Amplicon==A, grep("ratio", colnames(bigTable_filtered()))]
      observe <- as.data.frame(apply(t, 2, function(B) mean(as.numeric(B), na.rm = T)))
      colnames(observe) <- "Observation"
      observe$Sample <- row.names(observe)
      observe$Sample <- gsub(".ratio", "", observe$Sample)
      
      df_merge <- merge(expected(), observe, by = "Sample")
      df_merge <- df_merge[order(df_merge$Expectation),]
      
      x <- df_merge$Expectation
      y <- df_merge$Observation
      if(!any(is.na(y))) {
#        y0 <- min(y) #observed methylation levels at 0%
#        y1 <- max(y) #observed methylation levels at 100%
        y0 <- 0
        y1 <- 100
        b_uncorrected <- nls(y ~ ((b*y1-y0)*x+100*y0)/(b*x-x+100), start = list(b = 1), 
                             control = nls.control(maxiter = 500, warnOnly = TRUE), na.action = na.exclude) #the formula
        return(as.numeric(coef(b_uncorrected)))
      } else {
        return("NA")
      }
    })
    df_bs <- as.data.frame(cbind("Amplicon" = names(bs), "b_uncorrected" = bs),
                           stringsAsFactors = F)
    df_bs$b_uncorrected <- as.numeric(df_bs$b_uncorrected)
    df_bs$Amplicon <- factor(df_bs$Amplicon, levels = df_bs$Amplicon[order(df_bs$b_uncorrected)])
    plot_ly(df_bs, x = ~b_uncorrected, y = ~Amplicon,
      color = ~Amplicon, name = ~Amplicon,
      type = "bar"
    )%>%
      layout(title = "PCR bias values of all amplicons",
             xaxis = list(title = ""),
             yaxis = list(title = ""),
             shapes = list(type = "rect", 
                           fillcolor = "blue", line = list(color = "blue"), opacity = 0.3,
                           x0 = 0.5, x1 = 1.5, xref = "x",
                           y0 = -0.5, y1 = length(df_bs$Amplicon), yref = "y"),
             showlegend = F)
  })
  # Before corrected spikes-in
  output$spikesin_before <- renderPlotly({
    plot_ly(merge_oe(), x = ~Expectation, y = ~Observation,
            type = "scatter", mode = "lines+markers")%>%
      add_trace(x = c(0, 100), y = c(0, 100), type = "scatter", mode = "lines",
                line = list(dash = "dash", color = "grey"))%>%
      layout(title = "Before corrected of spike-in samples", showlegend = F,
             annotations = list(text = paste0("b = ", coef(b_uncorrected())), 
                                x = 80, y = 10, showarrow = F, font = list(size = 16)))
  })
  # After corrected spikes-in
  output$spikesin_after <- renderPlotly({
    y <- merge_oe()$Observation
    names(y) <- merge_oe()$Sample
#    y0 <- min(y)
#    y1 <- max(y)
    y0 <- 0
    y1 <- 100
    b <- b_uncorrected()
    y_corr<- (100*y0-100*y)/(coef(b)*y- coef(b)*y1+y0-y)
#    y0_corr <- min(y_corr)
#    y1_corr <- max(y_corr)
    y0_corr <- 0
    y1_corr <- 100
    x <- merge_oe()$Expectation
    names(x) <- merge_oe()$Sample
    y_corr <- y_corr[match(names(y_corr), names(x))]
    b_corrected <- nls(y_corr ~ ((b*y1_corr-y0_corr)*x+100*y0_corr)/(b*x-x+100), start = list(b=1),control=nls.control(maxiter = 500, warnOnly=TRUE),na.action=na.exclude) #the formula
  
    df_xy <- as.data.frame(cbind("Expectation" = x, "Corrected" = y_corr))
    plot_ly(df_xy, x = ~Expectation, y = ~Corrected,
            type = "scatter", mode = "lines+markers")%>%
      add_trace(x = c(0, 100), y = c(0, 100), type = "scatter", mode = "lines",
                line = list(dash = "dash", color = "grey"))%>%
      layout(title = "After corrected of spike-in samples", showlegend = F,
             yaxis = list(title = "Observation"),
             annotations = list(text = paste0("b = ", coef(b_corrected)), 
                                x = 80, y = 10, showarrow = F, font = list(size = 16)))
  })
  # Before-after plot for samples
  output$beforeafter_plot <- renderPlotly({
    t <- bigTable_filtered()[bigTable_filtered()$Amplicon == input$pcr_amplicon, grep("ratio", colnames(bigTable_filtered()))]
    observe <- as.data.frame(apply(t, 2, function(x) mean(as.numeric(x), na.rm = T)))
    colnames(observe) <- "Observation"
    observe$Sample <- row.names(observe)
    observe$Sample <- gsub(".ratio", "", observe$Sample)
    
    df_merge <- merge(expected(), observe, by = "Sample")
    df_merge <- df_merge[order(df_merge$Expectation),]
    
    x <- df_merge$Expectation
    y <- df_merge$Observation
#    y0 <- min(y) #observed methylation levels at 0%
#    y1 <- max(y) #observed methylation levels at 100%
    y0 <- 0
    y1 <- 100
    b_uncorrected <- nls(y ~ ((b*y1-y0)*x+100*y0)/(b*x-x+100), start = list(b=1),control=nls.control(maxiter = 500, warnOnly=TRUE),na.action=na.exclude) #the formula
    
    samples <- setdiff(observe$Sample, df_merge$Sample)
    
    y <- observe[observe$Sample %in% samples, "Observation"]
    names(y) <- observe[observe$Sample %in% samples, "Sample"]
    
    y0 <- min(y, na.rm = T)
    y1 <- max(y, na.rm = T)
    b <- b_uncorrected 
    
    y_corr<- (100*y0-100*y)/(coef(b)*y- coef(b)*y1+y0-y) #correcting the data using eq 3
    
    observe_correct <- data.frame(cbind("Before" = y[match(names(y), names(y_corr))], "After" = y_corr))
    observe_correct$Sample <- row.names(observe_correct)
    
    # Line plot
    plot_ly(observe_correct, x = ~Sample, y = ~Before, 
            type = "scatter", mode = "lines+markers", name = "Before")%>%
      add_trace(y = ~After, name = "After", mode = "lines+markers")%>%
      layout(title = "Non spike-in samples", xaxis = list(title = ""),
             yaxis = list(title = "Methylation (%"))
  })
  
  df_all_corrected_amps <- reactive({
    As <- unique(bigTable_filtered()$Amplicon)
    l <- lapply(As, function(A) {
      t <- bigTable_filtered()[bigTable_filtered()$Amplicon==A, grep("ratio", colnames(bigTable_filtered()))]
      observe <- as.data.frame(apply(t, 2, function(x) mean(as.numeric(x), na.rm = T)))
      colnames(observe) <- "Observation"
      observe$Sample <- row.names(observe)
      observe$Sample <- gsub(".ratio", "", observe$Sample)
      df_merge <- merge(expected(), observe, by = "Sample")
      y <- df_merge$Observation
      names(y) <- df_merge$Sample
      x <- df_merge$Expectation
      names(x) <- df_merge$Sample
      if(!any(is.na(y))) {
#        y0 <- min(y)
#        y1 <- max(y)
        y0 <- 0
        y1 <- 100
        b_uncorrected <- nls(y ~ ((b*y1-y0)*x+100*y0)/(b*x-x+100), start = list(b = 1), 
                             control = nls.control(maxiter = 500, warnOnly = TRUE), na.action = na.exclude) #the formula
        b <- b_uncorrected
        y <- observe$Observation
        names(y) <- observe$Sample
        y0 <- min(y, na.rm = T)
        y1 <- max(y, na.rm = T)
        y_corr<- (100*y0-100*y)/(coef(b)*y- coef(b)*y1+y0-y) #correcting the data using eq 3
        y_corr <- as.data.frame(y_corr)
        y_corr$Sample <- row.names(y_corr)
        y <- as.data.frame(y)
        y$Sample <- row.names(y)
        before_after <- merge(y, y_corr, by = "Sample")
        before_after$Amplicon <- A
        return(before_after)
      } else {
        return(NA)
      }
    })
    l1 <- do.call(rbind, l)
    l1 <- na.omit(l1)
    l2 <- reshape(l1, idvar = "Amplicon", v.names = c("y", "y_corr"), timevar = "Sample", direction = "wide")
    colnames(l2)[grep("y\\.", colnames(l2))] <- paste0(colnames(l2)[grep("y\\.", colnames(l2))], ".before")
    colnames(l2)[grep("before", colnames(l2))] <- gsub("y\\.", "", colnames(l2)[grep("before", colnames(l2))])
    colnames(l2)[grep("y_", colnames(l2))] <- paste0(colnames(l2)[grep("y_", colnames(l2))], ".after")
    colnames(l2)[grep("after", colnames(l2))] <- gsub("y_corr.", "", colnames(l2)[grep("after", colnames(l2))])
    l2
  })
  output$pcr_guide_2 <- renderUI({
    req(input$expected_file)
    p("PCR bias correction is optional. When b values exceed a limitation of 0.5 to 1.5 
    (the blue area), bias correction should be performed")
  })
  output$download_corrected_amps_guide <- renderUI({
    req(input$expected_file)
    p("You can click here to download the corrected DNA methylation values of all amplicons")
  })
  output$download_corrected_amps <- renderUI({
    req(input$expected_file)
    downloadButton(outputId = "all_corrected_amps",
                   label = "Download corrected amplicons")
  })
  output$all_corrected_amps <- downloadHandler(
    filename = "all_corrected_amps.tsv",
    content = function(con) {
      write.table(df_all_corrected_amps(), sep ="\t", row.names = F, con)
    }
  )
  ##### For CpGs #####
  output$pcr_guide_3 <- renderUI({
    req(input$expected_file)
    p("PCR bias correction is optional. When b values exceed a limitation of 0.5 to 1.5 
    (the blue area), bias correction should be performed")
  })
  df_obs <- reactive({
    df_obs <- bigTable_filtered()
    df_obs$Gpos_format <- sprintf("%03d", df_obs$Gpos)
    As_CpG <- paste(df_obs$Amplicon, df_obs$Gpos_format, sep="_")
    df_obs$Amplicon_CpG <- As_CpG
    df_obs
  })
  output$pcr_amplicon_CpGs <- renderUI({
    req(input$expected_file)
    selectInput(inputId = "pcr_amplicon_CpGs",
                label = "Select a CpG site",
                choices = df_obs()$Amplicon_CpG)
  })
  output$plot_b_uncorrected_CpGs <- renderPlotly({
    As_CpG <- df_obs()$Amplicon_CpG
    bs_CpG <- sapply(As_CpG, function(A) {
      t_CpG <- df_obs()[df_obs()$Amplicon_CpG==A, grep("ratio", colnames(df_obs()))]
      observe_CpG <- as.data.frame(apply(t_CpG, 2, as.numeric))
      colnames(observe_CpG) <- "Observation"
      observe_CpG$Sample <- row.names(observe_CpG)
      observe_CpG$Sample <- gsub(".ratio", "", observe_CpG$Sample)
      
      df_merge_CpG <- merge(expected(), observe_CpG, by = "Sample")
      df_merge_CpG <- df_merge_CpG[order(df_merge_CpG$Expectation),]
      
      x <- df_merge_CpG$Expectation
      y <- df_merge_CpG$Observation
      if(!any(is.na(y))) {
#        y0 <- min(y) #observed methylation levels at 0%
#        y1 <- max(y) #observed methylation levels at 100%
        y0 <- 0
        y1 <- 100
        b_uncorrected_CpG <- nls(y ~ ((b*y1-y0)*x+100*y0)/(b*x-x+100), start = list(b = 1), 
                                 control = nls.control(maxiter = 500, warnOnly = TRUE), na.action = na.exclude) #the formula
        return(as.numeric(coef(b_uncorrected_CpG)))
      } else {
        return("NA")
      }
    })
    df_bs_CpG <- as.data.frame(cbind("Amplicon_CpG" = names(bs_CpG), "b_uncorrected" = bs_CpG),
                               stringsAsFactors = F)
    df_bs_CpG$b_uncorrected_CpG <- as.numeric(df_bs_CpG$b_uncorrected)
    df_bs_CpG$Amplicon_CpG <- factor(df_bs_CpG$Amplicon, levels = df_bs_CpG$Amplicon[order(df_bs_CpG$b_uncorrected_CpG)])
    plot_ly(df_bs_CpG, x = ~b_uncorrected_CpG, y = ~Amplicon_CpG,
            color = ~Amplicon_CpG, name = ~Amplicon_CpG,
            type = "bar"
    )%>%
      layout(title = "PCR bias values of all CpG sites",
             xaxis = list(title = ""),
             yaxis = list(title = ""),
             shapes = list(type = "rect", 
                           fillcolor = "blue", line = list(color = "blue"), opacity = 0.3,
                           x0 = 0.5, x1 = 1.5, xref = "x",
                           y0 = -0.5, y1 = length(df_bs_CpG$Amplicon_CpG), yref = "y"),
             showlegend = F)
  })
  observed_CpGs <- reactive({
    table <- df_obs()[df_obs()$Amplicon_CpG == input$pcr_amplicon_CpGs, 
                      grep("ratio", colnames(df_obs()))]
    table <- as.data.frame(apply(table, 2, as.numeric))
    colnames(table) <- "Observation"
    table$Sample <- row.names(table)
    table$Sample <- gsub(".ratio", "", table$Sample)
    table
  })
  merge_oe_CpGs <- reactive({
    df_merge <- merge(expected(), observed_CpGs(), by = "Sample")
    df_merge <- df_merge[order(df_merge$Expectation),]
    df_merge
  })
  b_uncorrected_CpGs <- reactive({
    x <- merge_oe_CpGs()$Expectation
    y <- merge_oe_CpGs()$Observation
#    y0 <- min(y)
#    y1 <- max(y)
    y0 <- 0
    y1 <- 100
    nls(y ~ ((b*y1-y0)*x+100*y0)/(b*x-x+100), start = list(b=1),control=nls.control(maxiter = 500, warnOnly=TRUE),na.action=na.exclude)
  })
  # Before corrected spikes-in
  output$spikesin_before_CpGs <- renderPlotly({
    req(merge_oe_CpGs()$Expectation)
    req(b_uncorrected_CpGs())
    plot_ly(merge_oe_CpGs(), x = ~Expectation, y = ~Observation,
            type = "scatter", mode = "lines+markers")%>%
      add_trace(x = c(0, 100), y = c(0, 100), type = "scatter", mode = "lines",
                line = list(dash = "dash", color = "grey"))%>%
      layout(title = "Before corrected of spike-in samples", showlegend = F,
             annotations = list(text = paste0("b = ", coef(b_uncorrected_CpGs())), 
                                x = 80, y = 10, showarrow = F, font = list(size = 16)))
  })
  # After corrected spikes-in
  output$spikesin_after_CpGs <- renderPlotly({
    req(merge_oe_CpGs()$Expectation)
    req(b_uncorrected_CpGs())
    y <- merge_oe_CpGs()$Observation
    names(y) <- merge_oe_CpGs()$Sample
#    y0 <- min(y)
#    y1 <- max(y)
    y0 <- 0
    y1 <- 100
    b <- b_uncorrected_CpGs()
    y_corr<- (100*y0-100*y)/(coef(b)*y- coef(b)*y1+y0-y)
#    y0_corr <- min(y_corr)
#    y1_corr <- max(y_corr)
    y0_corr <- 0
    y1_corr <- 100
    x <- merge_oe_CpGs()$Expectation
    names(x) <- merge_oe_CpGs()$Sample
    y_corr <- y_corr[match(names(y_corr), names(x))]
    b_corrected <- nls(y_corr ~ ((b*y1_corr-y0_corr)*x+100*y0_corr)/(b*x-x+100), start = list(b=1),control=nls.control(maxiter = 500, warnOnly=TRUE),na.action=na.exclude) #the formula
    
    df_xy <- as.data.frame(cbind("Expectation" = x, "Corrected" = y_corr))
    plot_ly(df_xy, x = ~Expectation, y = ~Corrected,
            type = "scatter", mode = "lines+markers")%>%
      add_trace(x = c(0, 100), y = c(0, 100), type = "scatter", mode = "lines",
                line = list(dash = "dash", color = "grey"))%>%
      layout(title = "After corrected of spike-in samples", showlegend = F,
             annotations = list(text = paste0("b = ", coef(b_corrected)), 
                                x = 80, y = 10, showarrow = F, font = list(size = 16)))
  })
  # Before-after plot for CpGs
  output$beforeafter_plot_CpGs <- renderPlotly({
    samples <- setdiff(observed_CpGs()$Sample, merge_oe_CpGs()$Sample)
    y <- observed_CpGs()[observed_CpGs()$Sample %in% samples, "Observation"]
    names(y) <- observed_CpGs()[observed_CpGs()$Sample %in% samples, "Sample"]
    y0 <- min(y, na.rm = T)
    y1 <- max(y, na.rm = T)
    b <- b_uncorrected_CpGs()
    
    y_corr<- (100*y0-100*y)/(coef(b)*y- coef(b)*y1+y0-y) #correcting the data using eq 3
    
    observe_correct <- data.frame(cbind("Before" = y[match(names(y), names(y_corr))], "After" = y_corr))
    observe_correct$Sample <- row.names(observe_correct)
    
    # Line plot
    plot_ly(observe_correct, x = ~Sample, y = ~Before, 
            type = "scatter", mode = "lines+markers", name = "Before")%>%
      add_trace(y = ~After, name = "After", mode = "lines+markers")%>%
      layout(title = "Non spike-in samples", xaxis = list(title = ""),
             yaxis = list(title = "Methylation (%"))
  })
  df_all_corrected_CpGs <- reactive({
    As_CpG <- df_obs()$Amplicon_CpG
    l <- lapply(As_CpG, function(A) {
      t <- df_obs()[df_obs()$Amplicon_CpG==A, grep("ratio", colnames(df_obs()))]
      observe <- as.data.frame(apply(t, 2, as.numeric))
      colnames(observe) <- "Observation"
      observe$Sample <- row.names(observe)
      observe$Sample <- gsub(".ratio", "", observe$Sample)
      df_merge <- merge(expected(), observe, by = "Sample")
      y <- df_merge$Observation
      names(y) <- df_merge$Sample
      x <- df_merge$Expectation
      names(x) <- df_merge$Sample
      if(!any(is.na(y))) {
#        y0 <- min(y)
#        y1 <- max(y)
        y0 <- 0
        y1 <- 100
        b_uncorrected <- nls(y ~ ((b*y1-y0)*x+100*y0)/(b*x-x+100), start = list(b = 1), 
                             control = nls.control(maxiter = 500, warnOnly = TRUE), na.action = na.exclude) #the formula
        b <- b_uncorrected
        y <- observe$Observation
        names(y) <- observe$Sample
        y0 <- min(y, na.rm = T)
        y1 <- max(y, na.rm = T)
        y_corr<- (100*y0-100*y)/(coef(b)*y- coef(b)*y1+y0-y) #correcting the data using eq 3
        y_corr <- as.data.frame(y_corr)
        y_corr$Sample <- row.names(y_corr)
        y <- as.data.frame(y)
        y$Sample <- row.names(y)
        before_after <- merge(y, y_corr, by = "Sample")
        before_after$Amplicon <- A
        return(before_after)
      } else {
        return(NA)
      }
    })
    
    l1 <- do.call(rbind, l)
    l1 <- na.omit(l1)
    l2 <- reshape(l1, idvar = "Amplicon", v.names = c("y", "y_corr"), timevar = "Sample", direction = "wide")
    
    colnames(l2)[grep("y\\.", colnames(l2))] <- paste0(colnames(l2)[grep("y\\.", colnames(l2))], ".before")
    colnames(l2)[grep("before", colnames(l2))] <- gsub("y\\.", "", colnames(l2)[grep("before", colnames(l2))])
    colnames(l2)[grep("y_", colnames(l2))] <- paste0(colnames(l2)[grep("y_", colnames(l2))], ".after")
    colnames(l2)[grep("after", colnames(l2))] <- gsub("y_corr.", "", colnames(l2)[grep("after", colnames(l2))])
    l2
  })
  output$download_corrected_CpGs_guide <- renderUI({
    req(input$expected_file)
    p("You can click here to download the corrected DNA methylation values of all CpG sites")
  })
  output$download_corrected_CpGs <- renderUI({
    req(input$expected_file)
    downloadButton(outputId = "all_corrected_CpGs",
                   label = "Download corrected CpG sites")
  })
  output$all_corrected_CpGs <- downloadHandler(
    filename = "all_corrected_CpGs.tsv",
    content = function(con) {
      write.table(df_all_corrected_CpGs(), sep ="\t", row.names = F, con)
    }
  )
  ##### DNA methylation after correction #####
  df_after_plot <- reactive({
    df_after <- df_all_corrected_CpGs()[, grep("Amplicon|.after", colnames(df_all_corrected_CpGs()))]
    colnames(df_after) <- gsub(".after", "", colnames(df_after))
    df_after_plot <- reshape(df_after, idvar = "Amplicon", varying = grep("Amplicon", colnames(df_after), invert = T, value = T), 
                             v.names = "Methylation", times = grep("Amplicon", colnames(df_after), invert = T, value = T),
                             direction = "long")
    colnames(df_after_plot)[colnames(df_after_plot) == "time"] <- "Sample"
    
    df_after_plot
  })
  ##### By sample #####
  output$after_sample_select <- renderUI({
    req(df_after_plot())
    selectInput(inputId = "after_sample_select",
                label = "Select a sample to plot",
                choices = unique(df_after_plot()$Sample))
  })
  output$after_sample_plot <- renderPlotly({
    plot_ly(df_after_plot()[df_after_plot()$Sample == input$after_sample_select,], x = ~Amplicon, y = ~Methylation,
            type = "scatter", mode = "lines+markers", 
            hovertemplate = paste0("<b>Amplicon:</b> %{x}<br>",
                                   "<b>Methylation:</b> %{y}<br>",
                                   "<extra></extra>"),
            height = 700)%>%
      layout(title = paste0("Sample: ", input$after_sample_select),
             xaxis = list(title = ""),
             yaxis = list(range = c(0, 105)))
  })
  ##### By CpG site #####
  output$after_single_all <- renderUI({
    req(df_after_plot())
    radioButtons(inputId = "after_single_all", 
                 label = "", 
                 choices = c("All CpG sites" = 2, "Single CpG site" = 1))
  })
  output$after_CpG_select <- renderUI({
    req(input$after_single_all == 1)
    selectInput(inputId = "after_CpG_select", 
                label = "Select a CpG site", 
                choices = unique(df_after_plot()$Amplicon))
  })
  #Plot for methylation in 1 amplicon
  output$after_CpG_plot <- renderPlotly({
    req(input$after_single_all)
    if(input$after_single_all == 1) {
      df <- df_after_plot()[df_after_plot()$Amplicon == input$after_CpG_select,]
      plot_ly(df, x = ~Amplicon, y = ~Methylation, 
              type = "box", boxmean = T, hoverinfo = "y")%>%
        layout(title = paste0("Amplicon: ", input$after_CpG_select),
               xaxis = list(title = "", showticklabels = F))
    } else {
      df <- df_after_plot()
      plot_ly(df, x = ~Amplicon, y = ~Methylation, 
              type = "box", boxmean = T, color = ~Amplicon, height = 700)%>%
        layout(title = "DNA methylation of all CpGs",
               xaxis = list(title = ""),
               showlegend = F)
    }
  })
  ##### By group #####
  # Select an amplicon for plotting
  output$after_group_mode <- renderUI({
    req(dataset_group())
    radioButtons(inputId = "after_group_mode",
                 label = "Select a data visualization",
                 choices = c("MDS plot" = 2, "Box plot" = 1),
                 selected = 2)
  })
  output$after_group_select <- renderUI({
    req(input$file1)
    selectInput(inputId = "after_group_select",
                label = "Select a CpG site to plot DNA methylation of the groups of samples", 
                choices = unique(df_after_plot()$Amplicon), 
                multiple = F)
  })
  # Plotting DNA methylation for selected amplicon by group
  output$after_group_single <- renderPlotly({
    req(dataset_group())
    a <- df_after_plot()
    a$Group <- dataset_group()$Group[match(a$Sample, dataset_group()$Sample)]
    a2 <- a
    plot_ly(a2[a2$Amplicon==input$after_group_select,],x= ~Group ,y = ~Methylation, 
            type = "box", boxmean = T, color = ~Group,
            hoverinfo = "y")%>%
      layout(title = paste0("Amplicon: ", input$after_group_select),
             xaxis = list(title = ""),
             showlegend = F)
  })
  # Plotting DNA methylation by group
  output$after_group_sum <- renderPlotly({
    req(dataset_group())
    a <- df_after_plot()
    a$Group <- dataset_group()$Group[match(a$Sample, dataset_group()$Sample)]
    a1 <- a
    plot_ly(a1,x= ~Group ,y = ~Methylation, 
            type = "box", boxmean = T, color = ~Group,
            hoverinfo = "y")%>%
      layout(title = "All CpG sites: DNA methylation in sample groups",
             xaxis = list(title = ""),
             showlegend = F)
  })
  # Plotting DNA methylation for all amplicon by group
  output$after_group_all <- renderPlotly({
    req(dataset_group())
    a <- df_after_plot()
    a$Group <- dataset_group()$Group[match(a$Sample, dataset_group()$Sample)]
    a3 <- a
    a4 <- na.omit(a3)
    plot_ly(a4, x = ~Amplicon, y = ~Methylation, color = ~Group, 
            type = "box", boxmean = T) %>%
      layout(boxmode = "group",
             title = "Box plots for all CpG sites",
             xaxis = list(title = ""),
             showlegend = F)
  })
  ##### After t test 1 #####
  output$after_group_ttest <- renderPrint({
    req(input$after_group_select)
    a <- df_after_plot()[df_after_plot()$Amplicon == input$after_group_select,]
    a$Group <- dataset_group()$Group[match(a$Sample, dataset_group()$Sample)]
    a <- na.omit(a)
    a$idnew <- paste(a$Amplicon, a$Sample, sep = "_")
    if (length(unique(a$Group)) == 2) {
      t_test <- t.test(a$Methylation ~ a$Group, var.equal = F, alternative = "two.sided")
      cat(paste("Two-sided t test:",
                paste("p-value =", t_test$p.value),
                sep = "\n"))
    } else {
      subset <- a[,c("idnew", "Methylation", "Group")]
      df_wide <- reshape(subset, idvar = "idnew", timevar = "Group", direction = "wide")
      colnames(df_wide) <- gsub("Methylation.", "", colnames(df_wide))
      p.mat <- multi.ttest(df_wide[,-1], var.equal = F, alternative = "two.sided")
      upper <- p.mat
      upper[upper.tri(p.mat)] <- ""
      upper <- as.data.frame(upper)
      cat("Two-sided t test:", "\n")
      print(upper)
    }
  })
  ##### After t test 2 #####
  output$after_group_ttest2 <- renderPrint({
    req(dataset_group())
    a <- df_after_plot()
    a$Group <- dataset_group()$Group[match(a$Sample, dataset_group()$Sample)]
    a <- na.omit(a)
    a$idnew <- paste(a$Amplicon, a$Sample, sep = "_")
    if (length(unique(a$Group)) == 2) {
      t_test <- t.test(a$Methylation ~ a$Group, var.equal = F, alternative = "two.sided")
      cat(paste("Two-sided t test:",
                paste("p-value =", t_test$p.value),
                sep = "\n"))
    } else {
      subset <- a[,c("idnew", "Methylation", "Group")]
      df_wide <- reshape(subset, idvar = "idnew", timevar = "Group", direction = "wide")
      colnames(df_wide) <- gsub("Methylation.", "", colnames(df_wide))
      p.mat <- multi.ttest(df_wide[,-1], var.equal = F, alternative = "two.sided")
      upper <- p.mat
      upper[upper.tri(p.mat)] <- ""
      upper <- as.data.frame(upper)
      cat("Two-sided t test:", "\n")
      print(upper)
    }
  })
  ##### MDS plot after #####
  df_mds_group_after <- reactive({
    df <- df_all_corrected_CpGs()
    df1 <- df[, grep("after", colnames(df))]
    dm <- 1 - abs(cor(df1, use = "complete.obs"))
    fit <- cmdscale(dm, eig = TRUE, k = 2)
    df_plot <- as.data.frame(fit$points)
    df_plot$Sample <- row.names(df_plot)
    df_plot$Sample <- gsub(".after", "", df_plot$Sample)
    df_mds_merge <- merge(df_plot, dataset_group(), by = "Sample")
    df_mds_merge
  })
  output$mds_group_after <- renderPlotly({
    plot_ly(df_mds_group_after(), x = ~V1, y = ~V2, text = ~Sample, color = ~Group,
            type = "scatter", mode = "markers", 
            marker = list(size = 12), 
            hovertemplate = paste0("<b>Sample:</b> %{text}<br>",
                                   "<b>Group:</b> %{fullData.name}<br>",
                                   "<extra></extra>"),
            height = 600)%>%
      add_text(textposition = "top right", showlegend = F)%>%
      layout(title = "MDS plot",
             xaxis = list(title = "PC1", 
                          zeroline = F, range = c(min(df_mds_group_after()$V1) - 0.2,
                                                  max(df_mds_group_after()$V1) + 0.2)),
             yaxis = list(title = "PC2", zeroline = F),
             showlegend = T)
  })
}

shinyApp(ui = ui, server = server)
