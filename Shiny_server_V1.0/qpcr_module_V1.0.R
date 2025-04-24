qpcr_ui <- function(id) {
  ns <- NS(id)
  fluidPage(
    titlePanel("RT‑qPCR 一键分析"),
    sidebarLayout(
      sidebarPanel(
        width = 3,
        fileInput(ns("raw"),    "上传 Raw Results", accept = c(".xls", ".xlsx")),
        fileInput(ns("layout"), "上传通用加样表", accept = c(".xls", ".xlsx")),
        downloadButton(ns("download_example"), "下载范例加样表", class = "btn-info", style = "margin-top: -15px; margin-bottom: 30px;"),
        textInput(ns("ref"), "参考样本 (默认第一个孔)", value = ""),
        actionButton(ns("run"), "开始分析"),
        br(), br(),
        downloadButton(ns("dl_excel"), "下载 Prism 分析表格"),
        br(), br(),
        actionButton(ns("back_home"), "返回模块选择", class = "btn-secondary")
      ),
      mainPanel(
        plotOutput(ns("plt"), width = "1200px", height = "500px")
      )
    )
  )
}

qpcr_server <- function(id, user_status) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    output$download_example <- downloadHandler(
      filename = function() "example_layout.xlsx",
      content = function(file) file.copy("example_layout.xlsx", file)
    )
    
    layout_df <- reactive({
      req(input$layout)
      Sample_data <- readxl::read_excel(input$layout$datapath, range = "A2:L9", col_names = FALSE)
      Gene_data <- readxl::read_excel(input$layout$datapath, range = "M2:X9", col_names = FALSE)
      
      Sample_long <- Sample_data %>%
        pivot_longer(cols = everything(), names_to = "col", values_to = "Sample") %>%
        mutate(Well = as.character(row_number()))
      
      Gene_long <- Gene_data %>%
        pivot_longer(cols = everything(), names_to = "col", values_to = "Gene") %>%
        mutate(Well = as.character(row_number()))
      
      sample_order <- unique(Sample_long$Sample) %>% na.omit()
      list(Sample_long = Sample_long, Gene_long = Gene_long, sample_order = sample_order)
    })
    
    result_store <- reactiveVal(NULL)
    
    observeEvent(input$run, {
      req(input$raw, layout_df())
      data <- readxl::read_excel(input$raw$datapath, sheet = "Results", range = "A48:Q144", col_names = TRUE) %>%
        mutate(Well = as.character(Well))
      
      layout <- layout_df()
      Gene_long <- layout$Gene_long
      Sample_long <- layout$Sample_long
      sample_order <- layout$sample_order
      
      data <- data %>%
        left_join(Gene_long, by = "Well") %>%
        left_join(Sample_long, by = "Well") %>%
        mutate(`Target Name` = Gene, `Sample Name` = Sample) %>%
        dplyr::select(-Gene, -Sample)
      
      refsample <- if (input$ref != "") input$ref else na.omit(Sample_long$Sample)[1]
      
      colnames(data)[colnames(data) == "Sample Name"] <- "Sample"
      colnames(data)[colnames(data) == "Target Name"] <- "Detector"
      colnames(data)[colnames(data) == "CT"] <- "Ct"
      
      data <- data %>%
        filter(Task != "", !is.na(Ct), Ct != "Undetermined") %>%
        mutate(Ct = as.numeric(Ct)) %>%
        group_by(Sample, Detector) %>%
        mutate(replicate = row_number()) %>%
        ungroup()
      
      data$Detector <- tolower(data$Detector)
      data$Detector[data$Detector %in% c("actin", "actb")] <- "actin"
      
      actin_mean <- data %>%
        filter(Detector == "actin") %>%
        group_by(Sample) %>% summarise(actin_mean = mean(Ct, na.rm = TRUE), .groups = "drop")
      
      samples <- unique(data$Sample)
      detectors <- setdiff(unique(data$Detector), "actin")
      
      result_list <- list()
      
      for (s in samples) for (d in detectors) {
        ct <- data %>% filter(Sample == s, Detector == d) %>% pull(Ct)
        if (all(is.na(ct))) next
        delta <- ct - actin_mean$actin_mean[actin_mean$Sample == s]
        ct_ref <- data %>% filter(Sample == refsample, Detector == d) %>% pull(Ct)
        if (all(is.na(ct_ref))) next
        delta_ref <- ct_ref - actin_mean$actin_mean[actin_mean$Sample == refsample]
        ddct <- delta - mean(delta_ref, na.rm = TRUE)
        result_list[[length(result_list) + 1]] <- data.frame(Sample = s, Detector = d, ddCt = ddct)
      }
      
      result_store(bind_rows(result_list) %>% mutate(`2^-ddCt` = 2^-ddCt) %>%
                     group_by(Sample, Detector) %>% mutate(replicate = row_number()) %>% ungroup())
    })
    
    output$plt <- renderPlot({
      df <- result_store()
      req(df)
      sample_order <- layout_df()$sample_order
      
      df_summary <- df %>%
        group_by(Sample, Detector) %>%
        summarise(mean = mean(`2^-ddCt`), sd = sd(`2^-ddCt`), .groups = "drop") %>%
        mutate(Sample = factor(Sample, levels = sample_order))
      
      ggplot(df_summary, aes(Sample, mean, fill = Detector)) +
        geom_bar(stat = "identity", position = position_dodge(0.8), width = 0.7, color = "black") +
        geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), width = 0.25, position = position_dodge(0.8)) +
        facet_wrap(~Detector, scales = "free_y") +
        theme_pubr(base_size = 14) +
        labs(y = "2^-ΔΔCt", x = "Sample")
    })
    
    output$dl_excel <- downloadHandler(
      filename = function() paste0("qPCR_prism_", Sys.Date(), ".xlsx"),
      content = function(file) {
        df <- result_store() %>% filter(!is.na(`2^-ddCt`))
        sample_order <- layout_df()$sample_order
        df <- df %>%
          mutate(Sample = factor(Sample, levels = sample_order)) %>%
          arrange(Detector, Sample, replicate) %>%
          mutate(SampleRep = paste0(Sample, "_rep", replicate))
        
        prism_table <- df %>%
          dplyr::select(Detector, SampleRep, `2^-ddCt`) %>%
          pivot_wider(names_from = SampleRep, values_from = `2^-ddCt`) %>%
          arrange(Detector)
        
        wb <- openxlsx::createWorkbook()
        openxlsx::addWorksheet(wb, "Prism_Format")
        openxlsx::writeData(wb, "Prism_Format", prism_table)
        openxlsx::saveWorkbook(wb, file, overwrite = TRUE)
      }
    )
    
    observeEvent(input$back_home, {
      user_status$current_page <- NULL
    })
  })
}