qpcr_ui <- function(id) {
  ns <- NS(id)
  fluidPage(
    titlePanel("RT‑qPCR 一键分析"),
    sidebarLayout(
      sidebarPanel(
        width = 3,
        fileInput(ns("layout"), "上传通用加样表", accept = c(".xls", ".xlsx")),
        downloadButton(ns("download_example"), "下载范例加样表", class = "btn-info", style = "margin-top: -15px; margin-bottom: 30px;"),
        fileInput(ns("raw"),    "上传 Raw Results", accept = c(".xls", ".xlsx")),
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
    
    output$download_example <- downloadHandler(
      filename = function() {
        "example_layout.xlsx"
      },
      content = function(file) {
        file.copy("example_layout.xlsx", file)
      }
    )
    
    layout_df <- reactive({
      req(input$layout)
      Samp <- readxl::read_excel(input$layout$datapath, range = "A2:L9", col_names = FALSE) %>%
        pivot_longer(everything(), values_to = "Sample")
      Gene <- readxl::read_excel(input$layout$datapath, range = "M2:X9", col_names = FALSE) %>%
        pivot_longer(everything(), values_to = "Gene")
      tibble(Sample = Samp$Sample, Gene = Gene$Gene)
    })
    
    res <- eventReactive(input$run, {
      req(input$raw, layout_df())
      raw <- readxl::read_excel(input$raw$datapath, sheet = "Results", range = "A48:Q144", col_names = TRUE)
      flat <- layout_df()
      raw <- raw %>% mutate(`Sample Name` = flat$Sample, `Target Name` = flat$Gene)
      colnames(raw)[colnames(raw) == "Sample Name"] <- "Sample"
      colnames(raw)[colnames(raw) == "Target Name"] <- "Detector"
      colnames(raw)[colnames(raw) == "CT"] <- "Ct"
      raw <- raw %>% filter(Task != "") %>% group_by(Sample, Detector) %>% mutate(replicate = row_number())
      raw$Ct <- as.numeric(raw$Ct)
      raw$Detector <- stringr::str_replace_all(raw$Detector, regex("^actin$", ignore_case = TRUE), "actin")
      
      actin_mean <- raw %>% filter(Detector == "actin") %>%
        group_by(Sample) %>% summarise(actin_mean = mean(Ct, na.rm = TRUE))
      
      samples <- unique(raw$Sample)
      detectors <- setdiff(unique(raw$Detector), "actin")
      ref_sample <- if (input$ref != "") input$ref else samples[1]
      
      out <- list()
      for (s in samples) for (d in detectors) {
        ct <- raw %>% filter(Sample == s, Detector == d) %>% pull(Ct)
        if (all(is.na(ct))) next
        delta <- ct - actin_mean$actin_mean[actin_mean$Sample == s]
        ct_ref <- raw %>% filter(Sample == ref_sample, Detector == d) %>% pull(Ct)
        if (all(is.na(ct_ref))) next
        delta_ref <- ct_ref - actin_mean$actin_mean[actin_mean$Sample == ref_sample]
        out[[length(out) + 1]] <- data.frame(Sample = s, Detector = d, ddCt = delta - mean(delta_ref, na.rm = TRUE))
      }
      
      bind_rows(out) %>%
        mutate(`2^-ddCt` = 2^-ddCt) %>%
        group_by(Sample, Detector) %>%
        mutate(replicate = row_number()) %>%
        ungroup()
    })
    
    output$plt <- renderPlot({
      df <- res()
      sample_order <- layout_df()$Sample %>% unique()
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
        df <- res() %>% filter(!is.na(`2^-ddCt`))
        sample_order <- layout_df()$Sample %>% unique()
        df <- df %>%
          mutate(Sample = factor(Sample, levels = sample_order)) %>%
          arrange(Detector, Sample, replicate) %>%
          mutate(SampleRep = paste0(Sample, "_rep", replicate))
        prism_table <- df %>%
          select(Detector, SampleRep, `2^-ddCt`) %>%
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
