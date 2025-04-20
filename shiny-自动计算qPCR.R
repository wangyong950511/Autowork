# ---- 必备包 ----
library(shiny)
library(readxl)
library(openxlsx)
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(ggpubr)

ui <- fluidPage(
  tags$head(
    tags$style(HTML("
        .container-fluid {
          max-width: 1600px;
          margin-left: auto;
          margin-right: auto;
        }
      
        /* 控制 sidebar 和 mainpanel 宽度比例 */
        .row > .col-sm-3 {
          flex: 0 0 360px;
          max-width: 360px;
        }
      
        .row > .col-sm-9 {
          flex: 1;
          max-width: calc(100% - 360px);
        }
      
        .main-panel {
          min-height: 600px;
        }
      "))
  ),
  
  titlePanel("RT‑qPCR 一键分析"),
  
  sidebarLayout(
    sidebarPanel(
      width = 3,
      fileInput("layout", "上传通用加样表", accept = c(".xls", ".xlsx")),
      fileInput("raw",    "上传 Raw Results", accept = c(".xls", ".xlsx")),
      textInput("ref", "参考样本 (默认第一个孔)", value = ""),
      actionButton("run", "开始分析"),
      br(), br(),
      downloadButton("dl_excel", "下载 Prism 分析表格")
    ),
    mainPanel(class = "main-panel",
              plotOutput("plt", width = "1200px", height = "500px")
    )
  )
)

server <- function(input, output, session){
  
  layout_df <- reactive({
    req(input$layout)
    Samp <- read_excel(input$layout$datapath, range = "A2:L9", col_names = FALSE) %>%
      pivot_longer(everything(), values_to = "Sample")
    Gene <- read_excel(input$layout$datapath, range = "M2:X9", col_names = FALSE) %>%
      pivot_longer(everything(), values_to = "Gene")
    tibble(Sample = Samp$Sample, Gene = Gene$Gene)
  })
  
  res <- eventReactive(input$run, {
    req(input$raw, layout_df())
    
    raw <- read_excel(input$raw$datapath, sheet = "Results",
                      range = "A48:Q144", col_names = TRUE)
    
    flat <- layout_df()
    raw <- raw %>% mutate(`Sample Name` = flat$Sample,
                          `Target Name` = flat$Gene)
    
    names(raw)[names(raw)=="Sample Name"] <- "Sample"
    names(raw)[names(raw)=="Target Name"] <- "Detector"
    names(raw)[names(raw)=="CT"]          <- "Ct"
    
    raw <- raw %>%
      filter(Task != "") %>%
      group_by(Sample, Detector) %>%
      mutate(replicate = row_number())
    
    raw$Ct <- as.numeric(raw$Ct)
    raw$Detector <- str_replace_all(raw$Detector, regex("^actin$", ignore_case = TRUE), "actin")
    
    actin_mean <- raw %>%
      filter(Detector=="actin") %>%
      group_by(Sample) %>%
      summarise(actin_mean = mean(Ct, na.rm=TRUE))
    
    samples   <- unique(raw$Sample)
    detectors <- setdiff(unique(raw$Detector), "actin")
    ref_sample <- if(input$ref != "") input$ref else samples[1]
    
    out <- list()
    for(s in samples) for(d in detectors){
      ct <- raw %>% filter(Sample==s, Detector==d) %>% pull(Ct)
      if(all(is.na(ct))) next
      delta <- ct - actin_mean$actin_mean[actin_mean$Sample==s]
      
      ct_ref <- raw %>% filter(Sample==ref_sample, Detector==d) %>% pull(Ct)
      if(all(is.na(ct_ref))) next
      delta_ref <- ct_ref - actin_mean$actin_mean[actin_mean$Sample==ref_sample]
      
      out[[length(out)+1]] <- data.frame(Sample=s, Detector=d,
                                         ddCt = delta - mean(delta_ref, na.rm=TRUE),
                                         stringsAsFactors = FALSE)
    }
    
    bind_rows(out) %>%
      mutate(`2^-ddCt` = 2^-ddCt) %>%
      group_by(Sample, Detector) %>%
      mutate(replicate = row_number()) %>%
      ungroup()
  })
  
  output$plt <- renderPlot({
    df <- res()
    
    # 提取样本顺序
    sample_order <- read_excel(input$layout$datapath, range = "A2:L9", col_names = FALSE) %>%
      as.matrix() %>%
      t() %>%  # 转置，确保从左到右、从上到下
      as.vector() %>%
      .[!is.na(.)] %>%
      unique()
    
    df_summary <- df %>%
      group_by(Sample, Detector) %>%
      summarise(mean = mean(`2^-ddCt`), sd = sd(`2^-ddCt`), .groups = "drop") %>%
      mutate(Sample = factor(Sample, levels = sample_order))
    
    ggplot(df_summary, aes(Sample, mean, fill = Detector))+
      geom_bar(stat="identity", position=position_dodge(0.8), width=0.7, color="black")+
      geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=0.25,
                    position=position_dodge(0.8))+
      facet_wrap(~Detector, scales="free_y")+
      theme_pubr(base_size = 14)+
      labs(y="2^-ΔΔCt", x = "Sample")
  })
  
  output$dl_excel <- downloadHandler(
    filename = function(){ paste0("qPCR_prism_", Sys.Date(), ".xlsx") },
    content  = function(file){
      df <- res()
      df <- df %>% filter(!is.na(`2^-ddCt`))
      
      # --- 关键：从加样表中读取样本顺序 ---
      sample_order <- read_excel(input$layout$datapath, range = "A2:L9", col_names = FALSE) %>%
        as.matrix() %>%
        t() %>%
        as.vector() %>%
        .[!is.na(.)] %>%
        unique()
      
      # 应用样本排序
      df <- df %>%
        mutate(Sample = factor(Sample, levels = sample_order)) %>%
        arrange(Detector, Sample, replicate) %>%
        mutate(SampleRep = paste0(Sample, "_rep", replicate))
      
      # Prism 格式：每个基因为一行，每个 SampleRep 为列
      prism_table <- df %>%
        select(Detector, SampleRep, `2^-ddCt`) %>%
        pivot_wider(names_from = SampleRep, values_from = `2^-ddCt`) %>%
        arrange(Detector)
      
      # 写入 Excel
      wb <- createWorkbook()
      addWorksheet(wb, "Prism_Format")
      writeData(wb, "Prism_Format", prism_table)
      saveWorkbook(wb, file, overwrite = TRUE)
    }
  )
}

shinyApp(ui, server)
