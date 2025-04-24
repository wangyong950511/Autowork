library(shiny)
library(Biostrings)
library(xml2)
library(dplyr)



primer_ui <- function(id) {
  ns <- NS(id)
  fluidPage(
    tags$head(tags$style(HTML("
      .loader {
        border: 8px solid #f3f3f3;
        border-top: 8px solid #3498db;
        border-radius: 50%;
        width: 40px;
        height: 40px;
        animation: spin 1s linear infinite;
        margin: 20px auto;
        display: none;
      }
      @keyframes spin {
        0% { transform: rotate(0deg); }
        100% { transform: rotate(360deg); }
      }
    "))),
    
    titlePanel("引物特异性验证（BLAST）"),
    
    sidebarLayout(
      sidebarPanel(
        width = 4,
        textInput(ns("forward_primer"), "正向引物序列", value = "ACCGGCCTCTTCTTCGTCT"),
        textInput(ns("reverse_primer"), "反向引物序列", value = "AACTGCCTGTGTTGTCGATCT"),
        actionButton(ns("run_blast"), "运行 BLAST", class = "btn-primary"),
        br(), br(),
        div(id = ns("loading"), class = "loader"),
        actionButton(ns("back_home3"), "返回模块选择", class = "btn-secondary")
      ),
      mainPanel(
        fluidRow(
          column(6, tableOutput(ns("result_table"))),
          column(6, tableOutput(ns("unique_genes")))
        )
      )
    )
  )
}



primer_server <- function(id, user_status) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    # 加载图标控制
    observeEvent(input$run_blast, {
      session$sendCustomMessage("show_loader", ns("loading"))
    })
    
    observeEvent(input$run_blast, {
      req(input$forward_primer, input$reverse_primer)
      
      # 工作目录和环境配置
      workdir <- "/home/drwang/R/Autoanalyse/Premier/workdata"
      blastdb_path <- "/home/drwang/Documents/refseq_rna"
      blast_exe <- "/home/drwang/Downloads/ncbi-blast-2.15.0+/bin"
      Sys.setenv(PATH = paste(Sys.getenv("PATH"), blast_exe, sep = ":"))
      Sys.setenv(BLASTDB = blastdb_path)
      setwd(workdir)
      
      session$sendCustomMessage("show_loader", ns("loading"))   # 显示加载图标
      
      # 写入 FASTA
      writeLines(paste0(">forward\n", input$forward_primer), "forward.fasta")
      writeLines(paste0(">reverse\n", input$reverse_primer), "reverse.fasta")
      
      # 运行命令
      fwd_cmd <- "blastn -query forward.fasta -db refseq_rna -out forward.xml -outfmt 5 -task blastn-short -num_threads 19"
      rev_cmd <- "blastn -query reverse.fasta -db refseq_rna -out reverse.xml -outfmt 5 -task blastn-short -num_threads 19"
      fwd_msg <- system(fwd_cmd, intern = TRUE, ignore.stderr = FALSE)
      rev_msg <- system(rev_cmd, intern = TRUE, ignore.stderr = FALSE)
      
      # 日志显示
      output$blast_log <- renderPrint({
        cat(">>> 正向引物命令输出:\n", paste(fwd_msg, collapse = "\n"), "\n\n")
        cat(">>> 反向引物命令输出:\n", paste(rev_msg, collapse = "\n"))
      })
      
      # 解析函数
      parse_blast_results <- function(xml_file) {
        xml_data <- read_xml(xml_file)
        hits <- xml_find_all(xml_data, "//Hit")
        tibble(
          HitID = sapply(hits, function(hit) sub(".*ref\\|([^|]+)\\|.*", "\\1", xml_text(xml_find_first(hit, ".//Hit_id")))),
          GeneName = sapply(hits, function(hit) sub(".*\\(([^)]+)\\).*", "\\1", xml_text(xml_find_first(hit, ".//Hit_def"))))
        )
      }
      
      # 如果 BLAST 成功运行则继续解析
      if (file.exists("forward.xml") && file.exists("reverse.xml")) {
        fwd_df <- parse_blast_results("forward.xml")
        rev_df <- parse_blast_results("reverse.xml")
        result <- inner_join(fwd_df, rev_df, by = "HitID", suffix = c("_fwd", "_rev")) %>%
          select(HitID, GeneName = GeneName_rev) %>%
          distinct()
        
        output$unique_genes <- renderTable({
          unique_genes <- unique(result$GeneName)
          data.frame(序号 = seq_along(unique_genes), 基因名 = unique_genes)
        })
        
        output$result_table <- renderTable({
          result
        })
      } else {
        output$result_table <- renderTable({
          data.frame(错误 = "BLAST XML 文件不存在，请检查数据库路径或权限")
        })
      }
      
      session$sendCustomMessage("hide_loader", ns("loading"))   # 隐藏加载图标
      
    })
    
    observeEvent(input$back_home3, {
      user_status$current_page <- NULL
    })
  })
}