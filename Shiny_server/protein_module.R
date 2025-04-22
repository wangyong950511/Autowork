protein_ui <- function(id) {
  ns <- NS(id)
  fluidPage(
    titlePanel("蛋白定量分析工具"),
    sidebarLayout(
      sidebarPanel(
        width = 3,
        fileInput(ns("file"), "上传Excel文件 (包含Result sheet)", accept = ".xlsx"),
        numericInput(ns("volume"), "总加样体积 (µL)", value = 200),
        checkboxInput(ns("new_std"), "重新计算标准曲线", value = TRUE),
        conditionalPanel(
          condition = paste0("!input['", ns("new_std"), "']"),
          numericInput(ns("x_factor"), "既往标准曲线系数 (斜率)", value = 28.192847),
          numericInput(ns("constant"), "既往标准曲线常数 (截距)", value = -4.170041)
        ),
        numericInput(ns("n_samples"), "样本数", value = 9, min = 1),
        uiOutput(ns("sample_inputs")),
        actionButton(ns("analyze"), "开始分析"),
        br(), br(),
        actionButton(ns("back_home2"), "返回模块选择", class = "btn-secondary")
      ),
      mainPanel(
        width = 9,
        fluidRow(
          column(6, plotOutput(ns("stdcurve_plot"), width = "600px", height = "500px")),
          column(6, plotOutput(ns("bar_plot"), width = "600px", height = "500px"))
        ),
        fluidRow(
          column(6, tableOutput(ns("concentration_table"))),  # 默认受 column(6) 控制宽度
          column(6, plotOutput(ns("volume_plot"), width = "600px", height = "500px"))
        ),
        verbatimTextOutput(ns("summary_text"))
      )
    )
  )
}

protein_server <- function(id, user_status) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    output$sample_inputs <- renderUI({
      n <- input$n_samples
      cols <- 3
      rows <- ceiling(n / cols)
      padded <- c(1:n, rep(NA, rows * cols - n))
      matrix_inputs <- matrix(padded, nrow = rows, ncol = cols, byrow = TRUE)
      
      tagList(
        lapply(1:rows, function(i) {
          fluidRow(
            lapply(1:cols, function(j) {
              idx <- matrix_inputs[i, j]
              if (!is.na(idx)) {
                column(4, textInput(ns(paste0("sample_", idx)), paste0(idx, "."), value = as.character(idx)))
              }
            })
          )
        })
      )
    })
    
    observeEvent(input$analyze, {
      req(input$file)
      
      group_labels <- sapply(1:input$n_samples, function(i) {
        input[[paste0("sample_", i)]]
      })
      
      data <- readxl::read_excel(input$file$datapath, sheet = "Result sheet", col_names = TRUE, range = "B43:M51")
      
      if (input$new_std) {
        x <- as.numeric(data[1, 1:8])
        y <- c(0, 0.05, 0.1, 0.2, 0.4, 0.6, 0.8, 1)
        fit <- lm(y ~ x)
        r_squared <- summary(fit)$r.squared
        coef_values <- coef(fit)
        constant <- as.numeric(coef_values[1])
        x_factor <- as.numeric(coef_values[2])
        data <- data[-1, ]
      } else {
        constant <- input$constant
        x_factor <- input$x_factor
        r_squared <- NA
        fit <- NULL
      }
      
      group_size <- 3
      num_groups <- ncol(data) / group_size
      data_split <- lapply(1:num_groups, function(i) {
        cols <- ((i - 1) * group_size + 1):(i * group_size)
        tmp <- data[, cols]
        colnames(tmp) <- c("1", "2", "3")
        tmp
      })
      
      combined_data <- do.call(rbind, data_split)
      data_clean <- combined_data[complete.cases(combined_data), ]
      compute_new_value <- function(x) (x_factor * x + constant) * 10
      new_data <- as.data.frame(apply(data_clean, 2, compute_new_value))
      
      row_means <- apply(new_data, 1, mean)
      row_sd <- apply(new_data, 1, sd)
      summary_data <- data.frame(new_data, Mean = row_means, SD = row_sd)
      summary_data$Sample <- as.factor(1:nrow(summary_data))
      
      bad_row <- which(summary_data$SD > summary_data$Mean * 0.3)
      output_text1 <- if (length(bad_row) > 0) {
        paste(paste(bad_row, collapse = ","), "Need to be redone")
      } else {
        "All samples meet the requirements"
      }
      
      min_value <- min(summary_data$Mean)
      summary_data$Protein_volum <- as.integer((input$volume * min_value) / summary_data$Mean)
      summary_data$RIPA_volum <- input$volume - summary_data$Protein_volum
      
      out_data <- summary_data[, c("Mean", "RIPA_volum", "Protein_volum")]
      out_data$group <- as.factor(1:nrow(out_data))
      data_long <- tidyr::gather(out_data, key = "variable", value = "value", Protein_volum, RIPA_volum)
      data_long$variable <- factor(data_long$variable, levels = c("RIPA_volum", "Protein_volum"))
      data_long$group <- factor(data_long$group, levels = rev(levels(data_long$group)))
      group_labels_rev <- rev(group_labels)
      output_text2 <- paste("Protein concentration=", round(min_value, 1), "ug/ul, Protein volume=", input$volume, "µL")
      title_text <- paste("Add Rules -", output_text1, "\n", output_text2)
      
      output$stdcurve_plot <- renderPlot({
        if (!is.null(fit)) {
          plot(x, y, main = "Standard Curve", xlab = "Absorbance", ylab = "C")
          abline(fit, col = "red")
          text(x = min(x), y = max(y), labels = paste("R² =", round(r_squared, 4)), pos = 4, col = "blue")
        }
      })
      
      output$bar_plot <- renderPlot({
        ggplot(summary_data, aes(x = Sample, y = Mean)) +
          geom_bar(stat = "identity", fill = "blue", alpha = 0.5) +
          geom_errorbar(aes(ymin = Mean - SD, ymax = Mean + SD), width = 0.2, color = "black") +
          labs(title = "Mean and SD Bar Plot", x = "Sample", y = "Value") +
          theme_minimal() +
          theme(axis.text.x = element_text(angle = 45, hjust = 1))
      })
      
      output$volume_plot <- renderPlot({
        ggplot(data_long, aes(x = value, y = group, fill = variable)) +
          geom_bar(stat = "identity", position = "stack", width = 0.7) +
          geom_text(aes(label = value), position = position_stack(vjust = 0.5), size = 6, color = "black") +
          labs(title = title_text, x = "Volum", y = "Sample") +
          scale_fill_manual(values = c("Protein_volum" = "red", "RIPA_volum" = "yellow")) +
          scale_y_discrete(labels = group_labels_rev) +
          theme_minimal() +
          theme(axis.text.x = element_text(angle = 45, hjust = 1))
      })
      
      output$concentration_table <- renderTable({
        summary_data %>%
          select(Sample, Mean, SD) %>%
          mutate(across(c(Mean, SD), round, 2)) %>%
          rename(`样本` = Sample, `浓度(μg/μL)` = Mean, `标准差` = SD)
      })
      
      output$summary_text <- renderPrint({
        cat(output_text1, "\n", output_text2)
      })
    })
    
    observeEvent(input$back_home2, {
      user_status$current_page <- NULL
    })
  })
}