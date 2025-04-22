# app.R
library(shiny)
library(readxl)
library(openxlsx)
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(ggpubr)

# 加载模块
source("qpcr_module.R")
source("protein_module.R")
source("primer_module.R")

# 加载用户数据库
users_db <- read.csv("~/R/RNA/data/NASH/网页客户端/users.csv", stringsAsFactors = FALSE)

# UI
ui <- fluidPage(
  tags$head(
    # 加载动画控制器（全局）
    tags$script(HTML("
      Shiny.addCustomMessageHandler('show_loader', function(id) {
        document.getElementById(id).style.display = 'block';
      });
      Shiny.addCustomMessageHandler('hide_loader', function(id) {
        document.getElementById(id).style.display = 'none';
      });
    "))
  ),
  
  uiOutput("page_ui")  # 原有主UI渲染点
)

# Server
server <- function(input, output, session) {
  user_status <- reactiveValues(logged_in = FALSE, current_page = NULL, username = NULL)
  
  
  write_log <- function(user, action) {
    log_file <- "~/R/Autoanalyse/Autorun/user_log.csv"
    time <- Sys.time()
    log_entry <- data.frame(Time = time, User = user, Action = action)
    if (!file.exists(log_file)) {
      write.csv(log_entry, log_file, row.names = FALSE)
    } else {
      write.table(log_entry, log_file, row.names = FALSE, col.names = FALSE,
                  sep = ",", append = TRUE)
    }
  }
  
  # 登录页面 UI
  login_ui <- fluidPage(
    tags$head(tags$style(HTML("
    .center-box {
      display: flex;
      justify-content: center;
      align-items: center;
      height: 80vh;
    }
    .login-panel {
      width: 400px;
    }
  "))),
    div(class = "center-box",
        div(class = "login-panel",
            wellPanel(
              h3("用户登录", align = "center"),
              textInput("username", "用户名", value = "test"),
              passwordInput("password", "密码", value = "1234"),
              actionButton("login_btn", "登录", class = "btn btn-primary", width = "100%")
            )
        )
    )
  )
  
  # 模块选择 UI
  select_ui <- fluidPage(
    tags$head(tags$style(HTML("
    .select-center {
      display: flex;
      justify-content: center;
      align-items: center;
      height: 80vh;
      flex-direction: column;
    }
    .select-buttons {
      width: auto;
      display: flex;
      gap: 30px;
      justify-content: center;
      margin-bottom: 40px;
    }
    .module-button {
      width: 220px;
      height: 60px;
      font-size: 18px;
    }
    .select-title {
      margin-bottom: 60px;
    }
  "))),
    
    div(class = "select-center",
        div(class = "select-title",
            h2("模块选择")
        ),
        div(class = "select-buttons",
            actionButton("goto_qpcr", "进入 RT-qPCR 模块", class = "btn btn-primary module-button"),
            actionButton("goto_protein", "进入蛋白分析模块", class = "btn btn-success module-button"),
            actionButton("goto_primer", "引物验证", class = "btn-warning module-button")
        ),
        actionButton("logout", "退出登录", class = "btn btn-danger module-button")
    )
  )
  
  # 页面渲染
  output$page_ui <- renderUI({
    if (!user_status$logged_in) {
      login_ui
    } else if (is.null(user_status$current_page)) {
      select_ui
    } else if (user_status$current_page == "qpcr") {
      qpcr_ui("qpcr")
    } else if (user_status$current_page == "protein") {
      protein_ui("protein")
    } else if (user_status$current_page == "primer") {
      primer_ui("primer")
    }
  })
  
  # 登录逻辑
  observeEvent(input$login_btn, {
    user_row <- users_db[users_db$username == input$username & users_db$password == input$password, ]
    if (nrow(user_row) == 1) {
      user_status$logged_in <- TRUE
      user_status$username <- input$username
      write_log(input$username, "登录成功")
    } else {
      showModal(modalDialog("用户名或密码错误", easyClose = TRUE))
      write_log(input$username, "登录失败")
    }
  })
  
  # 页面跳转逻辑
  observeEvent(input$goto_qpcr, { 
    user_status$current_page <- "qpcr"
    write_log(user_status$username, "进入 RT-qPCR 模块")
  })
  observeEvent(input$goto_protein, { 
    user_status$current_page <- "protein"
    write_log(user_status$username, "进入蛋白分析模块")
  })
  observeEvent(input$goto_primer, { 
    user_status$current_page <- "primer"
    write_log(user_status$username, "进入引物验证模块")
  })
  observeEvent(input$logout, {
    user_status$logged_in <- FALSE
    user_status$current_page <- NULL
  })
  
  # 各模块服务器逻辑
  qpcr_server("qpcr", user_status)
  protein_server("protein", user_status)
  primer_server("primer", user_status)
}

shinyApp(ui, server)