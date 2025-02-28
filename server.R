library(shiny)
library(tidyverse)
library(broom)
library(car)
library(rhandsontable)
library(bslib)
library(ggplot2)
library(rsconnect)
library(shinythemes)
library(reshape2)
library(DT)
library(purrr)
library(EnvStats)
library(performance)
library(scales)




statistics <-  function(data, x, y) {
  data %>%
    reframe(
      model = list(lm(!!sym(y) ~ !!sym(x), data = cur_data())),
      Coef_deter = list(data.frame(R2 = performance::r2(model[[1]])$R2)),
      Regression = list(model[[1]] %>% tidy()), # Regression model
      Levene = list(leveneTest(!!sym(y) ~ as.factor(!!sym(x))) %>% tidy()), # Homogeneity test
      Durbin_Watson = list(durbinWatsonTest(model[[1]], max.lag = 1) %>% tidy()), # Autocorrelation test
      Shapiro_Wilk = list(shapiro.test(model[[1]]$residuals) %>% tidy()), # Normality test
      Anova = list(anovaPE(model[[1]]) %>% tidy()) # ANOVA test
    ) %>%
    select(-model)
}

limits <- function(data, x, y) {
  data %>%
    reframe(
      reg = list(lm(!!sym(y) ~ !!sym(x), data = cur_data())),
      anova_result = list(anova(reg[[1]])),
      dp = sqrt(anova_result[[1]]$`Mean Sq`[2]),
      b = reg[[1]]$coefficients[2], # Sensitivity
      # LOD
      LOD = 3.3 * dp / b,
      # LOQ
      LOQ = 10 * dp / b
    ) %>% 
    select(LOD,LOQ)
}


calculate_coefficients_compare <- function(data,  x, y) {
  data %>%
    reframe(
      model = list(lm(!!sym(y) ~ !!sym(x), data = cur_data())), # Perform regression
      model_tidy = list(tidy(model[[1]])),
      df = df.residual(model[[1]]), # Set degrees of freedom 
      slope = model_tidy[[1]]$estimate[2], # Get slope and format
      slope.std.error =  model_tidy[[1]]$std.error[2], # Get slope error and format
      t = qt(1 - 0.025, df), # Calculate t critical value 
      CI.slope = qt(1 - 0.025, df) * as.numeric(slope.std.error), # Calculate CI and format
      upr.slope = as.numeric(slope) + as.numeric(CI.slope), # Upper limit and format
      lwr.slope = as.numeric(slope) - as.numeric(CI.slope), # Lower limit and format
      intercept = model_tidy[[1]]$estimate[1], # Format intercept
      intercept.std.error = model_tidy[[1]]$std.error[1], # Format intercept error
      CI.intercept = qt(1 - 0.025, df) * as.numeric(intercept.std.error), # Calculate CI and format
      upr.intercept = as.numeric(intercept) + as.numeric(CI.intercept), # Upper limit and format
      lwr.intercept = as.numeric(intercept) - as.numeric(CI.intercept) # Lower limit and format
    )  %>%
    select(-model, -model_tidy, -t, -df, -slope.std.error, -intercept.std.error) # Remove undesired vectors
}

calculate_found <- function(data) {
  if(nrow(data) < 2) stop("Not enough data.")
  model <- lm(Concentration ~ Response, data = data)
  summary_model <- summary(model)
  sd_rec_IP <- summary_model$coefficients[, "Std. Error"][1]
  
  model_r <- data %>% 
    filter(Days == 1) %>%
    summarise(
      model = list(lm(Concentration ~ Response, data = cur_data()))
    )
  model_r_summary <- summary(model_r$model[[1]])
  sd_rec_r <- model_r_summary$coefficients[1, "Std. Error"]
  
  found <- -model$coefficients[1]
  RSD_r <- 100 * sd_rec_r / found
  RSD_IP <- 100 * sd_rec_IP / found
   results <- data.frame(
    `Found` = found,
    `± sd Repeatability` = sd_rec_r,
    `RSD (%) Repeatability` = RSD_r,
    `± sd Intermediate Precision` = sd_rec_IP,
    `RSD (%) Intermediate Precision` = RSD_IP,
    check.names = FALSE # Impede que o R renomeie automaticamente as colunas
  )
  return(results)
}


server <- function(input, output, session) {
  

  
  #selectivity
  
  
  # ReactiveVal to store table data
  DF_estd <- reactiveVal(data.frame(
    concentration = numeric(5),
    response = numeric(5),
    stringsAsFactors = FALSE
  ))
  
  # Render initial table
  output$table_estd <- renderRHandsontable({
    rhandsontable(DF_estd())
  })
  
  # upgrade table data 
  observe({
    if (!is.null(input$table_estd)) {
      DF_estd(hot_to_r(input$table_estd))
    }
  })
  
  # ReactiveVal to store table data
  DF_std <- reactiveVal(data.frame(
    concentration = numeric(5),
    response = numeric(5),
    stringsAsFactors = FALSE
  ))
  
  #  Render initial table
  output$table_std <- renderRHandsontable({
    rhandsontable(DF_std())
  })
  
  # upgrade table data 
  observe({
    if (!is.null(input$table_std)) {
      DF_std(hot_to_r(input$table_std))
    }
  })

  
  # Calculate outliers to External standard
  outliers_estd <- reactive({
    data_selected <- DF_estd()
    colnames(data_selected) <- c("x", "estd")
    model <- lm(estd~x, data = data_selected)
    result_outliers_estd <-  outlierTest(model)
    return(result_outliers_estd)
  })
 
  # Render result
  output$outliers_estd <- renderPrint({
    outliers_estd()
  })
  
  # test status
  output$p_ouliers_estd <- renderUI({
    req(outliers_estd)
    tryCatch({
      result <- outliers_estd()
    p <- result$bonf.p
    obs <- names(result$rstudent)

    if (p > 0.05) {
      box(
        title = "Passed",
        status = "success", # verde
        solidHeader = TRUE,
        paste("The largest residual observation (", obs,") is not an outlier" )
      )
    } else {
      box(
        title = "Failed",
        status = "danger", # vermelho
        solidHeader = TRUE,
        paste("The largest residual observation (", obs,") is an outlier"))
    }
    }, error = function(e) {
      # Retorna uma mensagem de erro amigável
      message <- data.frame(Message = "No data to analyze.")
      return(message)
    })

  })
  # Calculate outliers to standard additions
  outliers_std <- reactive({
    data_selected <- DF_std()
    colnames(data_selected) <- c("x", "std")
    model <- lm(std~x, data = data_selected)
    result_outliers_std <-  outlierTest(model)
    return(result_outliers_std)
  })
  
  output$outliers_std <- renderPrint({
    outliers_std()
  })
  
  output$p_ouliers_std <- renderUI({
    req(outliers_std)
    tryCatch({
      result <- outliers_std()
    p <- result$bonf.p
    obs <- names(result$rstudent)
    
    if (p > 0.05) {
      box(
        title = "Passed", 
        status = "success", # verde
        solidHeader = TRUE,
        paste("The largest residual observation (", obs,") is not an outlier" )
      )
    } else {
      box(
        title = "Failed", 
        status = "danger", # vermelho
        solidHeader = TRUE,
        paste("The largest residual observation (", obs,") is an outlier"))
    }
    }, error = function(e) {
      # Retorna uma mensagem de erro amigável
      message <- data.frame(Message = "No data to analyze.")
      return(message)
    })

    })

  
    # Create youden calibration data table input

    DF_youden <- reactiveVal(data.frame(
      volume = numeric(5),
      response = numeric(5),
      stringsAsFactors = FALSE
    ))
    
    output$table_youden <- renderRHandsontable({
      rhandsontable(DF_youden())
    })
    
  
    observe({
      if (!is.null(input$table_youden)) {
        DF_youden(hot_to_r(input$table_youden))
      }
    })
  
    # Calculate Youden coefficients
  youden_result_coef <- reactive({
    youden_data <- DF_youden()
    result <- calculate_coefficients_compare(youden_data, "volume", "response")
    return(result)
  })
  
  estd_result_coef <- reactive({
    estd_data <- DF_estd()
    colnames(estd_data) <- c("x", "estd")
    result <- calculate_coefficients_compare(estd_data, "x", "estd")
    return(result)
  })
  
  std_result_coef <- reactive({
    
    std_data <- DF_std()
    colnames(std_data) <- c("x", "std")
    result <- calculate_coefficients_compare(std_data, "x", "std")
    return(result)
  })
  
  output$coefficients_youden<- renderPrint({
    youden_result_coef()
  })
  
  output$coefficients_estd<- renderPrint({
    estd_result_coef()
  })
  
  output$coefficients_std<- renderPrint({
    std_result_coef()
  })
  
  output$p_slopes_compare <- renderUI({
    req(estd_result_coef)
    req(std_result_coef)
    tryCatch({
      estd <- estd_result_coef()
    std <- std_result_coef()
    
    upper1 <- estd[3]
    lower1 <- estd[4]
    upper2 <- std[3]
    lower2 <- std[4]
    
    
    if (is.na(upper1)){
      return(
        box(
          "No data to analyze."
        )
      )
      
    } 
    # Verifica a sobreposição
    if (upper1 < lower2 || upper2 < lower1) {
      box(
        title = "Failed", 
        status = "danger", # vermelho
        solidHeader = TRUE,
        paste("There is proportional matrix effect"))
    } else {
      box(
        title = "Passed", 
        status = "success", # verde
        solidHeader = TRUE,
        paste("There is no proportional matrix effect" )
      )
    }
    })
    
    
  })
  
  
  output$p_intercepts_compare <- renderUI({
    req(estd_result_coef)
    req(youden_result_coef)
    tryCatch({
      estd <- estd_result_coef()
      youden <- youden_result_coef()
      
      upper1 <- estd[6]
      lower1 <- estd[7]
      upper2 <- youden[6]
      lower2 <- youden[7]
      

      if (upper1 == 0){
        return(
          box(
            "No data to analyze."
          )
        )
        
      }
      
    
    
    # Verifica a sobreposição
    if (upper1 < lower2 || upper2 < lower1) {
      box(
        title = "Failed", 
        status = "danger", # vermelho
        solidHeader = TRUE,
        paste("There is constant matrix effect"))
    } else {
      box(
        title = "Passed", 
        status = "success", # verde
        solidHeader = TRUE,
        paste("There is no constant matrix effect" )
      )
    }
    })
    
    
  })
  calibrations <- reactive({
    estd_data <- DF_estd()  # Obtenha os dados para 'estd'
    colnames(estd_data) <- c("x", "value")  # Renomeia as colunas para ter 'x' e 'value'
    
    std_data <- DF_std()  # Obtenha os dados para 'std'
    colnames(std_data) <- c("x", "value")  # Renomeia as colunas para ter 'x' e 'value'
    
    # Adiciona uma coluna para identificar o grupo (External Standard ou Standard Additions)
    estd_data$group <- "External Standard"
    std_data$group <- "Standard Additions"
    
    # Junta os dois dataframes com as mesmas colunas
    rbind(estd_data, std_data)  # Retorna o dataframe combinado
  })
  
  coefficients_overlap_plot <- reactive({
  
    calibrations <- calibrations()
    # Plot comparison of matrix-effect coefficients
    slopes_graph <- ggplot(calibrations) +
      geom_point(aes(x = x, y = value, color = group)) + # Plot points
      geom_smooth(aes(x = x, y = value, color = group), se = TRUE, method = "lm") + # Linhas de ajuste
      labs(x = input$x_label, y = input$y_label) +
      scale_color_discrete(labels = c("External Standard", "Standard Additions"))
    
    
    # Apply the user-selected theme
    selected_theme <- switch(input$theme,
                             "theme_grey" = theme_grey(15),
                             "theme_bw" = theme_bw(15),
                             "theme_classic" = theme_classic(15),
                             "theme_minimal" = theme_minimal(15),
                             "theme_void" = theme_void(15)
    )
    
    slopes_graph  <- slopes_graph + selected_theme + theme(legend.title = element_blank(),
                                                           legend.position = "top")
    
    # return the graph
    return(slopes_graph)
    
  })
  
  output$coefficients_overlap_plot <- renderPlot({
    coefficients_overlap_plot()
  })
  
  #Regression analysis
  #select calibration data
  # ReactiveVal para armazenar os dados da tabela
  DF_linearity <- reactiveVal(data.frame(
    concentration = numeric(5),
    response = numeric(5),

    stringsAsFactors = FALSE
  ))
  
  # Renderiza a tabela inicial
  output$table_linearity <- renderRHandsontable({
    rhandsontable(DF_linearity())
  })
  
  # Atualiza os dados da tabela automaticamente
  observe({
    if (!is.null(input$table_linearity)) {
      DF_linearity(hot_to_r(input$table_linearity))
    }
  })
  
  regression_plot <- reactive({
    data <- DF_linearity()
    
    # Fit the linear model
    model <- lm(response~concentration, data = data)
    
    # Create prediction data with prediction intervals
    prediction_data <- cbind(data, predict(model, interval = "prediction"))
    
    # Create the plot with prediction intervals
    plot <- ggplot(prediction_data, aes(x = concentration, y = response)) +
      geom_point() +  # Plot the points with selected color
      geom_smooth(method = "lm", se = TRUE) + 
      geom_line(aes(y = upr, color = "Prediction Interval", linetype = "PI"), size = 1) +  # Upper prediction interval
      geom_line(aes(y = lwr, color = "Prediction Interval", linetype = "PI"), size = 1) +  # Lower prediction interval
      scale_color_manual(values = c("Prediction Interval" = "black")) +  # Custom colors
      scale_linetype_manual(values = c("PI" = "dashed")) +  # Dashed line for PI
      labs(x = input$x_label, y = input$y_label)
    
    # Apply the user-selected theme
    selected_theme <- switch(input$theme,
                             "theme_grey" = theme_grey(15),
                             "theme_bw" = theme_bw(15),
                             "theme_classic" = theme_classic(15),
                             "theme_minimal" = theme_minimal(15),
                             "theme_void" = theme_void(15)
    )
    
    # Add the theme to the plot
    plot <- plot + selected_theme + theme(legend.position='none')
    
    print(plot)
  })
  
  # Regression plot with prediction intervals
  output$regression_plot <- renderPlot({
    regression_plot()
  })
  
  residual_plot <- reactive({
    data <- DF_linearity()
    # Fit the linear model
    model <- lm(response~concentration, data = data)
    df_reg <- fortify(model)
    r1 <-ggplot(df_reg, aes(y = .resid, x = .fitted)) +
      geom_point() +
      xlab(expression("Fitted Values"))+
      ylab(expression(Residuals))+
      geom_hline(yintercept = 0, linetype= "dashed")
    
    # Apply the user-selected theme
    selected_theme <- switch(input$theme,
                             "theme_grey" = theme_grey(15),
                             "theme_bw" = theme_bw(15),
                             "theme_classic" = theme_classic(15),
                             "theme_minimal" = theme_minimal(15),
                             "theme_void" = theme_void(15)
    )
    r1 <- r1 + selected_theme + theme(legend.position='none')
    return(r1)
  })
  
  output$residual_plot <- renderPlot({
    residual_plot()
  })
  
  
  result <- reactive({
    data <- DF_linearity()
    
    statistics(data, "concentration",  "response")
  })
  
  regression <- reactive({
    tryCatch({
      reg <- result()$Regression[[1]]
    reg <- reg %>%
      mutate(across(where(~ !all(is.character(.))), ~ sprintf("%.2e", .))) 
    return(reg)
    }, error = function(e) {
      # Retorna uma mensagem de erro amigável
      message <- data.frame(Message = "No data to analyze.")
      return(message)
    })
    
  })
  
  output$regression <- renderTable({
    regression()
  })
  
  coef.determination <- reactive({
    tryCatch({
       r2 <- result()$Coef_deter[[1]]
    }, error = function(e) {
      # Retorna uma mensagem de erro amigável
      message <- data.frame(Message = "No data to analyze.")
      return(message)
    })
   
  })
  
  output$coef.determination <- renderTable({
    coef.determination()
  })
  
  levene <- reactive({
    tryCatch({
      result()$Levene[[1]]
    }, error = function(e) {
      # Retorna uma mensagem de erro amigável
      message <- data.frame(Message = "No data to analyze.")
      return(message)
    })
    
  })
  
  output$levene <- renderTable(levene())
  
  
  output$p_levene <- renderUI({
    tryCatch({
      p <- result()$Levene[[1]][2]
    
    if (p > 0.05) {
      box(
        title = "Passed", 
        status = "success", # verde
        solidHeader = TRUE,
        paste("The data are homogeneous" )
      )
    } else {
      box(
        title = "Failed", 
        status = "danger", # vermelho
        solidHeader = TRUE,
        paste("The data are not homogeneous"))
    }
    }, error = function(e) {
      # Retorna uma mensagem de erro amigável
      message <- data.frame(Message = "No data to analyze.")
      return(message)
    })
    
  })
  
  durbin_watson <- reactive({
    tryCatch({
      result()$Durbin_Watson[[1]]
    }, error = function(e) {
      # Retorna uma mensagem de erro amigável
      message <- data.frame(Message = "No data to analyze.")
      return(message)
    })
    
  })
  
  output$durbin_watson <- renderTable(durbin_watson() )
  
  output$p_dw <- renderUI({
    tryCatch({
      p <- result()$Durbin_Watson[[1]][2]
    
    if (p > 0.05) {
      box(
        title = "Passed", 
        status = "success", # verde
        solidHeader = TRUE,
        paste("The errors are independent" )
      )
    } else {
      box(
        title = "Failed", 
        status = "danger", # vermelho
        solidHeader = TRUE,
        paste("The errors are not independent"))
    }
    }, error = function(e) {
      # Retorna uma mensagem de erro amigável
      message <- data.frame(Message = "No data to analyze.")
      return(message)
    })
    
  })
  
  shapiro_wilk <- reactive({
    tryCatch({
      result()$Shapiro_Wilk[[1]]
    }, error = function(e) {
      # Retorna uma mensagem de erro amigável
      message <- data.frame(Message = "No data to analyze.")
      return(message)
    })
    
  })
  
  
  output$shapiro_wilk <- renderTable(shapiro_wilk())

  
  output$p_normality <- renderUI({
    tryCatch({
       p <- result()$Shapiro_Wilk[[1]][2]
    
  if (p > 0.05) {
    box(
      title = "Passed", 
      status = "success", # verde
      solidHeader = TRUE,
      paste("Data follow normal distribution" )
    )
  } else {
    box(
      title = "Failed", 
      status = "danger", # vermelho
      solidHeader = TRUE,
      paste("Data does not follow normal distribution"))
  }
    }, error = function(e) {
      # Retorna uma mensagem de erro amigável
      message <- data.frame(Message = "No data to analyze.")
      return(message)
    })
   
})
  
  anova <- reactive({
    tryCatch({
      aov <- result()$Anova[[1]]
    aov <- aov %>%
      mutate(across(where(~ !all(is.character(.))), ~ sprintf("%.2e", .))) 
    return(aov)
    }, error = function(e) {
      # Retorna uma mensagem de erro amigável
      message <- data.frame(Message = "No data to analyze.")
      return(message)
    })
    
  })
  
  output$anova <- renderTable({
    anova()
  })
  

  output$p_anova_reg <- renderUI({
    tryCatch({
      p <- result()$Anova[[1]][6][1,]

    if (p < 0.05) {
      box(
        title = "Passed",
        status = "success", # verde
        solidHeader = TRUE,
        paste("The regression is significant" )
      )
    } else {
      box(
        title = "Failed",
        status = "danger", # vermelho
        solidHeader = TRUE,
        paste("The regression is not significant"))
    }
    }, error = function(e) {
      # Retorna uma mensagem de erro amigável
      message <- data.frame(Message = "No data to analyze.")
      return(message)
    })
    
  })
  
  output$p_anova_lof <- renderUI({
    tryCatch({
      p <- result()$Anova[[1]][6][2,]
      
      if (p > 0.05) {
        box(
          title = "Passed", 
          status = "success", # verde
          solidHeader = TRUE,
          paste("There is no lack of fit" )
        )
      } else {
        box(
          title = "Failed", 
          status = "danger", # vermelho
          solidHeader = TRUE,
          paste("There is lack of fit"))
      }
    }, error = function(e) {
      # Retorna uma mensagem de erro amigável
      message <- data.frame(Message = "No data to analyze.")
      return(message)
    })
   
  })
  
  limits_result <- reactive({
    data <- DF_linearity()
    tryCatch({
      limits(data, "concentration", "response")
    }, error = function(e) {
      # Retorna uma mensagem de erro amigável
      message <- data.frame(Message = "No data to analyze.")
      return(message)
    })
    
  })
  
  output$limits <- renderTable({
    limits_result()
    })
  
  # Accuracy 
  
  selected_option <- reactiveVal("none")
  
  observeEvent(input$acc_estd,{
    selected_option("option1")
  })
  
  observeEvent(input$acc_std,{
    selected_option("option2")
  })
  
  output$acc_results <- renderUI({
    
    if (selected_option() == "option1") {
      fluidRow(
        box(
          h4("Enter Responses to Non-fortified Sample"),
          rHandsontableOutput("table_Q0_estd"),
          
          h4("Enter Responses to Q1 fortified Sample"),
          numericInput(
            inputId = "Q1_value_estd", 
            label = "Enter Q1 correspondent concentration:", 
            value = 0,       # Valor inicial
            min = 0,         # Valor mínimo
            max = Inf,       # Valor máximo
            step = 1         # Incremento
          ),
          rHandsontableOutput("table_Q1_estd"),
          
          h4("Enter Responses to Q2 fortified Sample"),
          numericInput(
            inputId = "Q2_value_estd", 
            label = "Enter Q2 correspondent concentration:", 
            value = 0,       # Valor inicial
            min = 0,         # Valor mínimo
            max = Inf,       # Valor máximo
            step = 1         # Incremento
          ),
          rHandsontableOutput("table_Q2_estd"),
          
          h4("Enter Responses to Q3 fortified Sample"),
          numericInput(
            inputId = "Q3_value_estd", 
            label = "Enter Q3 correspondent concentration:", 
            value = 0,       # Valor inicial
            min = 0,         # Valor mínimo
            max = Inf,       # Valor máximo
            step = 1         # Incremento
          ),
          rHandsontableOutput("table_Q3_estd")
          
        ),
        
        box(
          h4("Non-fortified results"),
          verbatimTextOutput("Q0_results_estd"),
          h4("Q1 results"),
          verbatimTextOutput("Q1_results_estd"),
          h4("Q2 results"),
          verbatimTextOutput("Q2_results_estd"),
          h4("Q3 results"),
          verbatimTextOutput("Q3_results_estd")
        )
      )
      
      
      
      
      
    } else if (selected_option() == "option2") {
      fluidRow(
          box(
            h4("Enter Std. Additions to Non-fortified Sample"),
            rHandsontableOutput("table_Q0"),
            
            h4("Enter Std. Additions to Q1 fortified Sample"),
            numericInput(
              inputId = "Q1_value", 
              label = "Enter Q1 correspondent concentration:", 
              value = 0,       # Valor inicial
              min = 0,         # Valor mínimo
              max = Inf,       # Valor máximo
              step = 1         # Incremento
            ),
            rHandsontableOutput("table_Q1"),
          
            h4("Enter Std. Additions to Q2 fortified Sample"),
              numericInput(
                inputId = "Q2_value", 
                label = "Enter Q2 correspondent concentration:", 
                value = 0,       # Valor inicial
                min = 0,         # Valor mínimo
                max = Inf,       # Valor máximo
                step = 1         # Incremento
              ),
              rHandsontableOutput("table_Q2"),
            
             h4("Enter Std. Additions to Q3 fortified Sample"),
              numericInput(
                inputId = "Q3_value", 
                label = "Enter Q3 correspondent concentration:", 
                value = 0,       # Valor inicial
                min = 0,         # Valor mínimo
                max = Inf,       # Valor máximo
                step = 1         # Incremento
              ),
              rHandsontableOutput("table_Q3")
          ),
          
          box(
            h4("Non-fortified results"),
            verbatimTextOutput("Q0_results"),
            h4("Q1 results"),
            verbatimTextOutput("Q1_results"),
            h4("Q2 results"),
            verbatimTextOutput("Q2_results"), 
            h4("Q3 results"),
            verbatimTextOutput("Q3_results")
          )
          
       
       
      )
    } 
  })
  
  # Q0 data
  
  DF_Q0_estd <- reactiveVal(data.frame(
    Days = numeric(5),
    Response = numeric(5),
    stringsAsFactors = FALSE
  ))
  
  # Renderiza a tabela inicial
  output$table_Q0_estd <- renderRHandsontable({
    rhandsontable(DF_Q0_estd())
  })
  
  # Atualiza os dados da tabela automaticamente
  observe({
    if (!is.null(input$table_Q0_estd)) {
      DF_Q0_estd(hot_to_r(input$table_Q0_estd))
    }
  })
  
  
  # Q1 data
  
  DF_Q1_estd <- reactiveVal(data.frame(
    Days = numeric(5),
    Response = numeric(5),
    stringsAsFactors = FALSE
  ))
  
  # Renderiza a tabela inicial
  output$table_Q1_estd <- renderRHandsontable({
    rhandsontable(DF_Q1_estd())
  })
  
  # Atualiza os dados da tabela automaticamente
  observe({
    if (!is.null(input$table_Q1_estd)) {
      DF_Q1_estd(hot_to_r(input$table_Q1_estd))
    }
  })
  
  #Q2 data
  
  DF_Q2_estd <- reactiveVal(data.frame(
    Days = numeric(5),
    Response = numeric(5),
    stringsAsFactors = FALSE
  ))
  
  # Renderiza a tabela inicial
  output$table_Q2_estd <- renderRHandsontable({
    rhandsontable(DF_Q2_estd())
  })
  
  # Atualiza os dados da tabela automaticamente
  observe({
    if (!is.null(input$table_Q2_estd)) {
      DF_Q2_estd(hot_to_r(input$table_Q2_estd))
    }
  })
  
  #Q3 data
  
  DF_Q3_estd <- reactiveVal(data.frame(
    Days = numeric(5),
    Response = numeric(5),
    stringsAsFactors = FALSE
  ))
  
  # Renderiza a tabela inicial
  output$table_Q3_estd <- renderRHandsontable({
    rhandsontable(DF_Q3_estd())
  })
  
  # Atualiza os dados da tabela automaticamente
  observe({
    if (!is.null(input$table_Q3_estd)) {
      DF_Q3_estd(hot_to_r(input$table_Q3_estd))
    }
  })




  # Precision

  Q0_results_estd <- reactive({
    req(DF_Q0_estd())
    req(DF_linearity())
    data <- DF_Q0_estd()
    
    data_regression <- DF_linearity()
    
    tryCatch({
      # Fit the linear model
      model <- lm(response ~ concentration, data = data_regression)
      
      b <- model[["coefficients"]][["concentration"]]
      a <- model[["coefficients"]][["(Intercept)"]]
      
      # Use mutate to handle the condition for each row
      data <- data %>%
        mutate(
          predicted.concentration = if_else(Response == 0, 0, (Response - a) / b)
        )
      
      results <- data %>%
        summarise(
          predicted.concentration.mean = mean(predicted.concentration),
          sd_predicted = sd(predicted.concentration)
        ) %>%
        select(predicted.concentration.mean, sd_predicted)  # Remove a coluna do modelo linear para não aparecer no output
      
      # Ajustando nomes das colunas
      colnames(results) <- c("Found", "± sd")
      
      return(results)
    }, error = function(e) {
      data.frame(Message = "No calibration data in linearity section.")
    })
  })
  
  
  output$Q0_results_estd <- renderPrint({
    Q0_results_estd()
  })


  Q1_results_estd <- reactive({
    req(Q0_results_estd())
    req(DF_linearity())
    
    data <- DF_Q1_estd()
    sample <- Q0_results_estd()[["Found"]]
    concentration <-  as.numeric(input$Q1_value_estd)
    data_regression <- DF_linearity()
    
    tryCatch({
      
      # Fit the linear model
      model <- lm(response~concentration, data = data_regression)
      
      b <- model[["coefficients"]][["concentration"]]
      a <- model[["coefficients"]][["(Intercept)"]]
      
      data <- data %>%
        mutate(
          predicted.concentration = (Response-a)/b
        )
      
      results <- data %>%
        summarise(
          reg_precision = list(lm(predicted.concentration  ~ Days, data = cur_data())),
          MSw = stats::anova(reg_precision[[1]])$`Mean Sq`[2],
          MSb = stats::anova(reg_precision[[1]])$`Mean Sq`[1],
          n = 2 + stats::anova(reg_precision[[1]])$Df[2],
          s_r = sqrt(MSw),
          s_b = sqrt(abs(MSb - MSw / n)),
          s_I = sqrt(s_r^2 + s_b^2),
          predicted.concentration.mean = mean(predicted.concentration),
          predicted_minus_sample = predicted.concentration.mean - sample,
          mean_r = mean(predicted.concentration[Days == 1])- sample,
          repeatability =  100 * s_r / mean_r,
          intermediate.precision = 100*s_I / predicted.concentration.mean,
          rec = 100 * predicted_minus_sample / concentration
        ) %>%
        select(predicted.concentration.mean, s_r, repeatability, s_I, intermediate.precision, rec)  # Remove a coluna do modelo linear para não aparecer no output
      
      # Ajustando nomes das colunas
      colnames(results) <- c("Found", "± sd Repeatability", "RSD Repeatability (%)", "± sd Intermediate Precision"," RSD Intermediate Precision (%)", "Recovery (%)")
     
      return(results)
    }, error = function(e) {
      data.frame(Message = "No calibration data in linearity section.")
    })
  })
  
  output$Q1_results_estd <- renderPrint({
    Q1_results_estd()
  })

  
  Q2_results_estd <- reactive({
    req(Q0_results_estd())
    req(DF_linearity())
    
    data <- DF_Q2_estd()
    sample <- Q0_results_estd()[["Found"]]
    concentration <-  as.numeric(input$Q2_value_estd)
    data_regression <- DF_linearity()
    
    tryCatch({
      
      # Fit the linear model
      model <- lm(response~concentration, data = data_regression)
      
      b <- model[["coefficients"]][["concentration"]]
      a <- model[["coefficients"]][["(Intercept)"]]
      
      data <- data %>%
        mutate(
          predicted.concentration = (Response-a)/b
        )
      
      results <- data %>%
        summarise(
          reg_precision = list(lm(predicted.concentration  ~ Days, data = cur_data())),
          MSw = stats::anova(reg_precision[[1]])$`Mean Sq`[2],
          MSb = stats::anova(reg_precision[[1]])$`Mean Sq`[1],
          n = 2 + stats::anova(reg_precision[[1]])$Df[2],
          s_r = sqrt(MSw),
          s_b = sqrt(abs(MSb - MSw / n)),
          s_I = sqrt(s_r^2 + s_b^2),
          predicted.concentration.mean = mean(predicted.concentration),
          predicted_minus_sample = predicted.concentration.mean - sample,
          mean_r = mean(predicted.concentration[Days == 1])- sample,
          repeatability =  100 * s_r / mean_r,
          intermediate.precision = 100*s_I / predicted.concentration.mean,
          rec = 100 * predicted_minus_sample / concentration
        ) %>%
        select(predicted.concentration.mean, s_r, repeatability, s_I, intermediate.precision, rec)  # Remove a coluna do modelo linear para não aparecer no output
      
      # Ajustando nomes das colunas
      colnames(results) <- c("Found", "± sd Repeatability", "RSD Repeatability (%)", "± sd Intermediate Precision"," RSD Intermediate Precision (%)", "Recovery (%)")
      
      return(results)
    }, error = function(e) {
      data.frame(Message = "No calibration data in linearity section.")
    })
  })
  
  output$Q2_results_estd <- renderPrint({
    Q2_results_estd()
  })
  
  
  Q3_results_estd <- reactive({
    req(Q0_results_estd())
    req(DF_linearity())
    
    data <- DF_Q3_estd()
    sample <- Q0_results_estd()[["Found"]]
    concentration <-  as.numeric(input$Q3_value_estd)
    data_regression <- DF_linearity()
    
    tryCatch({
      
      # Fit the linear model
      model <- lm(response~concentration, data = data_regression)
      
      b <- model[["coefficients"]][["concentration"]]
      a <- model[["coefficients"]][["(Intercept)"]]
      
      data <- data %>%
        mutate(
          predicted.concentration = (Response-a)/b
        )
      
      results <- data %>%
        summarise(
          reg_precision = list(lm(predicted.concentration  ~ Days, data = cur_data())),
          MSw = stats::anova(reg_precision[[1]])$`Mean Sq`[2],
          MSb = stats::anova(reg_precision[[1]])$`Mean Sq`[1],
          n = 2 + stats::anova(reg_precision[[1]])$Df[2],
          s_r = sqrt(MSw),
          s_b = sqrt(abs(MSb - MSw / n)),
          s_I = sqrt(s_r^2 + s_b^2),
          predicted.concentration.mean = mean(predicted.concentration),
          predicted_minus_sample = predicted.concentration.mean - sample,
          mean_r = mean(predicted.concentration[Days == 1])- sample,
          repeatability =  100 * s_r / mean_r,
          intermediate.precision = 100*s_I / predicted.concentration.mean,
          rec = 100 * predicted_minus_sample / concentration
        ) %>%
        select(predicted.concentration.mean, s_r, repeatability, s_I, intermediate.precision, rec)  # Remove a coluna do modelo linear para não aparecer no output
      
      # Ajustando nomes das colunas
      colnames(results) <- c("Found", "± sd Repeatability", "RSD Repeatability (%)", "± sd Intermediate Precision"," RSD Intermediate Precision (%)", "Recovery (%)")
      
      return(results)
    }, error = function(e) {
      data.frame(Message = "No calibration data in linearity section.")
    })
  })
  
  output$Q3_results_estd <- renderPrint({
    Q3_results_estd()
  })
  
  
  # Q0 data
  
  DF_Q0 <- reactiveVal(data.frame(
    Days = numeric(5),
    Concentration = numeric(5),
    Response = numeric(5),
    stringsAsFactors = FALSE
  ))
  
  # Renderiza a tabela inicial
  output$table_Q0 <- renderRHandsontable({
    rhandsontable(DF_Q0())
  })
  
  # Atualiza os dados da tabela automaticamente
  observe({
    if (!is.null(input$table_Q0)) {
      DF_Q0(hot_to_r(input$table_Q0))
    }
  })
  
  found_Q0 <- reactive({
    data <- DF_Q0()
    tryCatch({
      model <- lm(Concentration~Response, data)
      summary_model <- summary(model)
      sd_rec <- summary_model[["coefficients"]][, "Std. Error"][1]
      rec <- -model[["coefficients"]][1]
      results<-data.frame(rec, sd_rec)
      colnames(results) <- c("Found", "±sd")
      return(results)
    }, error = function(e) {
      # Retorna uma mensagem de erro amigável
      message <- data.frame(Message = "No data to analyze.")
      return(message)
    })
  })
  
  output$Q0_results <- renderPrint({
    found_Q0()
  })
  
  # Q1 data
  
  DF_Q1 <- reactiveVal(data.frame(
    Days = numeric(5),
    Concentration = numeric(5),
    Response = numeric(5),
    stringsAsFactors = FALSE
  ))
  
  
  
  # Renderiza a tabela inicial
  output$table_Q1 <- renderRHandsontable({
    rhandsontable(DF_Q1())
  })
  
  # Atualiza os dados da tabela automaticamente
  observe({
    if (!is.null(input$table_Q1)) {
      DF_Q1(hot_to_r(input$table_Q1))
    }
  })

  found_Q1 <- reactive({
    req(found_Q0)
    data <- DF_Q1()
    concentration <- as.numeric(input$Q1_value)
    sample <- found_Q0()$Found
    
    result <- tryCatch({
      calculate_found(data) %>% 
        mutate(`Recovery (%)` = 100 * (`Found` - sample) / concentration)
    }, error = function(e) {
      data.frame(Message = "No data to analyze.")
    })
    
   return(result)
  })
  
  output$Q1_results <- renderPrint({
      found_Q1()
  })
  
  # Q2 data
  DF_Q2 <- reactiveVal(data.frame(
    Days = numeric(5),
    Concentration = numeric(5),
    Response = numeric(5),
    stringsAsFactors = FALSE
  ))
  
  # Renderiza a tabela inicial
  output$table_Q2 <- renderRHandsontable({
    rhandsontable(DF_Q2())
  })
  
  # Atualiza os dados da tabela automaticamente
  observe({
    if (!is.null(input$table_Q2)) {
      DF_Q2(hot_to_r(input$table_Q2))
    }
  })
  
  
  
  found_Q2 <- reactive({
    req(found_Q0)
    data <- DF_Q2()
    concentration <- as.numeric(input$Q2_value)
    sample <- found_Q0()$Found
    
    result <- tryCatch({
      calculate_found(data) %>% 
        mutate(`Recovery (%)` = 100 * (`Found` - sample) / concentration)
    }, error = function(e) {
      data.frame(Message = "No data to analyze.")
    })
    
    return(result)
  })
  
  output$Q2_results <-renderPrint({
    found_Q2()
  })
  
  
  # Q3 data
  
  DF_Q3 <- reactiveVal(data.frame(
    Days = numeric(5),
    Concentration = numeric(5),
    Response = numeric(5),
    stringsAsFactors = FALSE
  ))
  
  # Renderiza a tabela inicial
  output$table_Q3 <- renderRHandsontable({
    rhandsontable(DF_Q3())
  })
  
  # Atualiza os dados da tabela automaticamente
  observe({
    if (!is.null(input$table_Q3)) {
      DF_Q3(hot_to_r(input$table_Q3))
    }
  })
  found_Q3 <- reactive({
    req(found_Q0)
    data <- DF_Q3()
    concentration <- as.numeric(input$Q3_value)
    sample <- found_Q0()$Found
    
    result <- tryCatch({
    calculate_found(data) %>% 
        mutate(`Recovery (%)` = 100 * (`Found` - sample) / concentration)
  }, error = function(e) {
    data.frame(Message = "No data to analyze.")
  })
    return(result)
  })
  
  output$Q3_results <- renderPrint({
    found_Q3()
  })
  
  
  

  output$report <- downloadHandler(
    
    filename = function() {
      paste('validation_report', Sys.Date(), switch(
        input$format,
        PDF = 'pdf',
        HTML = 'html',
        Word = 'docx'
      ), sep = '.')
    },
    
    content = function(file) {
      # Define o caminho temporário para o arquivo .Rmd
      tempReport <- file.path(tempdir(), "report.Rmd")
      
      # Copia o arquivo Rmd para o diretório temporário
      file.copy("report.Rmd", tempReport, overwrite = TRUE)
      
      # Set up parameters to pass to Rmd document
      params <- list(
                     outliers_estd = outliers_estd(),
                     outliers_std = outliers_std(),
                     coefficients_estd = estd_result_coef(),
                     coefficients_std = std_result_coef(),
                     coefficients_overlap_plot = coefficients_overlap_plot(),
                     regression_plot = regression_plot(),
                     residual_plot = residual_plot(),
                     regression = regression(),
                     coef.determination = coef.determination(),
                     levene = levene(),
                     durbin_watson = durbin_watson (),
                     shapiro_wilk = shapiro_wilk(),
                     anova = anova(),
                     limits_result = limits_result(),
                     ext.std.precision = rsd_precision(),
                     ext.std.trueness = trueness(),
                     std.add.Q0 = found_Q0(),
                     std.add.Q1 = found_Q1(),
                     std.add.Q2 = found_Q2(),
                     std.add.Q3 = found_Q3()
      )
      
      
      # Knit the document, passing in the `params` list, and eval it in a
      # child of the global environment (this isolates the code in the document
      # from the code in this app).
      # Renderiza o documento Rmd
      rmarkdown::render(tempReport, output_file = file,
                        params = params,
                        output_format <- switch(input$format,
                                                PDF = rmarkdown::pdf_document(),
                                                HTML = rmarkdown::html_document(),
                                                Word = rmarkdown::word_document()),
                        envir = new.env(parent = globalenv())
      )
    }
  )
}