#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    https://shiny.posit.co/
#
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

set.seed(123)  # Para reprodutibilidade

concentration <- rep(c(1, 2, 3, 4, 5, 6, 7, 8), each = 2)  # 8 concentrações em duplicata
response1 <- concentration * 0.95 + rnorm(16, mean = 0, sd = 0.05)  # Simula resposta linear com erro aleatório
response2 <- concentration * 1.05 + rnorm(16, mean = 0, sd = 0.05)

data1 <- data.frame(Concentrations = concentration, `External standard`= response1, `Standard additions`= response2)
example_data <- data1

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


calculate_coefficients_compare <- function(data, calibration_type, x, y) {
  data %>%
    group_by(!!sym(calibration_type)) %>%
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
    select(-model, -model_tidy, -t, -df, -slope, -slope.std.error, -intercept, -intercept.std.error) # Remove undesired vectors
}

outliers <- function(data, x, y) {
  data %>% 
    reframe(
      model = list(lm(!!sym(y) ~ !!sym(x), data = cur_data())),
      outliers_r = list(outlierTest(model[[1]])),
      largest.rstudent.obs = names(outliers_r[[1]]$rstudent[1]),
      rstudent = outliers_r[[1]]$rstudent[1],
      p.value = outliers_r[[1]]$p,
      bonf.p = outliers_r[[1]]$bonf.p
    ) %>% 
    select(largest.rstudent.obs,rstudent, p.value, bonf.p)
}

ui <- navbarPage(
  theme = shinytheme("sandstone"),
  title = "Analytical Methods Validation",
  
  tabPanel("Data Upload",
           
           h4("The data below provide an example of how you should structure your table. They are also used as an example to demonstrate the other validation parameters. Once you upload your validation table, these data will be replaced by your own. In the other sections, you should modify the required parameters according to your data."),
           fileInput("dataset_file", "Choose your CSV file with validation experiments data", accept = ".csv"),
           dataTableOutput("your_data")
           
           
  ),
  
  
  tabPanel("Selectivity",
           sidebarLayout(
             sidebarPanel(
               
               uiOutput("xcol"),
               uiOutput("estd_col"),
               uiOutput("std_col"),
               
               selectInput("theme", "Plot Theme:", choices = list(
                 "theme_grey" = "theme_grey",
                 "theme_bw" = "theme_bw",
                 "theme_classic" = "theme_classic",
                 "theme_minimal" = "theme_minimal",
                 "theme_void" = "theme_void"
               ), selected = "theme_classic"),
               
               textInput("x_label", "X-axis label:", value = "Concentration"),
               textInput("y_label", "Y-axis label:", value = "Response")
             ),
             
             mainPanel(
               h3("Results to Matrix Effect"),
               h4("Bonferroni outliers test to external standard"),
               tableOutput("outliers_estd"),
               h4("Bonferroni outliers test to standard additions"),
               tableOutput("outliers_std"),
               h4("Plot overlap"),
               h5("The linear regression performed to different methods plotted together"),
               plotOutput("coefficients_overlap_plot", width = "70%"),
               h4("Coefficiet overlap"),
               dataTableOutput("coefficients_overlap")
             )
           )
  ),
  
  tabPanel(
    "Regression Analysis and Linearity",
    sidebarLayout(
      sidebarPanel(
        uiOutput("x_reg_col"),
        uiOutput("y_reg_col"),
        
      ),
      mainPanel(
        h3("Results"),
        fluidRow(
          column(6, 
                 h4("Linear Regression"),
                 plotOutput("regression_plot", height = "4in"),
                 h6("The regression line in blue, confidence band in gray and prediction band in dashed black line")
          ),
          column(6, 
                 h4("Residual graph"),
                 plotOutput("residual_plot", height = "4in")
          )
        ),
        h4("Coefficient of determination"),
        tableOutput("coef.determination"),
        h4("Regressions terms"),
        tableOutput("regression"),
        h4("Levene homogeneity test"),
        tableOutput("levene"),
        h4("Durbin-Watson test to auto-cooreleted errors"),
        tableOutput("durbin_watson"),
        h4("Shapiro-Wilk normality test"),
        tableOutput("shapiro_wilk"),
        h4("ANOVA of regression and lack-of-fit"),
        tableOutput("anova")
      )
    )
  ),
  
  tabPanel(
    "LOD and LOQ",
    mainPanel(
      h3("LOD and LOQ with calibration curve method"),
      h4("The formula used is LOD = 3.3s/b, and LOQ = 10s/b, where s is the residual standard deviation of a regression line and b the inclination."),
      actionButton("limits","Calculate LOD and LOQ"),
      tableOutput("limits"),
      h5("The values are in terms of concentration")
      
      
    )
  ),
  tabPanel(
    "Download Report",
    mainPanel(
      radioButtons('format', 'Document format', c('PDF', 'HTML', 'Word'),
                   inline = TRUE),
      downloadButton("report", "Generate report")
    )
  )
 
)



server <- function(input, output, session) {
  
  dataset <- reactive({
    req(input$dataset_file)
    data <- read.csv(input$dataset_file$datapath, header = TRUE, sep = ",", row.names = NULL)
    rownames(data) <- make.unique(as.character(seq_len(nrow(data))))
    return(data)
  })
  
  # Reactive para carregar os dados do usuário ou mostrar o exemplo
  dataset <- reactive({
    # Se o arquivo não for carregado, retorna o dataset de exemplo
    if (is.null(input$dataset_file)) {
      return(example_data)
    }
    
    # Se o arquivo for carregado, lê o arquivo CSV
    read.csv(input$dataset_file$datapath)
  })
  
  # Renderiza a tabela de dados
  output$your_data <- renderDataTable({
    datatable(dataset(), caption = htmltools::tags$caption(
      style = 'caption-side: top; text-align: center; font-weight: bold; font-size: 18px;color: black;',
      "Table structure"
    ))
  })
  
  
#selectivity
  
  # Popula as opções para seleção de X e Y com base nos dados do dataset carregado
    output$xcol <- renderUI({
      req(dataset())
      selectInput("xcol", "Select x variable", names(dataset()))
    })
    
    output$estd_col <- renderUI({
      req(dataset())
      selectInput("estd_col", "Select external standard column", names(dataset()), selected = "External.standard")
    })
    
    output$std_col <- renderUI({
      req(dataset())
      selectInput("std_col", "Select standard additions column", names(dataset()), selected = "Standard.additions")
    })
    
  outliers_estd <- reactive({
    data_selected <- dataset() %>%
    select(all_of(c(input$xcol,input$estd_col)))
    colnames(data_selected) <- c("x", "estd")
    
    result_outliers_estd <-  outliers(data_selected, "x", "estd")
    colnames(result_outliers_estd) <- c("Largest rstudent observation", "rstudent", "p-value", "Bonferroni p")
    return(result_outliers_estd)
  })
  
  output$outliers_estd <- renderTable({
    outliers_estd() 
  })
  
  outliers_std <- reactive({
      data_selected <- dataset() %>%
      select(all_of(c(input$xcol,input$std_col)))
      colnames(data_selected) <- c("x", "std")
      
      result_outliers_std<-  outliers(data_selected, "x", "std")
      colnames(result_outliers_std) <- c("Largest rstudent observation", "rstudent", "p-value", "Bonferroni p")
      return(result_outliers_std)
  })
  
  output$outliers_std <- renderTable({
    outliers_std()
  })
  
  
  coefficients_overlap <- reactive({
    
    # Cria o dataframe 'data_selected' com as colunas escolhidas pelo usuário
    data_selected <- dataset() %>%
      select(all_of(c(input$xcol, input$std_col, input$estd_col)))
    
    # Reorganiza os dados
    data_selected_melt <- melt(data_selected, id.vars = input$xcol)
    colnames(data_selected_melt) <- c("x", "calibration_type", "y")
    
    
    # Calcula os coeficientes
    result_coefficients <- calculate_coefficients_compare(data_selected_melt, "calibration_type", "x", "y")
    colnames(result_coefficients) <- c("Calibration type", "Slope CI", "Upper slope value", "Lower slope value", "Intercept CI", "Upper intercept value", "Lower intercept value")
    
    result_coefficients <- result_coefficients %>%
      mutate(`Calibration type` = as.character(`Calibration type`))  # Converte para texto
    
    
    
    # Formata os números em notação científica com 3 casas decimais
    result_coefficients_formatted <- result_coefficients %>%
      mutate(across(where(~ !all(is.character(.))), ~ sprintf("%.2e", .)))  # Exclui colunas de texto
    
    return(result_coefficients_formatted)
   
  })
  
  output$coefficients_overlap <- renderDT({
    # Exibe a tabela com DT
    datatable(coefficients_overlap(), options = list(
      pageLength = 10,
      autoWidth = TRUE
    ))
  })
  
  
  coefficients_overlap_plot <- reactive({
    data_selected <- dataset() %>%
      select(all_of(c(input$xcol, input$std_col, input$estd_col)))
    
    data_selected_melt <- melt(data_selected, id.vars = input$xcol)
    colnames(data_selected_melt) <- c("x", "calibration_type", "y")
    
    
    # Plot comparison of matrix-effect coefficients
    slopes_graph <- ggplot(data_selected_melt) +
      aes(x = x, y = y, color = calibration_type) + # Define axis variables
      geom_point() + # Plot points
      geom_smooth(se = T, method = "lm") +
      labs(x = input$x_label, y = input$y_label)+
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
  
  output$x_reg_col <- renderUI({
    req(dataset())
    selectInput("x_reg_col", "Select x variable", names(dataset()))
  })
  output$y_reg_col <- renderUI({
    req(dataset())
    selectInput("y_reg_col", "Select y variable", names(dataset()), selected = "External.standard")
  })
  
  regression_plot <- reactive({
    data <- dataset()
    
    # Fit the linear model
    model <- lm(as.formula(paste(input$y_reg_col, "~", input$x_reg_col)), data = data)
    
    # Create prediction data with prediction intervals
    prediction_data <- cbind(data, predict(model, interval = "prediction"))
    
    # Create the plot with prediction intervals
    plot <- ggplot(prediction_data, aes_string(x = input$x_reg_col, y = input$y_reg_col)) +
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
    data <- dataset()
    model <- lm(as.formula(paste(input$y_reg_col, "~", input$x_reg_col)), data = data)
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
      statistics(dataset(), input$x_reg_col, input$y_reg_col)
    })
      
    regression <- reactive({
      reg <- result()$Regression[[1]]
      reg <- reg %>%
        mutate(across(where(~ !all(is.character(.))), ~ sprintf("%.2e", .))) 
      return(reg)
    })
    
    output$regression <- renderTable({
      regression()
    })
    
    coef.determination <- reactive({
      r2 <- result()$Coef_deter[[1]]
    })
    
    output$coef.determination <- renderTable({
      coef.determination()
    })
    
    levene <- reactive({
      result()$Levene[[1]]
    })
    
    output$levene <- renderTable(levene())
    
    durbin_watson <- reactive({
      result()$Durbin_Watson[[1]]
    })
    
    output$durbin_watson <- renderTable(durbin_watson() )
    
    shapiro_wilk <- reactive({
      result()$Shapiro_Wilk[[1]]
    })
    
    output$shapiro_wilk <- renderTable(shapiro_wilk())
    
    anova <- reactive({
      aov <- result()$Anova[[1]]
      aov <- aov %>%
        mutate(across(where(~ !all(is.character(.))), ~ sprintf("%.2e", .))) 
      return(aov)
    })
    
    output$anova <- renderTable({
      anova()
    })
  
  
    limits_result <- reactive({
      req(dataset(), input$x_reg_col, input$y_reg_col)
      limits(dataset(), input$x_reg_col, input$y_reg_col)
    })
    
   observeEvent(input$limits, {
      
    output$limits <- renderTable(limits_result())
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
      params <- list(dataset = dataset(),
                     outliers_estd = outliers_estd(),
                     outliers_std = outliers_std(),
                     coefficients_overlap = coefficients_overlap(),
                     coefficients_overlap_plot = coefficients_overlap_plot(),
                     regression_plot = regression_plot(),
                     residual_plot = residual_plot(),
                     regression = regression(),
                     coef.determination = coef.determination(),
                     levene = levene(),
                     durbin_watson = durbin_watson (),
                     shapiro_wilk = shapiro_wilk(),
                     anova = anova(),
                     limits_result = limits_result()
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

# Run the application 
shinyApp(ui = ui, server = server)
