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

statistics <- funstatistics <- function(data, x, y) {
  data %>%
    reframe(
      model = list(lm(!!sym(y) ~ !!sym(x), data = cur_data())),
      Regression = list(model[[1]] %>% tidy()), # Regression model
      Levene = list(leveneTest(!!sym(y) ~ as.factor(!!sym(x))) %>% tidy()), # Homogeneity test
      Durbin_Watson = list(durbinWatsonTest(model[[1]], max.lag = 1) %>% tidy()), # Autocorrelation test
      Shapiro_Wilk = list(shapiro.test(model[[1]]$residuals) %>% tidy()), # Normality test
      Anova = list(anova(model[[1]]) %>% tidy()) # ANOVA test
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

ui <- navbarPage(
  theme = shinytheme("sandstone"),
  
  title = "Analytical Calibration Statistics",
  
  tabPanel(
    "Regression Analysis",
    sidebarLayout(
      sidebarPanel(
        fileInput("dataset_file", "Choose CSV file", accept = ".csv"),
        
        uiOutput("xcol"),
        uiOutput("ycol"),
        
        selectInput("theme", "Plot Theme:", choices = list(
          "theme_grey" = "theme_grey",
          "theme_bw" = "theme_bw",
          "theme_classic" = "theme_classic",
          "theme_minimal" = "theme_minimal",
          "theme_void" = "theme_void",
          "theme_economist" = "theme_economist",
          "theme_fivethirtyeight" = "theme_fivethirtyeight",
          "theme_solarized" = "theme_solarized",
          "theme_stata" = "theme_stata",
          "theme_wsj" = "theme_wsj"
        ), selected = "theme_classic"),
        
        textInput("x_label", "X-axis label:", value = ""),
        textInput("y_label", "Y-axis label:", value = ""),
        
        
      ),
      mainPanel(
        h3("Regression Plot"),
        plotOutput("regression_plot", width = "5in", height = "4in"),  # Set dimensions for the plot
        plotOutput("residual_plot", width = "8in", height = "8in")  
      )
    )
  ),
  
  tabPanel(
    "Adjust-of-fit",
    mainPanel(
      h3("Results"),
      actionButton("analyze","Calculate adjust-of-fit statistics"),
      h4("Regressions terms"),
      tableOutput("regression"),
      h4("Levene test to homogeinedy"),
      tableOutput("levene"),
      h4("Durbin-Watson test to autocooreleted errors"),
      tableOutput("durbin_watson"),
      h4("Shapiro-Wilk normality test"),
      tableOutput("shapiro_wilk"),
      h4("ANOVA of regression"),
      tableOutput("anova")
    )
  ),
    tabPanel(
      "LOD and LOQ",
      mainPanel(
        h3("LOD and LOQ with calibration curve methods"),
        actionButton("limits","Calculate LOD and LOQ"),
        tableOutput("limits"),
      )
    )
)



server <- function(input, output, session) {
  # Função para ler o CSV enviado pelo usuário
  dataset <- reactive({
    req(input$dataset_file)
    
    # Lê o arquivo CSV e remove nomes de linhas duplicados
    data <- read.csv(input$dataset_file$datapath, header = TRUE, sep = ",", row.names = NULL)
    
    # Garante que os nomes das linhas sejam únicos
    rownames(data) <- make.unique(as.character(seq_len(nrow(data))))
    
    return(data)
  })
  
  # Popula as opções para seleção de X e Y com base nos dados do dataset carregado
  output$xcol <- renderUI({
    req(dataset())
    selectInput("xcol", "Select x variable", names(dataset()))
  })
  
  output$ycol <- renderUI({
    req(dataset())
    selectInput("ycol", "Select y variable", names(dataset()))
  })
  
  # Realiza as análises quando o botão "Analisar" for clicado
  observeEvent(input$analyze, {
    req(input$xcol, input$ycol)
    
    result <- statistics(dataset(), input$xcol, input$ycol)
    
    output$regression <- renderTable(result$Regression[[1]])
    output$levene <- renderTable(result$Levene[[1]])
    output$durbin_watson <- renderTable(result$Durbin_Watson[[1]])
    output$shapiro_wilk <- renderTable(result$Shapiro_Wilk[[1]])
    output$anova <- renderTable(result$Anova[[1]])
  })
  
  observeEvent(input$limits, {
    req(input$xcol, input$ycol)
    
    limits_result <- limits(dataset(), input$xcol, input$ycol)
    output$limits <- renderTable(limits_result)
  })
  
  # Regression plot with prediction intervals
  output$regression_plot <- renderPlot({
    data <- dataset()
    
    # Fit the linear model
    model <- lm(as.formula(paste(input$ycol, "~", input$xcol)), data = data)
    
    # Create prediction data with prediction intervals
    prediction_data <- cbind(data, predict(model, interval = "prediction"))
    
    # Create the plot with prediction intervals
    plot <- ggplot(prediction_data, aes_string(x = input$xcol, y = input$ycol)) +
      geom_point() +  # Plot the points with selected color
      geom_smooth(method = "lm", se = TRUE) +  # Regression line with CI
      geom_line(aes(y = fit, color = "Regression Line"), size = 1) +  # Regression line
      geom_line(aes(y = upr, color = "Prediction Interval", linetype = "PI"), size = 1) +  # Upper prediction interval
      geom_line(aes(y = lwr, color = "Prediction Interval", linetype = "PI"), size = 1) +  # Lower prediction interval
      scale_color_manual(values = c("Prediction Interval" = "black")) +  # Custom colors
      scale_linetype_manual(values = c("PI" = "dashed")) +  # Dashed line for PI
      labs(x = input$x_label, y = input$y_label, title = "Regression Plot with Prediction Intervals")
    
    # Apply the user-selected theme
    selected_theme <- switch(input$theme,
                             "theme_grey" = theme_grey(),
                             "theme_bw" = theme_bw(),
                             "theme_classic" = theme_classic(),
                             "theme_minimal" = theme_minimal(),
                             "theme_void" = theme_void(),
                             "theme_economist" = theme_economist(),
                             "theme_fivethirtyeight" = theme_fivethirtyeight(),
                             "theme_solarized" = theme_solarized(),
                             "theme_stata" = theme_stata(),
                             "theme_wsj" = theme_wsj()
    )
    
    # Add the theme to the plot
    plot <- plot + selected_theme
    
    print(plot)
  })
  
  output$residual_plot <- renderPlot({
    data <- dataset()
    model <- lm(as.formula(paste(input$ycol, "~", input$xcol)), data = data)
    r1<-par(mfrow = c(2, 2))
    plot(model)
    r1  # Return the plot object
  })
}


shinyApp(ui = ui, server = server)
