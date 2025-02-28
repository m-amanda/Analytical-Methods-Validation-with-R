library(shiny)
library(shinydashboard)
library(rhandsontable)
library(mathjaxr)




ui <- dashboardPage(
  dashboardHeader(title = "Analytical Methods Validation",
                  titleWidth = 300),
  ## Sidebar content
  dashboardSidebar(
    sidebarMenu(
      menuItem("Instructions", tabName = "instructions", icon = icon("circle-info")),
      menuItem("Selectivity", tabName = "Selectivity", icon = icon("filter")),
      menuItem("Linearity", tabName = "Linearity", icon = icon("chart-line")),
      menuItem("Accuracy", tabName = "Accuracy", icon = icon("screenshot", lib = "glyphicon")),
      menuItem("Download Report", tabName = "DownloadReport", icon = icon("download"))
      
    )
  ),
  ## Body content
  dashboardBody(
    tabItems(

      #Instructions Tab
      tabItem(tabName = "instructions",
              h3("Instructions"),
            
              h4("Selectivity"),
              h5("Enter the data to external standard, and standard addition in correspondent place. The first colunm is the concentration and the second the response. The Bonferroni outlier test will be performed automatically. If the Bonferroni p-value is less than 0.05, the corresponding observation could be removed from your dataset, and the analysis will be repeated. The comparison of coefficients will indicate whether the matrix contributes to systematic or constant errors in the analysis. If the linear coefficients of both calibrations do not overlap, there are constant errors, and the method must be optimized to remove them. If the angular coefficients do not overlap, the standard additions method can eliminate matrix effects."),
              
              h4("Linearity"),
              h5("Enter the data to choseen calibration method in correspondent place. The linear regression will display confidence and prediction intervals. At a minimum, the observations must fall within the prediction interval. The residual plot should not display trends, such as a curve or a horizontal triangle pattern, and the residuals should cluster around 0. If there are trends, another type of regression may be more appropriate."),
              
              h5("Several statistical tests are performed to assess linearity. Your data should be normally distributed (Shapiro-Wilk p-value > 0.05), without autocorrelated errors (Durbin-Watson p-value > 0.05), and show homogeneity of variance (Levene's test p-value > 0.05). Additionally, the regression should be statistically significant (ANOVA regression p-value < 0.05) and show no lack of fit (ANOVA lack-of-fit p-value > 0.05)."),
              
              h5("If your calibration fails any of these tests, consider using a different regression model or recalibrating the method."),
              
              h4("LOD and LOQ"),
              h5("The limits of detection (LOD) and quantification (LOQ) are calculated based on the data used in the linearity test. These limits are derived from the standard error of the regression and the slope (angular coefficient). The LOD and LOQ should be equal to or below the first concentration in your working range."),
              
              # LaTeX formulas for LOD and LOQ
              withMathJax(
                helpText(
                  "$$ LOD = 3.3s/b $$",
                  "$$ LOQ = 10s/b $$"
                )

              ),
              
              
              
              
              
              h4("Precision"),
              h5("There are two methods to estimate precision: the external standard method and the standard additions method. Data should be collected from at least six replicates on the first day to calculate repeatability, followed by data collection on three additional, distinct days. Data should be collected at Q1, Q2, and Q3 points of the calibration curve for fortified solutions, and from a non-fortified sample to determine the original concentration and calculate trueness. The standard deviation is calculated using ANOVA, as recommended by EURACHEM guidelines. For the standard additions method, the standard deviation is calculated using the reverse axis method, as described in the literature by Gonçalves, Jones, and Donati."),
              
              h4("Trueness"),
              h5("Trueness is calculated using the recovery method with the same data collected for precision. The recovery method is based on the difference in predicted concentrations between the fortified sample and the non-fortified sample, which is then divided by the reference concentration. The result is multiplied by 100% to express recovery as a percentage."),
              
              h4("Download Report"),
              h5("A report with results is available for download in PDF, DOCX, or HTML format."),
              
              h4("References"),
              h5(HTML(
                "INTERNATIONAL COUNCIL FOR HARMONISATION OF TECHNICAL REQUIREMENTS FOR PHARMACEUTICALS FOR HUMAN USE.
               ICH Q2(R2) Validation of Analytical Procedures, 2022.
               <a href='htt           #     INMETRO. Orientação Sobre Validação De Métodos Analíticos, 2020.
               <a href='http://www.inmetro.gov.br/Sidoq/Arquivos/Cgcre/DOQ/DOQ-Cgcre-8_08.pdf' target='_blank'>http://www.inmetro.gov.br/Sidoq/Arquivos/Cgcre/DOQ/DOQ-Cgcre-8_08.pdf</a>.<br><br>

               B. Magnusson; U. Örnemark. Eurachem Guide: The Fitness for Purpose of Analytical Methods – A Laboratory Guide to Method Validation and Related Topics, 2nd ed.; 2014. <a href='www.eurachem.org' target='_blank'>www.eurachem.org</a>.<br><br>

               Francisco Raposo; Damià Barceló. Assessment of Goodness-of-Fit for the Main Analytical Calibration Models: Guidelines and Case Studies.
             <a href='https://doi.org/10.1016/j.trac.2021.116373' target='_blank'>https://doi.org/10.1016/j.trac.2021.116373</a>.<br><br>
                
                Daniel A. Goncalves, Bradley T. Jones, George L. Donati,
                The reversed-axis method to estimate precision in standard additions analysis,
                Microchemical Journal, Volume 124, 2016, Pages 155-158, ISSN 0026-265X, <a href='https://doi.org/10.1016/j.microc.2015.08.006' target='_blank'>https://doi.org/10.1016/j.microc.2015.08.006</a>.<br><br>

                
                "
                
                
              ))
      ),

   
      
      # Selectivity Tab
      tabItem(tabName ="Selectivity",
            
                  fluidRow(
                    box(width=4,
                      h4("Enter External Standard data"),
                      rHandsontableOutput("table_estd"),
                      h4("Bonferroni outliers test to external standard"),
                      verbatimTextOutput("outliers_estd"),
                      uiOutput("p_ouliers_estd")
                            
                    ),
                    
                    box(width=4,
                      h4("Enter Standard Additions data"),
                      rHandsontableOutput("table_std"),
                      h4("Bonferroni outliers test to standard additions"),
                      verbatimTextOutput("outliers_std"),
                      uiOutput("p_ouliers_std")
                      

                    ),
                    
                    box(width=4,
                      h4("Enter Youden calibration data"),
                      rHandsontableOutput("table_youden")
                      
                      
                    )
                  ),

              fluidRow(
                box(            
                  h4("Plot overlap"),
                  h5("The linear regression performed to different methods plotted together"),
                  plotOutput("coefficients_overlap_plot", width = "70%"),
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
                
                box(
                  h4("Coefficiet overlap"),
                  h5("Coefficients and Confidence intervals"),
                  h5("External Standard"),
                  verbatimTextOutput("coefficients_estd"),
                  h5("Standard Additions"),
                  verbatimTextOutput("coefficients_std"),
                  h5("Youden"),
                  verbatimTextOutput("coefficients_youden"),
                  uiOutput("p_slopes_compare"),
                  uiOutput("p_intercepts_compare")
                )
            )
            
      ),
        
      # Linearity Tab
      tabItem(tabName = "Linearity",
          fluidRow(
            box(
                h4("Enter calibration data"),
                rHandsontableOutput("table_linearity"),
                h4("Linear Regression"),
                h5("The regression line in blue, confidence band in gray and prediction band in dashed black line"),
                plotOutput("regression_plot", height = "4in"),
                h4("Residual graph"),
                plotOutput("residual_plot", height = "4in")
            ),
            box(
                h4("Coefficient of determination"),
                tableOutput("coef.determination"),
                h4("Regressions terms"),
                tableOutput("regression"),
                h4("ANOVA of regression and lack-of-fit"),
                tableOutput("anova"),
                uiOutput("p_anova_reg"),
                uiOutput("p_anova_lof")
          ),
          box(
            h4("Levene homogeneity test"),
            tableOutput("levene"),
            uiOutput("p_levene")
          ),
          
          box(
            h4("Durbin-Watson test to auto-cooreleted errors"),
            tableOutput("durbin_watson"),
            uiOutput("p_dw")
          ),
          
          box(
            h4("Shapiro-Wilk normality test"), 
            tableOutput("shapiro_wilk"),
            uiOutput("p_normality")
        ),
        box(
          h4("LOD and LOQ"),
          tableOutput("limits"),
          h5("The values are in terms of concentration")
        )
          )
      ),
      
      tabItem(tabName = "Accuracy",
              h4("Choose the method to estimate precision and trueness"),
              actionButton("acc_estd", "Extenal Standard"),
              actionButton("acc_std", "Standard Additions"),
              br(), br(),
              uiOutput("acc_results"),
             
              ),
      
      tabItem(tabName = "DownloadReport",
              
                radioButtons('format', 'Document format', c('PDF', 'HTML', 'Word'),
                             inline = TRUE),
                downloadButton("report", "Generate report")
              
      )
      
          
      
      
    )
  )
)



      

