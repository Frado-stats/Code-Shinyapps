library(shiny)
library(MASS)
library(ggplot2)
library(dplyr)
library(tidyr)
library(shiny)
###########################################################################
probX2givenX1 <- function(X_values, thres_v2, rho){
  musY_X <- rho * X_values
  varY_X <- 1 - rho^2
  1 - pnorm(thres_v2, mean = musY_X, sd = sqrt(varY_X))
}
solve_rho <- function(n, frac_pos, target_cor,
                      mean_vec = c(0, 0),
                      interval = c(-0.99, 0.99),
                      tol = 1e-6) {
  
  f <- function(rho) {
    Sigma <- matrix(c(1, rho, rho, 1), 2)
    d <- as.data.frame(mvrnorm(n, mu = mean_vec, Sigma = Sigma))
    colnames(d) <- c("V1", "V2")
    
    cutoff <- sort(d$V2)[floor(n * (1 - frac_pos))]
    d$V2 <- as.numeric(d$V2 > cutoff)
    
    cor(d$V1, d$V2)
  }
  
  uniroot(function(rho) f(rho) - target_cor,
          interval = interval)$root
}
solve_b0 <- function(d, b1, b2, target_prev) {
  
  f <- function(b0) {
    lp <- b0 + log(b1)*d$V1 + log(b2)*d$V2
    mean(plogis(lp)) - target_prev
  }
  
  uniroot(f, c(-20, 20))$root
}
###########################################################################
ui <- fluidPage(
  
  titlePanel("A tool to explore the relations between Pre-test Probabilities and Sensitivity and Specificity"),
  
  tags$div(
    style = "padding: 10px 0px 15px 0px; color: #555;",
    "This tool simulates post-test probabilities using a logistic model with a continuous standard test (V1) and a dichotomised biomarker (V2).
   Pre-test probabilities are obtained by marginalising over V2. Likelihood ratios are then derived from pre- and post-test probabilities
   (rearranging: post-test odds = pre-test odds * LR), and subsequently translated into sensitivity and specificity 
    (inversing LR+ = sens / (1-spec), LR- = (1-sens) / spec)."
  ),
  sidebarLayout(
    sidebarPanel(
      sliderInput("corr", "Correlation between V1 and V2", 0, 0.9, 0.1, step = 0.05),
      sliderInput("event", "Prevalence", 0.05, 0.5, 0.2, step = 0.01),
      sliderInput("b1", "Odds ratio V1", 1.0, 3, 1.64, step = 0.05),
      sliderInput("b2", "Odds ratio V2", 1.0, 3, 2.13, step = 0.05),
      sliderInput("frac_pos", "Proportion positive results V2", 0.1, 0.9, 0.4, step = 0.01)),
    mainPanel(plotOutput("plot", height = "650px"))
  )
)
###########################################################################
server <- function(input, output) {
  
  sim_data <- reactive({
    
    set.seed(4)
    
    n <- 50000
    
    rho <- solve_rho(
      n = n,
      frac_pos = input$frac_pos,
      target_cor = input$corr
    )
    
    R <- matrix(c(1, rho, rho, 1), 2)
    
    d <- as.data.frame(mvrnorm(n, mu = c(0, 0), Sigma = R))
    colnames(d) <- c("V1", "V2")
    
    cutoff <- sort(d$V2)[floor(n * (1 - input$frac_pos))]
    d$V2 <- as.numeric(d$V2 > cutoff)
    
    d$probV2 <- probX2givenX1(d$V1, cutoff, rho)
    
    b0 <- solve_b0(d, input$b1, input$b2, input$event)
    
    plow  <- plogis(b0 + log(input$b1)*d$V1)
    phigh <- plogis(b0 + log(input$b1)*d$V1 + log(input$b2))
    
    d$riskref_true <- plow * (1 - d$probV2) + phigh * d$probV2
    
    odds_base <- d$riskref_true / (1 - d$riskref_true)
    odds_high <- phigh / (1 - phigh)
    odds_low  <- plow / (1 - plow)
    
    LR_pos <- odds_high / odds_base
    LR_neg <- odds_low / odds_base
    
    d$sens <- (LR_pos * (1 - LR_neg)) / (LR_pos - LR_neg)
    d$spec <- (1 - d$sens) / LR_neg
    
    d
  })
  
  output$plot <- renderPlot({
    
    d <- sim_data()
    
    plot_dat <- d |>
      select(riskref_true, sens, spec, probV2) |>
      pivot_longer(cols = c(sens, spec, probV2),
                   names_to = "metric",
                   values_to = "value") |>
      mutate(metric = recode(metric,
                             probV2 = "prob_pos"))
    
    ggplot(plot_dat,
           aes(riskref_true, value, colour = metric)) +
      geom_line(linewidth = 1) +
      ylim(0, 1) +
      theme_minimal(base_size = 16) +
      scale_colour_manual(
        values = c(
          "prob_pos" = "red",
          "sens" = "green",
          "spec" = "blue"
        ),
        labels = c(
          "prob_pos" = "Probability of positive test result",
          "sens" = "Sensitivity",
          "spec" = "Specificity"
        )
      ) +
      theme(
        legend.position = c(0.98, 0.98),
        legend.justification = c(1, 1),
        axis.title = element_text(size = 18),
        axis.text  = element_text(size = 14),
        legend.text = element_text(size = 14)
      ) +
      labs(
        x = "Pre-test probability",
        y = "Value",
        colour = NULL
      )
  })
}
###########################################################################
shinyApp(ui, server)