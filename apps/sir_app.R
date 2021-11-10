library(shiny)
library(deSolve)

# Data
create_setup <- function(T_max, by, I0) {
  t <- as.numeric(seq(0.0, T_max, by = by))
  N <- length(t)
  data_list <- list(
    N = N,
    t = t,
    I0 = I0,
    pop_size = 1000,
    I_data = rep(1.0, N) # dummy
  )
  return(data_list)
}

# ODE system
sir_rhs <- function(t, y, p) {
  infection_rate <- p[1] * y[2] * y[1] / p[3]
  recovery_rate <- p[2] * y[2]
  d1 <- -infection_rate
  d2 <- infection_rate - recovery_rate
  out <- list(c(d1, d2))
  return(out)
}

# Define UI for miles per gallon app ----
ui <- pageWithSidebar(

  # App title ----
  headerPanel("SIR app"),

  # Sidebar panel for inputs ----
  sidebarPanel(
    sliderInput(
      "T_max",
      "Length of time interval",
      15,
      min = 1,
      max = 30,
    ),
    sliderInput(
      "by",
      "by",
      0.1,
      min = 0.01,
      max = 1,
    ),
    sliderInput(
      "I0",
      "Initial number of infected",
      3,
      min = 1,
      max = 30,
      step = 1
    ),
    selectInput("solver", "Solver:", c("RK45" = "ode45", "BDF" = "bdf")),
    sliderInput(
      "beta",
      "beta (infection rate)",
      3.0,
      min = 0.01,
      max = 10,
      step = 0.01
    ),
    sliderInput(
      "gamma",
      "gamma (recovery rate)",
      0.3,
      min = 0.01,
      max = 5,
      step = 0.01
    ),
  ),

  # Main panel for displaying outputs ----
  mainPanel(

    # Output: Plot of the requested variable against mpg ----
    plotOutput("mainPlot"),
    div(
      class = "footer",
      includeHTML("footer.html")
    )
  )
)


# Define server logic to plot various variables against mpg ----
server <- function(input, output) {
  output$mainPlot <- renderPlot({
    dat <- create_setup(input$T_max, input$by, input$I0)
    y_sol <- ode(
      y = c(dat$pop_size - dat$I0, dat$I0),
      times = dat$t,
      func = sir_rhs,
      parms = c(input$beta, input$gamma, dat$pop_size),
      method = input$solver
    )
    plot(y_sol[, 1], y_sol[, 3], "l",
      pch = 20, ylim = c(0, 1.1 * dat$pop_size),
      main = "Infected", xlab = "t", ylab = "I(t)",
      col = "firebrick3", lwd = 2
    )
    lines(c(-10, 1.1 * input$T_max), c(dat$pop_size, dat$pop_size),
      lty = 2,
      col = "steelblue"
    )
  })
}

shinyApp(ui, server)
