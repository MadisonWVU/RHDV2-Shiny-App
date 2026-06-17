#################### Make Shiny App to Predict RHDV2 in WV ####################
library(shiny)
library(ggplot2)
library(sf)
library(spdep)

# setwd("C:/Users/madis/Desktop/RHDV2_R/WV_predictions")

#load county data
counties <- read_sf("WV_counties.shp")
county_names <- sort(counties$NAME)

#neighborhood matrix
nb <- poly2nb(counties, queen = TRUE, snap = 1e-6)
W  <- nb2mat(nb, zero.policy = TRUE, style = "C")
n  <- nrow(W)

#CAR parameters from Miller et al. (2026)
rho  <- 0.9897
tau2 <- 1.1750

pm    <- rho * (diag(rowSums(W)) - W) + (1 - rho) * diag(n)
sigma <- tau2 * solve(pm)

################################### UI ########################################
#fluid page is container for user interface
ui <- fluidPage(
  titlePanel("RHDV2 Case Prediction in West Virginia"),
  
  #split into map on right and input dropdowns on left
  sidebarLayout(
    sidebarPanel(
      h4("Observed Cases"),
      #adds and removes rows of interest
      uiOutput("county_entries"),
      actionButton("add_row",    "Add County"),
      actionButton("remove_row", "Remove Last"),
      hr(),
      #run the calculation and map
      actionButton("run", "Run Prediction", class = "btn-primary"),
      hr(),
      #will show intercept value and total cases
      verbatimTextOutput("summary_text")
    ),
    mainPanel(
      plotOutput("map_plot", height = "600px")
    )
  )
)

######################## Server ###############################################
server <- function(input, output, session) {
  
  #how many counties are showing starting at 1
  n_rows <- reactiveVal(1)
  
  #changes when click buttons to add counties, maxed at 55
  observeEvent(input$add_row,    { n_rows(min(n_rows() + 1, 55)) })
  observeEvent(input$remove_row, { n_rows(max(n_rows() - 1, 1))  })
  
  #fills in county_entries from UI, changes and reruns when n_rows changes
  output$county_entries <- renderUI({
    lapply(seq_len(n_rows()), function(i) {
      fluidRow(
        #county drop down list placed in shiny grid cell #7
        column(7, selectInput(paste0("county_", i),
                              label = if (i == 1) "County" else NULL,
                              choices = c("-- select --" = "", county_names),
                              selected = isolate(input[[paste0("county_", i)]]))),
        #tracks cases placed in shiny grid cell #5
        column(5, numericInput(paste0("cases_", i),
                               label = if (i == 1) "Cases" else NULL,
                               #preserves current selection and adds to it
                               value = isolate(input[[paste0("cases_", i)]]) %||% 1,
                               min = 1, step = 1))
      )
    })
  })

  #collects all county/cases and inputs into a single data frame   
  get_obs <- reactive({
    input$run
    isolate({
      #loops through each row and reads case input
      rows <- lapply(seq_len(n_rows()), function(i) {
        co <- input[[paste0("county_", i)]]
        ca <- input[[paste0("cases_",  i)]]
        #checks that a county and case number are selected        
         if (!is.null(co) && nzchar(co) && !is.null(ca) && !is.na(ca) && ca >= 1)
          data.frame(county = co, cases = as.numeric(ca), stringsAsFactors = FALSE)
      })
      #drops rows user left blank, keeps remaining into a df
      do.call(rbind, Filter(Negate(is.null), rows))
    })
  })
 
  #fills in summary text for UI including total cases and intercept 
  output$summary_text <- renderPrint({
    obs <- get_obs()
    
    #placeholder message if no cases have been added yet
    if (is.null(obs) || nrow(obs) == 0) { cat("No observations yet."); return() }
    
    #prevents same county from being entered twice, collapses them
    obs <- aggregate(cases ~ county, data = obs, FUN = sum)
    total <- sum(obs$cases)
    cat("Observed counties:\n")
    
    #prints the added cases per county 
    for (i in seq_len(nrow(obs)))
      cat(sprintf("  %s: %d case(s)\n", obs$county[i], obs$cases[i]))
    cat(sprintf("\nTotal cases:     %d\n", total))
    cat(sprintf("Intercept (log): %.4f\n", log(total / n)))
  })
  
  #placeholder for the map in the UI
  output$map_plot <- renderPlot({
    input$run
    isolate({
      obs <- get_obs()
      
      if (is.null(obs) || nrow(obs) == 0) {
        plot(st_geometry(counties), main = "Select counties and click Run Prediction")
        return()
      }
      
      #aggregate duplicates by summing cases
      obs <- aggregate(cases ~ county, data = obs, FUN = sum)
      total_cases <- sum(obs$cases)
      
      #match observation to county shapefile by row indices
      obs_idx <- match(obs$county, counties$NAME)
      
      #if county name doesn't match, return NA
      obs     <- obs[!is.na(obs_idx), ]
      obs_idx <- obs_idx[!is.na(obs_idx)]
      
      #intercept total cases divided by total counties on the log scale
      xb_i <- rep(log(total_cases / n), n)
      
      #create vector of 55 zeros to hold predicted intensity for every county
      ir <- numeric(n)
      
      #fills in positive counties with actual log case counts
      ir[obs_idx] <- log(obs$cases)
      
      #gets index for counties with no cases
      unobs_idx <- setdiff(seq_len(n), obs_idx)
      
      #use covariance from sigma and matrix math to determine how far spread
      ir[unobs_idx] <-
        xb_i[unobs_idx] +
        sigma[unobs_idx, obs_idx, drop = FALSE] %*%
        solve(sigma[obs_idx, obs_idx, drop = FALSE]) %*%
        matrix(ir[obs_idx] - xb_i[obs_idx], ncol = 1)
      
      #exponentiate everything back to count scale
      counties$predicted <- exp(ir)
      counties$observed  <- seq_len(n) %in% obs_idx
      
      #plot with user input
      ggplot(counties) +
        geom_sf(aes(fill = predicted), color = "grey30", linewidth = 0.3) +
        geom_sf(data = counties[counties$observed, ],
                fill = NA, color = "red", linewidth = 0.8) +
        geom_sf_text(aes(label = round(predicted, 2)), size = 5, color = "black") +
        scale_fill_distiller(name = "Predicted\ncases", palette = "YlOrRd", direction = 1) +
        theme_void() +
        theme(legend.title = element_text(face= "bold", size = 16),
              legend.text  = element_text(size = 14)) +
        labs(title = sprintf("Total observed: %d  |  Log intercept: %.3f",
                             total_cases, log(total_cases / n)))
    })
  })
}

`%||%` <- function(a, b) if (!is.null(a)) a else b

shinyApp(ui, server)
