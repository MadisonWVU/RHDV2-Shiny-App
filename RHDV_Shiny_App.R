#################### Make Shiny App to Predict RHDV2 in WV ####################
library(shiny)
library(ggplot2)
library(sf)
library(spdep)
library(plotly)

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

#shared color palette so map fill and the legend match
ylorrd_pal <- colorRampPalette(c("#FFFFB2","#FECC5C","#FD8D3C","#F03B20","#BD0026"))

#converts predicted values to hex colors using YlOrRd palette
val_to_color <- function(vals) {
  breaks <- seq(min(vals, na.rm = TRUE), max(vals, na.rm = TRUE), length.out = 100)
  idx    <- findInterval(vals, breaks, all.inside = TRUE)
  ylorrd_pal(100)[idx]
}

#builds the plotly map from county sf object
#fill_vals: hex color per county, observed_idx: row indices of observed counties
make_plotly_map <- function(counties_sf, fill_vals, observed_idx, title_text) {
  
  fig <- plot_ly()
  
  #loops through each county and adds a filled polygon trace
  for (i in seq_len(nrow(counties_sf))) {
    geom   <- st_geometry(counties_sf)[[i]]
    coords <- st_coordinates(geom)
    xs     <- coords[, 1]
    ys     <- coords[, 2]
    
    fill_color <- fill_vals[i]
    
    #observed counties get red border, unobserved get grey
    line_color <- if (i %in% observed_idx) "red"  else "grey30"
    line_width <- if (i %in% observed_idx) 2      else 0.5
    
    #hover always shows county name regardless of whether prediction has been run
    fig <- add_trace(fig,
                     type       = "scatter",
                     mode       = "lines",
                     x          = c(xs, xs[1]),
                     y          = c(ys, ys[1]),
                     fill       = "toself",
                     fillcolor  = fill_color,
                     line       = list(color = line_color, width = line_width),
                     name       = counties_sf$NAME[i],
                     hoverinfo  = "name",
                     showlegend = FALSE
    )
  }
  
  #get county centroids for placing predicted count labels
  centroids <- st_coordinates(st_centroid(counties_sf))
  
  has_predictions <- !is.null(counties_sf$predicted) && any(!is.na(counties_sf$predicted))
  
  #only show predicted count labels if predictions have been run
  fig <- add_annotations(fig,
                         x         = centroids[, 1],
                         y         = centroids[, 2],
                         text      = if (has_predictions)
                           as.character(round(counties_sf$predicted, 2))
                         else
                           rep("", nrow(counties_sf)),
                         showarrow = FALSE,
                         font      = list(size = 11, color = "black"),
                         hoverinfo = "skip"
  )
  
  #once predictions exist, add an invisible "dummy" trace for the legend
  if (has_predictions) {
    rng <- range(counties_sf$predicted, na.rm = TRUE)
    
    #build a plotly-style colorscale from the same palette used for the fill
    stops      <- seq(0, 1, length.out = 100)
    pal_colors <- ylorrd_pal(100)
    
    fig <- add_trace(fig,
                     type       = "scatter",
                     mode       = "markers",
                     x          = centroids[1, 1],
                     y          = centroids[1, 2],
                     marker     = list(
                       color      = rng,
                       colorscale = Map(function(s, c) list(s, c), stops, pal_colors),
                       cmin       = rng[1],
                       cmax       = rng[2],
                       showscale  = TRUE,
                       size       = 0.0001,
                       opacity    = 0,
                       colorbar   = list(
                         title = list(text = "Predicted\nCases", side = "right"),
                         len   = 0.6,
                         thickness = 18
                       )
                     ),
                     showlegend = FALSE,
                     hoverinfo  = "none"
    )
  }
  
  layout(fig,
         title         = list(text = title_text, font = list(size = 16)),
         xaxis         = list(visible = FALSE),
         yaxis         = list(visible = FALSE, scaleanchor = "x"),
         margin        = list(l = 0, r = 0, t = 40, b = 0),
         paper_bgcolor = "white",
         plot_bgcolor  = "white"
  )
}

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
      #will show total cases
      verbatimTextOutput("summary_text")
    ),
    mainPanel(
      plotlyOutput("map_plot", height = "600px")
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
  
  #fills in summary text for UI including total cases
  output$summary_text <- renderPrint({
    obs <- get_obs()
    
    #placeholder message if no cases have been added yet
    if (is.null(obs) || nrow(obs) == 0) { cat("No observations yet."); return() }
    
    #prevents same county from being entered twice, collapses them
    obs   <- aggregate(cases ~ county, data = obs, FUN = sum)
    total <- sum(obs$cases)
    cat("Observed counties:\n")
    
    #prints the added cases per county
    for (i in seq_len(nrow(obs)))
      cat(sprintf("  %s: %d case(s)\n", obs$county[i], obs$cases[i]))
    cat(sprintf("\nTotal cases: %d\n", total))
  })
  
  #builds and renders the plotly map
  output$map_plot <- renderPlotly({
    input$run
    isolate({
      obs <- get_obs()
      
      #show empty grey map with hover names before any input
      if (is.null(obs) || nrow(obs) == 0) {
        counties$predicted <- NA
        fill_colors <- rep("lightgrey", nrow(counties))
        return(make_plotly_map(counties, fill_colors, integer(0),
                               "Select counties and click Run Prediction"))
      }
      
      #aggregate duplicates by summing cases
      obs         <- aggregate(cases ~ county, data = obs, FUN = sum)
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
      
      #map predicted values to YlOrRd color scale
      fill_colors <- val_to_color(counties$predicted)
      make_plotly_map(counties, fill_colors, obs_idx,
                      sprintf("Total observed: %d", total_cases))
    })
  })
}

`%||%` <- function(a, b) if (!is.null(a)) a else b

shinyApp(ui, server)
