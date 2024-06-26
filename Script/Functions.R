
########### Functions ############
##################################

#------------- Get HYCOM data

#The funcion builds a character string pointing to the url where certain 
#variables (currents u and v for a certain day) are stored. This is done by 
#concatenating the data for the day with the common parts of the url that are 
#fixed between different urls. Data on the url is downloaded at the end of the 
#script with the function 'curl_download()', and stored as a netcdf file on disk.



get.hycom <- function(limits, time, vars=c('water_u', 'water_v'), include_latlon=TRUE,
                      filename='', download.file=TRUE, dir = getwd(), depLevels=NULL) {
  
  ## Set the base URL based on the start date. If the ending date exceeds the
  ## period for this experiment, then print a warning and truncate the output
  ## early.
  expts = data.frame(
    start=c(as.Date('1992-10-02'), as.Date('1995-08-01'),
            as.Date('2013-01-01'), as.Date('2013-08-21'),
            as.Date('2014-04-05'), as.Date('2016-04-18')),
    end=c(as.Date('1995-07-31'), as.Date('2012-12-31'),
          as.Date('2013-08-20'), as.Date('2014-04-04'),
          as.Date('2016-04-17'), Sys.Date() + 1),
    url=c('http://ncss.hycom.org/thredds/ncss/GLBu0.08/expt_19.0/',
          'http://ncss.hycom.org/thredds/ncss/GLBu0.08/expt_19.1/',
          'http://ncss.hycom.org/thredds/ncss/GLBu0.08/expt_90.9?',
          'http://ncss.hycom.org/thredds/ncss/GLBu0.08/expt_91.0?',
          'http://ncss.hycom.org/thredds/ncss/GLBu0.08/expt_91.1?',
          'http://ncss.hycom.org/thredds/ncss/GLBu0.08/expt_91.2?'))
  
  if(time[1] < expts$start[1])
    stop('Data begins at %s and is not available at %s.',
         strftime(expts$start[1], '%d %b %Y'),
         strftime(time[1], '%d %b %Y'))
  if(time[1] > expts$end[nrow(expts)])
    stop('Data ends at %s and is not available at %s.',
         strftime(expts$end[nrow(expts)], '%d %b %Y'),
         strftime(time[1], '%d %b %Y'))
  for(i in seq(nrow(expts))) {
    if((time[1] >= expts$start[i]) & (time[1] <= expts$end[i]))
      url = expts$url[i]
  }
  
  if(any(grep('19', url))) url = sprintf('%s%s?', url, as.numeric(format(time, '%Y')))
  
  
  
  ## Add the variables.
  for(var in vars)
    url = sprintf('%svar=%s&', url, var)
  
  
  
  ## Add the spatial domain.
  url = sprintf('%snorth=%f&west=%f&east=%f&south=%f&horizStride=1&',
                url, limits[[4]], limits[[1]], limits[[2]], limits[[3]])
  # north, west, east, south
  
  
  
  ## Add the time domain.
  if(length(time) == 2){
    url = sprintf('%stime_start=%s%%3A00%%3A00Z&time_end=%s%%3A00%%3A00Z&timeStride=1&',
                  url, strftime(time[1], '%Y-%m-%dT00'),
                  strftime(time[2], '%Y-%m-%dT00'))
  } else if(length(time) == 1){
    url = sprintf('%stime_start=%s%%3A00%%3A00Z&time_end=%s%%3A00%%3A00Z&timeStride=1&',
                  url, strftime(time[1], '%Y-%m-%dT00'),
                  strftime(time[1], '%Y-%m-%dT00'))
  }
  
  
  
  ## Add the lat-lon points if requested.
  if(include_latlon)
    url = sprintf('%saddLatLon=true&', url)
  
  
  
  ## Finish the URL.
  if (is.null(depLevels)){
    url = sprintf('%sdisableProjSubset=on&vertCoord=&accept=netcdf', url)
  } else{
    url = paste(url,'disableProjSubset=on&vertCoord=', depLevels, '&accept=netcdf', sep='')
  }
  
  print(url)
  
  
  
  ## Download the data if a filename was provided.
  if(filename != ''){
    if(download.file == TRUE){
      #download.file(url, filename, method = 'auto')
      curl_download(url, filename, quiet=FALSE)
    } else if(download.file == FALSE){
      system(sprintf('curl -o "%s" "%s"', filename, url))
    }
  }
  return(url)
}





#-------------  Calculate cost matrix from HYCOM data


calculate_costMatrix <- function(files, geo_points) {
  
  
  cost_list <- list()
  
  for(i in 1:length(files)){
    
    
    
    cat("File ", i, "/n")  # to keep track of analyzed files
    
    #Read u and v from the netcdf raster (.nc extension)
    
    # nc_open from the netCDF package allows for getting the names of bands on 
    #other specific properties of netCDF that are lacking in the raster package
    
    nc <- nc_open(files[i])
    if(length(which(abs(nc$var$water_u$dim[[1]]$vals) > 180)) >= 1){  #depending 
      #on their coords, some rasters need to be rotated
      rr_u <- rotate(raster(files[i], level = 1, varname = "water_u", band = 1))
      rr_v <- rotate(raster(files[i], level = 1, varname = "water_v", band = 1))
      
    } else {
      rr_u <- raster(files[i], level = 1, varname = "water_u", band = 1)
      
      rr_v <- raster(files[i], level = 1, varname = "water_v", band = 1)
      
    }
    
    # We already have currents u (rr_u) and v (rr_v). Optionally, we
    # can save these in a stack with:
    # stack_u <- stack(stack_u, rr_u)
    # stack_v <- stack(stack_v, rr_v)
    
    # Use rr_u and rr_v to calculate direction and speed
    
    
    sea_dir <- atan2(rr_u, rr_v) #from u and v components, calculatess direction
    rad2deg <- function(rad) {(rad * 180) / (pi)}
    sea_dir <- rad2deg(sea_dir)
    sea_dir[sea_dir < 0] <- 360 + sea_dir[sea_dir < 0]
    names(sea_dir) <- "direction"
    sea_spe <- sqrt( (rr_u * rr_u) + (rr_v * rr_v)) #and speed
    names(sea_spe) <- "speed"
    
    # Let´s make a stack with direction and speed, that is what rWind needs to calculate conductance
    day_stack <- stack(sea_dir, sea_spe)
    
    # We now use rWind to calculate the Conductance matrix
    conductance <- flow.dispersion(day_stack, type = "passive", output = "transitionLayer")
    
    # LetÂ´s save the conductance matrix of this day in a list. IÃÂ´ll use the
    #date as the name of the element within the list (each position in a list can get a name)
    #conduct_list[[i]] <- conductance
    #names(conduct_list)[i] <- files[i] %>% str_replace("_.nc", "")
    
    # And letÃÂ´s put the day_stack in the list
    #day_stacks_list[[i]] <- day_stack
    #names(day_stacks_list)[i] <- files[i] %>% str_replace("_.nc", "")
    
    cost_table <- costDistance(conductance, as.matrix(geo_points)) %>%    # this
      #is the function that does it
      as_tibble %>%     # instead of a matrix, I save it as a table
      set_names(rownames(geo_points))    # with corrent column names
    
    cost_list[[i]] <- cost_table   # move the result to the list
    names(cost_list)[i] <- files[i] %>% str_replace("_.nc", "")    # and give it a proper name
    
    
    nc_close(nc)
    
    
  }
  
  return(cost_list)
  
}






#-------------  Calculate conductance from HYCOM data and cost matrix




calculate_conductance <- function(files) {
  
  conduct_list <- list()
  cost_list <- list()
  
  for(i in 1:length(files)){
    
    
    
    cat("File ", i, "/n")  # to keep track of analyzed files
    
    #Read u and v from the netcdf raster (.nc extension)
    
    # nc_open from the netCDF package allows for getting the names of bands on 
    #other specific properties of netCDF that are lacking in the raster package
    
    nc <- nc_open(files[i])
    if(length(which(abs(nc$var$water_u$dim[[1]]$vals) > 180)) >= 1){  #depending on their coords, some rasters need to be rotated
      rr_u <- rotate(raster(files[i], level = 1, varname = "water_u", band = 1))
      rr_v <- rotate(raster(files[i], level = 1, varname = "water_v", band = 1))
      
    } else {
      rr_u <- raster(files[i], level = 1, varname = "water_u", band = 1)
      
      rr_v <- raster(files[i], level = 1, varname = "water_v", band = 1)
      
    }
    
    # We already have currents u (rr_u) and v (rr_v). Optionally, we
    # can save these in a stack with:
    # stack_u <- stack(stack_u, rr_u)
    # stack_v <- stack(stack_v, rr_v)
    
    # Use rr_u and rr_v to calculate direction and speed
    
    
    sea_dir <- atan2(rr_u, rr_v) #from u and v components, calculatess direction
    rad2deg <- function(rad) {(rad * 180) / (pi)}
    sea_dir <- rad2deg(sea_dir)
    sea_dir[sea_dir < 0] <- 360 + sea_dir[sea_dir < 0]
    names(sea_dir) <- "direction"
    sea_spe <- sqrt( (rr_u * rr_u) + (rr_v * rr_v)) #and speed
    names(sea_spe) <- "speed"
    
    # LetÃÂ´s make a stack with direction and speed, that is what rWind needs to calculate conductance
    day_stack <- stack(sea_dir, sea_spe)
    
    # We now use rWind to calculate the Conductance matrix
    conductance <- flow.dispersion(day_stack, type = "passive", output = "transitionLayer")
    
    # LetÃÂ´s save the conductance matrix of this day in a list. IÃÂ´ll use the date as the name of the element within the list (each position in a list can get a name)
    
    conduct_list[[i]] <- conductance
    names(conduct_list)[i] <- files[i] %>% str_replace("_.nc", "")
    
    # And letÃÂ´s put the day_stack in the list
    #day_stacks_list[[i]] <- day_stack
    #names(day_stacks_list)[i] <- files[i] %>% str_replace("_.nc", "")
    
    cost_table <- costDistance(conductance, as.matrix(Galapagos_points)) %>%    # this is the function that does it
      as_tibble %>%     # instead of a matrix, I save it as a table
      set_names(rownames(Galapagos_points))    # with corrent column names
    
    cost_list[[i]] <- cost_table   # move the result to the list
    names(cost_list)[i] <- files[i] %>% str_replace("_.nc", "")    # and give it a proper name
    
    
    nc_close(nc)
    
    
  }
  
  return(conduct_list)
  
}





### Same function but return day stack list


obtain_daystacklist <- function(files) {
  
  day_stacks_list<- list()
  conduct_list <- list()
  cost_list <- list()
  
  for(i in 1:length(files)){
    
    
    
    cat("File ", i, "/n")  # to keep track of analyzed files
    
    #Read u and v from the netcdf raster (.nc extension)
    
    # nc_open from the netCDF package allows for getting the names of bands on 
    #other specific properties of netCDF that are lacking in the raster package
    
    nc <- nc_open(files[i])
    if(length(which(abs(nc$var$water_u$dim[[1]]$vals) > 180)) >= 1){  #depending on their coords, some rasters need to be rotated
      rr_u <- rotate(raster(files[i], level = 1, varname = "water_u", band = 1))
      rr_v <- rotate(raster(files[i], level = 1, varname = "water_v", band = 1))
      
    } else {
      rr_u <- raster(files[i], level = 1, varname = "water_u", band = 1)
      
      rr_v <- raster(files[i], level = 1, varname = "water_v", band = 1)
      
    }
    
    # We already have currents u (rr_u) and v (rr_v). Optionally, we
    # can save these in a stack with:
    # stack_u <- stack(stack_u, rr_u)
    # stack_v <- stack(stack_v, rr_v)
    
    # Use rr_u and rr_v to calculate direction and speed
    
    
    sea_dir <- atan2(rr_u, rr_v) #from u and v components, calculatess direction
    rad2deg <- function(rad) {(rad * 180) / (pi)}
    sea_dir <- rad2deg(sea_dir)
    sea_dir[sea_dir < 0] <- 360 + sea_dir[sea_dir < 0]
    names(sea_dir) <- "direction"
    sea_spe <- sqrt( (rr_u * rr_u) + (rr_v * rr_v)) #and speed
    names(sea_spe) <- "speed"
    
    # LetÃÂ´s make a stack with direction and speed, that is what rWind needs to calculate conductance
    day_stack <- stack(sea_dir, sea_spe)
    
    # We now use rWind to calculate the Conductance matrix
    conductance <- flow.dispersion(day_stack, type = "passive", output = "transitionLayer")
    
    # LetÃÂ´s save the conductance matrix of this day in a list. IÃÂ´ll use the date as the name of the element within the list (each position in a list can get a name)
    
    #conduct_list[[i]] <- conductance
    #names(conduct_list)[i] <- files[i] %>% str_replace("_.nc", "")
    
    # And letÃÂ´s put the day_stack in the list
    day_stacks_list[[i]] <- day_stack
    names(day_stacks_list)[i] <- files[i] %>% str_replace("_.nc", "")
    
    cost_table <- costDistance(conductance, as.matrix(Galapagos_points)) %>%    # this is the function that does it
      as_tibble %>%     # instead of a matrix, I save it as a table
      set_names(rownames(Galapagos_points))    # with corrent column names
    
    cost_list[[i]] <- cost_table   # move the result to the list
    names(cost_list)[i] <- files[i] %>% str_replace("_.nc", "")    # and give it a proper name
    
    
    nc_close(nc)
    
    
  }
  
  return(day_stacks_list)
  
}



#------------- Calculate p-value matrix from shapiro test matrix


Calculate_Pvalue_matrix <- function(shapiro_matrix) {
  
  n_value <- rep(0, 1) # value that will include the objects '[i]' from each Shapiro Test
  
  pvalue_matrix <- matrix(nrow=46, ncol=46)
  
  
  for (i in 1:46) {
    for (j in 1:46) {
      n_value[i] <- shapiro_matrix[i,j]
      pvalue_matrix[i,j] <- n_value[[i]]$p.value
      
    }
  }
  
  return(pvalue_matrix)
  
}



#------------- Calculate percentage of data distributed normally from matrix normality


Calculate_percent_norm_from_normMatrix <- function(normality_matrix) {
  Table_Norm_test <- table(normality_matrix) / length(normality_matrix)
  percentage_normallity <- percent(Table_Norm_test[1], accuracy = 0.01)
  return(percentage_normallity)
}




Calculate_percent_norm_from_normMatrix <- function(normality_matrix) {
  Table_Norm_test <- table(normality_matrix_Canaries) / length(normality_matrix_Canaries)
  percentage_normallity <- percent(Table_Norm_test[1], accuracy = 0.01)
  return(percentage_normallity)
}



#-------------  Calculate matrix with the median of the 5% minimum values


Calculate_median_minimum5percent <- function(array, n_points) {
  
  median_minimums_matrix <- matrix(nrow=n_points, ncol=n_points) # matrix that will include the median of the minimum values
  
  for (i in 1:n_points) {
    for (j in 1:n_points) {
      sorted_values <- sort(array[i,j,]) # the values from each column and row are sorted in ascending order
      x <- round(length(sorted_values)*0.05,0) # 'x' indicates how many values are the 5% of the total number of values (excluding NAs values)
      median_minimums_matrix[i,j] <- median(sorted_values[1:x]) # taken 5% of the values and done the median of these values
    } }
  
  return(median_minimums_matrix)
}




#------------- Calculate correlation min_matrix and median_min matrix


Correlation_minMatrix_medianMinMatrix <- function(median_minimum_matrix, min_matrix) {
  
  r <- cor(c(min_matrix), c(median_minimums_matrix))
  R2 <- r^2
  
  return(R2)
}


#------------- Make matrix median min symmetric


Make_matrix_medianMin_symmetric <- function(median_min_mat) {
  idx<-median_min_mat<=t(median_min_mat) # matrix TRUE/FALSE 
  #if it's for the minimum
  
  median_min_mat.2<-median_min_mat*idx+t(median_min_mat)*(1L-idx)
  
  return(median_min_mat.2)
  
}


Make_minMatrix_symmetric <- function(min_matrix) {
  idx1<-min_matrix<=t(min_matrix) # matrix TRUE/FALSE if it's for the minimum
  min_matrix.2<-min_matrix*idx1+t(min_matrix)*(1L-idx1)
  
  return(min_matrix.2)
}




#------------- Transpose and smooth matrix, calculate Bray curtis distances and save matrix



Calculate_distmatrix_and_save <- function(matrix, file_path) {
  matrix <- t(matrix) # Transpose rows for columns
  rownames(matrix) <- c() # Delete names of the rows
  matrix[1, 1] <- NA # Delete the name of the first row
  
  colnames(matrix) <- as.character(unlist(matrix[1,])) # Put first row as header
  matrix = matrix[-1, ] # delete the first row
  
  rownames(matrix) <- matrix[,1] # Put first column as 'row header'
  matrix <- matrix[order(as.numeric(row.names(matrix)), decreasing=FALSE),] # Put first column as 'row header'
  matrix <- matrix[,-1]  # Delete the first column
  
  matrix_df <- as.data.frame(matrix) # Transform to 'data.frame'
  
  matrix_df[] <- lapply(matrix_df, function(x) as.numeric(as.character(x))) # all columns as numeric
  
  floristic_distance_bray <- vegdist(x = matrix_df, method="bray", binary=TRUE, diag=FALSE, upper=FALSE, na.rm = FALSE) 
  
  write.table(as.matrix(floristic_distance_bray), file = file_path,
              append = FALSE, quote = TRUE, sep = " ",
              eol = "\n", na = "NA", dec = ".", row.names = TRUE,
              col.names = TRUE, qmethod = c("escape", "double"),
              fileEncoding = "")
  
}




#------------- Get weighted matrix of minimum connecting points from median_min_mat_island


compute_mat_min_weighted <- function(matrix_conn) {
  
  # Convert the dataframe to an igraph object
  g <- graph_from_adjacency_matrix(matrix_conn, mode = "directed", weighted = TRUE)
  
  
  # Extract the weighted edges
  edges_df <- as_data_frame(g, what = "edges")
  
  # Eliminate "_" in names
  
  edges_df$from <- sub("_.*", "", edges_df$from)
  edges_df$to <- sub("_.*", "", edges_df$to)
  
  
  # Calculate minimum weight for each unique 'from' to 'to' combination
  result <- aggregate(weight ~ from + to, data = edges_df, FUN = min)
  
  result <- edges_df %>%
    group_by(from, to) %>%
    summarize(min_weight = min(weight))
  
  # Create a graph from the data frame
  g_single_point <- graph_from_data_frame(result, directed = TRUE)
  
  # Get the adjacency matrix with weights
  adj_matrix <- as_adjacency_matrix(g_single_point, attr = "min_weight", sparse = FALSE)
  
  return(adj_matrix)
  
}



# Plot weighted directed network

plot_directed_network <- function(directed_matrix, centroids, node_size, range_edge_width) {
  
  # Convert the matrix to an igraph object
  g <- graph_from_adjacency_matrix(directed_matrix, mode = "directed", weighted = TRUE)
  
  # Handle zero weights to avoid division by zero
  #E(g)$weight[E(g)$weight == 0] <- NA
  
  # Create a data frame from the igraph object
  node_list <- data.frame(name = V(g)$name)
  
  # Merge the coordinates with the node list
  node_list <- merge(node_list, centroids, by.x = "name", by.y = "island", all.x = TRUE)
  
  # Ensure no missing coordinates
  node_list <- na.omit(node_list)
  
  # Ensure the graph vertices match the coordinates names
  V(g)$name <- as.character(V(g)$name)
  node_list$name <- as.character(node_list$name)
  
  # Convert the igraph object to a tidygraph object
  tg <- as_tbl_graph(g)
  
  # Add coordinates to the tidygraph object
  tg <- tg %>%
    activate(nodes) %>%
    left_join(node_list, by = c("name" = "name"))
  
  # Invert the weights: higher values become lower and vice versa
  E(tg)$inv_weight <- 1 / E(tg)$weight
  
  # Normalize the inverted weights for better visualization
  max_inv_weight <- max(E(tg)$inv_weight, na.rm = TRUE)
  E(tg)$norm_inv_weight <- E(tg)$inv_weight / max_inv_weight
  
  
  # Calculate in-degree and out-degree based on the inverted weights
  V(tg)$in_degree <- strength(tg, mode = "in", weights = E(tg)$inv_weight)
  V(tg)$out_degree <- strength(tg, mode = "out", weights = E(tg)$inv_weight)
  
  
  # Plot the graph using ggraph with spatial coordinates and color nodes based on in-degree
  netplot_indegree <- ggraph(tg, layout = 'manual', x = V(tg)$long, y = V(tg)$lat) + 
    geom_edge_link(aes(width = norm_inv_weight, edge_alpha = norm_inv_weight), arrow = arrow(length = unit(3, 'mm')), end_cap = circle(3, 'mm'), show.legend = FALSE) +
    geom_node_point(aes(color = in_degree), size = node_size) +
    scale_edge_width(range = range_edge_width, name = "Normalized Connectivity") + # Adjust the range for better visualization
    scale_color_viridis(name = "In-Degree", option = "turbo") + # Use viridis color scale
    geom_node_text(aes(label = name), repel = TRUE, 
                   nudge_y = 0.1, 
                   nudge_x = -0.04, size = 5, segment.size = 0.15) +
    theme_classic() + 
    labs(x = "long",
         y = "lat",
         edge_width = "Normalized connectivity",
         edge_alpha = NULL)+
    guides(color = guide_colorbar(barwidth = 8, barheight = 1,
                                  title.position = "top", label.position = "bottom",
                                  title.hjust = 0.5,
                                  ticks = FALSE,
                                  breaks = quantile(V(tg)$in_degree, probs = c(0, 0.5, 1)))) + # Customize legend breaks and labels
    theme(legend.position = "bottom")+
    my_theme
  
  return(netplot_indegree)
  
}



# Plot weighted undirected network



plot_undirected_network <- function(directed_matrix, centroids, node_size) {
  
  # Assuming your adjacency matrix is named 'single_point_min_Galap'
  # Convert the matrix to an igraph object
  g <- graph_from_adjacency_matrix(directed_matrix, mode = "undirected", weighted = TRUE)
  
  # Handle zero weights to avoid division by zero
  #E(g)$weight[E(g)$weight == 0] <- NA
  
  # Create a data frame from the igraph object
  node_list <- data.frame(name = V(g)$name)
  
  # Merge the coordinates with the node list
  node_list <- merge(node_list, centroids, by.x = "name", by.y = "island", all.x = TRUE)
  
  # Ensure no missing coordinates
  node_list <- na.omit(node_list)
  
  # Ensure the graph vertices match the coordinates names
  V(g)$name <- as.character(V(g)$name)
  node_list$name <- as.character(node_list$name)
  
  # Convert the igraph object to a tidygraph object
  tg <- as_tbl_graph(g)
  
  # Add coordinates to the tidygraph object
  tg <- tg %>%
    activate(nodes) %>%
    left_join(node_list, by = c("name" = "name"))
  
  # Invert the weights: higher values become lower and vice versa
  E(tg)$inv_weight <- 1 / E(tg)$weight
  
  # Normalize the inverted weights for better visualization
  max_inv_weight <- max(E(tg)$inv_weight, na.rm = TRUE)
  E(tg)$norm_inv_weight <- E(tg)$inv_weight / max_inv_weight
  
  # Calculate total degree
  V(tg)$degree <- strength(tg, mode = "all",
                           weights = E(tg)$inv_weight)
  
  # Plot the graph using ggraph with spatial coordinates and color nodes based on total degree
  netplot_degree <- ggraph(tg, layout = 'manual', x = V(tg)$long, y = V(tg)$lat) + 
    geom_edge_link(aes(width = norm_inv_weight), show.legend = FALSE) +
    geom_node_point(aes(color = degree), size = node_size) +
    scale_edge_width(range = range(sqrt(E(tg)$norm_inv_weight)), name = "Normalized Connectivity") + # Adjust the range for better visualization
    scale_color_viridis(name = "Total Degree", option = "turbo") + # Use viridis color scale
    geom_node_text(aes(label = name), repel = TRUE, 
                   nudge_y = 0.1, 
                   nudge_x = -0.04, size = 5, segment.size = 0.15) +
    theme_classic() + 
    labs(x = "Longitude",
         y = "Latitude",
         edge_width = "Normalized connectivity") +
    guides(color = guide_colorbar(barwidth = 8, barheight = 1,
                                  title.position = "top", label.position = "bottom",
                                  title.hjust = 0.5,
                                  ticks = FALSE,
                                  breaks = quantile(V(tg)$degree, probs = c(0, 0.5, 1)))) + # Customize legend breaks and labels
    theme(legend.position = "bottom") +
    my_theme
  
  return(netplot_degree)
  
}











################################################


#------------- Load floristic matrix, order rows and columns, and convert to distance matrix

#### table

load_order_convert_table_to_distanceMat <- function(filepath) {
  flor_matrix<-read.table(filepath)
  
  
  # Order dataframe's rows and columns
  flor_matrix1 <- flor_matrix[order(rownames(flor_matrix)), order(colnames(flor_matrix))]
  
  
  # convert to distance matrix
  distFlor<-as.dist(flor_matrix1)
  
  return(distFlor)
  
}

#### csv

load_order_convert_csv_to_distanceMat <- function(filepath) {
  
  # Load
  curr_matrix<-read.csv(filepath)
  
  # rownames
  rownames(curr_matrix) <- curr_matrix[,1]
  curr_matrix1<-curr_matrix[,-1]
  
  # Order 2
  
  
  curr_matrix2 <- curr_matrix1[order(rownames(curr_matrix1)), order(colnames(curr_matrix1))]
  
  # Convert to distance matrix
  distCurr<-as.dist(curr_matrix2)
  
  
  return(distCurr)
  
}



#------------- Run Procrustes test between distance matrices

Compute_Protest_distMat <- function(distMat_rot, distMat_base, k) {
  add1 <-  !(is.euclid(distMat_base))
  test1<-cmdscale(distMat_base, k = k, eig = TRUE, add = add1)
  ord1<-ordiplot(test1, type = "text", main = "PCoA Bray Distances")
  
  
  add2 <-  !(is.euclid(distMat_rot))
  test2<-cmdscale(distMat_rot, k = k, eig = TRUE, add = add2)
  ord2<-ordiplot(test2, type = "text", main = "PCoA Bray Distances")
  
  protest <- protest(ord1,ord2,  symmetric = TRUE, scores = "sites", permutations=100000)
  
  return(protest)
}




# --------------- Obtain a vector with the number of plants per island. 
#The input is a filtered floristic dataframe.

header.true <- function(df) {
  names(df) <- as.character(unlist(df[1,]))
  df[-1,]
}



Calculate_Nplants_perIsland <- function(df_flora) {
  df_header_names<-header.true(df_flora)
  
  island.count_sums<-colSums(df_header_names != 0)
  
  names(island.count_sums) <- NULL
  
  island.count_sums1<-island.count_sums[-1]
  
  return(island.count_sums1)
}



# --------------- Obtain a vector with the percentage of plants per island. 
#The input is a filtered floristic dataframe.



header.true <- function(df) {
  names(df) <- as.character(unlist(df[1,]))
  df[-1,]
}



fun_percentage <- function(i) {
  (i)/n_plants_total*100
}




Calculate_Percentage_plants_perIsland <- function(df_flora) {
  
  df_header_names<-header.true(df_flora)
  
  island.count_sums<-colSums(df_header_names != 0)
  
  island.count_sums1<-island.count_sums[-1] # eliminate first value which is the total counts
  
  island.count_sums1[order(names(island.count_sums1))] #order islands alphabetically
  
  names(island.count_sums1) <- NULL
  
  percentage_plants<-sapply(island.count_sums1, fun_percentage)
  
  percentage_plants_perIsland<-round(percentage_plants, digits=2)
  
  return(percentage_plants_perIsland)
  
}




# --------------- Calculate nodes (islands) degree

calculate_degrees_perIsland <- function(dist_matrix) {
  
  sim_matrix <- 1/dist_matrix # inverse values to obtain similarity matrix
  
  links_conn_curr<-dist2list(sim_matrix) # obtain edge list
  
  net_curr <- graph_from_data_frame(d=links_conn_curr, directed=F) #create graph from edgelist
  
  net_curr1 <- set_edge_attr(net_curr, "weight", value= links_conn_curr$value) # assign weight attribute
  
  degrees<-strength(
    net_curr1,
    vids = V(net_curr),
    loops = FALSE,
    weights = net_curr$value
  )
  
  degrees[order(names(degrees))]
  
  names(degrees) <- NULL
  
  return(degrees)
}


# --------------- Create dataframe island name - degree

Create_dataframe_IslName_perc_N_degree <- function(Archipelago_name_vector,island_names_vector, Percentage_plants_vector1, Percentage_plants_vector2, Nplants_vector1, Nplants_vector2,Degree_vector1, Degree_vector2, area_vector, age_vector) {
  
  df_names_perc_N_deg<-data.frame(Archipelago_name_vector,island_names_vector, Percentage_plants_vector1, Percentage_plants_vector2, Nplants_vector1, Nplants_vector2,Degree_vector1, Degree_vector2, area_vector, age_vector)
  
  df_names_perc_N_deg1 = rename(df_names_perc_N_deg,
                                Archipelago = Archipelago_name_vector,
                                island_names = island_names_vector,
                                Percentage_plants_tha = Percentage_plants_vector1,
                                Percentage_plants_litt = Percentage_plants_vector2,
                                Nplants_tha = Nplants_vector1,
                                Nplants_litt = Nplants_vector2,
                                Degree_curr = Degree_vector1,
                                Degree_geo = Degree_vector2,
                                area = area_vector,
                                age = age_vector)
  
  return(df_names_perc_N_deg1)
  
}


