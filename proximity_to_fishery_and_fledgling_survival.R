
# Penguin tracking, environmental data, 
# vessel proximity, survival analysis, and behavioral states
# 

# Clean memory and load libraries (each loaded once)
gc()

library(dplyr)
library(lubridate)
library(sf)
library(ggplot2)
library(patchwork)
library(terra) # using terra where possible, but raster still used for brick extraction
library(raster)
library(ncdf4)
library(maps)
library(mapdata)
library(matrixStats)
library(purrr)
library(survival)
library(boot)
library(ggsurvfit)
library(aniMotum)
library(RVAideMemoire)
library(EMbC)

library(glmmTMB)
library(pgirmess)
library(nlme)
library(car)
library(gganimate)

library(tidyr)




# 
# -----1. Load and process tracking data-----
# 

df <- read.csv("PenguinData/Kruger_etal_data.csv")
df$TimeStamp <- as.POSIXct(df$Loc..date, format = "%d-%m-%Y %H:%M:%S", tz = "GMT")
df$DayStamp <- as.POSIXct(substring(df$Loc..date, 1, 10), format = "%d-%m-%Y", tz = "GMT")

# Hinke et al. 2020 data
jh <- read.csv("PenguinData/Hinke_etal_data.csv")
pyn <- jh
pyn$TimeStamp <- as.POSIXct(pyn$LOC_DATE, format = "%m/%d/%Y %H:%M:%S", tz = "GMT")
pyn$DayStamp <- as.POSIXct(substring(pyn$LOC_DATE, 1, 10), format = "%m/%d/%Y", tz = "GMT")
pyn$year <- year(pyn$TimeStamp)
pyn <- filter(pyn, year != "2018")

# Prepare IDs
ids2017 <- pyn %>%
  group_by(SPP, PTT) %>%
  summarise(N = length(LATITUDE), .groups = "drop") %>%
  dplyr::select(species = SPP, id = PTT, N)

ids2025 <- df %>%
  group_by(Platform.ID.No.) %>%
  summarise(N = length(Longitude), .groups = "drop") %>%
  dplyr::select(id = Platform.ID.No., N) %>%
  mutate(species = "PYN")

all_ids <- rbind(ids2017, ids2025)

# Combine raw tracking data
jhdf <- data.frame(id = pyn$PTT, DayStamp = pyn$DayStamp, date = pyn$TimeStamp,
                   year = year(pyn$DayStamp), lon = pyn$LONGITUDE, lat = pyn$LATITUDE, lc = 1)
tsdf <- data.frame(id = df$Platform.ID.No., DayStamp = df$DayStamp, date = df$TimeStamp,
                   year = year(df$DayStamp), lon = df$Longitude, lat = df$Latitude, lc = df$Loc..quality)
adf <- rbind(jhdf, tsdf)
adf <- filter(adf, id != "165242")   # remove problem ID
table(adf$lc)
# Prepare tracks for aniMotum
tracks <- data.frame(id = adf$id, date = adf$date, lc = as.factor(adf$lc), lon = adf$lon, lat = adf$lat)
#tracks <- tracks[!duplicated(tracks[, c("id", "date")]), ]
tracks <- subset(tracks, id != "255211")   # only one day of data

table(year(tracks$date))

# Fit SSM
fit.crw <- fit_ssm(tracks, model = "crw", time.step = 1)

res.crw <- osar(fit.crw) # be patient, this take a while


res.crw$year <- year(res.crw$date)


# Residual plots 

# for acf 
conf.level <- 0.95
res.x<-res.crw$residual[res.crw$coord=="x"]
res.y<-res.crw$residual[res.crw$coord=="y"]

ciline.x <- qnorm((1 - conf.level)/2)/sqrt(length(res.x))
bacf.x <- acf(res.x, plot = FALSE)
bacfdf.x <- with(bacf.x, data.frame(lag, acf))

ciline.y <- qnorm((1 - conf.level)/2)/sqrt(length(res.y))
bacf.y <- acf(res.y, plot = FALSE)
bacfdf.y <- with(bacf.y, data.frame(lag, acf))

head(res.crw)
res.crw$year<-year(res.crw$date)
(ggplot(subset(res.crw,coord=="x" & year==2017),aes(date,residual))+
    geom_hline(yintercept = 0,linetype="dashed")+
    geom_smooth()+
    geom_point(alpha=0.2)+
    theme_bw()+ggtitle(label="a.")+
    
    ggplot(subset(res.crw,coord=="x" & year==2025),aes(date,residual))+
    geom_hline(yintercept = 0,linetype="dashed")+
    geom_smooth()+
    geom_point(alpha=0.2)+
    theme_bw()+ggtitle(label="b."))/
  
  (ggplot(subset(res.crw,coord=="y" & year==2017),aes(date,residual))+
     geom_hline(yintercept = 0,linetype="dashed")+
     geom_smooth()+
     geom_point(alpha=0.2)+
     theme_bw()+ggtitle(label="c.")+
     
     ggplot(subset(res.crw,coord=="y" & year==2025),aes(date,residual))+
     geom_hline(yintercept = 0,linetype="dashed")+
     geom_smooth()+
     geom_point(alpha=0.2)+
     theme_bw()+ggtitle(label="d."))/
  
  
  (ggplot(subset(res.crw,coord=="x"),aes(sample=residual))+
     geom_qq(alpha=0.2)+geom_qq_line(linetype="dashed")+
     theme_bw()+ggtitle(label="e.")+
     
     
     ggplot(subset(res.crw,coord=="y"),aes(sample=residual))+
     geom_qq(alpha=0.2)+geom_qq_line(linetype="dashed")+
     
     theme_bw()+ggtitle(label="f."))/
  
  (ggplot(data=bacfdf.x, mapping=aes(x=lag, y=acf)) +
     geom_bar(stat = "identity", position = "identity")+theme_bw()+
     ggtitle(label="g.")+
     
     ggplot(data=bacfdf.y, mapping=aes(x=lag, y=acf)) +
     geom_bar(stat = "identity", position = "identity")+theme_bw()+ggtitle(label="h."))

#ggsave("CRW_diagnostic_plots.png", width = 8, height = 8, dpi = 600)


# Extract  predicted tracks

track.predicted <- aniMotum::grab(fit.crw, what = "predicted", group = TRUE)
track.predicted$DayStamp <- as.POSIXct(paste(year(track.predicted$date),
                                             month(track.predicted$date),
                                             day(track.predicted$date),sep="-"), format = "%Y-%m-%d")

track.2017 <- filter(track.predicted, year(DayStamp) == 2017)
track.2025 <- filter(track.predicted, year(DayStamp) == 2025)

summary(track.predicted)


# 
# -----2. Environmental variables – Net Primary Productivity (npp)-----
# 
# function to extract values, handling missing dates
extract_env <- function(track_df, brick_obj, var_name) {
  extracted <- vector("list", nrow(track_df))
  for (i in seq_len(nrow(track_df))) {
    current_timestamp <- as.character(track_df$DayStamp[i])
    current_lat <- track_df$lat[i]
    current_lon <- track_df$lon[i]
    timestamp_index <- which(as.character(brick_obj@z$`Date/time`) == current_timestamp)
    if (length(timestamp_index) == 0) {
      extracted[[i]] <- NA_real_
      next
    }
    current_raster <- brick_obj[[timestamp_index]]
    sp <- SpatialPoints(matrix(c(current_lon, current_lat), ncol = 2))
    extracted[[i]] <- raster::extract(current_raster, sp, method = "bilinear", fun = "mean", na.rm = TRUE)
  }
  return(unlist(extracted))
}

# 2025 npp
ncpath <- "C:/AWRF/PenguinFledglings/Copernicus/"
ncname1 <- "cmems_mod_glo_bgc_my_0.25deg_P1D-m_1779109338425"
ncfname1 <- paste0(ncpath, ncname1, ".nc")
tmp_brick1 <- brick(ncfname1, varname = "nppv")
track.2025$npp <- extract_env(track.2025, tmp_brick1, "nppv")

summary(track.2025$npp)

# 2017 npp
ncname2 <- "cmems_mod_glo_bgc_my_0.25deg_P1D-m_1779109435194"
ncfname2 <- paste0(ncpath, ncname2, ".nc")
tmp_brick2 <- brick(ncfname2, varname = "nppv")
track.2017$npp <- extract_env(track.2017, tmp_brick2, "nppv")

summary(track.2017$npp)


# 
# ----- 3. Sea water velocity – uo (eastward) and vo (northward)------
# 

# 2025 uo
ncname_uo25 <- "cmems_mod_glo_phy_my_0.083deg_P1D-m_1779109823356"
ncfname_uo25 <- paste0("C:/AWRF/SWV/", ncname_uo25, ".nc")
brick_uo25 <- brick(ncfname_uo25, varname = "uo")
track.2025$cuo <- extract_env(track.2025, brick_uo25, "uo")

summary(track.2025$cuo)

# 2017 uo
ncname_uo17 <- "cmems_mod_glo_phy_my_0.083deg_P1D-m_1779109673336"
ncfname_uo17 <- paste0("C:/AWRF/SWV/", ncname_uo17, ".nc")
brick_uo17 <- brick(ncfname_uo17, varname = "uo")
track.2017$cuo <- extract_env(track.2017, brick_uo17, "uo")

summary(track.2017$cuo)

# 2025 vo
brick_vo25 <- brick(ncfname_uo25, varname = "vo")   # same file, different variable
track.2025$cvo <- extract_env(track.2025, brick_vo25, "vo")
summary(track.2025$cvo)


# 2017 vo
brick_vo17 <- brick(ncfname_uo17, varname = "vo")
track.2017$cvo <- extract_env(track.2017, brick_vo17, "vo")
summary(track.2017$cvo)
# 
# -----4. Merge tracks and calculate water resistance-----
# 

tracks <- rbind(track.2017, track.2025)

# Bird speed and direction (from SSM output, u and v are in tracks)
tracks$bird.dir <- atan2(tracks$u, tracks$v)
tracks$current.dir <- atan2(tracks$cuo, tracks$cvo)
tracks$bird.speed <- sqrt(tracks$u^2 + tracks$v^2)
tracks$current.speed <- sqrt(tracks$cuo^2 + tracks$cvo^2)
tracks$relative.angle <- tracks$bird.dir - tracks$current.dir
tracks$relative.angle2 <- (tracks$bird.dir - tracks$current.dir + pi) %% (2 * pi) - pi
tracks$current_effect <- tracks$current.speed * cos(tracks$relative.angle2)
tracks$year <- year(tracks$date)

all_ids$id<-as.character(all_ids$id)

# Combine with species info
all_tracks <- tracks %>%
  dplyr::select(id, lon, lat, npp, ce = current_effect, cs = current.speed, date, year) %>%
  left_join(all_ids, by = "id")





# 
# -----5. Proximity to vessels -----
# 

# process fishery data 


gc()

# Function to extract ship name from zip filename
get_ship_name_from_zip <- function(zip_path) {
  # Extract just the filename without path
  zip_name <- basename(zip_path)
  # Remove the date range part and .zip extension
  # Pattern: "SHIP_NAME (YYYY-MM-DD - YYYY-MM-DD).zip"
  ship_name <- gsub("\\s*\\([0-9]{4}-[0-9]{2}-[0-9]{2}\\s*-\\s*[0-9]{4}-[0-9]{2}-[0-9]{2}\\)\\.zip$", "", zip_name)
  return(ship_name)
}

# Function to process a single zip file and return its data
process_ship_zip <- function(zip_path, years_to_keep = 2012:2025) {
  
  ship_name <- get_ship_name_from_zip(zip_path)
  cat("Processing:", ship_name, "\n")
  
  # Get list of CSV files inside the zip
  csv_files <- unzip(zip_path, list = TRUE)
  csv_files <- csv_files[grepl("\\.csv$", csv_files$Name, ignore.case = TRUE), ]
  
  if(nrow(csv_files) == 0) {
    cat("  No CSV files found in", ship_name, "\n")
    return(NULL)
  }
  
  # Read all CSVs from the zip
  ship_data <- list()
  
  for(i in 1:nrow(csv_files)) {
    # Read specific file from zip
    temp_data <- tryCatch({
      read.csv(unz(zip_path, csv_files$Name[i]), stringsAsFactors = FALSE)
    }, error = function(e) {
      cat("    Error reading", csv_files$Name[i], ":", e$message, "\n")
      return(NULL)
    })
    
    if(!is.null(temp_data) && nrow(temp_data) > 0) {
      # Check if required columns exist (adjust column names as needed)
      # Common GFW column names: timestamp, lon, lat, or longitude, latitude
      col_names <- tolower(names(temp_data))
      
      if(any(grepl("timestamp", col_names)) && 
         any(grepl("lon", col_names) | grepl("longitude", col_names)) &&
         any(grepl("lat", col_names) | grepl("latitude", col_names))) {
        
        # Standardize column names
        temp_data <- temp_data %>%
          rename_with(~"timestamp", .cols = matches("^timestamp$|^time$|^date$", ignore.case = TRUE)) %>%
          rename_with(~"lon", .cols = matches("^lon$|^longitude$|^long$", ignore.case = TRUE)) %>%
          rename_with(~"lat", .cols = matches("^lat$|^latitude$", ignore.case = TRUE)) %>%
          mutate(
            ship_id = ship_name,
            timestamp = as.POSIXct(timestamp, format = "%Y-%m-%d T %H:%M:%S", tz = "GMT"),
            year = year(timestamp),
            date_only = as.Date(timestamp)
          ) %>%
          filter(year %in% years_to_keep) %>%
         dplyr::select(ship_id, timestamp, year, date_only, lon, lat, everything())
        
        ship_data[[length(ship_data) + 1]] <- temp_data
      }
    }
  }
  
  if(length(ship_data) == 0) {
    cat("  No valid data for", ship_name, "\n")
    return(NULL)
  }
  
  result <- bind_rows(ship_data)
  cat("  Loaded", nrow(result), "records for", ship_name, "\n")
  
  return(result)
}

# Main function: Process all zip files year by year to avoid memory issues
process_vessels_by_year <- function(zip_folder_path, output_folder = "C:/AWRF/PenguinFledglings/GFWData/processed_yearly", 
                                    years = 2012:2025, chunk_size = 5) {
  
  # Create output directory
  if(!dir.exists(output_folder)) dir.create(output_folder, recursive = TRUE)
  
  # Get all zip files
  zip_files <- list.files(zip_folder_path, pattern = "\\.zip$", full.names = TRUE, recursive = FALSE)
  cat("Found", length(zip_files), "zip files to process\n")
  
  # Process each year separately
  for(year in years) {
    cat("\n========== Processing year", year, "==========\n")
    
    yearly_file <- file.path(output_folder, paste0("vessels_", year, ".rds"))
    
    # Skip if already processed
    if(file.exists(yearly_file)) {
      cat("Year", year, "already processed, skipping...\n")
      next
    }
    
    # Process ships in chunks to avoid memory issues
    all_yearly_data <- list()
    
    for(i in seq(1, length(zip_files), by = chunk_size)) {
      chunk_end <- min(i + chunk_size - 1, length(zip_files))
      cat("Processing ships", i, "to", chunk_end, "of", length(zip_files), "\n")
      
      chunk_files <- zip_files[i:chunk_end]
      
      chunk_data <- map_dfr(chunk_files, function(zip_path) {
        ship_data <- process_ship_zip(zip_path, years_to_keep = year)
        return(ship_data)
      })
      
      if(!is.null(chunk_data) && nrow(chunk_data) > 0) {
        all_yearly_data[[length(all_yearly_data) + 1]] <- chunk_data
      }
      
      # Clear memory
      gc()
    }
    
    if(length(all_yearly_data) > 0) {
      yearly_combined <- bind_rows(all_yearly_data)
      
      # Save as RDS (more efficient than CSV)
      saveRDS(yearly_combined, yearly_file)
      cat("Saved", nrow(yearly_combined), "records for", year, "\n")
      
      # Also save a spatial version for quick access
      yearly_sf <- st_as_sf(yearly_combined, coords = c("lon", "lat"), crs = 4326, remove = FALSE)
      saveRDS(yearly_sf, file.path(output_folder, paste0("vessels_", year, "_sf.rds")))
      
    } else {
      cat("No data for year", year, "\n")
    }
    
    # Clear memory
    rm(all_yearly_data, yearly_combined)
    gc()
  }
  
  cat("\nProcessing complete! Files saved to:", output_folder, "\n")
}

# Run the extraction (adjust path as needed)
process_vessels_by_year(
  zip_folder_path = "C:/AWRF/PenguinFledglings/GFWData/",
  output_folder = "C:/AWRF/PenguinFledglings/GFWData/processed_yearly/",
  years = c(2017,2025),
  chunk_size = 5  # Process 5 ships at a time
)



# Function to load all yearly vessel data
load_all_vessel_data <- function(processed_folder = "C:/AWRF/PenguinFledglings/GFWData/processed_yearly/") {
  
  # Find all RDS files
  vessel_files <- list.files(processed_folder, pattern = "^vessels_[0-9]{4}\\.rds$", full.names = TRUE)
  
  if(length(vessel_files) == 0) {
    cat("No vessel RDS files found in", processed_folder, "\n")
    return(NULL)
  }
  
  cat("Found", length(vessel_files), "yearly files\n")
  
  # Load and combine all years
  all_vessels <- map_dfr(vessel_files, function(file) {
    year <- gsub(".*vessels_([0-9]{4})\\.rds", "\\1", basename(file))
    cat("Loading", year, "...\n")
    data <- readRDS(file)
    return(data)
  })
  
  cat("Total records loaded:", nrow(all_vessels), "\n")
  cat("Date range:", min(all_vessels$timestamp), "to", max(all_vessels$timestamp), "\n")
  cat("Ships:", length(unique(all_vessels$ship_id)), "\n")
  
  return(all_vessels)
}

# Load all data (caution: may be large!)
all_vessels <- load_all_vessel_data("C:/AWRF/PenguinFledglings/GFWData/processed_yearly/")



all_vessels<-all_vessels[!duplicated(
  all_vessels[, c("timestamp","ship_id")]), ]

fishery<-all_vessels%>%
  na.omit()%>%
  group_by(ship_id)%>%
  filter(n() >= 10) %>%
  ungroup()


head(fishery)

# proximity function 

calculate_proximity_continuous <- function(predator_data, vessel_path, 
                                           threshold_km = 50, 
                                           time_window_hours = 1) {
  if (!inherits(predator_data$date, "POSIXct")) {
    predator_data$date <- as.POSIXct(predator_data$date, origin = "1970-01-01", tz = "GMT")
  }
  results <- predator_data
  results$encounter <- FALSE
  results$min_distance_km <- NA
  results$vessels_nearby <- 0
  results$closest_vessel_id <- NA
  
  
  for (yr in unique(predator_data$year)) {
    cat("\nProcessing year", yr, "...")
    vessel_file <- file.path(vessel_path, paste0("vessels_", yr, ".rds"))
    if (!file.exists(vessel_file)) {
      cat(" No vessel data\n")
      next
    }
    vessels <- readRDS(vessel_file)
    if (!inherits(vessels$timestamp, "POSIXct")) 
      vessels$timestamp <- as.POSIXct(vessels$timestamp, tz = "GMT")
    
    cat(" Loaded", nrow(vessels), "vessel positions\n")
    year_indices <- which(predator_data$year == yr)
    n_pred <- length(year_indices)
    cat("  Checking", n_pred, "predator positions...\n")
    encounters_found <- 0
    
    for (i in seq_len(n_pred)) {
      if (i %% 100 == 0) 
        cat("    Progress:", round(i/n_pred*100), "% (", encounters_found, "encounters)\n")
      idx <- year_indices[i]
      pred_time <- predator_data$date[idx]
      pred_lon <- predator_data$lon[idx]
      pred_lat <- predator_data$lat[idx]
      time_start <- pred_time - hours(time_window_hours)
      time_end <- pred_time + hours(time_window_hours)
      
      vessels_time <- vessels %>%
        filter(timestamp >= time_start, timestamp <= time_end)
      
      if (nrow(vessels_time) == 0) {
        results$min_distance_km[idx] <- 500   # large finite value as placeholder
        results$closest_vessel_distance_km[idx] <- 500
        next
      }
      
      # Haversine distance
      lon1 <- pred_lon * pi/180
      lat1 <- pred_lat * pi/180
      lon2 <- vessels_time$lon * pi/180
      lat2 <- vessels_time$lat * pi/180
      dlon <- lon2 - lon1
      dlat <- lat2 - lat1
      a <- sin(dlat/2)^2 + cos(lat1) * cos(lat2) * sin(dlon/2)^2
      c <- 2 * asin(sqrt(a))
      dist_km <- 6371 * c
      
      min_dist <- min(dist_km)
      closest <- vessels_time$ship_id[which.min(dist_km)]
      vessels_close <- sum(dist_km <= threshold_km)
      
      results$min_distance_km[idx] <- min_dist
      results$closest_vessel_distance_km[idx] <- min_dist
      if (vessels_close > 0) {
        encounters_found <- encounters_found + 1
        results$encounter[idx] <- TRUE
        results$vessels_nearby[idx] <- vessels_close
        results$closest_vessel_id[idx] <- closest
      }
    }
    cat("  ✅ Year", yr, "complete:", encounters_found, "encounters (", 
        round(encounters_found/n_pred*100, 2), "%)\n")
  }
  return(results)
}

# Run proximity analysis
vessel_path <- "C:/AWRF/PenguinFledglings/GFWData/processed_yearly/"
all_proximity <- calculate_proximity_continuous(all_tracks, vessel_path, 50, 1)

# Save full results
saveRDS(all_proximity, "chinstrap_proximity_final.rds")

# 
# -----6. Summaries (species, yearly, individual)-----
# 

species_summary <- all_proximity %>%
  group_by(species) %>%
  summarise(
    positions = n(),
    encounters = sum(encounter),
    encounter_rate = round(encounters / positions * 100, 4),
    mean_distance_km = round(mean(min_distance_km, na.rm = TRUE), 2),
    median_distance_km = round(median(min_distance_km, na.rm = TRUE), 2),
    .groups = "drop"
  ) %>%
  arrange(desc(encounter_rate))

print(species_summary)

id_summary <- all_proximity %>%
  group_by(id, species) %>%
  summarise(
    positions = n(),
    encounters = sum(encounter),
    encounter_rate = round(encounters / positions * 100, 4),
    mean_distance_km = round(mean(min_distance_km, na.rm = TRUE), 2),
    .groups = "drop"
  ) %>%
  arrange(desc(encounter_rate))


summary(id_summary)

yearly_summary <- all_proximity %>%
  group_by(year) %>%
  summarise(
    positions = n(),
    encounters = sum(encounter),
    encounter_rate = round(encounters / positions * 100, 2),
    .groups = "drop"
  )

# Save summaries
output_path <- "C:/AWRF/PenguinFledglings/Proximity_Results_Final/"
if(!dir.exists(output_path)) dir.create(output_path, recursive = TRUE)
write.csv(species_summary, file.path(output_path, "species_summary_final.csv"), row.names = FALSE)
write.csv(yearly_summary, file.path(output_path, "yearly_summary_final.csv"), row.names = FALSE)
write.csv(id_summary, file.path(output_path, "id_summary_final.csv"), row.names = FALSE)



# year comparison


# Compute individual-level metrics # this serves for the plot and the next analysis
last_positions <- all_proximity %>%
  group_by(id, year, species) %>%
  arrange(date) %>%
  summarise(
    deployment_date = min(date),
    deployment_doy = yday(min(date)),
    last_date = max(date),
    tracking_duration_days = as.numeric(difftime(max(date), min(date), units = "days")),
    first_encounter_date = min(date[encounter == TRUE], na.rm = TRUE),
    first_encounter_doy = yday(first_encounter_date),
    days_to_first_encounter = as.numeric(difftime(first_encounter_date, deployment_date, units = "days")),
    encountered_fisheries = any(encounter == TRUE, na.rm = TRUE),
    n_locations = n(),
    n_encounters = sum(encounter == TRUE, na.rm = TRUE),
    encounter_rate = n_encounters / n_locations * 100,
    min_dist = min(min_distance_km, na.rm = TRUE),
    npp = mean(npp, na.rm = TRUE),
    ce = mean(ce, na.rm = TRUE),
    cs = mean(cs, na.rm = TRUE),
    .groups = "drop"
  )




last_positions<-last_positions %>%
  group_by(year)%>%
  mutate(time_proportion=tracking_duration_days/max(tracking_duration_days))

ggplot(last_positions, aes(time_proportion)) +
  geom_histogram(aes(y = after_stat(density)), 
                 bins = 30, 
                 fill = "steelblue", 
                 colour = "black", 
                 alpha = 0.7) +
  geom_density(colour = "red3", linewidth = 1.2) +
  geom_vline(xintercept = 0.5, linetype = "dashed", colour = "red3", linewidth = 1) +
  theme_bw() +
  labs(
    x = "Proportion of maximum tracking duration",
    y = "Density",
    title = "Distribution of tracking durations",
    subtitle = "Dashed line = 50% threshold"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, colour = "grey40")
  )
  

# this justifies using 0.5 as the transmission threshold

#ggsave("tracking_duration_histogram.png", width = 5, height = 5, dpi = 600)




summary(last_positions$tracking_duration_days)
summary(last_positions$tracking_duration_days[last_positions$year=="2017"])
summary(last_positions$tracking_duration_days[last_positions$year=="2025"])


perm.anova(tracking_duration_days~year,  data=last_positions, 
           nperm = 999,
           progress = TRUE)


summary(last_positions$npp)


perm.anova(npp~year,  data=last_positions, 
           nperm = 999,
           progress = TRUE)

summary(last_positions$ce)

perm.anova(ce~year,  data=last_positions, 
           nperm = 999,
           progress = TRUE)

summary(last_positions$encounter_rate)

length(all_proximity$lat[all_proximity$encounter==TRUE])/
  length(all_proximity$lat)


length(all_proximity$lat[all_proximity$encounter==TRUE & all_proximity$year=="2017"])/
  length(all_proximity$lat[all_proximity$year=="2017"])

length(all_proximity$lat[all_proximity$encounter==TRUE & all_proximity$year=="2025"])/
  length(all_proximity$lat[all_proximity$year=="2025"])


rates_over_zero<-filter(last_positions,encounter_rate>0)



chp<-filter(last_positions,species=="PYN")



perm.anova(encounter_rate~year,  data=last_positions, 
           nperm = 1999,
           progress = TRUE)

perm.anova(encounter_rate~year,  data=rates_over_zero, 
           nperm = 1999,
           progress = TRUE)

perm.anova(encounter_rate~year,  data=chp, 
           nperm = 999,
           progress = TRUE)



summary(rates_over_zero$days_to_first_encounter)

summary(rates_over_zero$days_to_first_encounter[rates_over_zero$year=="2017"])

summary(rates_over_zero$days_to_first_encounter[rates_over_zero$year=="2025"])


perm.anova(days_to_first_encounter~year,  data=rates_over_zero, 
           nperm = 1999,
           progress = TRUE)


perm.anova(days_to_first_encounter~year,  data=filter(rates_over_zero,species=="PYN"), 
           nperm = 1999,
           progress = TRUE)



summary(rates_over_zero$min_dist)

perm.anova(min_dist~year,  data=rates_over_zero, 
           nperm = 1999,
           progress = TRUE)


perm.anova(min_dist~year,  data=filter(rates_over_zero,species=="PYN"), 
           nperm = 1999,
           progress = TRUE)



ggplot(last_positions,aes(as.factor(year),time_proportion))+
  geom_boxplot(fill="grey50")+
  theme_bw()+
  xlab("year")+ylab("proportion of days")+
  ggtitle(label="a. Tracking duration (proxy for survival)")+
  
  ggplot(last_positions,aes(as.factor(year),npp))+
  geom_boxplot(fill="grey50")+
  theme_bw()+
  xlab("year")+ylab("mg/m3")+
  ggtitle(label="b. net primary productivity")+
  
  
  ggplot(last_positions,aes(as.factor(year),ce))+
  geom_boxplot(fill="grey50")+
  theme_bw()+
  xlab("year")+ylab("")+
  ggtitle(label="c. current effect")+  
  
  
  
  ggplot(filter(last_positions,encounter_rate>0),aes(as.factor(year),encounter_rate))+
  geom_boxplot(fill="grey50")+
  theme_bw()+scale_y_log10()+
  xlab("year")+ylab("%")+
  ggtitle(label="d. Encounter rate (excluding zeros)")+
  
  ggplot(last_positions,aes(as.factor(year),days_to_first_encounter))+
  geom_boxplot(fill="grey50")+
  theme_bw()+
  xlab("year")+ylab("days")+
  ggtitle(label="e. first encounter with vessels")+
  
  ggplot(filter(last_positions,encounter_rate>0),aes(as.factor(year),min_dist))+
  geom_boxplot(fill="grey50")+
  theme_bw()+
  xlab("year")+ylab("km")+
  ggtitle(label="f. minimum distance to vessels")



#ggsave("boxplot_per_year.pdf", width = 10, height = 8, dpi = 600)



# SPATIAL HOTSPOTS

# Load coastline data 
antarctica <- st_as_sf(vect("ADD/add_coastline_medium_res_polygon_v7_10.shp"))
antarctica <- st_transform(antarctica, 4326)
plot(antarctica)

head(all_proximity)

last_positions.a <- all_proximity %>%
  group_by(id, year, species) %>%
  arrange(date) %>%
  summarise(
    last_date = max(date),
    last_lon = last(lon),
    last_lat = last(lat),
    tracking_duration_days = as.numeric(difftime(max(date), min(date), units = "days")),
    n_locations = n()
  ) %>%
  ungroup()

# Colony reference points
kopaitic <- data.frame(lat = -63.3, lon = -57.9)
shirreff <- data.frame(lat = -62.458889, lon = -60.788611)
cierva   <- data.frame(lat = -64.143,   lon = -60.984)
adbay    <- data.frame(lat = -62.175,   lon = -58.446)

# Combined map using patchwork
(ggplot() +
    geom_point(data = filter(all_proximity, year == "2017" & species != "PYD"),
               aes(x = lon, y = lat), size = 1, colour = "grey50", alpha = 0.1) +
    geom_point(data = filter(all_proximity, encounter == TRUE & year == "2017" & species != "PYD"),
               aes(x = lon, y = lat), size = 3, alpha = 0.05, colour = "red3") +
    geom_sf(data = antarctica, fill = "gray90", color = "gray70", size = 0.2) +
    coord_sf(xlim = c(-65, -55), ylim = c(-65, -62)) +
    labs(title = "a. Penguin-Vessel Encounters year 2017",
         x = "Longitude", y = "Latitude") +
    theme_bw() +
    geom_point(data = filter(last_positions.a, year == "2017" & species != "PYD"), 
               aes(x = last_lon, y = last_lat), size = 3) +
    geom_point(data = shirreff, aes(lon, lat), shape = 18, size = 5, colour = "blue3")+
    geom_point(data = cierva, aes(lon, lat), shape = 18, size = 5, colour = "blue3")) /
  
  (ggplot() +
     geom_point(data = filter(all_proximity, year == "2025"),
                aes(x = lon, y = lat), size = 1, colour = "grey50", alpha = 0.1) +
     geom_point(data = filter(all_proximity, encounter == TRUE & year == "2025"),
                aes(x = lon, y = lat), size = 3, alpha = 0.05, colour = "red3") +
     geom_sf(data = antarctica, fill = "gray90", color = "gray70", size = 0.2) +
     coord_sf(xlim = c(-64, -35), ylim = c(-65, -58)) +
     labs(title = "b. Penguin-Vessel Encounters year 2025",
          x = "Longitude", y = "Latitude") +
     theme_bw() +
     geom_point(data = filter(last_positions.a, year == "2025"), 
                aes(x = last_lon, y = last_lat), size = 3) +
     geom_point(data = kopaitic, aes(lon, lat), shape = 18, size = 5, colour = "blue3"))

#ggsave("c:/AWRF/PenguinFledglings/encounters_map_chins_gentoo.png", width = 8, height = 8, dpi = 600)


# ----- 7. Survival analysis -----

# Compute individual-level metrics (already done earlier, but included here for completeness)
last_positions <- all_proximity %>%
  group_by(id, year, species) %>%
  arrange(date) %>%
  summarise(
    deployment_date = min(date),
    deployment_doy = yday(min(date)),
    last_date = max(date),
    tracking_duration_days = as.numeric(difftime(max(date), min(date), units = "days")),
    first_encounter_date = min(date[encounter == TRUE], na.rm = TRUE),
    first_encounter_doy = yday(first_encounter_date),
    days_to_first_encounter = as.numeric(difftime(first_encounter_date, deployment_date, units = "days")),
    encountered_fisheries = any(encounter == TRUE, na.rm = TRUE),
    n_locations = n(),
    n_encounters = sum(encounter == TRUE, na.rm = TRUE),
    encounter_rate = n_encounters / n_locations * 100,
    min_dist = min(min_distance_km, na.rm = TRUE),
    npp = mean(npp, na.rm = TRUE),
    ce = mean(ce, na.rm = TRUE),
    cs = mean(cs, na.rm = TRUE),
    .groups = "drop"
  )

# Survival status: 1 = "died" (tracking duration below 60th percentile)
surv_data <- last_positions %>%
  mutate(status = ifelse(tracking_duration_days < quantile(tracking_duration_days, 0.5, na.rm = TRUE), 1, 0))

# Fit Kaplan‑Meier curves
km_fit <- survfit(Surv(tracking_duration_days, status) ~ year, data = surv_data)

# Plot KM curves (corrected colours for two years)
ggsurvfit(km_fit, linetype_aes = TRUE) +
  add_confidence_interval() +
  labs(x = "Days since fledging", y = "Survival Probability") +
  theme_bw() + xlim(0, 40) +
  scale_colour_manual(values = c("red2", "green3"), name = "Year") +
  scale_fill_manual(values = c("red2", "green3"), name = "Year") +
  scale_linetype_manual(values = c("solid", "dotted"), name = "Year") +
  theme(legend.position = "bottom")

#ggsave("km_fit_survival.pdf", width = 8, height = 6, dpi = 600)

# Jackknife sensitivity analysis (leave‑one‑out)

# Function to generate full survival curve data (only one definition, correct)
get_surv_curve <- function(data) {
  fit <- survfit(Surv(tracking_duration_days, status) ~ year, data = data)
  surv_summary <- summary(fit)
  data.frame(
    time = surv_summary$time,
    surv = surv_summary$surv,
    year = gsub("year=", "", surv_summary$strata)
  )
}

# Original curve
original_curve <- get_surv_curve(surv_data)

# Jackknife resampling (leave-one-out)
jackknife_curves <- purrr::map(unique(surv_data$id), ~{
  get_surv_curve(surv_data[surv_data$id != .x, ])
}) %>% 
  bind_rows(.id = "iteration")

# Mean jackknife curve with 95% CI
mean_jackknife <- jackknife_curves %>%
  group_by(year, time) %>%
  summarise(
    mean_surv = mean(surv),
    lower_ci = quantile(surv, 0.025, na.rm = TRUE),
    upper_ci = quantile(surv, 0.975, na.rm = TRUE),
    .groups = "drop"
  )

# Plot jackknife sensitivity with confidence ribbon
ggplot() +
  geom_ribbon(data = mean_jackknife,
              aes(x = time, ymin = lower_ci, ymax = upper_ci, fill = year), alpha = 0.2) +
  geom_line(data = original_curve, 
            aes(x = time, y = surv, colour = year, linetype = year), 
            size = 1) +
  labs(title = "Survival Curve Sensitivity Analysis",
       x = "Days since fledging",
       y = "Survival Probability",
       colour = "Year", fill = "Year", linetype = "Year") +
  theme_bw() + xlim(0, 30) +
  theme(legend.position = "bottom")

ggsave("leave_one_out_jacknife_survival.pdf", width = 6, height = 6, dpi = 600)


# 
# -----8. Behavioral states using EMbC-----
# 

embc_data <- all_proximity %>%
  dplyr::select(date, lon, lat, id) %>%
  mutate(datetime = as.POSIXct(date, format = "%Y-%m-%d %H:%M:%S"))

result <- stbc(embc_data, stdv = c(0.1, 5 * pi/180), spdLim = 20, smth = 0, 
               maxItr = 200, info = 1)
all_proximity$behavior <- result@A
all_proximity$speed <- result@X[, 1]
all_proximity$turn <- result@X[, 2]

# Assign behavior names (example mapping; adjust as needed)
all_proximity <- all_proximity %>%
  mutate(behavior_name = case_when(
    behavior == 1 ~ "Commuting / Transit",
    behavior == 2 ~ "Foraging (ARS)",
    behavior == 3 ~ "Commuting / Transit",
    behavior == 4 ~ "Exploratory",
    behavior == 5 ~ "Stationary / Resting",
    TRUE ~ "Unknown"
  ))


table(all_proximity$species)

# Filter to the two main behaviors for analysis
behaviors <- all_proximity %>%
  filter(behavior_name %in% c("Commuting / Transit", "Foraging (ARS)")) %>%
  mutate(
    species = dplyr::recode(species, "PYN" = "CHP", "PYD" = "ADP","PYP" = "GEP"),
    bebin = ifelse(behavior_name == "Foraging (ARS)", 1, 0), # behavior binary
    encounter_category = as.factor(ifelse(encounter, "TRUE", "FALSE"))
  )

table(behaviors$species)

gep<-filter(behaviors,species=="GEP")
chp<-filter(behaviors,species=="CHP")
adp<-filter(behaviors,species=="ADP")

table(adp$id,adp$encounter_category)


# GLM and Pemrutations 
set.seed(123)

# Function for permutation LRT on a single term
test_term_perm <- function(data, formula_full, term_to_test, n_perm = 1999) {
  set.seed(123)
  
  # Create reduced formula (remove the term)
  terms_full <- attr(terms(formula_full), "term.labels")
  terms_red <- terms_full[terms_full != term_to_test]
  
  if(length(terms_red) == 0) {
    formula_red <- as.formula("bebin ~ 1")
  } else {
    formula_red <- as.formula(paste("bebin ~", paste(terms_red, collapse = " + ")))
  }
  
  # Fit models
  fit_full <- glm(formula_full, data = data, family = binomial)
  fit_red <- glm(formula_red, data = data, family = binomial)
  
  obs_lrt <- deviance(fit_red) - deviance(fit_full)
  
  # Permutation
  perm_lrt <- replicate(n_perm, {
    data_perm <- data
    data_perm$bebin <- sample(data_perm$bebin)
    
    fit_full_perm <- glm(formula_full, data = data_perm, family = binomial)
    fit_red_perm <- glm(formula_red, data = data_perm, family = binomial)
    
    deviance(fit_red_perm) - deviance(fit_full_perm)
  })
  
  p_value <- mean(perm_lrt >= obs_lrt, na.rm = TRUE)
  
  return(list(term = term_to_test, obs_lrt = obs_lrt, p_value = p_value))
}

# Test all three terms for each species
test_terms <- c("npp", "encounter", "npp:encounter")


model_adp<-glm(bebin ~ npp * encounter,data=adp,family = "binomial")

car::Anova(model_adp)


model_chp<-glm(bebin ~ npp * encounter,data=chp,family = "binomial")

car::Anova(model_chp)
summary(model_chp)

# odds rate

coef_chp <- coef(model_chp)
odds_ratios <- exp(coef_chp)
ci <- exp(confint(model_chp))

# Combine into a table
results_chp <- data.frame(
  Term = names(coef_chp),
  Coefficient = coef_chp,
  Odds_Ratio = odds_ratios,
  CI_2.5 = ci[, 1],
  CI_97.5 = ci[, 2],
  P_value = summary(model_chp)$coefficients[, 4]
)

print(results_chp)


model_gep<-glm(bebin ~ npp * encounter,data=gep,family = "binomial")

car::Anova(model_gep)


coef_gep <- coef(model_gep)
odds_ratios <- exp(coef_gep)
ci <- exp(confint(model_gep))

# Combine into a table
results_gep <- data.frame(
  Term = names(coef_gep),
  Coefficient = coef_gep,
  Odds_Ratio = odds_ratios,
  CI_2.5 = ci[, 1],
  CI_97.5 = ci[, 2],
  P_value = summary(model_gep)$coefficients[, 4]
)

print(results_gep)


# Adélie
results_adp <- lapply(test_terms, function(term) {
  test_term_perm(adp, bebin ~ npp * encounter, term, n_perm = 1999)
})

# Chinstrap
results_chp <- lapply(test_terms, function(term) {
  test_term_perm(chp, bebin ~ npp * encounter, term, n_perm = 1999)
})

# Gentoo
results_gep <- lapply(test_terms, function(term) {
  test_term_perm(gep, bebin ~ npp * encounter, term, n_perm = 1999)
})

# Format results
make_results_table <- function(results, species_name) {
  data.frame(
    species = species_name,
    term = sapply(results, function(x) x$term),
    p_value = sapply(results, function(x) x$p_value)
  )
}

all_results <- rbind(
  make_results_table(results_adp, "Adélie"),
  make_results_table(results_chp, "Chinstrap"),
  make_results_table(results_gep, "Gentoo")
)

print(all_results)


behaviors<-behaviors%>%
  mutate(spkf=paste(species,encounter_category))
 # plot 

ggplot(filter(behaviors,species!="ADP"),aes(npp,bebin,
                                            colour=spkf,linetype=spkf))+
  stat_smooth(method="glm", method.args = list(family = "binomial"),
              fullrange = T)+
  
  #facet_wrap(species~.,scales="free")+
  xlab("Net primary productivity (mg/m3)")+
  ylab("Foraging probability")+
  scale_colour_manual(values=c("steelblue","red3","green3","orange2"),name="")+
  scale_linetype_manual(values=c("dashed","solid","dotted","dotdash"),name="")+
  theme_bw()+theme(legend.position = "bottom")


#ggsave("glm_plots.pdf", width = 5, height = 5, dpi = 600)

# 
# ----- 9. Animation -----
# summary()

fishery2017<-readRDS("GFWData/processed_yearly/vessels_2017.rds")%>%
  filter(timestamp >= (min(all_proximity$date[all_proximity$year=="2017"]))&
           timestamp <= (max(all_proximity$date[all_proximity$year=="2017"])))


fishery2025<-readRDS("GFWData/processed_yearly/vessels_2025.rds")%>%
  filter(timestamp >= (min(all_proximity$date[all_proximity$year=="2025"]))&
           timestamp <= (max(all_proximity$date[all_proximity$year=="2025"])))

summary(fishery2017$timestamp)

summary(fishery2025$timestamp)


fishery<-fishery2017%>%
  rbind(fishery2025)%>%
  mutate(type=c("vessel"),year=year(timestamp))%>%
  select(type,id=ship_id,lon,lat,date=timestamp)

head(fishery)


summary(all_proximity$date)


merged<-all_proximity%>%
  mutate(type=c("Penguin"))%>%
  select(type,id,lon,lat,date)%>%
  rbind(fishery)%>%
  mutate(year = year(date))%>%
  mutate(day=as.POSIXct(substring(date,first=1, last=10), 
                        format = "%Y-%m-%d"))


an17<-(ggplot() +
         
         
         geom_sf(data = antarctica, 
                 fill = "gray90", color = "black") +
         geom_point(data=filter(merged,year=="2017"),
                    aes(x = lon, y=lat,
                        colour = type,shape=type),size=2)+
         coord_sf(xlim =c(-65,-52.5),
                  ylim = c(-65,-60)) +
         
         
         geom_point(data=shirreff,aes(lon,lat),shape=18,size=5,colour="blue3")+
         geom_point(data=cierva,aes(lon,lat),shape=18,size=5,colour="blue3")+
         geom_point(data=adbay,aes(lon,lat),shape=18,size=5,colour="blue3")+
         
         xlab("Longitude")+ylab("Latitude")+ggtitle(label="year 2017")+
         theme_bw()+
         
         scale_colour_manual(values=c("red2","grey30"),name="")+
         scale_shape_manual(values=c(16,17),name=""))

an25<-(ggplot() +
         
         
         geom_sf(data = antarctica, 
                 fill = "gray90", color = "black") +
         geom_point(data=filter(merged,year=="2025"),
                    aes(x = lon, y=lat,
                        colour = type,shape=type),size=2)+
         coord_sf(xlim =c(-65,-40),
                  ylim = c(-65,-58)) +
         
         geom_point(data=kopaitic,aes(lon,lat),shape=18,size=5,colour="blue3")+
         
         
         xlab("Longitude")+ylab("Latitude")+ggtitle(label="year 2025")+
         theme_bw()+
         
         scale_colour_manual(values=c("red2","grey30"),name="")+
         scale_shape_manual(values=c(16,17),name="")  
)



anim17<-an17+transition_time(day)


animation17<-gganimate::animate(anim17,detail=100,res=100,
                                width = 20,height= 15,units='cm',
                                duration=30)
animation17

anim_save(file="Supplementary_Video_S1.gif",animation=animation17)


anim25<-an25+transition_time(day)
animation25<-gganimate::animate(anim25,detail=100,res=100,
                                width = 20,height= 15,units='cm',
                                duration=30)
animation25

anim_save(file="Supplementary_Video_S2.gif",animation=animation25)



# ============================================================================
# End of script
# ============================================================================