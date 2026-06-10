# Penguin tracking, environmental data, 
# vessel proximity, survival analysis, and behavioral states
# 

# Clean memory and load libraries 
gc()

library(dplyr)
library(lubridate)
library(sf)
library(ggplot2)
library(patchwork)
library(terra) 
library(raster)
library(ncdf4)
library(maps)
library(mapdata)
library(matrixStats)
library(purrr)
library(aniMotum)
library(momentuHMM)
library(amt)
library(brms)
library(bayesplot)
library(nlme)
library(car)
library(gganimate)
library(tidyr)
library(EMbC)
library(geosphere)
library(ggnewscale)

# 
# -----Load and process tracking data-----
# 

df <- read.csv("ArgosData_2025_07_24_16_53_35.csv")
df$TimeStamp <- as.POSIXct(df$Loc..date, format = "%d-%m-%Y %H:%M:%S", tz = "GMT")
df$DayStamp <- as.POSIXct(substring(df$Loc..date, 1, 10), format = "%d-%m-%Y", tz = "GMT")
head(df)
df <- (df[!duplicated(df[, c("Platform.ID.No.", "TimeStamp")]), ]) # eliminate duplicates

# Hinke et al. 2020 data
jh <- read.csv("tracks.csv")
pyn <- jh
pyn$TimeStamp <- as.POSIXct(pyn$LOC_DATE, format = "%m/%d/%Y %H:%M:%S", tz = "GMT")
pyn$DayStamp <- as.POSIXct(substring(pyn$LOC_DATE, 1, 10), format = "%m/%d/%Y", tz = "GMT")
pyn$year <- year(pyn$TimeStamp)
pyn <- filter(pyn, year != "2018" & SPP=="PYN")

# Prepare IDs
ids2017 <- pyn %>%
  group_by(SPP, PTT) %>%
  summarise(N = length(LATITUDE), .groups = "drop") %>%
  select(species = SPP, id = PTT, N)

ids2025 <- df %>%
  group_by(Platform.ID.No.) %>%
  summarise(N = length(Longitude), .groups = "drop") %>%
  select(id = Platform.ID.No., N) %>%
  mutate(species = "PYN")

all_ids <- rbind(ids2017, ids2025)

# Combine raw tracking data
jhdf <- data.frame(id = pyn$PTT, DayStamp = pyn$DayStamp, date = pyn$TimeStamp,
                   year = year(pyn$DayStamp), lon = pyn$LONGITUDE, lat = pyn$LATITUDE, lc = 1)
tsdf <- data.frame(id = df$Platform.ID.No., DayStamp = df$DayStamp, date = df$TimeStamp,
                   year = year(df$DayStamp), lon = df$Longitude, lat = df$Latitude, lc = df$Loc..quality)
adf <- rbind(jhdf, tsdf)
adf <- filter(adf, id != "165242")   # remove problem ID

# Prepare tracks for aniMotum
tracks <- data.frame(id = adf$id, date = adf$date, lc = as.factor(adf$lc), lon = adf$lon, lat = adf$lat)
tracks <- tracks[!duplicated(tracks[, c("id", "date")]), ]
tracks <- subset(tracks, id != "255211")   # only one day of data

# Fit SSM
fit.crw <- fit_ssm(tracks, model = "crw", time.step = 0.5,
                   vmax = 15,spdf=TRUE)

aniMotum::map(fit.crw,by.id=F,what="predicted")

fit_rerouted <- route_path(
  fit.crw,               
  what = "predicted",    
  dist = 1000            # Buffer distance in meters
)

aniMotum::map(fit_rerouted,by.id=F,what="predicted")

res.crw <- osar(fit_rerouted)
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

ggsave("supplemental_figure_S1.png", 
       width = 8, height = 8, dpi = 600)

# Extract  predicted tracks

track.predicted <- aniMotum::grab(fit_rerouted, what = "rerouted",group = TRUE)
track.predicted$DayStamp <- as.POSIXct(paste(year(track.predicted$date),
                                             month(track.predicted$date),
                                             day(track.predicted$date),sep="-"), format = "%Y-%m-%d")


track.predicted$DayStamp <- as.POSIXct(paste(year(track.predicted$date),
                                             month(track.predicted$date),
                                             day(track.predicted$date),sep="-"), format = "%Y-%m-%d")

track.2017 <- filter(track.predicted, year(DayStamp) == 2017)
track.2025 <- filter(track.predicted, year(DayStamp) == 2025)

summary(track.2017)

# Load coastline data 
antarctica <- st_as_sf(vect("ADD/add_coastline_medium_res_polygon_v7_10.shp/add_coastline_medium_res_polygon_v7_10.shp"))
antarctica <- st_transform(antarctica, 4326)
plot(antarctica)

# Colony reference points
kopaitic <- data.frame(lat = -63.3, lon = -57.9)
shirreff <- data.frame(lat = -62.458889, lon = -60.788611)
cierva   <- data.frame(lat = -64.143,   lon = -60.984)

ggplot()+
  geom_point(data=track.2025,aes(lon,lat),colour="steelblue")+
  geom_point(data=track.2017,aes(lon,lat),colour="darkred")+   
  geom_sf(data = antarctica, fill = "gray90", color = "gray70", size = 0.2) +
  coord_sf(xlim = c(-65, -10), ylim = c(-65, -54))+
  geom_point(data = shirreff, aes(lon, lat), shape = 18, size = 5, colour = "black")+
  geom_point(data = cierva, aes(lon, lat), shape = 18, size = 5, colour = "black")+
  geom_point(data = kopaitic, aes(lon, lat), shape = 18, size = 5, colour = "black")+
  theme_bw()

# 
# Environmental variables 
# 
# function to safely extract values, handling missing dates
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

ncpath <- "C:/Fledge/Copernicus/"

# 2017 npp
ncname2 <- "cmems_mod_glo_bgc_my_0.25deg_P1D-m_1779109435194"
ncfname2 <- paste0(ncpath, ncname2, ".nc")
tmp_brick2 <- brick(ncfname2, varname = "nppv")
track.2017$npp <- extract_env(track.2017, tmp_brick2, "nppv")

summary(track.2017$npp)

# 2025 npp
ncpath <- "C:/AWRF/CHL/"
ncname1 <- "cmems_mod_glo_bgc_my_0.25deg_P1D-m_1779109338425"
ncfname1 <- paste0(ncpath, ncname1, ".nc")
tmp_brick1 <- brick(ncfname1, varname = "nppv")
track.2025$npp <- extract_env(track.2025, tmp_brick1, "nppv")

summary(track.2025$npp)

# 
# Sea water velocity – uo (eastward) and vo (northward)
# 

# 2025 uo
ncname_uo25 <- "cmems_mod_glo_phy_my_0.083deg_P1D-m_1779109823356"
ncfname_uo25 <- paste0(ncpath, ncname_uo25, ".nc")
brick_uo25 <- brick(ncfname_uo25, varname = "uo")
track.2025$cuo <- extract_env(track.2025, brick_uo25, "uo")

# 2017 uo
ncname_uo17 <- "cmems_mod_glo_phy_my_0.083deg_P1D-m_1779109673336"
ncfname_uo17 <- paste0(ncpath, ncname_uo17, ".nc")
brick_uo17 <- brick(ncfname_uo17, varname = "uo")
track.2017$cuo <- extract_env(track.2017, brick_uo17, "uo")

# 2025 vo
brick_vo25 <- brick(ncfname_uo25, varname = "vo")   # same file, different variable
track.2025$cvo <- extract_env(track.2025, brick_vo25, "vo")

# 2017 vo
brick_vo17 <- brick(ncfname_uo17, varname = "vo")
track.2017$cvo <- extract_env(track.2017, brick_vo17, "vo")

# 
# Merge tracks and calculate water speed
# 

tracks <- rbind(track.2017, track.2025)

tracks$current.dir <- atan2(tracks$cuo, tracks$cvo)
summary(tracks$current.dir)

tracks$current.speed <- sqrt(tracks$cuo^2 + tracks$cvo^2)
summary(tracks$current.speed)

tracks$year <- year(tracks$date)

# Combine with species info
all_tracks <- tracks %>%
  select(id, lon, lat, npp, cd = current.dir, 
         cs = current.speed, date, year) %>%
  mutate(timestamp=date)


# Proximity to vessels 


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
  results$closest_vessel_distance_km <- NA
  
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

# process fishing data before applying the proximity function

# fisheries processing 

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
      # Check if required columns exist 
      
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
          select(ship_id, timestamp, year, date_only, lon, lat, everything())
        
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
process_vessels_by_year <- function(zip_folder_path, output_folder = "C:/AWRF/Fledge/GFWData/processed_yearly/", 
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

# Function to load all yearly vessel data
load_all_vessel_data <- function(processed_folder = "C:/AWRF/Fledge/GFWData/processed_yearly") {
  
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

# Load all data 
all_vessels <- load_all_vessel_data("C:/AWRF/Fledge/GFWData/processed_yearly")%>%
  filter(lat<(-50))%>%
  filter(lon<(-10) & lon>(-80))

all_vessels<-all_vessels[!duplicated(
  all_vessels[, c("timestamp","ship_id")]), ]

all_vessels<-all_vessels%>%
  na.omit()%>%
  group_by(ship_id)%>%
  filter(n() >= 10) %>%
  ungroup()

summary(all_vessels$speed)

ggplot(filter(all_vessels,speed<6),aes(lon,lat))+geom_point()



# now run proximity analysis
vessel_path <- "C:/AWRF/Fledge/GFWData/"
all_proximity <- calculate_proximity_continuous(all_tracks, vessel_path, 50, 1)

# day since fledging
all_proximity <- all_proximity %>%
  mutate(day_month_date=as.Date(paste("2000", 
                                      month(date), 
                                      day(date), 
                                      sep = "-")))%>%
  group_by(id) %>%
  arrange(day_month_date) %>%
  mutate(
    # Day 1 = first day of tracking for that individual
    day_since_fledging = as.numeric(day_month_date - min(day_month_date)) + 1
  ) %>%
  ungroup()

# Save full results
saveRDS(all_proximity, "all_proximity_chinstrap.rds")

table(all_proximity$colony_name)

# Compute individual-level metrics 
last_positions <- all_proximity %>%
  group_by(id, year,colony_name) %>%
  arrange(date) %>%
  summarise(
    deployment_date = min(date),
    deployment_doy = yday(min(date)),
    last_date = max(date),
    last_lon = last(lon),
    last_lat = last(lat),
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
    cs = mean(cs, na.rm = TRUE),
    .groups = "drop"
  )

last_positions<-last_positions %>%
  group_by(year)%>%
  mutate(time_proportion=tracking_duration_days/max(tracking_duration_days),
         cskm=cs*3.6)

# Save full results
saveRDS(last_positions, "last_positions.rds")

summary(last_positions$tracking_duration_days)
summary(last_positions$tracking_duration_days[last_positions$year=="2017"])

summary(last_positions$tracking_duration_days[last_positions$year=="2025"])

median(last_positions$tracking_duration_days[last_positions$year=="2025"])/
  median(last_positions$tracking_duration_days[last_positions$year=="2017"])

summary(last_positions$npp[last_positions$year=="2017"])

summary(last_positions$npp[last_positions$year=="2025"])

summary(last_positions$npp[last_positions$year=="2017"])

summary(last_positions$npp[last_positions$year=="2025"])

summary(last_positions$cskm[last_positions$year=="2017"])

summary(last_positions$cskm[last_positions$year=="2025"])

# Summaries (species, yearly, individual)

all_proximity %>%
  summarise(
    positions = n(),
    encounters = sum(encounter),
    encounter_rate = round(encounters / positions * 100, 4),
    mean_distance_km = round(mean(min_distance_km, na.rm = TRUE), 2),
    .groups = "drop"
  ) 

yearly_summary <- all_proximity %>%
  group_by(year) %>%
  summarise(
    positions = n(),
    encounters = sum(encounter),
    encounter_rate = round(encounters / positions * 100, 2),
    .groups = "drop"
  )

yearly_summary

length(last_positions$id[last_positions$encounter_rate>0])/25

summary(last_positions$encounter_rate[last_positions$year=="2017"])
summary(last_positions$encounter_rate[last_positions$year=="2025"])

summary(last_positions$days_to_first_encounter[last_positions$year=="2017" & 
                                                 last_positions$days_to_first_encounter !="Inf"])

summary(last_positions$days_to_first_encounter[last_positions$year=="2025"& 
                                                 last_positions$days_to_first_encounter !="Inf"])

summary(last_positions$min_dist[last_positions$year=="2017" & 
                                  last_positions$days_to_first_encounter !="Inf"])

summary(last_positions$min_dist[last_positions$year=="2025"& 
                                  last_positions$days_to_first_encounter !="Inf"])

# Save summaries
output_path <- "C:/AWRF/Proximity_Results_Final/"
if(!dir.exists(output_path)) dir.create(output_path, recursive = TRUE)

write.csv(yearly_summary, file.path(output_path, "yearly_summary_final.csv"), row.names = FALSE)


#Identify individuals that had at least one encounter
encountered_ids <- all_proximity %>%
  filter(encounter == TRUE) %>%
  distinct(id) %>%
  pull(id)

if(!"event" %in% colnames(all_proximity)) {
  # Calculate median tracking duration per year
  median_duration <- all_proximity %>%
    group_by(year, id) %>%
    summarise(max_day = max(day_since_fledging, na.rm = TRUE), .groups = "drop") %>%
    group_by(year) %>%
    summarise(median_duration = median(max_day, na.rm = TRUE), .groups = "drop")
  
  all_proximity <- all_proximity %>%
    left_join(median_duration, by = "year") %>%
    group_by(id) %>%
    mutate(
      max_day_id = max(day_since_fledging, na.rm = TRUE),
      event = ifelse(max_day_id < median_duration, 1, 0)
    ) %>%
    ungroup() %>%
    select(-median_duration, -max_day_id)
}

### distance from colony 


# Colony reference points (for later)
kopaitic <- data.frame(lat = -63.3, lon = -57.9)
shirreff <- data.frame(lat = -62.458889, lon = -60.788611)
cierva   <- data.frame(lat = -64.143,   lon = -60.984)


# Define colony locations as sf points
colony_locations <- data.frame(
  colony_name = c("shirreff", "cierva", "kopaitic"),
  colony_lon = c(-60.788611, -60.984, -57.9),
  colony_lat = c(-62.458889, -64.143, -63.3)
)

# Convert colonies to sf
colonies_sf <- colony_locations %>%
  st_as_sf(coords = c("colony_lon", "colony_lat"), crs = 4326)

# AUTOMATIC COLONY ASSIGNMENT USING SF 

# Get first position for each penguin
first_positions <- all_proximity %>%
  group_by(id) %>%
  arrange(date) %>%
  slice(1) %>%
  ungroup() %>%
  select(id, first_lon = lon, first_lat = lat, first_date = date)

# Convert first positions to sf
first_positions_sf <- first_positions %>%
  st_as_sf(coords = c("first_lon", "first_lat"), crs = 4326)

# Calculate distances to each colony and assign closest
penguin_colony_assignment <- data.frame()

for(i in 1:nrow(first_positions_sf)) {
  # Calculate distances to all colonies
  distances <- st_distance(first_positions_sf[i,], colonies_sf)
  min_dist_idx <- which.min(distances)
  
  penguin_colony_assignment <- bind_rows(penguin_colony_assignment, 
                                         data.frame(
                                           id = first_positions_sf$id[i],
                                           colony_name = colony_locations$colony_name[min_dist_idx],
                                           distance_to_nearest_colony = 
                                             as.numeric(distances[min_dist_idx]) / 1000,
                                           first_lon = first_positions$first_lon[i],
                                           first_lat = first_positions$first_lat[i],
                                           first_date = first_positions$first_date[i]
                                         )
  )
}

print("=== COLONY ASSIGNMENT ===")
print(penguin_colony_assignment)

# add colony info to all_proximity
all_proximity <- all_proximity %>%
  left_join(penguin_colony_assignment %>% select(id, colony_name), by = "id")

all_proximity <- all_proximity %>%
  left_join(colony_locations, by = "colony_name")

# calculate distances using geosphere 
all_proximity <- all_proximity %>%
  mutate(
    distance_from_colony = distm(
      cbind(lon, lat),
      cbind(colony_lon, colony_lat),
      fun = distHaversine
    )[,1] / 1000,
    date = as.Date(date)
  )


# Behavioral states using HMMI 
# 

head(all_proximity)

# Ensure proper column structure
df_clean <- all_proximity %>%
  # Rename critical columns
  rename(ID = id) %>%  # momentuHMM requires 'ID' column
  arrange(ID, timestamp) %>%
  group_by(ID) %>%
  mutate(
    time = as.numeric(difftime(timestamp, min(timestamp), units = "hours"))) %>%
  ungroup()

head(df_clean)

summary(df_clean$cs)

# Create properly structured dataframe

# create  dataframe with only required elements
hmm_data <- prepData(
  data.frame(
    ID = df_clean$ID,
    x = df_clean$lon,
    y = df_clean$lat,
    time = as.numeric(df_clean$time),
    mindist=df_clean$min_distance_km,
    encounter=df_clean$encounter,
    npp=df_clean$npp,cs=df_clean$cs
  ),
  type = "LL"
)

ggplot(hmm_data,aes(time))+geom_histogram()

# calculate steps
hmm_data_fixed <- hmm_data %>%
  group_by(ID) %>%
  mutate(
    # Calculate step length in meters using proper haversine formula
    step_m = c(NA, distHaversine(
      cbind(x[1:(n()-1)], y[1:(n()-1)]),
      cbind(x[2:n()], y[2:n()])
    )),
    # Convert to km (more interpretable)
    step_km = step_m / 1000,
    # Calculate turning angle (in radians, between -pi and pi)
    angle = c(NA, atan2(
      y[2:n()] - y[1:(n()-1)],
      x[2:n()] - x[1:(n()-1)]
    )),
    # Ensure angle is wrapped correctly
    angle = ifelse(is.na(angle), 0, 
                   atan2(sin(angle), cos(angle)))  # Wrap to [-pi, pi]
  ) %>%
  ungroup() %>%
  # Remove rows with missing steps or angles
  filter(!is.na(step_km), !is.na(angle), step_km > 0)

# Check the distribution of step lengths
summary(hmm_data_fixed$step_km)
hist(hmm_data_fixed$step_km, breaks = 50, main = "Step lengths (km)")

# Check reasonable values for penguins (should be 0.5-10 km per 30 min)
quantile(hmm_data_fixed$step_km, probs = c(0.1, 0.25, 0.5, 0.75, 0.9), na.rm = TRUE)

# Remove steps > 20 km (penguins swim max ~2-3 km per 30 minutes)
hmm_data_clean <- hmm_data_fixed #%>%
#filter(step_km < 20, step_km > 0)  # 20 km is still generous

# Check the cleaned distribution
summary(hmm_data_clean$step_km)
hist(hmm_data_clean$step_km, breaks = 50, main = "Step lengths (km) - Cleaned")

# Now the values should be reasonable
quantile(hmm_data_clean$step_km, probs = c(0.1, 0.25, 0.5, 0.75, 0.9), na.rm = TRUE)

# Prepare cleaned data for HMM
hmm_data_ready <- data.frame(
  ID = hmm_data_clean$ID,
  x = hmm_data_clean$x,
  y = hmm_data_clean$y,
  time = hmm_data_clean$time,
  step = hmm_data_clean$step_km,
  angle = hmm_data_clean$angle,
  encounter = hmm_data_clean$encounter,
  npp = hmm_data_clean$npp,
  cs = hmm_data_clean$cs,
  mindist=hmm_data_clean$mindist
)

hmm_prep <- prepData(
  data.frame(
    ID = hmm_data_ready$ID,
    x = hmm_data_ready$x,
    y = hmm_data_ready$y,
    time = hmm_data_ready$time,
    mindist=scale(hmm_data_ready$mindist),
    encounter=ifelse(hmm_data_ready$encounter==TRUE,1,0),
    npp=scale(hmm_data_ready$npp),cs=scale(hmm_data_ready$cs)
  ),
  type = "LL",
  coordNames = c("x","y"),
  covNames = c("mindist","encounter","npp","cs")
)



# Realistic starting values based on actual data

step_quantiles <- quantile(hmm_data_ready$step, probs = c(0.33, 0.67), na.rm = TRUE)
angle_quantiles <- quantile(hmm_data_ready$angle, probs = c(0.33, 0.67), na.rm = TRUE)

step_quantiles
angle_quantiles

sd(hmm_data_ready$step)
sd(hmm_data_ready$angle)

summary(hmm_prep$npp)
summary(hmm_prep$cs)
summary(hmm_prep$encounter)
summary(hmm_prep$mindist)

hmm_prep$cs[is.na(hmm_prep$cs)]<-0

fit_hmm_final <- fitHMM(
  data = (hmm_prep),
  nbStates = 2,
  dist = list(step = "gamma", angle = "vm"),
  Par0 = list(
    step = c(
      0.7,  # mean step state 1 (foraging - shorter)
      1.2,  # mean step state 2 (transit - longer)
      0.1,  # SD state 1 (relatively small)
      0.3   # SD state 2 (more variable)
    ),
    angle = c(0.01, 0.5)   # concentration: low (foraging), high (transit)
  ),
  formula = ~ (encounter+ npp + cs+mindist),  
  estAngleMean = list(angle = FALSE),
  stateNames = c("Foraging", "Transit")
)

getPar(fit_hmm_final)

#estimated parameters
print(fit_hmm_final)

# Check if any parameters look extreme
fit_hmm_final$mle

# Check if standard errors are huge or missing
fit_hmm_final$se

# Decode most probable state sequence
states <- viterbi(fit_hmm_final)

# Add to data
hmm_data_ready$state <- states
head(hmm_data_ready)

table(as.numeric(hmm_data_ready$ID))
hmm_data_ready$idnum<-as.numeric(hmm_data_ready$ID)

table(hmm_data_ready$ID,hmm_data_ready$idnum)

hmm_data_ready$year<-ifelse(as.numeric(hmm_data_ready$ID)<=8,"2017","2025")

table(hmm_data_ready$year)

subs<-st_as_sf(vect("FAO/FAO_Major_Fishing_Areas/FAO_Major_Fishing_Areas.shp"))
subs<-st_transform(subs,4326)
subs<-subs%>%
  filter(F_AREA=="48" | F_AREA=="47" |F_AREA=="41")

ggplot()+
  geom_sf(data = subs,aes(fill=F_SUBAREA), size = 0.2) +
  geom_sf_text(data = subs, aes(label = F_SUBAREA), 
               size = 3, check_overlap = TRUE) +
  coord_sf(xlim = c(-80, 0), ylim = c(-68, -50))

(ggplot() +
    geom_point(data=filter(hmm_data_ready,year=="2017"), 
               aes(x = x, y = y),colour="grey50",alpha=0.01) +
    geom_point(data=filter(hmm_data_ready,year=="2017" & state=="1"), 
               aes(x = x, y = y),colour="darkgreen",alpha=0.05) +
    
    theme_bw() +
    labs(title = "Behavioral states over time")+
    geom_sf(data = subs, fill = "transparent", color = "gray30", size = 0.2) +
    geom_sf(data = antarctica, fill = "gray90", color = "gray70", size = 0.2) +
    
    coord_sf(xlim = c(-65, -10), ylim = c(-65, -54)) +
    labs(title = "a. year 2017",
         x = "Longitude", y = "Latitude") +
    theme_bw() +
    geom_point(data = filter(last_positions,year=="2017"), 
               aes(x = last_lon, y = last_lat), size = 3) +
    geom_point(data = shirreff, aes(lon, lat), shape = 18, size = 5, colour = "blue3")+
    geom_point(data = cierva, aes(lon, lat), shape = 18, size = 5, colour = "blue3"))/
  
  (ggplot() +
     geom_point(data=filter(hmm_data_ready,year=="2025"), 
                aes(x = x, y = y),colour="grey50",alpha=0.01) +
     geom_point(data=filter(hmm_data_ready,year=="2025" & state=="1"), 
                aes(x = x, y = y),colour="darkgreen",alpha=0.05) +
     
     theme_bw() +
     labs(title = "Behavioral states over time")+
     geom_sf(data = subs, fill = "transparent", color = "gray30", size = 0.2) +
     geom_sf(data = antarctica, fill = "gray90", color = "gray70", size = 0.2) +
     
     
     coord_sf(xlim = c(-65, -10), ylim = c(-65, -54)) +
     labs(title = "b. year 2025",
          x = "Longitude", y = "Latitude") +
     theme_bw() +
     geom_point(data = filter(last_positions,year=="2025"), 
                aes(x = last_lon, y = last_lat), size = 3) +
     geom_point(data = kopaitic, aes(lon, lat), shape = 18, size = 5, colour = "blue3")
  )

ggsave(filename = "C:/AWRF/map_states_years.png",dpi = 600,
       height = 20,width=20,units="cm")

# ids, colonies and first transmission:

data.frame(all_proximity%>%
             group_by(year,id,colony_name)%>%
             dplyr::summarise(start=min(date)))
data.frame(all_proximity%>%
             group_by(year,id,colony_name)%>%
             dplyr::summarise(end=max(date)))

# get state probabilities
state_probs <- stateProbs(fit_hmm_final)

# Add to original data
hmm_data_ready$p_state1 <- state_probs[, 1]
hmm_data_ready$p_state2 <- state_probs[, 2]

head(hmm_data_ready)

head(data.frame(all_proximity))
all_proximity$year<-as.character(all_proximity$year)

prox_states<-all_proximity%>%
  left_join(hmm_data_ready,by=c("id"="ID", 
                                "lon"="x",
                                "lat"="y",
                                "npp","cs","year","encounter"))%>%
  select(id,year,timestamp,time,lon,lat,npp,cs,year,encounter,min_distance_km,mindist,
         state,p_state1,p_state2,step,day_since_fledging,distance_from_colony,
         colony_name,colony_lon,colony_lat)%>%
  na.omit()

plot(prox_states$min_distance_km,prox_states$mindist) # should match perfectly

plot(prox_states$timestamp,prox_states$time)

ind_states<-prox_states%>%
  group_by(id, year) %>%
  arrange(timestamp) %>%
  summarise(
    deployment_date = min(timestamp),
    deployment_doy = yday(min(timestamp)),
    last_date = max(timestamp),
    last_lon = last(lon),
    last_lat = last(lat),
    tracking_duration_days = as.numeric(difftime(max(timestamp), min(timestamp), units = "days")),
    first_encounter_date = min(timestamp[encounter == TRUE], na.rm = TRUE),
    first_encounter_doy = yday(first_encounter_date),
    days_to_first_encounter = as.numeric(difftime(first_encounter_date, deployment_date, units = "days")),
    encountered_fisheries = any(encounter == TRUE, na.rm = TRUE),
    n_locations = n(),
    n_encounters = sum(encounter == TRUE, na.rm = TRUE),
    encounter_rate = n_encounters / n_locations * 100,
    min_dist = min(min_distance_km, na.rm = TRUE),
    npp = mean(npp, na.rm = TRUE),
    cs = mean(cs, na.rm = TRUE),
    prob_foraging= mean(p_state1, na.rm = TRUE),
    se_foraging=sd(p_state1, na.rm = TRUE)/sqrt(n_locations),
    .groups = "drop"
  )

data.frame(ind_states)

# foraging probability

ggplot(ind_states,aes(npp,prob_foraging,colour=encountered_fisheries,
                      shape=encountered_fisheries,
                      linetype=encountered_fisheries))+
  geom_errorbar(aes(ymin=prob_foraging-se_foraging,
                    ymax=prob_foraging+se_foraging))+
  stat_smooth(method="glm", method.args = list(family = "binomial"),
              formula=y~poly(x,2),se=F)+
  geom_label(aes(label=year),size=3)+
  theme_bw()+theme(legend.position = "none")+
  scale_colour_manual(values=c("steelblue","darkred"))+
  scale_linetype_manual(values=c("dashed","solid"))+
  xlab("productivity")+ylab("probability of foraging")+
  ggtitle(label="a.")+
  
  ggplot(filter(ind_states,year=="2017"),
         aes(npp,tracking_duration_days))+
  
  geom_smooth(method="gam", 
              formula=y~s(x,k=4),se=F)+
  geom_point(aes(colour=encountered_fisheries,
                 shape=encountered_fisheries),size=3)+
  theme_bw()+
  scale_y_log10()+
  scale_shape_manual(values=c("circle","triangle"),name="encounter")+
  scale_colour_manual(values=c("steelblue","darkred"),name="encounter")+
  xlab("productivity")+ylab("tracking duration")+
  ggtitle(label="b.")+
  
  ggplot(filter(ind_states,year!="2017"),
         aes(npp,tracking_duration_days))+
  
  geom_smooth(method="gam", 
              formula=y~s(x,k=4),se=F)+
  geom_point(aes(colour=encountered_fisheries,
                 shape=encountered_fisheries),size=3)+
  theme_bw()+
  scale_y_log10()+
  theme(legend.position = "none")+
  scale_colour_manual(values=c("steelblue","darkred"),name="encounter")+
  
  xlab("productivity")+ylab("tracking duration")+
  ggtitle(label="c.")+
  
  ggplot(filter(ind_states,encountered_fisheries==TRUE),
         aes(days_to_first_encounter,tracking_duration_days))+
  
  geom_smooth(method="gam", 
              formula=y~s(x,k=4),se=F)+
  geom_point(size=3,aes(colour=as.factor(year),shape=as.factor(year)))+
  theme_bw()+
  scale_colour_manual(values=c("darkred","steelblue"),name="year")+
  scale_shape_manual(values=c("circle","triangle"),name="year")+
  
  scale_y_log10()+xlab("days to first encounter")+
  ylab("tracking duration")+
  ggtitle(label="d.")

ggsave(filename = "C:/AWRF/states_and_duration.png",dpi = 600,
       height = 15,width=20,units="cm")



# CATCH AND EFFORT ANALYSIS


hauls_sf<-st_as_sf(all_vessels,coords=c("lon","lat"),crs = st_crs(4326))

subs<-st_as_sf(vect("FAO/FAO_Major_Fishing_Areas/FAO_48_subareas.shp"))
subs<-st_transform(subs,4326)

plot(subs)

hauls_sf$subarea <- sapply(st_intersects(hauls_sf, subs), function(x) {
  if(length(x) > 0) subs$F_SUBAREA[x[1]] else NA
})

head(hauls_sf)
table(hauls_sf$subarea)

hauls_year<-data.frame(hauls_sf)%>%
  group_by(subarea,year)%>%
  summarise(effort=n())%>%
  mutate(year=as.integer(year))

catch<-read.csv("Catch (tonnes) - Krill - Area 48.csv")%>%
  select(year=CCAMLR.Season,Area.481,Area.482)%>%
  pivot_longer(cols=c("Area.481","Area.482"))%>%
  select(year,subarea=name,catch=value)

catch$subarea<-ifelse(catch$subarea=="Area.481","48.1","48.2")
catch$catch[is.na(catch$catch)]<-0

head(catch)
head(hauls_year)

catch_effort<-hauls_year%>%
  left_join(catch,by=c("year","subarea"))%>%
  na.omit()

ggplot(catch_effort,aes(effort,catch))+
  geom_smooth(method="lm",formula=y~log(x))+
  geom_point(aes(colour=subarea,shape=subarea))+
  xlab("number of points (effort proxy)")+
  ylab("catch (ton)")+
  theme_bw()+theme(legend.position = "inside",
                   legend.position.inside = c(0.9,0.15))+
  scale_color_manual(values=c("darkred","red3"))

cor(log10(catch_effort$effort),catch_effort$catch)

lm1<-lm(catch~log10(effort),data=catch_effort)
summary(lm1)

ggsave(filename = "catch_and_effort.png",dpi = 600,
       height = 10,width=13,units="cm")

# 
# ANIMATIONS


gc()

fishery2017<-readRDS("GFWData/processed_yearly/vessels_2017.rds")
fishery2025<-readRDS("GFWData/processed_yearly/vessels_2025.rds")

fishery2017<-fishery2017%>%
  filter(
    date_only<(as.Date("2017-03-31"
    )))%>%
  mutate(date=date_only)

fishery2025<-fishery2025%>%
  filter(date_only<(as.Date(max(all_proximity$date)))&
           date_only>(as.Date(min(all_proximity$date[all_proximity$year=="2025"])))&
           lon<(-35))%>%
  mutate(date=date_only)

head(all_proximity)

penguin2017<-all_proximity%>%
  filter(date<(as.Date("2017-03-31")))

penguin2025<-all_proximity%>%
  filter(year=="2025" & lon<(-36))

an17<-(ggplot() +
         
         geom_point(data=penguin2017,
                    aes(x = lon, y=lat,
                        colour = colony_name,shape=colony_name),size=2)+
         geom_point(data=fishery2017,
                    aes(x = lon, y=lat),size=2,colour="grey50")+
         
         geom_sf(data = antarctica, 
                 fill = "gray90", color = "black") +
         coord_sf(xlim =c(-65,-52.5),
                  ylim = c(-65,-60)) +
         
         
         geom_point(data=shirreff,aes(lon,lat),shape=18,size=5,colour="steelblue")+
         geom_point(data=cierva,aes(lon,lat),shape=18,size=5,colour="darkred")+
         
         xlab("Longitude")+ylab("Latitude")+
         theme_bw()+
         
         scale_colour_manual(values=c("darkred","steelblue"),name="")+
         scale_shape_manual(values=c(15,16,17,18),name=""))

an17

an25<-(ggplot() +
         
         geom_point(data=penguin2025,
                    aes(x = lon, y=lat),size=2,colour="darkred")+
         geom_point(data=fishery2025,
                    aes(x = lon, y=lat),size=2,colour="grey50")+
         geom_sf(data = antarctica, 
                 fill = "gray90", color = "black") +
         coord_sf(xlim =c(-65,-40),
                  ylim = c(-65,-58)) +
         
         geom_point(data=kopaitic,aes(lon,lat),shape=18,size=5,colour="darkred")+
         
         
         xlab("Longitude")+ylab("Latitude")+
         theme_bw()+
         
         scale_colour_manual(values=c("steelblue","darkred"),name="")+
         scale_shape_manual(values=c(15,16,17,18),name=""))

an25

anim17<-an17+transition_time(date)+
  labs(title = "a. year 2017: {format(frame_time, '%d %b')}")

animation17<-gganimate::animate(anim17,detail=100,res=100,
                                width = 20,height= 15,units='cm',
                                duration=20)
animation17

anim_save(file="c:/AWRF/Supplementary_Video_S1.gif",animation=animation17)

anim25<-an25+transition_time(date)+
  labs(title = "a. year 2025: {format(frame_time, '%d %b')}")

animation25<-gganimate::animate(anim25,detail=100,res=100,
                                width = 20,height= 15,units='cm',
                                duration=30)
animation25

anim_save(file="c:/AWRF/Supplementary_Video_S2.gif",animation=animation25)



# FISHERY & SEA ICE ANIMATION 

# Read sea ice data
ncpath <- "C:/AWRF/Fledge/SeaIce_and_Fishery/"
ncname1 <- "cmems_mod_glo_phy_my_0.083deg_P1D-m_1779888444246"
ncfname1 <- paste0(ncpath, ncname1, ".nc")
tmp_brick1 <- brick(ncfname1, varname = "siconc")

# Extract coordinates and values - convert to data frame
points_df <- as.data.frame(rasterToPoints(tmp_brick1))
colnames(points_df) <- c("lon", "lat", paste0("sic_", 1:nlayers(tmp_brick1)))

# Get dates from the brick
dates <- getZ(tmp_brick1)  # This gets the time dimension
if(is.null(dates)) {
  # If no dates in the brick, create sequential dates
  dates <- seq(as.Date("2024-12-23"), by = "day", length.out = nlayers(tmp_brick1))
}

# Reshape to long format
sic_df <- points_df %>%
  pivot_longer(cols = starts_with("sic_"),
               names_to = "layer",
               values_to = "sic") %>%
  mutate(layer_num = as.numeric(gsub("sic_", "", layer)),
         date = dates[layer_num]) %>%
  select(lon, lat, date, sic) %>%
  filter(!is.na(sic))  # Remove NA values

# Check the result
head(sic_df)

sic_df<-sic_df%>%
  filter(lat>(-65) & lat<(-58))%>%
  filter(lon>(-65) & lon<(-42.5))%>%
  mutate(dates=as.Date(date))

# Create 0.25° grid cells by rounding coordinates
sic_df_resampled <- sic_df %>%
  mutate(lon_025 = round(lon / 0.25) * 0.25,
         lat_025 = round(lat / 0.25) * 0.25) %>%
  group_by(lon_025, lat_025, date) %>%
  summarise(sic = max(sic, na.rm = TRUE), .groups = 'drop') %>%
  rename(lon = lon_025, lat = lat_025) %>%
  filter(is.finite(sic))  # Remove any infinite values

head(sic_df_resampled)

# read- subarea file

# Load coastline data 

subs<-st_as_sf(vect("C:/GIS/FAO_Major_Fishing_Areas/FAO_48_subareas.shp"))
subs<-st_transform(subs,4326)

plot(subs)

# load fishing vessels data 

# Run the extraction (adjust path as needed)
process_vessels_by_year(
  zip_folder_path = "C:/AWRF/Fledge/SeaIce_and_Fishery/",
  output_folder = "C:/AWRF/Fledge/SeaIce_and_Fishery/processed_yearly/",
  years = c(2024,2025),
  chunk_size = 5  # Process 5 ships at a time
)

# Load all data 
gerlache_vessels <- load_all_vessel_data("C:/AWRF/Fledge/SeaIce_and_Fishery/processed_yearly/")

head(gerlache_vessels)

gv<-gerlache_vessels%>%
  mutate(date=as.Date(date_only))

sic<-sic_df_resampled%>%
  mutate(date=as.Date(date))

summary(sic$sic)

table(sic$date)

ggplot(sic, aes(as.factor(date),sic))+geom_boxplot()

# Get unique dates
unique_dates <- unique(sic$date) %>% sort()

# create frames one by one and save directly
create_animation_sequential <- function(sic, gv, antarctica, subs, 
                                        shirreff, cierva, kopaitic, 
                                        output_file = "animation.gif") {
  
  # Create temporary directory for frames
  temp_dir <- tempdir()
  frame_files <- c()
  
  # Progress bar
  pb <- txtProgressBar(min = 0, max = length(unique_dates), style = 3)
  
  for(i in seq_along(unique_dates)) {
    current_date <- unique_dates[i]
    
    # Filter data for current date
    sic_date <- sic %>% filter(date == current_date)
    gv_date <- gv %>% filter(date == current_date)
    
    # Create plot
    p <- ggplot() +
      geom_raster(data = sic_date, aes(x = lon, y = lat, fill = sic)) +
      scale_fill_gradient(low = "white", high = "blue", 
                          na.value = "transparent",
                          limits = c(0, 1),
                          name = "Sea Ice\nConcentration") +
      geom_sf(data = antarctica, fill = "gray90", color = "gray70", size = 0.2) +
      geom_sf(data = subs, fill = "transparent", color = "gray20", size = 0.2) +
      coord_sf(xlim = c(-65, -40), ylim = c(-65, -58)) +
      labs(title = paste("Date:", current_date),
           x = "Longitude", y = "Latitude") +
      theme_bw() +
      geom_point(data = gv_date, aes(x = lon, y = lat, colour = ship_id, shape = ship_id)) +
      geom_point(data = shirreff, aes(lon, lat), shape = 18, size = 5, colour = "blue3") +
      geom_point(data = cierva, aes(lon, lat), shape = 18, size = 5, colour = "blue3") +
      geom_point(data = kopaitic, aes(lon, lat), shape = 18, size = 5, colour = "blue3")
    
    # Save frame as PNG
    frame_file <- file.path(temp_dir, sprintf("frame_%04d.png", i))
    ggsave(frame_file, p, width = 20, height = 15, units = "cm", dpi = 100)
    frame_files <- c(frame_files, frame_file)
    
    setTxtProgressBar(pb, i)
  }
  
  close(pb)
  
  # Combine frames into GIF 
  library(magick)
  img_list <- lapply(frame_files, image_read)
  animation <- image_animate(image_join(img_list), fps = 5)
  image_write(animation, output_file)
  
  # Clean up
  unlink(frame_files)
  
  return(output_file)
}

# Run the animation 
animation_file <- create_animation_sequential(
  sic = sic,
  gv = gv,
  antarctica = antarctica,
  subs=subs,
  shirreff = shirreff,
  cierva = cierva,
  kopaitic = kopaitic,
  output_file = "sea_ice_animation.gif"
)
