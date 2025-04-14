

library(tidyverse)
library(raster)
library(sf)
library(magick)
library(here)
library(tictoc)
library(randomForest)

source_directory <- "D:/SERDP/lst_gv_models/santa_clara/"

loadRaster <- function(filename, var_name, mask)
{
  print(paste("Loading file", filename))
  # Get year out of filename
  year <- as.numeric(str_extract(basename(filename), "[0-9][0-9][0-9][0-9]"))
  # Load data from raster
  new_data <- terra::rast(filename)
  # Extract values within mask to a dataframe
  new_df <- as.data.frame(terra::extract(new_data, mask, xy=TRUE))
  # Convert to long format (each band is a week)
  new_df_long <- new_df %>%
    pivot_longer(2:(ncol(new_df)-2), names_to="period_str", values_to=var_name) %>%
    mutate(period = as.numeric(substr(period_str, 1, nchar(period_str)-20)),
           year = year) %>%
    dplyr::select(-period_str)
  return(new_df_long)
}

# Mask to riparian area
riparian_mask <- st_read(paste(source_directory, "extent/scr_mainstem_v2.kml", sep="")) %>%
  st_transform("EPSG:32611")
# Mask to clip off some non-riparian bits (agriculture especially)
extra_clipping_mask <- st_read(paste(source_directory, "extent/scr_selective_removal.shp", sep="")) %>%
  st_transform("EPSG:32611")

gv_files <- list.files(paste(source_directory, "GV/", sep=""), 
                       pattern = "GV_phenology_Santa_Clara_[0-9][0-9][0-9][0-9].tif$",
                       full.names = TRUE)
# remove the data from 2024 because air temperature data weren't available
gv_files <- gv_files[-length(gv_files)]
lst_files <- list.files(paste(source_directory, "LST/", sep=""), 
                       pattern = "santa_clara_lst_rel_[0-9][0-9][0-9][0-9]_phenoseries.tif$",
                       full.names = TRUE)

# SPEI data
spei <- merge(merge(merge(merge(read_csv(here::here("spei_01.csv")),
                                read_csv(here::here("spei_03.csv")), 
                                by = "dates"),
                          read_csv(here::here("spei_06.csv")), 
                          by = "dates"),
                    read_csv(here::here("spei_12.csv")), 
                    by = "dates"),
              read_csv(here::here("spei_24.csv")), 
              by = "dates") %>%
  mutate(date = as.Date(dates)) %>%
  mutate(year = lubridate::year(date),
         month = lubridate::month(date),
         day = lubridate::day(date)) %>%
  dplyr::select(-dates) %>%
  mutate(period = month*2-1)
# SPEI has half the period we're using for the satellite imagery - fill in the gaps
spei <- rbind(spei, spei %>% mutate(period = period-1))
spei_annual_summer <- spei %>% 
  filter(month %in% c(5,6,7,8,9)) %>%
  group_by(year) %>%
  summarize(spei01 = mean(spei01),
            spei03 = mean(spei03),
            spei06 = mean(spei06),
            spei12 = mean(spei12),
            spei24 = mean(spei24))
spei_plot <- ggplot(spei) + 
  geom_line(aes(x=year+month/12, y=spei01), col="red") + 
  geom_line(aes(x=year+month/12, y=spei03), col="orange") + 
  geom_line(aes(x=year+month/12, y=spei06), col="goldenrod") + 
  geom_line(aes(x=year+month/12, y=spei12), col="forestgreen") + 
  geom_line(aes(x=year+month/12, y=spei24), col="blue") + 
  theme_bw() + 
  xlab("Year") + 
  ylab("SPEI") + 
  scale_x_continuous(limits=c(2010,2020)) + 
  ggtitle("SPEI Timeseries")
spei_plot
ggsave(here::here("output_plots/spei_plot.png"),
       spei_plot, width=5, height=3)
# Check for mutual correlation among pairs of SPEI variables
summary(lm(data=spei, spei01 ~ spei03))$adj.r.squared
summary(lm(data=spei, spei01 ~ spei06))$adj.r.squared
summary(lm(data=spei, spei01 ~ spei12))$adj.r.squared
summary(lm(data=spei, spei01 ~ spei24))$adj.r.squared
summary(lm(data=spei, spei03 ~ spei06))$adj.r.squared
summary(lm(data=spei, spei03 ~ spei12))$adj.r.squared
summary(lm(data=spei, spei03 ~ spei24))$adj.r.squared
summary(lm(data=spei, spei06 ~ spei12))$adj.r.squared
summary(lm(data=spei, spei06 ~ spei24))$adj.r.squared
summary(lm(data=spei, spei12 ~ spei24))$adj.r.squared
# From the above, there are strong correlations between 1+3, 3+6, 6+12, and 12+24, and moderate between 1+6 and 3+12
# BUT it's probably safe to retain both 1+12 or 3+24

# Generate some GIF timelapses from the LST and NDVI datasets
#   Prevent GDAL from writing extra .xml files for .png outputs
#terra::setGDALconfig("GDAL_PAM_ENABLED", "FALSE")
# Generate raster
writeRasterToPng <- function(filename_in, min_value, max_value, folder, invert=FALSE){
  # Get year out of filename
  year <- str_extract(basename(filename_in), "[0-9][0-9][0-9][0-9]")
  print(paste("Working on year ", year, " for folder ", folder, sep=""))
  raster_in <- terra::rast(filename_in)
  raster_in <- mask(raster_in, riparian_mask)
  raster_in <- mask(raster_in, extra_clipping_mask, inverse=TRUE)
  # Rescale raster to RGB after selecting target bands
  raster_rescaled <- c(raster_in[[9]], raster_in[[14]], raster_in[[19]])
  raster_rescaled[raster_rescaled < min_value] <- min_value
  raster_rescaled[raster_rescaled > max_value] <- max_value
  raster_rescaled <- (raster_rescaled-min_value)/(max_value-min_value)*1
  if(invert)
    raster_rescaled <- 1 - raster_rescaled
  raster_rescaled[is.na(raster_rescaled)] <- 1
  # Convert to magick raster
  #   note - looks like changing text color (and box background color) isn't supported on my version of R. Annoying 
  #   https://stackoverflow.com/questions/75587728/text-color-not-honored-by-image-annotate-from-the-magick-package-in-r-of-imagema
  magick_img <- magick::image_read(as.array(raster_rescaled))
  magick_img <- image_annotate(magick_img, year, location="+900+450", font="Arial", size=20, col="#FFFFFF")
  # Generate output filename and path
  directory <- paste(dirname(filename_in)) # Get directory for source NDVI data
  file_basename <- basename(filename_in) # Get filename
  file_basename <- substr(file_basename, 1, (nchar(file_basename)-4)) # remove ".tif" from end of filename
  new_filename <- paste(directory, folder, file_basename, ".png", sep="")
  # Write new PNG
  image_write(magick_img, new_filename) # write a PNG 
  rm(magick_img)
  return(new_filename)
}
# Generate LST and gv png files
gv_PNG_files <- lapply(gv_files, writeRasterToPng, min_value=0.1, max_value=0.9, folder="/gv_png/")
lst_PNG_files <- lapply(lst_files, writeRasterToPng, min_value=-5, max_value=15, folder="/lst_png/", invert=TRUE)
# Generate animated gifs
gv_PNG_images <- image_join(lapply(gv_PNG_files[17:40], image_read))
gv_animated <- image_animate(gv_PNG_images, fps = 2)
image_write(image = gv_animated, path = here::here("output_plots/gv_animation.gif"))
lst_PNG_images <- image_join(lapply(lst_PNG_files[17:40], image_read))
lst_animated <- image_animate(lst_PNG_images, fps = 2)
image_write(image = lst_animated, path = here::here("output_plots/lst_animation.gif"))
rm(gv_PNG_files, gv_PNG_images, gv_animated, lst_PNG_files, lst_PNG_images, lst_animated)

# Load all files
all_gv <- lapply(gv_files, FUN=loadRaster, var_name="GV",  mask=riparian_mask) %>%
  bind_rows()
all_lst <- lapply(lst_files, FUN=loadRaster, var_name="LST",  mask=riparian_mask) %>%
  bind_rows()
all_data <- all_gv %>%
  mutate(LST = all_lst$LST) 
rm(all_gv, all_lst)
# Add SPEI data to satellite imagery
all_data <- all_data %>% merge(spei %>% 
                                 dplyr::select(year, period, spei01, spei03, spei06, spei12, spei24), 
                               by=c("year","period"))

# Compare relationships between LST and NDVI across SPEI integration bands
getModel <- function(target_period, formula, model_name, GV_threshold = 0)
{
  print(paste("  Processing the period", target_period))
  # Filter to data within period
  target_data <- all_data %>% 
    filter(period == target_period) %>%
    filter(GV > GV_threshold)
  new_model <- lm(data=target_data, formula)
  model_coef <- summary(new_model)$coefficients
  model_estimates <- model_coef[,1]
  output_names <- paste(model_name, rownames(model_coef), sep="_")
  output_df <- data.frame(value = (model_estimates),
                          variable = rownames(model_coef),
                          model = model_name)
  output_df <- rbind(output_df, 
                     data.frame(value = summary(new_model)$adj.r.squared,
                                variable = "r_sqd",
                                model = model_name)) %>%
    remove_rownames() %>%
    mutate(period = target_period)
  return(list(output_df, new_model))
}
getModelByWeek <- function(formula, model_name, GV_threshold = 0)
{
  print(paste("Processing the model", model_name))
  list_of_model_results <- lapply(unique(all_data$period) %>% sort(), 
                                  getModel, formula=formula, model_name=model_name, GV_threshold=GV_threshold)
  dataframe <- lapply(list_of_model_results,
                            function(new_list){
                              return(new_list[[1]])
                            }) %>% bind_rows()
  model_list <- lapply(list_of_model_results,
                       function(new_list){
                         return(new_list[[2]])
                       })
  return(list(dataframe, model_list))
}

# First, a plot showing the slope of each SPEI product vs. LST in each week of the year, when modeled alongside GV
LST_SPEI_intercomparison <- rbind(getModelByWeek("LST ~ GV + spei01", "GV_spei01")[[1]],
                                  getModelByWeek("LST ~ GV + spei03", "GV_spei03")[[1]],
                                  getModelByWeek("LST ~ GV + spei06", "GV_spei06")[[1]],
                                  getModelByWeek("LST ~ GV + spei12", "GV_spei12")[[1]],
                                  getModelByWeek("LST ~ GV + spei24", "GV_spei24")[[1]])
LST_SPEI_intercomparison_plot_rsqd <- ggplot(LST_SPEI_intercomparison %>% filter(variable == "r_sqd")) + 
  geom_line(aes(x=period/2, y=value, group=model, col=model)) + 
  theme_bw() + 
  geom_hline(yintercept=0, linetype="dashed", col="gray") + 
  xlab("Month") + 
  ylab("Slope of Relative LST vs. SPEI")
scale_x_continuous(limits=c(0,12), expand=c(0,0)) 
LST_SPEI_intercomparison_plot_rsqd
ggsave(here::here("output_plots/LST_SPEI_intercomparison_rsqd.png"),
       LST_SPEI_intercomparison_plot_rsqd, width=5, height=3)

# Next, a model showing how much of the variation in LST is explained by each SPEI product alone (no GV)
LST_SPEI_intercomparison_no_GV <- rbind(getModelByWeek("LST ~ spei01", "spei01")[[1]],
                                        getModelByWeek("LST ~ spei03", "spei03")[[1]],
                                        getModelByWeek("LST ~ spei06", "spei06")[[1]],
                                        getModelByWeek("LST ~ spei12", "spei12")[[1]],
                                        getModelByWeek("LST ~ spei24", "spei24")[[1]])
LST_SPEI_intercomparison_no_GV_plot_rsqd <- ggplot(LST_SPEI_intercomparison_no_GV %>% filter(variable == "r_sqd")) + 
  geom_line(aes(x=period/2, y=value, group=model, col=model)) + 
  theme_bw() + 
  geom_hline(yintercept=0, linetype="dashed", col="gray") + 
  xlab("Month") + 
  ylab("R^2 of Relative LST vs. SPEI")
scale_x_continuous(limits=c(0,12), expand=c(0,0)) 
LST_SPEI_intercomparison_no_GV_plot_rsqd
ggsave(here::here("output_plots/LST_SPEI_intercomparison_no_GV_plot_rsqd.png"),
       LST_SPEI_intercomparison_no_GV_plot_rsqd, width=5, height=3)

# Next, let's build models which combine spei24, spei3, GV, and year to predict LST
# Then we'll test what the effects are if we leave out only one of each of those variables
LST_models <- rbind(getModelByWeek("LST ~ GV + spei03 + spei24 + year", "full_model")[[1]],
                    getModelByWeek("LST ~ GV + spei03 + spei24", "no_year")[[1]],
                    getModelByWeek("LST ~ spei03 + spei24 + year", "no_GV")[[1]],
                    getModelByWeek("LST ~ GV + spei24 + year", "no_spei03")[[1]],
                    getModelByWeek("LST ~ GV + spei03 + year", "no_spei24")[[1]])
LST_model_intercomparison_r_sqd_plot <- ggplot(LST_models %>% filter(variable=="r_sqd")) + 
  geom_line(aes(x=period/2, y=value, col=model, group=model)) + 
  theme_bw() + 
  geom_hline(yintercept=0, linetype="dashed", col="gray") + 
  xlab("Month") + 
  ylab("R^2 of Model to Predict Relative LST") + 
  scale_color_manual(values=c("full_model"="black",
                              "no_GV"="forestgreen",
                              "no_spei03"="cyan",
                              "no_spei24"="blue",
                              "no_year"="orange"))
LST_model_intercomparison_r_sqd_plot
ggsave(here::here("output_plots/LST_model_intercomparison_r_sqd.png"),
       LST_model_intercomparison_r_sqd_plot, width=5, height=3)


# Generate some output plots demonstrating model coefficients
LST_spei_dependence <- ggplot(LST_models %>% filter(model=="full_model", variable %in% c("spei03","spei24"))) + 
  geom_line(aes(x=period/2, y=value, group=variable, col=variable)) + 
  theme_bw() + 
  geom_hline(yintercept=0, linetype="dashed", col="gray") + 
  xlab("Month") + 
  ylab("Slope of Relative LST vs. SPEI (°C)") + 
  theme(legend.position="none")
LST_spei_dependence
ggsave(here::here("output_plots/LST_spei_dependence.png"), LST_spei_dependence, width=5, height=3)
# Generate some output plots demonstrating model coefficients
LST_GV_dependence <- ggplot(LST_models %>% filter(model=="full_model", variable %in% c("GV"))) + 
  geom_line(aes(x=period/2, y=value), col="forestgreen") + 
  theme_bw() + 
  geom_hline(yintercept=0, linetype="dashed", col="gray") + 
  xlab("Month") + 
  ylab("Slope of Relative LST vs. GV (°C)")
LST_GV_dependence
ggsave(here::here("output_plots/LST_GV_dependence.png"), LST_GV_dependence, width=5, height=3)
# Generate some output plots demonstrating model coefficients
LST_year_dependence <- ggplot(LST_models %>% filter(model=="full_model", variable %in% c("year"))) + 
  geom_line(aes(x=period/2, y=value), col="orange") + 
  theme_bw() + 
  geom_hline(yintercept=0, linetype="dashed", col="gray") + 
  xlab("Month") + 
  ylab("Slope of Relative LST vs. Year (°C/year)")
LST_year_dependence
ggsave(here::here("output_plots/LST_year_dependence.png"), LST_year_dependence, width=5, height=3)

#LST_final_model <- lm(data=all_data,
#                      LST ~ GV + spei24 + spei03 + year)
LST_final_model_list <- getModelByWeek("LST ~ GV + spei03 + spei24 + year", "final_model")
LST_final_model_df <- LST_final_model_list[[1]]
LST_final_models <- LST_final_model_list[[2]]
# Set the names equal to each period
names(LST_final_models) <- unique(all_data$period) %>% sort()

# Add predicted and residual LST to data 
applyLSTModel <- function(target_period)
{
  # Filter to just data within the target seasonal period
  data_in_period <- all_data %>% 
    filter(period==target_period)
  # Get LST model for that period
  LST_model <- LST_final_models[[as.character(target_period)]]
  # Apply model to predict LST
  data_in_period$LST_predicted <- predict(LST_model, newdata = data_in_period)
  # Return output
  return(data_in_period)
}
all_data <- lapply(unique(all_data$period) %>% sort(), 
                   applyLSTModel) %>%
  bind_rows()
all_data$LST_predicted <- predict(LST_final_model, newdata=all_data)
all_data$LST_residual <- all_data$LST - all_data$LST_predicted
all_data$pixel_loc <- paste(all_data$x, all_data$y)

write_csv(all_data,
          here::here("all_data.csv"))

# Write output rasters with LST residual estimates
correctRaster <- function(GV_filename, LST_filename, mask)
{
  print(paste("Processing files", GV_filename, "and", LST_filename))
  # Get year out of filename
  year <- as.numeric(str_extract(basename(GV_filename), "[0-9][0-9][0-9][0-9]"))
  output_filepath_prediction <- paste(source_directory, "LST_residual/LST_prediction", year, ".tif", sep="")
  output_filepath_residual <- paste(source_directory, "LST_residual/LST_residual", year, ".tif", sep="")
  # Load data from rasters
  GV_rast <- terra::rast(GV_filename)
  LST_rast <- terra::rast(LST_filename)
  # Period names
  period_names <- unique(all_data$period) %>% sort()
  names(GV_rast) <- paste("period_", period_names, sep="")
  names(LST_rast) <- paste("period_", period_names, sep="")
  correctPeriod <- function(period){
    print(paste("  Operating on period", period))
    # Select the target model
    target_model <- LST_final_models[[as.character(period)]]
    # Subset raster data to target period
    GV_period_rast <- GV_rast[[paste("period_",period,sep="")]]
    LST_period_rast <- LST_rast[[paste("period_",period,sep="")]]
    # Generate single-valued rasters for year and SPEI values
    spei_values <- (spei %>% rbind(spei %>% mutate(period = as.numeric(period)-1))) %>%
      filter(year == !!year, period == as.numeric(!!period))
    spei03_period_rast <- terra::rast(GV_period_rast, vals=spei_values[1,]$spei03)
    spei24_period_rast <- terra::rast(GV_period_rast, vals=spei_values[1,]$spei24)
    year_rast <- terra::rast(GV_period_rast, vals=year)
    # Combine into one raster
    input_raster <- c(GV_period_rast, LST_period_rast, spei03_period_rast, spei24_period_rast, year_rast)
    names(input_raster) <- c("GV", "LST", "spei03", "spei24", "year")
    # Predict and return a raster
    return(terra::predict(model=target_model, object=input_raster))
  }
  
  new_raster <- lapply(period_names, 
                       correctPeriod) %>% 
    terra::rast() %>%
    terra::mask(riparian_mask)
  terra::writeRaster(new_raster, output_filepath_prediction, overwrite=TRUE)
  terra::writeRaster(LST_rast - new_raster, output_filepath_residual, overwrite=TRUE)
  return(LST_rast - new_raster)
}
lapply(1:length(gv_files), 
       function(ind){
         return(correctRaster(gv_files[[ind]], lst_files[[ind]], 
                              st_union(riparian_mask[,1], extra_clipping_mask[,1])[,1]))
       })



# Rewrite to a matrix format
tic()
GV_matrix <- all_data %>% 
  mutate(time = paste(year, period, sep="_")) %>%
  dplyr::select(pixel_loc, time, GV) %>%
  pivot_wider(names_from = time,
              values_from = GV) %>% 
  as.matrix()
write.table(GV_matrix, here::here("GV_matrix.csv"))
toc()
tic()
LST_matrix <- all_data %>% 
  mutate(time = paste(year, period, sep="_")) %>%
  dplyr::select(pixel_loc, time, LST) %>%
  pivot_wider(names_from = time,
              values_from = LST) %>% 
  as.matrix()
write.table(LST_matrix, here::here("LST_matrix.csv"))
toc()
tic()
LST_residual_matrix <- all_data %>% 
  mutate(time = paste(year, period, sep="_")) %>%
  dplyr::select(pixel_loc, time, LST_residual) %>%
  pivot_wider(names_from = time,
              values_from = LST_residual) %>% 
  as.matrix()
write.table(LST_residual_matrix, here::here("LST_residual_matrix.csv"))
toc()

# Simple metric for distance upstream - distance from origin in Northing and Easting
distances_upstream <- lapply(strsplit(GV_matrix[,1], " "),
                             function(str_tuple){
                               return(sqrt(as.numeric(str_tuple[[1]])^2 + as.numeric(str_tuple[[2]])^2))
                             }) %>%
  unlist()

# To just re-read in the existing, pre-saved matrices: 
GV_matrix <- read.table(here::here("GV_matrix.csv"))
LST_matrix <- read.table(here::here("LST_matrix.csv"))
LST_residual_matrix <- read.table(here::here("LST_residual_matrix.csv"))



# Run some simple regression tests predicting greenness based on past GV and LST data, plus SPEI
# Build test dataset from 2016
set.seed(1)
example_df <- as.data.frame((data.frame(GV = GV_matrix[,"X2016_10"],
                                        GV_change_early_this = GV_matrix[,"X2016_16"] - GV_matrix[,"X2016_10"],
                                        GV_change_late_last = GV_matrix[,"X2016_16"] - GV_matrix[,"X2015_16"],
                                        GV_final = GV_matrix[,"X2016_16"],
                                        GV_m2 = GV_matrix[,"X2016_8"],
                                        GV_m3 = GV_matrix[,"X2016_6"],
                                        GV_m18 = GV_matrix[,"X2015_16"],
                                        GV_m24 = GV_matrix[,"X2015_10"],
                                        LST_residual = LST_residual_matrix[,"X2016_10"],
                                        LST_residual_m24 = LST_residual_matrix[,"X2015_10"],
                                        LST = LST_matrix[,"X2016_10"],
                                        distance = distances_upstream,
                                        distance_sqd = distances_upstream^2))) %>% 
  # filter(GV >= 0.4, GV_m18 >= 0.4) %>%
  drop_na() %>%
  mutate(ind = 1:n())


# Simple function to get root mean square error between 
#   actual - numeric vector of true values
#   predicted - numeric vector of predicted values
rmse <- function(actual, predicted)
{
  squared_error <- (predicted-actual)^2
  
  return(sqrt(mean(squared_error, na.rm=TRUE)))
}


# Write a function to build and test a model
buildNewModel <- function(target_formula, target_variable, dataset)
{
  # Divide dataset randomly in half for test and training
  training_df <- dataset[sample(1:nrow(dataset), nrow(dataset)/2),]
  validation_df <- dataset[-training_df$ind,]
  
  # Predicted variable name
  predicted_variable <- paste(target_variable, "_predicted", sep="")
  # Add new column to training df for prediction variable
  new_data <- data.frame(var = rep(0, nrow(training_df)))
  names(new_data) <- predicted_variable
  training_df <- cbind(training_df, new_data)
  # Add new column to validation df for prediction variable
  new_data <- data.frame(var = rep(0, nrow(validation_df)))
  names(new_data) <- predicted_variable
  validation_df <- cbind(validation_df, new_data)
  
  # Build the model
  new_model <- randomForest(as.formula(target_formula), data=training_df, ntree=1000, mtry=5)
  # Predict results for training and validation datasets
  training_df[,predicted_variable] <- predict(new_model, newdata=training_df)
  validation_df[,predicted_variable] <- predict(new_model, newdata=validation_df)
  # Explain the model
  #   Variable importance
  print(importance(new_model))
  #   Visualize distribution of values
  validation_plot <- ggplot(validation_df) + 
    geom_density_2d_filled(aes(x=!!ensym(predicted_variable), y=!!ensym(target_variable))) + 
    scale_x_continuous(limits=c(quantile(validation_df[,predicted_variable],0.05),quantile(validation_df[,predicted_variable],0.95))) + 
    scale_y_continuous(limits=c(quantile(validation_df[,target_variable],0.05),quantile(validation_df[,target_variable],0.95)))
  print(validation_plot)
  # Summarize model residuals
  #summary(lm(data=validation_df %>% mutate(residual = !!as.name(target_variable) - !!as.name(predicted_variable)),
  #           paste(residual, " ~ ", predicted_variable, sep="")))
  print(paste("RMSE Training: ", rmse(training_df[,predicted_variable], training_df[,target_variable])), sep="")
  print(paste("RMSE Validation: ", rmse(validation_df[,predicted_variable], validation_df[,target_variable])), sep="")
  print(paste("R^2 Training: ", summary(lm(data=training_df, paste(predicted_variable, " ~ ", target_variable, sep="")))$adj.r.squared), sep="")
  print(paste("R^2 Validation: ", summary(lm(data=validation_df, paste(predicted_variable, " ~ ", target_variable, sep="")))$adj.r.squared), sep="")
  print(paste("Slope Validation: ", summary(lm(data=validation_df, paste(predicted_variable, " ~ ", target_variable, sep="")))$coefficients[2,1]), sep="")
  
  return(list(new_model, training_df, validation_df))
}

explainModel <- function(model_training_list, predicted_variable, target_variable)
{
  #   Variable importance
  print(importance(model_training_list[[1]])/sum(importance(model_training_list[[1]]))*100)
  #   Visualize distribution of values
  prediction_plot <- ggplot(model_training_list[[3]]) + 
    geom_density_2d_filled(aes(x=!!ensym(predicted_variable), y=!!ensym(target_variable))) + 
    scale_x_continuous(limits=c(quantile(model_training_list[[3]][,predicted_variable],0.1),quantile(model_training_list[[3]][,predicted_variable],0.9))) + 
    scale_y_continuous(limits=c(quantile(model_training_list[[3]][,target_variable],0.1),quantile(model_training_list[[3]][,target_variable],0.9))) + 
    theme_bw() + 
    xlab("GV Predicted") + 
    ylab("GV Actual")
  print('survived')
  print(prediction_plot)
  # Summarize model residuals
  #summary(lm(data=model_training_list %>% mutate(residual = !!as.name(target_variable) - !!as.name(predicted_variable)),
  #           paste(residual, " ~ ", predicted_variable, sep="")))
  print(rmse(model_training_list[[2]][,predicted_variable], model_training_list[[2]][,target_variable]))
  print(rmse(model_training_list[[3]][,predicted_variable], model_training_list[[3]][,target_variable]))
  print(summary(lm(data=model_training_list[[2]], paste(predicted_variable, " ~ ", target_variable, sep="")))$adj.r.squared)
  print(summary(lm(data=model_training_list[[3]], paste(predicted_variable, " ~ ", target_variable, sep="")))$adj.r.squared)
  print(summary(lm(data=model_training_list[[3]], paste(predicted_variable, " ~ ", target_variable, sep="")))$coefficients[2,1])
}


# Model change relative to last year
etm_model_GV_early_this <- buildNewModel("GV_change_early_this ~ GV+GV_m24+GV_m18+LST_residual+LST_residual_m24+LST+distance+distance_sqd",
                                     dataset = example_df, target_variable="GV_change_early_this")
# Model change relative to start of this year
etm_model_GV_late_last <- buildNewModel("GV_change_late_last ~ GV+GV_m24+GV_m18+LST_residual+LST_residual_m24+LST+distance+distance_sqd",
                                     dataset = example_df, target_variable="GV_change_late_last")
# Model overall fraction at end of year
etm_model_GV_final <- buildNewModel("GV_final ~ GV+GV_m24+GV_m18+LST_residual+LST_residual_m24+LST+distance+distance_sqd",
                                     dataset = example_df, target_variable="GV_final")









# Visualize Mortality in Drought
all_data_2016 <- all_data %>% filter(year==2016)
woodland_sites_2016 <- (all_data_2016 %>% 
                          filter(period == 10, GV > 0.5))$pixel_loc %>% sort()
woodland_data_2016 <- all_data_2016 %>% filter(pixel_loc %in% woodland_sites_2016)
woodland_declined_2016 <- woodland_data_2016


tic()
# State Model
temp_storage_file <- here::here("temporary_data_with_shifts_incomplete.csv") # backup in-progress work 
# Add data with +1 month lag
all_data_with_shifts <- all_data %>% arrange(year, period, pixel_loc)
back_shift_p2 <- all_data %>%
  dplyr::select("year","period","pixel_loc","GV","LST","LST_residual","spei03","spei24") %>% 
  arrange(year, period, pixel_loc)
names(back_shift_p2) <- c("year","period","pixel_loc","GV_p2","LST_p2","LST_residual_p2","spei03_p2","spei24_p2")
back_shift_p2$period <- back_shift_p2$period - 2
all_data_with_shifts <- merge(all_data_with_shifts,
                              back_shift_p2, 
                              by=c("year", "period", "pixel_loc"))
write_csv(all_data_with_shifts, temp_storage_file)
toc()
tic()
# Add data with +12 month lag
back_shift_p12 <- all_data %>%
  dplyr::select("year","period","pixel_loc","GV","LST","LST_residual","spei03","spei24") %>% 
  arrange(year, period, pixel_loc)
names(back_shift_p12) <- c("year","period","pixel_loc","GV_p12","LST_p12","LST_residual_p12","spei03_p12","spei24_p12")
back_shift_p12$year <- back_shift_p12$year - 1
all_data_with_shifts <- merge(all_data_with_shifts,
                              back_shift_p12, 
                              by=c("year", "period", "pixel_loc"))
write_csv(all_data_with_shifts, temp_storage_file)
toc()
tic()
# Add data with -1 month lag
back_shift_m2 <- all_data %>%
  dplyr::select("year","period","pixel_loc","GV","LST","LST_residual","spei03","spei24") %>% 
  arrange(year, period, pixel_loc)
names(back_shift_m2) <- c("year","period","pixel_loc","GV_m2","LST_m2","LST_residual_m2","spei03_m2","spei24_m2")
back_shift_m2$period <- back_shift_m2$period + 2
all_data_with_shifts <- merge(all_data_with_shifts,
                              back_shift_m2, 
                              by=c("year", "period", "pixel_loc"))
write_csv(all_data_with_shifts, temp_storage_file)
toc()
tic()
# Add data with -2 month lag
back_shift_m4 <- all_data %>%
  dplyr::select("year","period","pixel_loc","GV","LST","LST_residual","spei03","spei24") %>% 
  arrange(year, period, pixel_loc)
names(back_shift_m4) <- c("year","period","pixel_loc","GV_m4","LST_m4","LST_residual_m4","spei03_m4","spei24_m4")
back_shift_m4$period <- back_shift_m4$period + 4
all_data_with_shifts <- merge(all_data_with_shifts,
                              back_shift_m4, 
                              by=c("year", "period", "pixel_loc"))
write_csv(all_data_with_shifts, temp_storage_file)
toc()
tic()
# Add data with -3 month lag
back_shift_m6 <- all_data %>%
  dplyr::select("year","period","pixel_loc","GV","LST","LST_residual","spei03","spei24") %>% 
  arrange(year, period, pixel_loc)
names(back_shift_m6) <- c("year","period","pixel_loc","GV_m6","LST_m6","LST_residual_m6","spei03_m6","spei24_m6")
back_shift_m6$period <- back_shift_m6$period + 6
all_data_with_shifts <- merge(all_data_with_shifts,
                              back_shift_m6, 
                              by=c("year", "period", "pixel_loc"))
write_csv(all_data_with_shifts, temp_storage_file)
toc()
tic()
# Add data with -12 month lag
back_shift_m12 <- all_data %>%
  dplyr::select("year","period","pixel_loc","GV","LST","LST_residual","spei03","spei24") %>% 
  arrange(year, period, pixel_loc)
names(back_shift_m12) <- c("year","period","pixel_loc","GV_m12","LST_m12","LST_residual_m12","spei03_m12","spei24_m12")
back_shift_m12$year <- back_shift_m12$year + 1
all_data_with_shifts <- merge(all_data_with_shifts,
                              back_shift_m12, 
                              by=c("year", "period", "pixel_loc"))
write_csv(all_data_with_shifts, temp_storage_file)
toc()
tic()
# Add data with -24 month lag
back_shift_m24 <- all_data %>%
  dplyr::select("year","period","pixel_loc","GV","LST","LST_residual","spei03","spei24") %>% 
  arrange(year, period, pixel_loc)
names(back_shift_m24) <- c("year","period","pixel_loc","GV_m24","LST_m24","LST_residual_m24","spei03_m24","spei24_m24")
back_shift_m24$year <- back_shift_m24$year + 2
all_data_with_shifts <- merge(all_data_with_shifts,
                              back_shift_m24, 
                              by=c("year", "period", "pixel_loc"))
write_csv(all_data_with_shifts, here::here("all_data_with_shifts_complete.csv"))
toc()
if(file.exists(temp_storage_file))
{
  file.remove(temp_storage_file)
}


growing_season_start <- 8
growing_season_t2_start <- 12
growing_season_t3_start <- 16
growing_season_end <- 20
all_data %>% 
  group_by(year, pixel_loc) %>%
  summarize(GV_mean_growing = mean(GV*(period > growing_season_start)*(period <= growing_season_end), na.rm=TRUE),
            GV_mean_T1 = mean(GV*(period > growing_season_start)*(period <= growing_season_t2_start), na.rm=TRUE),
            GV_mean_T2 = mean(GV*(period > growing_season_t2_start)*(period <= growing_season_t3_start), na.rm=TRUE),
            GV_mean_T3 = mean(GV*(period > growing_season_t3_start)*(period <= growing_season_end), na.rm=TRUE))



# Generate a model to predict GV transitions
generateStateModel <- function(formula, model_name){
  return(lapply((min(all_data_with_shifts$period) : max(all_data_with_shifts$period)), function(target_period){
  print(paste("Working on week ", target_period))
  new_model <- summary(lm(data=as.data.frame(scale(all_data_with_shifts %>%
                                                     filter(period==target_period,
                                                            GV_m12 > 0.6) %>%
                                                     mutate(GV_change_monthly = GV-GV_m2) %>%
                                                     dplyr::select(-pixel_loc))),
                          formula))
  print(new_model)
  coef <- new_model$coefficients
  effects <- t(coef[,1])
  effect_names <- colnames(effects)
  
  model_df <- data.frame(effects)
  names(model_df) <- effect_names
  model_df$r_sqd <- new_model$adj.r.squared
  model_df$period <- target_period
  model_df$model_name <- model_name
  
  return(model_df)}) %>%
  bind_rows())
}
GV_state_model_no_lst <- generateStateModel("GV_p2 ~ spei03_m2 + spei24_m2 + year + GV_m2",
                                            "GV_model_no_LST")
# Visualize model effects
GV_model_plot <- ggplot(LST_models %>% filter(model=="full_model", variable %in% c("year"))) + 
  geom_line(aes(x=period/2, y=value), col="orange") + 
  theme_bw() + 
  geom_hline(yintercept=0, linetype="dashed", col="gray") + 
  xlab("Month") + 
  ylab("Slope of Relative LST vs. Year (°C/year)")
LST_year_dependence
ggsave(here::here("output_plots/LST_year_dependence.png"), LST_year_dependence, width=5, height=3)

GV_state_model <- generateStateModel("GV_p2 ~ spei03_m2 + spei24_m2 + year + GV_m2 + LST_residual_m2",
                                            "GV_model")
# Visualize model effects
GV_model_plot <- ggplot(GV_state_model %>% filter(model=="full_model", variable %in% c("year"))) + 
  geom_line(aes(x=period/2, y=value), col="orange") + 
  theme_bw() + 
  geom_hline(yintercept=0, linetype="dashed", col="gray") + 
  xlab("Month") + 
  ylab("Slope of Relative LST vs. Year (°C/year)")
LST_year_dependence
ggsave(here::here("output_plots/LST_year_dependence.png"), LST_year_dependence, width=5, height=3)


































# Get list of pixels which are initially highly vegetated
# Filter to upper 70th percentile of summer NDVI in the pre-drought period
dense_vegetation_pixels <- (all_data %>%
  filter(year < 2012, week > 22, week < 32) %>%
  mutate(pixel_loc = paste(x, y)) %>%
  group_by(pixel_loc) %>%
  summarize(NDVI_mean = mean(NDVI)) %>%
  filter(NDVI_mean > quantile(NDVI_mean, 0.7)))$pixel_loc

# Add model to full dataset
all_data_with_model <- merge(all_data, before_drought_weekly_patterns, by=c("week")) %>%
  mutate(LST_predicted = LST_NDVI_intercept + LST_NDVI_slope*NDVI) %>%
  mutate(LST_residual = LST - LST_predicted) %>%
  mutate(pixel_loc = paste(x, y))
all_data_with_model <- merge(all_data_with_model, before_drought_spatial_stats, by=c("week","pixel_loc"))


# Create a dataset including values from previous week and year
current_year <- all_data_with_model %>% filter(year > 2005) 
previous_year <- all_data_with_model %>% filter(year < 2022) %>% dplyr::select(pixel_loc, year, week, NDVI, LST, LST_residual)
names(previous_year)[4:6] <- paste(names(previous_year)[4:6], "_previous_year", sep="")
previous_year$year <- previous_year$year + 1
all_data_with_previous_year <- merge(current_year, previous_year, by=c("year","week","pixel_loc"))

previous_month <- all_data_with_model %>% filter(week < 48) %>% dplyr::select(pixel_loc, year, week, NDVI, NDVI_mean, NDVI_std, LST, LST_residual)
names(previous_month)[4:8] <- paste(names(previous_month)[4:8], "_previous_month", sep="")
previous_month$week <- previous_month$week + 4
all_data_with_previous_month_and_year <- merge(all_data_with_previous_year, previous_month, by=c("year","week","pixel_loc"))

next_month <- all_data_with_model %>% filter(week > 4) %>% dplyr::select(pixel_loc, year, week, NDVI, NDVI_mean, NDVI_std, LST, LST_residual)
names(next_month)[4:8] <- paste(names(next_month)[4:8], "_previous_month", sep="")
next_month$week <- next_month$week - 4
all_data_with_previous_month_and_year <- merge(all_data_with_previous_year, next_month, by=c("year","week","pixel_loc"))

next_month <- all_data_with_model %>% filter(week > 4) %>% dplyr::select(pixel_loc, year, week, NDVI, NDVI_mean, NDVI_std, LST, LST_residual)
names(next_month)[4:8] <- paste(names(next_month)[4:8], "_next_month", sep="")
next_month$week <- next_month$week - 4
all_data_time_shifted <- merge(all_data_with_previous_month_and_year, next_month, by=c("year","week","pixel_loc"))

# Summarize interannual change across timeseries
annual_change_summary <- all_data_with_previous_month_and_year %>% 
  group_by(year) %>% 
  filter(week > 22, week < 32) %>%
  summarize(woodland_area = sum(NDVI > woodland_NDVI_threshold),
            LST_NDVI_annualized_slope = summary(lm(LST ~ NDVI))$coefficients[2,1],
            LST_NDVI_annualized_intercept = summary(lm(LST ~ NDVI))$coefficients[1,1],
            NDVI_change = mean(NDVI - NDVI_previous_year),
            NDVI_in_season_change = mean(NDVI - NDVI_previous_month),
            NDVI_residual = mean((NDVI - NDVI_mean)/NDVI_std), 
            NDVI = mean(NDVI),
            LST_change = mean(LST - LST_previous_year),
            LST_in_season_change = mean(LST - LST_previous_month),
            LST_residual_change = mean(LST_residual - LST_residual_previous_year),
            LST_residual_in_season_change = mean(LST_residual - LST_residual_previous_month),
            LST = mean(LST),
            LST_residual = mean(LST_residual)) %>%
  merge(spei_annual_summer %>% filter(year >= 2005, year <= 2022), by="year")
annual_change_summary_woodland <- all_data_with_previous_month_and_year %>% 
  group_by(year) %>% 
  filter(week > 22, week < 32) %>%
  filter(NDVI_previous_year > woodland_NDVI_threshold) %>%
  summarize(LST_NDVI_annualized_slope = summary(lm(LST ~ NDVI))$coefficients[2,1],
            LST_NDVI_annualized_intercept = summary(lm(LST ~ NDVI))$coefficients[1,1],
            NDVI_change = mean(NDVI - NDVI_previous_year),
            NDVI_in_season_change = mean(NDVI - NDVI_previous_month),
            NDVI_residual = mean((NDVI - NDVI_mean)/NDVI_std), 
            NDVI = mean(NDVI),
            LST_change = mean(LST - LST_previous_year),
            LST_in_season_change = mean(LST - LST_previous_month),
            LST_residual_change = mean(LST_residual - LST_residual_previous_year),
            LST_residual_in_season_change = mean(LST_residual - LST_residual_previous_month),
            LST = mean(LST),
            LST_residual = mean(LST_residual)) %>%
  merge(spei_annual_summer %>% filter(year >= 2005, year <= 2022), by="year")
ggplot(annual_change_summary) + 
  geom_line(aes(x=year, y=NDVI))

annual_NDVI_plot <- ggplot(annual_change_summary) + 
  geom_line(aes(x=year, y=NDVI), col="forestgreen") + 
  geom_hline(yintercept=0, linetype="dashed", col="gray") + 
  xlab("Year") + 
  ylab("Average Monthly NDVI") + 
  ggtitle("Annual Average NDVI") + 
  theme_bw() 
annual_NDVI_plot
ggsave("C:/Users/grad/Documents/lst_predicting_NDVI/annual_NDVI_plot.png", annual_NDVI_plot, width=5, height=3)

annual_NDVI_woodland_plot <- ggplot(annual_change_summary_woodland) + 
  geom_line(aes(x=year, y=NDVI), col="forestgreen") + 
  geom_hline(yintercept=0, linetype="dashed", col="gray") + 
  xlab("Year") + 
  ylab("Average Monthly NDVI") + 
  ggtitle("Annual Average NDVI in Woodland") + 
  theme_bw() 
annual_NDVI_woodland_plot
ggsave("C:/Users/grad/Documents/lst_predicting_NDVI/annual_NDVI_woodland_plot.png", annual_NDVI_woodland_plot, width=5, height=3)

annual_NDVI_residual_plot <- ggplot(annual_change_summary_woodland) + 
  geom_line(aes(x=year, y=NDVI_residual), col="forestgreen") + 
  geom_hline(yintercept=0, linetype="dashed", col="gray") + 
  xlab("Year") + 
  ylab("Average Monthly NDVI") + 
  ggtitle("Annual Average NDVI Residual in Woodland") + 
  theme_bw() 
annual_NDVI_residual_plot
ggsave("C:/Users/grad/Documents/lst_predicting_NDVI/annual_NDVI_residual_plot.png", annual_NDVI_residual_plot, width=5, height=3)

annual_NDVI_interannual_change_plot <- ggplot(annual_change_summary_woodland) + 
  geom_line(aes(x=year, y=NDVI_change), col="forestgreen") + 
  geom_hline(yintercept=0, linetype="dashed", col="gray") + 
  xlab("Year") + 
  ylab("Average NDVI Change") + 
  ggtitle("Average Interannual NDVI Change in Woodland") + 
  theme_bw() 
annual_NDVI_interannual_change_plot
ggsave("C:/Users/grad/Documents/lst_predicting_NDVI/annual_NDVI_interannual_change_plot.png", annual_NDVI_interannual_change_plot, width=5, height=3)

annual_NDVI_in_season_change_plot <- ggplot(annual_change_summary_woodland) + 
  geom_line(aes(x=year, y=NDVI_in_season_change), col="forestgreen") + 
  geom_hline(yintercept=0, linetype="dashed", col="gray") + 
  xlab("Year") + 
  ylab("Average NDVI Change") + 
  ggtitle("Average Within-Season NDVI Change in Woodland") + 
  theme_bw() 
annual_NDVI_in_season_change_plot
ggsave("C:/Users/grad/Documents/lst_predicting_NDVI/annual_NDVI_in_season_change_plot.png", annual_NDVI_in_season_change_plot, width=5, height=3)


annual_LST_plot <- ggplot(annual_change_summary) + 
  geom_line(aes(x=year, y=LST), col="blue") + 
  geom_hline(yintercept=0, linetype="dashed", col="gray") + 
  xlab("Year") + 
  ylab("Average Monthly LST") + 
  ggtitle("Annual Average LST") + 
  theme_bw() 
annual_LST_plot
ggsave("C:/Users/grad/Documents/lst_predicting_NDVI/annual_LST_plot.png", annual_LST_plot, width=5, height=3)

annual_LST_woodland_plot <- ggplot(annual_change_summary_woodland) + 
  geom_line(aes(x=year, y=LST), col="blue") + 
  geom_hline(yintercept=0, linetype="dashed", col="gray") + 
  xlab("Year") + 
  ylab("Average Monthly LST") + 
  ggtitle("Annual Average LST in Woodland") + 
  theme_bw() 
annual_LST_woodland_plot
ggsave("C:/Users/grad/Documents/lst_predicting_NDVI/annual_LST_woodland_plot.png", annual_LST_woodland_plot, width=5, height=3)

annual_LST_residual_plot <- ggplot(annual_change_summary_woodland) + 
  geom_line(aes(x=year, y=LST_residual), col="blue") + 
  geom_hline(yintercept=0, linetype="dashed", col="gray") + 
  xlab("Year") + 
  ylab("Average Monthly LST Residual") + 
  ggtitle("Annual Average LST Residual in Woodland") + 
  theme_bw() 
annual_LST_residual_plot
ggsave("C:/Users/grad/Documents/lst_predicting_NDVI/annual_LST_residual_plot.png", annual_LST_residual_plot, width=5, height=3)

annual_LST_interannual_change_plot <- ggplot(annual_change_summary_woodland) + 
  geom_line(aes(x=year, y=LST_change), col="blue") + 
  geom_hline(yintercept=0, linetype="dashed", col="gray") + 
  xlab("Year") + 
  ylab("Average Monthly LST Change") + 
  ggtitle("Average Interannual LST Change in Woodland") + 
  theme_bw() 
annual_LST_interannual_change_plot
ggsave("C:/Users/grad/Documents/lst_predicting_NDVI/annual_LST_interannual_change_plot.png", annual_LST_interannual_change_plot, width=5, height=3)

annual_LST_in_season_change_plot <- ggplot(annual_change_summary_woodland) + 
  geom_line(aes(x=year, y=LST_in_season_change), col="blue") + 
  geom_hline(yintercept=0, linetype="dashed", col="gray") + 
  xlab("Year") + 
  ylab("Average Monthly LST Change") + 
  ggtitle("Average Within-Season LST Change in Woodland") + 
  theme_bw() 
annual_LST_in_season_change_plot
ggsave("C:/Users/grad/Documents/lst_predicting_NDVI/annual_LST_in_season_change_plot.png", annual_LST_in_season_change_plot, width=5, height=3)


# How many sites were reduced in extent or increased in extent vs. the previous year? 
woodland_NDVI_threshold <- 0.6

sites_reduced <- all_data_with_previous_month_and_year %>%
  filter(week > 22, week < 32,
         NDVI_previous_year > woodland_NDVI_threshold) %>% 
  filter((NDVI_previous_year - NDVI) > 1.96*NDVI_std) %>%
  mutate(change = "Mortality")

sites_increased <- all_data_with_previous_month_and_year %>%
  filter(week > 22, week < 32,
         NDVI_previous_year > woodland_NDVI_threshold) %>% 
  filter((NDVI_previous_year - NDVI) < -1.96*NDVI_std)  %>%
  mutate(change = "Growth")

sites_steady <- all_data_with_previous_month_and_year %>%
  filter(week > 22, week < 32,
         NDVI_previous_year > woodland_NDVI_threshold) %>% 
  filter((NDVI_previous_year - NDVI) < 1.96*NDVI_std,
         (NDVI_previous_year - NDVI) > -1.96*NDVI_std) %>%
  mutate(change = "Static")

# How many sites reduced in extent not just since last year, but since last month?
sites_reduced_suddenly <- sites_reduced %>% 
  filter(((NDVI_previous_month-NDVI_mean_previous_month)-(NDVI-NDVI_mean)) > 1.96*((NDVI_std+NDVI_std_previous_month)/2),
         NDVI_previous_year > woodland_NDVI_threshold) %>%
  mutate(change = "Sudden Mortality")

ggplot() + 
  geom_line(data=sites_reduced %>% group_by(year) %>% tally(), 
            aes(x=year, y=n), col="red") + 
  geom_line(data=sites_increased %>% group_by(year) %>% tally(), 
            aes(x=year, y=n), col="blue")

# Are spots which reduced suddenly, which reduced, and which grew different in LAST MONTH's residual LST? 
ggplot(rbind(sites_reduced, sites_reduced_suddenly, sites_increased) %>% filter(year>2010, year<2019)) + 
  geom_density(aes(x=LST_residual_previous_month, group=change, col=change)) 
t.test(sites_increased$LST_residual_previous_month, sites_reduced_suddenly$LST_residual_previous_month)
t.test(sites_reduced$LST_residual_previous_month, sites_reduced_suddenly$LST_residual_previous_month)

# Are spots which reduced suddenly, which reduced, and which grew different in LAST MONTH's overall relative LST? 
ggplot(rbind(sites_reduced, sites_reduced_suddenly, sites_increased) %>% filter(year>2010, year<2019)) + 
  geom_density(aes(x=LST_previous_month, group=change, col=change)) 
t.test(sites_increased$LST_previous_month, sites_reduced_suddenly$LST_previous_month)
t.test(sites_reduced$LST_previous_month, sites_reduced_suddenly$LST_previous_month)









# Some tests looking at mortality specifically in 2016

sudden_mortality_2016 <- all_data_with_previous_month_and_year %>%
  filter(week > 21, week < 40,
         NDVI_previous_year > woodland_NDVI_threshold) %>% 
  filter((NDVI_previous_year - NDVI) > 1.96*NDVI_std) %>%
  mutate(change = "Mortality") %>%
  filter(((NDVI_previous_month-NDVI_mean_previous_month)-(NDVI-NDVI_mean)) > ((NDVI_std+NDVI_std_previous_month)/2),
         NDVI_previous_year > woodland_NDVI_threshold) %>%
  mutate(change = "Sudden Mortality") %>%
  filter(year==2016) %>%
  group_by(pixel_loc) %>%
  arrange(week) %>%
  filter(row_number()==1,
         week > 22)
sudden_mortality_2016_rast <- sudden_mortality_2016 %>% 
  ungroup() %>% 
  dplyr::select(x,y,week,NDVI,NDVI_previous_month) %>%
  terra::rast()
plot(sudden_mortality_2016_rast[[1]])

sudden_mortality_2016_all_data <- sudden_mortality_2016 %>%
  mutate(mortality_week = week,
         mortality_NDVI = NDVI,
         mortality_NDVI_interannual_decrease = NDVI_previous_year - NDVI,
         mortality_NDVI_seasonal_decrease = NDVI_previous_month - NDVI,
         mortality_LST = LST,
         mortality_LST_residual = LST_residual) %>%
  ungroup() %>% 
  dplyr::select(mortality_week,
                mortality_NDVI,
                mortality_NDVI_interannual_decrease,
                mortality_NDVI_seasonal_decrease,
                mortality_LST,
                mortality_LST_residual,
                pixel_loc) %>%
  merge(all_data_with_previous_month_and_year, by=c("pixel_loc"))

ggplot() + 
  geom_density_2d_filled(data=sudden_mortality_2016_all_data %>%
                           filter(year==2016, week==(mortality_week-4)),
                         aes(x=mortality_NDVI_interannual_decrease, y=LST_residual))


getMortalityPrediction <- function(week_offset){
  linmod <- summary(lm(data=sudden_mortality_2016_all_data %>%
                         filter(year==2016, week==(mortality_week-week_offset)) %>%
                         mutate(NDVI),
                       mortality_NDVI_interannual_decrease ~ LST+LST_residual+LST_previous_year+LST_previous_month+mortality_week))
  return( data.frame(week_offset = week_offset,
                     intercept = linmod$coefficients[1,1],
                     LST_slope = linmod$coefficients[2,1],
                     LST_residual_slope = linmod$coefficients[3,1],
                     LST_previous_year = linmod$coefficients[4,1],
                     LST_previous_month = linmod$coefficients[5,1],
                     mortality_week = linmod$coefficients[6,1],
                     r_sqd = linmod$adj.r.squared) )
}

prediction_information <- bind_rows(lapply(-20:20, getMortalityPrediction))

ggplot(prediction_information) + 
  geom_smooth(aes(x=week_offset, y=r_sqd)) + 
  geom_hline(yintercept=0, linetype="dashed", col="gray") +
  scale_x_continuous(limits=c(-20, 0), expand=c(0,0)) + 
  scale_y_continuous(limits=c(0, 0.8), expand=c(0,0)) +
  xlab("Weeks in Advance") + 
  ylab("Strength of Prediction (R^2)") + 
  theme_bw() + 
  ggtitle("Strength of LST Prediction of Impending Mortality")

ggplot(prediction_information) + 
  geom_smooth(aes(x=week_offset, y=LST_residual_slope)) + 
  geom_hline(yintercept=0, linetype="dashed", col="gray") +
  scale_x_continuous(limits=c(-20, 0), expand=c(0,0)) + 
  scale_y_continuous(expand=c(0,0)) + 
  xlab("Weeks in Advance") + 
  ylab("Strength of Prediction (R^2)") + 
  theme_bw() + 
  ggtitle("Strength of LST Prediction of Impending Mortality")




# Compare all data from one year
target_year <- 2016
late_2016_decline <- all_data_with_previous_month_and_year %>%
  filter(year == target_year, week == 37) %>% 
  mutate(NDVI_norm_late = NDVI - NDVI_mean,
         interannual_decline = NDVI - NDVI_previous_year) %>%
  dplyr::select(pixel_loc, NDVI_norm_late, interannual_decline)
early_season_NDVI_norm <- (all_data_with_previous_month_and_year %>%
  filter(year == target_year, week == 22) %>% 
  mutate(NDVI_norm = NDVI - NDVI_mean))$NDVI_norm
late_2016_decline$seasonal_decline <- late_2016_decline$NDVI_norm_late - early_season_NDVI_norm

all_data_with_decline <- all_data_with_previous_month_and_year %>%
  merge(late_2016_decline, by="pixel_loc")
summary(lm(data = all_data_with_decline %>% filter(year==target_year, week==24),
           interannual_decline ~ LST+LST_residual+LST_previous_year+LST_previous_month+NDVI_previous_year))

# Regression we want to test...
#    for each year, get the interannual and the seasonal decline in greenness, plus the date when seasonal decline onset (how?)
#    then, combine with previous 1 year of advance data + current year (inc. SPEI)
# Predict each type of decline based on 
#   LST_residual*NDVI + previous_month(LST_residual*NDVI) + previous_year(LST_residual*NDVI)
#   test strength with lag before decline occurs (by week within year, or by years before)
getData
