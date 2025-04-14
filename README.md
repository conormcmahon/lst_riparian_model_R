# LST-GV Riparian Model
R-based model to predict riparian vegetation cover and relative temperature.

First commit. Currently, primary functions are contained in the file 'lst_gv_model.R' and the R project 'lst_riparian_model_R.Rproj'

This will currently:
1. Load GV and relative LST data (these need to be generated beforehand with Earth Engine), plus SPEI (from SPEI drought monitor)
2. Generate animated .gif images showing the change year-to-year in seasonality for LST and GV
3. Create a simple linear model for LST ~ GV + SPEI + Year
4. Predict LST from the model and then evaluate monthly residuals away from predictions -> indicator of drought stress
5. Model GV loss (interannually and from start to end of season) during drought based on initial GV, LST, and LST residual.

# Animated Green Vegetation Fraction
Bands are: 
Red = Week 16 GV
Green = Week 28 GV
Blue = Week 40 GV

![Alt Text](https://github.com/conormcmahon/lst_riparian_model_R/blob/main/gv_animation.gif)

# Animated Relative Land Surface Temperature (vs. Air Temperature)
Bands are: 
Red = Week 16 GV
Green = Week 28 GV
Blue = Week 40 GV

![Alt Text](https://github.com/conormcmahon/lst_riparian_model_R/blob/main/lst_animation.gif)

# R^2 Performance of different models, depending on which parameters are EXCLUDED: 

![Alt Text](https://github.com/conormcmahon/lst_riparian_model_R/blob/main/output_plots/LST_model_intercomparison_r_sqd.png)

# SLOPE of Relationship Between LST and SPEI for Two Time Intervals

![Alt Text](https://github.com/conormcmahon/lst_riparian_model_R/blob/main/output_plots/LST_spei_dependence.png)

This probably suggests increased dependence on short-term water availability during the early part of the season (rains of that winter influence spring LST) vs. dependence on long-term water availability late in the growing season. 

# SLOPE of Relationship Between LST and GV

![Alt Text](https://github.com/conormcmahon/lst_riparian_model_R/blob/main/output_plots/LST_GV_dependence.png)

This probably suggests increased importance of GV in determining LST during periods when the deciduous trees on the river are leaf-on, and also when there is greater radiative forcing availabile to drive evapotranspiration. 
