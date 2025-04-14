# LST-GV Riparian Model
R-based model to predict riparian vegetation cover and relative temperature.

First commit. Currently, primary functions are contained in the file 'lst_gv_model.R' and the R project 'lst_riparian_model_R.Rproj'

This will currently:
1. Load GV and relative LST data (these need to be generated beforehand with Earth Engine), plus SPEI (from SPEI drought monitor)
2. Generate animated .gif images showing the change year-to-year in seasonality for LST and GV
3. Create a simple linear model for LST ~ GV + SPEI + Year
4. Predict LST from the model and then evaluate monthly residuals away from predictions -> indicator of drought stress
5. Model GV loss (interannually and from start to end of season) during drought based on initial GV, LST, and LST residual.