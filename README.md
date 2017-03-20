# vis-qrodm
Plot quantile regression and outlier diagnostics models in high-D space

Key notes: 

When visualizing QR model:

- orange: quantile regression model on 0.1 quantile
- green: quantile regression model on 0.5 quantile
- purple: quantile regression model on 0.9 quantile

When visualizing outlier diagnostic models:

- green: normal points
- purple: outliers


## R file

### plot-qr-model: contains code for plotting QR model

* linear-3d-plane.R : plotting linear regression model with 2 predictors

* linear-4d-cuboid.R : plotting linear regression model with 3 predictors

* non-linear-qr1.R : plotting non-linear regression model of elliptic paraboloid or hyperbolic paraboloid surface

* non-linear-qr2.R : plotting non-linear regression model of hyperbolic paraboloid surface

* non-linear-qr3.R : plotting non-linear regression model of elliptic paraboloid surface

### visualization of outlier diagnostic models

* simulating-points.R : generating points in space

* non-linear-diagnostic.R : visualize outlier diagnostic models for non-linear quantile regression (use non-linear curve SSlogis for example)

* linear-diagnostic.R : visualize outlier diagnostic models (we focus on two models: generalized cook distance and q-function distance) for linear quantile regression


## Figures file

### QR-model

This file contains figures of QR-models in high space

### Outlier-Diagnostic

This file contains figures of outlier diagnostic models for quantile regression


## Data file

### linear-2D

This file contains data used in plotting outlier diagnostic models for linear quantile regression with 1 predictor (only simulated data (sim_data in code), original data not included, we combined the original data in our code)

### linear-3D

This file contains data used in plotting outlier diagnostic models for linear quantile regression with 2 predictor (only simulated data (sim_data in code), original data not included, we combined the original data in our code)

### linear-4D

This file contains data used in plotting outlier diagnostic models for linear quantile regression with 3 predictor (only simulated data (sim_data in code), original data not included, we combined the original data in our code)

### non-linear

Will be added

## Markdown file

Our paper



