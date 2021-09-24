## ULHAdominica
# Unified Landslide Hazard Assessment: an Island of Dominica case study.

R code to reproduce the analysis in "Unified landslide hazard assessment using hurdle models: a case study in the Island of Dominica" by Erin Bryce, Luigi Lombardo and Daniela Castro-Camilo.

The article proposes a method to simultaneously address two main concepts in landslide hazard assessment: susceptibility (or probability of occurrence) and size. Our model, which can estimate where and how large a given landslide might be, corresponds to a flexible Bernoulli-log-Gaussian hurdle model. To cope with the high spatial resolution of the data, we use a Markovian representation of the Matérn covariance function based on the stochastic partial differential equation (SPDE) approach. Assuming Gaussian priors, our model can be integrated into the class of latent Gaussian models, for which inference is conveniently performed based on the integrated nested Laplace approximation method. We test our modelling approach in Dominica, where Hurricane Maria (September 2017) induced thousands of shallow flow-like landslides passing over the island.

Here is a description of the files available:
- ULHAcode.R: R code to run the two parts of the hurdle model (Bernoulli for landslide susceptibility and Gaussian for landslide log-size).
- ULHAlandslides.csv: the original data. Additional data for analysis is created in ULHAcode.R.
- ULHAcoords.csv: locations (in degrees) of the SU centroids.
