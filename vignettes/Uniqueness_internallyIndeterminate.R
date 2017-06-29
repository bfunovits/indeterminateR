## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)
library(indeterminateR)
library(tidyverse)

## ----Run the example-----------------------------------------------------
results <- evaluate_my_model(deep_params = 
                               list(s = 1.5, 
                                    b = 0.99, 
                                    k = 1/1.5, 
                                    f_y = -1.3/1.5, 
                                    f_p = 1.3, 
                                    f_r = 0.6),
                             my_model = "Lubik_Marzo_26_0", 
                             verbose = TRUE)

## ------------------------------------------------------------------------
t(results)

