###### custom function to add filename when reading data from excel spreadsheets

required <- list("tidyverse", 
                 "readxl", 
                 "ggbeeswarm",
                 "ggdark",
                 "colorspace",
                 "ggdist",
                 "lme4",
                 "lmerTest",
                 "MuMIn",
                 "ggh4x",
                 "glmmTMB",
                 "emmeans",
                 "ungeviz")

lapply(required, library, character.only = T)

library(tidyverse)
library(readxl)
library(ggbeeswarm)
library(ggdark) #desaturate()
library(colorspace) #lighten()
library(ggdist) # rainclouds
library(lme4)
library(lmerTest)
library(MuMIn) # dredge
library(ggh4x) # facet nested wrap
library(glmmTMB)
library(gt)
library(gtsummary)
library(patchwork)


## Functions

read_plus <- function(flnm, sheet, skip) {
  read_excel(flnm ,sheet=sheet, skip=skip, na=c("na","NA", "-")) %>% 
    mutate(filename = flnm)
  
}

#####

##### force bind, bind rows of data when column names do not match #####

force_bind = function(df1, df2) {
  colnames(df2) = colnames(df1)
  bind_rows(df1, df2)
}

#####


##### DHARMa_check function to simulate residuals and plot - for mixed models

DHARMa_check <- function(model){
  sim <- DHARMa::simulateResiduals(model) 
  plot(sim, asFactor=T)
}


##### binned plot to check overdispersion in binomial model

bin_plot <- function(model){
  x <- predict(model)
  y <- resid(model)
  arm::binnedplot(x,y)
}



#' Nudge any other positioning function from ggplot2
#' @param x Nudge in X direction
#' @param y Nudge in Y direction
#' @param position Any positioning operator from ggplot2 like \link{position_jitter}[ggplot2]
#' @return A combination positioning operator which first runs the original positioning,
#' then adds a nudge on top of that
#' @export
position_nudge_any <- function(x = 0, y = 0, position) {
  ggproto(NULL, PositionNudgeAny,
          nudge = ggplot2::position_nudge(x, y),
          position = position
  )
}

# Beeswarm nudge====

#' Internal class doing the actual nudging on top of the other operation
#' @keywords internal
PositionNudgeAny <- ggplot2::ggproto("PositionNudgeAny", ggplot2::Position,
                                     nudge = NULL,
                                     nudge_params = NULL,
                                     position = NULL,
                                     position_params = NULL,
                                     
                                     setup_params = function(self, data) {
                                       list(nudge = self$nudge,
                                            nudge_params = self$nudge$setup_params(data),
                                            position = self$position,
                                            position_params = self$position$setup_params(data))
                                     },
                                     
                                     setup_data = function(self, data, params) {
                                       data <- params$position$setup_data(data, params$position_params)
                                       params$nudge$setup_data(data, params$nudge_params)
                                     },
                                     
                                     compute_layer = function(self, data, params, layout) {
                                       data <- params$position$compute_layer(data, params$position_params, layout)
                                       params$nudge$compute_layer(data, params$nudge_params, layout)
                                     }
)