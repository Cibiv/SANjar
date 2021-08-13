library(data.table)
devtools::load_all(".")

message("*** Saving parameterset lt47.model to data/lt47.model.RData")
lt47.model <- san_model(rates=data.table(Tmax = c(3, 6, 11, 40),
                                         r_S = c(0, 0.6, 0.939572775547456, 1.67598732913722),
                                         r_0 = c(0.35, NA, NA, NA),
                                         r_R = c(0.15, 0.15, 1.13951204627429, NA),
                                         r_A = c(NA, NA, 0, 1.69289943912353),
                                         r_N = c(NA, NA, NA, 0.708099388708624),
                                         r_D = c(NA, NA, NA, 0.072323413732001)),
                        s0 = 30000)
save(lt47.model, file="data/lt47.model.RData", version=2)
