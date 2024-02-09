library(data.table)
devtools::load_all(".")

message("*** Saving parameterset esk2020.model to data/esk2020.model.RData")
esk2020.model <- san_model(rates=data.table(Tmax = c(3, 6, 11, 40),
                                         r_S = c(0, 0.6, 0.939572775547456, 1.67598732913722),
                                         r_0 = c(0.35, 0.0, 0.0, 0.0),
                                         r_R = c(0.15, 0.15, 1.13951204627429, 0.0),
                                         r_A = c(0.0, 0.0, 0, 1.69289943912353),
                                         r_N = c(0.0, 0.0, 0.0, 0.708099388708624),
                                         r_D = c(0.0, 0.0, 0.0, 0.072323413732001)),
                        s0 = 30000)
save(esk2020.model, file="data/esk2020.model.RData", version=2)
