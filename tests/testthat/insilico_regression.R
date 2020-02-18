rm(list=ls())

library(InSilicoVA)

data("RandomVA2", envir = environment())
RandomVA2 <- get("RandomVA2", envir = environment())
data <- RandomVA2
rm(RandomVA2)

# run without sub-population
regress <- insilico( data, subpop = NULL,
                  Nsim = 400, burnin = 200, thin = 10 , seed = 13,
                  auto.length = FALSE, warning.write = TRUE)

save.image("insilico_regression.RData")
