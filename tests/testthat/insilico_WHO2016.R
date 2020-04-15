data(RandomVA5)
fit1 <- insilico(RandomVA5, data.type="WHO2016", datacheck = FALSE,
                 Nsim = 400, burnin = 100, thin = 10, seed = 1,
                 auto.length = FALSE)
fit1.indiv <- get.indiv(fit1)
n <- nrow(RandomVA5)
save.image('insilico_WHO2016.RData')
