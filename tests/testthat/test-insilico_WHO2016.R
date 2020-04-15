context("test-insilico_WHO2016")

test_that("successful run with RandomVA5", {
    load('insilico_WHO2016.RData')
    expect_true(class(fit1) == "insilico")
    expect_true(class(summary(fit1) )== "insilico_summary")
    expect_true(nrow(fit1.indiv[[1]]) == n)
    expect_true(ncol(fit1.indiv[[1]]) == 61)
})
