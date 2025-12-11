test_that("NAs in inputs are rejected", {
  Sim <- data.frame(Q=c(1, 2, NA, 4), year=1:4,ID=1)
  Obs <-data.frame(Q= c(1, 2,  3, 4), P=1:4,year="1:4",ID=1)
  expect_error(BRRAT(Sim, Obs), regexp = "NA", ignore.case = TRUE)
})

test_that("Obs <=0 is rejected",{
	Sim <- data.frame(Q=1:10,year=1:10,ID=1)
	Obs <- data.frame(Q=-1:8,P=1:10,year=1:10,ID=1)
	expect_error(BRRAT(Sim, Obs), regexp = "Obs|positive", ignore.case = TRUE)
})

test_that("Inputs are data frame",{
  Sim <- 1:10
  Obs <- 1:10
  expect_error(BRRAT(Sim, Obs),regexp = "data.frame", ignore.case = TRUE)
})

test_that("data frame has requried headers",{
  Sim <- data.frame(S=1:10,year=1:10,ID=1)
  Obs <- data.frame(A=-1:8,P=1:10,year=1:10,ID=1)
  expect_error(BRRAT(Sim, Obs), regexp = "columns", ignore.case = TRUE)
})
