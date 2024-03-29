context("test primarySuppression()")

test_that("primarySuppression works", {
  skip_on_cran()

  # load problem without suppressions
  problem <- sdc_testproblem(with_supps = FALSE)

  p1 <- primarySuppression(problem, type = "freq", maxN = 2)
  p1.sdc <- get.problemInstance(p1@problemInstance, "sdcStatus")
  expect_is(p1, "sdcProblem")
  expect_equal(sum(p1.sdc == "u"), 1)
  expect_equal(which(p1.sdc == "u"), 6)

  expect_error(primarySuppression(problem, type = "nk", n = 2, k = 80))
  p2 <- primarySuppression(problem, type = "nk", n = 2, k = 75, numVarName = "val")
  p2.sdc <- get.problemInstance(p2@problemInstance, "sdcStatus")
  expect_is(p2, "sdcProblem")
  expect_equal(sum(p2.sdc == "u"), 2)
  expect_equal(which(p2.sdc == "u"), c(6, 14))

  expect_error(primarySuppression(problem, type = "p", p = 80))
  expect_error(primarySuppression(problem, type = "p", p = 0, numVarName = "val"))
  expect_error(primarySuppression(problem, type = "p", p = 100, numVarName = "val"))
  p3 <- primarySuppression(problem, type = "p", p = 70, numVarName = "val")
  p3.sdc <- get.problemInstance(p3@problemInstance, "sdcStatus")
  expect_is(p3, "sdcProblem")
  expect_equal(sum(p3.sdc == "u"), 2)
  expect_equal(which(p3.sdc == "u"), c(6, 14))

  expect_error(primarySuppression(problem, type = "pq", pq = c(80, 90)))
  expect_error(primarySuppression(problem, type = "pq", pq = c(85, 80), numVarName = "val"))
  expect_error(primarySuppression(problem, type = "pq", pq = c(85, 180), numVarName = "val"))
  expect_error(primarySuppression(problem, type = "pq", pq = c(110, 120), numVarName = "val"))
  p4 <- primarySuppression(problem, type = "pq", pq = c(60, 80), numVarName = "val")
  p4.sdc <- get.problemInstance(p4@problemInstance, "sdcStatus")
  expect_is(p4, "sdcProblem")
  expect_equal(sum(p4.sdc == "u"), 2)
  expect_equal(which(p4.sdc == "u"), c(6, 14))

  problem@dataObj@rawData$val[5] <- -5
  expect_error(primarySuppression(problem, type = "nk", n = 2, k = 80, numVarName = "val"))
  expect_error(primarySuppression(problem, type = "p", p = 70, numVarName = "val"))
  expect_error(primarySuppression(problem, type = "pq", pq = c(60, 80), numVarName = "val"))


  # use weights
  rm(list = ls())
  utils::data("microdata1", package = "sdcTable")
  microdata1$samp_weights <- sample(rnorm(nrow(microdata1), mean = 10))
  dim_region <- sdcHierarchies::hier_create(
    root = "Total",
    nodes = LETTERS[1:5]
  )
  dim_gender <- sdcHierarchies::hier_create(
    root = "Total",
    nodes = c("male", "female")
  )
  dimList <- list(
    region = dim_region,
    gender = dim_gender
  )
  problem <- makeProblem(
    data = microdata1,
    dimList = dimList,
    freqVarInd = NULL,
    numVarInd = 3,
    sampWeightInd = 4)

  # no numVarInd or numVarName
  expect_error(
    primarySuppression(problem, type = "nk", n = 2, k = 80)
  )
  # variable does not exist
  expect_error(
    primarySuppression(problem, type = "nk", n = 2, k = 80, numVarName = "samp_weight")
  )

  # n < 2
  expect_error(
    primarySuppression(problem, type = "nk", n = 1, k = 80, numVarName = "samp_weights")
  )

  res <- primarySuppression(
    object = problem,
    type = "nk",
    n = 2,
    k = 80,
    numVarName = "val")

  res <- sdcProb2df(res, addDups = FALSE, addNumVars = TRUE)
  expect_equal(res[sdcStatus == "s", .N], 15)
  expect_equal(res[sdcStatus == "z", .N], 3)


  res <- primarySuppression(
    object = problem,
    type = "nk",
    n = 2,
    k = 80,
    allowZeros = TRUE,
    numVarName = "val"
  )
  res <- sdcProb2df(res, addDups = FALSE)

  expect_error(primarySuppression(
    object = problem,
    type = "p",
    p = -5,
    numVarName = "samp_weights"
  ))

  expect_error(primarySuppression(
    object = problem,
    type = "p",
    p = c(-5, 5),
    numVarName = "samp_weights"
  ))

  # print method
  print(res)
})
