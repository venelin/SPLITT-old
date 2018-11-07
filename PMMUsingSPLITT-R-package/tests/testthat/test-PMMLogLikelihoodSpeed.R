library(testthat)
context("Test MiniBenchmark(1000, 10) passes without errors")


library(PMMUsingSPLITT)

test_that("MiniBenchmark runs without errors", {
          expect_output(df <- PMMUsingSPLITT::MiniBenchmark(1000, 10), "Measuring calculation times...")
          expect_s3_class(df, "data.frame")
          })