testthat::context("digest_args")
library(ssvQC)

test_get_args.dots = function(a = 1, b = 2, ...){
    length(get_args(...))
}

test_get_args = function(a = 1, b = 2){
    length(get_args())
}

test_digest_args.dots = function(a = 1, b = 2, ...){
    digest_args(...)
}

test_digest_args = function(a = 1, b = 2){
    digest_args()
}

test_that("get_args", {
    expect_equal(test_get_args(), 2)
    expect_equal(test_get_args.dots(), 2)
    expect_equal(test_get_args.dots(extra = 3), 3)
    expect_equal(test_get_args.dots(extra = 3, extra2 = 4), 4)
})

test_that("digest_args", {
    expect_true(test_digest_args.dots() == test_digest_args())
    expect_false(test_digest_args.dots(extra = 3) == test_digest_args())
    expect_false(test_digest_args.dots(extra = 3) == test_digest_args.dots(extra = 4))
    expect_false(test_digest_args.dots(extra = 3) == test_digest_args.dots(extra2 = 3))
    expect_true(test_digest_args.dots(extra1 = 3, extra2 = 4) == test_digest_args.dots(extra2 = 4, extra1 = 3))
})
