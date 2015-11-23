context("Functionality that probably belongs in the S4Vectors package")

test_that("combine,DataFrame,DataFrame-method works", {
  expect_identical(combine(S4Vectors::DataFrame(), S4Vectors::DataFrame()),
                   S4Vectors::DataFrame())
  # Based on ?BiocGenerics::combine
  x <- data.frame(x = 1:5,
                  y = factor(letters[1:5], levels = letters[1:8]),
                  row.names = letters[1:5])
  X <- S4Vectors::DataFrame(x)
  y <- data.frame(z = 3:7,
                  y = factor(letters[3:7], levels = letters[1:8]),
                  row.names = letters[3:7])
  Y <- S4Vectors::DataFrame(y)
  expect_identical(combine(X, Y), DataFrame(combine(x, y)))
  expect_identical(combine(X, DataFrame()), X)
  expect_identical(combine(DataFrame(), Y), Y)
  YY <- Y
  row.names(YY) <- NULL
  expect_error(combine(X, YY), "'row.names' of 'x' and 'y' must be non-NULL")

  w <- data.frame(w = 4:8,
                  y = factor(letters[4:8], levels = letters[1:8]),
                  row.names = letters[4:8])
  W <- DataFrame(w)
  expect_identical(combine(W, X, Y), DataFrame(combine(w, x, y)))

  # y is converted to 'factor' with different levels
  df1 <- data.frame(x = 1:5, y = letters[1:5], row.names = letters[1:5])
  DF1 <- S4Vectors::DataFrame(df1)
  df2 <- data.frame(z = 3:7, y = letters[3:7], row.names = letters[3:7])
  DF2 <- S4Vectors::DataFrame(df2)
  expect_warning(combine(DF1, DF2)) # Fails
  # solution 1: ensure identical levels
  y1 <- factor(letters[1:5], levels = letters[1:7])
  y2 <- factor(letters[3:7], levels = letters[1:7])
  df1 <- data.frame(x = 1:5, y = y1, row.names = letters[1:5])
  DF1 <- S4Vectors::DataFrame(df1)
  df2 <- data.frame(z = 3:7, y = y2, row.names = letters[3:7])
  DF2 <- S4Vectors::DataFrame(df2)
  expect_identical(combine(DF1, DF2), DataFrame(combine(df1, df2)))
  # solution 2: force column to be 'character'
  df1 <- data.frame(x = 1:5, y = I(letters[1:5]), row.names = letters[1:5])
  DF1 <- S4Vectors::DataFrame(df1)
  df2 <- data.frame(z = 3:7, y = I(letters[3:7]), row.names = letters[3:7])
  DF2 <- S4Vectors::DataFrame(df2)
  expect_identical(combine(DF1, DF2), DataFrame(combine(df1, df2)))
})

test_that("combine,SimpleList,SimpleList-method works", {
  x <- SimpleList(C = matrix(10:1, ncol = 2),
                  D = S4Vectors::DataFrame(x = 1:5,
                                           y = letters[1:5],
                                           row.names = LETTERS[1:5]))
  expect_identical(combine(x, SimpleList()), x)
  expect_identical(combine(SimpleList(), x), x)
  expect_error(combine(x[1], x[2]),
               "'SimpleList' objects have different element names:")
  expect_error(combine(x, x[1]),
               "'SimpleList' objects have different number of elements:")
  y <- SimpleList(C = matrix(c(10:7, 5:2), ncol = 2),
                  D = S4Vectors::DataFrame(x = 1:3,
                                           y = letters[1:3],
                                           row.names = LETTERS[1:3]))
  expect_error(combine(x, y), "matricies must have dimnames for 'combine'")
  xx <- SimpleList(C = matrix(10:1, ncol = 2,
                              dimnames = list(LETTERS[5:1], letters[1:2])),
                  D = S4Vectors::DataFrame(x = 1:5,
                                           y = letters[1:5],
                                           row.names = LETTERS[1:5]))
  yy <- SimpleList(C = matrix(c(10:7, 5:2), ncol = 2,
                             dimnames = list(LETTERS[5:2], letters[1:2])),
                  D = S4Vectors::DataFrame(x = 1:3,
                                           y = letters[1:3],
                                           row.names = LETTERS[1:3]))
  expect_identical(combine(xx, yy), xx)
})
