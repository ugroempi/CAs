library(CAs)

test_that("bestCA-based constructions record their ingredient CAs", {
  composed <- compositCA(2, 3, 6)
  composed_ingredients <- attr(composed, "ingredients")

  expect_named(composed_ingredients, c("D1", "D2"))
  expect_named(composed_ingredients$D1, c("constr", "t", "k", "v"))
  expect_named(composed_ingredients$D2, c("constr", "t", "k", "v"))
  expect_identical(composed_ingredients$D1[c("t", "k", "v")],
                   list(t=2, k=3, v=2))
  expect_identical(composed_ingredients$D2[c("t", "k", "v")],
                   list(t=2, k=3, v=3))
  expect_true(composed_ingredients$D1$constr %in% names(Ns(2, 3, 2)))
  expect_true(composed_ingredients$D2$constr %in% names(Ns(2, 3, 3)))

  doubled <- CK_doublingCA(8, 2)
  doubled_ingredients <- attr(doubled, "ingredients")

  expect_named(doubled_ingredients, c("D3", "D2"))
  expect_named(doubled_ingredients$D3, c("constr", "t", "k", "v"))
  expect_named(doubled_ingredients$D2, c("constr", "t", "k", "v"))
  expect_identical(doubled_ingredients$D3[c("t", "k", "v")],
                   list(t=3, k=4, v=2))
  expect_identical(doubled_ingredients$D2[c("t", "k", "v")],
                   list(t=2, k=4, v=2))
})
