a <- tar_read(multivariable_blca)

b <- a |> 
  mutate(
    base_level = map2(fit_anova, names, function(x, y) {
      x <- x$xlevels
      y <- enframe(y)
      x <- lapply(x, \(x) x[1]) |> 
        enframe()
      full_join(y, x, by = c("value" = "name"))
    })) |> 
  unnest(cols = base)
