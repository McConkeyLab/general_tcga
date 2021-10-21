library(targets)
library(tidymodels)
library(censored)

a <- tar_read(surv_tidy_blca)

names <- c("follow_up_time", "death", "sex", "age", "path_stage")

rec <- recipe(a) |>   
  step_select(names) |> 
  step_dummy(all_nominal()) 
  

prepped <- prep(rec)

baked <- bake(prepped, new_data = NULL)


fit_mod <- function(df) {
  proportional_hazards() |> 
    fit(Surv(follow_up_time, death) ~ ., data = df)
}


b <- tibble(data = list(baked)) |> 
  mutate(cox_fit = map(data, fit_mod),
         tidy_fit = map(cox_fit, tidy, conf.int = TRUE, exponentiate = TRUE))


function (mod, type = c("II", "III", 2, 3), test.statistic = c("LR", 
                                                               "Wald"), ...) 
{
  type <- as.character(type)
  type <- match.arg(type)
  test.statistic <- match.arg(test.statistic)
  if (length((mod$rscore) > 0) && (test.statistic == "LR")) {
    warning("LR tests unavailable with robust variances\n  Wald tests substituted")
    test.statistic <- "Wald"
  }
  names <- term.names(mod)
  clusters <- grepl("cluster\\(", names)
  strata <- grepl("strata\\(", names)
  if ((any(clusters) || any(strata)) && test.statistic == 
      "LR") {
    warning("LR tests not supported for models with clusters or strata\n Wald tests substituted")
    test.statistic <- "Wald"
  }
  switch(type, 
         II = switch(test.statistic, 
                     LR = Anova.II.LR.coxph(mod), 
                     Wald = Anova.default(mod, type = "II", test.statistic = "Chisq", 
                                          vcov. = vcov(mod, complete = FALSE))), 
         III = switch(test.statistic, 
                      LR = Anova.III.LR.coxph(mod), 
                      Wald = Anova.default(mod, type = "III", test.statistic = "Chisq", 
                                           vcov. = vcov(mod, complete = FALSE))), 
         `2` = switch(test.statistic, 
                      LR = Anova.II.LR.coxph(mod), 
                      Wald = Anova.default(mod, type = "II", test.statistic = "Chisq", 
                                           vcov. = vcov(mod, complete = FALSE))), 
         `3` = switch(test.statistic, 
                      LR = Anova.III.LR.coxph(mod), 
                      Wald = Anova.default(mod, type = "III", test.statistic = "Chisq", 
                                           vcov. = vcov(mod, complete = FALSE))))
}
