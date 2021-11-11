library(targets)
library(survival)


temp <- tar_read(combined_survs) |> 
  mutate(sex = "all")
a <- tar_read(combined_survs) |> 
  bind_rows(temp)


my_survfit <- function(df) {
  survfit(Surv(follow_up_time, death) ~ b_bin, data = df) |> 
    surv_pvalue(data = df)
}


pvals <- a |> 
  group_by(project, sex) |> 
  nest() |> 
  mutate(surv_diff_p = map(data, my_survfit)) |> 
  unnest(surv_diff_p)

main <- dplyr::filter(a, project %in% c("hnsc", "lihc", "luad", "skcm")) |> 
  mutate(sex = if_else(sex == "all", "All", sex),
         project = toupper(project))


# By b_bin, sex, project
b <- survfit(Surv(follow_up_time, death) ~ project + sex + b_bin, data = main)



c <- ggsurvplot(b, color = "b_bin", size = 0.2, 
                censor.size = 1, 
                break.time.by = 365.25, 
                xscale = "d_y")

d <- c$plot

d + facet_grid(project~sex) + 
  theme(legend.position = "none") + 
  bladdr::theme_tufte(10) + 
  coord_cartesian(xlim = c(0, 3653)) +
  scale_color_manual(values = c("red", "blue")) +
  guides(color = guide_legend(title = "B-cell Gene Signature"))


# By sex, project
b <- survfit(Surv(follow_up_time, death) ~ project + sex, data = a)

c <- ggsurvplot(b, color = "sex", size = 0.2, censor.size = 1)

d <- c$plot

d + facet_grid(rows = "project") + 
  theme(legend.position = "none") + 
  bladdr::theme_tufte(10) + 
  coord_cartesian(xlim = c(0, 6000)) +
  scale_color_manual(values = c("black", "blue", "red"))
# walk using autoplot &c from here

