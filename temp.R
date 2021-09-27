library(bladdr)

df <- tar_read(combined_hrs)

df_summ <- df |>
  group_by(group, term) |> 
  summarize(mean = mean(estimate), 
            sd = sd(estimate))

df_plot <- df |> 
  mutate(project = toupper(project),
         group = case_when(group == "all" ~ "All",
                           group == "female" ~ "F",
                           group == "male" ~ "M"),
         group = factor(group, levels = c("F", "M", "All")),
         term = case_when(term == "b_binHi" ~ "B-cell Signature",
                          term == "cd8_binHi" ~ "CD8+ T-cell Signature",
                          term == "imm_binHi" ~ "Pan-Immune Signature"))

ggplot(df_plot, aes(x = estimate, y = group, color = group)) +
  geom_vline(xintercept = 1, alpha = 0.5) + 
  scale_color_manual(values = c("#F8B7CD", "#0671B7", "black")) + 
  geom_linerange(aes(xmin = conf.low, xmax = conf.high)) + 
  geom_point() + 
  facet_grid(project~term) + 
  bladdr:::theme_tufte(10) + 
  labs(x = "Hazard Ratio") + 
  theme(legend.position = "none",
        axis.title.y = element_blank())

df_2d <- df |> 
  select(-c(std.error, statistic, p.value)) |> 
  pivot_wider(names_from = term, values_from = c(estimate, conf.low, conf.high)) |> 
  mutate(alpha = case_when(conf.high_cd8_binHi < estimate_b_binHi ~ 0.8,
                           conf.low_cd8_binHi > estimate_b_binHi ~ 0.8,
                           conf.high_b_binHi < estimate_cd8_binHi ~ 0.8,
                           conf.low_b_binHi > estimate_cd8_binHi ~ 0.8,
                           T ~ 0.1))

ggplot(df_2d, aes(x = estimate_b_binHi, y = estimate_cd8_binHi, shape = group, color = project)) +
  geom_abline(slope = 1, intercept = 0, alpha = 0.1, size = 2) +
  geom_hline(yintercept = 1, alpha = 0.1, size = 2) + 
  geom_vline(xintercept = 1, alpha = 0.1, size = 2) + 
  scale_alpha_continuous(limits = c(0,1)) + 
  geom_errorbar(aes(xmin = conf.low_b_binHi, xmax = conf.high_b_binHi, alpha = alpha)) +
  geom_errorbar(aes(ymin = conf.low_cd8_binHi, ymax = conf.high_cd8_binHi, alpha = alpha)) + 
  geom_point(size = 2, aes(alpha = alpha)) +
  coord_cartesian(xlim = c(0, 2), ylim = c(0, 2)) + 
  bladdr:::theme_tufte(10)

ggplot(df_2d, aes(x = log2(estimate_b_binHi), y = log2(estimate_cd8_binHi), shape = group, color = project)) + 
  geom_abline(slope = 1, intercept = 0, alpha = 0.1, size = 2) +
  geom_hline(yintercept = 0, alpha = 0.1, size = 2) + 
  geom_vline(xintercept = 0, alpha = 0.1, size = 2) + 
  geom_point(size = 2) +
  bladdr:::theme_tufte(10) + 
  coord_cartesian(xlim = c(-1.5, 1.5), y = c(-1.5, 1.5))
