library(targets)
library(tidyverse)

a <- tar_read(gsva_coldata) 

all <- tar_read(gsva_coldata) |> 
  mutate(sex = "all") |> 
  bind_rows(a) |> 
  pivot_longer(cols = c(b_cell, cd8, exp_immune), names_to = "signature") |> 
  mutate(project = toupper(project),
         signature = case_when(signature == "b_cell" ~ "B-cell",
                               signature == "cd8" ~ "CD8+ T-cell",
                               signature == "exp_immune" ~ "IFNg"),
         sex = if_else(sex == "all", "All", sex),
         sex = factor(sex, levels = c("F", "M", "All"))) # To  put 'all' on top

linetypes <- c(all = "dashed", M = "solid", `F` = "solid")
colors <- c(all = "gray30", M = "#0671B7", `F` = "#F8B7CD")

ggplot(all, aes(x = value, color = sex, linetype = sex)) +
  scale_linetype_manual(values = linetypes) +
  scale_color_manual(values = colors) +
  geom_density(size = 0.5) + 
  facet_grid(project~signature) + 
  labs(x = "Gene Signature Score", y = "Density") +
  bladdr::theme_tufte(10) + 
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank())
