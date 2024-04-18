library(tidyverse)

df <- read_csv("data/Artificial_proteins_Knot_type.csv") %>%
  mutate(Major_knot = recode(Major_knot, 'TooManyCrossings' = 'N/A', '0_1' = 'N/A')) %>%
  inner_join(read_csv("data/minimal.csv"), by="ID") %>%
  filter(Tool!="Real") %>%
  arrange(Major_knot)

df2 <- starwars <- df %>% 
  group_by(Tool, Major_knot) %>%
  summarise(n = n()) %>%
  pivot_wider(names_from = Major_knot, values_from = n)

df2[, c(1,1+order(colnames(df2)[-1]))]
