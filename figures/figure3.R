library(tidyverse)
library(patchwork)  # for 2x2 plot layout

my_theme <- theme_bw # replacedtheme_minimal

# plot A ------------------------------------------------------------------

evodiff <- read_csv("data/Evodiff_perplexity.csv") %>% 
  inner_join(read_csv("data/Artificial_proteins_pLDTT.csv")) %>%
  inner_join(read_csv("data/minimal.csv"))

plotA <- ggplot(evodiff %>% filter(pLDDT > 1), 
       aes(x = pLDDT, y = Perplexity, color = as.factor(Label))) + 
  geom_point(alpha = 0.5) +  # Making points partially transparent
  scale_color_manual(values = c("0" = "orange", "1" = "blue"),  # Specify colors
                     labels = c("No", "Yes"),  # Custom labels for the legend
                     name = "Knotted?") +  # Custom legend title
  xlab("pLDDT") + 
  ylab("Perplexity") + 
  ggtitle("EvoDiff: Perplexity vs. pLDDT") +
  my_theme()

# plot B ------------------------------------------------------------------

RFdiffusion <- read_csv("data/RFdiffusion_perplexity.csv") %>% 
  inner_join(read_csv("data/Artificial_proteins_pLDTT.csv")) %>%
  inner_join(read_csv("data/minimal.csv"))

plotB <- ggplot(RFdiffusion %>% filter(pLDDT > 1), 
       aes(x = pLDDT, y = Perplexity, color = as.factor(Label))) + 
  geom_point(alpha = 0.5) +  # Making points partially transparent
  scale_color_manual(values = c("0" = "#FF7F00", "1" = "blue"),  # Specify colors
                     labels = c("No", "Yes"),  # Custom labels for the legend
                     name = "Knotted?") +  # Custom legend title
  xlab("pLDDT") + 
  ylab("Perplexity") + 
  ggtitle("RFdiffusion + MPNN: Perplexity vs. pLDDT") +
  my_theme()
plotB

# plot C ------------------------------------------------------------------

blast <- read_csv("data/Artificial_proteins_Blast.csv") %>%
  inner_join(read_csv("data/minimal.csv")) %>%
  mutate(`Match_e-value` = as.numeric(`Match_e-value`),  # Convert to numeric
         match = !is.na(`Match_e-value`) & `Match_e-value` < 1)  # Determine match condition

blast_summary <- blast %>%
  group_by(interaction(Label, Tool)) %>%
  summarise(Percentage = mean(match, na.rm = TRUE) * 100) %>%
  ungroup()
blast_summary
names(blast_summary)[1] <- "Group"
blast_summary$Group <- c("Evodiff,\n unknotted", "Evodiff,\n knotted", 
                         "RFdiff+MPNN,\n unknotted", "RFdiff+MPNN,\n knotted")
blast_summary$Tool <- c("Evodiff", "Evodiff", "RFdiffusion", "RFdiffusion")

plotC <- ggplot(blast_summary, aes(x = Group, y = Percentage, fill = as.factor(Tool))) + 
  geom_bar(stat = "identity", position = "dodge") + 
  labs(title = "Percentage of BLAST matches to real protein sequences",
       x = "Tool and Knot Status Combination", 
       y = "Percentage of BLAST Matches (e-value < 1)", fill = "Label") + 
  ylim(c(0,100)) +
  theme_minimal() + 
  scale_fill_manual(values = c("RFdiffusion" = "#57db5f", "Evodiff" = "#db5f57")) +
  guides(fill = FALSE)  # to remove legend

# plot D ------------------------------------------------------------------

foldseek <- read_csv("data/Artificial_proteins_Foldseek.csv") %>%
  inner_join(read_csv("data/minimal.csv")) %>% 
  mutate(`AF_match_score` = as.numeric(`AF_match_score`),  # Convert to numeric
         match = !is.na(`AF_match_score`) & `AF_match_score` > 0)  # Determine match condition

foldseek_summary <- foldseek %>%
  group_by(interaction(Label, Tool)) %>%
  summarise(Percentage = mean(match, na.rm = TRUE) * 100) %>%
  ungroup()
foldseek_summary
names(foldseek_summary)[1] <- "Group"
foldseek_summary$Group <- c("Evodiff,\n unknotted", "Evodiff,\n knotted", 
                         "RFdiff+MPNN,\n unknotted", "RFdiff+MPNN,\n knotted")
foldseek_summary$Tool <- c("Evodiff", "Evodiff", "RFdiffusion", "RFdiffusion")

plotD <- ggplot(foldseek_summary, aes(x = Group, y = Percentage, fill = as.factor(Tool))) + 
  geom_bar(stat = "identity", position = "dodge") + 
  labs(title = "Percentage of FoldSeek matches to AlphaFold structures",
       x = "Tool and Knot Status Combination", 
       y = "Percentage of FoldSeek Matches (match-score > 0)", fill = "Label") + 
  theme_minimal() + 
  scale_fill_manual(values = c("RFdiffusion" = "#57db5f", "Evodiff" = "#db5f57")) +
  guides(fill = FALSE)  # to remove legend

# Combined plot -----------------------------------------------------------

combined_plot <- plotA + plotB + plotC + plotD +
  plot_layout(ncol = 2, nrow = 2, byrow = FALSE) +  # Arrange in a 2x2 grid
  plot_annotation(tag_levels = 'A')  # Add labels (A-D) to each plot

ggsave("figures/figure3.pdf", combined_plot, width = 12, height = 12)
