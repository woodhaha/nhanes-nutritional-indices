# Flow diagram using ggplot2
library(ggplot2)
library(dplyr)

total <- 44066
any_cancer <- 3586
no_cancer <- total - any_cancer
excluded_small <- 109
final <- 2942
gi <- 353; br <- 440; fr <- 415; pu <- 516; os <- 1218

# Build nodes
nodes <- data.frame(
  id = c("total", "no_cancer", "any_cancer", "excluded", "final",
         "gi", "br", "fr", "pu", "os"),
  label = c(
    sprintf("NHANES III + Continuous NHANES\nParticipants with nutrition & mortality data\n(N = %s)", format(total, big.mark=",")),
    sprintf("No cancer diagnosis\n(n = %s)", format(no_cancer, big.mark=",")),
    sprintf("Any cancer diagnosis\n(n = %s)", format(any_cancer, big.mark=",")),
    sprintf("Excluded: lung cancer / hematologic\nmalignancies (n = %d)", excluded_small),
    sprintf("Final analytic sample\n(N = %s)", format(final, big.mark=",")),
    sprintf("GI\nn = %d", gi),
    sprintf("Breast\nn = %d", br),
    sprintf("FemaleRepro\nn = %d", fr),
    sprintf("Prostate/\nUrinary\nn = %d", pu),
    sprintf("Other\nSolid\nn = %d", os)
  ),
  x = c(3, 1.5, 4.5, 6, 3, 1, 2, 3, 4, 5),
  y = c(9, 7.5, 7.5, 7.5, 4.5, 2.5, 2.5, 2.5, 2.5, 2.5)
)

# Connections
edges <- data.frame(
  from = c("total", "any_cancer", "any_cancer", "any_cancer", "final", "final", "final", "final", "final"),
  to = c("any_cancer", "excluded", "final", "no_cancer", "gi", "br", "fr", "pu", "os")
)

# Fix: total goes to no_cancer + any_cancer
# any_cancer goes to excluded + final

p <- ggplot() +
  # Boxes
  geom_rect(data = nodes[nodes$id %in% c("total","any_cancer","final"),],
            aes(xmin = x - 1.8, xmax = x + 1.8, ymin = y - 0.5, ymax = y + 0.5),
            fill = "#E8F0FE", color = "black", linewidth = 0.3) +
  geom_rect(data = nodes[nodes$id %in% c("excluded","no_cancer"),],
            aes(xmin = x - 1.3, xmax = x + 1.3, ymin = y - 0.4, ymax = y + 0.4),
            fill = "#FCE8E6", color = "black", linewidth = 0.3) +
  geom_rect(data = nodes[nodes$id %in% c("gi","br","fr","pu","os"),],
            aes(xmin = x - 0.7, xmax = x + 0.7, ymin = y - 0.5, ymax = y + 0.5),
            fill = "white", color = "black", linewidth = 0.3) +
  # Labels
  geom_text(data = nodes, aes(x = x, y = y, label = label),
            size = c(3.2, 2.8, 2.8, 2.8, 3.2, rep(2.8, 5)), lineheight = 0.9) +
  # Arrows (simplified with segments)
  # total -> no_cancer (left)
  geom_segment(aes(x = 2.5, y = 8.5, xend = 1.1, yend = 7.9), linewidth = 0.3, arrow = arrow(length = unit(2, "mm"))) +
  # total -> any_cancer (right)
  geom_segment(aes(x = 3.8, y = 8.5, xend = 4.5, yend = 8.0), linewidth = 0.3, arrow = arrow(length = unit(2, "mm"))) +
  # any_cancer -> excluded
  geom_segment(aes(x = 4.9, y = 7.1, xend = 5.8, yend = 7.9), linewidth = 0.3, arrow = arrow(length = unit(2, "mm"))) +
  # any_cancer -> final (skip level for direct)
  geom_segment(aes(x = 4.5, y = 7.0, xend = 3.3, yend = 5.0), linewidth = 0.3, arrow = arrow(length = unit(2, "mm"))) +
  # excluded -> final (diagonal down-left)
  geom_segment(aes(x = 5.8, y = 7.1, xend = 3.7, yend = 5.0), linewidth = 0.3, arrow = arrow(length = unit(2, "mm"))) +
  # final -> subgroups
  geom_segment(aes(x = 2.2, y = 4.0, xend = 1.2, yend = 3.0), linewidth = 0.3, arrow = arrow(length = unit(2, "mm"))) +
  geom_segment(aes(x = 2.6, y = 4.0, xend = 2.0, yend = 3.0), linewidth = 0.3, arrow = arrow(length = unit(2, "mm"))) +
  geom_segment(aes(x = 3.0, y = 4.0, xend = 3.0, yend = 3.0), linewidth = 0.3, arrow = arrow(length = unit(2, "mm"))) +
  geom_segment(aes(x = 3.4, y = 4.0, xend = 4.0, yend = 3.0), linewidth = 0.3, arrow = arrow(length = unit(2, "mm"))) +
  geom_segment(aes(x = 3.8, y = 4.0, xend = 5.0, yend = 3.0), linewidth = 0.3, arrow = arrow(length = unit(2, "mm"))) +
  coord_fixed(xlim = c(0, 7), ylim = c(1, 10)) +
  theme_void() +
  labs(title = "Figure 1: Study Flow Diagram",
       subtitle = sprintf("NHANES III (1988-1994) and Continuous NHANES (2005-2016), mortality follow-up through 2019"),
       caption = "GI: gastrointestinal. FemaleRepro: cervix, ovary, uterine. ProstateUrinary: prostate, bladder, kidney, testicular. OtherSolid: melanoma, thyroid, brain, bone, head/neck, skin.") +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
    plot.subtitle = element_text(hjust = 0.5, size = 9, color = "grey40"),
    plot.caption = element_text(hjust = 0, size = 7, color = "grey50")
  )

ggsave("figures/flow_diagram.png", p, width = 8, height = 7, dpi = 200)
cat("Flow diagram saved\n")
