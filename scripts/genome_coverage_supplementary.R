library(ggplot2)
library(patchwork)

depth <- read.table("depth.txt")

colnames(depth) <- c("genome","position","depth")

# gráfico de cobertura ao longo do genoma
p1 <- ggplot(depth, aes(x=position, y=depth)) +
  geom_line(color="black") +
  theme_classic() +
  labs(x="Genome position (bp)",
       y="Sequencing depth")

# histograma de profundidade
p2 <- ggplot(depth, aes(x=depth)) +
  geom_histogram(bins=50, fill="grey70", color="black") +
  theme_classic() +
  labs(x="Sequencing depth",
       y="Frequency")

fig <- p1 / p2

ggsave("genome_coverage_supplementary.svg",
       fig,
       width=7,
       height=6)
