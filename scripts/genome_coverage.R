

library(ggplot2)

# ler fasta
lines <- readLines("aligned.fasta")

# identificar headers
headers <- grep("^>", lines)

# sequência 1
seq1 <- paste(lines[(headers[1]+1):(headers[2]-1)], collapse="")

# sequência 2
seq2 <- paste(lines[(headers[2]+1):length(lines)], collapse="")

seq1 <- toupper(seq1)
seq2 <- toupper(seq2)

s1 <- strsplit(seq1,"")[[1]]
s2 <- strsplit(seq2,"")[[1]]

# parâmetros sliding window
window_size <- 200
step_size <- 20

positions <- seq(1, length(s1) - window_size, by = step_size)

identity <- numeric(length(positions))

for(i in seq_along(positions)){
  
  start <- positions[i]
  end <- start + window_size - 1
  
  w1 <- s1[start:end]
  w2 <- s2[start:end]
  
  identity[i] <- sum(w1 == w2) / window_size * 100
}

data <- data.frame(
  position = positions,
  identity = identity
)

# ORFs (coordenadas típicas de hepadnavírus)
orfs <- data.frame(
  gene = c("Pol","Surface","Precore","Core","X"),
  start = c(2132,1,1618,1714,1211),
  end   = c(1466,672,1713,2280,1618)
)

genome_length <- length(s1)

orfs_expanded <- data.frame()

for(i in 1:nrow(orfs)){
  
  if(orfs$start[i] > orfs$end[i]){
    
    part1 <- data.frame(
      gene = orfs$gene[i],
      start = orfs$start[i],
      end = genome_length
    )
    
    part2 <- data.frame(
      gene = orfs$gene[i],
      start = 1,
      end = orfs$end[i]
    )
    
    orfs_expanded <- rbind(orfs_expanded, part1, part2)
    
  } else {
    
    orfs_expanded <- rbind(orfs_expanded, orfs[i,])
    
  }
}

# gráfico
p <- ggplot() +
  
  geom_line(
    data = data,
    aes(x = position, y = identity),
    linewidth = 1
  ) +
  
  geom_rect(
    data = orfs_expanded,
    aes(
      xmin = start,
      xmax = end,
      ymin = -5,
      ymax = 0,
      fill = gene
    ),
    color = "black"
  ) +
  
  labs(
    x = "Genome position (nt)",
    y = "Sequence identity (%)"
  ) +
  
  ylim(-5,100) +
  
  theme_classic()

print(p)

# salvar figura
ggsave(
  "identity_orf_plot.svg",
  p,
  width = 8,
  height = 4
)

