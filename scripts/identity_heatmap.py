

from Bio import SeqIO, pairwise2
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

fasta = "genome_dataset.fasta"

records = list(SeqIO.parse(fasta, "fasta"))

names = [r.id for r in records]

matrix = []

for i in range(len(records)):
    row = []
    for j in range(len(records)):
        seq1 = str(records[i].seq)
        seq2 = str(records[j].seq)

        align = pairwise2.align.globalxx(seq1, seq2, one_alignment_only=True)[0]

        matches = align.score
        length = max(len(seq1), len(seq2))

        identity = (matches / length) * 100

        row.append(identity)

    matrix.append(row)

df = pd.DataFrame(matrix, index=names, columns=names)

df.to_csv("identity_matrix.csv")

plt.figure(figsize=(12,10))
sns.heatmap(df, cmap="viridis")

plt.title("Genome Identity Heatmap")
plt.tight_layout()

plt.savefig("identity_heatmap.png", dpi=300)

