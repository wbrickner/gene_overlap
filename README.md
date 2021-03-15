# Usage

```
gene_overlap ./genes.csv ./rna.csv ./output.csv
```

# Performance

A naive `O(n^2)` implementation takes 40 minutes to complete an anlysis of ~7,000 genes and ~12,000,000 RNA reads.
This tool analyzes the same dataset in ~5 seconds, and has a time complexity of O(g log r).
