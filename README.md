# Contraction Hierarchies – Shortest Path Engine

A Python implementation of the **Contraction Hierarchies (CH)** algorithm for computing shortest paths in large directed graphs. CH is a widely used speed-up technique in real-world routing engines such as OSRM and GraphHopper.

## What is Contraction Hierarchies?

Standard Dijkstra's algorithm becomes slow on large graphs (e.g. city road networks with millions of nodes). Contraction Hierarchies solve this with two phases:

1. **Preprocessing**: Nodes are contracted one by one in order of importance. When a node is removed, shortcut edges are added to preserve shortest path distances.
2. **Query**: A bidirectional Dijkstra search runs on the augmented graph, only expanding nodes of increasing rank — dramatically reducing the search space.

## Features

- DIMACS graph file parser (directed & undirected)
- Standard Dijkstra (reference implementation)
- Distance-limited Dijkstra (used during contraction)
- Full node contraction with shortcut detection
- Bidirectional CH-Dijkstra with shortcut unpacking
- Benchmark script with timing and expanded-node statistics

## Results on `rome99.gr` (Rome Road Network)

| Metric | Value |
|---|---|
| Nodes | 3,353 |
| Original edges | 8,859 |
| Shortcuts inserted | ~10,000+ |
| Processing time | < 10s |
| Avg expanded nodes – Dijkstra | ~1,600 |
| Avg expanded nodes – CH Bidir | ~130 |
| **Speedup** | **~15x** |

> Results may vary slightly depending on hardware.

## Getting Started

**Requirements:** Python 3.10+, no external libraries needed.

```bash
git clone https://github.com/YOUR_USERNAME/contraction-hierarchies.git
cd contraction-hierarchies
```

**Run the benchmark** (requires `rome99.gr` from [DIMACS](http://www.diag.uniroma1.it/challenge9/download.shtml)):

```bash
python ch_framework.py
```

**Use in your own code:**

```python
from ch_framework import read_file, contraction, bidirectional_dijkstra

graph = read_file('rome99.gr', directed=True)
g_ch, rank = contraction(graph)

distance, expanded, path = bidirectional_dijkstra(g_ch, rank, s=1, t=100)
print(f"Shortest distance: {distance}")
print(f"Path: {path}")
```

## Project Structure

```
contraction-hierarchies/
├── ch_framework.py   # Full implementation
└── README.md
```

## How It Works

### Preprocessing
Each node `v` is contracted by:
1. Finding all pairs `(u, w)` where `u → v → w`
2. Running a limited Dijkstra from `u` (excluding `v`) to check if a shorter path exists
3. If not, inserting a shortcut edge `u → w`

Nodes are contracted in order of minimum degree (in + out), breaking ties by node ID.

### Query
A bidirectional Dijkstra runs simultaneously from `s` (forward, upward in rank) and `t` (backward, upward in rank). The search stops when `min(L_up, L_down) > best_distance_found`. Shortcut edges in the result path are recursively unpacked to return the true node sequence.

## References

- Geisberger et al. (2008) – *Contraction Hierarchies: Faster and Simpler Hierarchical Routing in Road Networks*
- [DIMACS Shortest Path Challenge](http://www.diag.uniroma1.it/challenge9/download.shtml)
