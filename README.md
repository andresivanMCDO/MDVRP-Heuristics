# Multi-Depot Vehicle Routing Problem — Heuristic Framework

A self-contained **Julia** implementation of the three-stage MDVRP heuristic framework proposed by Geetha, Vanathi & Poonthalir (2012), validated against 25 benchmark instances from the literature. The routing stage compares the **Nearest Neighbour Heuristic (NNH)** against an exact **MILP** solver (HiGHS via JuMP) to quantify the cost of using a heuristic instead of an optimal router.

---

## Table of Contents

- [Problem Description](#problem-description)
- [Solution Framework](#solution-framework)
  - [Clustering](#1-clustering)
  - [Grouping](#2-grouping)
  - [Routing](#3-routing)
- [Instances](#instances)
- [Results](#results)
  - [Route Comparisons](#route-comparisons)
  - [Cost Comparison](#cost-comparison)
  - [Computational Time](#computational-time)
- [Implementation](#implementation)
- [References](#references)

---

## Problem Description

The **Multi-Depot Vehicle Routing Problem (MDVRP)** is an NP-hard combinatorial optimization problem. Given a graph $G = (V, E)$ with a set of customers $V_c = \{v_1, \ldots, v_n\}$ and depots $V_d = \{v_{n+1}, \ldots, v_{n+t}\}$, the goal is to construct a set of minimum-cost vehicle routes such that:

- Every customer is visited exactly once.
- Each route starts and ends at the same depot.
- The total demand on each route does not exceed the vehicle capacity $Q$.
- The total duration of each route does not exceed the bound $D$.

The objective minimizes total travel cost across all vehicles at all depots:

$$\min Z = \sum_{p=1}^{t} \sum_{q=1}^{k_p} \sum_{i=1}^{n} \sum_{j=1}^{n} c_{ij}\, x_{ijqp}$$

where $x_{ijqp} = 1$ if vehicle $q$ at depot $p$ travels from $v_i$ to $v_j$, and $c_{ij}$ is the Euclidean distance between them.

---

## Solution Framework

The framework operates in three sequential stages as described in Geetha et al. (2012).

### 1. Clustering

Each customer is assigned to its nearest depot by Euclidean distance, subject to the depot's aggregate capacity $Q_d = k_s \cdot Q$. If the nearest depot is full, the customer falls back to the next nearest depot with available capacity.

| Instance group | Clustering |
|---|---|
| p01–p07 | ![p01–p07 clusters](Clustering/p0107Clusters.png) |
| p08–p11 | ![p08–p11 clusters](Clustering/p0811Clusters.png) |
| p12–p23 | ![p12–p23 clusters](Clustering/p1223Clusters.png) |
| pr01–pr10 | ![pr01–pr10 clusters](Clustering/pr0110Clusters.png) |
| pfbo | ![pfbo clusters](Clustering/pfboClusters.png) |

### 2. Grouping

Within each cluster, customers are partitioned among the depot's $k_j$ vehicles using an **angular sweep** (Nurcahyo et al., 2003). The depot is treated as the origin of a polar coordinate system. A ray sweeps around it in angular order; when adding the next customer would exceed capacity $Q$, a new vehicle slice begins.

| Instance group | Grouping |
|---|---|
| p01–p07 | ![p01–p07 grouping](Grouping/p0107Grouping.png) |
| p08–p11 | ![p08–p11 grouping](Grouping/p0811Grouping.png) |
| p12–p23 | ![p12–p23 grouping](Grouping/p1223Grouping.png) |
| pr01–pr10 | ![pr01–pr10 grouping](Grouping/pr0110Grouping.png) |
| pfbo | ![pfbo grouping](Grouping/pfboGrouping.png) |

### 3. Routing

Each vehicle's assignment is a TSP on $n_v + 1$ nodes (customers + home depot). Two methods are compared:

**Nearest Neighbour Heuristic (NNH)** — from the depot, always move to the closest unvisited customer, then return. Implemented via `TravelingSalesmanHeuristics.jl`.

**MILP (exact TSP)** — Miller–Tucker–Zemlin formulation solved with HiGHS via JuMP:
- Binary arc variables $x_{ij} \in \{0,1\}$.
- MTZ subtour elimination: $u_i - u_j + n \cdot x_{ij} \leq n - 1$ for $i, j \geq 2$.

---

## Instances

The benchmark set includes 33 MDVRP instances from four sources. 25 were tested in this work.

| Instances | Source |
|-----------|--------|
| p01–p07 | Christofides & Eilon (1969) |
| p08–p11 | Gillett & Johnson (1976) |
| p12–p23 | Chao, Golden & Wasil (1993) |
| pr01–pr10 | Cordeau, Gendreau & Laporte (1997) |

Marker color in the plots below encodes customer demand — darker means higher demand.

| Instance group | Visualisation |
|---|---|
| p01–p07 | ![p01–p07 instances](Instances/p0107Instances.png) |
| p08–p11 | ![p08–p11 instances](Instances/p0811Instances.png) |
| p12–p23 | ![p12–p23 instances](Instances/p1223Instances.png) |
| pr01–pr10 | ![pr01–pr10 instances](Instances/pr0110Instances.png) |
| pfbo | ![pfbo instance](Instances/pfboInstance.png) |

---

## Results

### Route Comparisons

For p01–p07, p12–p23, and pr01–pr10 each plot shows NNH on the left and MILP on the right. Instances p08–p11 are NNH only — their per-vehicle subproblems are large enough to make the exact solver computationally prohibitive.

<details>
<summary><strong>p01–p07 (click to expand)</strong></summary>

| | |
|---|---|
| ![p01](Routing/p01_routing.png) | ![p02](Routing/p02_routing.png) |
| p01 — NNH vs MILP | p02 — NNH vs MILP |
| ![p03](Routing/p03_routing.png) | ![p04](Routing/p04_routing.png) |
| p03 — NNH vs MILP | p04 — NNH vs MILP |
| ![p05](Routing/p05_routing.png) | ![p06](Routing/p06_routing.png) |
| p05 — NNH vs MILP | p06 — NNH vs MILP |
| ![p07](Routing/p07_routing.png) | |
| p07 — NNH vs MILP | |

</details>

<details>
<summary><strong>p08–p11 — NNH only (click to expand)</strong></summary>

| | |
|---|---|
| ![p08](Routing/p08_routing.png) | ![p09](Routing/p09_routing.png) |
| p08 — NNH | p09 — NNH |
| ![p10](Routing/p10_routing.png) | ![p11](Routing/p11_routing.png) |
| p10 — NNH | p11 — NNH |

</details>

<details>
<summary><strong>p12–p23 (click to expand)</strong></summary>

| | |
|---|---|
| ![p12](Routing/p12_routing.png) | ![p13](Routing/p13_routing.png) |
| p12 — NNH vs MILP | p13 — NNH vs MILP |
| ![p14](Routing/p14_routing.png) | ![p15](Routing/p15_routing.png) |
| p14 — NNH vs MILP | p15 — NNH vs MILP |
| ![p16](Routing/p16_routing.png) | ![p17](Routing/p17_routing.png) |
| p16 — NNH vs MILP | p17 — NNH vs MILP |
| ![p18](Routing/p18_routing.png) | ![p19](Routing/p19_routing.png) |
| p18 — NNH vs MILP | p19 — NNH vs MILP |
| ![p20](Routing/p20_routing.png) | ![p21](Routing/p21_routing.png) |
| p20 — NNH vs MILP | p21 — NNH vs MILP |
| ![p22](Routing/p22_routing.png) | ![p23](Routing/p23_routing.png) |
| p22 — NNH vs MILP | p23 — NNH vs MILP |

</details>

<details>
<summary><strong>pr01–pr10 and pfbo (click to expand)</strong></summary>

| | |
|---|---|
| ![pr01](Routing/pr01_routing.png) | ![pr02](Routing/pr02_routing.png) |
| pr01 — NNH vs MILP | pr02 — NNH vs MILP |
| ![pr03](Routing/pr03_routing.png) | ![pr04](Routing/pr04_routing.png) |
| pr03 — NNH vs MILP | pr04 — NNH vs MILP |
| ![pr05](Routing/pr05_routing.png) | ![pr06](Routing/pr06_routing.png) |
| pr05 — NNH vs MILP | pr06 — NNH vs MILP |
| ![pr07](Routing/pr07_routing.png) | ![pr08](Routing/pr08_routing.png) |
| pr07 — NNH vs MILP | pr08 — NNH vs MILP |
| ![pr09](Routing/pr09_routing.png) | ![pr10](Routing/pr10_routing.png) |
| pr09 — NNH vs MILP | pr10 — NNH vs MILP |
| ![pfbo](Routing/pfbo_routing.png) | |
| pfbo — NNH vs MILP | |

</details>

---

### Cost Comparison

`—` indicates MILP was not run for that instance group.

| Instance | NNH Cost | MILP Cost | Gap (%) | Mean gap/vehicle (%) |
|----------|--------:|--------:|--------:|--------------------:|
| p01 | 690.37 | 687.26 | 0.45 | 0.32 |
| p02 | 541.11 | 540.76 | 0.07 | 0.05 |
| p03 | 767.97 | 755.27 | 1.68 | 1.21 |
| p04 | 1206.93 | 1206.78 | 0.01 | 0.01 |
| p05 | 826.83 | 810.89 | 1.97 | 2.14 |
| p06 | 959.46 | 957.08 | 0.25 | 0.22 |
| p07 | 1129.75 | 1128.36 | 0.12 | 0.18 |
| p08 | 5112.18 | — | — | — |
| p09 | 4609.38 | — | — | — |
| p10 | 4573.16 | — | — | — |
| p11 | 4562.21 | — | — | — |
| p12 | 1365.69 | 1365.69 | 0.00 | 0.00 |
| p13 | 1365.69 | 1365.69 | 0.00 | 0.00 |
| p14 | 1365.69 | 1365.69 | 0.00 | 0.00 |
| p15 | 2731.37 | 2731.37 | 0.00 | 0.00 |
| p16 | 2731.37 | 2731.37 | 0.00 | 0.00 |
| p17 | 2731.37 | 2731.37 | 0.00 | 0.00 |
| p18 | 4097.06 | 4097.06 | 0.00 | 0.00 |
| p19 | 4097.06 | 4097.06 | 0.00 | 0.00 |
| p20 | 4097.06 | 4097.06 | 0.00 | 0.00 |
| p21 | 6145.59 | 6145.59 | 0.00 | 0.00 |
| p22 | 6145.59 | 6145.59 | 0.00 | 0.00 |
| p23 | 6145.59 | 6145.59 | 0.00 | 0.00 |
| pr01 | 904.48 | 891.54 | 1.45 | 0.85 |
| pr02 | 1671.89 | 1663.39 | 0.51 | 0.48 |
| pr03 | 1985.72 | 1968.90 | 0.85 | 0.82 |
| pr04 | 2512.70 | 2479.01 | 1.36 | 1.22 |
| pr05 | 3299.10 | 3279.50 | 0.60 | 0.53 |
| pr06 | 3537.08 | 3515.60 | 0.61 | 0.51 |

**Key observations:**
- Instances p12–p23 have a regular grid-like structure around each depot — NNH is already optimal (0% gap) for all of them.
- The maximum overall gap across all instances with MILP results is **1.97%** (p05).
- On p05, only 2 of 9 vehicles show any divergence, both under 1.5% per vehicle, demonstrating the gap is highly localized.

---

### Computational Time

The efficiency ratio is $t_{\text{MILP}} / t_{\text{NNH}}$.

| Instance | NNH (s) | MILP (s) | Ratio |
|----------|--------:|---------:|------:|
| p01 | 0.0001 | 0.2885 | 4,373× |
| p02 | 0.0001 | 0.6728 | 10,385× |
| p03 | 0.0001 | 0.7028 | 8,877× |
| p04 | 0.0001 | 2.6327 | 22,009× |
| p05 | 0.0001 | 9.5797 | 116,460× |
| p06 | 0.0001 | 1.6316 | 15,991× |
| p07 | 0.0001 | 2.1024 | 25,266× |
| p08 | 0.000275 | — | — |
| p09 | 0.000278 | — | — |
| p10 | 0.000259 | — | — |
| p11 | 0.000272 | — | — |
| p12 | 0.0001 | 0.4636 | 6,567× |
| p17 | 0.0001 | 0.9403 | 10,978× |
| p21 | 0.0002 | 2.0762 | 13,290× |
| p23 | 0.0002 | 2.0856 | 12,097× |
| pr01 | 0.0001 | 1.1443 | 19,195× |
| pr02 | 0.0001 | 5.4612 | 56,513× |
| pr03 | 0.0001 | 7.8457 | 64,937× |
| pr04 | 0.0002 | 36.6888 | 234,613× |
| pr05 | 0.0002 | 51.2414 | 277,306× |
| pr06 | 0.0002 | 63.1172 | 350,291× |

NNH is **4,000× to 350,000× faster** than the MILP solver. MILP runtime scales steeply with instance size — pr06 takes over a minute; NNH solves it in under a millisecond. NNH solve time rises slightly from ~0.0001 s on p01–p07 to ~0.00027 s on p08–p11, reflecting the larger number of customers per vehicle in that group.

---

## Implementation

The codebase is organized into two files:

**`MDVRP_Module.jl`** — core module

| Function | Description |
|----------|-------------|
| `parse_mdvrp(filepath)` | Reads an instance file; returns an `MDVRPInstance` struct with customer/depot data, coordinates, and the full $C$ matrix |
| `instance_to_dataframe(inst)` | Converts the instance to a `DataFrame` for plotting |
| `create_static_df(inst)` | Builds a customer-to-depot distance ranking table |
| `create_dynamic_df(inst, static_df, df)` | Runs clustering; returns each customer's assigned depot |
| `create_grouping_df(inst, dynamic_df)` | Runs the angular sweep; returns each customer's assigned vehicle |
| `plot_instance / plot_clustering / plot_groups` | Visualization functions |

**`MDVRP.ipynb`** — orchestration notebook: loads instances, calls module functions, runs NNH and MILP routing, produces tables and figures.

**Dependencies:**

```julia
using Distances, DataFrames, Plots, Colors, Random
using JuMP, HiGHS, TravelingSalesmanHeuristics
```

---

## References

1. N. Christofides and S. Eilon, "An algorithm for the vehicle-dispatching problem," *Oper. Res. Q.*, vol. 20, no. 3, pp. 309–318, 1969.
2. B. Gillett and J. Johnson, "Multi-terminal vehicle-dispatch algorithm," *Omega*, vol. 4, no. 6, pp. 711–718, 1976.
3. I. Chao, B. Golden, and E. Wasil, "A new heuristic for the multi-depot vehicle routing problem that improves upon best-known solutions," *Am. J. Math. Manag. Sci.*, vol. 13, no. 3, pp. 371–406, 1993.
4. J.-F. Cordeau, M. Gendreau, and G. Laporte, "A tabu search heuristic for periodic and multi-depot vehicle routing problems," *Networks*, vol. 30, no. 2, pp. 105–119, 1997.
5. S. Geetha, P. T. Vanathi, and G. Poonthalir, "Metaheuristic approach for the multi-depot vehicle routing problem," *Applied Artificial Intelligence*, vol. 26, no. 9, pp. 878–901, 2012. https://doi.org/10.1080/08839514.2012.727344
6. J. K. Lenstra and A. H. G. Rinnooy Kan, "Complexity of vehicle routing and scheduling problems," *Networks*, vol. 11, no. 2, pp. 221–227, 1981.
7. G. W. Nurcahyo, R. A. Alias, S. Shamsuddin, and M. Sap, "Sweep algorithm in vehicle routing problem for public transport," *Jurnal Antarabangsa Teknologi Maklumat*, vol. 2, pp. 51–64, 2003.
8. T. Vidal, T. G. Crainic, M. Gendreau, N. Lahrichi, and W. Rei, "A hybrid genetic algorithm for multi-depot and periodic vehicle routing problems," *Oper. Res.*, vol. 60, no. 3, pp. 611–624, 2012.
