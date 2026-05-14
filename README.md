# Hidden Markov Model & Viterbi Algorithm — Gene Finding

A computational biology assignment implementing a **Hidden Markov Model (HMM)** to identify gene structure (exons, introns, and splice sites) in a DNA sequence, using the **Viterbi algorithm** for optimal state decoding.

---

## Problem Description

Given a DNA query sequence, the goal is to determine the most likely annotation of each nucleotide as belonging to an **Exon**, **Intron**, or **5' splice site** using an HMM with five hidden states:

| State | Symbol | Description |
|-------|--------|-------------|
| Start | `s` | Initial state (no emission) |
| Exon | `E` | Protein-coding region |
| 5' Splice Site | `5` | Boundary between exon and intron |
| Intron | `I` | Non-coding region |
| End | `e` | Terminal state (no emission) |

### Query Sequence

```
CTTCATGTGAAAGCAGACGTAAGTCA
```

---

## HMM Parameters

### State Transition Probabilities

|       | s   | E   | 5   | I   | e   |
|-------|-----|-----|-----|-----|-----|
| **s** | 0.0 | 1.0 | 0.0 | 0.0 | 0.0 |
| **E** | 0.0 | 0.9 | 0.1 | 0.0 | 0.0 |
| **5** | 0.0 | 0.0 | 0.0 | 1.0 | 0.0 |
| **I** | 0.0 | 0.0 | 0.0 | 0.9 | 0.1 |
| **e** | 0.0 | 0.0 | 0.0 | 0.0 | 0.0 |

### Emission Probabilities

|       | A    | C    | G    | T    |
|-------|------|------|------|------|
| **s** | 0.00 | 0.00 | 0.00 | 0.00 |
| **E** | 0.25 | 0.25 | 0.25 | 0.25 |
| **5** | 0.05 | 0.00 | 0.95 | 0.00 |
| **I** | 0.40 | 0.10 | 0.10 | 0.40 |
| **e** | 0.00 | 0.00 | 0.00 | 0.00 |

---

## Implementation Details

The notebook is organized into three main stages:

### 1. Manual Path Evaluation

Several candidate state paths are manually constructed and scored using `get_log_prob_for_state_path()` to demonstrate how different exon–intron boundaries affect the overall log-probability:

| Path | Splice Position | Log Probability |
|------|----------------|-----------------|
| `EEEEEE5IIIIIIIIIIIIIIIIIII` | 7 | -43.90 |
| `EEEEEEEE5IIIIIIIIIIIIIIIII` | 9 | -43.45 |
| `EEEEEEEEEEEE5IIIIIIIIIIIII` | 13 | -43.94 |
| `EEEEEEEEEEEEEEE5IIIIIIIIII` | 16 | -42.58 |
| `EEEEEEEEEEEEEEEEEE5IIIIIII` | 19 | -41.22 |
| `EEEEEEEEEEEEEEEEEEEEEE5III` | 23 | -41.71 |
| `EEEEEEEEEEEEEEEEEEEEEEEEEE` | None (all exon) | -40.98 |

### 2. Viterbi Algorithm (Dynamic Programming)

The algorithm efficiently finds the optimal state path using two matrices:

- **`viterbi_value_matrix`** — Stores the maximum log-probability of reaching each state at each position
- **`viterbi_trace_matrix`** — Stores backpointers to the previous best state for traceback

**Key function: `calculate_prob_for_a_node(current_state, t)`**

For each cell `(state, time_step)`, the function:
1. Checks if the current state can emit the observed nucleotide
2. Iterates over all possible previous states
3. Computes: `V[prev][t-1] + log(transition) + log(emission)`
4. Returns the maximum value and the index of the best previous state

### 3. Traceback

Starting from the highest-probability state in the last column, the algorithm traces back through `viterbi_trace_matrix` to reconstruct the optimal hidden state path.

---

## Results

```
Query:      CTTCATGTGAAAGCAGACGTAAGTCA
Best path:  EEEEEEEEEEEEEEEEEEEEEEEEEE
Log prob:   -38.6777
```

The Viterbi algorithm determines that the most probable decoding assigns all nucleotides to the **Exon (E)** state for this particular sequence and model parameterization.

---

## Requirements

- Python 3.x
- NumPy

```bash
pip install numpy
```

---

## Usage

1. Clone the repository:
   ```bash
   git clone https://github.com/<your-username>/HMM-Viterbi.git
   cd HMM-Viterbi
   ```

2. Open and run the Jupyter Notebook:
   ```bash
   jupyter notebook HMM_Viterbi.ipynb
   ```

3. All cells can be executed sequentially — no external data files are required.

---

## Project Structure

```
HMM-Viterbi/
├── HMM_Viterbi.ipynb   # Main notebook with HMM setup, Viterbi algorithm, and results
└── README.md            # Project documentation
```

---

## References

- Durbin, R., Eddy, S., Krogh, A., & Mitchison, G. (1998). *Biological Sequence Analysis: Probabilistic Models of Proteins and Nucleic Acids*. Cambridge University Press.
- Viterbi, A. (1967). Error bounds for convolutional codes and an asymptotically optimum decoding algorithm. *IEEE Transactions on Information Theory*, 13(2), 260–269.
