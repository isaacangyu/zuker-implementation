# Zuker Algorithm Implementation for RNA Secondary Structure Prediction

## Motivation

**Formulation:** Given an RNA sequence \( S \in \{A, U, G, C\}^* \), find the non-crossing secondary structure $ P $ with Minimal Free Energy (MFE).

### Advantages over Nussinov Algorithm
- Minimizes free energy instead of maximizing the number of base pairs.
- Uses experimentally determined free energies rather than treating all base pairs equally.

## Recurrence Relations

The optimal RNA secondary structure is computed using a 3D dynamic programming approach.

### General Recurrence (W Matrix)
Let $ W_{i,j} $ represent the MFE of the substring $ S[i..j] $.

$$
W_{i,j} = \min \begin{cases}
W_{i,j-1} & \text{S[j] unpaired} \\
\min_{i \leq k < j - m} W_{i,k-1} + V_{k,j} & \text{otherwise}
\end{cases}
$$

The final MFE is $ W_{1,n} $.

### Paired Recurrence (V Matrix)
Let $ V_{i,j} $ represent the MFE of the substring $ S[i..j] $ where bases $ i $ and $ j $ are paired.

$$
V_{i,j} = \min \begin{cases}
\text{hairpin}(i,j) \\
\text{stacking}(i,j) + V_{i+1,j-1} \\
\min_{i<i'<j'<j} \text{internal}(i,j,i',j') \\
\min_{i<k<j} \text{multiloop}(i+1,k) + \text{multiloop}(k+1,j-1) + a
\end{cases}
$$

- **Hairpin:** Energy for a hairpin loop of size $ j - i - 1 $.
- **Stacking:** Energy for stacking pair plus the MFE of the enclosed region.
- **Internal/Bulge:** Minimum energy over possible internal loops. A bulge loop is a special case of an internal loop where $i'=i+1$ or $j'=j-1$. 
- **Multiloop:** Initialize a split for multiloop

### Multi-loop Recurrence (WM Matrix)
Let $ WM_{i,j} $ represent the MFE of a multi-loop in $ S[i..j] $.

$$
WM_{i,j} = \min \begin{cases}
WM_{i,j-1} + c & \text{$ S[j] $ unpaired} \\
WM_{i+1,j} + c & \text{$ S[i] $ unpaired} \\
V_{i,j} + b & \text{closed} \\
\min_{i < k < j} WM_{i,k} + WM_{k+1,j} & \text{split}
\end{cases}
$$

Where $ a, b, c $ are multiloop parameters (offset, helix penalty, unpaired nucleotide penalty).

## Complexity

- **Time Complexity:** $ O(n^3) $. The algorithm uses three DP tables of size $ O(n^2) $, and filling each cell involves $ O(n) $ operations for finding minimums in the V and WM matrices. By restricting the internal loop size to 30, the computation for internal loops is bounded to $ O(n) $ per cell.
- **Space Complexity:** $ O(n^2) $, due to storing the three DP tables (V, W, WM) of size $ n \times n $.

## Architecture

- `lookup.py`: Handles energy lookup tables for stacking, hairpin, bulge, internal, and multiloop energies. Parses data files and creates JSON caches.
- `stack.json`, `bulge.json`, `hairpin.json`, `internal.json`: JSON files containing precomputed energy values.
- `zuker.py`: Main implementation of the Zuker algorithm, including DP tables (V, W, WM) and backtracing for structure prediction.
- `validation.py`: Scripts for validating the implementation against known structures and datasets.

## Validation

### Dataset: ArchiveII
- Source: [Hugging Face ArchiveII Dataset](https://huggingface.co/datasets/multimolecule/archiveii)
- Provides RNA sequences and corresponding dot-bracket notations for comparison.

Validation was performed with the following steps:

1. Call the Zuker algorithm on all sequences in the dataset and retrieve a list of the predicted base pairs (a list of tuples representing the indices that paired).

2. Represent the base pairs in the actual solution as a list of the predicted base pairs.

3. Compare the two lists and compute the following metrics:
- **F1 Score:** Harmonic mean of precision and recall.
- **Precision:** Fraction of predicted base pairs that are correct.
- **Recall:** Fraction of true base pairs that are predicted.

## Results
We created histograms representing the metrics:
!(zuker-implementation/validation/hist_scores.png)

Precision mean: ___
Precision median: ___

Recall mean: ___
Recall median: ___

F1 Score mean: ___
F1 Score median: ___

## Future Steps

- **Custom Lookup Tables:**
  - Hairpin: Factor in closing nucleotides.
  - Internal: Same as above.
  - Multiloop: Vectorize $ a, b, c $ to account for every inner loop.

## Citations

- Lookup table data files: [RNAstructure] (https://rna.urmc.rochester.edu/Overview/index.html)
- Algorithm: [Zuker Algorithm Slides](https://math.mit.edu/classes/18.417/Slides/rna-prediction-zuker.pdf)
- Dataset: [ArchiveII on Hugging Face](https://huggingface.co/datasets/multimolecule/archiveii)