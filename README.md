# Gröbner Bases Applied to Algebraic Cryptography — Code & Experiments

This repository contains the complete code accompanying the Master’s thesis
**“Gröbner Bases Applied to Algebraic Cryptography”** by Frederick Spreafichi (University of Trieste, MSc in Mathematics).

The project investigates the effectiveness and role of Gröbner bases in algebraic cryptography, with experiments on classic HFE over GF(2) and on ANEMOI in P_CICO mode (using the ANEMOI generation code from the official repository: https://github.com/isec-tugraz/six-worlds-anemoi.git).

1. computation of a DRL Gröbner basis via **F4/F5**;

2. change of monomial order (DRL → LEX) using **FGLM**;

3. root finding (typically, factorization of a univariate polynomial in the LEX basis and back-substitution).

The codebase uses Julia, SageMath, and Singular to run Gröbner-basis attacks on the polynomial systems arising from algebraic models of HFE and ANEMOI. It supports both local execution via scripts and HPC execution via Slurm (sbatch) job files. The original experimental results reported in the thesis were obtained on Demetra, the HPC cluster of the Department of Mathematics, Informatics and Geosciences (DMIG), University of Trieste.

---

## Pipeline overview

1. **HFE instance generation (SageMath)**
Generate HFE(𝑛,𝐷) instances over GF(2), with optional planted solution and field equations.
Output: data/hfe_instances/HFE_n{n}_D{D}.in (+ metadata log).

2. **System diagnosis (SageMath)**
Sanity checks and quick invariants (variable/eq counts, degrees, sparsity, optional random evaluation tests).
Output: summary report.

3. **F4 Gröbner basis (AlgebraicSolving.jl)**
Compute a DRL Gröbner basis using the linear-algebra (Macaulay) approach with configurable pairing/degree buckets.
Output: DRL GB, per-degree logs (matrix sizes, density, new reducers), timings.

4. **F5 / signature-based (AlgebraicSolving.jl)**
Optional runs with signature criteria (syzygy/rewritable) to compare against F4 on the same instances.
Output: DRL GB, zero-reduction statistics, timings.

5. **FGLM (SageMath backend: Singular)**
Change of ordering DRL → LEX, compute multiplication matrices, standard monomials, and (when in shape) the univariate.
Output: reduced LEX GB.

6. **Solution retrieval (SageMath)**
Factor the univariate, back-substitute, and verify solutions against the original system.
Output: solutions.

Execution modes:
Local: run the scripts in sequence on a workstation/laptop.
HPC (Slurm): submit the provided sbatch files (single-node, threaded) for long runs and automatic logging.

---

### 1. **HFE instance generation**

**Script:** `scripts/construct_HFE_fast.sage` 

**What it does:**
Fast, bounded, fail-safe generator of classical HFE instances over GF(2) (targets degree 𝐷, enforces Jacobian rank ≈ 𝑛−1, emits best-so-far on timeout).

**Outputs:** a pipeline .in file — variables, field id 2, 𝑛 public equations, then n field equations — plus public and secret logs.

**Example usage**
```bash
# Generate an HFE(n=80, D=96) instance with a fixed seed; write .in and logs
sage scripts/construct_HFE_fast.sage 80 96 data/ --seed 42
# Produces:
#   data/HFE_n80_D96.in
#   logs/HFE_n80_D96_genlog.txt
#   logs/HFE_n80_D96_secret.txt
```




