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

**Script:** `construct_HFE_fast.sage` 

**What it does:**
Fast, bounded, fail-safe generator of classical HFE instances over GF(2) (targets degree 𝐷, enforces Jacobian rank ≈ 𝑛−1, emits best-so-far on timeout).

**Outputs:** a pipeline `.in` file — variables, field id 2, 𝑛 public equations, then n field equations — plus public and secret logs.

**Example usage**
```bash
# Generate an HFE(n=80, D=96) instance with a fixed seed; write .in and logs
sage scripts/construct_HFE_fast.sage 80 96 data/ --seed 42
# Produces:
#   data/HFE_n80_D96.in
#   logs/HFE_n80_D96_genlog.txt
#   logs/HFE_n80_D96_secret.txt
```

### 2. **System diagnosis**

**Script:** `system_diagnosis.sage`

**What it does:**  
Reads a pipeline `.in` polynomial system, builds the ring/ideal, and reports quick diagnostics (vars/eqs, degrees, sparsity, homogeneity). Optionally computes **Krull dimension** and **degree via Hilbert polynomial** within user-specified wall-time budgets; prints to stdout and saves a log.

**Outputs:** `logs/<basename>_diagnostics.log` (e.g., `logs/HFE_n80_D96_diagnostics.log`), plus the same summary on stdout.

**Example usage**
```bash
# Diagnose an input system with 60s budget for dimension and 120s for Hilbert polynomial
sage scripts/system_diagnosis.sage data/HFE_n80_D96.in --dim=60 --hilbert=120

# Produces:
#   logs/HFE_n80_D96_diagnostics.log
```

### 3. **F4 Gröbner basis (DRL) — streaming parser**

**Script:** `solve_F4_from_file_parsing.jl`

**What it does:**  
Parses a pipeline `.in` file (vars, prime `p`, expanded polynomials), builds a DRL ring over `GF(p)` via `AbstractAlgebra/Nemo`, and runs `AlgebraicSolving.groebner_basis` (F4). Writes a timestamped log and result file, and verifies that the input generators reduce to `0` modulo the computed basis.

**Outputs:**  
- `logs/<basename>_F4_<yyyymmdd_HHMMSS>.log`  
- `results/<basename>_F4_<yyyymmdd_HHMMSS>.txt`

**Example usage**
```bash
# Basic run, 16 threads requested to F4; environment knobs use safe defaults
julia scripts/solve_F4_from_file_parsing.jl data/HFE_n80_D96.in 16

# (Optional) tune F4 via environment variables:
#   F4_THREADS     — threads used inside F4 (default: CLI nthreads)
#   F4_MAX_PAIRS   — cap on S-pairs in the queue (default: 6000)
#   F4_INITIAL_HTS — initial hash table size log2 (default: 19)
#   F4_LA_OPTION   — linear algebra backend/strategy (default: 22)
F4_THREADS=32 F4_MAX_PAIRS=8000 F4_INITIAL_HTS=20 F4_LA_OPTION=22 \
julia scripts/solve_F4_from_file_parsing.jl data/HFE_n80_D96.in 32

# Produces:
#   logs/HFE_n80_D96_F4_<timestamp>.log
#   results/HFE_n80_D96_F4_<timestamp>.txt
```

### 4. **F5 Gröbner basis (DRL) — single-threaded, small systems**

**Script:** `solve_F5_from_file.jl`

**What it does:**  
Reads a pipeline `.in` system (variables, prime `p`, expanded polynomials), builds a DRL ring over `GF(p)`, and computes a Gröbner basis using **F5** (`sig_groebner_basis` from AlgebraicSolving.jl).  
**Single-threaded by design** and intended for **small/medium, sanity-check experiments** — **not** for large instances or long HPC runs.

**Outputs:**  
- `logs/<basename>_F5_<yyyymmdd_HHMMSS>.log`  
- `results/<basename>_F5_<yyyymmdd_HHMMSS>.txt`

**Example usage**
```bash
# Run F5 on a small instance (single-threaded)
julia scripts/solve_F5_from_file.jl data/HFE_n10_D16.in

# Produces:
#   logs/HFE_n10_D16_F5_<timestamp>.log
#   results/HFE_n10_D16_F5_<timestamp>.txt
```

### 5. **DRL → LEX conversion (FGLM with std fallback)**

**Script:** `fglm_adjusted.sage`

**What it does:**  
Converts a **DRL** Gröbner basis (from F4/F5) over `GF(p)` to **LEX**. Tries a **safe FGLM** first, then (on failure/timeout or if preferred) falls back to **Singular’s** Buchberger in **pure LEX** (`algorithm='singular:std'`). Optional post-reduction of the LEX basis via `std`.

**Outputs:**  
- `results/<basename>_LEX.txt` — the (possibly reduced) LEX Gröbner basis  
- `logs/<basename>_FGLM_adjusted.log` — full conversion log (method used, degrees, shape checks)

**Example usage**
```bash
# Preferred: try FGLM (unlimited), fallback to std(lex) if needed; skip dim() check
sage scripts/fglm_adjusted.sage results/HFE_n80_D96_F4_20250101_120000.txt \
  --assume-zerodim --fglm-timeout=0 --reduce=auto --reduce-timeout=300

# Strict budget: allow 60s to confirm dim=0, 1800s for FGLM; always reduce (fail on timeout)
sage scripts/fglm_adjusted.sage results/HFE_n80_D96_F4_20250101_120000.txt \
  --dim-timeout=60 --fglm-timeout=1800 --reduce=always --reduce-timeout=600

# Force direct std(lex), skipping FGLM
sage scripts/fglm_adjusted.sage results/HFE_n80_D96_F4_20250101_120000.txt \
  --prefer=std --fglm-timeout=0 --reduce=never

# CLI flags
--assume-zerodim
    Treat the ideal as zero-dimensional and skip dimension().

--dim-timeout=SECONDS
    Wall-time budget for dimension(); default: 60. Use 0 to skip.

--fglm-timeout=SECONDS
    Wall-time budget for FGLM; default: 0 (unlimited).

--reduce=never|auto|always
    Post-reduce the produced LEX basis via std. Default: never.
      never  -> write the raw LEX basis.
      auto   -> try to reduce; on timeout keep the raw basis.
      always -> enforce reduction; on timeout the job fails.

--reduce-timeout=SECONDS
    Wall-time budget for std reduction; default: 120.

--prefer=auto|fglm|std
    Pipeline policy; default: auto.
      auto -> try FGLM; on failure/timeout fall back to std(lex).
      fglm -> only attempt FGLM (no fallback).
      std  -> skip FGLM and compute std(lex) directly.
```