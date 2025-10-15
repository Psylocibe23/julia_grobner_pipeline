# Gröbner Bases Applied to Algebraic Cryptography — Code & Experiments

This repository contains the complete code accompanying the Master’s thesis
**“Gröbner Bases Applied to Algebraic Cryptography”** by Frederick Spreafichi (University of Trieste, MSc in Mathematics).

The project investigates the effectiveness and role of Gröbner bases in algebraic cryptography, with experiments on classic HFE over `GF(2)` and on ANEMOI in P_CICO mode (using the ANEMOI generation code from the official repository: https://github.com/isec-tugraz/six-worlds-anemoi.git).

1. computation of a DRL Gröbner basis via **F4/F5**;

2. change of monomial order (DRL → LEX) using **FGLM**;

3. root finding (typically, factorization of a univariate polynomial in the LEX basis and back-substitution).

The codebase uses Julia, SageMath, and Singular to run Gröbner-basis attacks on the polynomial systems arising from algebraic models of HFE and ANEMOI. It supports both local execution via scripts and HPC execution via Slurm (sbatch) job files. The original experimental results reported in the thesis were obtained on Demetra, the HPC cluster of the Department of Mathematics, Informatics and Geosciences (DMIG), University of Trieste.

---

## Pipeline overview

1. **HFE instance generation (SageMath)**
Generate HFE(𝑛,𝐷) instances over `GF(2)`, with optional planted solution and field equations.
Output: data/hfe_instances/HFE_n_D.in (+ metadata log).

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
Fast, bounded, fail-safe generator of classical HFE instances over `GF(2)` (targets degree 𝐷, enforces Jacobian rank ≈ 𝑛−1, emits best-so-far on timeout).

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

### 6. **Solution retrieval from LEX basis**

**Script:** `extract_solutions_from_lex.sage`

**What it does:**  
Given a Gröbner basis in **LEX** order over GF(p), extracts all solutions in the base field. It prefers a fast triangular route: solve a univariate in the last variable if present; otherwise branch over the last variable’s values (cheap over GF(2)) and ascend with pruning. Verifies and deduplicates all candidates; includes an optional small-n fallback via `I.variety(GF(p))`.

**Outputs:**  
- `results/<stem>_LEX_sols.txt` — one solution per line (e.g., `{x0: 1, x1: 0, ...}`)  
- `logs/<stem>_LEX_extract.log` — verbose run log

**Example usage**
```bash
# Fast triangular extraction, keep default safeguards
sage scripts/extract_solutions_from_lex.sage results/HFE_n80_D96_F4_20250101_120000_LEX.txt

# Return only the first solution, cap exploration nodes, allow elimination fallback
sage scripts/extract_solutions_from_lex.sage results/HFE_n10_D16_F4_20250101_120000_LEX.txt \
  --first-only --node-limit 200000 --elim

# Enumerate at most 8 solutions, disable small-n fallback
sage scripts/extract_solutions_from_lex.sage results/HFE_n24_D48_F4_20250101_120000_LEX.txt \
  --max-solutions 8 --no-fallback
```

### 7. **HFE solution validity checker**

**Script:** `test_hfe_solution_validity.sage`

**What it does:**  
Given a pipeline `.in` file (variables, field, public equations, optional field equations), a solutions file (one dict-like solution per line), and optionally a generation/secret log, verifies for each candidate:  
(1) satisfies all public equations;  
(2) satisfies field equations `x_i^2 + x_i` if present;  
(3) equals the planted secret if available;  
(4) S-equivalence to the secret via `A_S x + b_S` over `GF(2)`. Optionally prints the Jacobian rank at each solution.

**Outputs:**  
- Prints a per-solution validation report to stdout  
- No files written by default (intended for quick checks in the pipeline)

**Example usage**
```bash
# Basic check: public + (optional) field equations
sage scripts/test_hfe_solution_validity.sage data/HFE_n80_D96.in results/HFE_n80_D96_F4_20250101_120000_LEX_sols.txt

# Also load secret & affine maps from generation log, and print Jacobian rank
sage scripts/test_hfe_solution_validity.sage data/HFE_n80_D96.in \
  results/HFE_n80_D96_F4_20250101_120000_LEX_sols.txt logs/HFE_n80_D96_secret.txt --rank

# Provide a secret explicitly on the CLI (comma/space-separated), no log file
sage scripts/test_hfe_solution_validity.sage data/HFE_n10_D16.in \
  results/HFE_n10_D16_F4_20250101_120000_LEX_sols.txt --secret 1,0,1,0,0,1 --rank
```

### 8. **ANEMOI P_CICO emitter (l = 1)**

**Script:** `emit_anemoi_pcico.sage`

**What it does:**  
Emits **ANEMOI** polynomial systems in **P_CICO** mode with **l = 1** over **GF(p)** (odd prime). You specify the prime `p`, S-box exponent `alpha` (odd, `gcd(alpha, p−1)=1`), and the number of rounds `r`. The script builds the model (ordering 1|2|3), and writes a pipeline `.in` file ready for the F4/F5 → FGLM → solution pipeline.

**Outputs:**  
- `data/ANEMOI_p<p>_r<r>_a<alpha>_PCICO[ _noFLL ].in` (unless `--outfile` is provided)  
  - Line 1: comma-separated variable names (in the DRL tie-break order used by the model)  
  - Line 2: the prime `p`  
  - Line 3+: one polynomial per line (Sage textual form; coefficients in `[0..p-1]`)

**Example usage**
```bash
# Minimal: prime by value, odd alpha, r rounds
sage scripts/emit_anemoi_pcico.sage --p 251 --alpha 3 --rounds 2

# Choose variable ordering and omit the final linear layer constraint
sage scripts/emit_anemoi_pcico.sage --p 101 --alpha 5 --rounds 3 --ordering 2 --no-final-ll

# Hex literal for p and explicit output path
sage scripts/emit_anemoi_pcico.sage --p 0xFFFFFFFF00000001 --alpha 3 --rounds 2 \
  --outfile data/ANEMOI_custom.in

# CLI flags
--p <prime|0xHEX|NAME>
    Field prime. Accepts a decimal integer, a 0x... hex literal, or a symbol defined in scripts/constants.py.

--alpha <odd>
    S-box exponent; must satisfy gcd(alpha, p-1) = 1.

--rounds <r>
    Number of rounds (P_CICO, l = 1).

--outfile <path>
    Optional output file. Default: data/ANEMOI_p<p>_r<r>_a<alpha>_PCICO[ _noFLL ].in

--ordering {1|2|3}
    Variable ordering in the P_CICO model. Default: 1 (recommended).

--no-final-ll
    Omit the final linear layer constraint. Default: include it.
```

### 9. **ANEMOI algebraic-closure analysis & solution enumeration**

**Script:** `extract_anemoi_solutions.sage`

**What it does:**  
Given a **LEX** Gröbner basis over GF(p), locates the univariate **eliminant** in the last variable, reports algebraic-closure solution counts (total with multiplicity = `deg(f_last)`, distinct = `deg(squarefree(f_last))`), provides a factor-degree breakdown over GF(p), and (optionally) **enumerates concrete solutions** in minimal extensions GF(p^d) for small factor degrees.

**Outputs:**  
- `results/<stem>_AC_report.txt` — counts + factor-degree breakdown  
- `results/<stem>_AC_solutions.txt` — optional: enumerated solutions over GF(p^d)  
- `logs/<stem>_AC_extract.log` — detailed log (shape checks, timing, enumeration summary)

**Example usage**
```bash
# Produce algebraic-closure counts and factor breakdown (no enumeration limit changes)
sage scripts/extract_anemoi_solutions.sage results/ANEMOI_p251_r2_a3_PCICO_LEX.txt

# Enumerate roots only up to degree 4 in the minimal extension, stop after first solution
sage scripts/extract_anemoi_solutions.sage results/ANEMOI_p251_r2_a3_PCICO_LEX.txt \
  --enum-maxdeg 4 --first-only

# Enumerate up to 20 solutions overall; try elimination if no univariate is present
sage scripts/extract_anemoi_solutions.sage results/ANEMOI_p251_r3_a5_PCICO_LEX.txt \
  --max-solutions 20 --elim

# Skip enumeration entirely; only the algebraic-closure count report is produced
sage scripts/extract_anemoi_solutions.sage results/ANEMOI_p101_r3_a5_PCICO_LEX.txt \
  --skip-enum
```