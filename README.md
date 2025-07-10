# Algebraic Cryptanalysis Pipeline: Gröbner Bases and Field Extension Reduction

This repository implements a **modular pipeline for solving systems of polynomial equations over finite fields**, with a focus on cryptographic applications (e.g., HFE, ANEMOI). The pipeline includes automatic support for field extensions, efficient Gröbner basis computation (F4/F5), basis conversion (FGLM), solution extraction, and verification utilities.

---

## Pipeline Overview

1. **Field Extension Expansion** (SageMath):  
   Transforms systems over extension fields 𝐹<sub>pⁿ</sub> into equivalent systems over the base field 𝐹<sub>p</sub>, introducing extra variables/equations as needed.

2. **System Diagnostics & Metadata Extraction** (SageMath):  
   Computes structural properties of the polynomial system and ideal (e.g., dimension, degree, homogeneity, etc.).

3. **Grobner Basis Computation (F4/F5)** (Julia/AlgebraicSolving.jl):  
   Computes a Gröbner basis in degree reverse lexicographic (DRL) order.  
   - F4 (matrix-based, multi-threaded)  
   - F5 (signature-based, single-threaded)

4. **FGLM Basis Conversion** (SageMath/Singular):  
   Converts the DRL basis to lexicographic (LEX) order using the FGLM algorithm.

5. **Solution Extraction** (SageMath):  
   Extracts and verifies all solutions of the system from the lex basis.

6. **Mapping Solutions Back to Extension Field** (SageMath):  
   For systems originally over 𝐹<sub>pⁿ</sub>, reconstructs extension field solutions from base field coordinates.

7. **Correctness & Expansion Verification** (SageMath):  
   Checks correctness of field extension expansion and solution mappings.

All steps produce logs and human-readable output files, organized in the `logs` and `results` folders.

---

## Script Reference

Below, each script is described with its role, how to invoke it, and the expected input/output files.

---

### 1. **Field Extension Expansion**

**Script:** `scripts/expand_field_extension_to_base.sage`  
**What it does:**  
- Converts a polynomial system over 𝐹<sub>pⁿ</sub> to an equivalent system over 𝐹<sub>p</sub>.
- Writes new system to a `.in` file in `data/`.

**How to use:**
```sh
sage scripts/expand_field_extension_to_base.sage data/<input_system>.in
```

**Output:** `data/<input_system>_expanded.in`

### 2. **System Diagnostics** 
**Script:** `scripts/system_diagnosis.sage`  
**What it does:**  
- Prints and logs properties of the polynomial system and the corresponding ideal:
    - Number of variables/equations
    - Degree statistics
    - Sparsity
    - Homogeneity
    - Krull dimension
    - Zero-dimensionality
    - Quadratic/Boolean check

**How to use:**
```sh
sage sage scripts/system_diagnosis.sage data/<system>.in
```

**Output:** `logs/<system>_DIAGNOSIS.log`

### 3. **Grobner Basis Computation (F4/F5)**

**Scripts:** 
- `scripts/solve_F4_from_file.jl`  
- `scripts/solve_F5_from_file.jl`

**What they do:**  
- Compute a Grobner basis in DRL (degrevlex) order using F4 (multi-threaded) or F5 (signature-based).
- Logs the computation and saves basis in a results file.

**How to use:**
```sh
julia scripts/solve_F4_from_file.jl data/<system>.in <num_threads>

julia scripts/solve_F5_from_file.jl data/<system>.in
```

**Outputs:** 
- `results/<system>_F4_<timestamp>.txt`
- `results/<system>_F5_<timestamp>.txt`

### 4. **FGLM Conversion (DRL → LEX Order)**

**Script:** `scripts/convert_to_lex_fglm.sage`  

**What it does:**  
- Converts a DRL Groebner basis (from Julia) to LEX order using the FGLM algorithm (via Singular).
- Saves new basis and logs statistics/timings.

**How to use:**
```sh
sage scripts/convert_to_lex_fglm.sage results/<system>_F4_<timestamp>.txt
```

**Outputs:** 
- `results/<system>_F4_<timestamp>_LEX.txt`
- `logs/<system>_F4_<timestamp>_FGLM.log`

### 5. **Solution Extraction**

**Script:** `scripts/extract_solutions_from_lex.sage`  

**What it does:**  
- Extracts all solutions from the LEX Groebner basis file (using variety() or brute-force if needed).
- Verifies correctness and logs stats.

**How to use:**
```sh
sage scripts/extract_solutions_from_lex.sage results/<system>_F4_<timestamp>_LEX.txt
```

**Outputs:** 
- `results/<system>_F4_<timestamp>_LEX_sols.txt`
- `logs/<system>_F4_<timestamp>_LEX_SOLUTIONS.log`

### 6. **Mapping Solutions to Extension Field**

**Script:** `scripts/map_base_field_solutions_to_extension.sage`  

**What it does:**  
- For systems originally defined over 𝐹<sub>pⁿ</sub>, maps solutions in base field coordinates back to original extension field variables.

**How to use:**
```sh
sage scripts/map_base_field_solutions_to_extension.sage results/<system>_F4_<timestamp>_LEX_sols.txt data/<original_system>.in
```

**Output:** 
- `results/<system>_F4_<timestamp>_LEX_sols_mapped.txt`

### 7. **Expansion Correctness Verification**

**Script:** `scripts/check_expansion_correctness.sage`  

**What it does:**  
- Verifies that the expanded system over 𝐹<sub>p</sub> is equivalent to the original over 𝐹<sub>pⁿ</sub>, by checking assignments (exhaustively for small systems).

**How to use:**
```sh
sage scripts/check_expansion_correctness.sage data/<original_system>.in data/<original_system>_expanded.in
```

**Output:** 
- Prints verification result to console.