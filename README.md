# Algebraic Cryptanalysis Pipeline: Gröbner Bases and Field Extension Reduction

This repository implements a **modular pipeline for solving systems of polynomial equations over finite fields**, with a focus on cryptographic applications (e.g., HFE, ANEMOI). The pipeline includes automatic support for field extensions, efficient Gröbner basis computation (F4/F5), basis conversion (FGLM), solution extraction, and verification utilities.

---

## Pipeline Overview

1. **Field Extension Expansion** (SageMath):  
   Transforms systems over extension fields 𝐹<sub>pⁿ</sub> into equivalent systems over the base field 𝐹<sub>p</sub>, introducing extra variables/equations as needed.

2. **System Diagnostics & Metadata Extraction** (SageMath):  
   Computes structural properties of the polynomial system and ideal (e.g., dimension, degree, homogeneity, etc.).

3. **Groebner Basis Computation (F4/F5)** (Julia/AlgebraicSolving.jl):  
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
- Converts a polynomial system over \(\mathbb{F}_{p^n}\) to an equivalent system over 𝐹<sub>p</sub>.
- Writes new system to a `.in` file in `data/`.

**How to use:**
```sh
sage scripts/expand_field_extension_to_base.sage data/<input_system>.in


