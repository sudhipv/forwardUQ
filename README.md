# forwardUQ

Collection of MATLAB prototypes for forward uncertainty quantification (UQ) of a Poisson diffusion problem with random material properties. Each folder contains a self‑contained workflow that assembles finite‐element operators, propagates lognormal/gaussian uncertainty via stochastic spectral finite element methods (SSFEM), and compares those solutions to Monte Carlo surrogates.


<img width="711" height="890" alt="forward_uq2" src="https://github.com/user-attachments/assets/bd7ebb96-d7aa-4825-884c-dcf7b0c71489" />


<img width="1010" height="817" alt="forward_UQ" src="https://github.com/user-attachments/assets/69004248-6b9f-4795-8f6b-b1dd2b27d8ac" />

## Repository Overview
- `Intrusive` — Intrusive SSFEM for random variables. Key scripts: `ssfem_Intrusive.m`, `AssembleMatrix_ss.m`, `AssembleVector_ss.m`, `Poisson_ss.m`. Includes pre/post-processing data (`*.mat`, `.fig`) plus gmsh meshes under `mac/` and `source/`.
- `NISP` — Non-intrusive spectral projection (Monte Carlo and quadrature variants). Core drivers are `ssfem_NISP.m`, `ssfem_NISP_quadrature.m`, and `Normalized_ssfem_NISP_quadrature.m`; helper PDE definitions live beside gmsh utilities (`ExtractGmsh.m`, `square_refined.m`).
- `NISP_Process` — NISP workflow for lognormal random processes driven by Karhunen–Loève expansions. Uses `NISP_process.m`, `GetKLEterms.m`, `GetLambdaSymmetric_2D.m`, and quadrature data from UQTk.
- `NISP_analytical` — Analytical/non-intrusive reference solutions for verification (`NISP_analytical.m`, `NISP_quadrature.m`, `PC_std_quadrature.m`) plus saved Monte Carlo statistics (`MCS_*.mat`).
- `ssfem_LNP_Intrusive` — Intrusive SSFEM for lognormal random processes (`ssfem_Intrusive.m`, `ssfem_Intrusive_Normalized.m`, `GetKLEterms.m`, `Cijk.m`). Contains meshes (`mac/`, `source/`) and comparison figures.
- `ssfem_MCS_Process` — Monte Carlo sampling of random processes (`sfem_MCS_Process.m`) using the same assembly utilities as the intrusive lognormal process workflow.
- `ssfem_MCS_RV` — Monte Carlo sampling for single random variables (`ssfem_MCS.m`) that mirrors the non-intrusive RV setup.

Shared helper files include:
- `AssembleMatrix.m`, `AssembleMatrix_ss.m`, `AssembleVector.m`, and `AssembleVector_ss.m` – generic stiffness/load assemblers for deterministic vs. stochastic forms.
- `Poisson.m`, `Poisson_ss.m`, `Poisson_ss_Normalized.m` – PDE weak forms with configurable random diffusivity.
- `Cijk.m` – generates Hermite polynomial triple products (`Cijk*.mat`) and norm-squared data via UQTk (`gen_mi`).
- `GetKLEterms.m`, `GetLambdaSymmetric_2D.m`, `GetKLEterms_Normalized.m` – build KLE eigenpairs and evaluate lognormal scaling terms.
- Mesh utilities such as `ExtractGmsh.m`, `square.m`, `square_refined.m`, and gmsh sources under each `{module}/mac` or `{module}/source` directory.

## Requirements
- MATLAB R2017a or newer (scripts rely on base language plus `pdesurf` from PDE Toolbox).
- [Gmsh](https://gmsh.info/) to generate `.msh` files from the `.geo` geometries stored under `mac/` or `source/`.
- [UQTk](https://uqtk.org/) (tested with v3.0.4) for `gen_mi` and `generate_quad`. Update `uqtk_path` inside the scripts to match your install.
- Access to the `misc/cijk_ord_dim` directory referenced by several scripts for precomputed `Cijk*.mat`, `norm_squared*.mat`, and `mindex*.mat/.dat` files. Recreate these files with `Cijk.m` and UQTk if they are not present.

## Typical Workflow
1. **Prepare a mesh**
   - Edit or reuse the `.geo` files in `*/source` or `*/mac`.
   - Run Gmsh to export a version 2 ASCII `.msh`.
   - Convert the mesh with `ExtractGmsh.m`, which writes `points.txt`, `edges.txt`, and `triangles.txt`. Scripts such as `square_refined.m` can also be run directly in MATLAB to load a hard-coded mesh into `p`, `e`, `t`.
2. **Generate spectral data (if needed)**
   - Use `Cijk.m` to compute multiplication tensors for the desired polynomial order and dimensionality. The script writes `CijkXXXXX.mat` and `norm_squaredXXXXX.mat` files whose suffixes encode order/dimension (e.g., `Cijk030003.mat`).
   - Create multi-index tables with UQTk’s `gen_mi` and place the resulting `mindex*.dat` alongside the solver that needs it.
   - For random processes, run `GetLambdaSymmetric_2D.m` to sample the covariance eigenpairs and reuse them through `GetKLEterms.m`.
3. **Run a solver**
   - **Intrusive RV** – from `Intrusive/`, set `ord_in`, `ord_out`, `dim`, and `num_spectral` in `ssfem_Intrusive.m`, ensure the correct `Cijk`/`norm_squared` paths, then run the script to obtain polynomial chaos coefficients, means, and standard deviations.
   - **Non-intrusive RV (NISP)** – use `ssfem_NISP.m` for Monte Carlo projection or `ssfem_NISP_quadrature.m` / `Normalized_ssfem_NISP_quadrature.m` for Gauss–Hermite quadrature. Configure sample counts or quadrature orders and Hermite coefficients (`kappa(*)`).
   - **Intrusive/NISP lognormal process** – use the drivers in `ssfem_LNP_Intrusive/` or `NISP_Process/`. These scripts expect `lambda_2D.mat`, `mindex*.dat`, and `norm_squared*.mat`, and they call UQTk’s `generate_quad` internally; update the `uqtk_path` variable and normalize inputs to match your workspace.
   - **Monte Carlo baselines** – `ssfem_MCS_RV/ssfem_MCS.m` and `ssfem_MCS_Process/sfem_MCS_Process.m` sample random inputs directly to compare against the spectral solvers.
4. **Post-process**
   - Most scripts plot `pdesurf` surfaces for the mean field, higher-order coefficients, and standard deviations while also saving `.fig`, `.mat`, or KDE data (`MCS_50000_*.mat`, `compare_*.fig`, `ssfem_intrsv.mat`, etc.).
   - `NISP_analytical/` holds scripts for closed-form coefficient evaluations and PDF reconstruction at specific nodes to benchmark the numerical schemes.

## Usage Notes
- Add the desired module folder to the MATLAB path before invoking any driver script so that shared assemblers and PDE definitions resolve correctly.
- Several scripts contain hard-coded `cd` statements (for example, to `/Users/sudhipv/...`) around UQTk calls or result exports; replace those paths with locations on your machine to avoid runtime failures.
- `ExtractGmsh.m` assumes a 2D mesh saved in Gmsh 2 ASCII format where edges precede triangle elements. Confirm that your `.msh` files follow the same layout or update the reader accordingly.
- Normalized vs. non-normalized polynomial chaos: the “Normalized” scripts expect Hermite norms loaded from `norm_squared*.mat`. Be sure that the file suffix (`_O11_D1`, `030003`, etc.) matches the `(order, dimension)` pair configured at the top of the solver.
- The Monte Carlo and NISP scripts often loop over ~10k–50k samples. Adjust `n`, `n_sample`, or `num_qd` for quicker smoke tests when validating new setups.

## Reference

> **[Scalable Domain Decomposition Methods for Nonlinear and Time-Dependent Stochastic Systems](https://doi.org/10.22215/etd/2023-15817)**

**Authors:** Vasudevan, Padillath and Sharma, Sudhi  
**Institution:** Carleton University (2023)  
**DOI:** [10.22215/etd/2023-15817](https://doi.org/10.22215/etd/2023-15817)

<details>
<summary><b>Click to expand BibTeX citation</b></summary>

```bibtex
@phdthesis{vasudevan2023scalable,
  title={Scalable Domain Decomposition Methods for Nonlinear and Time-Dependent Stochastic Systems},
  author={Vasudevan, Padillath and Sharma, Sudhi},
  year={2023},
  school={Carleton University},
  doi={10.22215/etd/2023-15817}
}
\```
</details>



## Questions?
Contact : Sudhi Sharma P V  
Email: sudhisharmapadillath@gmail.com
