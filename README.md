# Solid_Hydrogen2

Repository for Solid Hydrogen analysis tools and notebooks. Reproduces thesis work by Camilla Roselli: https://amslaurea.unibo.it/id/eprint/37654/

Purpose
- Produce, inspect and prepare reconstructed Solid Hydrogen event samples for analysis and unfolding of neutrino energy.
- Provide interactive event displays and plotting utilities.

Repository layout
- `Generate_file_lists/` — scripts to create file lists for processing (input lists).
- `Snap_file_processing/` — data-processing scripts that produce the processed ROOT trees.
- `Jupyter_analysis/` — analysis notebooks and plotting scripts (event displays, performance studies, pre-unfolding histograms).
- `Unfolding/` — scripts and notebooks to prepare and run unfolding of neutrino energy spectra.
- `Solid_Hydrogen2.sh` — top-level helper script to run subfolders' scripts.

Quick start
1. Ensure dependencies are installed (see Dependencies).
2. Set the input ROOT file path inside the notebooks (they use `Rough_tree` from `Processed_tree.root`) or run the processing scripts in `Snap_file_processing/` to create processed trees.
3. Run analyses or notebooks in `Jupyter_analysis/` (examples: `SH_neutron_reco_performance.ipynb`, `SH_Event_display.ipynb`, `SH_pre_unfolding.ipynb`).

Run the top-level script (interactive):
```bash
bash Solid_Hydrogen2.sh
```
A batch script that is submittable to the CNAF machines is also provided and can be run as:
```bash
condor_submit Solid_Hydrogen2_batch.sh
```

Typical notebook usage
- Open the notebook, set `df = RDataFrame("Rough_tree", "<path-to>/Processed_tree.root")` and `Out_dir` to your desired output folder, then run cells in order.
- `SH_Event_display.ipynb` lets you set `event_index` to visualize a single event (saves HTML with `fig.write_html(...)`).

Dependencies
- ROOT with PyROOT (ROOT 6+)
- Python packages: `numpy`, `plotly` (and optionally `kaleido` to export images), plus any local helper modules (e.g., `df_functions_utils.py`).
- Build tools and compilers if using the C++ code under `src/` and `include/`.
- A compiled version of sandreco (legacy version) https://github.com/DUNE/sandreco

Notes
- The notebooks rely on a local helper module `df_functions_utils.py` for dataframe extension and selection flags. Keep that module on `PYTHONPATH` or in the same folder as the notebooks.
- Outputs (PDF/PNG plots and ROOT histogram files) are written under `Plots/` subfolders by default—ensure those directories exist or set `Out_dir` to an existing path.

Contact / Maintainer
- Owner: `Trenchcoat95` (repository `Solid_Hydrogen2`)

License
- (Add license information here if applicable)
