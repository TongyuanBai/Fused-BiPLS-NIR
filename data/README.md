# Data

The six benchmark datasets are publicly available and are not redistributed in this repository. Download them from the sources below and place the files in this folder (`data/`) under the file names listed. The R and Python scripts read the files from here; another folder can be used by setting the environment variable `FUSED_BIPLS_DATA_DIR`.

| Dataset | File expected in `data/` | Content read by the scripts | Source |
|---|---|---|---|
| Corn | `corn.mat` | `m5spec` (80 × 700, 1100–2498 nm, m5 instrument); response: first column of `propvals` (moisture) | Eigenvector Research, [NIR of corn samples](https://eigenvector.com/resources/data-sets/nir-of-corn-samples-for-standardization-benchmarking/) |
| Diesel | `Filtered_Diesel_Spec_Data.csv` | Column 1: density (ASTM D4052); columns 2–402: spectra at 750–1550 nm (395 samples with a density value) | Southwest Research Institute, distributed by Eigenvector Research, [NIR spectra of diesel fuels](https://eigenvector.com/resources/data-sets/near-infrared-spectra-of-diesel-fuels/) |
| Meat | `Tecator.xlsx` | Column 1: fat content; columns 2–101: absorbance at 850–1050 nm (215 samples) | Tecator data, R package [`fsemipar`](https://search.r-project.org/CRAN/refmans/fsemipar/html/Tecator.html); also [StatLib](http://lib.stat.cmu.edu/datasets/tecator) |
| Milk | `milk.mat` | `X` (67 × 1557, 1000–7224 cm⁻¹); response: `protein` | |
| Soil | `NIRsoil.mat` | `soil` (108 × 1050, 400–2500 nm); response: first column of `soilref` (soil organic matter) | University of Copenhagen, [chemometrics datasets](https://ucphchemometrics.com/datasets/) (NIR soil); Rinnan and Rinnan, *Soil Biology and Biochemistry* 39 (2007) 1664–1673 |
| Tablet | `NIRdata_tablets.MAT` | Columns 4–407 of `Matrix` (310 × 404, 7400–10,500 cm⁻¹); response: column 1 of `Matrix` (active ingredient) | University of Copenhagen, [chemometrics datasets](https://ucphchemometrics.com/datasets/) (Tablets); Dyrby et al., *Applied Spectroscopy* 56 (2002) 579–585 |

File names are matched without regard to case, so `NIRdata_tablets.mat` and `NIRdata_tablets.MAT` are both accepted.
