# Code for "Southern Annular Mode dynamics, projections and impacts in a changing climate" by Purich et al.

This code repository allows for re-creation of figures from Purich et al., "Southern Annular Mode dynamics, projections and impacts in a changing climate".

Most script require data files stored separately [on Zenodo.org](https://doi.org/10.5281/zenodo.17364274). Download all data from that repository and store it inside one dedicated folder.

In general, required python packages are `xarray`, `pandas`, `matplotlib`, `numpy`, `cartopy`, `scipy`.

To reproduce the figures, follow below steps:

## All figures

Set variable `data_path` in each script to the path of the data files downloaded [from Zenodo.org](https://doi.org/10.5281/zenodo.17364274). In all our scripts, `data_path` defaults to `../data/`. So one idea would be to store all data from Zenodo in a folder called "data" under the same folder as you store this code repository (that is, not inside this code folder).

## Figure 1

The script for this figure is `Fig1.ipynb`.

Run all cells, or the script with `jupyter execute Fig1.ipynb`. Resulting figure will be `Fig1.pdf`.

## Figure 2

The main script for this figure is `Fig2.py`. 

Run `python Fig2.py` and resulting figure will be `Fig2.pdf`.

Requires `aostools` available [in this directory](https://github.com/SAMworkshop2024/aostools) or [this repository](https://github.com/mjucker/aostools).

## Figure 3

The main script for this figure is `Fig3.py`.

Run `python Fig3.py` and resulting figure will be `Fig3.pdf`.

Requires `aostools` available [in this directory](https://github.com/SAMworkshop2024/aostools) or [this repository](https://github.com/mjucker/aostools).

## Figure 4

The script for this figure is `Fig4.ipynb`.

Run  all cells. Resulting figure will be `Fig4.pdf`.


## Figure 5

The script for this figure is `Fig5.ipynb`.

Run all cells or exectue notebook. Resulting figures will be `Fig5ac.pdf` and `Fig5bd.pdf`.

## Figure 6

The script for this figure is `Fig6.py`.

Run `python Fig6.py` and resulting figure will be `Fig6.pdf`.

Requires `sacpy` and `cmaps` packages. `sacpy` is included in the code base and does not need to be installed, but it is available via `pip` or `conda`. `cmaps` is available via `pip` or `conda`.

