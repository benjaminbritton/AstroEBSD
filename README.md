# AstroEBSD
AstroEBSD - indexing tool for EBSD patterns

The main repository can be found https://github.com/benjaminbritton/AstroEBSD/

The paper describing this work can be found: 

Thomas Benjamin Britton, Vivian Tong, Jim Hikey, Alex Foden, and Angus Wilkinson "AstroEBSD: exploring new space in pattern indexing with methods launched from an astronomical approach" Journal of Applied Crystallography (2018) https://doi.org/10.1107/S1600576718010373 (& preprint https://arxiv.org/abs/1804.02602 )


The RTM module is presented in:

Alex Foden, David Collins, Angus Wilkinson, and Thomas Benjamin Britton "Indexing electron backscatter diffraction patterns with a refined template matching approach" Ultramicroscopy (2019) https://doi.org/10.1016/j.ultramic.2019.112845 (& preprint https://arxiv.org/abs/1807.11313 )


...with PCA as described in:

Angus Wilkinson, David Collins, Yevhen Zayachuk, Rajesh Korla, Arantxa Vilalta-Clemente "Applications of multivariate statistical methods and simulation libraries to analysis of electron backscatter diffraction and transmission Kikuchi diffraction datasets" (https://doi.org/10.1016/j.ultramic.2018.09.011)

...and with EDS for phase identification:

Thomas McAuliffe, Alex Foden, Chris Bilsland, Dafni Daskalaki Mountanou, and Thomas Benjamin Britton "Advancing characterisation with statistics from correlative electron diffraction and X-ray spectroscopy, in the scanning electron microscope" Ultramicroscopy (2020) (https://doi.org/10.1016/j.ultramic.2020.112944 & preprint http://arxiv.org/abs/1908.04084 )

An example (large, full resolution) correlative EBSD/EDS dataset for RTM and PCA is available at https://doi.org/10.5281/zenodo.3617455.
A compressed version (900 MB) is available at https://doi.org/10.5281/zenodo.3737987.

======

Master has been updated 25/09/2019

Changes include:

cyan/yellow colouring fix

new background correction options

mac/windows loader fixed

With help from Alex Foden & Tom McAuliffe

======

Master updated 12/02/2020 - Tom McAuliffe

Change include:

New phase data structures: phase files, cifnames and masterpatterns (.bin files) specified. 
Phase_Builder_RTM function updated accordingly to handle these.

Integration of the RTM and PCA modules.
- RTM example deck for indexing of a full map with a single candidate phase.
- PCA example decks for clustering of combined EBSD/EDS datasets, and EBSD-based phase ID and indexing.
- PCA postanalysis decks for chemical analysis (as function of phase) and RC-EBSP comparison to template simulations.
- Available with (adjustable) spatial weighting kernel.

Gitignores for phase files / folders fixed




