# CHESS <img src='man/figures/logo.png' align="right" height="139" />
<!-- badges: start -->
[![Travis build
status](https://travis-ci.com/T-Heide/CHESS.cpp.svg?token=cqqKUEunNazVxSmjTbrG&branch=r_package)](https://travis-ci.com/T-Heide/CHESS.cpp)
[![codecov](https://codecov.io/gh/T-Heide/CHESS.cpp/branch/r_package/graph/badge.svg?token=JQ7YMKP381)](https://codecov.io/gh/T-Heide/CHESS.cpp)


<!-- badges: end -->

Cancer Heterogeneity with Spatial Simulations (CHESS)

-----

#### Installation

``` r
# install.packages("devtools")
devtools::install_github("T-Heide/CHESS.cpp", ref="r_package")
```
-----

#### Getting started



To run a tumour simulation that simulates a spatial tumour growth, followed by bulk and/or single cell sampling data generation, use the functions:

```SimulateTumourSample()``` and/or ```SimulateTumourTree()```

The functions output sample VAF and/or phylogenetic tree files, respectively.

The package implements three strategies of a tumour growth parameter inference (using the ABC SMC algorithm) by approximating a target tumour, using its:

```
1. Bulk samples (punch or needle biopsies) - ABCSMCwithBulkSamples()
2. Single cell sample phylogenetic trees - ABCSMCwithTreeSampleBL() and ABCSMCwithTreeSampleBT() (using Branch Lengths or Branching Times as summary statistics)
3. Whole tumour bulk sample - ABCSMCwithWholeTumour()
```

Depending on the strategy, a user would need to provide either target tumour bulk sample VAFs (list of R data.frames where each row should correspond to a unique mutation with the following columns: clone (Clone type label), alt (Number of reads), depth (Sequencing depth), id (Unique ID)), an array of whole tumour sample VAFs or single cell sampling phylogenetic trees. Alternatively, a user can provide a set of parameters (please refer to the package documentation for the details of each input parameter format) to simulate a synthetic target tumour to then recover these input parameters.  

The functions output sequence of files containing sets of inferred parameters corresponding to each SMC round (that can then be used to construct the posterior distributions of each parameter).


-----

#### Help and support

[![](https://img.shields.io/badge/GitHub%20Pages-https://t--heide.github.io/CHESS.cpp/-informational)](https://t-heide.github.io/CHESS.cpp/)

-----

#### Copyright

Copyright (C) 2018, Timon Heide (timon.heide@icr.ac.uk)
                  & Kate Chkhaidze (Kate.Chkhaidze@icr.ac.uk)

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 3 of the License, or
(at your option) any later version.
This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program. If not, see http://www.gnu.org/licenses/.

-----

#### Contacts

Timon Heide, _Institute of Cancer Research, London, UK_.

[![](https://img.shields.io/badge/Email-timon.heide@icr.ac.uk-informational.svg?style=social)](mailto:timon.heide@icr.ac.uk)
[![](https://img.shields.io/badge/Github-T--Heide-informational.svg?style=social&logo=GitHub)](https://github.com/T-Heide)


Kate Chkhaidze, _Institute of Cancer Research, London, UK_.

[![](https://img.shields.io/badge/Email-kate.chkhaidze@icr.ac.uk-informational.svg?style=social)](mailto:kate.chkhaidze@icr.ac.uk)
[![](https://img.shields.io/badge/Github-kchkhaidze-informational.svg?style=social&logo=GitHub)](https://github.com/kchkhaidze)
