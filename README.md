# 1. Introduction

This repository contains the source code of our space ACA on data nodes, space ACA on internal nodes and time ACAs.

# 2. Contents

- Space ACA on data nodes
  - [MCK implementation](https://github.com/ruiyang00/aca_dlis_review/blob/master/attack/attack.h): This file contains main logic of the attack (MCK implementation)
  - [Source code](https://github.com/ruiyang00/aca_dlis_review/blob/master/src/benchmark/space_aca_dn.cpp)
    - White-box: setting=1
    - Gray-box: setting=2
      - python [script](https://github.com/ruiyang00/aca_dlis_review/blob/master/scripts/populate_graybox_dataset.py) to generate dataset to 
 constrcut a "substitude" tree (described in paper).
  - Scripts to regenerate the result in Figure 2
    - [White-box](https://github.com/ruiyang00/aca_dlis_review/blob/master/scripts/run_space_aca_dn_whitebox.sh)
    - [Gray-box](https://github.com/ruiyang00/aca_dlis_review/blob/master/scripts/run_space_aca_dn_graybox.sh)
- Space ACA on internal nodes
  - Source code
    - [Our attack](https://github.com/ruiyang00/aca_dlis_review/tree/master/src/benchmark/space_aca_in.cpp)
    - [SZEGP Attack](https://github.com/ruiyang00/aca_dlis_review/tree/master/src/benchmark/lis.cpp)
  - Script to regenerate the result in Table 2
    - [Black-box](https://github.com/ruiyang00/aca_dlis_review/tree/master/scripts/run_space_aca_in_blackbox.sh) 
- Time ACA on ALEX
  - Source code
    - [White-box](https://github.com/ruiyang00/aca_dlis_review/tree/master/src/benchmark/time_aca_whitebox.cpp)
    - [Gray-box](https://github.com/ruiyang00/aca_dlis_review/tree/master/src/benchmark/time_aca_graybox.cpp)
      - python [script](https://github.com/ruiyang00/aca_dlis_review/blob/master/scripts/populate_graybox_dataset.py) to generate dataset to constrcut a "substitude" tree. 
    - [Black-box](https://github.com/ruiyang00/aca_dlis_review/tree/master/src/benchmark/time_aca_blackbox.cpp)
  - Scripts to regenerate the results in Figure 5 ( Note: Before mounting any time ACA attacks to modified ALEX (forcing catastrohpic expansion instead of catastrohpic split) or vanilla ALEX, we need to comment/uncomment line 1316 to 1321 in [alex.h](https://github.com/ruiyang00/aca_dlis_review/blob/master/src/core/alex.h), as described in paper)
    - [White-box](https://github.com/ruiyang00/aca_dlis_review/tree/master/scripts/run_time_aca_whitebox.sh)
    - [Gray-box](https://github.com/ruiyang00/aca_dlis_review/tree/master/scripts/run_time_aca_graybox.sh)
    - [Black-box](https://github.com/ruiyang00/aca_dlis_review/tree/master/scripts/run_time_aca_blackbox.sh)
    - [Legitmate workload](https://github.com/ruiyang00/aca_dlis_review/tree/master/scripts/run_time_aca_legit.sh)

- [ALEX](https://github.com/ruiyang00/aca_dlis_review/tree/master/src/core): ALEX's source code.

# 3. Dependencies for Attacks
- Space ACA on data nodes: [OR-Tools](https://github.com/google/or-tools) and [KDE](https://scikit-learn.org/stable/install.html)
- Space ACA on internal nodes: No dependencies
- Time ACA on ALEX: [KDE](https://scikit-learn.org/stable/install.html)

# 4. Datasets
Datasets (Longitudes, Longlat, Lognormal, YCSB) used in this paper could be found in [ALEX](https://github.com/microsoft/ALEX) repo.


