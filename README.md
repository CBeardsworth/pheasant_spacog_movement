## Repository for: Spatial cognitive ability is associated with speed of movement in the early exploration of an environment

_Christine E. Beardsworth, Mark A. Whiteside, Lucy A. Capstick, Philippa R. Laker, Ellis JG Langley, Ran Nathan, Yotam Orchan, Sivan Toledo, Jayden O. van Horik, Joah R. Madden_



For any questions about the code please contact Christine at c.e.beardsworth@gmail.com

To use any data contained in this repository contact Joah at j.r.madden@exeter.ac.uk for permission.

### code/

The data for the following code sections are found in `data/`, apart from the full filtered movement dataset (used in `1_cleaning and run HMM.R`) as it is too large to be included here. Figures created using this code are stored in `figs/`

- `1_cleaning and run HMM.R`: Some extra cleaning of pheasant movement data and subsetting so that only birds that have 7 days of data (of at least 6 hours) are present in the dataset. The beginning section R code is an example only as the data is not available on this repository. The data can be retrieved from Christine or Joah. However, the code shows how the subsequent files that are included here are created and how we ran Hidden Markov models. 
- `2_Choose HMM and describe states and HRs.R`: This code shows how we chose the best HMM to describe state transitions. We also use the same dataset to create HRs that are shown in Fig 2 of the manuscript. 
- `3_Statistics.R`: This code performs the statistics used in the manuscript as well as some figures.
- `4_ESM.R`: This code produces the figures found in the supplementary material;.


