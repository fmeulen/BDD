# Bdd: Julia code for Bayesian decompounding of discrete distributions

## Description

This Julia code is accompanying the paper *Decompounding discrete distributions: a non-parametric Bayesian approach* by Shota Gugushvili, Ester Mariucci and Frank van der Meulen. The following datasets are analysed (for details see the paper):

- In Section 3.1: synthetic data with a uniform base distribution. The data under settings (a), (b) and (c) are stored as 'testdat_A1.csv', 'testdat_A2.csv' and 'testdat_A3.csv', respectively.
- In Section 3.2: synthetic data with a geometric base distribution. The data under settings (a) and (b) are stored as 'testdat_B1.csv' and 'testdat_B2.csv', respectively.
- In Section 4.2: the plant data, stored as 'plantpopulation.csv'. 

All datasets are contained in the folder named 'data'. 

To compute Bayesian estimates, as well as the truncated estimate from Buchmann and Grubel (2004), put all the files (the files in the scr folder and the files in the data folder) into one directory and create a subdirectory named 'out'. Next run the script 'cpp_discr.jl'. Within this file, one can choose the dataset to be analysed by setting 'data_choice' to any of the scenarios available (which are "generated1", "generated2", "testdata_A1", "testdata_A2", "testdata_A3", "testdata_B1", "testdata_B2", "testdata_C", "horsekicks", "plantpopulation"; details are seen in 'setdata.jl'). 

The code for the Monte-Carlo study in Section 3.3 of the paper is contained in 'cpp_discr_mc.jl'. 

Julia dependencies: RCall, Distributions, DataFrames, DelimitedFiles, SpecialFunctions, TimerOutputs, Random

## Citation

If you use this software in your work, consider citing it using the following BibTeX code:

```
@misc{bdd,
  author       = {Gugushvili, Shota and
                  Mariucci, Ester and
                  van der Meulen, Frank},
  title        = {{Bdd: Julia code for Bayesian decompounding of 
                   discrete distributions}},
  month        = {mar},
  year         = {2019},
  doi          = {10.5281/zenodo.2598802},
  url          = {https://doi.org/10.5281/zenodo.2598802}
}
```


