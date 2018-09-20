# denovo_example

This repository contains code providing a basic example of an application of Kwonsang Lee's Denovo package.
This code uses the lalonde dataset. There is a file that contains code that matches observations within the lalonde dataset
and a file that uses the matched data in the implementation of the denovo method.

# Files

**lalonde.csv**: The lalonde dataset stored in csv format

**match_lalonde_denovo.R**: The code used to generate a dataset in the format required by denovo.
This code contains parameters in the begining that can be changed to facilitate use with other datasets. Generates
a csv representing the dataframe required by denovo.

**run_denovo.R**: The code implementing the denovo method from the user perspective. Downloads the package, installs it,
runs the method, then saves all generated objects.
