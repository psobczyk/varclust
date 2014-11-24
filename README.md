# README #


### What is this repository for? ###

* varclust is a R package for clustering variables
* Current version is 0.9.20

### How do I get set up? ###

* Install package using devtools package
* Download the package and install it manually from R console
* No configuration is needed
* Dependencies will install automatically
* Read the vignette to get familiar with basic usage

### Who do I talk to? ###
* If help provided in the package does not solve your problem please contact Piotr.Sobczyk[malpka]pwr.edu.pl

### To do ###

#### Iterative estimation:

1. remove estimating subspace in mlcc.bic - it is now done in mlcc.reps. Right now the 
  way of computing sigma is not coherent with mlcc.kmeans

#### Initialization:

1. implement hierarchical clustering (similarly to clustofvar)

#### Testing

1. write tests for main aux functions

#### Documentation

1. update mlcc.kmeans variables names to meet the dot standard
2. rewrite documentation with @inheritParams
