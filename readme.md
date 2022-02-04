# Dynamical Dirichlet Mixture Model

In this repo, we share the source codes written in Matlab back in 2007 and the corresponding report from the public [official link](https://infoscience.epfl.ch/record/146114/files/) to Github with the
hope that this model could find more applications.

* For technical details, please refer the [technical report](./le-idiap-rr-07-02.pdf):
```
Le Chen, David Barber, Jean-Marc Odobez, "Dynamical Dirichlet Mixture Model", IDIAP technical report, IDIAP-07-02, 2007, April.
```
* Please refer [contents.m](./source_codes/contents.m) for a brief overview of the codes.
* For the impatient, run [Script_Test_EstHMMDM.m](./source_codes/Script_Test_EstHMMDM.m) directly to have a flavor.

# Migrating to Octave  (not working yet)
* [Octave](https://www.gnu.org/software/octave/index) is a free software counterpart of Matlab. They
	are sufficiently similar but not the same. To run some codes, you need to install a statistics
	package. For example on Ubuntu, you need to run:
```
sudo apt-get install octave-statistics
```
* Then you need to load the package at the beginning:	`pkg load statistics`. For example, start Octave in the
	directory [source_codes](./source_codes) and run:
```
octave:> pkg load statistics
octave:> run("Script_Test_EstHMMDM.m")
```
* The script is still not running well on Octave unfortunately. It works well under Matlab as it was
	originally written under Matlab.

# Others
* Any comments and suggestions are extremely welcome.
* Citations from [Google Scholar](https://scholar.google.com/scholar?oi=bibs&hl=en&cites=13445517119047258916&as_sdt=100)
* Similar model can be found here:
	*  [Dynamic Mixture Models for Multiple Time-Series](https://dl.acm.org/doi/10.5555/1625275.1625744).

# To do
* [ ] Translate this package to Python.
