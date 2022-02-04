# Dynamical Dirichlet Mixture Model

In this repo, we share the source codes written in Matlab back in 2007 and the corresponding report from the public [official link](https://infoscience.epfl.ch/record/146114/files/) to Github with the
hope that this model could find more applications.

* For technical details, please refer the [technical report](./le-idiap-rr-07-02.pdf):
```
Le Chen, David Barber, Jean-Marc Odobez, "Dynamical Dirichlet Mixture Model", IDIAP technical report, IDIAP-07-02, 2007, April.
```
* Please refer [contents.m](./source_codes/contents.m) for a brief overview of the codes.
* For the impatient, run [Script_Test_EstHMMDM.m](./source_codes/Script_Test_EstHMMDM.m) directly to have a flavor.

# For Octave
* If you want to run the codes under [Octave](https://www.gnu.org/software/octave/index), you need
	to install the statistics package. For Ubuntu, you need to run
	```
sudo apt-get install octave-statistics
	```
* You need to load the package at the beginning:	`pkg load statistics`. For example, start Octave in the
	directory [source_codes](./source_codes) and run:
		```
octave:> pkg load statistics
octave:> run("Script_Test_EstHMMDM.m")
		```


# Others
* Any comments and suggestions are extremely welcome.
* Citations from [Google Scholar](https://scholar.google.com/scholar?oi=bibs&hl=en&cites=13445517119047258916&as_sdt=100)
