paleochrono
===========

A statistical and physical model to optimize chronologies of paleoclimatic sites.


What this manual is and is not?
-------------------------------

This manual is a documentation on how to use the paleochrono software.  
It is _not_ a description of the paleochrono principles. Please read to the scientific articles
describing paleochrono for that purpose:\
Parrenin, F., Bouchet, M., Buizert, C., Capron, E., Corrick, E., Drysdale, R., Kawamura, K., Landais, A., Mulvaney, R., Oyabu, I., and Rasmussen, S. O.: 
The Paleochrono-1.1 probabilistic model to derive a common age model for several paleoclimatic sites using absolute and relative dating constraints, 
_Geoscientific Model Development_, 17, 8735–8750, https://doi.org/10.5194/gmd-17-8735-2024, 2024.
It is _not_ an operating system or python documentation.
Please use your operating system or python documentation instead.


Where can I get help on paleochrono?
------------------------------------

A mailing list has been set up on Google Groups:  
https://groups.google.com/forum/?hl=en#!forum/icechrono  
You just need a google account to access this mailing list.
You can also directly email to Frédéric Parrenin: frederic.parrenin@univ-grenoble-alpes.fr

How to download paleochrono?
----------------------------

Go here:  
https://github.com/parrenin/paleochrono/  
and click on the donwload button.  
In the downloaded folder, you will find the following files:
- README.md		: is the current documentation of paleochrono.
- LICENCE		: is the paleochrono licence file.
- paleochrono.py		: is the main paleochrono program that you will run.
- pccfg.py, pcmath.py, pcsite.py and pcsitepair.py	: are python modules used by paleochrono.py
- Clean.py		: is a python script to clean a dating experiment directory
- AICC2023-Hulu		: is an example experiment directory: it contains all the necessary
numerical settings, prior information and observations for the different ice cores in the AICC2023
dating experiment and for the MSD and MSL Hulu speleothem. It takes a few minutes to run an a recent
computer.
- AICC2023-LowRes : is the experiment directory used to create the AICC2023 ice core chronology.
It uses a lower resolution than the official AICC2023 experiment so that it runs faster.
- Hulu_MSL : is an experiment directory used in the Paleochrono-1.1 article for only the Hulu/MSL speleothem.
It is a simple example that you can follow when you date only one record from a simple archive.

What do I need to run paleochrono?
----------------------------------

paleochrono is a scientific python3 software, therefore you need a scipy distribution.  
paleochrono is developed and tested using the anaconda distribution, therefore we recommend it.  
Anaconda can be downloaded here (use the python3 version):  
https://www.anaconda.com/download 

Paleochrono probably works on other scipy distributions, provided they contain the following python
modules:  
- sys
- os
- time
- math
- numpy
- matplotlib
- warnings
- scipy
- yaml
- gc
- multiprocessing
- warnings
- pickle
- iosacal (optional, for C14 calibration)


How to run paleochrono?
---------------------

Assuming you use anaconda, you can go in the spyder IDE and type the following commands in the
ipython interpreter:

```
cd path-to-palechrono
run paleochron.py exp_directory/
```

where `path-to-paleochrono` is the directory containing paleochrono and `exp_directory` is the name of
your experiment directory. 
The `AICC2023-Hulu` experiment directory is provided for you convenience.
It is an AICC2023 + Hulu experiment.
It takes a few minutes to run on a recent computer.


What are the outputs of a run:
------------------------------

If the run went correctly, it has created output files.

In the main directory, you have the following output file:
- `output.txt`: contains the size of the residuals and variables vectors, the initial and final
values of the cost function, the program execution times and the memory usage.
- `obs_residuals.pdf` : is a distribution diagram for the observations

In each site directory, you have the following output files:
- `output.txt`: is the main output file. It gives you the posterior estimates and
uncertainties of the three input variables (accumulation, LID and thinning) and of the output
 variables (ice age, air age, Δdepth, etc.). The header in the file tells you which column is which
 variable.
- `restart.txt`: is a restart file, which can be used to start an optimization experiment from the
result of a previous optimization experiment, for a faster convergence.
- `deposition.pdf`: is the deposition rate figure
- `deposition_log.pdf`: is the deposition rate figure with a log scale
- `age.pdf`: is the age for a non-ice-core archive
- `ice_age.pdf` and `air_age.pdf`: are the ice and air age figures for an ice core
- `delta_depth.pdf`: is the Δdepth figure for an ice core
- `ice_layer_thickness` and `air_layer_thickness.pdf`: are the ice and air layer thickness figures
for an ice core
- `lock_in_depth.pdf`: is the Lock-In Depth figure for an ice core
- `thinning.pdf`		: is the thinning figure for an ice core
- `obs_residuals.pdf` : is a distribution diagram for the observations related to this site

In each site-pair directory, you have the following output files:
- `synchro.pdf`: is the stratigraphic links figure for a pair of non-ice-cores archives
- `ice_synchro.pdf`: is the ice stratigraphic links figure when there is one ice core in the pair
- `air_synchro.pdf`: is the air stratigraphic links figure when there is one ice core in the pair
- `air_air_synchro.pdf`		: is the air-air stratigraphic links figure for a pair of ice cores
- `air_ice_synchro.pdf`		: is the air-ice stratigraphic links figure for a pair of ice cores
- `ice_air_synchro.pdf`		: is the ice-air stratigraphic links figure for a pair of ice cores
- `ice_ice_synchro.pdf`		: is the ice-ice stratigraphic links figure for a pair of ice cores
- `residuals.pdf` : is a distribution diagram for the observations of this site pair


How to clean an experiment directory after a run?
-------------------------------------------------

If your run was successful, it has produced output files and figure files.
To clean it from the results files, you can run the following command in ipython:

```
run Clean.py exp_directory/
```


What is the structure of an experiment directory?
-------------------------------------------------

You can have a look at the provided `AICC20123-Hulu` directory.
You need to specify your prior scenarios for deposition rate (in all cases) and LID and thinning
(for an ice core) and your age observations.
Note that the .txt files should be space or tab separated values (do not use commas).

You have five general files:
- `parameters.yml`: contains general parameters for the
experiment
- `parameters_all_sites.yml`: defines site parameters that are the same
for all sites (there are overidded by site specific parameters).
- `parameters_covariance_observations_all_sites.py`: defines the covariance of the
observations that are the same for all sites  (there are overidded by site specific parameters).
- `parameters_covariance_observations_all_site_pairs.py`: defines the covariance for the
observations that are the same for all site pairs  (there are overidded by site pair specific
parameters).

Then you have one directory per site, which contains:
- `parameters.yml`: all the site specific parameters
- `parameters_covariance_observations.py`: this file allows to define the correlation of site
 specific observations
- `deposition.txt`: depth / background accu (in ice-equivalent) / unc. (opt., relative)
- `age_horizons.txt`: depth / age / unc. (absolute) for dated horizons for a non-ice-core
- `age_horizons_C14.txt`: depth / age / unc. (absolute) for C14-dated horizons for a non-ice-core
- `age_intervals.txt`: epth\_top / depth\_bottom / duration / unc. (absolute) for intervals for a non-ice-core
- `ice_age_horizons.txt`: depth / age / unc. (absolute) for ice dated horizons for an ice core
- `air_age_horizons.txt`: depth / age / unc. (absolute) for air dated horizons for an ice core
- `ice_age_intervals.txt`: depth\_top / depth\_bottom / duration / unc. (absolute) for ice intervals for an
ice core
- `air_age_intervals.txt`: depth\_top / depth\_bottom / duration / unc. (absolute) for air intervals for an
ice core
- `density.txt`: depth / relative density for an ice core
- `lock_in_depth.txt`: depth / Lock-in-Depth / unc. (opt, relative) for an ice core
- `thinning.txt`: depth / thinning function / unc. (opt, relative) for an ice core
- `delta_depths.txt`: depth / Delta-depth / unc. (absolute) for an ice core

Then you have one directory per site pair, which contains:
- `parameters_covariance_observations.py`: this file allows to define the correlation of site pair
specific observations
- `
- `ice_depth.txt`: depth1 / depth2 / age unc. (absolute) for ice-ice stratigraphic links observations
- `air_depth.txt`: depth1 / depth2 / age unc. (absolute) for air-air stratigraphic links observations
- `iceair_depth.txt`: depth1 / depth2 / age unc. (absolute) for ice-air stratigraphic links observations
- `airice_depth.txt`: depth1 / depth2 / age unc. (absolute) for air-ice stratigraphic links observations

A few things you need to know to use paleochrono:
1) You can use whatever units you want but they need to be consistent. For example, if you use meters for the depths and years for the dated horizons, you need to use meters per years for the accumulation rates. 
2) The site specific parameters override the general parameters for all sites. In the very same way, the site-pair specific parameters override the general parameters for all site-pairs.
3) The standard deviations defined in the parameters-Covariance*.py override the standard deviation defined in the observation or prior files.
4) Most of these files are optional. If there is no file for an certain type of observations, that means that there is no observation of this type. If a covariance matrix is not defined for a prior or an observation type, that means that the correlation matrix is supposed to be equal to identity and that the standard deviation is given in the prior or observation file.


What is the structure of the general `parameters.yml` file?
--------------------------------------------------------

It contains the list of sites, the optimization method to be used and some settings for the figures.
It is where you define the names of your sites.
Have a look at the file `AICC2023-Hulu/parameters.yml`, it is commented.


What is the structure of a site `parameters.yml` file?
---------------------------------------------------------

It defines age at the top of the core, the unthinned depth at the top of the core, the age equation grid, the correction functions grids and the type of representation of the prior accu scenario (linear or staircase). You can also define other parameters that are used to defined the covariance matrices of the priors.
Have a look at the files `AICC2023-Hulu/EDC/parameters.yml`, it is commented.

How do I use the C14 calibration?
---------------------------------
First, you need to install iosacal, which is a python software for C14 calibration.
isoacal is included in conda-forge, so you have to add the conda-forge repo.
If you are using Anaconda, you can add conda-forge like this:
```
conda config --append channels conda-forge
```
Then you can install iosacal like this:
```
conda install iosacal
```
To change the C14 calibration curve, you can set the `c14_cal` key in the main parameters.yml file.
Please have a look at the iosacal doc to have a list of C14 calibration curves.

How to set up the `parameters-CovarianceObservations.py` file?
--------------------------------------------------------------

You need to know a little bit of python to do that.
Feel free to send an email on the mailing list if you need assistance.

For site specific observations, you set up the correlation matrices in the file `parameters-CovarianceObservations.py` in the site directory.
- `self.icemarkers_correlation`     : for ice dated horizons
- `self.airmarkers_correlation`     : for air dated horizons
- `self.iceintervals_correlation`   : for ice dated intervals
- `self.airintervals_correlation`   : for air dated intervals
- `self.Ddepth_correlation`         : for Delta-depth observations

For site pair specific observations (stratigraphic links), you set up the correlation matrices in the file `parameters-CovarianceObservations.py` in the site pair directory.
- `self.iceicemarkers_correlation`  : for ice-ice stratigraphic links
- `self.airairmarkers_correlation`  : for air-air stratigraphic links
- `self.iceairmarkers_correlation`  : for ice-air stratigraphic links
- `self.airicemarkers_correlation`  : for air-ice stratigraphic links

Let us take a concrete example and assume we want a correlation matrix for ice dated horizons with ones in the diagonal and with a constant correlation factor k outside the diagonal, you can write:

```
self.icemarkers_correlation=k*np.ones((np.shape(self.icemarkers_correlation)))+(1-k)*np.diag(np.ones(np.shape(self.icemarkers_correlation)[0]))
```

Don't forget that if you find the use of python and the Paleochrono internal variables too difficult, you can define your correlation matrices outside paleochrono and import them here by using for example the `np.loadtxt` function.


What to do if something goes wrong?
-----------------------------------

Some errors can be eliminated by restarting the kernel in spyder (under "Console">"Restart kernel").
If the problem persist, please post an email to the author or on the mailing list with the error message appearing on the command line.
