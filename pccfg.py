#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Feb  7 11:52:46 2019
Global variables used across all modules.
@author: parrenif
"""
import sys
import yaml

##Default Parameters
LIST_SITES = []
OPT_METHOD = 'none'  #leastsq, leastsq-parallel, none
NB_NODES = 6         #Number of nodes for the leastsq-parallel mode
DATADIR='./'
COLOR_OBS = 'r'       #color for the observations
COLOR_OPT = 'k'       #color for the posterior scenario
COLOR_MOD = 'b'       #color for the prior scenario
COLOR_CI = '0.8'      #color for the confidence intervals
COLOR_SIGMA = 'm'     #color for the uncertainty
COLOR_DI = 'g'        #color for the dated intervals
SHOW_INITIAL = False  #always put to False for now
COLOR_INIT = 'c'      #always put to 'c' for now
SCALE_AGECI = 10.     #scaling of the confidence interval in the ice and air age figures
SHOW_FIGURES = False  #whether to show or not the figures at the end of the run
SHOW_AIRLAYERTHICK = False #whether to show the air layer thickness figure (buggy on anaconda)

###Reading parameters directory
DATADIR = sys.argv[1]
if DATADIR[-1] != '/':
    DATADIR = DATADIR+'/'
print('Parameters directory is: ', DATADIR)
#os.chdir(DATADIR)

with open(DATADIR+"parameters.yml", "r") as f:
    data = yaml.load(f)
    globals().update(data)