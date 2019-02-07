#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Feb  7 11:52:46 2019
Global variables used across all modules.
@author: parrenif
"""
import sys

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
print 'Parameters directory is: ', DATADIR
#os.chdir(DATADIR)

execfile(DATADIR+'/parameters.py')

try:
    LIST_SITES = list_drillings
except NameError:
    pass
try:
    OPT_METHOD = opt_method
except NameError:
    pass
try:
    NB_NODES = nb_nodes
except NameError:
    pass
try:
    COLOR_OBS = color_obs
except NameError:
    pass
try:
    COLOR_OPT = color_opt
except NameError:
    pass
try:
    COLOR_MOD = color_mod
except NameError:
    pass
try:
    COLOR_CI = color_ci
except NameError:
    pass
try:
    COLOR_SIGMA = color_sigma
except NameError:
    pass
try:
    COLOR_DI = color_di
except NameError:
    pass
try:
    SHOW_INITIAL = show_initial
except NameError:
    pass
try:
    COLOR_INIT = color_init
except NameError:
    pass
try:
    SCALE_AGECI = scale_ageci
except NameError:
    pass
try:
    SHOW_FIGURES = show_figures
except NameError:
    pass
try:
    SHOW_AIRLAYERTHICK = show_airlayerthick
except NameError:
    pass