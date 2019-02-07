#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Feb  7 16:55:54 2019
Module for the SitePair class.
@author: parrenif
"""

import os
import numpy as np
import matplotlib.pyplot as mpl
from matplotlib.backends.backend_pdf import PdfPages
from scipy.linalg import lu_factor, lu_solve
from scipy.linalg import cholesky
import pccfg

class SitePair(object):
    """Class for a pair of sites."""

    def __init__(self, site1, site2):
        self.site1 = site1
        self.site2 = site2
        self.label = self.site1.label+'-'+self.site2.label


#TODO: allow to have either dlabel1+'-'dlabel2 or dlbel2+'-'dlabel1 as directory
        filename = pccfg.DATADIR+self.site1.label+'-'+self.site2.label+'/ice_depth.txt'
        if os.path.isfile(filename) and open(filename).read():
            readarray = np.loadtxt(filename)
            self.iceicemarkers_depth1 = readarray[:, 0]
            self.iceicemarkers_depth2 = readarray[:, 1]
            self.iceicemarkers_sigma = readarray[:, 2]
        else:
            self.iceicemarkers_depth1 = np.array([])
            self.iceicemarkers_depth2 = np.array([])
            self.iceicemarkers_sigma = np.array([])

        filename = pccfg.DATADIR+self.site1.label+'-'+self.site2.label+'/air_depth.txt'
        if os.path.isfile(filename) and open(filename).read():
            readarray = np.loadtxt(filename)
            self.airairmarkers_depth1 = readarray[:, 0]
            self.airairmarkers_depth2 = readarray[:, 1]
            self.airairmarkers_sigma = readarray[:, 2]
        else:
            self.airairmarkers_depth1 = np.array([])
            self.airairmarkers_depth2 = np.array([])
            self.airairmarkers_sigma = np.array([])

        filename = pccfg.DATADIR+self.site1.label+'-'+self.site2.label+'/iceair_depth.txt'
        if os.path.isfile(filename) and open(filename).read():
            readarray = np.loadtxt(filename)
            self.iceairmarkers_depth1 = readarray[:, 0]
            self.iceairmarkers_depth2 = readarray[:, 1]
            self.iceairmarkers_sigma = readarray[:, 2]
        else:
            self.iceairmarkers_depth1 = np.array([])
            self.iceairmarkers_depth2 = np.array([])
            self.iceairmarkers_sigma = np.array([])

        filename = pccfg.DATADIR+self.site1.label+'-'+self.site2.label+'/airice_depth.txt'
        if os.path.isfile(filename) and open(filename).read():
            readarray = np.loadtxt(filename)
            self.airicemarkers_depth1 = readarray[:, 0]
            self.airicemarkers_depth2 = readarray[:, 1]
            self.airicemarkers_sigma = readarray[:, 2]
        else:
            self.airicemarkers_depth1 = np.array([])
            self.airicemarkers_depth2 = np.array([])
            self.airicemarkers_sigma = np.array([])


        self.iceicemarkers_correlation = np.diag(np.ones(np.size(self.iceicemarkers_depth1)))
        self.airairmarkers_correlation = np.diag(np.ones(np.size(self.airairmarkers_depth1)))
        self.iceairmarkers_correlation = np.diag(np.ones(np.size(self.iceairmarkers_depth1)))
        self.airicemarkers_correlation = np.diag(np.ones(np.size(self.airicemarkers_depth1)))
        filename = pccfg.DATADIR+'/parameters-CovarianceObservations-AllSitePairs.py'
        if os.path.isfile(filename):
            execfile(filename)
        filename = pccfg.DATADIR+self.label+'/parameters-CovarianceObservations.py'
        if os.path.isfile(filename):
            execfile(filename)
        if np.size(self.iceicemarkers_depth1) > 0:
            self.iceicemarkers_chol = cholesky(self.iceicemarkers_correlation)
            self.iceicemarkers_lu_piv = lu_factor(self.iceicemarkers_chol)
        if np.size(self.airairmarkers_depth1) > 0:
            self.airairmarkers_chol = cholesky(self.airairmarkers_correlation)
            self.airairmarkers_lu_piv = lu_factor(self.airairmarkers_chol)
        if np.size(self.iceairmarkers_depth1) > 0:
            self.iceairmarkers_chol = cholesky(self.iceairmarkers_correlation)
            self.iceairmarkers_lu_piv = lu_factor(self.iceairmarkers_chol)
        if np.size(self.airicemarkers_depth1) > 0:
            self.airicemarkers_chol = cholesky(self.airicemarkers_correlation)
            self.airicemarkers_lu_piv = lu_factor(self.airicemarkers_chol)


    def residuals(self):
        """Calculate the residual terms of a pair of sites."""

        resi_iceice = (self.site1.fct_age(self.iceicemarkers_depth1)-\
                       self.site2.fct_age(self.iceicemarkers_depth2))/self.iceicemarkers_sigma
        if np.size(self.iceicemarkers_depth1) > 0:
            resi_iceice = lu_solve(self.iceicemarkers_lu_piv, resi_iceice)
        resi_airair = (self.site1.fct_airage(self.airairmarkers_depth1)-\
                       self.site2.fct_airage(self.airairmarkers_depth2))/self.airairmarkers_sigma
        if np.size(self.airairmarkers_depth1) > 0:
            resi_airair = lu_solve(self.airairmarkers_lu_piv, resi_airair)
        resi_iceair = (self.site1.fct_age(self.iceairmarkers_depth1)-\
                       self.site2.fct_airage(self.iceairmarkers_depth2))/self.iceairmarkers_sigma
        if np.size(self.iceairmarkers_depth1) > 0:
            resi_iceair = lu_solve(self.iceairmarkers_lu_piv, resi_iceair)
        resi_airice = (self.site1.fct_airage(self.airicemarkers_depth1)-\
                       self.site2.fct_age(self.airicemarkers_depth2))/self.airicemarkers_sigma
        if np.size(self.airicemarkers_depth1) > 0:
            resi_airice = lu_solve(self.airicemarkers_lu_piv, resi_airice)
        resi = np.concatenate((resi_iceice, resi_airair, resi_iceair, resi_airice))

        return resi


    def figures(self):
        """Build the figures related to a pair of sites."""

        if not os.path.isdir(pccfg.DATADIR+self.label):
            os.mkdir(pccfg.DATADIR+self.label)


        mpl.figure(self.label+' ice-ice')
        mpl.xlabel(self.site1.label+' ice age (yr b1950)')
        mpl.ylabel(self.site2.label+' ice age (yr b1950)')
        if np.size(self.iceicemarkers_depth1) > 0:
            if pccfg.SHOW_INITIAL:
                mpl.errorbar(self.site1.fct_age_init(self.iceicemarkers_depth1),
                             self.site2.fct_age_init(self.iceicemarkers_depth2), color=pccfg.COLOR_INIT,
                             xerr=self.iceicemarkers_sigma, linestyle='', marker='o', markersize=2,
                             label="Initial")
            mpl.errorbar(self.site1.fct_age_model(self.iceicemarkers_depth1),
                         self.site2.fct_age_model(self.iceicemarkers_depth2), color=pccfg.COLOR_MOD,
                         xerr=self.iceicemarkers_sigma, linestyle='', marker='o', markersize=2,
                         label="Prior")
            mpl.errorbar(self.site1.fct_age(self.iceicemarkers_depth1),
                         self.site2.fct_age(self.iceicemarkers_depth2), color=pccfg.COLOR_OPT,
                         xerr=self.iceicemarkers_sigma, linestyle='', marker='o', markersize=2,
                         label="Posterior")
        x_low, x_up, y_low, y_up = mpl.axis()
        x_low = self.site1.age_top
        y_low = self.site2.age_top
        mpl.axis((x_low, x_up, y_low, y_up))
        rangefig = np.array([max(x_low, y_low), min(x_up, y_up)])
        mpl.plot(rangefig, rangefig, color=pccfg.COLOR_OBS, label='perfect agreement')
        mpl.legend(loc="best")
        printed_page = PdfPages(pccfg.DATADIR+self.label+'/ice-ice.pdf')
        printed_page.savefig(mpl.figure(self.label+' ice-ice'))
        printed_page.close()
        if not pccfg.SHOW_FIGURES:
            mpl.close()

        mpl.figure(self.label+' air-air')
        mpl.xlabel(self.site1.label+' air age (yr b1950)')
        mpl.ylabel(self.site2.label+' air age (yr b1950)')
        if np.size(self.airairmarkers_depth1) > 0:
            if pccfg.SHOW_INITIAL:
                mpl.errorbar(self.site1.fct_airage_init(self.airairmarkers_depth1),
                             self.site2.fct_airage_init(self.airairmarkers_depth2),
                             color=pccfg.COLOR_INIT, xerr=self.airairmarkers_sigma, linestyle='',
                             marker='o', markersize=2, label="Initial")
            mpl.errorbar(self.site1.fct_airage_model(self.airairmarkers_depth1),
                         self.site2.fct_airage_model(self.airairmarkers_depth2), color=pccfg.COLOR_MOD,
                         xerr=self.airairmarkers_sigma, linestyle='', marker='o', markersize=2,
                         label="Prior")
            mpl.errorbar(self.site1.fct_airage(self.airairmarkers_depth1),
                         self.site2.fct_airage(self.airairmarkers_depth2), color=pccfg.COLOR_OPT,
                         xerr=self.airairmarkers_sigma, linestyle='', marker='o', markersize=2,
                         label="Posterior")
        x_low, x_up, y_low, y_up = mpl.axis()
        x_low = self.site1.age_top
        y_low = self.site2.age_top
        mpl.axis((x_low, x_up, y_low, y_up))
        rangefig = np.array([max(x_low, y_low), min(x_up, y_up)])
        mpl.plot(rangefig, rangefig, color=pccfg.COLOR_OBS, label='perfect agreement')
        mpl.legend(loc="best")
        printed_page = PdfPages(pccfg.DATADIR+self.label+'/air-air.pdf')
        printed_page.savefig(mpl.figure(self.label+' air-air'))
        printed_page.close()
        if not pccfg.SHOW_FIGURES:
            mpl.close()

        mpl.figure(self.label+' ice-air')
        mpl.xlabel(self.site1.label+' ice age (yr b1950)')
        mpl.ylabel(self.site2.label+' air age (yr b1950)')
        if np.size(self.iceairmarkers_depth1) > 0:
            if pccfg.SHOW_INITIAL:
                mpl.errorbar(self.site1.fct_age_init(self.iceairmarkers_depth1),
                             self.site2.fct_airage_init(self.iceairmarkers_depth2),
                             color=pccfg.COLOR_INIT, xerr=self.iceairmarkers_sigma, linestyle='',
                             marker='o', markersize=2, label="Initial")
            mpl.errorbar(self.site1.fct_age_model(self.iceairmarkers_depth1),
                         self.site2.fct_airage_model(self.iceairmarkers_depth2), color=pccfg.COLOR_MOD,
                         xerr=self.iceairmarkers_sigma, linestyle='', marker='o', markersize=2,
                         label="Prior")
            mpl.errorbar(self.site1.fct_age(self.iceairmarkers_depth1),
                         self.site2.fct_airage(self.iceairmarkers_depth2), color=pccfg.COLOR_OPT,
                         xerr=self.iceairmarkers_sigma, linestyle='', marker='o', markersize=2,
                         label="Posterior")
        x_low, x_up, y_low, y_up = mpl.axis()
        x_low = self.site1.age_top
        y_low = self.site2.age_top
        mpl.axis((x_low, x_up, y_low, y_up))
        rangefig = np.array([max(x_low, y_low), min(x_up, y_up)])
        mpl.plot(rangefig, rangefig, color=pccfg.COLOR_OBS, label='perfect agreement')
        mpl.legend(loc="best")
        printed_page = PdfPages(pccfg.DATADIR+self.label+'/ice-air.pdf')
        printed_page.savefig(mpl.figure(self.label+' ice-air'))
        printed_page.close()
        if not pccfg.SHOW_FIGURES:
            mpl.close()

        mpl.figure(self.label+' air-ice')
        mpl.xlabel(self.site1.label+' air age (yr b1950)')
        mpl.ylabel(self.site2.label+' ice age (yr b1950)')
        if np.size(self.airicemarkers_depth1) > 0:
            if pccfg.SHOW_INITIAL:
                mpl.errorbar(self.site1.fct_airage_init(self.airicemarkers_depth1),
                             self.site2.fct_age_init(self.airicemarkers_depth2),
                             color=pccfg.COLOR_INIT, xerr=self.airicemarkers_sigma,
                             linestyle='', marker='o', markersize=2, label="Initial")
            mpl.errorbar(self.site1.fct_airage_model(self.airicemarkers_depth1),
                         self.site2.fct_age_model(self.airicemarkers_depth2), color=pccfg.COLOR_MOD,
                         xerr=self.airicemarkers_sigma, linestyle='', marker='o', markersize=2,
                         label="Prior")
            mpl.errorbar(self.site1.fct_airage(self.airicemarkers_depth1),
                         self.site2.fct_age(self.airicemarkers_depth2), color=pccfg.COLOR_OPT,
                         xerr=self.airicemarkers_sigma, linestyle='', marker='o', markersize=2,
                         label="Posterior")
        x_low, x_up, y_low, y_up = mpl.axis()
        x_low = self.site1.age_top
        y_low = self.site2.age_top
        mpl.axis((x_low, x_up, y_low, y_up))
        rangefig = np.array([max(x_low, y_low), min(x_up, y_up)])
        mpl.plot(rangefig, rangefig, color=pccfg.COLOR_OBS, label='perfect agreement')
        mpl.legend(loc="best")
        printed_page = PdfPages(pccfg.DATADIR+self.label+'/air-ice.pdf')
        printed_page.savefig(mpl.figure(self.label+' air-ice'))
        printed_page.close()
        if not pccfg.SHOW_FIGURES:
            mpl.close()
