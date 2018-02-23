# -*- coding: utf-8 -*-
"""
Created on Tue Oct 14 18:32:46 2014

@author: noore
"""
import sys,pdb
sys.path.append("../")
from scripts.max_min_driving_force import KeggFile2ModelList, MaxMinDrivingForce
from component_contribution.component_contribution_trainer import ComponentContribution
from scripts.html_writer import HtmlWriter
import logging
import numpy as np
import matplotlib.pyplot as plt


# REACTION_FNAME = '../examples/cterm_butanol.txt'
# HTML_FNAME = '../res/mdf_cterm_butanol.html'
REACTION_FNAME = '../examples/cterm_022318_P1.txt'
HTML_FNAME = '../examples/res_ctherm/mdf_cterm_test.html'

exp_data_file = "ctherm_exp_data.txt"
exp_data = []
f = open(exp_data_file, "r")
l = f.readline()
l = f.readline()
mets = [i for i in l.strip("\n").split(" ") if i != ""]
l = f.readline()
while l != "":
    data = [float(i) for i in l.strip("\n").split(" ") if i != ""]
    exp_data.append(dict(zip(mets, data)))
    l = f.readline()

f.close()

html_writer = HtmlWriter(HTML_FNAME)
pathways = KeggFile2ModelList(REACTION_FNAME)
p = pathways[0]
cc = ComponentContribution.init()

p['model'].add_thermo(cc)

p['bounds']['C00003'] = (1e-2,) * 2
p['bounds']['C00004'] = (0.3e-2,) * 2

cid = p['model'].cids
ref = True
ref_conc = {}
for d in exp_data:
    if not ref:
        for cpd, conc in ref_conc.iteritems():
            if cpd == "C00469":
                p['bounds'][cpd] = (d[cpd], ) * 2
            elif cpd == "C00004":
                p['bounds'][cpd] = (1e-2 * float(d[cpd]), ) * 2
            else:
                p['bounds'][cpd] = (conc * d[cpd])
        
    mdf = MaxMinDrivingForce(p['model'], p['fluxes'], p['bounds'],
                         pH=p['pH'], I=p['I'], T=p['T'],
                         html_writer=html_writer)

mdf_solution, dG_r_prime, param = mdf.Solve(uncertainty_factor=3.0)
plt.show()
