# -*- coding: utf-8 -*-
"""
Created on Tue Oct 14 18:32:46 2014

@author: noore
"""
import sys,pdb
#sys.path.append("../")
from scripts.max_min_driving_force import KeggFile2ModelList, MaxMinDrivingForce
from component_contribution.component_contribution_trainer import ComponentContribution
from scripts.html_writer import HtmlWriter
import logging
import numpy as np
import matplotlib.pyplot as plt


# REACTION_FNAME = '../examples/cterm_butanol.txt'
# HTML_FNAME = '../res/mdf_cterm_butanol.html'

HTML_FNAME = './examples/res_ctherm/'

REACTION_FNAMES = ['./examples/cterm_022318_P1.txt',
                   './examples/cterm_022318_P1_ATPonly.txt',
                   ]
saveDirs = ["ATPADPfree",
            "ATPADPfree_noGTP"
            ]
atp_ratio = [1, 5, 10, 20]

exp_data_file = "ctherm_exp_data.txt"
exp_data = []
f = open(exp_data_file, "r")
l = f.readline()
l = f.readline()
mets_constr = [i for i in l.strip("\n").split(" ") if i != ""]
l = f.readline()
while l != "":
    data = [float(i) for i in l.strip("\n").split(" ") if i != ""]
    exp_data.append(dict(zip(mets_constr, data)))
    l = f.readline()

f.close()



for k in range(len(saveDirs)):
    pathways = KeggFile2ModelList(REACTION_FNAMES[k])
    p = pathways[0]
    cc = ComponentContribution.init()

    p['model'].add_thermo(cc)

    p['bounds']['C00003'] = (1e-2,) * 2
    p['bounds']['C00004'] = (0.3e-2,) * 2

    p['bounds']['C00008'] = (1e-4,) * 2
    if k == 0:
        p['bounds']['C00035'] = (1e-4,) * 2

    cid = p['model'].cids

    for r in atp_ratio:

        # fix ATP/ADP ratio
        p['bounds']['C00002'] = (1e-4 * r, ) * 2
        if k == 0:
            # fix GTP/GDP ratio
            p['bounds']['C00044'] = (1e-4 * r, ) * 2

        ref = True
        ref_conc = {}
        t = 0

        for d in exp_data:
            if not ref:
                for cpd, conc in ref_conc.iteritems():
                    if cpd == "C00469":
                        # absolute EtOH conc
                        p['bounds'][cpd] = (d[cpd], ) * 2
                    elif cpd == "C00004":
                        # NADH/NAD+ ratio
                        p['bounds'][cpd] = (1e-2 * float(d[cpd]), ) * 2
                    else:
                        # other relative data
                        p['bounds'][cpd] = (conc * d[cpd], ) * 2

            html_writer = HtmlWriter(HTML_FNAME + saveDirs[k] + "/atp%d_t%d.html" %(r, t))
            mdf = MaxMinDrivingForce(p['model'], p['fluxes'], p['bounds'],
                                 pH=p['pH'], I=p['I'], T=p['T'],
                                 html_writer=html_writer)


            mdf_solution, dG_r_prime, param = mdf.Solve(uncertainty_factor=3.0)
            #plt.show()
            t += 1

            if ref:
                ref = False
                conc_list = [x[0] for x in param['concentrations'].tolist()]
                ref_conc = dict(zip(cid, conc_list))
                ref_conc = {k:v for k, v in ref_conc.iteritems() if k in mets_constr}
