# -*- coding: utf-8 -*-
"""
Created on Tue Oct 14 18:32:46 2014

@author: noore
"""
import json
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

HTML_FNAME = './examples/res_ctherm_ref_tca/'

REACTION_FNAMES = ['./examples/cterm_112418_TCA_ethanol.txt',
                   #'../examples/cterm_022318_P1_ATPonly.txt',
                   './examples/cterm_112418_TCA.txt',
                   #'../examples/cterm_022318_P3_PYK.txt',
                   #'../examples/cterm_022318_P4_PFK_PYK.txt',
                   #'../examples/cterm_022318_P5_MalShunt.txt',
                   #'../examples/cterm_022318_P6_H2.txt',
                   #'../examples/cterm_022318_P7_AdhNADP.txt',
                   ]
saveDirs = ["P1_TCA_ethanol_upreg",
            "P1_TCA_only_upreg",
            #"P2_ATP_PFK_ATPratio",
#            "P3_PYK",
#            "P4_PFK_PYK",
#            "P5_MalShunt",
#            "P6_H2",
#            "P7_AdhNADP",
            ]
atp_ratio = [12,13,14,15,16]#[25.51]#[1, 5, 10, 20]

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
ref_conc0 = {'C00469': 10,
            'C00004': 0.08,
            'C00024': 0.83,
            'C00002': 2.70,
            'C00008': 0.11,
            'C00020': 0.22,
            'C00354': 1.50,
            'C00092': 8.19,
            'C00085': 1.49,
            'C00103': 6.66,
            'C00118': 0.10,
            'C00197': 1.35,
            'C00011': 1.27,
            'C00074': 0.69,
            'C00005': 0.39,
            'C00003': 2.25,
            'C00006': 0.26,
            'C00010': 0.02,
            'C00022': 12.66,
            'C00103': 6.66
            }
'''ref_conc0 = {'C00469': 10,
            'C00004': 0.08,
            'C00024': 0.83,
            'C00002': 2.51,
            'C00008': 0.11,
            'C00020': 0.22,
            'C00354': 1.50,
            'C00092': 8.19,
            'C00074': 0.69,
            'C00005': 0.38,
            'C00022': 12.65,
            'C00103': 6.66
            }'''

'''ref_conc0 = {'C00469': 10,
            'C00004': 0.3,
            'C00024': 0.75,
            'C00002': 13.55,
            'C00008': 0.7,
            'C00020': 1.98,
            'C00354': 4.77,
            'C00092': 23.75,
            'C00074': 0.14,
            #'C00022': 0.074,
            }'''


for k in range(len(saveDirs)):
    pathways = KeggFile2ModelList(REACTION_FNAMES[k])
    p = pathways[0]
    cc = ComponentContribution.init()
    ratio = 0.03
    p['model'].add_thermo(cc)

    #p['bounds']['C00003'] = (2.25e-4,) * 2
    #p['bounds']['C00004'] = (0.8e-5,) * 2

    #p['bounds']['C00008'] = (1e-4,) * 2
    #if k == 0:
        #p['bounds']['C00035'] = (1e-4,) * 2

    cid = p['model'].cids

    ref_conc = {k: float(v) / 1000 for k,v in ref_conc0.iteritems()}
    for cpd, conc in ref_conc.iteritems():
        p['bounds'][cpd] = (conc, ) * 2
    for r in atp_ratio:

        dictCur = {}
        # fix initial ATP/ADP ratio
        #p['bounds']['C00002'] = (1e-4 * r, ) * 2

        # fix GTP/GDP ratio
        #p['bounds']['C00044'] = (1e-4 * r, ) * 2

        ref = True
#        ref_conc = {}
        t = 0

        print "ATP & GTP, time point %d" %t
            #else:
             #   print "ATP only, time point %d" %t
        p['bounds']['C00158'] = (1e-6 * r, ) * 2
        print p["bounds"]
        html_writer = HtmlWriter(HTML_FNAME + saveDirs[k] + "/atp%d_t%d.html" %(r, t))
        mdf = MaxMinDrivingForce(p['model'], p['fluxes'], p['bounds'],
        pH=p['pH'], I=p['I'], T=p['T'],
        html_writer=html_writer)


        mdf_solution, dG_r_prime, param = mdf.Solve(uncertainty_factor=3.0)
            #plt.show()
            #t += 1

        if ref:
            ref = False
            conc_list = [x[0] for x in param['concentrations'].tolist()]
            ref_conc_comput = dict(zip(cid, conc_list))
            for c in cid:
                if c in mets_constr and c not in ref_conc:
                    ref_conc[c] = ref_conc_comput[c]

                #print ref_conc

