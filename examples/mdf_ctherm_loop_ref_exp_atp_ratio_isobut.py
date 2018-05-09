# mdf_ctherm_loop_ref_exp_atp_ratio_isobut.py
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

HTML_FNAME = '../examples/res_ctherm_ref_exp_max10/'

REACTION_FNAMES = [#'../examples/cterm_022318_P1.txt',
                   #'../examples/cterm_022318_P1_ATPonly.txt',
                   '../examples/cterm_050818_P1_isobut.txt',
                   #'../examples/cterm_050818_P1_isobut_nadph.txt',
                   #'../examples/cterm_050818_P1_isobut_lump_nadh.txt',
                   #'../examples/cterm_050818_P1_isobut_lump_nadph.txt',
                   #'../examples/cterm_022318_P3_PYK.txt',
                   #'../examples/cterm_022318_P4_PFK_PYK.txt',
                   #'../examples/cterm_022318_P5_MalShunt.txt',
                   #'../examples/cterm_022318_P6_H2.txt',
                   #'../examples/cterm_022318_P7_AdhNADP.txt',
                   ]
saveDirs = [#"P1_ref_path",
            #"P1_ref_path_noGTP",
            #"P1_isobut",
            #"P1_isobut_nadph",
            #"P1_isobut_lump_nadh",
            #"P1_isobut_lump_nadph",
            "test",
#            "P3_PYK",
#            "P4_PFK_PYK",
#            "P5_MalShunt",
#            "P6_H2",
#            "P7_AdhNADP",
            ]
atp_ratio = [20] #[1, 5, 10, 20]

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

nadphRatio = 1
ref_conc0 = {'C14710': 10,  # replace ethanol by isobutanol
            'C00004': 0.3,
            'C00024': 0.75,
            'C00002': 13.55,
            'C00008': 0.7,
            'C00020': 1.98,
            'C00354': 4.77,
            'C00092': 23.75,
            'C00074': 0.14,
            'C00005': nadphRatio * 1e-4,
            #'C00022': 0.074,
            }

downstreamMets = ['C06010', 'C04272', 'C00141', 'C80061']

nadphRatio2Test = [0.1]#[0.1, 1, 10]
downstreamMaxConc = [1e-5]#[1e-5, 1e-6, 1e-4, 1e-3, 1e-2]

for ratio in nadphRatio2Test:
    for q in range(len(downstreamMaxConc)):
        for k in range(len(saveDirs)):
            if q == 0 or "lump" not in saveDirs[k]:
                pathways = KeggFile2ModelList(REACTION_FNAMES[k])
                p = pathways[0]
                cc = ComponentContribution.init()
                
                p['model'].add_thermo(cc)
                
                p['bounds']['C00003'] = (1e-4,) * 2
                p['bounds']['C00004'] = (0.3e-4,) * 2
                p['bounds']['C00006'] = (1e-4,) * 2
                
                p['bounds']['C00008'] = (1e-4,) * 2
             
                #if k == 0:
                p['bounds']['C00035'] = (1e-4,) * 2
                
                if k <= 1:
                    for m in downstreamMets:
                        p['bounds'][m] = (1e-6, downstreamMaxConc[q])

                cid = p['model'].cids
                
                ref_conc = {k: float(v) / 1000 for k,v in ref_conc0.iteritems()}
                
                for r in atp_ratio:
                    
                    # fix initial ATP/ADP ratio        
                    p['bounds']['C00002'] = (1e-4 * r, ) * 2
                    
                    # fix GTP/GDP ratio
                    p['bounds']['C00044'] = (1e-4 * r, ) * 2

                    ref = True
            #        ref_conc = {}
                    t = 0
                
                    for d in exp_data:
                        for cpd, conc in ref_conc.iteritems():
                            if cpd == "C14710":
                                # absolute EtOH conc
                                ## replace ethanol by isobutanol
                                #p['bounds'][cpd] = (d[cpd], ) * 2
                                p['bounds'][cpd] = (d['C00469'], ) * 2
                            elif cpd in ["C00004"]:
                                # NADH/NAD+ ratio
                                p['bounds'][cpd] = (1e-4 * float(d[cpd]), ) * 2
                            #elif cpd == "C00008" and t == 0:
                                #p['bounds'][cpd] = (conc, ) * 2
                            elif cpd in ['C00005']:
                                # varying NADPH/NADPratios
                                # p['bounds'][cpd] = (ratio * 1e-4 * float(d[cpd]), ) * 2
                                # fixedOrInit = 'Init'
                                # fixed NADPH/NADP ratio
                                p['bounds'][cpd] = (ratio * 1e-4, ) * 2
                                fixedOrInit = 'Fix'
                            elif cpd not in ["C00008"]:
                                # other relative data
                                p['bounds'][cpd] = (conc * d[cpd], ) * 2
                            
                            
                        
                        #if k == 0:
                        print "ATP & GTP, time point %d" %t
                        #else:
                         #   print "ATP only, time point %d" %t
                            
                        print p["bounds"]    
                        saveNameCur = HTML_FNAME + saveDirs[k] + "/atp%d_ratio%s%.0e" %(r, fixedOrInit, ratio)
                        if k <= 1:
                            saveNameCur += "_dsMax%.0e" %downstreamMaxConc[q]

                        html_writer = HtmlWriter(saveNameCur + "_t%d.html" %t)
                        mdf = MaxMinDrivingForce(p['model'], p['fluxes'], p['bounds'],
                                             pH=p['pH'], I=p['I'], T=p['T'],
                                             html_writer=html_writer)
                        
                        
                        mdf_solution, dG_r_prime, param = mdf.Solve(uncertainty_factor=3.0)
                        #plt.show()
                        t += 1
                        
                        if ref:
                            ref = False
                            conc_list = [x[0] for x in param['concentrations'].tolist()]
                            ref_conc_comput = dict(zip(cid, conc_list))
                            for c in cid:
                                if c in mets_constr and c not in ref_conc:
                                    ref_conc[c] = ref_conc_comput[c]
                            
                            print ref_conc
