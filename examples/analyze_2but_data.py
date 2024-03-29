import json
import matplotlib
import numpy as np
import matplotlib.cm as cm
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt

f = open("./examples/res_ctherm_ref_exp_max10_2but/mdf_2but.json", "r")
s = f.read()
data = json.loads(s)
f.close()

ibutConc = []
f = open("ctherm_exp_data.txt", "r")
f.readline()
f.readline()
line = f.readline()
while line != "":
	ibutConc.append(float(line.split(" ")[0]))
	line = f.readline()

nadphRatio2Test = [4.9]#, 0.3, 1, 3, 10]#[0.1, 0.2, 0.5, 1, 2, 5, 10]
downstreamMaxConc = [1e-6, 2e-6, 5e-6, 1e-5, 2e-5, 5e-5, 1e-4, 2e-4, 5e-4, 1e-3, 2e-3, 5e-3, 1e-2]

p, q = len(nadphRatio2Test), len(downstreamMaxConc)

dsRxnIDs = ['ALS', 'ALDC', '23BDH', 'DDHT', 'SADH']
# get the shadow prices of the downstream reactions as a function of isobutanol conc, NADPH/NADP ratios and max concentrations of downstream mets
keys = ["nadph", "maxConc", "ibutConc", "mdf", "atp", "time", "PFK", "GAPDH"] + dsRxnIDs
data2 = {}
for k1, v1 in data.iteritems():
	# loop over different simulated pathways
	data2[k1] = {k: [] for k in keys}
	KARI = "KARI_NADH" if k1 == "P1_isobut" else "KARI_NADPH"
	for k2, v2 in v1.iteritems():
		# loop over different nadph/nadp ratios
		for k3, v3 in v2.iteritems():
			# loop over maximum concentrations allowed for downstream mets
			for k4, v4 in v3.iteritems():
				# loop over ATP/ADP ratios
				for k5, v5 in v4.iteritems():
					# loop over different time points (thus different isobutanol concentrations)
					data2[k1]["nadph"].append(float(k2))
					data2[k1]["maxConc"].append(float(k3))
					data2[k1]["atp"].append(float(k4))
					data2[k1]["time"].append(float(k5))
					data2[k1]["ibutConc"].append(ibutConc[int(k5)])
					data2[k1]["mdf"].append(v5["mdf"])
					data2[k1]["PFK"].append(v5["reaction prices"]["PPI_PFK"])
					data2[k1]["GAPDH"].append(v5["reaction prices"]["GAPDH"])
					for dsR in dsRxnIDs:
						data2[k1][dsR].append(v5["reaction prices"][dsR])

data3 = {}
for k1, v1 in data.iteritems():
	data3[k1] = {}
	for t in range(10):
		data3[k1][t] = {k: [] for k in keys}
		for i in nadphRatio2Test:
			listTmp = {k: [] for k in keys}
			for j in downstreamMaxConc:
				index = [k for k in range(p * q * 10) if  data2[k1]["time"][k] == t and data2[k1]["nadph"][k] == i and data2[k1]["maxConc"][k] == j][0]
				for k in keys:
					listTmp[k].append(data2[k1][k][index])

			for k in keys:
				data3[k1][t][k].append(listTmp[k])

		for k in keys:
			data3[k1][t][k] = np.asarray(data3[k1][t][k])

k1 = "P1_2but"
#ext = (min(nadphRatio2Test), max(nadphRatio2Test), min(downstreamMaxConc), max(downstreamMaxConc))
##ext=(0,1,0,1)
keyToPlot = ["mdf", "GAPDH", "PFK"] + dsRxnIDs[-1:]#['CRT', 'BTNOLDH']
tToPlot = [0, 1, 3, 8]
#plt.figure()
#
maxCur = {k: -1000. for k in keys}
minCur = {k: 1000. for k in keys}
for k in keys:
	for tPlot in tToPlot:
		maxCur[k] = max(maxCur[k], data3[k1][tPlot][k].max())
		minCur[k] = min(minCur[k], data3[k1][tPlot][k].min())

	if maxCur[k] == minCur[k]:
		maxCur[k] += 1
		minCur[k] -= 1

plt.figure(figsize=(3, 2), dpi=300)
p = 0
for k in keyToPlot:
	#	levels = [minCur[k] + (float(j)/8) * (maxCur[k] - minCur[k]) for j in range(9)]
	for tPlot in tToPlot:
		p += 1
		plt.subplot(len(keyToPlot), 4, p)
		cs = plt.plot(np.ndarray.tolist(data3[k1][tPlot]["maxConc"] * 1000)[0], np.ndarray.tolist(data3[k1][tPlot][k])[0])
#		ext = (min(downstreamMaxConc), max(downstreamMaxConc), min(nadphRatio2Test), max(nadphRatio2Test))
		#ext=(0,1,0,1)

#		cs = plt.contourf(data3[k1][tPlot]["maxConc"] * 1000, data3[k1][tPlot]["nadph"], \
#					   data3[k1][tPlot][k], levels)#, colors=((0.2,0.2,1), (0.2,0.2,.9), (0.2,0.2,.75), (0.2,0.2,.5), (0,0,0), (.5,0.2,0.2), (.75,0.2,0.2), (.9,0.2,0.2), (1,0.2,0.2)))
#		if tPlot == 8:
##			CB = plt.colorbar(cs,fraction=0.046, pad=0.04)
#			CB = plt.colorbar(cs, extend='both', aspect = 5)
#			CB.set_ticks([levels[kL] for kL in [0, 2, 4, 6, 8]])

		ax = plt.gca()
		ax.set_xscale('log')
		if k == "mdf":
			ax.set_ylim(bottom=minCur[k], top=maxCur[k])
		else:
			ax.set_ylim(bottom=0, top=1)

		ax.set_xticks([1e-3, 1e-1, 1e+1])
#		ax.set_yscale('log')
#		ax.yaxis.minorTicks = []
		if tPlot != 0:
			ax.set_yticks([])
		else:
			if k == "mdf":
				ax.set_yticks([minCur[k], 0, maxCur[k]])
			else:
				ax.set_yticks([0, 0.5, 1])
		#elif k == "R5":
			#ax.set_ylabel('NADPH/NADP+ ratio')

		if k != keyToPlot[-1]:
			ax.set_xticks([])
		#elif tPlot == 0:
			#ax.set_xlabel('max. allowable conc. of downstream metabolites')

#		plt.minorticks_off()
		plt.rcParams.update({'font.size': 3.5})

plt.show()
