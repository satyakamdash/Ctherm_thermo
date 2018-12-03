import json
import matplotlib
import numpy as np
import matplotlib.cm as cm
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt

f = open("./examples/res_ctherm_ref_exp_max10_isobut/mdf_isobut.json", "r")
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

#nadphRatio2Test = [0.1, 0.3, 1, 3, 10]#[0.1, 0.2, 0.5, 1, 2, 5, 10]
#downstreamMaxConc = [1e-6, 3e-6, 1e-5, 3e-5, 1e-4, 1e-3, 1e-2] #[1e-6, 2e-6, 5e-6, 1e-5, 2e-5, 5e-5, 1e-4, 2e-4, 5e-4, 1e-3, 2e-3, 5e-3, 1e-2]
#nadphRatio2Test = [1, 10]#[0.1, 0.3, 1, 3, 10]#[0.1, 0.2, 0.5, 1, 2, 5, 10]
#downstreamMaxConc = [1e-4, 1e-2]#[1e-6, 3e-6, 1e-5, 3e-5, 1e-4, 1e-3, 1e-2] #[1e-6, 2e-6, 5e-6, 1e-5, 2e-5, 5e-5, 1e-4, 2e-4, 5e-4, 1e-3, 2e-3, 5e-3, 1e-2]
nadphRatio2Test = sorted([float(k) for k in data["P1_isobut"].keys()])
downstreamMaxConc = sorted([float(k) for k in data["P1_isobut"][str(nadphRatio2Test[0])].keys()])
p, q = len(nadphRatio2Test), len(downstreamMaxConc)

# get the shadow prices of the downstream reactions as a function of isobutanol conc, NADPH/NADP ratios and max concentrations of downstream mets
keys = ["nadph", "maxConc", "ibutConc", "mdf", "atp", "time", "PFK", "GAPDH", "R1", "R2", "R3", "R4", "R5"]
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
					data2[k1]["R1"].append(v5["reaction prices"]["ALS"])
					data2[k1]["R2"].append(v5["reaction prices"][KARI])
					data2[k1]["R3"].append(v5["reaction prices"]["DHAD"])
					data2[k1]["R4"].append(v5["reaction prices"]["KIVD"])
					data2[k1]["R5"].append(v5["reaction prices"]["ADH_IBT"])

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

k1 = "P1_isobut"
#for tPlot in range(9):
#	ext = (min(downstreamMaxConc), max(downstreamMaxConc), min(nadphRatio2Test), max(nadphRatio2Test))
#	#ext=(0,1,0,1)
#	plt.figure()
#	cs = plt.contourf(data3[k1][tPlot]["maxConc"] * 1000, data3[k1][tPlot]["nadph"], data3[k1][tPlot]["mdf"])
#	CB = plt.colorbar(cs, shrink=0.8, extend='both')
#
#	ax = plt.gca()
#	ax.set_xscale('log')
#	ax.set_yscale('log')
#
#	plt.show()

levels = [-1,-0.75, -0.5, -0.25, 0, 0.25, 0.5, 0.75, 1]
maxCur = {k: -1000. for k in keys}
minCur = {k: 1000. for k in keys}
for tPlot in range(9):
	for k in keys:
		maxCur[k] = max(maxCur[k], data3[k1][tPlot][k].max())
		minCur[k] = min(minCur[k], data3[k1][tPlot][k].min())

plt.figure(figsize=(3, 2), dpi=300)
p = 0
for k in ["mdf", "GAPDH", "PFK", "R5"]:
	levels = [minCur[k] + (float(j)/8) * (maxCur[k] - minCur[k]) for j in range(9)]
	for tPlot in [0, 1, 3, 8]:
		p += 1
		plt.subplot(4, 4, p)
		ext = (min(downstreamMaxConc), max(downstreamMaxConc), min(nadphRatio2Test), max(nadphRatio2Test))
		#ext=(0,1,0,1)
		if k == "mdf":
			level = [minCur[k]/8.0 * (8 - j) for j in range(9)] + [maxCur[k]/3.0 * j for j in range(1,4)]
			cs = plt.contourf(data3[k1][tPlot]["maxConc"] * 1000, data3[k1][tPlot]["nadph"], \
					   data3[k1][tPlot][k], level, colors=tuple([(1, 0, 0) for j in range(8)]+[ (0.,0.,0.4 + 0.3*j) for j in range(3)]))
					   #tuple([(1/7.0 * (7-j), 0, 0) for j in range(8)]+[ (0,0,1/3.0*j) for j in range(1,4)]))
		else:
			cs = plt.contourf(data3[k1][tPlot]["maxConc"] * 1000, data3[k1][tPlot]["nadph"], \
					   data3[k1][tPlot][k], levels)#, colors=((0.2,0.2,1), (0.2,0.2,.9), (0.2,0.2,.75), (0.2,0.2,.5), (0,0,0), (.5,0.2,0.2), (.75,0.2,0.2), (.9,0.2,0.2), (1,0.2,0.2)))
		if tPlot == 8:
#			CB = plt.colorbar(cs,fraction=0.046, pad=0.04)
			CB = plt.colorbar(cs, extend='both', aspect = 5)
			if k != "mdf":
				CB.set_ticks([levels[kL] for kL in [0, 2, 4, 6, 8]])

		ax = plt.gca()
		ax.set_xscale('log')
		ax.set_yscale('log')
		ax.yaxis.minorTicks = []
		if tPlot != 0:
			ax.set_yticks([])
		#elif k == "R5":
			#ax.set_ylabel('NADPH/NADP+ ratio')

		if k != "R5":
			ax.set_xticks([])
		#elif tPlot == 0:
			#ax.set_xlabel('max. allowable conc. of downstream metabolites')

		plt.minorticks_off()
		plt.rcParams.update({'font.size': 3.5})

plt.show()

