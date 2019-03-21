import os
import sys
import math
from decimal import Decimal
#usage: python extractStatDrawPerfGraph.py textfile.tex SCIPSDPoutfile1.out testset.solu folder-with-UG-outfiles-1 folder-with-UG-outfiles-2 ...

#example: python extractStatDrawPerfGraph_ug.py /home/gally/gally_Dissertation/Supplement/Tables/UGMISDP.tex /home/gally/gally_Dissertation/resultfiles/Parallel/SCIPSDP/check.SCIPSDPpaper.scipsdp.linux.x86_64.gnu.opt.msk.fb04549.1thread.out /home/gally/gally_SCIPSDP/check/testset/SCIPSDPpaper.solu /home/gally/gally_Dissertation/resultfiles/Parallel/Mosek8.1.0.54/1thread/ /home/gally/gally_Dissertation/resultfiles/Parallel/Mosek8.1.0.54/2thread/ /home/gally/gally_Dissertation/resultfiles/Parallel/Mosek8.1.0.54/4thread/ /home/gally/gally_Dissertation/resultfiles/Parallel/Mosek8.1.0.54/8thread/ /home/gally/gally_Dissertation/resultfiles/Parallel/Mosek8.1.0.54/16thread/ /home/gally/gally_Dissertation/resultfiles/Parallel/Mosek8.1.0.54/32thread/


#example: python extractStatDrawPerfGraph.py /local/gally/results/WarmstartResults20171123/warmstartcomparison_SCIPSDPRIPpaper_wosimplepreopt001.tex /local/gally/results/WarmstartResults20171123/check.SCIPSDPRIPpaper.scipsdp.linux.x86_64.gnu.opt.sdpa.moskito.disable_warmstart.out /local/gally/results/WarmstartResults20171123/check.SCIPSDPRIPpaper.scipsdp.linux.x86_64.gnu.opt.sdpa.moskito.warmstart.out /local/gally/results/WarmstartResults20171123/check.SCIPSDPRIPpaper.scipsdp.linux.x86_64.gnu.opt.sdpa.moskito.warmstartpreoptgap05.out /local/gally/results/WarmstartResults20171123/check.SCIPSDPRIPpaper.scipsdp.linux.x86_64.gnu.opt.sdpa.moskito.warmstartpreoptgap001.out /local/gally/results/WarmstartResults20171123/check.SCIPSDPRIPpaper.scipsdp.linux.x86_64.gnu.opt.sdpa.moskito.warmstartipfactor001_proj2_pdsame.out /local/gally/results/WarmstartResults20171123/check.SCIPSDPRIPpaper.scipsdp.linux.x86_64.gnu.opt.sdpa.moskito.warmstartipfactor05_proj2_pddiff.out /local/gally/results/WarmstartResults20171123/check.SCIPSDPRIPpaper.scipsdp.linux.x86_64.gnu.opt.sdpa.moskito.warmstartipfactor05_proj2_pdsame.out /local/gally/results/WarmstartResults20171123/check.SCIPSDPRIPpaper.scipsdp.linux.x86_64.gnu.opt.sdpa.moskito.warmstart_analcent_05_proj2.out /local/gally/results/WarmstartResults20171123/check.SCIPSDPRIPpaper.scipsdp.linux.x86_64.gnu.opt.sdpa.moskito.warmstartprojminev-1pdsame.out /local/gally/results/WarmstartResults20171123/check.SCIPSDPRIPpaper.scipsdp.linux.x86_64.gnu.opt.sdpa.moskito.warmstartipfactor05_proj4.out /local/gally/results/WarmstartResults20171123/check.SCIPSDPRIPpaper.scipsdp.linux.x86_64.gnu.opt.sdpa.moskito.warmstartroundinfonly.out

#example: python extractStatDrawPerfGraph.py /local/gally/results/WarmstartResults20171123/warmstartcomparison_SCIPSDPRIPpaper.tex /local/gally/results/WarmstartResults20171123/check.SCIPSDPRIPpaper.scipsdp.linux.x86_64.gnu.opt.sdpa.moskito.disable_warmstart.out /local/gally/results/WarmstartResults20171123/check.SCIPSDPRIPpaper.scipsdp.linux.x86_64.gnu.opt.sdpa.moskito.warmstart.out /local/gally/results/WarmstartResults20171123/check.SCIPSDPRIPpaper.scipsdp.linux.x86_64.gnu.opt.sdpa.moskito.warmstartpreoptgap05.out /local/gally/results/WarmstartResults20171123/check.SCIPSDPRIPpaper.scipsdp.linux.x86_64.gnu.opt.sdpa.moskito.warmstartpreoptgap01.out /local/gally/results/WarmstartResults20171123/check.SCIPSDPRIPpaper.scipsdp.linux.x86_64.gnu.opt.sdpa.moskito.warmstartpreoptgap001.out /local/gally/results/WarmstartResults20171123/check.SCIPSDPRIPpaper.scipsdp.linux.x86_64.gnu.opt.sdpa.moskito.warmstartpreoptgap0001.out /local/gally/results/WarmstartResults20171123/check.SCIPSDPRIPpaper.scipsdp.linux.x86_64.gnu.opt.sdpa.moskito.warmstartipfactor001_proj2_pddiff.out /local/gally/results/WarmstartResults20171123/check.SCIPSDPRIPpaper.scipsdp.linux.x86_64.gnu.opt.sdpa.moskito.warmstartipfactor001_proj2_pdsame.out /local/gally/results/WarmstartResults20171123/check.SCIPSDPRIPpaper.scipsdp.linux.x86_64.gnu.opt.sdpa.moskito.warmstartipfactor05_proj2_pddiff.out /local/gally/results/WarmstartResults20171123/check.SCIPSDPRIPpaper.scipsdp.linux.x86_64.gnu.opt.sdpa.moskito.warmstartipfactor05_proj2_pdsame.out /local/gally/results/WarmstartResults20171123/check.SCIPSDPRIPpaper.scipsdp.linux.x86_64.gnu.opt.sdpa.moskito.warmstart_analcent_05_proj2.out /local/gally/results/WarmstartResults20171123/check.SCIPSDPRIPpaper.scipsdp.linux.x86_64.gnu.opt.sdpa.moskito.warmstart_analcent_001_proj2.out /local/gally/results/WarmstartResults20171123/check.SCIPSDPRIPpaper.scipsdp.linux.x86_64.gnu.opt.sdpa.moskito.warmstartprojminev10.out /local/gally/results/WarmstartResults20171123/check.SCIPSDPRIPpaper.scipsdp.linux.x86_64.gnu.opt.sdpa.moskito.warmstartprojminev-1pddifferent.out /local/gally/results/WarmstartResults20171123/check.SCIPSDPRIPpaper.scipsdp.linux.x86_64.gnu.opt.sdpa.moskito.warmstartprojminev-1pdsame.out /local/gally/results/WarmstartResults20171123/check.SCIPSDPRIPpaper.scipsdp.linux.x86_64.gnu.opt.sdpa.moskito.warmstartipfactor001_proj4.out /local/gally/results/WarmstartResults20171123/check.SCIPSDPRIPpaper.scipsdp.linux.x86_64.gnu.opt.sdpa.moskito.warmstartipfactor05_proj4.out /local/gally/results/WarmstartResults20171123/check.SCIPSDPRIPpaper.scipsdp.linux.x86_64.gnu.opt.sdpa.moskito.warmstartroundinfonly.out

#TODO: irgendwas klappt mit unsolved nicht, Gesamtdurchschnitt ist viel hoeher als pro Testset, gefuehlt cardleastsquares zu hoch, rest zu gering
maxinstances = 321 # maximum number of instances
completeTable = 1 # make table listing all instances
overviewTable = 0 # make table for each setting listing all subsets
overviewBySubset = 0 # make table for each subset listing all settings
overviewTotal = 0 # make table for complete testset listing all settings
performancePlotOverview = 0 # draw a performance plot for the complete testset
performancePlotBySubset = 0 # draw a performance plot for each subset
racingWinnerBySubset = 0 # make table listing the winning parameters for each instanceset
racingWinnerOverall = 0 # make table listing the winning parameters for each instanceset
racingWinnerBarplotAll = 0 #make racing winner barplot for all settings
subsets = [[0,64],[65,133],[134,193]] # combine all instances from subsets[i][0] to subsets[i][1] to a single row and make a performance plot for each subset
#subsets = [[0,64],[65,133],[134,193], [194,319]] # combine all instances from subsets[i][0] to subsets[i][1] to a single row and make a performance plot for each subset
#subsets = [[0,125]]
#subsets = [[0,1],[2,3]]
#colors = ["blue", "red"]
colors = ["blue", "red", "green", "Apricot", "gray", "SkyBlue","yellow", "DarkOrchid", "Aquamarine", "Sepia", "black", "LimeGreen", "CadetBlue", "Lavender", "Tan", "SeaGreen", "Bittersweet", "BlueViolet", "Cerulean", "Brown", "Cyan", "ForestGreen", "Goldenrod",  "JungleGreen", "Mahogany", "Melon", "Mulberry", "OliveGreen", "OrangeRed", "Peach", "PineGreen", "ProcessBlue", "RawSienna", "RedOrange", "Rhodamine", "RoyalPurple", "Salmon", "SpringGreen", "TealBlue", "Turquoise", "VioletRed", "WildStrawberry", "YellowGreen", "BlueGreen", "BrickRed", "BurntOrange", "CarnationPink", "CornflowerBlue", "Dandelion", "Emerald", "Fuchsia", "GreenYellow", "Magenta", "Maroon", "MidnightBlue", "NavyBlue", "Orange", "Orchid", "Periwinkle", "Plum", "Purple", "RedViolet", "RoyalBlue", "RubineRed", "Thistle", "Violet", "YellowOrange"]  # colors to use for the performance plot, needs to have at least as many elements as settings were used
#colors = ["red", "red", "red", "red", "red", "red", "red", "red", "red", "red", "red", "red", "red", "red", "red", "red", "red", "red", "red", "red", "red", "red", "red", "red", "red", "red", "red", "red", "red", "red", "red", "red", "red", "red", "red", "red", "red", "red", "red", "red", "red", "red", "red", "red", "red", "red", "red", "red", "red", "red", "red", "red"]
#subsetnames = ["first half", "second half"]
subsetnames = ["cardinality constrained least squares", "min-$k$-partitioning", "truss topology", "RIP"] # names given to the subsets used in the tables, needs to have at least as many elements as there are subsets
timeshift = 10
nodeshift = 100
itershift = 1000
fixingsshift = 100
heurarithmetic = 0 #TODO: also do this for percentages or at least count number of solved instances
fixingsarithmetic = 0
itergeomonlysolved = 1
heurshift = 1
nodegeomonlysolved = 1 # in the tables for each setting only include solved instances into the geometric mean of the nodes, for the overviews only include instances which were solved by all settings
overviewlongtable = 0 # should overview tables be done as longtables
plotsbysettingssubsets = 1 # make seperate plots and tables for each subset in settingssubsets
performanceprofileNormed = 0 #should the performance plots be normalized by the fastest solver like in Dolan & More to get the relative time or should we just plot over the absolute time instead?
#settingsets = [[0,1,2,3,4,5], [6,7,8,9,10,11], [12,13,14,15,16,17], [18,19,20,21,22,23], [24,25,26,27,28,29], [30,31,32,33,34,35], [36,37,38,39,40,41], [42,43,44,45,46,47], [0,6,12,18,24,30,36,42], [1,7,13,19,25,31,37,43], [2,8,14,20,26,32,38,44], [3,9,15,21,27,33,39,45], [4,10,16,22,28,34,40,46], [5,11,17,23,29,35,41,47]]
#settingsets = [[0,2,3,4,5,11,17,23], [17,18,19,20,21,22], [23,24,25,26,27,28]]
#settingsets = [[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17]]
settingsets = [[0,1,2,3,4,5,6]]
#settingsets = [[0,1,2,3,4,5,6,7,8,9,10]]
#settingsets = [[0,1,2,3,4,5,6,7,8,9]]
#settingsets = [[0,1,2,3]]
usesettingnames = 1
#settingnames = ["SDPA", "DSDP", "MOSEK"]
#settingnames = [ "no warmstart", "simple warmstart", "preoptgap 0.5", "preoptgap 0.1", "preoptgap 0.01", "preoptgap 0.001", "0.01 id pddiff", "0.01 id pdsame", "0.5 id pddiff", "0.5 id pdsame", "0.5 analcent", "0.01 analcent", "proj minev 10", "proj minev auto pddiff", "proj minev auto pdsame", "roundingprob 0.01 id", "roundingprob 0.5 id", "roundingprob inf only"]
#settingnames = [ "no warmstart", "simple warmstart", "preoptgap 0.01", "preoptgap 0.5", "0.01 id pdsame", "0.5 id pddiff", "0.5 id pdsame", "0.5 analcent", "proj minev auto pdsame", "roundingprob 0.5 id", "roundingprob inf only"]
settingnames = ["SCIP-SDP", "UG-MISDP-1thread", "UG-MISDP-2threads", "UG-MISDP-4threads", "UG-MISDP-8threads", "UG-MISDP-16threads", "UG-MISDP-32threads"]
epsilon = 0.001
rel_epsilon = 0.05
captions=["","Complete results for \ugmisdp with 1 thread on shared memory environment of Intel Xeon E5-4650 CPUs running at \SI{2.70}{GHz} with \SI{512}{GB} of shared RAM","Complete results for \ugmisdp with 2 threads on shared memory environment of Intel Xeon E5-4650 CPUs running at \SI{2.70}{GHz} with \SI{512}{GB} of shared RAM","Complete results for \ugmisdp with 4 threads on shared memory environment of Intel Xeon E5-4650 CPUs running at \SI{2.70}{GHz} with \SI{512}{GB} of shared RAM","Complete results for \ugmisdp with 8 threads on shared memory environment of Intel Xeon E5-4650 CPUs running at \SI{2.70}{GHz} with \SI{512}{GB} of shared RAM","Complete results for \ugmisdp with 16 threads on shared memory environment of Intel Xeon E5-4650 CPUs running at \SI{2.70}{GHz} with \SI{512}{GB} of shared RAM","Complete results for \ugmisdp with 32 threads on shared memory environment of Intel Xeon E5-4650 CPUs running at \SI{2.70}{GHz} with \SI{512}{GB} of shared RAM"]
shortcaptions=["","Complete results for \ugmisdp with 1 thread","Complete results for \ugmisdp with 2 threads","Complete results for \ugmisdp with 4 threads","Complete results for \ugmisdp with 8 threads","Complete results for \ugmisdp with 16 threads","Complete results for \ugmisdp with 32 threads"]
labels=["","UGMISDP1thread","UGMISDP2threads","UGMISDP4threads","UGMISDP8threads","UGMISDP16threads","UGMISDP32threads"]

def readFile(arg,i):
	j=0;
	file=open(arg, "r")
	filestring = file.read()
	while filestring.find("@01") > -1:
		substring = filestring[filestring.find("@01"):filestring.find("@04")]
		assert(substring.find("@01 ") > -1)
		names[i][j] = os.path.basename(substring[substring.find("@01 ") + 4:].split()[0])
		if substring.find("Solving Time (sec) : ") == -1:
			times[i][j] = 3600
		else:
			time = substring[substring.find("Solving Time (sec) : ") + 21:].split()[0]
			times[i][j] = float(time)
			if times[i][j] > 3600.0:
				times[i][j] = 3600
		if substring.find("Solving Nodes      : ") == -1:
			nodes[i][j] = "-"
		else:
			node = substring[substring.find("Solving Nodes      : ") + 21:].split()[0]
			nodes[i][j] = int(node)
		if substring.find("sdpredcost       :") == -1:
			sdpredcostfixings[i][j] = "-"
		else:
			sdpredcostfixing = substring[substring.find("sdpredcost       :") + 18:].split()[3]
			sdpredcostfixings[i][j] = int(sdpredcostfixing)
		if substring.find("sdpfracdiving    :") == -1:
			fracdivefound[i][j] = "-"
		else:
			sdpfracdive = substring[substring.find("sdpfracdiving    :") + 18:].split()[3]
			fracdivefound[i][j] = int(sdpfracdive)
		if substring.find("Dual Bound         : ") == -1:
			dualresults[i][j] = -1e20
		else:
			dualresult = substring[substring.find("Dual Bound         : ") + 21:].split()[0]
			dualresults[i][j] = float(dualresult)
		if substring.find("Primal Bound       : ") == -1:
			primalresults[i][j] = 1e20
		else:
			primalresult = substring[substring.find("Primal Bound       : ") + 21:].split()[0]
			primalresults[i][j] = float(primalresult)
		if substring.find("Gap                : ") == -1:
			gaps[i][j] = "1e20"
		else:
			gap = substring[substring.find("Gap                : ") + 21:].split()[0]
			gaps[i][j] = gap
		if substring.find("SDP iterations:				") == -1:
			sdpiters[i][j] = "-"
		else:
			sdpiter = substring[substring.find("SDP iterations:") + 15:].split()[0]
			sdpiters[i][j] = int(sdpiter)
		if substring.find("Percentage penalty formulation used:	 ") == -1:
			penalties[i][j] = "-"
		else:
			penalty = substring[substring.find("Percentage penalty formulation used:") + 38:].split()[0]
			penalties[i][j] = float(penalty)
		if substring.find("Percentage unsolved even with penalty:	") == -1:
			unsolved[i][j] = "-"
		else:
			unsolv = substring[substring.find("Percentage unsolved even with penalty:") + 40:].split()[0]
			unsolved[i][j] = float(unsolv)
		filestring = filestring[filestring.find("@04")+36:]
		j = j+1
		assert( j < maxinstances)
	ninstances[i] = j
	file.close()

def readUG(i,j,folder,soluline):
	outfile = folder + soluline.split()[1] + "_LC0_T.status"
	file=open(outfile, "r")
	filestring = file.readlines()
	#find racing winner line
	l=0
	while True:
		if "Racing ramp-up finished" in filestring[l]:
			racingwinner[i][j] = (filestring[l].split()[len(filestring[l].split()) - 1])[:-1]
			break
		elif l == len(filestring) - 1:
			racingwinner[i][j] = -1
			break
		else:
			l=l+1
	#find solving time
	l=0
	while True:
		if "Total Time         :" in filestring[l]:
			times[i][j] = float(filestring[l].split()[3])
			break
		elif l == len(filestring) - 1:
			times[i][j] = 3600.0
			break
		else:
			l=l+1
	#find nodes
	while True:
		if "nodes (total)    :" in filestring[l]:
			nodes[i][j] = float(filestring[l].split()[3])
			break
		elif l == len(filestring) - 1:
			nodes[i][j] = "-"
			break
		else:
			l=l+1
	#find primal bound
	while True:
		if "Primal Bound     :" in filestring[l]:
			if filestring[l].split()[3] == "-":
				primalresults[i][j] = 1e20
			else:
				primalresults[i][j] = float(filestring[l].split()[3])
			break
		elif l == len(filestring) - 1:
			primalresults[i][j] = 1e20
			break
		else:
			l=l+1
	#find dual bound
	while True:
		if "Dual Bound       :" in filestring[l]:
			if filestring[l].split()[3] == "-":
				dualresults[i][j] = -1e20
			else:
				dualresults[i][j] = float(filestring[l].split()[3])
			break
		elif l == len(filestring) - 1:
			dualresults[i][j] = -1e20
			break
		else:
			l=l+1
	#find gap
	while True:
		if "Gap                :" in filestring[l]:
			gaps[i][j] = filestring[l].split()[2]
			break
		elif l == len(filestring) - 1:
			gaps[i][j] = "infinite"
			break
		else:
			l=l+1

	#compare with solu file
	if "=opt=" in soluline.split()[0]:
		if (dualresults[i][j] - float(soluline.split()[2])) / max([1, abs(float(soluline.split()[2]))]) > rel_epsilon:
			print("wrong dual bound for " + settingnames[i] + " on instance " + soluline.split()[1] + ": " + str(dualresults[i][j]) + " > " + soluline.split()[2] + " = opt result\n")
                        times[i][j] = 3600.0 ##### TODO change for full table
                        gaps[i][j] = "infinite" ##### TODO change for full table
		if (primalresults[i][j] - float(soluline.split()[2])) / max([1, abs(float(soluline.split()[2]))]) < -rel_epsilon:
			print("wrong primal bound for " + settingnames[i] + " on instance " + soluline.split()[1] + ": " + str(primalresults[i][j]) + " < " + soluline.split()[2] + " = opt result\n")
                        times[i][j] = 3600.0 ##### TODO change for full table
                        gaps[i][j] = "infinite" ##### TODO change for full table
	if "=best=" in soluline.split()[0]:
		if (dualresults[i][j] - float(soluline.split()[2])) / max([1, abs(float(soluline.split()[2]))]) > rel_epsilon:
			print("wrong dual bound  for " + settingnames[i] + " on instance " + soluline.split()[1] + ": " + str(dualresults[i][j]) + " > " + soluline.split()[2] + " = best known solution\n")
                        times[i][j] = 3600.0 ##### TODO change for full table
                        gaps[i][j] = "infinite" ##### TODO change for full table
	file.close

def instanceSolved(i, j):
	""" check if the duality gap is zero """
	if gaps[i][j] == "infinite" or ("Large" in gaps[i][j])  or "--" in gaps[i][j] or float(gaps[i][j]):
		return False
	return True

def allSolved(ind, settings):
	""" check if for this index all settings could solve the instance"""
	for i in settings:
		if not instanceSolved(i, ind):
			return False
	return True

def makeCompleteTable(arg, file, i):
	file.write("\\begin{longtable}[ht] {p{.3\\textwidth} | p{.09\\textwidth} p{.09\\textwidth} p{.056\\textwidth} p{.05\\textwidth} p{.05\\textwidth} p{.05\\textwidth} p{.05\\textwidth} p{.05\\textwidth} p{.05\\textwidth} p{.05\\textwidth}} \\toprule \n problem & dual bound &  primal bound & gap & nodes & fixings & fracdive & time & sdpiter & penalty & unsolved \% \\\ \midrule \n")
	for j in range(ninstances[i]):
		file.write("$" + names[i][j] + "$ & $ %.4g" % dualresults[i][j] + "$ & $ %.4g" % primalresults[i][j] + "$ & $" + str(gaps[i][j]) + "$ & $" + str(nodes[i][j]) + "$ & $" + str(sdpredcostfixings[i][j])+ "$ & $" + str(fracdivefound[i][j])  + "$ & $" + str(times[i][j]) + "$ & $" + str(sdpiters[i][j]) + "$ & $" + str(penalties[i][j]) + "$ & $" + str(unsolved[i][j]) + "$ \\\ \n")
	file.write("\\bottomrule \caption{Results for $" + str((os.path.basename(arg)).split(".")[7]) + " | " + str((os.path.basename(arg)).split(".")[9]) + "$} \end{longtable} \n")

def makeCompleteTableCaption(arg, i):
	file.write("\\newpage \n \\begin{scriptsize} \n  \\setlength{\\tabcolsep}{2pt} \n \\tablehead{\\toprule \n")
	file.write("problem & dbound &  pbound & gap & nodes & time & racingwinner ")
	file.write("\\\ \\midrule} \n \\tabletail{ \\midrule \\multicolumn{3}{@{}l}{continued on next page \\dots}\\\ \\bottomrule} \n \\tablelasttail{\\bottomrule} \n \\tablecaption[")
	file.write(shortcaptions[i])
	file.write("]{")
	file.write(captions[i])
	file.write("}\label{")
	file.write(labels[i])
	file.write("}\n")
	file.write("\\begin{xtabular*}{\\textwidth}{@{\extracolsep{\\fill}}lrrrrrr")
	file.write("@{}} \n ")
	for j in range(ninstances[0]):
		file.write(names[0][j].replace("_", "\\_").split(".")[0] + "& ")
		if dualresults[i][j] > -1e20:
			file.write("\\num{" + "%.2f" % float(dualresults[i][j]) + "} ")
		else:
			file.write(" $-\infty$ ")
		file.write("& ")
		if primalresults[i][j] < 1e20:
			file.write("\\num{" + "%.2f" % float(primalresults[i][j]) + "}")
		else:
			file.write(" $\infty$ ")
		file.write("& ")
		if gaps[i][j] != "infinite" and gaps[i][j] != "1e20":
			if gaps[i][j] != "-":
				file.write("\\SI{" + "%.2f" % float(gaps[i][j]) + "}{\percent} & ")
			else:
				file.write("-- & ")
		else:
			file.write("$\infty$ & ")
		if nodes[i][j] == "-":
			file.write("-- & ")
		else:
			file.write("\\num{" + "%.0f" % float(nodes[i][j]) + "} &")
		if times[i][j] == "-":
			file.write("-- & ")
		else:
			file.write("\\num{" + "%.1f" % float(times[i][j]) + "} &")
		if racingwinner[i][j] == -1:
			file.write("-- ")
		else:
			file.write("\\num{" + "%.0f" % float(racingwinner[i][j]) + "}")
		file.write("\\\ \n")
	file.write("  \\end{xtabular*} \n \\end{scriptsize} \n")

def makeOverviewTable(arg, file, i):
	if overviewlongtable:
		file.write("\\begin{longtable}[ht] {p{.3\\textwidth} | p{.05\\textwidth} p{.09\\textwidth} p{.09\\textwidth} p{.05\\textwidth} p{.05\\textwidth} p{.09\\textwidth} p{.05\\textwidth} p{.05\\textwidth}} \\toprule \n type & \# solved & nodes & fixings & fracdive & time & sdpiter & penalty & unsolved \% \\\ \midrule \n")
	else:
		file.write("\\begin{table}[ht]  \n \\caption{Results for $" + str((os.path.basename(arg)).split(".")[7]) + " | " + str((os.path.basename(arg)).split(".")[9]) + "$} \n \\begin{tabular*}{\\textwidth}{l | r r r r r r} \\toprule \n type & \# solved & nodes & time & sdpiter & penalty & unsolved \% \\\ \midrule \n")
	subsetnum = 0
	for (s,f) in subsets:
		assert(s >= 0)
		assert(f < ninstances[i])
		nsolved = 0
		nodegeom = 1.0
		fixingsgeom = 1.0
		if heurarithmetic:
			totalheur = 0.0
			nsolved = 0
		else:
			totalheur = 1.0
		timegeom = 1.0
		itergeom = 1.0
		totalpenalty = 0
		totalunsolved = 0
		j = s
		numinstances = f - s + 1
		while j <= f:
			if instanceSolved(i,j):
				nsolved = nsolved + 1
			if nodes[i][j] != "-" and (not nodegeomonlysolved or instanceSolved(i,j)):
				nodegeom = math.pow(nodegeom, float(j-s) / float(j-s+1)) * math.pow(float(nodes[i][j] + nodeshift), 1.0/float(j-s+1))
			if sdpredcostfixings[i][j] != "-" and instanceSolved(i,j):
				fixingsgeom = math.pow(fixingsgeom, float(j-s) / float(j-s+1)) * math.pow(float(sdpredcostfixings[i][j] + fixingsshift), 1.0/float(j-s+1))
			if heurarithmetic and fracdivefound[i][j] != "-":
				totalheur = totalheur + float(fracdivefound[i][j])
				nsolved = nsolved + 1
			elif fracdivefound[i][j] != "-":
				totalheur = math.pow(totalheur, float(j-s) / float(j-s+1)) * math.pow(float(fracdivefound[i][j] + heurshift), 1.0/float(j-s+1))
			timegeom = math.pow(timegeom, float(j-s) / float(j-s+1)) * math.pow(float(times[i][j] + timeshift), 1.0/float(j-s+1))
			if sdpiters[i][j] != "-":
				itergeom = math.pow(itergeom, float(j-s) / float(j-s+1)) * math.pow(float(sdpiters[i][j] + itershift), 1.0/float(j-s+1))
			if penalties[i][j] != "-":
				totalpenalty = totalpenalty + penalties[i][j]
			if unsolved[i][j] != "-":
				totalunsolved = totalunsolved + unsolved[i][j]
			j = j + 1
		nodegeom = nodegeom - nodeshift
		fixingsgeom = fixingsgeom - fixingsshift
		if fixingsgeom < epsilon:
			fixingsgeom = 0.0
		if heurarithmetic:
			avgheur = totalheur / nsolved
		else:
			avgheur = totalheur - heurshift
		timegeom = timegeom - timeshift
		itergeom = itergeom - itershift
		avgpenalty = float(totalpenalty) / float(numinstances)
		avgunsolved = float(totalunsolved) / float(numinstances)
		file.write(subsetnames[subsetnum] + "& $" + str(nsolved) + "$ & $ %.4g" % nodegeom + "$ & $ %.4g" % fixingsgeom + "$ & $ %.4g" % avgheur + "$ & $ %.4g" % timegeom + "$ & $ %.4g" % itergeom + "$ & $ %.4g" % avgpenalty + "$ & $ %.4g" % avgunsolved + "$ \\\ \n")
		subsetnum = subsetnum + 1
	#make one more line for the whole testset
	file.write("\midrule \n")
	nsolved = 0
	nodegeom = 1.0
	fixingsgeom = 1.0
	if heurarithmetic:
		totalheur = 0.0
		nsolved = 0
	else:
		totalheur = 1.0
	timegeom = 1.0
	itergeom = 1.0
	totalpenalty = 0
	totalunsolved = 0
	s = 0
	f = ninstances[i] - 1
	j = s
	numinstances = f - s + 1
	while j <= f:
		if instanceSolved(i,j) :
			nsolved = nsolved + 1
		if nodes[i][j] != "-" and (not nodegeomonlysolved or instanceSolved(i,j)):
			nodegeom = math.pow(nodegeom, float(j-s) / float(j-s+1)) * math.pow(float(nodes[i][j] + nodeshift), 1.0/float(j-s+1))
		if sdpredcostfixings[i][j] != "-" and instanceSolved(i,j):
			fixingsgeom = math.pow(fixingsgeom, float(j-s) / float(j-s+1)) * math.pow(float(sdpredcostfixings[i][j] + fixingsshift), 1.0/float(j-s+1))
		if heurarithmetic and fracdivefound[i][j] != "-":
			totalheur = totalheur + float(fracdivefound[i][j])
			nsolved = nsolved + 1
		elif fracdivefound[i][j] != "-":
			totalheur = math.pow(totalheur, float(j-s) / float(j-s+1)) * math.pow(float(fracdivefound[i][j] + heurshift), 1.0/float(j-s+1))
		timegeom = math.pow(timegeom, float(j-s) / float(j-s+1)) * math.pow(float(times[i][j] + timeshift), 1.0/float(j-s+1))
		if sdpiters[i][j] != "-":
			itergeom = math.pow(itergeom, float(j-s) / float(j-s+1)) * math.pow(float(sdpiters[i][j] + itershift), 1.0/float(j-s+1))
		if penalties[i][j] != "-":
			totalpenalty = totalpenalty + penalties[i][j]
		if unsolved[i][j] != "-":
			totalunsolved = totalunsolved + unsolved[i][j]
		j = j + 1
	nodegeom = nodegeom - nodeshift
	fixingsgeom = fixingsgeom - fixingsshift
	if fixingsgeom < epsilon:
		fixingsgeom = 0.0
	if heurarithmetic:
		avgheur = totalheur / nsolved
	else:
		avgheur = totalheur - heurshift
	timegeom = timegeom - timeshift
	itergeom = itergeom - itershift
	avgpenalty = float(totalpenalty) / float(numinstances)
	avgunsolved = float(totalunsolved) / float(numinstances)
	file.write("Total & $" + str(nsolved) + "$ & $ %.4g" % nodegeom + "$ & $ %.4g" % fixingsgeom + "$ & $ %.4g" % avgheur + "$ & $ %.4g" % timegeom + "$ & $ %.4g" % itergeom + "$ & $ %.4g" % avgpenalty + "$ & $ %.4g" % avgunsolved + "$ \\\ \n")
	if overviewlongtable:
		file.write("\\bottomrule \n \\caption{Results for $" + str((os.path.basename(arg)).split(".")[7]) + " | " + str((os.path.basename(arg)).split(".")[9]) + "$} \n \\end{longtable} \n")
	else:
		file.write("\\bottomrule \n \\end{tabular*} \n \\end{table} \n")

def makeOverviewBySubset(file, args, settings, s, f, caption, label, settingnames, fixingscol):
	if overviewlongtable:
		file.write("\\begin{longtable}[ht] {p{.5\\textwidth} | p{.15\\textwidth} p{.15\\textwidth} p{.15\\textwidth}} \\toprule \n settings & solved & nodes & time  \\\ \midrule \n")
	else:
		file.write("\\begin{table} \n \\begin{scriptsize} \\caption{" + caption + "} \n \\label{" + label + "} \n \\begin{tabular*}{\\textwidth}{@{}l@{\\;\\;\extracolsep{\\fill}}rrr")
		if fixingscol:
			file.write("r")
		file.write("@{}}\\toprule \n")
		file.write(" settings &  solved & nodes & time\\\ \midrule \n")
	i = 0
	ind = 0
	argind = -1
	for arg in args:
		if i in settings:
			nsolved = 0
			nodegeom = 1.0
			timegeom = 1.0
			nunaborted = 0
			nnodescounted = 0
			j = s
			numinstances = f - s + 1
			while j <= f:
				if nodes[i][j] != "-":
					nunaborted = nunaborted + 1
				if nodes[i][j] != "-" and (not nodegeomonlysolved or allSolved(j, settings)):
					nnodescounted = nnodescounted + 1
				j = j + 1

			j = s
			while j <= f:
				if instanceSolved(i,j) :
					nsolved = nsolved + 1
				if nodes[i][j] != "-" and (not nodegeomonlysolved or allSolved(j, settings)):
					nodegeom = nodegeom * math.pow(float(nodes[i][j] + nodeshift), 1.0/float(nnodescounted))
				timegeom = timegeom * math.pow(float(times[i][j] + timeshift), 1.0/float(f-s+1))
				j = j + 1

			nodegeom = nodegeom - nodeshift
			timegeom = timegeom - timeshift
			file.write(settingnames[ind] + " & \\num{" + "%.0f" % nsolved + "} & \\num{" + "%.2f" % nodegeom + "} & \\num{" + "%.2f" % timegeom + "} \\\ \n")
			ind = ind + 1
		i = i + 1
	if overviewlongtable:
		file.write("\\bottomrule \n \\caption{" + caption + "} \n \\label{" + label + "} \n \\end{longtable} \n")
	else:
		file.write("\\bottomrule \n \\end{tabular*} \n \end{scriptsize} \n \\end{table} \n")

def makeRacingWinnerBySubset(file, args, settings, s, f, caption, label, settingnames):
	file.write("\\begin{table} \n \\begin{scriptsize} \\caption{" + caption + "} \n \\label{" + label + "} \n \\begin{tabular*}{\\textwidth}{@{}l@{\\;\\;\extracolsep{\\fill}}rrrrrrrrrrrrrrrrrrrrrrrrrrrrrrrr")
	file.write("@{}}\\toprule \n")
	file.write(" settings & 1 & 2 & 3 & 4 & 5 & 6 & 7 & 8 & 9 & 10 & 11 & 12 & 13 & 14 & 15 & 16 & 17 & 18 & 19 & 20 & 21 & 22 & 23 & 24 & 25 & 26 & 27 & 28 & 29 & 30 & 31 & 32\\\ \midrule \n")
	i = 1 #ignore SCIP-SDP run
	ind = 1
	argind = -1
	for arg in args:
		if i in settings:
			chosenparam = [0 for x in range(32)]
			j = s
			numinstances = f - s + 1
			while j <= f:
				if racingwinner[i][j] != -1 :
					assert(int(racingwinner[i][j]) > 0)
					assert(int(racingwinner[i][j]) <= 32)
					chosenparam[int(racingwinner[i][j]) - 1] += 1
				j = j + 1
			file.write(settingnames[ind])
			for pos in range(32):
				file.write(" & " + str(chosenparam[pos]))
			file.write(" \\\ \n")
			ind = ind + 1
		i = i + 1
	file.write("\\bottomrule \n \\end{tabular*} \n \end{scriptsize} \n \\end{table} \n")

def makeRacingWinnerPlot(file, i, caption, label, settingnames):
        chosenparam = [[0 for x in range(32)] for y in range(len(subsets))]
        ind = 0
        for (s,f) in subsets:
		j = s
		while j <= f:
		        if racingwinner[i][j] != -1 :
			        assert(int(racingwinner[i][j]) > 0)
			        assert(int(racingwinner[i][j]) <= 32)
			        chosenparam[ind][int(racingwinner[i][j]) - 1] += 1
		        j = j + 1
		ind = ind + 1
        #draw bar plot
        file.write("\\begin{tikzpicture} \n \\begin{axis}[scale only axis, width=0.67*\\textwidth, ybar stacked, ymin=0, xmin=-1, xmax=33, legend style={at={(0.5,-0.20)}, anchor=north,legend columns=-1}, ylabel={\#racing winner}, symbolic x coords={0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33}, xtick=data, xticklabels={1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32}, xmin=0, xmax=33]\n")
        for ind in range(len(subsets)):
                file.write("\\addplot+[ybar] plot coordinates { ")
                for pos in range(32):
                        file.write("(" + str(pos+1) + "," + str(chosenparam[ind][pos]) + ") ")
                file.write("};\n")
        file.write("\\legend{")
        for ind in range(len(subsets)):
                file.write(subsetnames[ind])
                if ind < range(len(subsets)):
                        file.write(", ")
        file.write("}\n \\end{axis} \n \\end{tikzpicture}\n\n")

def makePerformancePlot(file, args, settings, s, f, ind):
	performance=[[0 for x in range(f-s+1)] for x in range(len(sys.argv))] #initialize performance matrix
	j = s
	assert(len(colors) >= len(args))
	while j <= f:
		best = float('inf')
		# compute optimum
		for arg in settings:
			if times[arg][j] < best:
				best = times[arg][j]
		# divide all values by the optimum
		for arg in settings:
			if instanceSolved(arg,j):
				if performanceprofileNormed:
					performance[arg][j-s] = times[arg][j] / best
				else:
					performance[arg][j-s] = times[arg][j]
			else:
				performance[arg][j-s] = 1e20
		j = j + 1
	# sort performance for each arg
	for arg in settings:
		(performance[arg]).sort()
	#write the plot
	i = 0
	color = 0
	file.write("\\begin{figure} \n \\centering \n \\begin{tikzpicture} \n \\begin{semilogxaxis}[xlabel={time},ylabel={\\# solved instances}, legend style={cells={anchor=east},legend pos=outer north east,}]")
	for arg in args:
		if i in settings:
			file.write("\\addplot[" + colors[color] + "]table[row sep=crcr]{x 	y \\\ \n")
			for j in range(f - s + 1):
				# if multiple instances are solved at the same (relative) time, only write the last one
				if (j == f -s or performance[i][j] < performance[i][j+1]) and performance[i][j] < 1e20:
					file.write(str(performance[i][j]) + " " + str(j + 1) + " \\\ \n")
			if usesettingnames and settingnames[i] != "":
				file.write("}; \n \\addlegendentry{" + settingnames[i] + "}")
			else:
				file.write("}; \n \\addlegendentry{$" + str((os.path.basename(arg)).split(".")[7]) + " | " + str((os.path.basename(arg)).split(".")[9]) + "$}")
			color = color + 1
		i = i + 1
	file.write("\\end{semilogxaxis} \n \\end{tikzpicture} \n")
	if ind > -1:
		file.write("\\caption{Performanceplot for " + subsetnames[ind] + "} \n")
	else:
		file.write("\\caption{Performanceplot for full testset} \n")
	file.write("\\end{figure} \n")


if __name__=="__main__":
	"""give any number of .out-files for the same testset, then loops over them and returns a .tex-file given as first argument with some tables and a performance graph """
	ninstances=[0 for x in range(len(sys.argv) -2)] #initialize ninstances matrix
	names=[[0 for x in range(maxinstances)] for x in range(len(sys.argv) -3)] #initialize names matrix
	dualresults=[[0 for x in range(maxinstances)] for x in range(len(sys.argv) -3)] #initialize dual results matrix
	primalresults=[[0 for x in range(maxinstances)] for x in range(len(sys.argv) -3)] #initialize primal results matrix
	gaps=[[0 for x in range(maxinstances)] for x in range(len(sys.argv) -3)]   #initialize gaps matrix
	times=[[0 for x in range(maxinstances)] for x in range(len(sys.argv) -3)]   #initialize solvingtime matrix
	nodes=[[0 for x in range(maxinstances)] for x in range(len(sys.argv) -3)]   #initialize nodes matrix
	sdpredcostfixings=[[0 for x in range(maxinstances)] for x in range(len(sys.argv) -3)]   #initialize redcostfixings matrix
	fracdivefound=[[0 for x in range(maxinstances)] for x in range(len(sys.argv) -3)]   #initialize fracdivefound matrix
	sdpiters=[[0 for x in range(maxinstances)] for x in range(len(sys.argv) -3)]   #initialize sdpiterations matrix
	penalties=[[0 for x in range(maxinstances)] for x in range(len(sys.argv) -3)]   #initialize penalty% matrix
	unsolved=[[0 for x in range(maxinstances)] for x in range(len(sys.argv) -3)]   #initialize unsolved% matrix
	racingwinner=[[0 for x in range(maxinstances)] for x in range(len(sys.argv) -3)]   #initialize racingwiner matrix
	sys.argv.remove("extractStatDrawPerfGraph_ug.py")				 #remove function call from list of files
	texfilename = sys.argv[0]
	file=open(texfilename, "w")
	file.write("\\documentclass[landscape]{article} \n \\usepackage{amsmath} \n \\usepackage{amsthm} \n \\usepackage{dsfont} \n \\usepackage[dvipsnames]{xcolor} \n \\usepackage{booktabs} \n \\usepackage{multirow} \n \\usepackage{mathtools} \n \\usepackage{longtable} \n \usepackage{tikz} \n \usepackage{pgfplots} \n  \usepackage[margin=1in,footskip=0.25in]{geometry} \n \\extrafloats{100} \n \\usepackage{sistyle} \n \\begin{document} \n")
	sys.argv.remove(sys.argv[0]) #remove results file
	assert(len(subsetnames) >= len(subsets))
	i=0
	#read SCIP-SDP results
	readFile(sys.argv[0],i)
	#if completeTable:
	#	makeCompleteTable(arg, file, i)
	if overviewTable:
		makeOverviewTable(arg, file, i)
	i=i+1
	#read UG results
	solufile = open(sys.argv[1], "r")
	solulines = solufile.readlines()
        sys.argv.remove(sys.argv[1]) #remove solu file
	for arg in sys.argv[1:]:#ignore SCIP-SDP-resultsfile
		j=0
		for soluline in solulines:
			readUG(i,j,arg,soluline)
			j=j+1
		if completeTable:
			#makeCompleteTable(arg, file, i)
			makeCompleteTableCaption(arg, i)
		i=i+1

	if plotsbysettingssubsets:
		for settings in settingsets:
			if overviewBySubset:
				ind = 0
				for (s,f) in subsets:
					makeOverviewBySubset(file, sys.argv, settings, s, f, subsetnames[ind], subsetnames[ind], settingnames, 0)
					ind = ind + 1
			if overviewTotal:
				makeOverviewBySubset(file, sys.argv, settings, 0, ninstances[0] - 1, "totalresults", "Overall", settingnames, 0)
			if performancePlotBySubset:
				ind = 0
				for (s,f) in subsets:
					makePerformancePlot(file, sys.argv, settings, s, f, ind)
					ind = ind + 1
			if performancePlotOverview:
				makePerformancePlot(file, sys.argv, settings, 0, ninstances[0] - 1, -1)
	else:
		if overviewBySubset:
			ind = 0
			for (s,f) in subsets:
				makeOverviewBySubset(file, sys.argv, range(len(sys.argv)), s, f, ind)
				ind = ind + 1
		if overviewTotal:
			makeOverviewBySubset(file, sys.argv, range(len(sys.argv)), 0, ninstances[0] - 1, -1)
		if performancePlotBySubset:
			ind = 0
			for (s,f) in subsets:
				makePerformancePlot(file, sys.argv, range(len(sys.argv)), s, f, ind)
				ind = ind + 1
		if performancePlotOverview:
			makePerformancePlot(file, sys.argv, range(len(sys.argv)), 0, ninstances[0] - 1, -1)

	if racingWinnerBySubset:
		ind = 0
		for (s,f) in subsets:
			makeRacingWinnerBySubset(file, sys.argv, settings, s, f, "racing wins " + subsetnames[ind], "racing wins " + subsetnames[ind], settingnames)
			ind = ind + 1
	if racingWinnerOverall:
		makeRacingWinnerBySubset(file, sys.argv, settings, 0, 193, "racing wins " + "overall", "racing wins " + "overall", settingnames)

        if racingWinnerBarplotAll:
                i = 1
                for arg in sys.argv[1:]:#ignore SCIP-SDP-resultsfile
                        makeRacingWinnerPlot(file, i, "racing wins for testrun number " + str(i), "racingwins" + str(i), settingnames)
                        i = i + 1

	file.write("\\end{document}")
	file.close()
