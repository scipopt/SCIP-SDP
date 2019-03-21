import os
import sys
import math
from decimal import Decimal
#usage: python extractStatDrawPerfGraph.py textfile.tex outfile1.out outfile2.out
#example: python extractStatDrawPerfGraph_roundingprob.py /home/gally/gally_Dissertation/resultfiles/WarmstartResults20180711/warmstartcomparison_SCIPSDPRIPpaper_roundingprob.tex  /home/gally/gally_Dissertation/resultfiles/WarmstartResults20180711/check.SCIPSDPRIP.scipsdp.linux.x86_64.gnu.opt.sdpa.moskito.warmstartipfactor05_proj4_slatercheck.out

#example python extractStatDrawPerfGraph_roundingprob.py /home/gally/gally_Dissertation/Supplement/Tables/roundingprobStats.tex /home/gally/gally_Dissertation/resultfiles/WarmstartResults20180711/check.SCIPSDPRIP.scipsdp.linux.x86_64.gnu.opt.sdpa.moskito.warmstartipfactor05_proj4_slatercheck.out /home/gally/gally_Dissertation/resultfiles/WarmstartResults20180711/check.SCIPSDPRIP.scipsdp.linux.x86_64.gnu.opt.sdpa.moskito.warmstartipfactor05_proj4_slatercheck_norr.out

#TODO: explicitly output number of solved relaxatiosn instead of summing up (which misses the failed checks)

maxinstances = 321 # maximum number of instances
completeTable = 1 # make table listing all instances
overviewTable = 0 # make table for each setting listing all subsets
overviewBySubset = 0 # make table for each subset listing all settings
overviewTotal = 0 # make table for complete testset listing all settings
roundingprobStatBySubset = 0 #make table for each subsets with statistics on rounding problem success and timing
roundingprobStatTotal = 0 #make table for whole testset with statistics on rounding problem success and timing
performancePlotOverview = 0 # draw a performance plot for the complete testset
performancePlotBySubset = 0 # draw a performance plot for each subset
subsets = [[0,64],[65,133],[134,193], [194,319]] # combine all instances from subsets[i][0] to subsets[i][1] to a single row and make a performance plot for each subset
#subsets = [[0,125]]
#subsets = [[0,2],[3,5],[6,8]]
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
nodegeomonlysolved = 0 # in the tables for each setting only include solved instances into the geometric mean of the nodes, for the overviews only include instances which were solved by all settings
overviewlongtable = 0 # should overview tables be done as longtables
plotsbysettingssubsets = 1 # make seperate plots and tables for each subset in settingssubsets
#settingsets = [[0,1,2,3,4,5], [6,7,8,9,10,11], [12,13,14,15,16,17], [18,19,20,21,22,23], [24,25,26,27,28,29], [30,31,32,33,34,35], [36,37,38,39,40,41], [42,43,44,45,46,47], [0,6,12,18,24,30,36,42], [1,7,13,19,25,31,37,43], [2,8,14,20,26,32,38,44], [3,9,15,21,27,33,39,45], [4,10,16,22,28,34,40,46], [5,11,17,23,29,35,41,47]]
#settingsets = [[0,2,3,4,5,11,17,23], [17,18,19,20,21,22], [23,24,25,26,27,28]]
settingsets = [[0]]
usesettingnames = 1
settingnames = [ "roundprob 0.5 id"]
epsilon = 0.001
captions=["Complete statistics for rounding problems", "Complete statistics for rounding problems with disabled randomized rounding"]
shortcaptions=["Complete statistics for rounding problems", "Complete statistics for rounding problems with disabled randomized rounding"]
labels=["roundingprobs", "roundingprobsNorandround"]


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
		if substring.find("Number of nodes that were successfully warmstarted using the rounding problems:") == -1:
			roundingsuccess[i][j] = "-"
		else:
			roundsucc = substring[substring.find("Number of nodes that were successfully warmstarted using the rounding problems:") + 80:].split()[0]
			roundingsuccess[i][j] = int(roundsucc)
		if substring.find("Number of nodes where the primal rounding problem failed:") == -1:
			roundingprimalfail[i][j] = "-"
		else:
			roundpfail = substring[substring.find("Number of nodes where the primal rounding problem failed:") + 58:].split()[0]
			roundingprimalfail[i][j] = int(roundpfail)
		if substring.find("Number of nodes where the dual rounding problem failed:") == -1:
			roundingdualfail[i][j] = "-"
		else:
			rounddfail = substring[substring.find("Number of nodes where the dual rounding problem failed:") + 56:].split()[0]
			roundingdualfail[i][j] = int(rounddfail)
		if substring.find("Time spent in rounding problems for warmstarting / detecting infeasibility:") == -1:
			roundingtime[i][j] = "-"
		else:
			roundtime = substring[substring.find("Time spent in rounding problems for warmstarting / detecting infeasibility:") + 76:].split()[0]
			roundingtime[i][j] = float(roundtime)
                totalsolved[i][j] = 0
		if substring.find("Number of infeasible nodes:") == -1:
			infeasible[i][j] = "-"
		else:
			infeas = substring[substring.find("Number of infeasible nodes:") + 28:].split()[0]
			infeasible[i][j] = int(infeas)
                        totalsolved[i][j] += infeasible[i][j]
		if substring.find("Number of nodes detected infeasible through primal rounding problem:") == -1:
			roundinginfeasible[i][j] = "-"
		else:
			roundinf = substring[substring.find("Number of nodes detected infeasible through primal rounding problem:") + 69:].split()[0]
			roundinginfeasible[i][j] = int(roundinf)
		if substring.find("Number of nodes where the optimal solution was determined by the rounding problem:") == -1:
			roundingoptimal[i][j] = "-"
		else:
			roundopt = substring[substring.find("Number of nodes where the optimal solution was determined by the rounding problem:") + 83:].split()[0]
			roundingoptimal[i][j] = int(roundopt)
		if substring.find("Number of nodes cut off through bounding by the rounding problem:") == -1:
			roundingcutoff[i][j] = "-"
		else:
			roundcutoff = substring[substring.find("Number of nodes cut off through bounding by the rounding problem:") + 66:].split()[0]
			roundingcutoff[i][j] = int(roundcutoff)
		if substring.find("Number of nodes with primal and dual slater holding:") != -1:
			pdholdingnodes = substring[substring.find("Number of nodes with primal and dual slater holding:") + 53:].split()[0]
			totalsolved[i][j] += int(pdholdingnodes)
		if substring.find("Number of nodes with either primal or dual slater not holding:") != -1:
			pdfailingnodes = substring[substring.find("Number of nodes with either primal or dual slater not holding:") + 63:].split()[0]
			totalsolved[i][j] += int(pdfailingnodes)
		filestring = filestring[filestring.find("@04")+36:]
		j = j+1
		assert( j < maxinstances)
	ninstances[i] = j
	file.close()

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
	file.write(" & \\multicolumn{5}{c}{all rounding problems} & \\multicolumn{2}{c}{infeasibility} \\\ \n ")
	file.write("problem & opt & pruned & warmstart & pfail & dfail & detected & total ")
	file.write("\\\ \\midrule} \n \\tabletail{ \\midrule \\multicolumn{3}{@{}l}{continued on next page \\dots}\\\ \\bottomrule} \n \\tablelasttail{\\bottomrule} \n \\tablecaption[")
	file.write(shortcaptions[i])
	file.write("]{")
	file.write(captions[i])
	file.write("}\label{")
	file.write(labels[i])
	file.write("}\n")
	file.write("\\begin{xtabular*}{\\textwidth}{@{\extracolsep{\\fill}}lrrrrrrr")
	file.write("@{}} \n ")
	for j in range(ninstances[i]):
		file.write(names[i][j].replace("_", "\\_").split(".")[0] + "& ")
		if roundingoptimal[i][j] == "-":
			file.write("-- & ")
		else:
			file.write("\\num{" + "%.0f" % float(roundingoptimal[i][j]) + "} &")
		if roundingcutoff[i][j] == "-":
			file.write("-- & ")
		else:
			file.write("\\num{" + "%.0f" % float(roundingcutoff[i][j]) + "} &")
		if roundingsuccess[i][j] == "-":
			file.write("-- & ")
		else:
			file.write("\\num{" + "%.0f" % float(roundingsuccess[i][j]) + "} &")
		if roundingprimalfail[i][j] == "-":
			file.write("-- & ")
		else:
			file.write("\\num{" + "%.0f" % float(roundingprimalfail[i][j]) + "} &")
		if roundingdualfail[i][j] == "-":
			file.write("-- & ")
		else:
			file.write("\\num{" + "%.0f" % float(roundingdualfail[i][j]) + "} &")
		if roundinginfeasible[i][j] == "-":
			file.write("-- & ")
		else:
			file.write("\\num{" + "%.0f" % float(roundinginfeasible[i][j]) + "} &")
		if infeasible[i][j] == "-":
			file.write("-- ")
		else:
			file.write("\\num{" + "%.0f" % float(infeasible[i][j]) + "}")
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
		file.write("\\begin{longtable}[ht] {p{.3\\textwidth} | p{.05\\textwidth} p{.09\\textwidth} ")
		if fixingscol:
			file.write("p{.09\\textwidth} ")
		file.write("p{.05\\textwidth} p{.05\\textwidth} p{.09\\textwidth} p{.09\\textwidth} p{.05\\textwidth}} \\toprule \n settings &  solved & nodes & ")
		if fixingscol:
			file.write("fixings & ")
		file.write("fracdive & time & sdpiter & penalty & unsolved \% \\\ \midrule \n")
	else:
		file.write("\\begin{table} \n \\begin{scriptsize} \\caption{" + caption + "} \n \\label{" + label + "} \n \\begin{tabular*}{\\textwidth}{@{}l@{\\;\\;\extracolsep{\\fill}}rrrrrrr")
		if fixingscol:
			file.write("r")
		file.write("@{}}\\toprule \n")
		file.write(" settings &  solved & nodes & ")
		if fixingscol:
			file.write("fixings & ")
		file.write("fracdive & time & sdpiter & penalty & unsolved \\\ \midrule \n")
	i = 0
	ind = 0
	argind = -1
	for arg in args:
		if i in settings:
			nsolved = 0
			nodegeom = 1.0
			if fixingsarithmetic:
				totalfixings = 0.0
			else:
				totalfixings = 1.0
			if heurarithmetic:
				totalheur = 0.0
				nsolved = 0
			else:
				totalheur = 1.0
			timegeom = 1.0
			itergeom = 1.0
			totalpenalty = 0
			totalunsolved = 0
                        nunaborted = 0
                        nnodescounted = 0
                        niterscounted = 0
			j = s
			numinstances = f - s + 1
			while j <= f:
                                if nodes[i][j] != "-":
                                        nunaborted = nunaborted + 1
				if nodes[i][j] != "-" and (not nodegeomonlysolved or allSolved(j, settings)):
                                        nnodescounted = nnodescounted + 1
                                if sdpiters[i][j] != "-" and (not itergeomonlysolved or allSolved(j, settings)):
                                        niterscounted = niterscounted + 1
                                j = j + 1
                        j = s
			while j <= f:
				if instanceSolved(i,j) :
					nsolved = nsolved + 1
				if nodes[i][j] != "-" and (not nodegeomonlysolved or allSolved(j, settings)):
					nodegeom = nodegeom * math.pow(float(nodes[i][j] + nodeshift), 1.0/float(nnodescounted))
				if fixingsarithmetic and sdpredcostfixings[i][j] != "-":
					totalfixings = totalfixings + float(sdpredcostfixings[i][j])
				elif sdpredcostfixings[i][j] != "-":
					totalfixings = totalfixings * math.pow(float(sdpredcostfixings[i][j] + fixingsshift), 1.0/float(nunaborted))
				if heurarithmetic and fracdivefound[i][j] != "-":
					totalheur = totalheur + float(fracdivefound[i][j])
				elif fracdivefound[i][j] != "-":
					totalheur = totalheur * math.pow(float(fracdivefound[i][j] + heurshift), 1.0/float(nunaborted))
				timegeom = timegeom * math.pow(float(times[i][j] + timeshift), 1.0/float(f-s+1))
				if sdpiters[i][j] != "-" and (not itergeomonlysolved or allSolved(j, settings)):
					itergeom = itergeom * math.pow(float(sdpiters[i][j] + itershift), 1.0/float(niterscounted))
				if penalties[i][j] != "-":
					totalpenalty = totalpenalty + penalties[i][j]
				if unsolved[i][j] != "-":
					totalunsolved = totalunsolved + unsolved[i][j]
				j = j + 1

			nodegeom = nodegeom - nodeshift
			if fixingsarithmetic:
				avgfixings = totalfixings / nunaborted
			else:
				avgfixings = totalfixings - fixingsshift
				if avgfixings < epsilon:
					avgfixings = 0.0
			if heurarithmetic:
				avgheur = totalheur / nunaborted
			else:
				avgheur = totalheur - heurshift
			timegeom = timegeom - timeshift
			itergeom = itergeom - itershift
			avgpenalty = float(totalpenalty) / float(nunaborted)
			avgunsolved = float(totalunsolved) / float(nunaborted)
			file.write(settingnames[ind] + " & \\num{" + "%.0f" % nsolved + "} & \\num{" + "%.2f" % nodegeom + "} & ")
			if fixingscol:
				file.write("\\num{" + "%.2f" % avgfixings + "} & ")
			file.write("\\num{" + "%.2f" % avgheur + "} & \\num{" + "%.2f" % timegeom + "} & \\num{" + "%.2f" % itergeom + "} & \\num{" + "%.2f" % avgpenalty + "} \% & \\num{" + "%.2f" % avgunsolved + "} \% \\\ \n")
			ind = ind + 1
		i = i + 1
	if overviewlongtable:
		file.write("\\bottomrule \n \\caption{" + caption + "} \n \\label{" + label + "} \n \\end{longtable} \n")
	else:
		file.write("\\bottomrule \n \\end{tabular*} \n \end{scriptsize} \n \\end{table} \n")

def makeRoundingStatBySubset(file, args, settings, s, f, caption, label, settingnames, fixingscol):
	if overviewlongtable:
		file.write("\\begin{longtable}[ht] {p{.3\\textwidth} | p{.05\\textwidth} p{.09\\textwidth} ")
		if fixingscol:
			file.write("p{.09\\textwidth} ")
		file.write("p{.05\\textwidth} p{.05\\textwidth} p{.09\\textwidth} p{.09\\textwidth} p{.05\\textwidth}} \\toprule \n settings &  solved & nodes & ")
		if fixingscol:
			file.write("fixings & ")
		file.write("fracdive & time & sdpiter & penalty & unsolved \% \\\ \midrule \n")
	else:
		file.write("\\begin{table} \n \\begin{scriptsize} \\caption{" + caption + "} \n \\label{" + label + "} \n \\begin{tabular*}{\\textwidth}{@{}l@{\\;\\;\extracolsep{\\fill}}rrrrrrrrrrrr")
		file.write("@{}}\\toprule \n")
		file.write(" settings &  solved & time & roundtime & sdpiter & ipsolved & \multicolumn{5}{c}{roundingprob} & \multicolumn{2}{c}{infeasibility} \\\  \n")
                file.write(" & & & & & & opt & cutoff & warmstart & pfail & dfail & detected & total \\\ \midrule \n")
	i = 0
	ind = 0
	argind = -1
	for arg in args:
		if i in settings:
			timegeom = 1.0
			itergeom = 1.0
                        nsolved = 0
                        nunaborted = 0
                        nnodescounted = 0
                        niterscounted = 0
                        roundingsuccessarith=0
                        roundingprimalfailarith=0
                        roundingdualfailarith=0
                        roundingtimegeom=1.0
                        infeasiblearith=0
                        roundinginfeasiblearith=0
                        totalrelaxarith=0
                        roundcutoffarith=0
                        roundoptarith=0
			j = s
			numinstances = f - s + 1
			while j <= f:
                                if nodes[i][j] != "-":
                                        nunaborted = nunaborted + 1
				if nodes[i][j] != "-" and (not nodegeomonlysolved or allSolved(j, settings)):
                                        nnodescounted = nnodescounted + 1
                                if sdpiters[i][j] != "-" and (not itergeomonlysolved or allSolved(j, settings)):
                                        niterscounted = niterscounted + 1
                                j = j + 1
                        j = s
			while j <= f:
				if instanceSolved(i,j) :
					nsolved = nsolved + 1
				timegeom = timegeom * math.pow(float(times[i][j] + timeshift), 1.0/float(f-s+1))
                                if nodes[i][j] != "-":
                                        roundingtimegeom = roundingtimegeom * math.pow(float(roundingtime[i][j] + timeshift), 1.0/float(nunaborted))
				if sdpiters[i][j] != "-" and (not itergeomonlysolved or allSolved(j, settings)):
					itergeom = itergeom * math.pow(float(sdpiters[i][j] + itershift), 1.0/float(niterscounted))
                                if nodes[i][j] != "-":
                                        roundingsuccessarith = roundingsuccessarith + roundingsuccess[i][j]
                                        roundingprimalfailarith = roundingprimalfailarith + roundingprimalfail[i][j]
                                        roundingdualfailarith = roundingdualfailarith + roundingdualfail[i][j]
                                        if (infeasible[i][j] != "-"):
                                                infeasiblearith = infeasiblearith + infeasible[i][j]
                                        roundinginfeasiblearith = roundinginfeasiblearith + roundinginfeasible[i][j]
                                        totalrelaxarith = totalrelaxarith + totalsolved[i][j]
                                        roundcutoffarith = roundcutoffarith + roundingcutoff[i][j]
                                        roundoptarith = roundoptarith + roundingoptimal[i][j]
				j = j + 1

			timegeom = timegeom - timeshift
                        roundingtimegeom = roundingtimegeom - timeshift
			itergeom = itergeom - itershift
			roundingsuccessarith = float(roundingsuccessarith) / float(nunaborted)
			roundingprimalfailarith = float(roundingprimalfailarith) / float(nunaborted)
                        roundingdualfailarith = float(roundingdualfailarith) / float(nunaborted)
                        infeasiblearith = float(infeasiblearith) / float(nunaborted)
                        roundinginfeasiblearith = float(roundinginfeasiblearith) / float(nunaborted)
                        totalrelaxarith = float(totalrelaxarith) / float(nunaborted)
                        roundcutoffarith = float(roundcutoffarith) / float(nunaborted)
                        roundoptarith = float(roundoptarith) / float(nunaborted)
			file.write(settingnames[ind] + " & \\num{" + "%.0f" % nsolved + "} & \\num{" + "%.2f" % timegeom + "} & \\num{" + "%.2f" % roundingtimegeom + "} & \\num{" + "%.2f" % itergeom + "} & \\num{" + "%.2f" % totalrelaxarith + "} & \\num{" + "%.2f" % roundoptarith + "} & \\num{" + "%.2f" % roundcutoffarith + "} & \\num{" + "%.2f" % roundingsuccessarith + "} & \\num{" + "%.2f" % roundingprimalfailarith + "} & \\num{" + "%.2f" % roundingdualfailarith + "} & \\num{" + "%.2f" % roundinginfeasiblearith + "} & \\num{" + "%.2f" % infeasiblearith + "}  \\\ \n")
			ind = ind + 1
		i = i + 1
	if overviewlongtable:
		file.write("\\bottomrule \n \\caption{" + caption + "} \n \\label{" + label + "} \n \\end{longtable} \n")
	else:
		file.write("\\bottomrule \n \\end{tabular*} \n \end{scriptsize} \n \\end{table} \n")

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
				performance[arg][j-s] = times[arg][j] / best
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
	names=[[0 for x in range(maxinstances)] for x in range(len(sys.argv) -2)] #initialize names matrix
	dualresults=[[0 for x in range(maxinstances)] for x in range(len(sys.argv) -2)] #initialize dual results matrix
	primalresults=[[0 for x in range(maxinstances)] for x in range(len(sys.argv) -2)] #initialize primal results matrix
	gaps=[[0 for x in range(maxinstances)] for x in range(len(sys.argv) -2)]   #initialize gaps matrix
	times=[[0 for x in range(maxinstances)] for x in range(len(sys.argv) -2)]   #initialize solvingtime matrix
	nodes=[[0 for x in range(maxinstances)] for x in range(len(sys.argv) -2)]   #initialize nodes matrix
	sdpredcostfixings=[[0 for x in range(maxinstances)] for x in range(len(sys.argv) -2)]   #initialize redcostfixings matrix
	fracdivefound=[[0 for x in range(maxinstances)] for x in range(len(sys.argv) -2)]   #initialize fracdivefound matrix
	sdpiters=[[0 for x in range(maxinstances)] for x in range(len(sys.argv) -2)]   #initialize sdpiterations matrix
	penalties=[[0 for x in range(maxinstances)] for x in range(len(sys.argv) -2)]   #initialize penalty% matrix
	unsolved=[[0 for x in range(maxinstances)] for x in range(len(sys.argv) -2)]   #initialize unsolved% matrix
        roundingsuccess=[[0 for x in range(maxinstances)] for x in range(len(sys.argv) -2)]
        roundingprimalfail=[[0 for x in range(maxinstances)] for x in range(len(sys.argv) -2)]
        roundingdualfail=[[0 for x in range(maxinstances)] for x in range(len(sys.argv) -2)]
        roundingtime=[[0 for x in range(maxinstances)] for x in range(len(sys.argv) -2)]
        infeasible=[[0 for x in range(maxinstances)] for x in range(len(sys.argv) -2)]
        roundinginfeasible=[[0 for x in range(maxinstances)] for x in range(len(sys.argv) -2)]
        roundingcutoff=[[0 for x in range(maxinstances)] for x in range(len(sys.argv) -2)]
        roundingoptimal=[[0 for x in range(maxinstances)] for x in range(len(sys.argv) -2)]
        totalsolved=[[0 for x in range(maxinstances)] for x in range(len(sys.argv) -2)]
	sys.argv.remove("extractStatDrawPerfGraph_roundingprob.py")				 #remove function call from list of files
	texfilename = sys.argv[0]
	file=open(texfilename, "w")
	file.write("\\documentclass[landscape]{article} \n \\usepackage{amsmath} \n \\usepackage{amsthm} \n \\usepackage{dsfont} \n \\usepackage[dvipsnames]{xcolor} \n \\usepackage{booktabs} \n \\usepackage{multirow} \n \\usepackage{mathtools} \n \\usepackage{longtable} \n \usepackage{tikz} \n \usepackage{pgfplots} \n  \usepackage[margin=1in,footskip=0.25in]{geometry} \n \\extrafloats{100} \n \\usepackage{sistyle} \n \\begin{document} \n")
	sys.argv.remove(sys.argv[0])
	assert(len(subsetnames) >= len(subsets))
	i=0
	for arg in sys.argv:
		readFile(arg,i)
		print(str(i + 1) + "/" + str(len(sys.argv)) + " files read")
		if completeTable:
			#makeCompleteTable(arg, file, i)
			makeCompleteTableCaption(arg, i)
		if overviewTable:
			makeOverviewTable(arg, file, i)
		i=i+1
	# if plotsbysettingssubsets:
	# 	for settings in settingsets:
	# 		if overviewBySubset:
	# 			ind = 0
	# 			for (s,f) in subsets:
	# 				makeOverviewBySubset(file, sys.argv, settings, s, f, subsetnames[ind], subsetnames[ind], settingnames, 0)
	# 				ind = ind + 1
	# 		if overviewTotal:
	# 			makeOverviewBySubset(file, sys.argv, settings, 0, ninstances[0] - 1, "totalresults", "Overall", settingnames, 0)
	# 		if performancePlotBySubset:
	# 			ind = 0
	# 			for (s,f) in subsets:
	# 				makePerformancePlot(file, sys.argv, settings, s, f, ind)
	# 				ind = ind + 1
	# 		if performancePlotOverview:
	# 			makePerformancePlot(file, sys.argv, settings, 0, ninstances[0] - 1, -1)
    #                     if roundingprobStatBySubset:
    #                             ind = 0
	# 			for (s,f) in subsets:
	# 				makeRoundingStatBySubset(file, sys.argv, settings, s, f, subsetnames[ind], subsetnames[ind], settingnames, 0)
	# 				ind = ind + 1
    #                     if roundingprobStatTotal:
    #                             makeRoundingStatBySubset(file, sys.argv, settings, 0, ninstances[0] - 1, "totalresults", "Overall", settingnames, 0)
	# else:
	# 	if overviewBySubset:
	# 		ind = 0
	# 		for (s,f) in subsets:
	# 			makeOverviewBySubset(file, sys.argv, range(len(sys.argv)), s, f, ind)
	# 			ind = ind + 1
	# 	if overviewTotal:
	# 		makeOverviewBySubset(file, sys.argv, range(len(sys.argv)), 0, ninstances[0] - 1, -1)
	# 	if performancePlotBySubset:
	# 		ind = 0
	# 		for (s,f) in subsets:
	# 			makePerformancePlot(file, sys.argv, range(len(sys.argv)), s, f, ind)
	# 			ind = ind + 1
	# 	if performancePlotOverview:
	# 		makePerformancePlot(file, sys.argv, range(len(sys.argv)), 0, ninstances[0] - 1, -1)
	file.write("\\end{document}")
	file.close()
