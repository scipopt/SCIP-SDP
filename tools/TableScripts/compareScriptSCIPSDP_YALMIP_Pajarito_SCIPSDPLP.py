import os
import sys
import math
from decimal import Decimal

#first argument should be the tex-file to write to, second should be SCIP-SDP-results file,  third SCIP-SDP-lpapproxsecond, fourth yalmip-results-file, fifth YALMIP CUTSDP, sixth Pajarito,
#example for a correct call: python compareScriptSCIPSDP_YALMIP_Pajarito_SCIPSDPLP.py /home/gally/gally_Dissertation/resultfiles/MISDPvergleich/MISDPvergleich.tex /home/gally/gally_Dissertation/resultfiles/MISDPvergleich/check.SCIPSDPpaper.scipsdp.linux.x86_64.gnu.opt.msk.fb04668.1thread.out /home/gally/gally_Dissertation/resultfiles/MISDPvergleich/check.SCIPSDPpaper.scipsdp.linux.x86_64.gnu.opt.msk.fb04668.lp_approx.out /home/gally/gally_Dissertation/resultfiles/MISDPvergleich/SCIPSDPpaper_MATLAB2017a_YALMIPBNB_R20180926_MOSEK81054.test.results /home/gally/gally_Dissertation/resultfiles/MISDPvergleich/SCIPSDPpaper_MATLAB2017a_YALMIPCUTSDP_R20180926_CPLEX1261.test.results /home/gally/gally_Dissertation/resultfiles/MISDPvergleich/SCIPSDPpaper_CBF_drives_Pajarito050_Mosek81054_CPLEX1261_tolerances.test.pajaritoresults

#nodes: python compareScriptSCIPSDP_YALMIP_Pajarito_SCIPSDPLP.py /home/gally/gally_Dissertation/Supplement/Tables/MISDP.tex /home/gally/gally_Dissertation/resultfiles/MISDPvergleich/check.SCIPSDPpaper.scipsdp.linux.x86_64.gnu.opt.msk.fb04668.1thread.out /home/gally/gally_Dissertation/resultfiles/MISDPvergleich/check.SCIPSDPpaper.scipsdp.linux.x86_64.gnu.opt.msk.fb04668.lp_approx.out /home/gally/gally_Dissertation/resultfiles/MISDPvergleich/SCIPSDPpaper_MATLAB2017a_YALMIPBNB_R20180926_MOSEK81054.test.results /home/gally/gally_Dissertation/resultfiles/MISDPvergleich/SCIPSDPpaper_MATLAB2017a_YALMIPBNB_R20180926_MOSEK81054_noRounder.test.results /home/gally/gally_Dissertation/resultfiles/MISDPvergleich/YALMIP_20180612/SCIPSDPpaper_MATLAB2017a_YALMIPBNB_R20180612_MOSEK81054.test.results /home/gally/gally_Dissertation/resultfiles/MISDPvergleich/SCIPSDPpaper_MATLAB2017a_YALMIPCUTSDP_R20180926_CPLEX1261.test.results /home/gally/gally_Dissertation/resultfiles/MISDPvergleich/YALMIP_20180612/SCIPSDPpaper_MATLAB2017a_YALMIPCUTSDP_R20180612_CPLEX1261.test.results /home/gally/gally_Dissertation/resultfiles/MISDPvergleich/SCIPSDPpaper_CBF_drives_Pajarito050_Mosek81054_CPLEX1261.test.pajaritoresults

maxinstances = 200 # maximum number of instances
completeTable = 1 # make table listing all instances
#overviewTable = 0 # make table for each setting listing all subsets
#overviewBySubset = 1 # make table for each subset listing all settings
#overviewTotal = 1 # make table for complete testset listing all settings
#performancePlotOverview = 1 # draw a performance plot for the complete testset
#performancePlotBySubset = 1 # draw a performance plot for each subset
subsets = [[0,64],[65,133],[134,193]] # combine all instances from subsets[i][0] to subsets[i][1] to a single row and make a performance plot for each subset
#colored performance plots:
colors = ["SkyBlue, thick", "blue, thick", "Apricot, thick", "red, thick", "gray, thick", "Sepia, thick", "green, thick", "PineGreen, thick", "Rhodamine, thick", "DarkOrchid, thick",  "Sepia, thick", "black", "LimeGreen", "CadetBlue", "Lavender", "Tan", "SeaGreen", "Bittersweet", "BlueViolet", "Cerulean", "Brown", "Cyan", "ForestGreen", "Goldenrod",  "JungleGreen", "Mahogany", "Melon", "Mulberry", "OliveGreen", "OrangeRed", "Peach", "yellow", "ProcessBlue", "RawSienna", "RedOrange", "Aquamarine", "RoyalPurple", "Salmon", "SpringGreen", "TealBlue", "Turquoise", "VioletRed", "WildStrawberry", "YellowGreen", "BlueGreen", "BrickRed", "BurntOrange", "CarnationPink", "CornflowerBlue", "Dandelion", "Emerald", "Fuchsia", "GreenYellow", "Magenta", "Maroon", "MidnightBlue", "NavyBlue", "Orange", "Orchid", "Periwinkle", "Plum", "Purple", "RedViolet", "RoyalBlue", "RubineRed", "Thistle", "Violet", "YellowOrange"]
#black and white for different solvers
#colors = ["black, solid, thick", "black, dashed, thick", "black, dotted, thick", "black, dash pattern=on 1pt off 2pt on 3pt off 2pt, thick", "gray, solid, thick", "gray, dashed, thick", "gray, dotted, thick", "gray, dash pattern=on 1pt off 2pt on 3pt off 2pt, thick", "black, dash pattern=on 1pt off 2pt on 1pt off 2pt on 3pt off 2pt on 3pt off 2pt, thick", "gray, dash pattern=on 1pt off 2pt on 1pt off 2pt on 3pt off 2pt on 3pt off 2pt, thick", "red", "red", "red", "red", "red", "red", "red", "red", "red", "red", "red", "red", "red", "red", "red", "red", "red", "red","red", "red", "red", "red", "red", "red","red", "red", "red", "red", "red", "red","red", "red", "red", "red", "red", "red","red", "red", "red", "red", "red", "red",]
#black and white for different settings
#colors = ["very thick, 3d", "very thick, 10d", "very thick, 4a", "10a, solid, very thick", "very thick, 8b", "red", "red", "red", "red", "red", "red", "red", "red", "red", "red", "red", "red", "red", "red", "red", "red", "red", "red","red", "red", "red", "red", "red", "red","red", "red", "red", "red", "red", "red","red", "red", "red", "red", "red", "red","red", "red", "red", "red", "red", "red",]
subsetnames = ["cardinality constrained least squares", "min-$k$-partitioning", "truss topology"] # names given to the subsets used in the tables, needs to have at least as many elements as there are subsets
timeshift = 10
nodeshift = 100
itershift = 1000
fixingsarithmetic = 1
fixingsshift = 100
heurarithmetic = 1
heurshift = 1
texfile = 1
nodegeomonlysolved = 1 # in the tables for each setting only include solved instances into the geometric mean of the nodes, for the overviews only include instances which were solved by all settings
itergeomonlysolved = 1 # in the tables for each setting only include solved instances into the geometric mean of the nodes, for the overviews only include instances which were solved by all settings
overviewlongtable = 0 # should overview tables be done as longtables
usesettingnames = 1
settingnames = ["YALMIP", "SCIP-SDP"]
epsilon = 0.001
# without randrounding
#completetablesettings = [0,2,3,4,5,6,7,8,9,10,16,17,18,19,20,21,22,28]
#completeTableCaptions = ["Complete results for DSDP with inference branching", "", "Complete results for DSDP with combined infeasibility/objective branching", "Complete results for DSDP with combined infeasibility/objective branching and dualfixing", "Complete results for DSDP with combined infeasibility/objective branching and fractional diving in all nodes with depth a multiple of 10", "Complete results for DSDP with combined infeasibility/objective branching and dual fixing and fractional diving in all nodes with depth a multiple of 10", "Complete results for DSDP with combined infeasibility/objective branching and without fractional diving", "Complete results for DSDP with combined infeasibility/objective branching and dual fixing and without fractional diving", "Complete results for DSDP with objective branching", "Complete results for DSDP with infeasibility branching", "Complete results for SDPA with inference branching", "", "", "", "", "", "Complete results for SDPA with combined infeasibility/objective branching", "Complete results for SDPA with combined infeasibility/objective branching and dualfixing", "Complete results for SDPA with combined infeasibility/objective branching and fractional diving in all nodes with depth a multiple of 10", "Complete results for SDPA with combined infeasibility/objective branching and dual fixing and fractional diving in all nodes with depth a multiple of 10", "Complete results for SDPA with combined infeasibility/objective branching and without fractional diving", "Complete results for SDPA with combined infeasibility/objective branching and dual fixing and without fractional diving", "Complete results for SDPA with objective branching", "", "", "", "", "", "Complete results for SDPA with infeasibility branching", "Complete results for SDPA with infeasibility branching and dualfixing", "Complete results for SDPA with infeasibility branching and fractional diving in all nodes with depth a multiple of 10", "Complete results for SDPA with infeasibility branching and dual fixing and fractional diving in all nodes with depth a multiple of 10", "Complete results for SDPA with infeasibility branching and without fractional diving", "Complete results for SDPA with infeasibility branching and dual fixing and without fractional diving"]
#completeTableLabels = ["CompleteDSDPinfer", "", "CompleteDSDPinfobj", "CompleteDSDPinfobjdualfix", "CompleteDSDPinfobjdive10", "CompleteDSDPinfobjdive10dualfix", "CompleteDSDPinfobjnodive", "CompleteDSDPinfobjnodivedualfix", "CompleteDSDPobj", "CompleteDSDPinf", "CompleteSDPAinfer", "", "", "", "", "", "CompleteSDPAinfobj", "CompleteSDPAinfobjdualfix", "CompleteSDPAinfobjdive10", "CompleteSDPAinfobjdive10dualfix", "CompleteSDPAinfobjnodive", "CompleteSDPAinfobjnodivedualfix", "CompleteSDPAobj", "", "", "", "", "", "CompleteSDPAinf", "CompleteSDPAinfdualfix", "CompleteSDPAinfdive10", "CompleteSDPAinfdive10dualfix", "CompleteSDPAinfnodive", "CompleteSDPAinfnodivedualfix"]
#with randrounding (new = 8,9,10,11,26,27,28,29)
completetablesettings = [0,1,2]
completeTableCaptions = ["Complete results and performance indicators for YALMIP", "Complete results and performance indicators for Pajarito", "Complete results and performance indicators for SCIP-SDP"]
completeTableShortCaptions = ["Complete results and performance indicators for YALMIP", "Complete results and performance indicators for Pajarito", "Complete results and performance indicators for SCIP-SDP"]
completeTableLabels = ["CompleteYALMIP", "CompletePajarito", "CompleteSCIPSDP"]
captions=["Complete results for \scipsdp with nonlinear branch-and-bound on 8-core Intel i7-4770 CPU with \SI{3.40}{GHz} and \SI{16}{GB} memory", "Complete results for \scipsdp with LP-based cutting plane approach on 8-core Intel i7-4770 CPU with \SI{3.40}{GHz} and \SI{16}{GB} memory","Complete results for \\yalmip-\\bnb on 8-core Intel i7-4770 CPU with \SI{3.40}{GHz} and \SI{16}{GB} memory","Complete results for \\yalmip-\\bnb without rounding on 8-core Intel i7-4770 CPU with \SI{3.40}{GHz} and \SI{16}{GB} memory","Complete results for \\yalmip-\\bnb R20180612 on 8-core Intel i7-4770 CPU with \SI{3.40}{GHz} and \SI{16}{GB} memory","Complete results for \\yalmip-\\cutsdp on 8-core Intel i7-4770 CPU with \SI{3.40}{GHz} and \SI{16}{GB} memory","Complete results for \\yalmip-\\cutsdp R20180612 on 8-core Intel i7-4770 CPU with \SI{3.40}{GHz} and \SI{16}{GB} memory","Complete results for \\pajarito on 8-core Intel i7-4770 CPU with \SI{3.40}{GHz} and \SI{16}{GB} memory"]
shortcaptions=["Complete results for \scipsdp with nonlinear branch-and-bound", "Complete results for \scipsdp with LP-based cutting plane approach", "Complete results for \\yalmip-\\bnb", "Complete results for \\yalmip-\\bnb without rounding", "Complete results for \\yalmip-\\bnb R20180612", "Complete results for \\yalmip-\\cutsdp","Complete results for \\yalmip-\\cutsdp R20180612", "Complete results for \\pajarito"]
labels=["SCIPSDP-NLBB", "SCIPSDP-LP", "YALMIPBNB", "YALMIPBNBnoRounder", "YALMIPBNB20180612", "YALMIPCUTSDP", "YALMIPCUTSDP20180612", "Pajarito"]

def readFile(arg,i):
	j=0;
	file=open(arg, "r")
	filestring = file.read()
	while filestring.find("@01") > -1:
		substring = filestring[filestring.find("@01"):filestring.find("@04")]
		assert(substring.find("@01 ") > -1)
		names[i][j] = os.path.basename(substring[substring.find("@01 ") + 3:].split()[0])
		if substring.find("Solving Time (sec) : ") == -1:
			times[i][j] = 3600.0
		else:
			time = substring[substring.find("Solving Time (sec) : ") + 21:].split()[0]
			times[i][j] = float(time)
			if times[i][j] > 3600.0:
				times[i][j] = 3600.0
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
		if substring.find("sdprand          :") == -1:
			randfound[i][j] = "-"
		else:
			rand = substring[substring.find("sdprand          :") + 18:].split()[3]
			randfound[i][j] = int(rand)
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
			gaps[i][j] = "infinite"
		else:
			gap = substring[substring.find("Gap                : ") + 21:].split()[0]
			gaps[i][j] = gap
		if substring.find("SDP iterations:") == -1:
			sdpiters[i][j] = "-"
		else:
			sdpiter = substring[substring.find("SDP iterations:") + 15:].split()[0]
			sdpiters[i][j] = int(sdpiter)
		if substring.find("Percentage penalty formulation used:") == -1:
			penalties[i][j] = "-"
		else:
			penalty = substring[substring.find("Percentage penalty formulation used:") + 36:].split()[0]
			penalties[i][j] = float(penalty)
		if substring.find("Percentage unsolved even with penalty:") == -1:
			unsolved[i][j] = "-"
		else:
			unsolv = substring[substring.find("Percentage unsolved even with penalty:") + 38:].split()[0]
			unsolved[i][j] = float(unsolv)
		filestring = filestring[filestring.find("@04")+3:]
		j = j+1
		assert( j < maxinstances)
	ninstances[i] = j
	assert(ninstances[i] == ninstances[0])
	file.close()

def readYALMIP(arg, i):
	j = 0;
	file=open(arg, "r")
	filestring = file.read()
	lines = filestring.split("\n")
	iterlines = iter(lines)
	#next(iterlines) #skip the header line
	for l in iterlines:
		if l == "":
			break #for some reason python has an empty element at the end of the list
		words = l.split()
		names[i][j] = os.path.basename(words[0].split()[0])
		if words[1] != "error":
			primalresults[i][j] = float(words[1])
			if float(words[2]) < 3600:
				times[i][j] = float(words[2])
			else:
				times[i][j] = 3600.0
			nodes[i][j] = int(words[3])
		else:
			primalresults[i][j] = "error"
                        times[i][j] = 3600.0
		j = j+1
		assert( j < maxinstances)
	ninstances[i] = j
	assert(ninstances[i] == ninstances[0])
	file.close()

def readPajarito(arg, i):
	j = 0;
	file=open(arg, "r")
	filestring = file.read()
	lines = filestring.split("\n")
	iterlines = iter(lines)
	for l in iterlines:
		if l == "":
			break #for some reason python has an empty element at the end of the list
		words = l.split()
		names[i][j] = os.path.basename(words[0].split()[0])
		if words[1] == "Optimal" or words[1] == "Suboptimal":
			primalresults[i][j] = float(words[2])
			if float(words[3]) < 3600:
				times[i][j] = float(words[3])
			else:
				times[i][j] = 3600.0
			nodes[i][j] = int(words[5])
		else:
			primalresults[i][j] = "error"
			times[i][j] = 3600.0
		j = j+1
		assert( j < maxinstances)
	ninstances[i] = j
	assert(ninstances[i] == ninstances[0])
	file.close()

def instanceSolved(i, j):
	""" check if the duality gap is zero """
	if times[i][j] > 3599.999:
		return False
	return True

def oneSolved(ind, settings):
	""" check if for this index at least one setting could solve the instance to optimality"""
	for i in settings:
		if instanceSolved(i, ind):
			return True
	return False

def allSolved(ind, settings):
	""" check if for this index all settings could solve the instance"""
	settings=[0,1,2] #############<- hardcoded since the others use different kinds of "nodes" (and otherwise we have too few instances)
	for i in settings:
		if not instanceSolved(i, ind):
			return False
	return True

def optimalSolution(ind, settings):
	""" if at least one of settings could solve instance ind to optimality this returns the optimal objective"""
	if not oneSolved(ind, settings):
		return -float("inf")
	for i in settings:
		if times[i][ind] <= 3599.999:
			return primalresults[i][ind]
	return -float("inf")

def BestPrimal(ind, settings):
	""" returns minimum of primalbounds of all settings for instance ind"""
	best = float("inf")
	for i in settings:
		if primalresults[i][ind] < best:
			best = primalresults[i][ind]
	return best

def BestDual(ind, settings):
	""" returns maximum of dualbounds of all settings for instance ind"""
	best = -float("inf")
	for i in settings:
		if dualresults[i][ind] > best and not dualresults[i][ind] == 1e+20:
			best = dualresults[i][ind]
	return best

def makeCompleteTableCaption(arg, i, SCIP, SDP):
	file.write("\\newpage \n \\begin{scriptsize} \n  \\setlength{\\tabcolsep}{2pt} \n \\tablehead{\\toprule \n")
	if SCIP and SDP:
		file.write("problem & dbound &  pbound & gap & nodes & time & iters & pen & uns")
	elif SCIP:
		file.write("problem & dbound &  pbound & gap & nodes & time")
	else:
		file.write("problem & pbound & nodes & time")
	file.write("\\\ \\midrule} \n \\tabletail{ \\midrule \\multicolumn{3}{@{}l}{continued on next page \\dots}\\\ \\bottomrule} \n \\tablelasttail{\\bottomrule} \n \\tablecaption[")
	file.write(shortcaptions[i])
	file.write("]{")
	file.write(captions[i])
	file.write("}\label{")
	file.write(labels[i])
	file.write("}\n")
	if SCIP and SDP:
		file.write("\\begin{xtabular*}{\\textwidth}{@{\extracolsep{\\fill}}lrrrrrrrr")
	elif SCIP:
		file.write("\\begin{xtabular*}{\\textwidth}{@{\extracolsep{\\fill}}lrrrrr")
	else:
		file.write("\\begin{xtabular*}{\\textwidth}{@{\extracolsep{\\fill}}lrrr")
	file.write("@{}} \n ")
	for j in range(ninstances[i]):
		file.write(names[i][j].replace("_", "\\_").split(".")[0] + "& ")
		if SCIP:
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
		if SCIP:
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
			file.write("--")
		else:
			file.write("\\num{" + "%.1f" % float(times[i][j]) + "}")
		if SDP:
			if sdpiters[i][j] == "-":
				file.write("& -- & ")
			else:
				file.write("& \\num{" + "%.0f" % float(sdpiters[i][j]) + "} &")
			if penalties[i][j] == "-":
				file.write("-- & ")
			else:
				file.write("\\SI{" + "%.2f" % float(penalties[i][j]) + "}{\percent} &")
			if unsolved[i][j] == "-":
				file.write("-- ")
			else:
				file.write("\\SI{" + "%.2f" % float(unsolved[i][j]) + "}{\percent}")
		file.write("\\\ \n")
	file.write("  \\end{xtabular*} \n \\end{scriptsize} \n")

def makeOverviewBySubset(file, args, settings, s, f, caption, label, settingnames):
	if overviewlongtable:
		file.write("\\begin{longtable}[ht] {p{.5\\textwidth} | p{.15\\textwidth} p{.15\\textwidth} p{.15\\textwidth} ")
		file.write("p{.15\\textwidth}} \\toprule \n settings & solved & nodes & time ")
		file.write("\\\ \midrule \n")
	else:
		file.write("\\begin{table} \n \\begin{scriptsize} \\caption{" + caption + "} \n \\label{" + label + "} \n \\begin{tabular*}{\\textwidth}{@{}l@{\\;\\;\extracolsep{\\fill}}rrr")
		file.write("@{}}\\toprule \n")
		file.write(" settings & solved & time & nodes")
		file.write("\\\ \midrule \n")
	i = 0
	ind = 0
	argind = -1
	nyalmiperrors = 0
	for arg in args:
		if i in settings:
			nsolved = 0
			unfailed = 0
			aborts = 0
			timegeom = 1.0
			nodegeom = 1.0
			nnodescounted = 0
			j = s
			numinstances = f - s + 1
			if i == 0:
				nyalmiperrors = 0
			while j <= f:
				if nodes[i][j] != "-" and (not nodegeomonlysolved or allSolved(j, settings)):
				    nnodescounted = nnodescounted + 1
				j = j + 1
			j = s
			while j <= f:
				if primalresults[0][j] == "error":
					if i == 0:
						nyalmiperrors = nyalmiperrors + 1
					j = j + 1
					continue #do not count instances that YALMIP solved incorrectly
				if instanceSolved(i,j) :
					nsolved = nsolved + 1
				timegeom = math.pow(timegeom, float(j-s) / float(j-s+1)) * math.pow(float(times[i][j] + timeshift), 1.0/float(j-s+1))
				if nodes[i][j] != "-" and (not nodegeomonlysolved or allSolved(j, settings)):
					nodegeom = nodegeom * math.pow(float(nodes[i][j] + nodeshift), 1.0/float(nnodescounted))
				j = j + 1
			timegeom = timegeom - timeshift
			nodegeom = nodegeom - nodeshift
			file.write(settingnames[ind] + " & \\num{" + "%.0f" % nsolved + "} & \\num{" + "%.1f" % timegeom + "} & \\num{" + "%.1f" % nodegeom + "}")
			file.write(" \\\ \n")
			ind = ind + 1
		i = i + 1
	print("yalmip errors: " + str(nyalmiperrors))
	if overviewlongtable:
		file.write("\\bottomrule \n \\caption{" + caption + " except for " + str(nyalmiperrors) + " instances with wrong results from YALMIP} \n \\label{" + label + "} \n \\end{longtable} \n")
	else:
		file.write("\\bottomrule \n \\end{tabular*} \n \end{scriptsize} \n \\end{table} \n")

def makePerformancePlot(file, args, settings, s, f, caption, label, settingnames):
	performance=[[0 for x in range(f-s+1)] for x in range(len(sys.argv))] #initialize performance matrix
	j = s
	assert(len(colors) >= len(args))
	nyalmiperrors = 0
	while j <= f:
		if primalresults[0][j] == "error":
			nyalmiperrors = nyalmiperrors + 1
			for arg in settings:
				performance[arg][j-s] = 1e20
			j = j + 1
			continue #do not count instances that YALMIP solved incorrectly
		best = float('inf')
		# compute optimum
		for arg in settings:
			if times[arg][j] < best:
				best = times[arg][j]
		# divide all values by the optimum
		for arg in settings:
			if instanceSolved(arg,j):
                                if (best > 0.00):
				        performance[arg][j-s] = times[arg][j]# / best
                                else:
				        performance[arg][j-s] = 1e20
			else:
				performance[arg][j-s] = 1e20
		j = j + 1
	# sort performance for each arg
	for arg in settings:
		(performance[arg]).sort()
	#write the plot
	i = 0
	color = 0
	file.write("\\begin{figure} \n \\begin{scriptsize} \n \\centering \n \\begin{tikzpicture} \n \\begin{semilogxaxis}[xlabel={factor of time of fastest setting},xmin=0.92,xticklabels={0,1,10},ylabel={\\# solved instances}, legend style={cells={anchor=west},legend pos=outer north east,}, x label style={at={(axis description cs:0.5,0.025)}},y label style={at={(axis description cs:0.075,0.5)}}]\n")
	for arg in args:
		if i in settings:
			file.write("\\addplot[" + colors[color] + "]table[row sep=crcr]{x 	y \\\ \n")
			for j in range(f - s + 1):
				# if multiple instances are solved at the same (relative) time, only write the last one
				if (j == f -s or performance[i][j] < performance[i][j+1]) and performance[i][j] < 1e20:
					file.write(str(performance[i][j]) + " " + str(j + 1) + " \\\ \n")
			if usesettingnames:
				file.write("}; \n \\addlegendentry{" + settingnames[color] + "}")
			else:
				file.write("}; \n \\addlegendentry{$" + str((os.path.basename(arg)).split(".")[7]) + " | " + str((os.path.basename(arg)).split(".")[9]) + "$}")
			color = color + 1
		i = i + 1
	file.write("\\end{semilogxaxis} \n \\end{tikzpicture} \n")
	file.write("\\caption{" + caption + " except for " + str(nyalmiperrors) + " instances with wrong results in YALMIP} \n")
	file.write("\\label{" + label + "}\n")
	file.write("\\end{scriptsize} \n \\end{figure} \n")


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
	randfound=[[0 for x in range(maxinstances)] for x in range(len(sys.argv) -2)]   #initialize randfound matrix
	sdpiters=[[0 for x in range(maxinstances)] for x in range(len(sys.argv) -2)]   #initialize sdpiterations matrix
	penalties=[[0 for x in range(maxinstances)] for x in range(len(sys.argv) -2)]   #initialize penalty% matrix
	unsolved=[[0 for x in range(maxinstances)] for x in range(len(sys.argv) -2)]   #initialize unsolved% matrix
	sys.argv.remove("compareScriptSCIPSDP_YALMIP_Pajarito_SCIPSDPLP.py")				 #remove function call from list of files
	texfilename = sys.argv[0]
	file=open(texfilename, "w")
	if texfile:
		file.write("\\documentclass[landscape]{article} \n \\usepackage{amsmath} \n \\usepackage{amsthm} \n \\usepackage{dsfont} \n \\usepackage[dvipsnames]{xcolor} \n \\usepackage{booktabs} \n \\usepackage{multirow} \n \\usepackage{mathtools} \n \\usepackage{xtab} \n \usepackage{tikz} \n \usepackage{pgfplots} \n  \usepackage[margin=1in,footskip=0.25in]{geometry} \n \\extrafloats{100} \n \\usepackage{sistyle} \n \\usepackage{tabularx} \n \\newcommand{\\setting}[1]{\\texttt{#1}} \n \\begin{document} \n")
	sys.argv.remove(sys.argv[0])
	assert(len(subsetnames) >= len(subsets))
	#if len(sys.argv) != 5:
		#print("please call this using sixth arguments, a tex-file for output, 2 YALMIP-results-files, 1 Pajarito output-file, 2 SCIP-SDP output-files")
		#assert(0)
	readFile(sys.argv[0],0)
	readFile(sys.argv[1],1)
	readYALMIP(sys.argv[2], 2)
	readYALMIP(sys.argv[3], 3)
	readYALMIP(sys.argv[4], 4)
	readYALMIP(sys.argv[5], 5)
	readYALMIP(sys.argv[6], 6)
	readPajarito(sys.argv[7], 7)

	if completeTable:
		makeCompleteTableCaption(sys.argv[0], 0, 1, 1)
		makeCompleteTableCaption(sys.argv[1], 1, 1, 0)
		makeCompleteTableCaption(sys.argv[2], 2, 0, 0)
		makeCompleteTableCaption(sys.argv[3], 3, 0, 0)
		makeCompleteTableCaption(sys.argv[4], 4, 0, 0)
		makeCompleteTableCaption(sys.argv[5], 5, 0, 0)
		makeCompleteTableCaption(sys.argv[6], 6, 0, 0)
		makeCompleteTableCaption(sys.argv[7], 7, 0, 0)

	else:
		makeOverviewBySubset(file, sys.argv, [0,1,2,3,4], 0, 64, "Results for SCIP-SDP, YALMIP and Pajarito for the \emph{cardinality constrained least squares} testset of 65 instances", "CLSTable", ["\setting{SCIP-SDP (NL-BB)}", "\setting{SCIP-SDP (Cut-LP)}", "\setting{YALMIP (BNB)}", "\setting{YALMIP (CUTSDP)}", "\setting{Pajarito}"])
		makeOverviewBySubset(file, sys.argv, [0,1,2,3,4], 65, 133, "Results for SCIP-SDP, YALMIP and Pajarito for the \emph{min-$k$-partitioning} testset of 69 instances", "MinKTable", ["\setting{SCIP-SDP (NL-BB)}", "\setting{SCIP-SDP (Cut-LP)}", "\setting{YALMIP (BNB)}", "\setting{YALMIP (CUTSDP)}", "\setting{Pajarito}"])
		makeOverviewBySubset(file, sys.argv, [0,1,2,3,4], 134, 193, "Results for SCIP-SDP, YALMIP and Pajarito for the \emph{truss topology} testset of 60 instances", "TTTable", ["\setting{SCIP-SDP (NL-BB)}", "\setting{SCIP-SDP (Cut-LP)}", "\setting{YALMIP (BNB)}", "\setting{YALMIP (CUTSDP)}", "\setting{Pajarito}"])
		makeOverviewBySubset(file, sys.argv, [0,1,2,3,4], 0, 193, "Results for SCIP-SDP, YALMIP and Pajarito for the \emph{complete testset} of 194 instances", "OverallTable", ["\setting{SCIP-SDP (NL-BB)}", "\setting{SCIP-SDP (Cut-LP)}", "\setting{YALMIP (BNB)}", "\setting{YALMIP (CUTSDP)}", "\setting{Pajarito}"])
		makePerformancePlot(file, sys.argv, [0,1,2,3,4], 0, 64, "Results for SCIP-SDP, YALMIP and Pajarito for the \emph{cardinality constrained least squares} testset of 65 instances", "CLSPlot", ["\setting{SCIP-SDP (NL-BB)}", "\setting{SCIP-SDP (Cut-LP)}", "\setting{YALMIP (BNB)}", "\setting{YALMIP (CUTSDP)}", "\setting{Pajarito}"])
		makePerformancePlot(file, sys.argv, [0,1,2,3,4], 65, 133, "Results for SCIP-SDP, YALMIP and Pajarito for the \emph{min-$k$-partitioning} testset of 69 instances", "MinKPlot", ["\setting{SCIP-SDP (NL-BB)}", "\setting{SCIP-SDP (Cut-LP)}", "\setting{YALMIP (BNB)}", "\setting{YALMIP (CUTSDP)}", "\setting{Pajarito}"])
		makePerformancePlot(file, sys.argv, [0,1,2,3,4], 134, 193, "Results for SCIP-SDP, YALMIP and Pajarito for the \emph{truss topology} testset of 60 instances", "TTPlot", ["\setting{SCIP-SDP (NL-BB)}", "\setting{SCIP-SDP (Cut-LP)}", "\setting{YALMIP (BNB)}", "\setting{YALMIP (CUTSDP)}", "\setting{Pajarito}"])
		makePerformancePlot(file, sys.argv, [0,1,2,3,4], 0, 193, "Results for SCIP-SDP, YALMIP and Pajarito for the \emph{complete testset} of 194 instances", "OverallPlot", ["\setting{SCIP-SDP (NL-BB)}", "\setting{SCIP-SDP (Cut-LP)}", "\setting{YALMIP (BNB)}", "\setting{YALMIP (CUTSDP)}", "\setting{Pajarito}"])

	if texfile:
		file.write("\\end{document}")
	file.close()
