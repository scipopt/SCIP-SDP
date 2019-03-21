import os
import sys
import math
from decimal import Decimal
from numpy import inf

#TODO: write tables for online supplement

#which tables/plots for the paper should be prepared?
texfile = 0 # make the output a compilable texfile or just a figure/table that can be included?
Lbounds = 0
Rbounds = 0
LRbounds = 0
Ltimes = 0
Rtimes = 0
LRtimes = 0
printFails = 0
completeTable = 0
CompleteRICtable = 0
completeResultsTable=1

MISDPfilename = "/workopt/gally/results/RIP-results/check.RIPMISDP.scipsdp.linux.x86_64.gnu.opt.sdpa.extra.branchinfobj_nofracdive.out"
Asp0708path = "/workopt/gally/results/SDPA-results/"
Asp07filenames = ["0+-115305A_07l.out", "0+-115305A_07r.out", "0+-115305B_07l.out", "0+-115305B_07r.out", "0+-115305C_07l.out", "0+-115305C_07r.out", "0+-125354A_07l.out", "0+-125354A_07r.out", "0+-125354B_07l.out", "0+-125354B_07r.out", "0+-125354C_07l.out", "0+-125354C_07r.out", "0+-130403A_07l.out", "0+-130403A_07r.out", "0+-130403B_07l.out", "0+-130403B_07r.out", "0+-130403C_07l.out", "0+-130403C_07r.out", "bina15305A_07l.out", "bina15305A_07r.out", "bina15305B_07l.out", "bina15305B_07r.out", "bina15305C_07l.out", "bina15305C_07r.out", "bina25354A_07l.out", "bina25354A_07r.out", "bina25354B_07l.out", "bina25354B_07r.out", "bina25354C_07l.out", "bina25354C_07r.out", "bina30403A_07l.out", "bina30403A_07r.out", "bina30403B_07l.out", "bina30403B_07r.out", "bina30403C_07l.out", "bina30403C_07r.out", "bern15305A_07l.out", "bern15305A_07r.out", "bern15305B_07l.out", "bern15305B_07r.out", "bern15305C_07l.out", "bern15305C_07r.out", "bern25354A_07l.out", "bern25354A_07r.out", "bern25354B_07l.out", "bern25354B_07r.out", "bern25354C_07l.out", "bern25354C_07r.out", "bern30403A_07l.out", "bern30403A_07r.out", "bern30403B_07l.out", "bern30403B_07r.out", "bern30403C_07l.out", "bern30403C_07r.out", "norm15305A_07l.out", "norm15305A_07r.out", "norm15305B_07l.out", "norm15305B_07r.out", "norm15305C_07l.out", "norm15305C_07r.out", "norm25354A_07l.out", "norm25354A_07r.out", "norm25354B_07l.out", "norm25354B_07r.out", "norm25354C_07l.out", "norm25354C_07r.out", "norm30403A_07l.out", "norm30403A_07r.out", "norm30403B_07l.out", "norm30403B_07r.out", "norm30403C_07l.out", "norm30403C_07r.out", "wish15305A_07l.out", "wish15305A_07r.out", "wish15305B_07l.out", "wish15305B_07r.out", "wish15305C_07l.out", "wish15305C_07r.out", "wish25354A_07l.out", "wish25354A_07r.out", "wish25354B_07l.out", "wish25354B_07r.out", "wish25354C_07l.out", "wish25354C_07r.out", "wish30403A_07l.out", "wish30403A_07r.out", "wish30403B_07l.out", "wish30403B_07r.out", "wish30403C_07l.out", "wish30403C_07r.out", "band30305A_07l.out", "band30305A_07r.out", "band30305B_07l.out", "band30305B_07r.out", "band30305C_07l.out", "band30305C_07r.out", "band35354A_07l.out", "band35354A_07r.out", "band35354B_07l.out", "band35354B_07r.out", "band35354C_07l.out", "band35354C_07r.out", "band40403A_07l.out", "band40403A_07r.out", "band40403B_07l.out", "band40403B_07r.out", "band40403C_07l.out", "band40403C_07r.out", "rnk130305A_07l.out", "rnk130305A_07r.out", "rnk130305B_07l.out", "rnk130305B_07r.out", "rnk130305C_07l.out", "rnk130305C_07r.out", "rnk135354A_07l.out", "rnk135354A_07r.out", "rnk135354B_07l.out", "rnk135354B_07r.out", "rnk135354C_07l.out", "rnk135354C_07r.out", "rnk140403A_07l.out", "rnk140403A_07r.out", "rnk140403B_07l.out", "rnk140403B_07r.out", "rnk140403C_07l.out", "rnk140403C_07r.out"]
Asp08basenames = ["0+-115305A_08l", "0+-115305A_08r", "0+-115305B_08l", "0+-115305B_08r", "0+-115305C_08l", "0+-115305C_08r", "0+-125354A_08l", "0+-125354A_08r", "0+-125354B_08l", "0+-125354B_08r", "0+-125354C_08l", "0+-125354C_08r", "0+-130403A_08l", "0+-130403A_08r", "0+-130403B_08l", "0+-130403B_08r", "0+-130403C_08l", "0+-130403C_08r", "bina15305A_08l", "bina15305A_08r", "bina15305B_08l", "bina15305B_08r", "bina15305C_08l", "bina15305C_08r", "bina25354A_08l", "bina25354A_08r", "bina25354B_08l", "bina25354B_08r", "bina25354C_08l", "bina25354C_08r", "bina30403A_08l", "bina30403A_08r", "bina30403B_08l", "bina30403B_08r", "bina30403C_08l", "bina30403C_08r", "bern15305A_08l", "bern15305A_08r", "bern15305B_08l", "bern15305B_08r", "bern15305C_08l", "bern15305C_08r", "bern25354A_08l", "bern25354A_08r", "bern25354B_08l", "bern25354B_08r", "bern25354C_08l", "bern25354C_08r", "bern30403A_08l", "bern30403A_08r", "bern30403B_08l", "bern30403B_08r", "bern30403C_08l", "bern30403C_08r", "norm15305A_08l", "norm15305A_08r", "norm15305B_08l", "norm15305B_08r", "norm15305C_08l", "norm15305C_08r", "norm25354A_08l", "norm25354A_08r", "norm25354B_08l", "norm25354B_08r", "norm25354C_08l", "norm25354C_08r", "norm30403A_08l", "norm30403A_08r", "norm30403B_08l", "norm30403B_08r", "norm30403C_08l", "norm30403C_08r", "wish15305A_08l", "wish15305A_08r", "wish15305B_08l", "wish15305B_08r", "wish15305C_08l", "wish15305C_08r", "wish25354A_08l", "wish25354A_08r", "wish25354B_08l", "wish25354B_08r", "wish25354C_08l", "wish25354C_08r", "wish30403A_08l", "wish30403A_08r", "wish30403B_08l", "wish30403B_08r", "wish30403C_08l", "wish30403C_08r", "band30305A_08l", "band30305A_08r", "band30305B_08l", "band30305B_08r", "band30305C_08l", "band30305C_08r", "band35354A_08l", "band35354A_08r", "band35354B_08l", "band35354B_08r", "band35354C_08l", "band35354C_08r", "band40403A_08l", "band40403A_08r", "band40403B_08l", "band40403B_08r", "band40403C_08l", "band40403C_08r", "rnk130305A_08l", "rnk130305A_08r", "rnk130305B_08l", "rnk130305B_08r", "rnk130305C_08l", "rnk130305C_08r", "rnk135354A_08l", "rnk135354A_08r", "rnk135354B_08l", "rnk135354B_08r", "rnk135354C_08l", "rnk135354C_08r", "rnk140403A_08l", "rnk140403A_08r", "rnk140403B_08l", "rnk140403B_08r", "rnk140403C_08l", "rnk140403C_08r"]
Asp08Datapath = "/workopt/gally/results/SDPA-results/"
nAsp08Steps = 15
epsilon = 0.0001
timeshift = 10



#argument should be the tex-file to write to

def readFileMISDP(arg):
	j=0;
	file=open(arg, "r")
	filestring = file.read()
	while filestring.find("@01") > -1:
		assert( j < 126)
		substring = filestring[filestring.find("@01"):filestring.find("@04")]
		assert(substring.find("@01 ") > -1)
		names[0][j] = os.path.basename(substring[substring.find("@01 ") + 3:].split()[0])
		if substring.find("Solving Time (sec) : ") == -1:
			times[0][j] = 14400
		else:
			time = substring[substring.find("Solving Time (sec) : ") + 21:].split()[0]
			times[0][j] = float(time)
			if times[0][j] > 14400.0:
				times[0][j] = 14400
		if substring.find("Solving Nodes      : ") == -1:
			nodes[0][j] = "-"
		else:
			node = substring[substring.find("Solving Nodes      : ") + 21:].split()[0]
			nodes[0][j] = int(node)
		if substring.find("Dual Bound         : ") == -1:
			dualresults[0][j] = -1e20
		else:
			dualresult = substring[substring.find("Dual Bound         : ") + 21:].split()[0]
			dualresults[0][j] = float(dualresult)
		if substring.find("Primal Bound       : ") == -1:
			primalresults[0][j] = 1e20
		else:
			primalresult = substring[substring.find("Primal Bound       : ") + 21:].split()[0]
			primalresults[0][j] = float(primalresult)
		if substring.find("Gap                : ") == -1:
			gaps[0][j] = "infinite"
		else:
			gap = substring[substring.find("Gap                : ") + 21:].split()[0]
			gaps[0][j] = gap
		if substring.find("SDP iterations:") == -1:
			sdpiters[0][j] = "-"
		else:
			sdpiter = substring[substring.find("SDP iterations:") + 15:].split()[0]
			sdpiters[0][j] = int(sdpiter)
		if j % 2 == 1:
			primalresults[0][j] *= -1
			dualresults[0][j] *= -1
		if substring.find("Percentage penalty formulation used:	 ") == -1:
			penalties[0][j] = "-"
		else:
			penalty = substring[substring.find("Percentage penalty formulation used:") + 38:].split()[0]
			penalties[0][j] = float(penalty)
		if substring.find("Percentage unsolved even with penalty:	") == -1:
			unsolved[0][j] = "-"
		else:
			unsolv = substring[substring.find("Percentage unsolved even with penalty:") + 40:].split()[0]
			unsolved[0][j] = float(unsolv)
		filestring = filestring[filestring.find("@04")+3:]
		j = j+1
	assert(j == 126)
	file.close()

def readFileAsp07(filename, j, rhs):
	file=open(Asp0708path + filename, "r")
	filestring = file.read()
	substring = filestring[:filestring.find("** Parameters **")]
	time = substring[substring.find("total time   = ") + 15:].split()[0]
	times[1][j] = float(time)
	primal = substring[substring.find("objValPrimal = ") + 15:].split()[0]
	primalresults[1][j] = float(primal)
	dual = substring[substring.find("objValDual   = ") + 15:].split()[0]
	dualresults[1][j] = float(dual)
	gap = substring[substring.find("relative gap = ") + 15:].split()[0]
	gaps[1][j] = float(gap)
	iters = substring[substring.find("   Iteration = ") + 15:].split()[0]
	sdpiters[1][j] = float(iters)
	returncode = substring[substring.find("phase.value  = ") + 15:].split()[0]
	if substring.find("phase.value  = ") == -1:
		Asp07unsolved[j] = -1
	else:
		returncode = substring[substring.find("phase.value  = ") + 15:].split()[0]
		if returncode != "pdOPT":
			Asp07unsolved[j] = 1
	if rhs:
		primalresults[1][j] *= -1
		dualresults[1][j] *= -1
	file.close()

def readFileAsp08(basename, datafilename, j, rhs):
	print("reading Asp08 instance " + str(j))
	for i in range(nAsp08Steps):
		filename = basename + "_" + str(i+1) + ".out"
		file=open(Asp0708path + filename , "r")
		filestring = file.read()
		substring = filestring[:filestring.find("** Parameters **")]
		if substring.find("total time   = ") == -1:
			Asp08times[j][i] = 14400
		else:
			time = substring[substring.find("total time   = ") + 15:].split()[0]
			Asp08times[j][i] = float(time)
		if substring.find("objValPrimal = ") == -1:
			Asp08primalresults[j][i] = "-"
		else:
			primal = substring[substring.find("objValPrimal = ") + 15:].split()[0]
			Asp08primalresults[j][i] = float(primal)
		if substring.find("objValDual   = ") == -1:
			Asp08dualresults[j][i] = "-"
		else:
			dual = substring[substring.find("objValDual   = ") + 15:].split()[0]
			Asp08dualresults[j][i] = float(dual)
		if substring.find("relative gap = ") == -1:
			Asp08gaps[j][i] = "-"
		else:
			gap = substring[substring.find("relative gap = ") + 15:].split()[0]
			Asp08gaps[j][i] = float(gap)
		if substring.find("   Iteration = ") == -1:
			Asp08sdpiters[j][i] = "-"
		else:
			iters = substring[substring.find("   Iteration = ") + 15:].split()[0]
			Asp08sdpiters[j][i] = float(iters)
		if substring.find("phase.value  = ") == -1:
			Asp08unsolved[j][i] = -1
		else:
			returncode = substring[substring.find("phase.value  = ") + 15:].split()[0]
			if returncode != "pdOPT":
				Asp08unsolved[j][i] = 1	#TODO maybe also reset results, iters etc. or check for gap below
		file.close()
	# compute best bound
	datafilename = basename + "_1" + ".dat-s"
	file=open(Asp08Datapath + datafilename , "r")
	fileString = file.read()
	kString = fileString[fileString.find("order k=") + 8:].split()[0]
	k = int(kString)
	if rhs:
		dualresults[2][j] = float("inf")
	else:
		dualresults[2][j] = float("-inf")
		alphaString = fileString[fileString.find("alpha =") + 7:].split()[0]
		alpha = float(alphaString)
	file.close()
	times[2][j] = 0
	for i in range(nAsp08Steps):
		times[2][j] = times[2][j] + Asp08times[j][i]
		if Asp08dualresults[j][i] == "-":
			continue
		datafilename = basename + "_" + str(i+1) + ".dat-s"
		file=open(Asp08Datapath + datafilename , "r")
		fileString = file.read()
		rhoString = fileString[fileString.find("rho =") + 5:].split()[0]
		rho = float(rhoString)
		file.close()
		if rhs:
			result = Asp08dualresults[j][i] + rho * k
			if result < dualresults[2][j]:
				dualresults[2][j] = result
				sdpiters[2][j] = Asp08sdpiters[j][i]
		else:
			result = alpha - Asp08dualresults[j][i] - rho * k
			if result > dualresults[2][j]:
				dualresults[2][j] = result
				sdpiters[2][j] = Asp08sdpiters[j][i]


def CompleteResultsTable():
	#RIP-MISDP
	file.write("\\newpage \n \\begin{scriptsize} \n  \\setlength{\\tabcolsep}{2pt} \n \\tablehead{\\toprule \n")
	file.write("problem & dbound &  pbound & gap & nodes & time & iters & pen & uns ")
	dive = 0
	rand = 0
	fix = 0
	if dive:
		file.write("& dive ")
	if rand:
		file.write("& rand ")
	if fix:
		file.write("& fix ")
	file.write("\\\ \\midrule} \n \\tabletail{ \\midrule \\multicolumn{3}{@{}l}{continued on next page \\dots}\\\ \\bottomrule} \n \\tablelasttail{\\bottomrule} \n \\tablecaption[")
	file.write("RIP-MISDP")
	file.write("]{")
	file.write("Full results for MISDP formulation of RIP with \sdpa~7.3.8 on a Linux cluster with Intel i3 CPUs with \SI{3.20}{GHz} and \SI{8}{GB} memory")
	file.write("}\label{")
	file.write("RIPMISDP")
	file.write("}\n")
	file.write("\\begin{xtabular*}{\\textwidth}{@{\extracolsep{\\fill}}lrrrrrrrr")
	if dive:
		file.write("r")
	if rand:
		file.write("r")
	if fix:
		file.write("r")
	file.write("@{}} \n ")
	for j in range(126):
		file.write(names[0][j].replace("_", "\\_").split(".")[0] + "& ")
		if dualresults[0][j] > -1e20:
			file.write("\\num{" + "%.2f" % float(dualresults[0][j]) + "} ")
		else:
			file.write(" $-\infty$ ")
		file.write("& ")
		if primalresults[0][j] < 1e20:
			file.write("\\num{" + "%.2f" % float(primalresults[0][j]) + "}")
		else:
			file.write(" $\infty$ ")
		file.write("& ")
		if gaps[0][j] != "infinite":
			if gaps[0][j] != "-":
				file.write("\\num{" + "%.2f" % float(gaps[0][j]) + "}\,\% & ")
			else:
				file.write("-- & ")
		else:
			file.write("$\infty$ & ")
		if nodes[0][j] == "-":
			file.write("-- & ")
		else:
			file.write("\\num{" + "%.0f" % float(nodes[0][j]) + "} &")
		if times[0][j] == "-":
			file.write("-- & ")
		else:
			file.write("\\num{" + "%.1f" % float(times[0][j]) + "} &")
		if sdpiters[0][j] == "-":
			file.write("-- & ")
		else:
			file.write("\\num{" + "%.0f" % float(sdpiters[0][j]) + "} &")
		if penalties[0][j] == "-":
			file.write("-- & ")
		else:
			file.write("\\SI{" + "%.2f" % float(penalties[0][j]) + "}{\percent} &")
		if unsolved[0][j] == "-":
			file.write("-- ")
		else:
			file.write("\\SI{" + "%.2f" % float(unsolved[0][j]) + "}{\percent}")
		if dive:
			if fracdivefound[0][j] == "-":
				file.write("& -- ")
			else:
				file.write("& \\num{" + "%.0f" % float(fracdivefound[0][j]) + "}")
		if rand:
			if randfound[0][j] == "-":
				file.write("& -- ")
			else:
				file.write("& \\num{" + "%.0f" % float(randfound[0][j]) + "}")
		if fix:
			if sdpredcostfixings[0][j] == "-":
				file.write("& -- ")
			else:
				file.write("& \\num{" + "%.0f" % float(sdpredcostfixings[0][j]) + "} ")
		file.write("\\\ \n")
	file.write("  \\end{xtabular*} \n \\end{scriptsize} \n")
	#A07
	file.write("\\newpage \n \\begin{scriptsize} \n  \\setlength{\\tabcolsep}{2pt} \n \\tablehead{\\toprule \n")
	file.write("problem & dbound &  pbound & gap & time & iters  ")
	dive = 0
	rand = 0
	fix = 0
	if dive:
		file.write("& dive ")
	if rand:
		file.write("& rand ")
	if fix:
		file.write("& fix ")
	file.write("\\\ \\midrule} \n \\tabletail{ \\midrule \\multicolumn{3}{@{}l}{continued on next page \\dots}\\\ \\bottomrule} \n \\tablelasttail{\\bottomrule} \n \\tablecaption[")
	file.write("RIP-(A07)")
	file.write("]{")
	file.write("Full results for (A07) with \sdpa~7.3.8 on a Linux cluster with Intel i3 CPUs with \SI{3.20}{GHz} and \SI{8}{GB} memory")
	file.write("}\label{")
	file.write("RIPA07")
	file.write("}\n")
	file.write("\\begin{xtabular*}{\\textwidth}{@{\extracolsep{\\fill}}lrrrrr")
	if dive:
		file.write("r")
	if rand:
		file.write("r")
	if fix:
		file.write("r")
	file.write("@{}} \n ")
	for j in range(126):
		file.write(names[0][j].replace("_", "\\_").split(".")[0] + "& ")
		if dualresults[1][j] > -1e20:
			file.write("\\num{" + "%.2f" % float(dualresults[1][j]) + "} ")
		else:
			file.write(" $-\infty$ ")
		file.write("& ")
		if primalresults[1][j] < 1e20:
			file.write("\\num{" + "%.2f" % float(primalresults[1][j]) + "}")
		else:
			file.write(" $\infty$ ")
		file.write("& ")
		if gaps[1][j] != "infinite":
			if gaps[1][j] != "-":
				file.write("\\num{" + "%.2f" % float(gaps[1][j]) + "}\,\% & ")
			else:
				file.write("-- & ")
		else:
			file.write("$\infty$ & ")
		if times[1][j] == "-":
			file.write("-- & ")
		else:
			file.write("\\num{" + "%.1f" % float(times[1][j]) + "} &")
		if sdpiters[1][j] == "-":
			file.write("-- ")
		else:
			file.write("\\num{" + "%.0f" % float(sdpiters[1][j]) + "}")
		file.write("\\\ \n")
	file.write("  \\end{xtabular*} \n \\end{scriptsize} \n")
	#A08
	file.write("\\newpage \n \\begin{scriptsize} \n  \\setlength{\\tabcolsep}{2pt} \n \\tablehead{\\toprule \n")
	file.write("problem & dbound &  pbound & gap & time & iters  ")
	dive = 0
	rand = 0
	fix = 0
	if dive:
		file.write("& dive ")
	if rand:
		file.write("& rand ")
	if fix:
		file.write("& fix ")
	file.write("\\\ \\midrule} \n \\tabletail{ \\midrule \\multicolumn{3}{@{}l}{continued on next page \\dots}\\\ \\bottomrule} \n \\tablelasttail{\\bottomrule} \n \\tablecaption[")
	file.write("RIP-(A08)")
	file.write("]{")
	file.write("Full results for (A08D) subproblems with \sdpa~7.3.8 on a Linux cluster with Intel i3 CPUs with \SI{3.20}{GHz} and \SI{8}{GB} memory")
	file.write("}\label{")
	file.write("RIPA08D")
	file.write("}\n")
	file.write("\\begin{xtabular*}{\\textwidth}{@{\extracolsep{\\fill}}lrrrrr")
	if dive:
		file.write("r")
	if rand:
		file.write("r")
	if fix:
		file.write("r")
	file.write("@{}} \n ")
	for j in range(126):
		for i in range(nAsp08Steps):
			file.write((names[0][j] + "_" + str(i+1)).replace("_", "\\_").split(".")[0] + "& ")
			if Asp08dualresults[j][i] == "-":
				file.write("-- ")
			elif Asp08dualresults[j][i] > -1e20:
				file.write("\\num{" + "%.2f" % float(Asp08dualresults[j][i]) + "} ")
			else:
				file.write(" $-\infty$ ")
			file.write("& ")
			if Asp08primalresults[j][i] == "-":
				file.write("-- ")
			elif Asp08primalresults[j][i] < 1e20:
				file.write("\\num{" + "%.2f" % float(Asp08primalresults[j][i]) + "}")
			else:
				file.write(" $\infty$ ")
			file.write("& ")
			if Asp08gaps[j][i] != "infinite":
				if Asp08gaps[j][i] != "-":
					file.write("\\num{" + "%.2f" % float(Asp08gaps[j][i]) + "}\,\% & ")
				else:
					file.write("-- & ")
			else:
				file.write("$\infty$ & ")
			if Asp08times[j][i] == "-":
				file.write("-- & ")
			else:
				file.write("\\num{" + "%.1f" % float(Asp08times[j][i]) + "} &")
			if Asp08sdpiters[j][i] == "-":
				file.write("-- ")
			else:
				file.write("\\num{" + "%.0f" % float(Asp08sdpiters[j][i]) + "}")
			file.write("\\\ \n")
	file.write("  \\end{xtabular*} \n \\end{scriptsize} \n")
	#A08-total
	file.write("\\newpage \n \\begin{scriptsize} \n  \\setlength{\\tabcolsep}{2pt} \n \\tablehead{\\toprule \n")
	file.write("problem & bound & time & iters  ")
	dive = 0
	rand = 0
	fix = 0
	if dive:
		file.write("& dive ")
	if rand:
		file.write("& rand ")
	if fix:
		file.write("& fix ")
	file.write("\\\ \\midrule} \n \\tabletail{ \\midrule \\multicolumn{3}{@{}l}{continued on next page \\dots}\\\ \\bottomrule} \n \\tablelasttail{\\bottomrule} \n \\tablecaption[")
	file.write("RIP-(A08)")
	file.write("]{")
	file.write("Full results for (A08) (as sum over all subproblems) with \sdpa~7.3.8 on a Linux cluster with Intel i3 CPUs with \SI{3.20}{GHz} and \SI{8}{GB} memory")
	file.write("}\label{")
	file.write("RIPA07")
	file.write("}\n")
	file.write("\\begin{xtabular*}{\\textwidth}{@{\extracolsep{\\fill}}lrrr")
	if dive:
		file.write("r")
	if rand:
		file.write("r")
	if fix:
		file.write("r")
	file.write("@{}} \n ")
	for j in range(126):
		file.write(names[0][j].replace("_", "\\_").split(".")[0] + "& ")
		if float(dualresults[2][j]) == inf:
			file.write(" $\infty$ ")
		elif dualresults[2][j] > -1e20:
			file.write("\\num{" + "%.2f" % float(dualresults[2][j]) + "} ")
		else:
			file.write(" $-\infty$ ")
		file.write("& ")
		if times[2][j] == "-":
			file.write("-- & ")
		else:
			file.write("\\num{" + "%.1f" % float(times[2][j]) + "} &")
		if sdpiters[2][j] == "-":
			file.write("-- ")
		else:
			file.write("\\num{" + "%.0f" % float(sdpiters[2][j]) + "}")
		file.write("\\\ \n")
	file.write("  \\end{xtabular*} \n \\end{scriptsize} \n")



def makeCompleteTableCaptionMISDP(shortcaption, caption, label, file, i):
	file.write("\\newpage \n \\begin{scriptsize} \n  \\setlength{\\tabcolsep}{2pt} \n \\tablehead{\\toprule \n")
	file.write("problem & dbound &  pbound & gap & nodes & time & iters & pen & uns & dive & rand & fix ")
	file.write("\\\ \\midrule} \n \\tabletail{ \\midrule \\multicolumn{3}{@{}l}{continued on next page \\dots}\\\ \\bottomrule} \n \\tablelasttail{\\bottomrule} \n \\tablecaption[")
	file.write(shortcaption)
	file.write("]{")
	file.write(caption)
	file.write("}\label{")
	file.write(label)
	file.write("}\n")
	file.write("\\begin{xtabular*}{0.48\\textwidth}{@{\extracolsep{\\fill}}lrrrrrr@{}} \n ")
	for j in range(ninstances[i]):
		file.write(names[i][j].replace("_", "\\_").split(".")[0] + "& ")
		if dualresults[i][j] > -1e20:
			file.write("\\num{" + "%.2f" % dualresults[i][j] + "} ")
		else:
			file.write(" $-\infty$ ")
		file.write("& ")
		if primalresults[i][j] < 1e20:
			file.write("\\num{" + "%.2f" % primalresults[i][j] + "}")
		else:
			file.write(" $\infty$ ")
		file.write("& $")
		if gaps[i][j] != "infinite":
			file.write(str(gaps[i][j]))
		else:
			file.write("\infty")
		file.write("$ & $" + str(nodes[i][j]) + "$ & $" + str(times[i][j]) + "$ & $" + str(sdpiters[i][j])  + "$ \\\ \n")
	file.write("  \\end{xtabular*} \n \\end{scriptsize} \n")



def LhsResultsTable(instancesets, instancesetnames, caption, label):
	file.write("\\begin{table} \n \\begin{scriptsize} \\caption{" + caption + "} \n \\label{" + label + "} \n \\begin{tabular*}{\\linewidth}{@{}l@{\\;\\;\extracolsep{\\fill}}rr")
	file.write("@{}}\\toprule \n")
	file.write(" matrices & \eqref{Asp07} & \eqref{Asp08} ")
	file.write("\\\ \midrule \n")
	totalDiff07 = 0
	totalnInstances07 = 0
	totalDiff08 = 0
	totalnInstances08 = 0
	i = 0
	for instancetype in instancesets:
		diff07 = 0
		nInstances07 = 0
		diff08 = 0
		nInstances08 = 0
		for j in instancetype:
			# we only take instances with alpha_k > 0, since otherwise the trivial bound is optimal and we cannot compute a relative gap
			if dualresults[1][j] != -1e20 and gaps[0][j] != "infinite" and float(gaps[0][j]) < epsilon and float(dualresults[0][j]) > epsilon:
				diff07 += (dualresults[0][j] - dualresults[1][j]) / dualresults[0][j]
				totalDiff07 += (float(dualresults[0][j]) - float(dualresults[1][j])) / float(dualresults[0][j])
				nInstances07 += 1
				totalnInstances07 += 1
			if dualresults[2][j] != -1e20 and gaps[0][j] != "infinite" and float(gaps[0][j]) < epsilon and float(dualresults[0][j]) > epsilon and dualresults[2][j] != -float("inf"):
				diff08 += (float(dualresults[0][j]) - float(dualresults[2][j])) / float(dualresults[0][j])
				nInstances08 += 1
				totalDiff08 += (float(dualresults[0][j]) - float(dualresults[2][j])) / float(dualresults[0][j])
				totalnInstances08 += 1
		if nInstances07 > 0:
			avgdiff07 = 100 * diff07 / nInstances07
		else:
			avgdiff07 = "-"
		if nInstances08 > 0:
			avgdiff08 = 100 * diff08 / nInstances08
		else:
			avgdiff08 = "-"
		file.write(instancesetnames[i])
		if avgdiff07 != "-":
			file.write(" & \\num{" + "%.2f" % avgdiff07 + "}\,\% &")
		else:
			file.write("& - &")
		if avgdiff08 != "-":
			file.write(" \\num{" + "%.2f" % avgdiff08 + "}\,\% \\\ \n")
		else:
			file.write("- \\\ \n")
		i += 1
	file.write("\\midrule \n")
	if totalnInstances07 > 0:
		totalavgdiff07 = 100 * totalDiff07 / totalnInstances07
	else:
		totalavgdiff07 = "-"
	if nInstances08 > 0:
		totalavgdiff08 = 100 * totalDiff08 / totalnInstances08
	else:
		totalavgdiff08 = "-"
	file.write("total & ")
	if avgdiff07 != "-":
		file.write(" \\num{" + "%.2f" % totalavgdiff07 + "}\,\% &")
	else:
		file.write("- &")
	if avgdiff08 != "-":
		file.write(" \\num{" + "%.2f" % totalavgdiff08 + "}\,\%\\\ \n")
	else:
		file.write("- \\\ \n")
	file.write("\\bottomrule \n \\end{tabular*} \n \end{scriptsize} \n \\end{table} \n")


def RhsResultsTable(instancesets, instancesetnames, caption, label):
	file.write("\\begin{table} \n \\begin{scriptsize} \\caption{" + caption + "} \n \\label{" + label + "} \n \\begin{tabular*}{\\linewidth}{@{}l@{\\;\\;\extracolsep{\\fill}}rr")
	file.write("@{}}\\toprule \n")
	file.write(" matrices & \eqref{Asp07} & \eqref{Asp08} ")
	file.write("\\\ \midrule \n")
	totalDiff07 = 0
	totalnInstances07 = 0
	totalDiff08 = 0
	totalnInstances08 = 0
	i = 0
	for instancetype in instancesets:
		diff07 = 0
		nInstances07 = 0
		diff08 = 0
		nInstances08 = 0
		for j in instancetype:
			# we only take instances with alpha_k > 0, since otherwise the trivial bound is optimal and we cannot compute a relative gap
			if dualresults[1][j] != -1e20 and gaps[0][j] != "infinite" and float(gaps[0][j]) < epsilon:
				diff07 += (float(dualresults[1][j]) - float(dualresults[0][j])) / float(dualresults[0][j])
				totalDiff07 += (float(dualresults[1][j]) - float(dualresults[0][j])) / float(dualresults[0][j])
				nInstances07 += 1
				totalnInstances07 += 1
			if dualresults[2][j] != -1e20 and gaps[0][j] != "infinite" and float(gaps[0][j]) < epsilon and dualresults[2][j] != float("inf"):
				diff08 += (float(dualresults[2][j]) - float(dualresults[0][j])) / float(dualresults[0][j])
				nInstances08 += 1
				totalDiff08 += (float(dualresults[2][j]) - float(dualresults[0][j])) / float(dualresults[0][j])
				totalnInstances08 += 1
		if nInstances07 > 0:
			avgdiff07 = 100 * diff07 / nInstances07
		else:
			avgdiff07 = "-"
		if nInstances08 > 0:
			avgdiff08 = 100 * diff08 / nInstances08
		else:
			avgdiff08 = "-"
		file.write(instancesetnames[i])
		if avgdiff07 != "-":
			file.write(" & \\num{" + "%.2f" % avgdiff07 + "}\,\% &")
		else:
			file.write("& - &")
		if avgdiff08 != "-":
			file.write(" \\num{" + "%.2f" % avgdiff08 + "}\,\% \\\ \n")
		else:
			file.write("- \\\ \n")
		i += 1
	file.write("\\midrule \n")
	if totalnInstances07 > 0:
		totalavgdiff07 = 100 * totalDiff07 / totalnInstances07
	else:
		totalavgdiff07 = "-"
	if nInstances08 > 0:
		totalavgdiff08 = 100 * totalDiff08 / totalnInstances08
	else:
		totalavgdiff08 = "-"
	file.write("total & ")
	if avgdiff07 != "-":
		file.write(" \\num{" + "%.2f" % totalavgdiff07 + "}\,\% &")
	else:
		file.write("- &")
	if avgdiff08 != "-":
		file.write(" \\num{" + "%.2f" % totalavgdiff08 + "}\,\%\\\ \n")
	else:
		file.write("- \\\ \n")
	file.write("\\bottomrule \n \\end{tabular*} \n \end{scriptsize} \n \\end{table} \n")


def LhsRhsResultsTable(instancesets, instancesetnames, caption, label):
	file.write("\\begin{table} \n \\begin{scriptsize} \\caption{" + caption + "} \n \\label{" + label + "} \n \\begin{tabular*}{\\linewidth}{@{}l@{\\;\\;\extracolsep{\\fill}}rrrr")
	file.write("@{}}\\toprule \n")
	file.write("  & \\multicolumn{2}{c}{$\\alpha_k$} & \\multicolumn{2}{c}{$\\beta_k$} \\\ \n")
	file.write("\\cmidrule(r){2-3} \\cmidrule(l){4-5} \n")
	file.write("matrices & \eqref{Asp07} & \eqref{Asp08} & \eqref{Asp07} & \eqref{Asp08} \\ \n")
	file.write("\\\ \midrule \n")
	totalDiff07l = 0
	totalnInstances07l = 0
	totalDiff07r = 0
	totalnInstances07r = 0
	totalDiff08l = 0
	totalnInstances08l = 0
	totalDiff08r = 0
	totalnInstances08r = 0
	i = 0
	for instancetype in instancesets:
		diff07l = 0
		nInstances07l = 0
		diff07r = 0
		nInstances07r = 0
		diff08l = 0
		nInstances08l = 0
		diff08r = 0
		nInstances08r = 0
		for j1,j2 in instancetype:
			# we only take instances with alpha_k > 0, since otherwise the trivial bound is optimal and we cannot compute a relative gap
			if dualresults[1][j1] != -1e20 and gaps[0][j1] != "infinite" and float(gaps[0][j1]) < epsilon and float(dualresults[0][j1]) > epsilon:
				diff07l += (float(dualresults[0][j1]) - float(dualresults[1][j1])) / float(dualresults[0][j1])
				totalDiff07l += (float(dualresults[0][j1]) - float(dualresults[1][j1])) / float(dualresults[0][j1])
				nInstances07l += 1
				totalnInstances07l += 1
			if dualresults[1][j2] != -1e20 and gaps[0][j2] != "infinite" and float(gaps[0][j2]) < epsilon and dualresults[0][j2] != float("inf"):
				diff07r += (float(dualresults[1][j2]) - float(dualresults[0][j2])) / float(dualresults[0][j2])
				totalDiff07r += (float(dualresults[1][j2]) - float(dualresults[0][j2])) / float(dualresults[0][j2])
				nInstances07r += 1
				totalnInstances07r += 1
			if dualresults[2][j1] != -1e20 and gaps[0][j1] != "infinite" and float(gaps[0][j1]) < epsilon and float(dualresults[0][j1]) > epsilon and dualresults[2][j1] != -float("inf"):
				diff08l += (float(dualresults[0][j1]) - float(dualresults[2][j1])) / float(dualresults[0][j1])
				nInstances08l += 1
				totalDiff08l += (float(dualresults[0][j1]) - float(dualresults[2][j1])) / float(dualresults[0][j1])
				totalnInstances08l += 1
			if dualresults[2][j2] != -1e20 and gaps[0][j2] != "infinite" and float(gaps[0][j2]) < epsilon and dualresults[2][j2] != float("inf"):
				diff08r += (float(dualresults[2][j2]) - float(dualresults[0][j2])) / float(dualresults[0][j2])
				nInstances08r += 1
				totalDiff08r += (float(dualresults[2][j2]) - float(dualresults[0][j2])) / float(dualresults[0][j2])
				totalnInstances08r += 1
		if nInstances07l > 0:
			avgdiff07l = 100 * diff07l / nInstances07l
		else:
			avgdiff07l = "-"
		if nInstances07r > 0:
			avgdiff07r = 100 * diff07r / nInstances07r
		else:
			avgdiff07r = "-"
		if nInstances08l > 0:
			avgdiff08l = 100 * diff08l / nInstances08l
		else:
			avgdiff08l = "-"
		if nInstances08r > 0:
			avgdiff08r = 100 * diff08r / nInstances08r
		else:
			avgdiff08r = "-"
		file.write(instancesetnames[i])
		if avgdiff07l != "-":
			file.write(" & \\num{" + "%.2f" % avgdiff07l + "}\,\% (" + str(nInstances07l) + ") &")
		else:
			file.write("& -- (0) &")
		if avgdiff08l != "-":
			file.write(" \\num{" + "%.2f" % avgdiff08l + "}\,\% (" + str(nInstances08l) + ") &")
		else:
			file.write("-- (0) &")
		if avgdiff07r != "-":
			file.write(" \\num{" + "%.2f" % avgdiff07r + "}\,\% (" + str(nInstances07r) + ") &")
		else:
			file.write(" -- (0) &")
		if avgdiff08r != "-":
			file.write(" \\num{" + "%.2f" % avgdiff08r + "}\,\% (" + str(nInstances08r) + ") \\\ \n")
		else:
			file.write("-- (0) \\\ \n")
		i += 1
	file.write("\\midrule \n")
	if totalnInstances07l > 0:
		totalavgdiff07l = 100 * totalDiff07l / totalnInstances07l
	else:
		totalavgdiff07l = "-"
	if totalnInstances07r > 0:
		totalavgdiff07r = 100 * totalDiff07r / totalnInstances07r
	else:
		totalavgdiff07r = "-"
	if nInstances08l > 0:
		totalavgdiff08l = 100 * totalDiff08l / totalnInstances08l
	else:
		totalavgdiff08l = "-"
	if nInstances08r > 0:
		totalavgdiff08r = 100 * totalDiff08r / totalnInstances08r
	else:
		totalavgdiff08r = "-"
	file.write("total & ")
	if avgdiff07l != "-- (0)":
		file.write(" \\num{" + "%.2f" % totalavgdiff07l + "}\,\% (" + str(totalnInstances07l) + ") &")
	else:
		file.write("-- (0) &")
	if avgdiff08l != "-":
		file.write(" \\num{" + "%.2f" % totalavgdiff08l + "}\,\% (" + str(totalnInstances08l) + ") &")
	else:
		file.write("-- (0) &")
	if avgdiff07r != "-":
		file.write(" \\num{" + "%.2f" % totalavgdiff07r + "}\,\% (" + str(totalnInstances07r) + ") &")
	else:
		file.write("-- (0) &")
	if avgdiff08r != "-":
		file.write(" \\num{" + "%.2f" % totalavgdiff08r + "}\,\% (" + str(totalnInstances08r) + ") \\\ \n")
	else:
		file.write("- \\\ \n")
	file.write("\\bottomrule \n \\end{tabular*} \n \end{scriptsize} \n \\end{table} \n")


def TimeTable(instancesets, instancesetnames, caption, label):
	file.write("\\begin{table} \n \\begin{scriptsize} \\caption{" + caption + "} \n \\label{" + label + "} \n \\begin{tabular*}{\\linewidth}{@{}l@{\\;\\;\extracolsep{\\fill}}rrr")
	file.write("@{}}\\toprule \n")
	file.write(" matrices & \eqref{MISP} & \eqref{Asp07} & \eqref{Asp08} ")
	file.write("\\\ \midrule \n")
	totalMISDPtime = 1.0
	totalMISDPnum = 0
	totalA07time = 1.0
	totalA07num = 0
	totalA08time = 1.0
	totalA08num = 0
	i = 0
	for instancetype in instancesets:
		MISDPtime = 1.0
		MISDPnum = 0
		A07time = 1.0
		A07num = 0
		A08time = 1.0
		A08num = 0
		for j in instancetype:
			MISDPtime = math.pow(MISDPtime, float(MISDPnum) / float(MISDPnum+1)) * math.pow(float(times[0][j] + timeshift), 1.0/float(MISDPnum+1))
			MISDPnum += 1
			totalMISDPtime = math.pow(totalMISDPtime, float(totalMISDPnum) / float(totalMISDPnum+1)) * math.pow(float(times[0][j] + timeshift), 1.0/float(totalMISDPnum+1))
			totalMISDPnum += 1
			A07time = math.pow(A07time, float(A07num) / float(A07num+1)) * math.pow(float(times[1][j] + timeshift), 1.0/float(A07num+1))
			A07num += 1
			totalA07time = math.pow(totalA07time, float(totalA07num) / float(totalA07num+1)) * math.pow(float(times[1][j] + timeshift), 1.0/float(totalA07num+1))
			totalA07num += 1
			A08time = math.pow(A08time, float(A08num) / float(A08num+1)) * math.pow(float(times[2][j] + timeshift), 1.0/float(A08num+1))
			A08num += 1
			totalA08time = math.pow(totalA08time, float(totalA08num) / float(totalA08num+1)) * math.pow(float(times[2][j] + timeshift), 1.0/float(totalA08num+1))
			totalA08num += 1
		MISDPtime = MISDPtime - timeshift
		A07time = A07time - timeshift
		A08time = A08time - timeshift
		file.write(instancesetnames[i] + "& \\num{%.1f" % MISDPtime + "} & " + "\\num{%.1f" % A07time + "} & " + "\\num{%.1f" % A08time + "} \\\ \n ")
		i += 1
	file.write("\\midrule \n")
	totalMISDPtime = totalMISDPtime - timeshift
	totalA07time = totalA07time - timeshift
	totalA08time = totalA08time - timeshift
	file.write("total & \\num{%.1f" % totalMISDPtime + "} & " + "\\num{%.1f" % totalA07time + "} & " + "\\num{%.1f" % totalA08time + "} \\\ \n ")
	file.write("\\bottomrule \n \\end{tabular*} \n \end{scriptsize} \n \\end{table} \n")

def LhsRhsTimeTable(instancesets, instancesetnames, caption, label):
	file.write("\\begin{table} \n \\begin{scriptsize} \\caption{" + caption + "} \n \\label{" + label + "} \n \\begin{tabular*}{\\linewidth}{@{}l@{\\;\\;\extracolsep{\\fill}}rrrrrr")
	file.write("@{}}\\toprule \n")
	file.write("  & \\multicolumn{3}{c}{$\\alpha_k$} & \\multicolumn{3}{c}{$\\beta_k$} \\\ \n")
	file.write("\\cmidrule(r){2-4} \\cmidrule(l){5-7} \n")
	file.write(" matrices & \eqref{MISDP} & \eqref{Asp07} & \eqref{Asp08} & \eqref{MISDP} & \eqref{Asp07} & \eqref{Asp08} ")
	file.write("\\\ \midrule \n")
	totalMISDPtimel = 1.0
	totalMISDPnuml = 0
	totalA07timel = 1.0
	totalA07numl = 0
	totalA08timel = 1.0
	totalA08numl = 0
	totalMISDPtimer = 1.0
	totalMISDPnumr = 0
	totalA07timer = 1.0
	totalA07numr = 0
	totalA08timer = 1.0
	totalA08numr = 0
	i = 0
	for instancetype in instancesets:
		MISDPtimel = 1.0
		MISDPnuml = 0
		A07timel = 1.0
		A07numl = 0
		A08timel = 1.0
		A08numl = 0
		MISDPtimer = 1.0
		MISDPnumr = 0
		A07timer = 1.0
		A07numr = 0
		A08timer = 1.0
		A08numr = 0
		for j1,j2 in instancetype:
			if j1 > -1:
				MISDPtimel = math.pow(MISDPtimel, float(MISDPnuml) / float(MISDPnuml+1)) * math.pow(float(times[0][j1] + timeshift), 1.0/float(MISDPnuml+1))
				MISDPnuml += 1
				totalMISDPtimel = math.pow(totalMISDPtimel, float(totalMISDPnuml) / float(totalMISDPnuml+1)) * math.pow(float(times[0][j1] + timeshift), 1.0/float(totalMISDPnuml+1))
				totalMISDPnuml += 1
				A07timel = math.pow(A07timel, float(A07numl) / float(A07numl+1)) * math.pow(float(times[1][j1] + timeshift), 1.0/float(A07numl+1))
				A07numl += 1
				totalA07timel = math.pow(totalA07timel, float(totalA07numl) / float(totalA07numl+1)) * math.pow(float(times[1][j1] + timeshift), 1.0/float(totalA07numl+1))
				totalA07numl += 1
				A08timel = math.pow(A08timel, float(A08numl) / float(A08numl+1)) * math.pow(float(times[2][j1] + timeshift), 1.0/float(A08numl+1))
				A08numl += 1
				totalA08timel = math.pow(totalA08timel, float(totalA08numl) / float(totalA08numl+1)) * math.pow(float(times[2][j1] + timeshift), 1.0/float(totalA08numl+1))
				totalA08numl += 1
			else:
				MISDPtimel = -1.0
			MISDPtimer = math.pow(MISDPtimer, float(MISDPnumr) / float(MISDPnumr+1)) * math.pow(float(times[0][j2] + timeshift), 1.0/float(MISDPnumr+1))
			MISDPnumr += 1
			totalMISDPtimer = math.pow(totalMISDPtimer, float(totalMISDPnumr) / float(totalMISDPnumr+1)) * math.pow(float(times[0][j2] + timeshift), 1.0/float(totalMISDPnumr+1))
			totalMISDPnumr += 1
			A07timer = math.pow(A07timer, float(A07numr) / float(A07numr+1)) * math.pow(float(times[1][j2] + timeshift), 1.0/float(A07numr+1))
			A07numr += 1
			totalA07timer = math.pow(totalA07timer, float(totalA07numr) / float(totalA07numr+1)) * math.pow(float(times[1][j2] + timeshift), 1.0/float(totalA07numr+1))
			totalA07numr += 1
			A08timer = math.pow(A08timer, float(A08numr) / float(A08numr+1)) * math.pow(float(times[2][j2] + timeshift), 1.0/float(A08numr+1))
			A08numr += 1
			totalA08timer = math.pow(totalA08timer, float(totalA08numr) / float(totalA08numr+1)) * math.pow(float(times[2][j2] + timeshift), 1.0/float(totalA08numr+1))
			totalA08numr += 1
		MISDPtimel = MISDPtimel - timeshift
		A07timel = A07timel - timeshift
		A08timel = A08timel - timeshift
		MISDPtimer = MISDPtimer - timeshift
		A07timer = A07timer - timeshift
		A08timer = A08timer - timeshift
		if MISDPtimel >= 0.0:
			file.write(instancesetnames[i] + "& \\num{%.1f" % MISDPtimel + "} & " + "\\num{%.1f" % A07timel + "} & " + "\\num{%.1f" % A08timel + "} ")
		else:
			file.write(instancesetnames[i] + "& -- & -- & -- ")
		file.write("& \\num{%.1f" % MISDPtimer + "} & " + "\\num{%.1f" % A07timer + "} & " + "\\num{%.1f" % A08timer + "} \\\ \n ")
		i += 1
	file.write("\\midrule \n")
	totalMISDPtimel = totalMISDPtimel - timeshift
	totalA07timel = totalA07timel - timeshift
	totalA08timel = totalA08timel - timeshift
	totalMISDPtimer = totalMISDPtimer - timeshift
	totalA07timer = totalA07timer - timeshift
	totalA08timer = totalA08timer - timeshift
	file.write("total & \\num{%.1f" % totalMISDPtimel + "} & " + "\\num{%.1f" % totalA07timel + "} & " + "\\num{%.1f" % totalA08timel + "} ")
	file.write("& \\num{%.1f" % totalMISDPtimer + "} & " + "\\num{%.1f" % totalA07timer + "} & " + "\\num{%.1f" % totalA08timer + "} \\\ \n ")
	file.write("\\bottomrule \n \\end{tabular*} \n \end{scriptsize} \n \\end{table} \n")

def RICtable(instances, caption, label, tabularx):
	if tabularx:
		file.write("\\newpage \n \\begin{scriptsize} \n  \\setlength{\\tabcolsep}{2pt} \n \\tablehead{\\toprule \n")
		file.write("matrix & $m$ & $n$ & $k$ & $\\alpha_k^2$ & $\\beta_k^2$ & $\\gamma_k^2$ & $\\delta_k^2$")
		file.write("\\\ \\midrule} \n \\tabletail{ \\midrule \\multicolumn{3}{@{}l}{continued on next page \\dots}\\\ \\bottomrule} \n \\tablelasttail{\\bottomrule} \n \\tablecaption[")
		file.write(caption)
		file.write("]{")
		file.write(caption)
		file.write("}\label{")
		file.write(label)
		file.write("}\n")
		file.write("\\begin{xtabular*}{\\linewidth}{@{\extracolsep{\\fill}}lrrrrrrr@{}} \n ")
	else:
		file.write("\\begin{table} \n \\begin{scriptsize} \\caption{" + caption + "} \n \\label{" + label + "} \n \\begin{tabular*}{\\linewidth}{@{}l@{\\;\\;\extracolsep{\\fill}}rrrrrrr @{}}\\toprule \n")
		file.write("matrix & $m$ & $n$ & $k$ & $\\alpha_k^2$ & $\\beta_k^2$ & $\\gamma_k^2$ & $\\delta_k^2$")
		file.write("\\\ \midrule \n")
	i = 0
	j = 0
	correct = [-0.00574, -0.01336, -0.01220]
	m = [15,15,15,25,25,25,30,30,30,15,15,15,25,25,25,30,30,30,30,30,30,35,35,35,40,40,40,30,30,30,35,35,35,40,40,40,15,15,15,25,25,25,30,30,30,15,15,15,25,25,25,30,30,30,15,15,15,25,25,25,30,30,30]
	n = [30,30,30,35,35,35,40,40,40,30,30,30,35,35,35,40,40,40,30,30,30,35,35,35,40,40,40,30,30,30,35,35,35,40,40,40,30,30,30,35,35,35,40,40,40,30,30,30,35,35,35,40,40,40,30,30,30,35,35,35,40,40,40]
	k = [5,5,5,4,4,4,3,3,3,5,5,5,4,4,4,3,3,3,5,5,5,4,4,4,3,3,3,5,5,5,4,4,4,3,3,3,5,5,5,4,4,4,3,3,3,5,5,5,4,4,4,3,3,3,5,5,5,4,4,4,3,3,3]
	for j1,j2 in instances:
		if primalresults[0][j1] == 1e20:
			alpha = correct[j]
			j += 1
		else:
			alpha = float(primalresults[0][j1])
		beta = float(primalresults[0][j2])
		if alpha > 0.0001:
			gamma = beta / alpha
		else:
			gamma = "-"
		if 1 - alpha > beta - 1:
			delta = 1-alpha
		else:
			delta = beta-1
		file.write(names[0][j1].split("_")[0] + " & \\num{%.0f" % float(m[j1/2]) + "} & \\num{%.0f" % float(n[j1/2]) + "} & \\num{%.0f" % float(k[j1/2]) + "} & \\num{%.2f" % float(alpha) + "} & \\num{%.2f" % float(beta) + "} & ")
		if gamma == "-":
			file.write("-- & ")
		else:
			file.write("\\num{%.2f" % float(gamma) + "} & ")
		file.write("\\num{%.2f" % float(delta) + "} \\\ \n")
		i += 1
	if tabularx:
		file.write("\\end{xtabular*} \n \\end{scriptsize} \n")
	else:
		file.write("\\bottomrule \n \\end{tabular*} \n \end{scriptsize} \n \\end{table} \n")

def printUnsolved(instanceset):
	memory07 = 0
	fail07 = 0
	memory08 = 0
	fail08 = 0
	totalunsolved08 = 0
	for j in instanceset:
		if Asp07unsolved[j] < 0:
			memory07 += 1
		elif Asp07unsolved[j] > 0:
			fail07 += 1
		allfailed = True
		for s in range(nAsp08Steps):
			if Asp08unsolved[j][s] < 0:
				memory08 += 1
			elif Asp08unsolved[j][s] > 0:
				fail08 += 1
			if Asp08unsolved[j][s] >= 0:
				allfailed = False
				continue
		if allfailed:
			totalunsolved08 += 1
	print("number of memory fails d'Aspremont 2007: " + str(memory07))
	print("number of solver fails d'Aspremont 2007: " + str(fail07))
	print("number of memory fails d'Aspremont 2008: " + str(memory08))
	print("number of solver fails d'Aspremont 2007: " + str(fail08))
	print("number of complete fails d'Aspremont 2007: " + str(totalunsolved08))


if __name__=="__main__":
	"""give any number of .out-files for the same testset, then loops over them and returns a .tex-file given as first argument with some tables and a performance graph """
# only one argument: the tex-File
	ninstances=[0 for x in range(len(sys.argv) -2)] #initialize ninstances matrix
	names=[[0 for x in range(126)] for x in range(3)] #initialize names matrix
	dualresults=[[0 for x in range(126)] for x in range(3)] #initialize dual results matrix
	primalresults=[[0 for x in range(126)] for x in range(3)] #initialize primal results matrix
	gaps=[[0 for x in range(126)] for x in range(3)]   #initialize gaps matrix
	times=[[0 for x in range(126)] for x in range(3)]   #initialize solvingtime matrix (for Asp08 gives time of best mu)
	sdpiters=[[0 for x in range(126)] for x in range(3)]   #initialize sdp iterations matrix (for Asp08 gives iters of best mu)
	nodes=[[0 for x in range(126)] for x in range(1)]   #initialize nodes matrix
	unsolved=[[0 for x in range(126)] for x in range(1)]   #initialize unsolved matrix
	penalties=[[0 for x in range(126)] for x in range(1)]   #initialize penalties matrix
	Asp08dualresults=[[0 for x in range(15)] for x in range(126)] #initialize dual results matrix for Asp08
	Asp08primalresults=[[0 for x in range(15)] for x in range(126)] #initialize primal results matrix for Asp08
	Asp08gaps=[[0 for x in range(15)] for x in range(126)]   #initialize gaps matrix for Asp08
	Asp08times=[[0 for x in range(15)] for x in range(126)]   #initialize solvingtime matrix for Asp08
	Asp08sdpiters=[[0 for x in range(15)] for x in range(126)]   #initialize sdp iterations matrix for Asp08
	Asp08unsolved=[[0 for x in range(15)] for x in range(126)]	#initialize unsolved matrix for Asp08 | -1 = memory limit +1 = not converged
	Asp07unsolved=[0 for x in range(126)]	#initialize unsolved matrix for Asp07

	#read the results
	readFileMISDP(MISDPfilename)

	for j in range(len(Asp07filenames)):
		readFileAsp07(Asp07filenames[j], j, j % 2)
		readFileAsp08(Asp08basenames[j], Asp08basenames[j], j, j % 2)

	texfilename = sys.argv[1]
	file=open(texfilename, "w")
	if completeResultsTable:
		file.write("\\documentclass[landscape]{article} \n \\usepackage{amsmath} \n \\usepackage{amsthm} \n \\usepackage{xspace} \n \\usepackage{dsfont} \n \\usepackage[dvipsnames]{xcolor} \n \\usepackage{booktabs} \n \\usepackage{multirow} \n \\usepackage{mathtools} \n \\usepackage{xtab} \n \usepackage{tikz} \n \usepackage{pgfplots} \n  \usepackage[margin=1in,footskip=0.25in]{geometry} \n \\extrafloats{100} \n \\usepackage[detect-weight,group-minimum-digits = 4]{siunitx} \n \\usepackage{tabularx} \n \\newcommand{\\setting}[1]{\\texttt{#1}} \n \\newcommand{\\sdpa}{\\textsc{SDPA}\\xspace} \n \\begin{document} \n")
		CompleteResultsTable()

	if texfile:
		file.write("\\documentclass[landscape]{article} \n \\usepackage{amsmath} \n \\usepackage{amsthm} \n \\usepackage{dsfont} \n \\usepackage[dvipsnames]{xcolor} \n \\usepackage{booktabs} \n \\usepackage{multirow} \n \\usepackage{mathtools} \n \\usepackage{xtab} \n \usepackage{tikz} \n \usepackage{pgfplots} \n  \usepackage[margin=1in,footskip=0.25in]{geometry} \n \\extrafloats{100} \n \\usepackage{sistyle} \n \\usepackage{tabularx} \n \\newcommand{\\setting}[1]{\\texttt{#1}} \n \\begin{document} \n")

	if Lbounds:
		LhsResultsTable([[54,56,58,60,62,64,66,68,70],[18,20,22,24,26,28,30,32,34],[90,92,94,96,98,100,102,104,106],[72,74,76,78,80,82,84,86,88],[36,38,40,42,44,46,48,50,52],[0,2,4,6,8,10,12,14,16]], ["$N(0,1)$", "binary", "band matrix", "$N(0,1/m)$", "$\\pm 1/\\sqrt{m}$", "$0, \\pm \\sqrt{3/m}$"], "Average gap of relaxations for RIC $\\alpha_k$", "lhsGap")

	if Rbounds:
		RhsResultsTable([[55,57,59,61,63,65,67,69,71],[19,21,23,25,27,29,31,33,35],[91,93,95,97,99,101,103,105,107],[109,111,113,115,117,119,121,123,125],[73,75,77,79,81,83,85,87,89],[37,39,41,43,45,47,49,51,53],[1,3,5,7,9,11,13,15,17]], ["$N(0,1)$", "binary", "band matrix", "rank 1", "$N(0,1/m)$", "$\\pm 1/\\sqrt{m}$", "$0, \\pm \\sqrt{3/m}$"], "Average gap of relaxations for RIC $\\beta_k$", "rhsGap")

	if LRbounds:
		LhsRhsResultsTable([[[54,55],[56,57],[58,59],[60,61],[62,63],[64,65],[66,67],[68,69],[70,71]],[[18,19],[20,21],[22,23],[24,25],[26,27],[28,29],[30,31],[32,33],[34,35]],[[90,91],[92,93],[94,95],[96,97],[98,99],[100,101],[102,103],[104,105],[106,107]],[[108,109],[110,111],[112,113],[114,115],[116,117],[118,119],[120,121],[122,123],[124,125]],[[72,73],[74,75],[76,77],[78,79],[80,81],[82,83],[84,85],[86,87],[88,89]],[[36,37],[38,39],[40,41],[42,43],[44,45],[46,47],[48,49],[50,51],[52,53]],[[0,1],[2,3],[4,5],[6,7],[8,9],[10,11],[12,13],[14,15],[16,17]]], ["$N(0,1)$", "binary", "band matrix", "rank 1", "$N(0,1/m)$", "$\\pm 1/\\sqrt{m}$", "$0, \\pm \\sqrt{3/m}$"], "Average gap of relaxations for RICs", "lhsRhsGap")

	if Ltimes:
		TimeTable([[54,56,58,60,62,64,66,68,70],[18,20,22,24,26,28,30,32,34],[90,92,94,96,98,100,102,104,106],[72,74,76,78,80,82,84,86,88],[36,38,40,42,44,46,48,50,52],[0,2,4,6,8,10,12,14,16]], ["$N(0,1)$", "binary", "band matrix", "$N(0,1/m)$", "$\\pm 1/\\sqrt{m}$", "$0, \\pm \\sqrt{3/m}$"], "Shifted geometric mean of solving times for left-hand side of RIP", "lhsTime")

	if Rtimes:
		TimeTable([[55,57,59,61,63,65,67,69,71],[19,21,23,25,27,29,31,33,35],[91,93,95,97,99,101,103,105,107],[109,111,113,115,117,119,121,123,125],[73,75,77,79,81,83,85,87,89],[37,39,41,43,45,47,49,51,53],[1,3,5,7,9,11,13,15,17]], ["$N(0,1)$", "binary", "band matrix", "rank 1", "$N(0,1/m)$", "$\\pm 1/\\sqrt{m}$", "$0, \\pm \\sqrt{3/m}$"], "Shifted geometric mean of solving times for right-hand side of RIP", "rhsTime")
	if LRtimes:
		LhsRhsTimeTable([[[54,55],[56,57],[58,59],[60,61],[62,63],[64,65],[66,67],[68,69],[70,71]],[[18,19],[20,21],[22,23],[24,25],[26,27],[28,29],[30,31],[32,33],[34,35]],[[90,91],[92,93],[94,95],[96,97],[98,99],[100,101],[102,103],[104,105],[106,107]],[[-1,109],[-1,111],[-1,113],[-1,115],[-1,117],[-1,119],[-1,121],[-1,123],[-1,125]],[[72,73],[74,75],[76,77],[78,79],[80,81],[82,83],[84,85],[86,87],[88,89]],[[36,37],[38,39],[40,41],[42,43],[44,45],[46,47],[48,49],[50,51],[52,53]],[[0,1],[2,3],[4,5],[6,7],[8,9],[10,11],[12,13],[14,15],[16,17]]], ["$N(0,1)$", "binary", "band matrix", "rank 1", "$N(0,1/m)$", "$\\pm 1/\\sqrt{m}$", "$0, \\pm \\sqrt{3/m}$"], "Shifted geometric mean of solving times for RICs", "lhsRhsTime")

	if printFails:
		printUnsolved([0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84,85,86,87,88,89,90,91,92,93,94,95,96,97,98,99,100,101,102,103,104,105,106,107,109,111,113,115,117,119,121,123,125])

	if CompleteRICtable:
		RICtable([[54,55],[56,57],[58,59],[60,61],[62,63],[64,65],[66,67],[68,69],[70,71],[18,19],[20,21],[22,23],[24,25],[26,27],[28,29],[30,31],[32,33],[34,35],[90,91],[92,93],[94,95],[96,97],[98,99],[100,101],[102,103],[104,105],[106,107],[108,109],[110,111],[112,113],[114,115],[116,117],[118,119],[120,121],[122,123],[124,125],[72,73],[74,75],[76,77],[78,79],[80,81],[82,83],[84,85],[86,87],[88,89],[36,37],[38,39],[40,41],[42,43],[44,45],[46,47],[48,49],[50,51],[52,53],[0,1],[2,3],[4,5],[6,7],[8,9],[10,11],[12,13],[14,15],[16,17]], "List of all created matrices", "MatrixList",1)

	if texfile  or completeResultsTable:
		file.write("\\end{document}")

	file.close()
