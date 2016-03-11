import os
import sys
import math
from decimal import Decimal

#TODO: write tables for online supplement

#which tables/plots for the paper should be prepared?
texfile = 0 # make the output a compilable texfile or just a figure/table that can be included?
Lbounds = 0
Rbounds = 0
Ltimes = 0
Rtimes = 1
completeTable = 0

MISDPfilename = "/local/gally/results/RIP-MISDP-Paper/160311/RIP-results/check.RIPMISDP.scipsdp.linux.x86_64.gnu.opt.sdpa.extra.branchinfobj_nofracdive.out"
Asp0708path = "/local/gally/results/RIP-MISDP-Paper/160311/SDPA-results/"
Asp07filenames = ["0+-115305A_07l.out", "0+-115305A_07r.out", "0+-115305B_07l.out", "0+-115305B_07r.out", "0+-115305C_07l.out", "0+-115305C_07r.out", "0+-125354A_07l.out", "0+-125354A_07r.out", "0+-125354B_07l.out", "0+-125354B_07r.out", "0+-125354C_07l.out", "0+-125354C_07r.out", "0+-130403A_07l.out", "0+-130403A_07r.out", "0+-130403B_07l.out", "0+-130403B_07r.out", "0+-130403C_07l.out", "0+-130403C_07r.out", "bina15305A_07l.out", "bina15305A_07r.out", "bina15305B_07l.out", "bina15305B_07r.out", "bina15305C_07l.out", "bina15305C_07r.out", "bina25354A_07l.out", "bina25354A_07r.out", "bina25354B_07l.out", "bina25354B_07r.out", "bina25354C_07l.out", "bina25354C_07r.out", "bina30403A_07l.out", "bina30403A_07r.out", "bina30403B_07l.out", "bina30403B_07r.out", "bina30403C_07l.out", "bina30403C_07r.out", "bern15305A_07l.out", "bern15305A_07r.out", "bern15305B_07l.out", "bern15305B_07r.out", "bern15305C_07l.out", "bern15305C_07r.out", "bern25354A_07l.out", "bern25354A_07r.out", "bern25354B_07l.out", "bern25354B_07r.out", "bern25354C_07l.out", "bern25354C_07r.out", "bern30403A_07l.out", "bern30403A_07r.out", "bern30403B_07l.out", "bern30403B_07r.out", "bern30403C_07l.out", "bern30403C_07r.out", "norm15305A_07l.out", "norm15305A_07r.out", "norm15305B_07l.out", "norm15305B_07r.out", "norm15305C_07l.out", "norm15305C_07r.out", "norm25354A_07l.out", "norm25354A_07r.out", "norm25354B_07l.out", "norm25354B_07r.out", "norm25354C_07l.out", "norm25354C_07r.out", "norm30403A_07l.out", "norm30403A_07r.out", "norm30403B_07l.out", "norm30403B_07r.out", "norm30403C_07l.out", "norm30403C_07r.out", "wish15305A_07l.out", "wish15305A_07r.out", "wish15305B_07l.out", "wish15305B_07r.out", "wish15305C_07l.out", "wish15305C_07r.out", "wish25354A_07l.out", "wish25354A_07r.out", "wish25354B_07l.out", "wish25354B_07r.out", "wish25354C_07l.out", "wish25354C_07r.out", "wish30403A_07l.out", "wish30403A_07r.out", "wish30403B_07l.out", "wish30403B_07r.out", "wish30403C_07l.out", "wish30403C_07r.out", "band30305A_07l.out", "band30305A_07r.out", "band30305B_07l.out", "band30305B_07r.out", "band30305C_07l.out", "band30305C_07r.out", "band35354A_07l.out", "band35354A_07r.out", "band35354B_07l.out", "band35354B_07r.out", "band35354C_07l.out", "band35354C_07r.out", "band40403A_07l.out", "band40403A_07r.out", "band40403B_07l.out", "band40403B_07r.out", "band40403C_07l.out", "band40403C_07r.out", "rnk130305A_07l.out", "rnk130305A_07r.out", "rnk130305B_07l.out", "rnk130305B_07r.out", "rnk130305C_07l.out", "rnk130305C_07r.out", "rnk135354A_07l.out", "rnk135354A_07r.out", "rnk135354B_07l.out", "rnk135354B_07r.out", "rnk135354C_07l.out", "rnk135354C_07r.out", "rnk140403A_07l.out", "rnk140403A_07r.out", "rnk140403B_07l.out", "rnk140403B_07r.out", "rnk140403C_07l.out", "rnk140403C_07r.out"]
Asp08basenames = ["0+-115305A_08l", "0+-115305A_08r", "0+-115305B_08l", "0+-115305B_08r", "0+-115305C_08l", "0+-115305C_08r", "0+-125354A_08l", "0+-125354A_08r", "0+-125354B_08l", "0+-125354B_08r", "0+-125354C_08l", "0+-125354C_08r", "0+-130403A_08l", "0+-130403A_08r", "0+-130403B_08l", "0+-130403B_08r", "0+-130403C_08l", "0+-130403C_08r", "bina15305A_08l", "bina15305A_08r", "bina15305B_08l", "bina15305B_08r", "bina15305C_08l", "bina15305C_08r", "bina25354A_08l", "bina25354A_08r", "bina25354B_08l", "bina25354B_08r", "bina25354C_08l", "bina25354C_08r", "bina30403A_08l", "bina30403A_08r", "bina30403B_08l", "bina30403B_08r", "bina30403C_08l", "bina30403C_08r", "bern15305A_08l", "bern15305A_08r", "bern15305B_08l", "bern15305B_08r", "bern15305C_08l", "bern15305C_08r", "bern25354A_08l", "bern25354A_08r", "bern25354B_08l", "bern25354B_08r", "bern25354C_08l", "bern25354C_08r", "bern30403A_08l", "bern30403A_08r", "bern30403B_08l", "bern30403B_08r", "bern30403C_08l", "bern30403C_08r", "norm15305A_08l", "norm15305A_08r", "norm15305B_08l", "norm15305B_08r", "norm15305C_08l", "norm15305C_08r", "norm25354A_08l", "norm25354A_08r", "norm25354B_08l", "norm25354B_08r", "norm25354C_08l", "norm25354C_08r", "norm30403A_08l", "norm30403A_08r", "norm30403B_08l", "norm30403B_08r", "norm30403C_08l", "norm30403C_08r", "wish15305A_08l", "wish15305A_08r", "wish15305B_08l", "wish15305B_08r", "wish15305C_08l", "wish15305C_08r", "wish25354A_08l", "wish25354A_08r", "wish25354B_08l", "wish25354B_08r", "wish25354C_08l", "wish25354C_08r", "wish30403A_08l", "wish30403A_08r", "wish30403B_08l", "wish30403B_08r", "wish30403C_08l", "wish30403C_08r", "band30305A_08l", "band30305A_08r", "band30305B_08l", "band30305B_08r", "band30305C_08l", "band30305C_08r", "band35354A_08l", "band35354A_08r", "band35354B_08l", "band35354B_08r", "band35354C_08l", "band35354C_08r", "band40403A_08l", "band40403A_08r", "band40403B_08l", "band40403B_08r", "band40403C_08l", "band40403C_08r", "rnk130305A_08l", "rnk130305A_08r", "rnk130305B_08l", "rnk130305B_08r", "rnk130305C_08l", "rnk130305C_08r", "rnk135354A_08l", "rnk135354A_08r", "rnk135354B_08l", "rnk135354B_08r", "rnk135354C_08l", "rnk135354C_08r", "rnk140403A_08l", "rnk140403A_08r", "rnk140403B_08l", "rnk140403B_08r", "rnk140403C_08l", "rnk140403C_08r"]
Asp08Datapath = "/local/gally/instances/RIP-MISDP-Paper/Asp08/"
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
			times[0][j] = 3600
		else:
			time = substring[substring.find("Solving Time (sec) : ") + 21:].split()[0]
			times[0][j] = float(time)
			if times[0][j] > 3600.0:
				times[0][j] = 3600
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
	if rhs:
		primalresults[1][j] *= -1
		dualresults[1][j] *= -1
	#if returncode != "pdOPT":
		#print(filename + " returned " + returncode) 
	file.close()

def readFileAsp08(basename, datafilename, j, rhs):
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
		returncode = substring[substring.find("phase.value  = ") + 15:].split()[0]
		#if returncode != "pdOPT":
			#print(filename + " returned " + returncode) 
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
	times[2][j] = 14400
	for i in range(nAsp08Steps):
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
				times[2][j] = Asp08times[j][i]
				sdpiters[2][j] = Asp08sdpiters[j][i]
		else:
			result = alpha - Asp08dualresults[j][i] - rho * k
			if result > dualresults[2][j]:
				dualresults[2][j] = result
				times[2][j] = Asp08times[j][i]
				sdpiters[2][j] = Asp08sdpiters[j][i]
				


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
	file.write("\\begin{table} \n \\begin{scriptsize} \\caption{" + caption + "} \n \\label{" + label + "} \n \\begin{tabular*}{0.48\\textwidth}{@{}l@{\\;\\;\extracolsep{\\fill}}rr")
	file.write("@{}}\\toprule \n")
	file.write(" type of matrix & d'Asp 07 & d'Asp 08 ")
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
	file.write("\\begin{table} \n \\begin{scriptsize} \\caption{" + caption + "} \n \\label{" + label + "} \n \\begin{tabular*}{0.48\\textwidth}{@{}l@{\\;\\;\extracolsep{\\fill}}rr")
	file.write("@{}}\\toprule \n")
	file.write(" type of matrix & d'Asp 07 & d'Asp 08 ")
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
				diff07 += (float(dualresults[1][j]) - float(dualresults[0][j])) / float(dualresults[1][j])
				totalDiff07 += (float(dualresults[1][j]) - float(dualresults[0][j])) / float(dualresults[1][j])
				nInstances07 += 1
				totalnInstances07 += 1
			if dualresults[2][j] != -1e20 and gaps[0][j] != "infinite" and float(gaps[0][j]) < epsilon and dualresults[2][j] != float("inf"):
				diff08 += (float(dualresults[2][j]) - float(dualresults[0][j])) / float(dualresults[1][j])
				nInstances08 += 1
				totalDiff08 += (float(dualresults[2][j]) - float(dualresults[0][j])) / float(dualresults[1][j])
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
		

def TimeTable(instancesets, instancesetnames, caption, label):
	file.write("\\begin{table} \n \\begin{scriptsize} \\caption{" + caption + "} \n \\label{" + label + "} \n \\begin{tabular*}{0.48\\textwidth}{@{}l@{\\;\\;\extracolsep{\\fill}}rrrr")
	file.write("@{}}\\toprule \n")
	file.write(" type of matrix & MISDP & d'Asp 07 & d'Asp 08 ")
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

if __name__=="__main__":
	"""give any number of .out-files for the same testset, then loops over them and returns a .tex-file given as first argument with some tables and a performance graph """
	ninstances=[0 for x in range(len(sys.argv) -2)] #initialize ninstances matrix
	names=[[0 for x in range(126)] for x in range(3)] #initialize names matrix
	dualresults=[[0 for x in range(126)] for x in range(3)] #initialize dual results matrix
	primalresults=[[0 for x in range(126)] for x in range(3)] #initialize primal results matrix
	gaps=[[0 for x in range(126)] for x in range(3)]   #initialize gaps matrix
	times=[[0 for x in range(126)] for x in range(3)]   #initialize solvingtime matrix (for Asp08 gives time of best mu)
	sdpiters=[[0 for x in range(126)] for x in range(3)]   #initialize sdp iterations matrix (for Asp08 gives iters of best mu)
	nodes=[[0 for x in range(126)] for x in range(1)]   #initialize nodes matrix
	Asp08dualresults=[[0 for x in range(15)] for x in range(126)] #initialize dual results matrix for Asp08 
	Asp08primalresults=[[0 for x in range(15)] for x in range(126)] #initialize primal results matrix for Asp08 
	Asp08gaps=[[0 for x in range(15)] for x in range(126)]   #initialize gaps matrix for Asp08 
	Asp08times=[[0 for x in range(15)] for x in range(126)]   #initialize solvingtime matrix for Asp08 
	Asp08sdpiters=[[0 for x in range(15)] for x in range(126)]   #initialize sdp iterations matrix for Asp08 

	#read the results
	readFileMISDP(MISDPfilename)

	for j in range(len(Asp07filenames)):
		readFileAsp07(Asp07filenames[j], j, j % 2)
		readFileAsp08(Asp08basenames[j], Asp08basenames[j], j, j % 2)

	texfilename = sys.argv[1]
	file=open(texfilename, "w")
	if texfile:
		file.write("\\documentclass[landscape]{article} \n \\usepackage{amsmath} \n \\usepackage{amsthm} \n \\usepackage{dsfont} \n \\usepackage[dvipsnames]{xcolor} \n \\usepackage{booktabs} \n \\usepackage{multirow} \n \\usepackage{mathtools} \n \\usepackage{xtab} \n \usepackage{tikz} \n \usepackage{pgfplots} \n  \usepackage[margin=1in,footskip=0.25in]{geometry} \n \\extrafloats{100} \n \\usepackage{sistyle} \n \\usepackage{tabularx} \n \\newcommand{\\setting}[1]{\\texttt{#1}} \n \\begin{document} \n")

	if Lbounds:
		LhsResultsTable([[54,56,58,60,62,64,66,68,70],[18,20,22,24,26,28,30,32,34],[90,92,94,96,98,100,102,104,106],[72,74,76,78,80,82,84,86,88],[36,38,40,42,44,46,48,50,52],[0,2,4,6,8,10,12,14,16]], ["$\\mathcal{N}(0,1)$", "binary", "band matrix", "$\\mathcal{N}(0,1/m)$", "$\\pm 1/\\sqrt{m}$", "$0, \\pm \\sqrt{3/m}$"], "Average gap of relaxations for RIC $\\alpha_k$", "lhsGap")

	if Rbounds:
		RhsResultsTable([[55,57,59,61,63,65,67,69,71],[19,21,23,25,27,29,31,33,35],[91,93,95,97,99,101,103,105,107],[109,111,113,115,117,119,121,123,125],[73,75,77,79,81,83,85,87,89],[37,39,41,43,45,47,49,51,53],[1,3,5,7,9,11,13,15,17]], ["$\\mathcal{N}(0,1)$", "binary", "band matrix", "rank 1", "$\\mathcal{N}(0,1/m)$", "$\\pm 1/\\sqrt{m}$", "$0, \\pm \\sqrt{3/m}$"], "Average gap of relaxations for RIC $\\beta_k$", "rhsGap")

	if Ltimes:
		TimeTable([[54,56,58,60,62,64,66,68,70],[18,20,22,24,26,28,30,32,34],[90,92,94,96,98,100,102,104,106],[72,74,76,78,80,82,84,86,88],[36,38,40,42,44,46,48,50,52],[0,2,4,6,8,10,12,14,16]], ["$\\mathcal{N}(0,1)$", "binary", "band matrix", "$\\mathcal{N}(0,1/m)$", "$\\pm 1/\\sqrt{m}$", "$0, \\pm \\sqrt{3/m}$"], "Solving times for left-hand side of RIP", "lhsTime")

	if Rtimes:
		TimeTable([[55,57,59,61,63,65,67,69,71],[19,21,23,25,27,29,31,33,35],[91,93,95,97,99,101,103,105,107],[109,111,113,115,117,119,121,123,125],[73,75,77,79,81,83,85,87,89],[37,39,41,43,45,47,49,51,53],[1,3,5,7,9,11,13,15,17]], ["$\\mathcal{N}(0,1)$", "binary", "band matrix", "rank 1", "$\\mathcal{N}(0,1/m)$", "$\\pm 1/\\sqrt{m}$", "$0, \\pm \\sqrt{3/m}$"], "Solving times for right-hand side of RIP", "rhsTime")


	if texfile:
		file.write("\\end{document}")
	file.close()
