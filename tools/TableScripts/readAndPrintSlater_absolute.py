import os
import sys
import math
from decimal import Decimal

#which tables/plots for the paper should be prepared?
texfile = 1 # make the output a compilable texfile or just a figure/table that can be included?
SlaterTableCLS = 0
SlaterTableMinK = 0
SlaterTableTT = 0
SlaterTableOverall = 0
SlaterSolvedTableCLS =0
SlaterSolvedTableMinK = 0
SlaterSolvedTableTT = 0
SlaterSolvedTableOverall = 0
completeTableSlater = 1
completeTableSolved = 1

SlaterTableShortCLS = 0
SlaterTableShortMinK = 0
SlaterTableShortTT = 0
SlaterTableShortOverall = 0
SlaterSolvedTableShortCLS = 0
SlaterSolvedTableShortMinK = 0
SlaterSolvedTableShortTT = 0
SlaterSolvedTableShortOverall = 0


printNonlyDiving = 0 #print the number of instances that could only be solved using the fractional diving heuristic
randrounding = 1 #if input includes randomized rounding instances

#first argument should be the tex-file to write to, rest should be the .out files
#example for a correct call: python readAndPrintSlater_absolute.py /home/gally/gally_Dissertation/resultfiles/SlaterResults180604/Slaterresults_complete.tex /home/gally/gally_Dissertation/resultfiles/SlaterResults180604/check.SCIPSDPpaper.scipsdp.linux.x86_64.gnu.opt.dsdp.moskito.1thread_slatercheck.out /home/gally/gally_Dissertation/resultfiles/SlaterResults180604/check.SCIPSDPpaper.scipsdp.linux.x86_64.gnu.opt.sdpa.moskito.1thread_slatercheck.out /home/gally/gally_Dissertation/resultfiles/SlaterResults180604/check.SCIPSDPpaper.scipsdp.linux.x86_64.gnu.opt.msk.moskito.1thread_slatercheck.out

maxinstances = 200 # maximum number of instances
epsilon = 0.001
# without randrounding
#completetablesettings = [0,2,3,4,5,6,7,8,9,10,16,17,18,19,20,21,22,28]
#completeTableCaptions = ["Complete results for DSDP with inference branching", "", "Complete results for DSDP with combined infeasibility/objective branching", "Complete results for DSDP with combined infeasibility/objective branching and dualfixing", "Complete results for DSDP with combined infeasibility/objective branching and fractional diving in all nodes with depth a multiple of 10", "Complete results for DSDP with combined infeasibility/objective branching and dual fixing and fractional diving in all nodes with depth a multiple of 10", "Complete results for DSDP with combined infeasibility/objective branching and without fractional diving", "Complete results for DSDP with combined infeasibility/objective branching and dual fixing and without fractional diving", "Complete results for DSDP with objective branching", "Complete results for DSDP with infeasibility branching", "Complete results for SDPA with inference branching", "", "", "", "", "", "Complete results for SDPA with combined infeasibility/objective branching", "Complete results for SDPA with combined infeasibility/objective branching and dualfixing", "Complete results for SDPA with combined infeasibility/objective branching and fractional diving in all nodes with depth a multiple of 10", "Complete results for SDPA with combined infeasibility/objective branching and dual fixing and fractional diving in all nodes with depth a multiple of 10", "Complete results for SDPA with combined infeasibility/objective branching and without fractional diving", "Complete results for SDPA with combined infeasibility/objective branching and dual fixing and without fractional diving", "Complete results for SDPA with objective branching", "", "", "", "", "", "Complete results for SDPA with infeasibility branching", "Complete results for SDPA with infeasibility branching and dualfixing", "Complete results for SDPA with infeasibility branching and fractional diving in all nodes with depth a multiple of 10", "Complete results for SDPA with infeasibility branching and dual fixing and fractional diving in all nodes with depth a multiple of 10", "Complete results for SDPA with infeasibility branching and without fractional diving", "Complete results for SDPA with infeasibility branching and dual fixing and without fractional diving"]
#completeTableLabels = ["CompleteDSDPinfer", "", "CompleteDSDPinfobj", "CompleteDSDPinfobjdualfix", "CompleteDSDPinfobjdive10", "CompleteDSDPinfobjdive10dualfix", "CompleteDSDPinfobjnodive", "CompleteDSDPinfobjnodivedualfix", "CompleteDSDPobj", "CompleteDSDPinf", "CompleteSDPAinfer", "", "", "", "", "", "CompleteSDPAinfobj", "CompleteSDPAinfobjdualfix", "CompleteSDPAinfobjdive10", "CompleteSDPAinfobjdive10dualfix", "CompleteSDPAinfobjnodive", "CompleteSDPAinfobjnodivedualfix", "CompleteSDPAobj", "", "", "", "", "", "CompleteSDPAinf", "CompleteSDPAinfdualfix", "CompleteSDPAinfdive10", "CompleteSDPAinfdive10dualfix", "CompleteSDPAinfnodive", "CompleteSDPAinfnodivedualfix"]
#with randrounding (new = 8,9,10,11,26,27,28,29)
completetablesettings = [0,1,2]
completeTableSlaterCaptions = ["Complete statistics of Slater condition for \\dsdp", "Complete statistics of Slater condition for \\sdpa","Complete statistics of Slater condition for \\mosek"]
completeTableSlaterShortCaptions = ["Slater statistics for \\dsdp", "Slater statistics for \\sdpa","Slater statistics for \\mosek"]
completeTableSlaterLabels = ["CompleteSlaterTestsDSDP", "CompleteSlaterTestsSDPA", "CompleteSlaterTestsMosek"]
completeTableSolvedHoldsCaptions = ["Complete statistics of solver fails with Slater condition holding for \\dsdp", "Complete statistics of solver fails with Slater condition holding for \\sdpa","Complete statistics of solver fails with Slater condition holding for \\mosek"]
completeTableSolvedHoldsShortCaptions = ["Solver fails with Slater for \\dsdp", "Solver fails with Slater for \\sdpa","Solver fails with Slater for \\mosek"]
completeTableSolvedHoldsLabels = ["CompleteSlaterSolvedHoldsDSDP", "CompleteSlaterSolvedHoldsSDPA", "CompleteSlaterSolvedHoldsMosek"]
completeTableSolvedFailsCaptions = ["Complete statistics of solver fails with Slater condition failing for \\dsdp", "Complete statistics of solver fails with Slater condition failing for \\sdpa","Complete statistics of solver fails with Slater condition failing for \\mosek"]
completeTableSolvedFailsShortCaptions = ["Solver fails without Slater for \\dsdp", "Solver fails without Slater for \\sdpa","Solver fails without Slater for \\mosek"]
completeTableSolvedFailsLabels = ["CompleteSlaterSolvedFailsDSDP", "CompleteSlaterSolvedFailsSDPA", "CompleteSlaterSolvedFailsMosek"]
completeTableSolvedInfeasibleCaptions = ["Complete statistics of solver fails with Slater condition showing infeasibility for \\dsdp", "Complete statistics of solver fails with Slater condition showing infeasibility for \\sdpa","Complete statistics of solver fails with Slater condition showing infeasibility for \\mosek"]
completeTableSolvedInfeasibleShortCaptions = ["Solver fails with infeasibility for \\dsdp", "Solver fails with infeasibility for \\sdpa","Solver fails with infeasibility for \\mosek"]
completeTableSolvedInfeasibleLabels = ["CompleteSlaterSolvedInfeasibleDSDP", "CompleteSlaterSolvedInfeasibleSDPA", "CompleteSlaterSolvedInfeasibleMosek"]

def readFile(arg,i):
	j=0;
	file=open(arg, "r")
	filestring = file.read()
	while filestring.find("@01") > -1:
		substring = filestring[filestring.find("@01"):filestring.find("@04")]
		assert(substring.find("@01 ") > -1)
		names[i][j] = os.path.basename(substring[substring.find("@01 ") + 3:].split()[0])
                if substring.find("Dual Slater   :") == -1:
                        dualholds[i][j] = "-"
                        dualholdspercent[i][j] = "-"
                        dualfails[i][j] = "-"
                        dualfailspercent[i][j] = "-"
                        dualinfeasible[i][j] = "-"
                        dualinfeasiblepercent[i][j] = "-"
                        dualnoinfo[i][j] = "-"
                        dualnoinfopercent[i][j] = "-"
                else:
                        dualline = substring[substring.find("Dual Slater   :"):]
                        duallineitems = dualline.split()
                        dualholds[i][j] = float(duallineitems[3])
                        dualfails[i][j] = float(duallineitems[4])
                        dualinfeasible[i][j] = float(duallineitems[5])
                        dualnoinfo[i][j] = float(duallineitems[6])
                        #compute percentages
                        dualholdspercent[i][j] = 100 * dualholds[i][j] / (dualholds[i][j] + dualfails[i][j] + dualinfeasible[i][j] + dualnoinfo[i][j])
                        dualfailspercent[i][j] = 100 * dualfails[i][j] / (dualholds[i][j] + dualfails[i][j] + dualinfeasible[i][j] + dualnoinfo[i][j])
                        dualinfeasiblepercent[i][j] = 100 * dualinfeasible[i][j] / (dualholds[i][j] + dualfails[i][j] + dualinfeasible[i][j] + dualnoinfo[i][j])
                        dualnoinfopercent[i][j] = 100 * dualnoinfo[i][j] / (dualholds[i][j] + dualfails[i][j] + dualinfeasible[i][j] + dualnoinfo[i][j])
                if substring.find("Primal Slater :") == -1:
                        primalholds[i][j] = "-"
                        primalholdspercent[i][j] = "-"
                        primalfails[i][j] = "-"
                        primalfailspercent[i][j] = "-"
                        primalnoinfo[i][j] = "-"
                        primalnoinfopercent[i][j] = "-"
                else:
                        primalline = substring[substring.find("Primal Slater :"):]
                        primallineitems = primalline.split()
                        primalholds[i][j] = float(primallineitems[3])
                        primalfails[i][j] = float(primallineitems[4])
                        primalnoinfo[i][j] = float(primallineitems[6])
                        #compute percentages
                        primalholdspercent[i][j] = 100 * primalholds[i][j] / (primalholds[i][j] + primalfails[i][j] + primalnoinfo[i][j])
                        primalfailspercent[i][j] = 100 * primalfails[i][j] / (primalholds[i][j] + primalfails[i][j] + primalnoinfo[i][j])
                        primalnoinfopercent[i][j] = 100 * primalnoinfo[i][j] / (primalholds[i][j] + primalfails[i][j] + primalnoinfo[i][j])
                if substring.find("Slater holds  :") == -1:
                        slaterfast[i][j] = "-"
                        slaterstable[i][j] = "-"
                        slaterpenalty[i][j] = "-"
                        slaterbound[i][j] = "-"
                        slaterunsolved[i][j] = "-"
                else:
                        slaterline = substring[substring.find("Slater holds  :"):]
                        slaterlineitems = slaterline.split()
                        if not slaterlineitems[7].isdigit(): # MOSEK and DSDP do not have "stable" column
                                slaterfast[i][j] = slaterlineitems[3]
                                slaterstable[i][j] = 0
                                slaterpenalty[i][j] = slaterlineitems[4]
                                slaterbound[i][j] = slaterlineitems[5]
                                slaterunsolved[i][j] = slaterlineitems[6]
                        else:
                                slaterfast[i][j] = slaterlineitems[3]
                                slaterstable[i][j] = slaterlineitems[4]
                                slaterpenalty[i][j] = slaterlineitems[5]
                                slaterbound[i][j] = slaterlineitems[6]
                                slaterunsolved[i][j] = slaterlineitems[7]
                if substring.find("Slater fails  :") == -1:
                        noslaterfast[i][j] = "-"
                        noslaterstable[i][j] = "-"
                        noslaterpenalty[i][j] = "-"
                        noslaterbound[i][j] = "-"
                        noslaterunsolved[i][j] = "-"
                else:
                        noslaterline = substring[substring.find("Slater fails  :"):]
                        noslaterlineitems = noslaterline.split()
                        if not noslaterlineitems[7].isdigit(): # MOSEK and DSDP do not have "stable" column
                                noslaterfast[i][j] = noslaterlineitems[3]
                                noslaterstable[i][j] = 0
                                noslaterpenalty[i][j] = noslaterlineitems[4]
                                noslaterbound[i][j] = noslaterlineitems[5]
                                noslaterunsolved[i][j] = noslaterlineitems[6]
                        else:
                                noslaterfast[i][j] = noslaterlineitems[3]
                                noslaterstable[i][j] = noslaterlineitems[4]
                                noslaterpenalty[i][j] = noslaterlineitems[5]
                                noslaterbound[i][j] = noslaterlineitems[6]
                                noslaterunsolved[i][j] = noslaterlineitems[7]
                if substring.find("Infeasible    :") == -1:
                        infeasiblefast[i][j] = "-"
                        infeasiblestable[i][j] = "-"
                        infeasiblepenalty[i][j] = "-"
                        infeasiblebound[i][j] = "-"
                        infeasibleunsolved[i][j] = "-"
                else:
                        infeasibleline = substring[substring.find("Infeasible    :"):]
                        infeasiblelineitems = infeasibleline.split()
                        if not infeasiblelineitems[6].isdigit(): # MOSEK and DSDP do not have "stable" column
                                infeasiblefast[i][j] = infeasiblelineitems[2]
                                infeasiblestable[i][j] = 0
                                infeasiblepenalty[i][j] = infeasiblelineitems[3]
                                infeasiblebound[i][j] = infeasiblelineitems[4]
                                infeasibleunsolved[i][j] = infeasiblelineitems[5]
                        else:
                                infeasiblefast[i][j] = infeasiblelineitems[2]
                                infeasiblestable[i][j] = infeasiblelineitems[3]
                                infeasiblepenalty[i][j] = infeasiblelineitems[4]
                                infeasiblebound[i][j] = infeasiblelineitems[5]
                                infeasibleunsolved[i][j] = infeasiblelineitems[6]
		filestring = filestring[filestring.find("@04")+3:]
		j = j+1
		assert( j < maxinstances)
	ninstances[i] = j
	file.close()


def makeCompleteTableSlater(shortcaption, caption, label, file, i):
	file.write("\\newpage \n \\begin{scriptsize} \n  \\setlength{\\tabcolsep}{2pt} \n \\tablehead{\\toprule \n")
	file.write("\phantom{abc} & \multicolumn{4}{c}{Dual Slater} & \multicolumn{3}{c}{Primal Slater} \\\ ")
	file.write("\\cmidrule{2-5} \\cmidrule{6-8} \n")
	file.write("problem & \\cmark & \\xmark & inf & ? & \\cmark & \\xmark  & ? ")
	file.write("\\\ \\midrule} \n \\tabletail{ \\midrule \\multicolumn{3}{@{}l}{continued on next page \\dots}\\\ \\bottomrule} \n \\tablelasttail{\\bottomrule} \n \\tablecaption[")
	file.write(shortcaption)
	file.write("]{")
	file.write(caption)
	file.write("}\label{")
	file.write(label)
	file.write("}\n")
	file.write("\\begin{xtabular*}{\\textwidth}{@{\extracolsep{\\fill}}lrrrrrrr@{}} \n ")
	for j in range(ninstances[i]):
		file.write(names[i][j].replace("_", "\\_").split(".")[0] + "& ")
		if dualholds[i][j] == "-":
			file.write("-- & ")
		else:
			file.write("\\num{" + "%.0f" % float(dualholds[i][j]) + "} &")
		if dualfails[i][j] == "-":
			file.write("-- & ")
		else:
			file.write("\\num{" + "%.0f" % float(dualfails[i][j]) + "} &")
		if dualinfeasible[i][j] == "-":
			file.write("-- & ")
		else:
			file.write("\\num{" + "%.0f" % float(dualinfeasible[i][j]) + "} &")
		if dualnoinfo[i][j] == "-":
			file.write("-- & ")
		else:
			file.write("\\num{" + "%.0f" % float(dualnoinfo[i][j]) + "} &")
		if primalholds[i][j] == "-":
			file.write("-- & ")
		else:
			file.write("\\num{" + "%.0f" % float(primalholds[i][j]) + "} &")
		if primalfails[i][j] == "-":
			file.write("-- & ")
		else:
			file.write("\\num{" + "%.0f" % float(primalfails[i][j]) + "} &")
		if primalnoinfo[i][j] == "-":
			file.write("-- \\\ \n")
		else:
			file.write("\\num{" + "%.0f" % float(primalnoinfo[i][j]) + "} \\\ \n")
	file.write("  \\end{xtabular*} \n \\end{scriptsize} \n")

def makeCompleteTableSolved(shortcaption1, caption1, label1, shortcaption2, caption2, label2, shortcaption3, caption3, label3, file, i, stableColumn):
	file.write("\\newpage \n \\begin{scriptsize} \n  \\setlength{\\tabcolsep}{2pt} \n \\tablehead{\\toprule \n")
	file.write("problem & number & fast & ")
	if stableColumn:
		file.write("stable & ")
	file.write("penalty & bound & unssucc ")
	file.write("\\\ \\midrule} \n \\tabletail{ \\midrule \\multicolumn{3}{@{}l}{continued on next page \\dots}\\\ \\bottomrule} \n \\tablelasttail{\\bottomrule} \n \\tablecaption[")
	file.write(shortcaption1)
	file.write("]{")
	file.write(caption1)
	file.write("}\label{")
	file.write(label1)
	file.write("}\n")
	file.write("\\begin{xtabular*}{\\textwidth}{@{\extracolsep{\\fill}}lrrrrr")
	if stableColumn:
		file.write("r")
	file.write("@{}} \n ")
	for j in range(ninstances[i]):
		file.write(names[i][j].replace("_", "\\_").split(".")[0] + "& ")
		if slater[i][j] == "-":
			file.write("-- & ")
		else:
			file.write("\\num{" + "%.0f" % float(slater[i][j]) + "} &")
		if slaterfast[i][j] == "-":
			file.write("-- & ")
		else:
			file.write("\\num{" + "%.0f" % float(slaterfast[i][j]) + "} &")
		if stableColumn:
			if slaterstable[i][j] == "-":
				file.write("-- & ")
			else:
				file.write("\\num{" + "%.0f" % float(slaterstable[i][j]) + "} &")
		if slaterpenalty[i][j] == "-":
			file.write("-- & ")
		else:
			file.write("\\num{" + "%.0f" % float(slaterpenalty[i][j]) + "} &")
		if slaterbound[i][j] == "-":
			file.write("-- & ")
		else:
			file.write("\\num{" + "%.0f" % float(slaterbound[i][j]) + "} &")
		if slaterunsolved[i][j] == "-":
			file.write("-- \\\ \n")
		else:
			file.write("\\num{" + "%.0f" % float(slaterunsolved[i][j]) + "} \\\ \n")
	file.write("  \\end{xtabular*} \n \\end{scriptsize} \n")

	file.write("\\newpage \n \\begin{scriptsize} \n  \\setlength{\\tabcolsep}{2pt} \n \\tablehead{\\toprule \n")
	file.write("problem & number & fast & ")
	if stableColumn:
		file.write("stable & ")
	file.write("penalty & bound & unssucc ")
	file.write("\\\ \\midrule} \n \\tabletail{ \\midrule \\multicolumn{3}{@{}l}{continued on next page \\dots}\\\ \\bottomrule} \n \\tablelasttail{\\bottomrule} \n \\tablecaption[")
	file.write(shortcaption2)
	file.write("]{")
	file.write(caption2)
	file.write("}\label{")
	file.write(label2)
	file.write("}\n")
	file.write("\\begin{xtabular*}{\\textwidth}{@{\extracolsep{\\fill}}lrrrrr")
	if stableColumn:
		file.write("r")
	file.write("@{}} \n ")
	for j in range(ninstances[i]):
		file.write(names[i][j].replace("_", "\\_").split(".")[0] + "& ")
		if noslater[i][j] == "-":
			file.write("-- & ")
		else:
			file.write("\\num{" + "%.0f" % float(noslater[i][j]) + "} &")
		if noslaterfast[i][j] == "-":
			file.write("-- & ")
		else:
			file.write("\\num{" + "%.0f" % float(noslaterfast[i][j]) + "} &")
		if stableColumn:
			if noslaterstable[i][j] == "-":
				file.write("-- & ")
			else:
				file.write("\\num{" + "%.0f" % float(noslaterstable[i][j]) + "} &")
		if noslaterpenalty[i][j] == "-":
			file.write("-- & ")
		else:
			file.write("\\num{" + "%.0f" % float(noslaterpenalty[i][j]) + "} &")
		if noslaterbound[i][j] == "-":
			file.write("-- & ")
		else:
			file.write("\\num{" + "%.0f" % float(noslaterbound[i][j]) + "} &")
		if noslaterunsolved[i][j] == "-":
			file.write("-- \\\ \n")
		else:
			file.write("\\num{" + "%.0f" % float(noslaterunsolved[i][j]) + "} \\\ \n")
	file.write("  \\end{xtabular*} \n \\end{scriptsize} \n")

	file.write("\\newpage \n \\begin{scriptsize} \n  \\setlength{\\tabcolsep}{2pt} \n \\tablehead{\\toprule \n")
	file.write("problem & number & fast & ")
	if stableColumn:
		file.write("stable & ")
	file.write("penalty & bound & unssucc ")
	file.write("\\\ \\midrule} \n \\tabletail{ \\midrule \\multicolumn{3}{@{}l}{continued on next page \\dots}\\\ \\bottomrule} \n \\tablelasttail{\\bottomrule} \n \\tablecaption[")
	file.write(shortcaption3)
	file.write("]{")
	file.write(caption3)
	file.write("}\label{")
	file.write(label3)
	file.write("}\n")
	file.write("\\begin{xtabular*}{\\textwidth}{@{\extracolsep{\\fill}}lrrrrr")
	if stableColumn:
		file.write("r")
	file.write("@{}} \n ")
	for j in range(ninstances[i]):
		file.write(names[i][j].replace("_", "\\_").split(".")[0] + "& ")
		if infeasible[i][j] == "-":
			file.write("-- & ")
		else:
			file.write("\\num{" + "%.0f" % float(infeasible[i][j]) + "} &")
		if infeasiblefast[i][j] == "-":
			file.write("-- & ")
		else:
			file.write("\\num{" + "%.0f" % float(infeasiblefast[i][j]) + "} &")
		if stableColumn:
			if infeasiblestable[i][j] == "-":
				file.write("-- & ")
			else:
				file.write("\\num{" + "%.0f" % float(infeasiblestable[i][j]) + "} &")
		if infeasiblepenalty[i][j] == "-":
			file.write("-- & ")
		else:
			file.write("\\num{" + "%.0f" % float(infeasiblepenalty[i][j]) + "} &")
		if infeasiblebound[i][j] == "-":
			file.write("-- & ")
		else:
			file.write("\\num{" + "%.0f" % float(infeasiblebound[i][j]) + "} &")
		if infeasibleunsolved[i][j] == "-":
			file.write("-- \\\ \n")
		else:
			file.write("\\num{" + "%.0f" % float(infeasibleunsolved[i][j]) + "} \\\ \n")
	file.write("  \\end{xtabular*} \n \\end{scriptsize} \n")

def makeOverviewBySubsetSlater(file, args, settings, s, f, caption, label, settingnames):
	file.write("\\begin{table} \n \\begin{scriptsize} \\caption{" + caption + "} \n \\label{" + label + "} \n \\begin{tabular*}{\\textwidth}{@{}l@{\\;\\;\extracolsep{\\fill}}rrrrrrr")
	file.write("@{}}\\toprule \n")
	file.write("\phantom{abc} & \multicolumn{4}{c}{Dual Slater} & \multicolumn{3}{c}{Primal Slater} \\\ ")
	file.write("\\cmidrule(r){2-5} \\cmidrule(l){6-8} \n")
	file.write("problem & \\cmark & \\xmark & inf & ? & \\cmark & \\xmark & ? ")
	file.write("\\\ \midrule \n")
	i = 0
	ind = 0
	argind = -1
	for arg in args:
		if i in settings:
			unfailed = 0
			pslaterh = 0.0
			pslaterf = 0.0
			pslatern = 0.0
			dslaterh = 0.0
			dslaterf = 0.0
			dslateri = 0.0
			dslatern = 0.0
			j = s
			numinstances = f - s + 1
			while j <= f:
				if primalholds[i][j] != "-":
					unfailed += 1
				if primalholds[i][j] != "-":
					pslaterh += float(primalholds[i][j])
				if primalfails[i][j] != "-":
					pslaterf += float(primalfails[i][j])
				if primalnoinfo[i][j] != "-":
					pslatern += float(primalnoinfo[i][j])
				if dualholds[i][j] != "-":
					dslaterh += float(dualholds[i][j])
				if dualfails[i][j] != "-":
					dslaterf += float(dualfails[i][j])
				if dualinfeasible[i][j] != "-":
					dslateri += float(dualinfeasible[i][j])
				if dualnoinfo[i][j] != "-":
					dslatern += float(dualnoinfo[i][j])
				j = j + 1
			if unfailed:
				avgpslaterh = float(pslaterh) / float(unfailed)
			else:
				avgpslaterh = "-"
			if unfailed:
				avgpslaterf = float(pslaterf) / float(unfailed)
			else:
				avgpslaterf = "-"
			if unfailed:
				avgpslatern = float(pslatern) / float(unfailed)
			else:
				avgpslatern = "-"
			if unfailed:
				avgdslaterh = float(dslaterh) / float(unfailed)
			else:
				avgpslatern = "-"
			if unfailed:
				avgdslaterf = float(dslaterf) / float(unfailed)
			else:
				avgdslaterf = "-"
			if unfailed:
				avgdslateri = float(dslateri) / float(unfailed)
			else:
				avgdslateri = "-"
			if unfailed:
				avgdslatern = float(dslatern) / float(unfailed)
			else:
				avgdslatern = "-"
			file.write(settingnames[ind] + " & ")
			if avgdslaterh != "-":
				file.write("\\num{" + "%.2f" % avgdslaterh + "}\,\% & ")
			else:
				file.write("- & ")
			if avgdslaterf != "-":
				file.write("\\num{" + "%.2f" % avgdslaterf + "}\,\% & ")
			else:
				file.write("- & ")
			if avgdslateri != "-":
				file.write("\\num{" + "%.2f" % avgdslateri + "}\,\% & ")
			else:
				file.write("- & ")
			if avgdslatern != "-":
				file.write("\\num{" + "%.2f" % avgdslatern + "}\,\% &")
			else:
				file.write("- &")
			if avgpslaterh != "-":
				file.write("\\num{" + "%.2f" % avgpslaterh + "}\,\% & ")
			else:
				file.write("- & ")
			if avgpslaterf != "-":
				file.write("\\num{" + "%.2f" % avgpslaterf + "}\,\% & ")
			else:
				file.write("- & ")
			if avgpslatern != "-":
				file.write("\\num{" + "%.2f" % avgpslatern + "}\,\% ")
			else:
				file.write("- ")
			file.write(" \\\ \n")
			ind = ind + 1
		i = i + 1
	file.write("\\bottomrule \n \\end{tabular*} \n \end{scriptsize} \n \\end{table} \n")

def makeOverviewBySubsetSolved(file, args, settings, s, f, caption, label, settingnames):
	#file.write("\\begin{table} \n \\begin{scriptsize} \\caption{" + caption + "} \n \\label{" + label + "} \n \\begin{tabular*}{\\textwidth}{@{}l@{\\;\\;\extracolsep{\\fill}}rrrrrrrrrr")
	file.write("\\begin{table} \n \\begin{scriptsize} \\caption{" + caption + "} \n \\label{" + label + "} \n \\begin{tabular*}{\\textwidth}{@{}l@{\\;\\;\extracolsep{\\fill}}rrrrrrrrrrrr")
	file.write("@{}}\\toprule \n")
	#file.write("\phantom{abc} & \multicolumn{5}{c}{with Slater} & \multicolumn{5}{c}{without Slater} \\\ ")
	#file.write("\\cmidrule{2-6} \\cmidrule{7-11} \n")
	#file.write("problem & fast & stable & penalty & bound & unsolved & fast & stable & penalty & bound & unsolved ")
	file.write("\phantom{abc} & \multicolumn{4}{c}{with Slater} & \multicolumn{4}{c}{without Slater} & \multicolumn{4}{c}{infeasible} \\\ ")
	file.write("\\cmidrule(r){2-5} \\cmidrule{6-9} \\cmidrule(l){10-13} \n")
	file.write("problem & default & penalty & bound & unsolved & default & penalty & bound & unsolved & default & penalty & bound & unsolved ")
	file.write("\\\ \midrule \n")
	i = 0
	ind = 0
	argind = -1
	for arg in args:
		if i in settings:
			nslater = 0
			slaterf = 0
			slaters = 0
			slaterp = 0
			slaterb = 0
			slateru = 0
			nnoslater = 0
			noslaterf = 0
			noslaters = 0
			noslaterp = 0
			noslaterb = 0
			noslateru = 0
			ninf = 0
			inff = 0
			infs = 0
			infp = 0
			infb = 0
			infu = 0
			j = s
			numinstances = f - s + 1
			while j <= f:
				if slaterfast[i][j] != "-":
					slaterf += int(slaterfast[i][j])
					nslater += int(slaterfast[i][j])
				if slaterstable[i][j] != "-":
					slaters += int(slaterstable[i][j])
					nslater += int(slaterstable[i][j])
				if slaterpenalty[i][j] != "-":
					slaterp += int(slaterpenalty[i][j])
					nslater += int(slaterpenalty[i][j])
				if slaterbound[i][j] != "-":
					slaterb += int(slaterbound[i][j])
					nslater += int(slaterbound[i][j])
				if slaterunsolved[i][j] != "-":
					slateru += int(slaterunsolved[i][j])
					nslater += int(slaterunsolved[i][j])
				if noslaterfast[i][j] != "-":
					noslaterf += int(noslaterfast[i][j])
					nnoslater += int(noslaterfast[i][j])
				if noslaterstable[i][j] != "-":
					noslaters += int(noslaterstable[i][j])
					nnoslater += int(noslaterstable[i][j])
				if noslaterpenalty[i][j] != "-":
					noslaterp += int(noslaterpenalty[i][j])
					nnoslater += int(noslaterpenalty[i][j])
				if noslaterbound[i][j] != "-":
					noslaterb += int(noslaterbound[i][j])
					nnoslater += int(noslaterbound[i][j])
				if noslaterunsolved[i][j] != "-":
					noslateru += int(noslaterunsolved[i][j])
					nnoslater += int(noslaterunsolved[i][j])
				if infeasiblefast[i][j] != "-":
					inff += int(infeasiblefast[i][j])
					ninf += int(infeasiblefast[i][j])
				if infeasiblestable[i][j] != "-":
					infs += int(infeasiblestable[i][j])
					ninf += int(infeasiblestable[i][j])
				if infeasiblepenalty[i][j] != "-":
					infp += int(infeasiblepenalty[i][j])
					ninf += int(infeasiblepenalty[i][j])
				if infeasiblebound[i][j] != "-":
					infb += int(infeasiblebound[i][j])
					ninf += int(infeasiblebound[i][j])
				if infeasibleunsolved[i][j] != "-":
					infu += int(infeasibleunsolved[i][j])
					ninf += int(infeasibleunsolved[i][j])
				j = j + 1
			if nslater:
				avgslaterf = float(slaterf) / float(nslater)
			else:
				avgslaterf = "-"
			if nslater:
				avgslaters = float(slaters) / float(nslater)
			else:
				avgslaters = "-"
			if nslater:
				avgslaterp = float(slaterp) / float(nslater)
			else:
				avgslaterp = "-"
			if nslater:
				avgslaterb = float(slaterb) / float(nslater)
			else:
				avgslaterb = "-"
			if nslater:
				avgslateru = float(slateru) / float(nslater)
			else:
				avgslateru = "-"
			if nnoslater:
				avgnoslaterf = float(noslaterf) / float(nnoslater)
			else:
				avgnoslaterf = "-"
			if nnoslater:
				avgnoslaters = float(noslaters) / float(nnoslater)
			else:
				avgnoslaters = "-"
			if nnoslater:
				avgnoslaterp = float(noslaterp) / float(nnoslater)
			else:
				avgnoslaterp = "-"
			if nnoslater:
				avgnoslaterb = float(noslaterb) / float(nnoslater)
			else:
				avgnoslaterb = "-"
			if nnoslater:
				avgnoslateru = float(noslateru) / float(nnoslater)
			else:
				avgnoslateru = "-"
			if ninf:
				avginff = float(inff) / float(ninf)
			else:
				avginff = "-"
			if ninf:
				avginfs = float(infs) / float(ninf)
			else:
				avginfs = "-"
			if ninf:
				avginfp = float(infp) / float(ninf)
			else:
				avginfp = "-"
			if ninf:
				avginfb = float(infb) / float(ninf)
			else:
				avginfb = "-"
			if ninf:
				avginfu = float(infu) / float(ninf)
			else:
				avginfu = "-"
			file.write(settingnames[ind] + " & ")
			#if avgslaterf != "-":
			#	file.write("\\num{" + "%.2f" % avgslaterf + "} & ")
			#else:
			#	file.write("- & ")
			#if avgslaters != "-":
			#	file.write("\\num{" + "%.2f" % avgslaters + "} & ")
			#else:
			#	file.write("- & ")
			if avgslaterf != "-":
				file.write("\\num{" + "%.2f" % float(avgslaterf + avgslaters) + "}\,\% & ")
			else:
				file.write("- & ")
			if avgslaterp != "-":
				file.write("\\num{" + "%.2f" % avgslaterp + "}\,\% & ")
			else:
				file.write("- & ")
			if avgslaterb != "-":
				file.write("\\num{" + "%.2f" % avgslaterb + "}\,\% & ")
			else:
				file.write("- & ")
			if avgslateru != "-":
				file.write("\\num{" + "%.2f" % avgslateru + "}\,\% & ")
			else:
				file.write("- & ")
			if avgnoslaterf != "-":
				file.write("\\num{" + "%.2f" % float(avgnoslaterf + avgnoslaters) + "}\,\% & ")
			else:
				file.write("- & ")
			if avgnoslaterp != "-":
				file.write("\\num{" + "%.2f" % avgnoslaterp + "}\,\% & ")
			else:
				file.write("- & ")
			if avgnoslaterb != "-":
				file.write("\\num{" + "%.2f" % avgnoslaterb + "}\,\% & ")
			else:
				file.write("- & ")
			if avgnoslateru != "-":
				file.write("\\num{" + "%.2f" % avgnoslateru + "}\,\% &")
			else:
				file.write("- &")
			if avginff != "-":
				file.write("\\num{" + "%.2f" % float(avginff + avginfs) + "}\,\% & ")
			else:
				file.write("- & ")
			if avginfp != "-":
				file.write("\\num{" + "%.2f" % avginfp + "}\,\% & ")
			else:
				file.write("- & ")
			if avginfb != "-":
				file.write("\\num{" + "%.2f" % avginfb + "}\,\% & ")
			else:
				file.write("- & ")
			if avginfu != "-":
				file.write("\\num{" + "%.2f" % avginfu + "}\,\%")
			else:
				file.write("-")
			file.write(" \\\ \n")
			ind = ind + 1
		i = i + 1
	file.write("\\bottomrule \n \\end{tabular*} \n \end{scriptsize} \n \\end{table} \n")


def makeOverviewBySubsetSolvedInfeasible(file, args, settings, s, f, caption1, caption2, caption3, label1, label2, label3, settingnames):
	#file.write("\\begin{table} \n \\begin{scriptsize} \\caption{" + caption + "} \n \\label{" + label + "} \n \\begin{tabular*}{\\textwidth}{@{}l@{\\;\\;\extracolsep{\\fill}}rrrrrrrrrr")
	file.write("\\begin{table} \n \\begin{scriptsize} \\caption{" + caption1 + "} \n \\label{" + label1 + "} \n \\begin{tabular*}{\\textwidth}{@{}l@{\\;\\;\extracolsep{\\fill}}rrrrr")
	file.write("@{}}\\toprule \n")
	#file.write("\phantom{abc} & \multicolumn{4}{c}{with Slater} \\\ ")
	#file.write("\\cmidrule{2-5}  \n")
	file.write("settings & number & default & penalty & bound & unsucc ")
	file.write("\\\ \midrule \n")
	i = 0
	ind = 0
	argind = -1
	for arg in args:
		if i in settings:
			nslater = 0
			slaterf = 0
			slaters = 0
			slaterp = 0
			slaterb = 0
			slateru = 0
			j = s
			numinstances = f - s + 1
			while j <= f:
				if slater[i][j] != "-":
					nslater += int(slater[i][j])
				if slaterfast[i][j] != "-":
					slaterf += int(slaterfast[i][j])
				if slaterstable[i][j] != "-":
					slaters += int(slaterstable[i][j])
				if slaterpenalty[i][j] != "-":
					slaterp += int(slaterpenalty[i][j])
				if slaterbound[i][j] != "-":
					slaterb += int(slaterbound[i][j])
				if slaterunsolved[i][j] != "-":
					slateru += int(slaterunsolved[i][j])
				j = j + 1
			if nslater:
				avgslaterf = 100 * float(slaterf) / float(nslater)
			else:
				avgslaterf = "-"
			if nslater:
				avgslaters = 100 * float(slaters) / float(nslater)
			else:
				avgslaters = "-"
			if nslater:
				avgslaterp = 100 * float(slaterp) / float(nslater)
			else:
				avgslaterp = "-"
			if nslater:
				avgslaterb = 100 * float(slaterb) / float(nslater)
			else:
				avgslaterb = "-"
			if nslater:
				avgslateru = 100 * float(slateru) / float(nslater)
			else:
				avgslateru = "-"
			file.write(settingnames[ind] + " & ")
			file.write("\\num{" + "%.0f" % float(nslater) + "} & ")
			if avgslaterf != "-":
				file.write("\\num{" + "%.2f" % float(avgslaterf + avgslaters) + "}\,\% & ")
			else:
				file.write("- & ")
			if avgslaterp != "-":
				file.write("\\num{" + "%.2f" % avgslaterp + "}\,\% & ")
			else:
				file.write("- & ")
			if avgslaterb != "-":
				file.write("\\num{" + "%.2f" % avgslaterb + "}\,\% & ")
			else:
				file.write("- & ")
			if avgslateru != "-":
				file.write("\\num{" + "%.2f" % avgslateru + "}\,\% ")
			else:
				file.write("- ")
			file.write(" \\\ \n")
			ind = ind + 1
		i = i + 1
	file.write("\\bottomrule \n \\end{tabular*} \n \end{scriptsize} \n \\end{table} \n")
	#noslater
	file.write("\\begin{table} \n \\begin{scriptsize} \\caption{" + caption2 + "} \n \\label{" + label2 + "} \n \\begin{tabular*}{\\textwidth}{@{}l@{\\;\\;\extracolsep{\\fill}}rrrrr")
	file.write("@{}}\\toprule \n")
	#file.write("\phantom{abc} &  \multicolumn{4}{c}{without Slater} \\\ ")
	#file.write("\\cmidrule{2-5}   \n")
	file.write("settings & number & default & penalty & bound & unsucc ")
	file.write("\\\ \midrule \n")
	i = 0
	ind = 0
	argind = -1
	for arg in args:
		if i in settings:
			nnoslater = 0
			noslaterf = 0
			noslaters = 0
			noslaterp = 0
			noslaterb = 0
			noslateru = 0
			j = s
			numinstances = f - s + 1
			while j <= f:
				if noslater[i][j] != "-":
					nnoslater += int(noslater[i][j])
				if noslaterfast[i][j] != "-":
					noslaterf += int(noslaterfast[i][j])
				if noslaterstable[i][j] != "-":
					noslaters += int(noslaterstable[i][j])
				if noslaterpenalty[i][j] != "-":
					noslaterp += int(noslaterpenalty[i][j])
				if noslaterbound[i][j] != "-":
					noslaterb += int(noslaterbound[i][j])
				if noslaterunsolved[i][j] != "-":
					noslateru += int(noslaterunsolved[i][j])
				j = j + 1
			if nnoslater:
				avgnoslaterf = 100 * float(noslaterf) / float(nnoslater)
			else:
				avgnoslaterf = "-"
			if nnoslater:
				avgnoslaters = 100 * float(noslaters) / float(nnoslater)
			else:
				avgnoslaters = "-"
			if nnoslater:
				avgnoslaterp = 100 * float(noslaterp) / float(nnoslater)
			else:
				avgnoslaterp = "-"
			if nnoslater:
				avgnoslaterb = 100 * float(noslaterb) / float(nnoslater)
			else:
				avgnoslaterb = "-"
			if nnoslater:
				avgnoslateru = 100 * float(noslateru) / float(nnoslater)
			else:
				avgnoslateru = "-"
			file.write(settingnames[ind] + " & ")
			file.write("\\num{" + "%.0f" % float(nnoslater) + "} & ")
			if avgnoslaterf != "-":
				file.write("\\num{" + "%.2f" % float(avgnoslaterf + avgnoslaters) + "}\,\% & ")
			else:
				file.write("- & ")
			if avgnoslaterp != "-":
				file.write("\\num{" + "%.2f" % avgnoslaterp + "}\,\% & ")
			else:
				file.write("- & ")
			if avgnoslaterb != "-":
				file.write("\\num{" + "%.2f" % avgnoslaterb + "}\,\% & ")
			else:
				file.write("- & ")
			if avgnoslateru != "-":
				file.write("\\num{" + "%.2f" % avgnoslateru + "}\,\%")
			else:
				file.write("-")
			file.write(" \\\ \n")
			ind = ind + 1
		i = i + 1
	file.write("\\bottomrule \n \\end{tabular*} \n \end{scriptsize} \n \\end{table} \n")
	#infeasible
	file.write("\\begin{table} \n \\begin{scriptsize} \\caption{" + caption3 + "} \n \\label{" + label3 + "} \n \\begin{tabular*}{\\textwidth}{@{}l@{\\;\\;\extracolsep{\\fill}}rrrrr")
	file.write("@{}}\\toprule \n")
	#file.write("\phantom{abc} & \multicolumn{4}{c}{infeasible} \\\ ")
	#file.write("\\cmidrule{2-5}   \n")
	file.write("settings & number & default & penalty & bound & unsucc")
	file.write("\\\ \midrule \n")
	i = 0
	ind = 0
	argind = -1
	for arg in args:
		if i in settings:
			ninf = 0
			inff = 0
			infs = 0
			infp = 0
			infb = 0
			infu = 0
			j = s
			numinstances = f - s + 1
			while j <= f:
				if infeasible[i][j] != "-":
					ninf += int(infeasible[i][j])
				if infeasiblefast[i][j] != "-":
					inff += int(infeasiblefast[i][j])
				if infeasiblestable[i][j] != "-":
					infs += int(infeasiblestable[i][j])
				if infeasiblepenalty[i][j] != "-":
					infp += int(infeasiblepenalty[i][j])
				if infeasiblebound[i][j] != "-":
					infb += int(infeasiblebound[i][j])
				if infeasibleunsolved[i][j] != "-":
					infu += int(infeasibleunsolved[i][j])
				j = j + 1
			if ninf:
				avginff = 100 * float(inff) / float(ninf)
			else:
				avginff = "-"
			if ninf:
				avginfs = 100 * float(infs) / float(ninf)
			else:
				avginfs = "-"
			if ninf:
				avginfp = 100 * float(infp) / float(ninf)
			else:
				avginfp = "-"
			if ninf:
				avginfb = 100 * float(infb) / float(ninf)
			else:
				avginfb = "-"
			if ninf:
				avginfu = 100 * float(infu) / float(ninf)
			else:
				avginfu = "-"
			file.write(settingnames[ind] + " & ")
			file.write("\\num{" + "%.0f" % float(ninf) + "} & ")
			if avginff != "-":
				file.write("\\num{" + "%.2f" % float(avginff + avginfs) + "}\,\% & ")
			else:
				file.write("- & ")
			if avginfp != "-":
				file.write("\\num{" + "%.2f" % avginfp + "}\,\% & ")
			else:
				file.write("- & ")
			if avginfb != "-":
				file.write("\\num{" + "%.2f" % avginfb + "}\,\% & ")
			else:
				file.write("- & ")
			if avginfu != "-":
				file.write("\\num{" + "%.2f" % avginfu + "}\,\%")
			else:
				file.write("-")
			file.write(" \\\ \n")
			ind = ind + 1
		i = i + 1
	file.write("\\bottomrule \n \\end{tabular*} \n \end{scriptsize} \n \\end{table} \n")


if __name__=="__main__":
	"""give any number of .out-files for the same testset, then loops over them and returns a .tex-file given as first argument with some tables and a performance graph """
	ninstances=[0 for x in range(len(sys.argv) -2)] #initialize ninstances matrix
	names=[[0 for x in range(maxinstances)] for x in range(len(sys.argv) -2)] #initialize names matrix
	primalholds=[[0 for x in range(maxinstances)] for x in range(len(sys.argv) -2)] #initialize  matrix
	primalholdspercent=[[0 for x in range(maxinstances)] for x in range(len(sys.argv) -2)] #initialize  matrix
	primalfails=[[0 for x in range(maxinstances)] for x in range(len(sys.argv) -2)] #initialize  matrix
	primalfailspercent=[[0 for x in range(maxinstances)] for x in range(len(sys.argv) -2)] #initialize  matrix
	primalnoinfo=[[0 for x in range(maxinstances)] for x in range(len(sys.argv) -2)] #initialize  matrix
	primalnoinfopercent=[[0 for x in range(maxinstances)] for x in range(len(sys.argv) -2)] #initialize  matrix
	dualholds=[[0 for x in range(maxinstances)] for x in range(len(sys.argv) -2)] #initialize  matrix
	dualholdspercent=[[0 for x in range(maxinstances)] for x in range(len(sys.argv) -2)] #initialize  matrix
	dualfails=[[0 for x in range(maxinstances)] for x in range(len(sys.argv) -2)] #initialize  matrix
	dualfailspercent=[[0 for x in range(maxinstances)] for x in range(len(sys.argv) -2)] #initialize  matrix
	dualinfeasible=[[0 for x in range(maxinstances)] for x in range(len(sys.argv) -2)] #initialize  matrix
	dualinfeasiblepercent=[[0 for x in range(maxinstances)] for x in range(len(sys.argv) -2)] #initialize  matrix
	dualnoinfo=[[0 for x in range(maxinstances)] for x in range(len(sys.argv) -2)] #initialize  matrix
	dualnoinfopercent=[[0 for x in range(maxinstances)] for x in range(len(sys.argv) -2)] #initialize  matrix
	slater=[[0 for x in range(maxinstances)] for x in range(len(sys.argv) -2)] #initialize  matrix
	slaterfast=[[0 for x in range(maxinstances)] for x in range(len(sys.argv) -2)] #initialize  matrix
	slaterstable=[[0 for x in range(maxinstances)] for x in range(len(sys.argv) -2)] #initialize  matrix
	slaterpenalty=[[0 for x in range(maxinstances)] for x in range(len(sys.argv) -2)] #initialize  matrix
	slaterbound=[[0 for x in range(maxinstances)] for x in range(len(sys.argv) -2)] #initialize  matrix
	slaterunsolved=[[0 for x in range(maxinstances)] for x in range(len(sys.argv) -2)] #initialize  matrix
	noslater=[[0 for x in range(maxinstances)] for x in range(len(sys.argv) -2)] #initialize  matrix
	noslaterfast=[[0 for x in range(maxinstances)] for x in range(len(sys.argv) -2)] #initialize  matrix
	noslaterstable=[[0 for x in range(maxinstances)] for x in range(len(sys.argv) -2)] #initialize  matrix
	noslaterpenalty=[[0 for x in range(maxinstances)] for x in range(len(sys.argv) -2)] #initialize  matrix
	noslaterbound=[[0 for x in range(maxinstances)] for x in range(len(sys.argv) -2)] #initialize  matrix
	noslaterunsolved=[[0 for x in range(maxinstances)] for x in range(len(sys.argv) -2)] #initialize  matrix
	infeasible=[[0 for x in range(maxinstances)] for x in range(len(sys.argv) -2)] #initialize  matrix
	infeasiblefast=[[0 for x in range(maxinstances)] for x in range(len(sys.argv) -2)] #initialize  matrix
	infeasiblestable=[[0 for x in range(maxinstances)] for x in range(len(sys.argv) -2)] #initialize  matrix
	infeasiblepenalty=[[0 for x in range(maxinstances)] for x in range(len(sys.argv) -2)] #initialize  matrix
	infeasiblebound=[[0 for x in range(maxinstances)] for x in range(len(sys.argv) -2)] #initialize  matrix
	infeasibleunsolved=[[0 for x in range(maxinstances)] for x in range(len(sys.argv) -2)] #initialize  matrix
	sys.argv.remove("readAndPrintSlater_absolute.py")				 #remove function call from list of files
	texfilename = sys.argv[0]
	file=open(texfilename, "w")
	if texfile:
		file.write("\\documentclass[landscape]{article} \n \\usepackage{amsmath} \n \\usepackage{amsthm} \n \\usepackage{dsfont} \n \\usepackage[dvipsnames]{xcolor} \n \\usepackage{booktabs} \n \\usepackage{multirow} \n \\usepackage{mathtools} \n \\usepackage{xtab} \n \usepackage{tikz} \n \usepackage{pgfplots} \n  \usepackage[margin=1in,footskip=0.25in]{geometry} \n \\extrafloats{100} \n \\usepackage{sistyle} \n \\usepackage{tabularx} \n \\usepackage{pifont} \n \\newcommand{\\setting}[1]{\\texttt{#1}} \n   \\newcommand{\\xmark}{\\ding{55}} \n \\newcommand{\\cmark}{\\ding{51}} \n \\begin{document} \n")
	sys.argv.remove(sys.argv[0])
	i=0
	for arg in sys.argv:
		readFile(arg,i)
		print(str(i + 1) + "/" + str(len(sys.argv)) + " files read")
		i=i+1



	#Complete Table with all results for all settings
	if completeTableSlater:
		#for i in completetablesettings:
			makeCompleteTableSlater(completeTableSlaterShortCaptions[2], completeTableSlaterCaptions[2], completeTableSlaterLabels[2], file, 2)

	#Complete Table with all results for all settings
	if completeTableSolved:
		for i in completetablesettings:
			if not i == 2:
				stableColumn = False
			else:
				stableColumn = True
			makeCompleteTableSolved(completeTableSolvedHoldsShortCaptions[i], completeTableSolvedHoldsCaptions[i], completeTableSolvedHoldsLabels[i], completeTableSolvedFailsShortCaptions[i], completeTableSolvedFailsCaptions[i], completeTableSolvedFailsLabels[i], completeTableSolvedInfeasibleShortCaptions[i], completeTableSolvedInfeasibleCaptions[i], completeTableSolvedInfeasibleLabels[i], file, i, stableColumn)

	if SlaterTableCLS:
		makeOverviewBySubsetSlater(file, sys.argv, [0,1,2,3,4,5], 0, 64, "Statistics of Slater condition for the \emph{cardinality constrained least squares} testset of 65 instances", "CLSslater", ["\setting{DSDP-nodive}", "\setting{DSDP-frac10-fix}", "\setting{DSDP-nodive-rand10-fix}", "\setting{SDPA-nodive}", "\setting{SDPA-frac10-fix}", "\setting{SDPA-nodive-rand10-fix}"])

	if SlaterTableMinK:
		makeOverviewBySubsetSlater(file, sys.argv, [0,1,2,3,4,5], 65, 133, "Statistics of Slater condition for the \emph{min-$k$-partitioning} testset of 69 instances", "MinKslater", ["\setting{DSDP-nodive}", "\setting{DSDP-frac10-fix}", "\setting{DSDP-nodive-rand10-fix}", "\setting{SDPA-nodive}", "\setting{SDPA-frac10-fix}", "\setting{SDPA-nodive-rand10-fix}"])

	if SlaterTableTT:
		makeOverviewBySubsetSlater(file, sys.argv, [0,1,2,3,4,5], 134, 193, "Statistics of Slater condition for the \emph{truss topology} testset of 60 instances", "TTslater", ["\setting{DSDP-nodive}", "\setting{DSDP-frac10-fix}", "\setting{DSDP-nodive-rand10-fix}", "\setting{SDPA-nodive}", "\setting{SDPA-frac10-fix}", "\setting{SDPA-nodive-rand10-fix}"])

	if SlaterTableOverall:
		makeOverviewBySubsetSlater(file, sys.argv, [0,1,2,3,4,5], 0, 193, "Statistics of Slater condition for the \emph{complete testset} of 194 instances", "Overallslater", ["\setting{DSDP-nodive}", "\setting{DSDP-frac10-fix}", "\setting{DSDP-nodive-rand10-fix}", "\setting{SDPA-nodive}", "\setting{SDPA-frac10-fix}", "\setting{SDPA-nodive-rand10-fix}"])

	if SlaterSolvedTableCLS:
		#makeOverviewBySubsetSolved(file, sys.argv, [0,1,2,3,4,5], 0, 64, "Statistics of solver fails depending on Slater condition for the \emph{cardinality constrained least squares} testset of 65 instances", "CLSslaterSolved", ["\setting{DSDP-nodive}", "\setting{DSDP-frac10-fix}", "\setting{DSDP-nodive-rand10-fix}", "\setting{SDPA-nodive}", "\setting{SDPA-frac10-fix}", "\setting{SDPA-nodive-rand10-fix}"])
		makeOverviewBySubsetSolvedInfeasible(file, sys.argv, [0,1,2,3,4,5], 0, 64, "Statistics of solver fails when Slater condition holds for the \emph{cardinality constrained least squares} testset of 65 instances", "Statistics of solver fails when Slater condition does not hold for the \emph{cardinality constrained least squares} testset of 65 instances", "Statistics of solver fails for infeasible subproblems for the \emph{cardinality constrained least squares} testset of 65 instances", "CLSslaterHoldsSolved", "CLSslaterFailsSolved", "CLSslaterInfSolved", ["\setting{DSDP-nodive}", "\setting{DSDP-frac10-fix}", "\setting{DSDP-nodive-rand10-fix}", "\setting{SDPA-nodive}", "\setting{SDPA-frac10-fix}", "\setting{SDPA-nodive-rand10-fix}"])

	if SlaterSolvedTableMinK:
		#makeOverviewBySubsetSolved(file, sys.argv, [0,1,2,3,4,5], 65, 133, "Statistics of solver fails depending on Slater condition for the \emph{min-$k$-partitioning} testset of 69 instances", "MinKslaterSolved", ["\setting{DSDP-nodive}", "\setting{DSDP-frac10-fix}", "\setting{DSDP-nodive-rand10-fix}", "\setting{SDPA-nodive}", "\setting{SDPA-frac10-fix}", "\setting{SDPA-nodive-rand10-fix}"])
		makeOverviewBySubsetSolvedInfeasible(file, sys.argv, [0,1,2,3,4,5], 65, 133, "Statistics of solver fails when Slater condition holds for the \emph{min-$k$-partitioning} testset of 69 instances", "Statistics of solver fails when Slater condition does not hold for the \emph{min-$k$-partitioning} testset of 69 instances", "Statistics of solver fails for infeasible subproblems for the \emph{min-$k$-partitioning} testset of 69 instances", "MinKslaterHoldsSolved", "MinKslaterFailsSolved", "MinKslaterInfSolved", ["\setting{DSDP-nodive}", "\setting{DSDP-frac10-fix}", "\setting{DSDP-nodive-rand10-fix}", "\setting{SDPA-nodive}", "\setting{SDPA-frac10-fix}", "\setting{SDPA-nodive-rand10-fix}"])

	if SlaterSolvedTableTT:
		#makeOverviewBySubsetSolved(file, sys.argv, [0,1,2,3,4,5], 134, 193, "Statistics of solver fails depending on Slater condition for the \emph{truss topology} testset of 60 instances", "TTslaterSolved", ["\setting{DSDP-nodive}", "\setting{DSDP-frac10-fix}", "\setting{DSDP-nodive-rand10-fix}", "\setting{SDPA-nodive}", "\setting{SDPA-frac10-fix}", "\setting{SDPA-nodive-rand10-fix}"])
		makeOverviewBySubsetSolvedInfeasible(file, sys.argv, [0,1,2,3,4,5], 134, 193, "Statistics of solver fails when Slater condition holds for the \emph{truss topology} testset of 60 instances", "Statistics of solver fails when Slater condition does not hold for the \emph{truss topology} testset of 60 instances", "Statistics of solver fails for infeasible subproblems for the \emph{truss topology} testset of 60 instances", "TTslaterHoldsSolved", "TTslaterFailsSolved", "TTslaterInfSolved", ["\setting{DSDP-nodive}", "\setting{DSDP-frac10-fix}", "\setting{DSDP-nodive-rand10-fix}", "\setting{SDPA-nodive}", "\setting{SDPA-frac10-fix}", "\setting{SDPA-nodive-rand10-fix}"])

	if SlaterSolvedTableOverall:
		#makeOverviewBySubsetSolved(file, sys.argv, [0,1,2,3,4,5], 0, 193, "Statistics of solver fails depending on Slater condition for the \emph{complete testset} of 194 instances", "OverallslaterSolved", ["\setting{DSDP-nodive}", "\setting{DSDP-frac10-fix}", "\setting{DSDP-nodive-rand10-fix}", "\setting{SDPA-nodive}", "\setting{SDPA-frac10-fix}", "\setting{SDPA-nodive-rand10-fix}"])
		makeOverviewBySubsetSolvedInfeasible(file, sys.argv, [0,1,2,3,4,5], 0, 193, "Statistics of solver fails when the primal and dual Slater condition holds for the \emph{complete testset} of 194 instances", "Statistics of solver fails when either the primal or dual Slater condition does not hold for the \emph{complete testset} of 194 instances", "Statistics of solver fails for infeasible subproblems for the \emph{complete testset} of 194 instances", "OverallslaterHoldsSolved", "OverallslaterFailsSolved", "OverallslaterInfSolved", ["\setting{DSDP-nodive}", "\setting{DSDP-frac10-fix}", "\setting{DSDP-nodive-rand10-fix}", "\setting{SDPA-nodive}", "\setting{SDPA-frac10-fix}", "\setting{SDPA-nodive-rand10-fix}"])


	if SlaterTableShortCLS:
		makeOverviewBySubsetSlater(file, sys.argv, [0,1,2], 0, 64, "Statistics of Slater condition for the \emph{cardinality constrained least squares} testset of 65 instances", "CLSslater", ["\setting{DSDP}", "\setting{SDPA}", "\setting{MOSEK}"])

	if SlaterTableShortMinK:
		makeOverviewBySubsetSlater(file, sys.argv, [0,1,2], 65, 133, "Statistics of Slater condition for the \emph{min-$k$-partitioning} testset of 69 instances", "MinKslater", ["\setting{DSDP}", "\setting{SDPA}", "\setting{MOSEK}"])

	if SlaterTableShortTT:
		makeOverviewBySubsetSlater(file, sys.argv, [0,1,2], 134, 193, "Statistics of Slater condition for the \emph{truss topology} testset of 60 instances", "TTslater", ["\setting{DSDP}", "\setting{SDPA}", "\setting{MOSEK}"])

	if SlaterTableShortOverall:
		makeOverviewBySubsetSlater(file, sys.argv, [0,1,2], 0, 193, "Statistics of Slater condition for the \emph{complete testset} of 194 instances", "Overallslater", ["\setting{DSDP}", "\setting{SDPA}", "\setting{MOSEK}"])

	if SlaterSolvedTableShortCLS:
		#makeOverviewBySubsetSolved(file, sys.argv, [0,1,2,3,4,5], 0, 64, "Statistics of solver fails depending on Slater condition for the \emph{cardinality constrained least squares} testset of 65 instances", "CLSslaterSolved", ["\setting{DSDP-nodive}", "\setting{DSDP-frac10-fix}", "\setting{DSDP-nodive-rand10-fix}", "\setting{SDPA-nodive}", "\setting{SDPA-frac10-fix}", "\setting{SDPA-nodive-rand10-fix}"])
		makeOverviewBySubsetSolvedInfeasible(file, sys.argv, [0,1,2], 0, 64, "Statistics of solver fails when Slater condition holds for the \emph{cardinality constrained least squares} testset of 65 instances", "Statistics of solver fails when Slater condition does not hold for the \emph{cardinality constrained least squares} testset of 65 instances", "Statistics of solver fails for infeasible subproblems for the \emph{cardinality constrained least squares} testset of 65 instances", "CLSslaterHoldsSolved", "CLSslaterFailsSolved", "CLSslaterInfSolved", ["\setting{DSDP}", "\setting{SDPA}", "\setting{MOSEK}"])

	if SlaterSolvedTableShortMinK:
		#makeOverviewBySubsetSolved(file, sys.argv, [0,1,2,3,4,5], 65, 133, "Statistics of solver fails depending on Slater condition for the \emph{min-$k$-partitioning} testset of 69 instances", "MinKslaterSolved", ["\setting{DSDP-nodive}", "\setting{DSDP-frac10-fix}", "\setting{DSDP-nodive-rand10-fix}", "\setting{SDPA-nodive}", "\setting{SDPA-frac10-fix}", "\setting{SDPA-nodive-rand10-fix}"])
		makeOverviewBySubsetSolvedInfeasible(file, sys.argv, [0,1,2], 65, 133, "Statistics of solver fails when Slater condition holds for the \emph{min-$k$-partitioning} testset of 69 instances", "Statistics of solver fails when Slater condition does not hold for the \emph{min-$k$-partitioning} testset of 69 instances", "Statistics of solver fails for infeasible subproblems for the \emph{min-$k$-partitioning} testset of 69 instances", "MinKslaterHoldsSolved", "MinKslaterFailsSolved", "MinKslaterInfSolved", ["\setting{DSDP}", "\setting{SDPA}", "\setting{MOSEK}"])

	if SlaterSolvedTableShortTT:
		#makeOverviewBySubsetSolved(file, sys.argv, [0,1,2,3,4,5], 134, 193, "Statistics of solver fails depending on Slater condition for the \emph{truss topology} testset of 60 instances", "TTslaterSolved", ["\setting{DSDP-nodive}", "\setting{DSDP-frac10-fix}", "\setting{DSDP-nodive-rand10-fix}", "\setting{SDPA-nodive}", "\setting{SDPA-frac10-fix}", "\setting{SDPA-nodive-rand10-fix}"])
		makeOverviewBySubsetSolvedInfeasible(file, sys.argv, [0,1,2], 134, 193, "Statistics of solver fails when Slater condition holds for the \emph{truss topology} testset of 60 instances", "Statistics of solver fails when Slater condition does not hold for the \emph{truss topology} testset of 60 instances", "Statistics of solver fails for infeasible subproblems for the \emph{truss topology} testset of 60 instances", "TTslaterHoldsSolved", "TTslaterFailsSolved", "TTslaterInfSolved", ["\setting{DSDP}", "\setting{SDPA}", "\setting{MOSEK}"])

	if SlaterSolvedTableShortOverall:
		#makeOverviewBySubsetSolved(file, sys.argv, [0,1,2,3,4,5], 0, 193, "Statistics of solver fails depending on Slater condition for the \emph{complete testset} of 194 instances", "OverallslaterSolved", ["\setting{DSDP-nodive}", "\setting{DSDP-frac10-fix}", "\setting{DSDP-nodive-rand10-fix}", "\setting{SDPA-nodive}", "\setting{SDPA-frac10-fix}", "\setting{SDPA-nodive-rand10-fix}"])
		makeOverviewBySubsetSolvedInfeasible(file, sys.argv, [0,1,2], 0, 193, "Statistics of solver fails when the primal and dual Slater condition holds for the \emph{complete testset} of 194 instances", "Statistics of solver fails when either the primal or dual Slater condition does not hold for the \emph{complete testset} of 194 instances", "Statistics of solver fails for infeasible subproblems for the \emph{complete testset} of 194 instances", "OverallslaterHoldsSolved", "OverallslaterFailsSolved", "OverallslaterInfSolved", ["\setting{DSDP}", "\setting{SDPA}", "\setting{MOSEK}"])


	if texfile:
		file.write("\\end{document}")
	file.close()
