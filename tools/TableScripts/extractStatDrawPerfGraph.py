import os
import sys
import math
from decimal import Decimal
#usage: python extractStatDrawPerfGraph.py textfile.tex outfile1.out outfile2.out

#example: python extractStatDrawPerfGraph.py /home/gally/gally_Dissertation/resultfiles/WarmstartResults20180711/warmstartcomparison_SCIPSDPRIP_wosimplepreopt001.tex /home/gally/gally_Dissertation/resultfiles/WarmstartResults20180711/check.SCIPSDPRIP.scipsdp.linux.x86_64.gnu.opt.sdpa.moskito.disable_warmstart.out /home/gally/gally_Dissertation/resultfiles/WarmstartResults20180711/check.SCIPSDPRIP.scipsdp.linux.x86_64.gnu.opt.sdpa.moskito.warmstartpreoptgap05.out /home/gally/gally_Dissertation/resultfiles/WarmstartResults20180711/check.SCIPSDPRIP.scipsdp.linux.x86_64.gnu.opt.sdpa.moskito.warmstartipfactor001_proj2_pdsame.out /home/gally/gally_Dissertation/resultfiles/WarmstartResults20180711/check.SCIPSDPRIP.scipsdp.linux.x86_64.gnu.opt.sdpa.moskito.warmstartipfactor05_proj2_pddiff.out /home/gally/gally_Dissertation/resultfiles/WarmstartResults20180711/check.SCIPSDPRIP.scipsdp.linux.x86_64.gnu.opt.sdpa.moskito.warmstartipfactor05_proj2_pdsame.out /home/gally/gally_Dissertation/resultfiles/WarmstartResults20180711/check.SCIPSDPRIP.scipsdp.linux.x86_64.gnu.opt.sdpa.moskito.warmstart_analcent_05_proj2.out /home/gally/gally_Dissertation/resultfiles/WarmstartResults20180711/check.SCIPSDPRIP.scipsdp.linux.x86_64.gnu.opt.sdpa.moskito.warmstartprojminev-1pdsame.out /home/gally/gally_Dissertation/resultfiles/WarmstartResults20180711/check.SCIPSDPRIP.scipsdp.linux.x86_64.gnu.opt.sdpa.moskito.warmstartipfactor05_proj4.out /home/gally/gally_Dissertation/resultfiles/WarmstartResults20180711/check.SCIPSDPRIP.scipsdp.linux.x86_64.gnu.opt.sdpa.moskito.warmstartroundinfonly.out



#example: python extractStatDrawPerfGraph.py /home/gally/gally_Dissertation/resultfiles/WarmstartResults20180711/warmstartcomparison_SCIPSDPRIPpaper.tex /home/gally/gally_Dissertation/resultfiles/WarmstartResults20180711/check.SCIPSDPRIP.scipsdp.linux.x86_64.gnu.opt.sdpa.moskito.disable_warmstart.out /home/gally/gally_Dissertation/resultfiles/WarmstartResults20180711/check.SCIPSDPRIP.scipsdp.linux.x86_64.gnu.opt.sdpa.moskito.warmstart.out /home/gally/gally_Dissertation/resultfiles/WarmstartResults20180711/check.SCIPSDPRIP.scipsdp.linux.x86_64.gnu.opt.sdpa.moskito.warmstartpreoptgap001.out /home/gally/gally_Dissertation/resultfiles/WarmstartResults20180711/check.SCIPSDPRIP.scipsdp.linux.x86_64.gnu.opt.sdpa.moskito.warmstartpreoptgap05.out /home/gally/gally_Dissertation/resultfiles/WarmstartResults20180711/check.SCIPSDPRIP.scipsdp.linux.x86_64.gnu.opt.sdpa.moskito.warmstartipfactor001_proj2_pdsame.out /home/gally/gally_Dissertation/resultfiles/WarmstartResults20180711/check.SCIPSDPRIP.scipsdp.linux.x86_64.gnu.opt.sdpa.moskito.warmstartipfactor05_proj2_pddiff.out /home/gally/gally_Dissertation/resultfiles/WarmstartResults20180711/check.SCIPSDPRIP.scipsdp.linux.x86_64.gnu.opt.sdpa.moskito.warmstartipfactor05_proj2_pdsame.out /home/gally/gally_Dissertation/resultfiles/WarmstartResults20180711/check.SCIPSDPRIP.scipsdp.linux.x86_64.gnu.opt.sdpa.moskito.warmstart_analcent_05_proj2.out /home/gally/gally_Dissertation/resultfiles/WarmstartResults20180711/check.SCIPSDPRIP.scipsdp.linux.x86_64.gnu.opt.sdpa.moskito.warmstartprojminev-1pdsame.out /home/gally/gally_Dissertation/resultfiles/WarmstartResults20180711/check.SCIPSDPRIP.scipsdp.linux.x86_64.gnu.opt.sdpa.moskito.warmstartipfactor05_proj4.out /home/gally/gally_Dissertation/resultfiles/WarmstartResults20180711/check.SCIPSDPRIP.scipsdp.linux.x86_64.gnu.opt.sdpa.moskito.warmstartroundinfonly.out

#example: python extractStatDrawPerfGraph.py /local/gally/results/WarmstartResults20171123/warmstartcomparison_SCIPSDPRIPpaper.tex /local/gally/results/WarmstartResults20171123/check.SCIPSDPRIPpaper.scipsdp.linux.x86_64.gnu.opt.sdpa.moskito.disable_warmstart.out /local/gally/results/WarmstartResults20171123/check.SCIPSDPRIPpaper.scipsdp.linux.x86_64.gnu.opt.sdpa.moskito.warmstart.out /local/gally/results/WarmstartResults20171123/check.SCIPSDPRIPpaper.scipsdp.linux.x86_64.gnu.opt.sdpa.moskito.warmstartpreoptgap05.out /local/gally/results/WarmstartResults20171123/check.SCIPSDPRIPpaper.scipsdp.linux.x86_64.gnu.opt.sdpa.moskito.warmstartpreoptgap01.out /local/gally/results/WarmstartResults20171123/check.SCIPSDPRIPpaper.scipsdp.linux.x86_64.gnu.opt.sdpa.moskito.warmstartpreoptgap001.out /local/gally/results/WarmstartResults20171123/check.SCIPSDPRIPpaper.scipsdp.linux.x86_64.gnu.opt.sdpa.moskito.warmstartpreoptgap0001.out /local/gally/results/WarmstartResults20171123/check.SCIPSDPRIPpaper.scipsdp.linux.x86_64.gnu.opt.sdpa.moskito.warmstartipfactor001_proj2_pddiff.out /local/gally/results/WarmstartResults20171123/check.SCIPSDPRIPpaper.scipsdp.linux.x86_64.gnu.opt.sdpa.moskito.warmstartipfactor001_proj2_pdsame.out /local/gally/results/WarmstartResults20171123/check.SCIPSDPRIPpaper.scipsdp.linux.x86_64.gnu.opt.sdpa.moskito.warmstartipfactor05_proj2_pddiff.out /local/gally/results/WarmstartResults20171123/check.SCIPSDPRIPpaper.scipsdp.linux.x86_64.gnu.opt.sdpa.moskito.warmstartipfactor05_proj2_pdsame.out /local/gally/results/WarmstartResults20171123/check.SCIPSDPRIPpaper.scipsdp.linux.x86_64.gnu.opt.sdpa.moskito.warmstart_analcent_05_proj2.out /local/gally/results/WarmstartResults20171123/check.SCIPSDPRIPpaper.scipsdp.linux.x86_64.gnu.opt.sdpa.moskito.warmstart_analcent_001_proj2.out /local/gally/results/WarmstartResults20171123/check.SCIPSDPRIPpaper.scipsdp.linux.x86_64.gnu.opt.sdpa.moskito.warmstartprojminev10.out /local/gally/results/WarmstartResults20171123/check.SCIPSDPRIPpaper.scipsdp.linux.x86_64.gnu.opt.sdpa.moskito.warmstartprojminev-1pddifferent.out /local/gally/results/WarmstartResults20171123/check.SCIPSDPRIPpaper.scipsdp.linux.x86_64.gnu.opt.sdpa.moskito.warmstartprojminev-1pdsame.out /local/gally/results/WarmstartResults20171123/check.SCIPSDPRIPpaper.scipsdp.linux.x86_64.gnu.opt.sdpa.moskito.warmstartipfactor001_proj4.out /local/gally/results/WarmstartResults20171123/check.SCIPSDPRIPpaper.scipsdp.linux.x86_64.gnu.opt.sdpa.moskito.warmstartipfactor05_proj4.out /local/gally/results/WarmstartResults20171123/check.SCIPSDPRIPpaper.scipsdp.linux.x86_64.gnu.opt.sdpa.moskito.warmstartroundinfonly.out

#example: python extractStatDrawPerfGraph.py /home/gally/Documents/KonferenzenVortraegePoster/UGworkshop2018/SCIPSDPcomparison.tex /home/gally/gally_Dissertation/resultfiles/Parallel/SCIPSDP/check.SCIPSDPpaper.scipsdp.linux.x86_64.gnu.opt.msk.fb04549.1thread.out /home/gally/gally_Dissertation/resultfiles/MISDPvergleich/check.SCIPSDPpaper.scipsdp.linux.x86_64.gnu.opt.msk.fb04668.lp_approx.out

#example: python extractStatDrawPerfGraph.py /home/gally/gally_Dissertation/Supplement/Tables/Heuristics.tex /home/gally/gally_Dissertation/resultfiles/Heuristiken/check.SCIPSDPpaper.scipsdp.linux.x86_64.gnu.opt.msk.moskito.no_rand.out /home/gally/gally_Dissertation/resultfiles/Heuristiken/check.SCIPSDPpaper.scipsdp.linux.x86_64.gnu.opt.msk.moskito.rand_10.out /home/gally/gally_Dissertation/resultfiles/Solververgleich/check.SCIPSDPpaper.scipsdp.linux.x86_64.gnu.opt.msk.moskito.1thread.out /home/gally/gally_Dissertation/resultfiles/Heuristiken/check.SCIPSDPpaper.scipsdp.linux.x86_64.gnu.opt.msk.moskito.rand_10rounds.out /home/gally/gally_Dissertation/resultfiles/Heuristiken/check.SCIPSDPpaper.scipsdp.linux.x86_64.gnu.opt.msk.moskito.no_rand_fracdive_root.out /home/gally/gally_Dissertation/resultfiles/Heuristiken/check.SCIPSDPpaper.scipsdp.linux.x86_64.gnu.opt.msk.moskito.no_rand_fracdive_10.out /home/gally/gally_Dissertation/resultfiles/Heuristiken/check.SCIPSDPpaper.scipsdp.linux.x86_64.gnu.opt.msk.moskito.fracdive_root.out /home/gally/gally_Dissertation/resultfiles/Heuristiken/check.SCIPSDPpaper.scipsdp.linux.x86_64.gnu.opt.msk.moskito.fracdive_10.out

#example python extractStatDrawPerfGraph.py /home/gally/gally_Dissertation/resultfiles/Solververgleich/Solververgleich.tex /home/gally/gally_Dissertation/resultfiles/Solververgleich/check.SCIPSDPpaper.scipsdp.linux.x86_64.gnu.opt.dsdp.moskito.1thread.out /home/gally/gally_Dissertation/resultfiles/Solververgleich/check.SCIPSDPpaper.scipsdp.linux.x86_64.gnu.opt.sdpa.moskito.disable_warmstart.out /home/gally/gally_Dissertation/resultfiles/Solververgleich/check.SCIPSDPpaper.scipsdp.linux.x86_64.gnu.opt.msk.moskito.1thread.out

#example python extractStatDrawPerfGraph.py /home/gally/gally_Dissertation/resultfiles/Branchingrules/BranchRulesvergleich.tex /home/gally/gally_Dissertation/resultfiles/Branchingrules/check.SCIPSDPpaper.scipsdp.linux.x86_64.gnu.opt.msk.moskito.branchmostinf.out /home/gally/gally_Dissertation/resultfiles/Branchingrules/check.SCIPSDPpaper.scipsdp.linux.x86_64.gnu.opt.msk.moskito.branchmostfrac.out /home/gally/gally_Dissertation/resultfiles/Branchingrules/check.SCIPSDPpaper.scipsdp.linux.x86_64.gnu.opt.msk.moskito.branchobj.out /home/gally/gally_Dissertation/resultfiles/Solververgleich/check.SCIPSDPpaper.scipsdp.linux.x86_64.gnu.opt.msk.moskito.1thread.out

#example python extractStatDrawPerfGraph.py /home/gally/gally_Dissertation/resultfiles/obbt/obbtVergleichFracdive.tex /home/gally/gally_Dissertation/resultfiles/Solververgleich/check.SCIPSDPpaper.scipsdp.linux.x86_64.gnu.opt.msk.moskito.1thread.out /home/gally/gally_Dissertation/resultfiles/obbt/check.SCIPSDPpaper.scipsdp.linux.x86_64.gnu.opt.msk.moskito.1thread_obbt_bin.out /home/gally/gally_Dissertation/resultfiles/Heuristiken/check.SCIPSDPpaper.scipsdp.linux.x86_64.gnu.opt.msk.moskito.fracdive_root.out /home/gally/gally_Dissertation/resultfiles/obbt/check.SCIPSDPpaper.scipsdp.linux.x86_64.gnu.opt.msk.moskito.1thread_obbt_bin_fracdive.out

#example python extractStatDrawPerfGraph.py /home/gally/gally_Dissertation/resultfiles/DualFixing/DualFixingVergleichDomReds.tex /home/gally/gally_Dissertation/resultfiles/Solververgleich/check.SCIPSDPpaper.scipsdp.linux.x86_64.gnu.opt.msk.moskito.1thread.out /home/gally/gally_Dissertation/resultfiles/DualFixing/check.SCIPSDPpaper.scipsdp.linux.x86_64.gnu.opt.msk.moskito.1thread_nofixing.out

#example python extractStatDrawPerfGraph.py /home/gally/gally_Dissertation/Supplement/Tables/DSDPSDPA.tex /home/gally/gally_Dissertation/resultfiles/Solververgleich/check.SCIPSDPpaper.scipsdp.linux.x86_64.gnu.opt.dsdp.moskito.1thread.out /home/gally/gally_Dissertation/resultfiles/Solververgleich/check.SCIPSDPpaper.scipsdp.linux.x86_64.gnu.opt.sdpa.moskito.disable_warmstart.out

#example python extractStatDrawPerfGraph.py /home/gally/gally_Dissertation/Supplement/Tables/OBBT.tex /home/gally/gally_Dissertation/resultfiles/obbt/check.SCIPSDPpaper.scipsdp.linux.x86_64.gnu.opt.msk.moskito.1thread_obbt_bin.out /home/gally/gally_Dissertation/resultfiles/obbt/check.SCIPSDPpaper.scipsdp.linux.x86_64.gnu.opt.msk.moskito.1thread_obbt_bin_fracdive.out

#example python extractStatDrawPerfGraph.py /home/gally/gally_Dissertation/Supplement/Tables/branching.tex /home/gally/gally_Dissertation/resultfiles/Branchingrules/check.SCIPSDPpaper.scipsdp.linux.x86_64.gnu.opt.msk.moskito.branchmostinf.out /home/gally/gally_Dissertation/resultfiles/Branchingrules/check.SCIPSDPpaper.scipsdp.linux.x86_64.gnu.opt.msk.moskito.branchmostfrac.out /home/gally/gally_Dissertation/resultfiles/Branchingrules/check.SCIPSDPpaper.scipsdp.linux.x86_64.gnu.opt.msk.moskito.branchobj.out

#example python extractStatDrawPerfGraph.py /home/gally/gally_Dissertation/Supplement/Tables/warmstarts_roundingtime.tex /home/gally/gally_Dissertation/resultfiles/WarmstartResults20180711/check.SCIPSDPRIP.scipsdp.linux.x86_64.gnu.opt.sdpa.moskito.warmstart.out /home/gally/gally_Dissertation/resultfiles/WarmstartResults20180711/check.SCIPSDPRIP.scipsdp.linux.x86_64.gnu.opt.sdpa.moskito.warmstartpreoptgap001.out /home/gally/gally_Dissertation/resultfiles/WarmstartResults20180711/check.SCIPSDPRIP.scipsdp.linux.x86_64.gnu.opt.sdpa.moskito.warmstartpreoptgap05.out /home/gally/gally_Dissertation/resultfiles/WarmstartResults20180711/check.SCIPSDPRIP.scipsdp.linux.x86_64.gnu.opt.sdpa.moskito.warmstartipfactor001_proj2_pdsame.out /home/gally/gally_Dissertation/resultfiles/WarmstartResults20180711/check.SCIPSDPRIP.scipsdp.linux.x86_64.gnu.opt.sdpa.moskito.warmstartipfactor05_proj2_pddiff.out /home/gally/gally_Dissertation/resultfiles/WarmstartResults20180711/check.SCIPSDPRIP.scipsdp.linux.x86_64.gnu.opt.sdpa.moskito.warmstartipfactor05_proj2_pdsame.out /home/gally/gally_Dissertation/resultfiles/WarmstartResults20180711/check.SCIPSDPRIP.scipsdp.linux.x86_64.gnu.opt.sdpa.moskito.warmstart_analcent_05_proj2.out /home/gally/gally_Dissertation/resultfiles/WarmstartResults20180711/check.SCIPSDPRIP.scipsdp.linux.x86_64.gnu.opt.sdpa.moskito.warmstartprojminev-1pdsame.out /home/gally/gally_Dissertation/resultfiles/WarmstartResults20180711/check.SCIPSDPRIP.scipsdp.linux.x86_64.gnu.opt.sdpa.moskito.warmstartipfactor05_proj4.out /home/gally/gally_Dissertation/resultfiles/WarmstartResults20180711/check.SCIPSDPRIP.scipsdp.linux.x86_64.gnu.opt.sdpa.moskito.warmstartroundinfonly.out

#example python extractStatDrawPerfGraph.py /home/gally/gally_Dissertation/Supplement/Tables/mosek_parallel.tex /home/gally/gally_Dissertation/resultfiles/Parallel/SCIPSDP/check.SCIPSDPpaper.scipsdp.linux.x86_64.gnu.opt.msk.fb04549.1thread.out /home/gally/gally_Dissertation/resultfiles/Parallel/SCIPSDP/check.SCIPSDPpaper.scipsdp.linux.x86_64.gnu.opt.msk.fb04549.2thread.out /home/gally/gally_Dissertation/resultfiles/Parallel/SCIPSDP/check.SCIPSDPpaper.scipsdp.linux.x86_64.gnu.opt.msk.fb04549.4thread.out

#TODO: irgendwas klappt mit unsolved nicht, Gesamtdurchschnitt ist viel hoeher als pro Testset, gefuehlt cardleastsquares zu hoch, rest zu gering
maxinstances = 321 # maximum number of instances
completeTable = 0 # make table listing all instances
overviewTable = 0 # make table for each setting listing all subsets
overviewBySubset = 1 # make table for each subset listing all settings
overviewTotal = 1 # make table for complete testset listing all settings
performancePlotOverview = 0 # draw a performance plot for the complete testset
performancePlotBySubset = 0 # draw a performance plot for each subset
heurcolumns = 0 #add columns for fracdivefoundbest, randfoundbest and primaldualintegral
redcostcolumn = 0
roundingtimecolumn = 0
obbtcolumn = 0 #add column for number of domain reductions of sdp-obbt
bestrow = 1 #add additional row for optimal auto settings (currenlty only supports nsolved, time and nodes)
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
fixingsshift = 1
heurshift = 1
heurarithmetic = 0 #TODO: also do this for percentages or at least count number of solved instances
fixingsarithmetic = 0
itergeomonlysolved = 1
nodegeomonlysolved = 1 # in the tables for each setting only include solved instances into the geometric mean of the nodes, for the overviews only include instances which were solved by all settings
abortscolumn = 0 # should an additional column for number of aborts be added to the overview by subsets tables
primaldualintegralTmax = 1 #should the primal dual integral be nomralized with Tmax? (otherwise it is normalized with the solving time of each instance)
Tmax = 3600
primaldualintegralarithmetic = 0 #should arithmetic mean be used for primaldualintegral?
primaldualintegralshift = 1
performanceprofileNormed = 0 #should the performance plots be normalized by the fastest solver like in Dolan & More to get the relative time or should we just plot over the absolute time instead?
overviewlongtable = 0 # should overview tables be done as longtables
plotsbysettingssubsets = 1 # make seperate plots and tables for each subset in settingssubsets
#settingsets = [[0,1,2,3,4,5], [6,7,8,9,10,11], [12,13,14,15,16,17], [18,19,20,21,22,23], [24,25,26,27,28,29], [30,31,32,33,34,35], [36,37,38,39,40,41], [42,43,44,45,46,47], [0,6,12,18,24,30,36,42], [1,7,13,19,25,31,37,43], [2,8,14,20,26,32,38,44], [3,9,15,21,27,33,39,45], [4,10,16,22,28,34,40,46], [5,11,17,23,29,35,41,47]]
#settingsets = [[0,2,3,4,5,11,17,23], [17,18,19,20,21,22], [23,24,25,26,27,28]]
#settingsets = [[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17]]
#settingsets = [[0,1,2,3,4,5,6,7,8]]
#settingsets = [[0,1,2,3,4,5,6,7,8,9,10]]
#settingsets = [[0,1,2,3,4,5,6,7,8,9]]
#settingsets = [[0,1,2,3]]
#settingsets = [[0,1,2,3,4,5,6,7]]
settingsets = [[0,1]]
usesettingnames = 1
#settingnames = ["\setting{noSDPheur}","\setting{rand\\_10\\_2\\_no\\_dive}","\setting{rand\\_1\\_2\\_no\\_dive}","\setting{rand\\_1\\_10\\_no\\_dive}", "\setting{no\\_rand\\_frac\\_root}", "\setting{no\\_rand\\_frac\\_10}", "\setting{rand\\_1\\_2\\_frac\\_root}", "\setting{rand\\_1\\_2\\_frac\\_10}"]
#settingnames = ["DSDP", "SDPA", "MOSEK"]
#settingnames = ["\setting{no\_obbt}", "\setting{obbt}", "\setting{fracdive}", "\setting{obbt+fracdive}"]
settingnames = ["NL-BB", "LP-cuts"]
#settingnames = ["\setting{default}", "\setting{no dualfixing}"]
#settingnames = [ "no warmstart", "simple warmstart", "preoptgap 0.5", "preoptgap 0.1", "preoptgap 0.01", "preoptgap 0.001", "0.01 id pddiff", "0.01 id pdsame", "0.5 id pddiff", "0.5 id pdsame", "0.5 analcent", "0.01 analcent", "proj minev 10", "proj minev auto pddiff", "proj minev auto pdsame", "roundingprob 0.01 id", "roundingprob 0.5 id", "roundingprob inf only"]
#settingnames = [ "no warmstart", "simple warmstart", "preoptgap 0.01", "preoptgap 0.5", "0.01 id pdsame", "0.5 id pddiff", "0.5 id pdsame", "0.5 analcent", "proj minev auto pdsame", "roundingprob 0.5 id", "roundingprob inf only"]
#settingnames = [ "no warmstart", "preoptgap 0.5", "0.01 id pdsame", "0.5 id pddiff", "0.5 id pdsame", "0.5 analcent", "proj minev auto pdsame", "roundingprob 0.5 id", "roundingprob inf only"]
epsilon = 0.001
captions=["Complete results for \scipsdp with \mosek running on 1 thread on shared memory environment of Intel Xeon E5-4650 CPUs running at \SI{2.70}{GHz} with \SI{512}{GB} of shared RAM", "Complete results for \scipsdp with \mosek running on 2 threads on shared memory environment of Intel Xeon E5-4650 CPUs running at \SI{2.70}{GHz} with \SI{512}{GB} of shared RAM", "Complete results for \scipsdp with \mosek running on 4 threads on shared memory environment of Intel Xeon E5-4650 CPUs running at \SI{2.70}{GHz} with \SI{512}{GB} of shared RAM"]
shortcaptions=["Complete results for \scipsdp with \mosek running on 1 thread", "Complete results for \scipsdp with \mosek running on 2 threads", "Complete results for \scipsdp with \mosek running on 4 threads"]
labels=["MOSEK1thread", "MOSEK2threads", "MOSEK4threads"]


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
			aborted[i][j] = 1
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
		if substring.find("sdp-obbt         :") == -1:
			sdpobbtfixings[i][j] = "-"
		else:
			sdpobbtfixing = substring[substring.find("sdp-obbt         :") + 18:].split()[3]
			sdpobbtfixings[i][j] = int(sdpobbtfixing)
		if substring.find("sdpfracdiving    :") == -1:
			fracdivefoundbest[i][j] = "-"
		else:
			sdpfracdive = substring[substring.find("sdpfracdiving    :") + 18:].split()[4]
			fracdivefoundbest[i][j] = int(sdpfracdive)
		if substring.find("sdprand          :") == -1:
			randfoundbest[i][j] = "-"
		else:
			sdprand = substring[substring.find("sdprand          :") + 18:].split()[4]
			randfoundbest[i][j] = int(sdprand)
		if primaldualintegralTmax:
			if substring.find("Avg. Gap         :") == -1:
				primaldualintegral[i][j] = "-"
			else:
				pdintegral = substring[substring.find("Avg. Gap         :") + 18:].split()[2][1:len(substring[substring.find("Avg. Gap         :") + 18:].split()[2])]
				primaldualintegral[i][j] = float(pdintegral)/3600.0
		else:
			if substring.find("Avg. Gap         :") == -1:
				primaldualintegral[i][j] = "-"
			else:
				pdintegral = substring[substring.find("Avg. Gap         :") + 18:].split()[0]
				primaldualintegral[i][j] = float(pdintegral)
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
			if substring.find("Relaxators         :       ") and substring[substring.find("Relaxators         :       "):].find("  SDP              :") == -1:
				sdpiters[i][j] = "-"
			else:
				sdpline = substring[substring.find("Relaxators         :       "):][substring[substring.find("Relaxators         :       "):].find("  SDP              :"):]
				sdpiters[i][j] = sdpline.split()[4]
		else:
			sdpiter = substring[substring.find("SDP iterations:") + 15:].split()[0]
			sdpiters[i][j] = int(sdpiter)
		if substring.find("Percentage penalty formulation used:	 ") == -1:
			if substring.find("     MOSEK 8.1.0.54:   ") == -1:
				penalties[i][j] = "-"
			else:
					mosekline = substring[substring.find("     MOSEK 8.1.0.54:   "):]
					penalties[i][j] = mosekline.split()[4]
		else:
			penalty = substring[substring.find("Percentage penalty formulation used:") + 38:].split()[0]
			penalties[i][j] = float(penalty)
		if substring.find("Percentage unsolved even with penalty:	") == -1:
			if substring.find("     MOSEK 8.1.0.54:   ") == -1:
				unsolved[i][j] = "-"
			else:
					mosekline = substring[substring.find("     MOSEK 8.1.0.54:   "):]
					unsolved[i][j] = mosekline.split()[6]
		else:
			unsolv = substring[substring.find("Percentage unsolved even with penalty:") + 40:].split()[0]
			unsolved[i][j] = float(unsolv)
		if substring.find("Time spent in rounding problems for warmstarting / detecting infeasibility:") == -1:
			roundingtime[i][j] = "-"
		else:
			roundtime = substring[substring.find("Time spent in rounding problems for warmstarting / detecting infeasibility:") + 76:].split()[0]
			roundingtime[i][j] = float(roundtime)
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

def bestTime(ind, settings):
	""" returns best time of any solver in settings for instance ind"""
	best = Tmax
	for i in settings:
		if times[i][ind] < best:
			best = times[i][ind]
	return best

def bestNodes(ind, settings):
	""" returns smallest number of nodes of any solver in settings for instance ind"""
	best = float("inf")
	for i in settings:
		if nodes[i][ind] < best:
			best = nodes[i][ind]
	return best

def makeCompleteTableCaption(arg, i):
	file.write("\\newpage \n \\begin{scriptsize} \n  \\setlength{\\tabcolsep}{2pt} \n \\tablehead{\\toprule \n")
	file.write("problem & dbound &  pbound & gap & nodes & time & iters ")
	if redcostcolumn or obbtcolumn:
		file.write("& domreds ")
	if heurcolumns:
		file.write("& rrbf & fdbf & pdintegral ")
	else:
		file.write("& pen & uns")
	if roundingtimecolumn:
		file.write("& roundtime")
	file.write("\\\ \\midrule} \n \\tabletail{ \\midrule \\multicolumn{3}{@{}l}{continued on next page \\dots}\\\ \\bottomrule} \n \\tablelasttail{\\bottomrule} \n \\tablecaption[")
	file.write(shortcaptions[i])
	file.write("]{")
	file.write(captions[i])
	file.write("}\label{")
	file.write(labels[i])
	file.write("}\n")
	file.write("\\begin{xtabular*}{\\textwidth}{@{\extracolsep{\\fill}}lrrrrrrrr")
	if heurcolumns:
		file.write("rr")
	if redcostcolumn or obbtcolumn:
		file.write("r")
	if roundingtimecolumn:
		file.write("r")
	file.write("@{}} \n ")
	for j in range(ninstances[i]):
		file.write(names[i][j].replace("_", "\\_").split(".")[0] + "& ")
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
		if sdpiters[i][j] == "-":
			file.write("-- & ")
		else:
			file.write("\\num{" + "%.0f" % float(sdpiters[i][j]) + "} &")
		if roundingtimecolumn:
			if roundingtime[i][j] == "-":
				file.write("-- & ")
			else:
				file.write("\\num{" + "%.1f" % float(roundingtime[i][j]) + "} &")
		if redcostcolumn:
			if sdpredcostfixings[i][j] == "-":
				file.write("-- & ")
			else:
				file.write("\\num{" + "%.0f" % float(sdpredcostfixings[i][j]) + "} &")
		if obbtcolumn:
			if sdpobbtfixings[i][j] == "-":
				file.write("-- & ")
			else:
				file.write("\\num{" + "%.0f" % float(sdpobbtfixings[i][j]) + "} &")
		if heurcolumns:
			if randfoundbest[i][j] == "-":
				file.write("-- & ")
			else:
				file.write("\\num{" + "%.0f" % float(randfoundbest[i][j]) + "} &")
			if fracdivefoundbest[i][j] == "-":
				file.write("-- & ")
			else:
				file.write("\\num{" + "%.0f" % float(fracdivefoundbest[i][j]) + "} &")
			if primaldualintegral[i][j] == "-":
				file.write("-- & ")
			else:
				file.write("\\SI{" + "%.2f" % float(primaldualintegral[i][j]) + "}{\percent} &")
		else:
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

def makeCompleteTable(arg, file, i):
	if heurcolumns:
		file.write("\\begin{longtable}[ht] {p{.25\\textwidth} | p{.09\\textwidth} p{.09\\textwidth} p{.056\\textwidth} p{.05\\textwidth} p{.05\\textwidth} p{.05\\textwidth} p{.05\\textwidth} p{.05\\textwidth} p{.05\\textwidth}} \\toprule \n problem & dual bound &  primal bound & gap & nodes & fixings & fracdive & time & sdpiter & randbestfound & divebestfound & pdintegral \% \\\ \midrule \n")
		for j in range(ninstances[i]):
			file.write("$" + names[i][j] + "$ & $ %.4g" % dualresults[i][j] + "$ & $ %.4g" % primalresults[i][j] + "$ & $" + str(gaps[i][j]) + "$ & $" + str(nodes[i][j]) + "$ & $" + str(times[i][j]) + "$ & $" + str(sdpiters[i][j]) + "$ & $" + str(randfoundbest[i][j]) + "$ & $" + str(fracdivefoundbest[i][j]) + "$ & $" + str(primaldualintegral[i][j]) + "$ \\\ \n")
		file.write("\\bottomrule \caption{Results for $" + str((os.path.basename(arg)).split(".")[7]) + " | " + str((os.path.basename(arg)).split(".")[9]) + "$} \end{longtable} \n")
	else:
		file.write("\\begin{longtable}[ht] {p{.3\\textwidth} | p{.09\\textwidth} p{.09\\textwidth} p{.056\\textwidth} p{.05\\textwidth} p{.05\\textwidth} p{.05\\textwidth} p{.05\\textwidth} p{.05\\textwidth}} \\toprule \n problem & dual bound &  primal bound & gap & nodes & fixings & fracdive & time & sdpiter & penalty & unsolved \% \\\ \midrule \n")
		for j in range(ninstances[i]):
			file.write("$" + names[i][j] + "$ & $ %.4g" % dualresults[i][j] + "$ & $ %.4g" % primalresults[i][j] + "$ & $" + str(gaps[i][j]) + "$ & $" + str(nodes[i][j]) + "$ & $" + str(times[i][j]) + "$ & $" + str(sdpiters[i][j]) + "$ & $" + str(penalties[i][j]) + "$ & $" + str(unsolved[i][j]) + "$ \\\ \n")
		file.write("\\bottomrule \caption{Results for $" + str((os.path.basename(arg)).split(".")[7]) + " | " + str((os.path.basename(arg)).split(".")[9]) + "$} \end{longtable} \n")

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
			if heurarithmetic and fracdivefoundbest[i][j] != "-":
				totalheur = totalheur + float(fracdivefoundbest[i][j])
				nsolved = nsolved + 1
			elif fracdivefoundbest[i][j] != "-":
				totalheur = math.pow(totalheur, float(j-s) / float(j-s+1)) * math.pow(float(fracdivefoundbest[i][j] + heurshift), 1.0/float(j-s+1))
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
		if heurarithmetic and fracdivefoundbest[i][j] != "-":
			totalheur = totalheur + float(fracdivefoundbest[i][j])
			nsolved = nsolved + 1
		elif fracdivefoundbest[i][j] != "-":
			totalheur = math.pow(totalheur, float(j-s) / float(j-s+1)) * math.pow(float(fracdivefoundbest[i][j] + heurshift), 1.0/float(j-s+1))
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
        assert(not(redcostcolumn and obbtcolumn))
	if overviewlongtable:
		file.write("\\begin{longtable}[ht] {p{.3\\textwidth} | p{.05\\textwidth} p{.09\\textwidth} ")
		if fixingscol:
			file.write("p{.09\\textwidth} ")
		file.write("p{.05\\textwidth} p{.05\\textwidth} p{.09\\textwidth} p{.09\\textwidth} p{.05\\textwidth}} \\toprule \n settings &  solved & nodes & ")
		if fixingscol:
			file.write("fixings & ")
		file.write("fracdive & time & sdpiter & penalty & unsolved \% \\\ \midrule \n")
	else:
		if heurcolumns:
			file.write("\\begin{table} \n \\begin{scriptsize} \\caption{" + caption + "} \n \\label{" + label + "} \n \\begin{tabular*}{\\textwidth}{@{}l@{\\;\\;\extracolsep{\\fill}}rrrrrrr")
			if abortscolumn:
				file.write("r")
			file.write("@{}}\\toprule \n")
			file.write(" settings &  solved & ")
			if abortscolumn:
				file.write("aborts & ")
			file.write("nodes & ")
			file.write("time & sdpiter & randbestfound & divebestfound & pdintegral \\\ \midrule \n")
		else:
			file.write("\\begin{table} \n \\begin{scriptsize} \\caption{" + caption + "} \n \\label{" + label + "} \n \\begin{tabular*}{\\textwidth}{@{}l@{\\;\\;\extracolsep{\\fill}}rrrrrr")
			if redcostcolumn or obbtcolumn:
				file.write("r")
			if abortscolumn:
				file.write("r")
			file.write("@{}}\\toprule \n")
			file.write(" settings &  solved & ")
			if abortscolumn:
				file.write("aborts & ")
			file.write("nodes & ")
			if redcostcolumn or obbtcolumn:
				file.write("domreds & ")
			file.write("time & sdpiter & penalty & unsolved \\\ \midrule \n")
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
				totalredcostfixings = 1.0
                        totalobbtfixings = 1.0
			fracdivebestgeom = 1.0
			randbestgeom = 1.0
			if primaldualintegralarithmetic:
				primaldualintegralsum = 0
			else:
				primaldualintegralgeom = 1.0
			timegeom = 1.0
			itergeom = 1.0
			totalpenalty = 0
			totalunsolved = 0
                        nunaborted = 0
                        nnodescounted = 0
                        niterscounted = 0
			naborts = 0 #note that this and nunaborted could be merged
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
					totalredcostfixings = totalredcostfixings * math.pow(float(sdpredcostfixings[i][j] + fixingsshift), 1.0/float(nunaborted))
                                if sdpobbtfixings[i][j] != "-":
					totalobbtfixings = totalobbtfixings * math.pow(float(sdpobbtfixings[i][j] + fixingsshift), 1.0/float(nunaborted))
				if fracdivefoundbest[i][j] != "-":
					fracdivebestgeom = fracdivebestgeom * math.pow(float(fracdivefoundbest[i][j] + heurshift), 1.0/float(nunaborted))
				if randfoundbest[i][j] != "-":
					randbestgeom = randbestgeom * math.pow(float(randfoundbest[i][j] + heurshift), 1.0/float(nunaborted))
				if primaldualintegral[i][j] != "-":
					if primaldualintegralarithmetic:
						primaldualintegralsum += primaldualintegral[i][j];
					else:
						primaldualintegralgeom = primaldualintegralgeom * math.pow(float(primaldualintegral[i][j] + primaldualintegralshift), 1.0/float(nunaborted))
				timegeom = timegeom * math.pow(float(times[i][j] + timeshift), 1.0/float(f-s+1))
				if sdpiters[i][j] != "-" and (not itergeomonlysolved or allSolved(j, settings)):
					itergeom = itergeom * math.pow(float(sdpiters[i][j]) + itershift, 1.0/float(niterscounted))
				if penalties[i][j] != "-":
					totalpenalty = totalpenalty + float(penalties[i][j])
				if unsolved[i][j] != "-":
					totalunsolved = totalunsolved + float(unsolved[i][j])
				naborts += aborted[i][j]
				j = j + 1

			nodegeom = nodegeom - nodeshift
			if fixingsarithmetic:
				avgfixings = totalfixings / nunaborted
			else:
				avgredcostfixings = totalredcostfixings - fixingsshift
                        avgobbtfixings = totalobbtfixings - fixingsshift
			fracdivebestgeom = fracdivebestgeom - heurshift
			randbestgeom = randbestgeom - heurshift
			if primaldualintegralarithmetic:
				avgpdintegral = float(primaldualintegralsum) / float(nunaborted)
			else:
				avgpdintegral = primaldualintegralgeom - primaldualintegralshift
			timegeom = timegeom - timeshift
			itergeom = itergeom - itershift
			avgpenalty = float(totalpenalty) / float(nunaborted)
			avgunsolved = float(totalunsolved) / float(nunaborted)
			if heurcolumns:
				file.write(settingnames[ind] + " & \\num{" + "%.0f" % nsolved + "} & ")
				if abortscolumn:
					file.write("\\num{" + "%.0f" % naborts + "} & ")
				file.write("\\num{" + "%.2f" % nodegeom + "} & ")
				file.write("\\num{" + "%.2f" % timegeom + "} & \\num{" + "%.2f" % itergeom + "} & \\num{" + "%.2f" % randbestgeom + "} & \\num{" + "%.2f" % fracdivebestgeom + "} & \\num{" + "%.2f" % avgpdintegral + "} \% \\\ \n")
			else:
				file.write(settingnames[ind] + " & \\num{" + "%.0f" % nsolved + "} & ")
				if abortscolumn:
					file.write("\\num{" + "%.0f" % naborts + "} & ")
				file.write("\\num{" + "%.2f" % nodegeom + "} & ")
				if redcostcolumn:
					file.write("\\num{" + "%.2f" % avgredcostfixings + "} & ")
                                elif obbtcolumn:
                                        file.write("\\num{" + "%.2f" % avgobbtfixings + "} & ")
				file.write("\\num{" + "%.2f" % timegeom + "} & \\num{" + "%.2f" % itergeom + "} & \\num{" + "%.2f" % avgpenalty + "} \% & \\num{" + "%.2f" % avgunsolved + "} \% \\\ \n")
			ind = ind + 1
		i = i + 1
	if bestrow:
		nsolved = 0
		nodegeom = 1.0
		timegeom = 1.0
		nnodescounted = 0
		j = s
		numinstances = f - s + 1
		assert(nodegeomonlysolved) #otherwise while loop for finding number of instances to average over needs to be adjusted
		while j <= f:
 			if allSolved(j, settings):
				nnodescounted = nnodescounted + 1
			j = j + 1
		j = s
		while j <= f:
			if bestTime(j, settings) < Tmax :
				nsolved = nsolved + 1
			if allSolved(j, settings):
				nodegeom = nodegeom * math.pow(float(bestNodes(j, settings) + nodeshift), 1.0/float(nnodescounted))
			timegeom = timegeom * math.pow(float(bestTime(j, settings) + timeshift), 1.0/float(f-s+1))
			j = j + 1

		nodegeom = nodegeom - nodeshift
		timegeom = timegeom - timeshift
		if heurcolumns:
			file.write("best" + " & \\num{" + "%.0f" % nsolved + "} & ")
			if abortscolumn:
				file.write("-- & ")
			file.write("\\num{" + "%.2f" % nodegeom + "} & ")
			file.write("\\num{" + "%.2f" % timegeom + "} & -- & -- & -- & -- \\\ \n")
		else:
			file.write("best" + " & \\num{" + "%.0f" % nsolved + "} & ")
			if abortscolumn:
				file.write("-- & ")
			file.write("\\num{" + "%.2f" % nodegeom + "} & ")
			if redcostcolumn:
				file.write("-- & ")
			elif obbtcolumn:
				file.write("-- & ")
			file.write("\\num{" + "%.2f" % timegeom + "} & -- & -- & -- \\\ \n")
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
	names=[[0 for x in range(maxinstances)] for x in range(len(sys.argv) -2)] #initialize names matrix
	dualresults=[[0 for x in range(maxinstances)] for x in range(len(sys.argv) -2)] #initialize dual results matrix
	primalresults=[[0 for x in range(maxinstances)] for x in range(len(sys.argv) -2)] #initialize primal results matrix
	gaps=[[0 for x in range(maxinstances)] for x in range(len(sys.argv) -2)]   #initialize gaps matrix
	times=[[0 for x in range(maxinstances)] for x in range(len(sys.argv) -2)]   #initialize solvingtime matrix
	nodes=[[0 for x in range(maxinstances)] for x in range(len(sys.argv) -2)]   #initialize nodes matrix
	sdpredcostfixings=[[0 for x in range(maxinstances)] for x in range(len(sys.argv) -2)]   #initialize redcostfixings matrix
        sdpobbtfixings=[[0 for x in range(maxinstances)] for x in range(len(sys.argv) -2)]   #initialize obbtfixings matrix
	fracdivefoundbest=[[0 for x in range(maxinstances)] for x in range(len(sys.argv) -2)]   #initialize fracdivefoundbest matrix
	randfoundbest=[[0 for x in range(maxinstances)] for x in range(len(sys.argv) -2)]   #initialize randfoundbest matrix
	primaldualintegral=[[0 for x in range(maxinstances)] for x in range(len(sys.argv) -2)]   #initialize primaldualintegral matrix
	sdpiters=[[0 for x in range(maxinstances)] for x in range(len(sys.argv) -2)]   #initialize sdpiterations matrix
	penalties=[[0 for x in range(maxinstances)] for x in range(len(sys.argv) -2)]   #initialize penalty% matrix
	unsolved=[[0 for x in range(maxinstances)] for x in range(len(sys.argv) -2)]   #initialize unsolved% matrix
	aborted=[[0 for x in range(maxinstances)] for x in range(len(sys.argv) -2)]   #initialize unsolved% matrix
	roundingtime=[[0 for x in range(maxinstances)] for x in range(len(sys.argv) -2)]   #initialize roundingtime matrix
	sys.argv.remove("extractStatDrawPerfGraph.py")				 #remove function call from list of files
	texfilename = sys.argv[0]
	assert(texfilename.split(".")[1]=="tex")
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
			makeOverviewTable(sys.argv[i], i)
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
	file.write("\\end{document}")
	file.close()
