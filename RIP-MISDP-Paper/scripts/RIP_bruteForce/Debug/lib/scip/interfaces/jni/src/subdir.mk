################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
C_SRCS += \
../lib/scip/interfaces/jni/src/JniScip.c \
../lib/scip/interfaces/jni/src/JniScipBranch.c \
../lib/scip/interfaces/jni/src/JniScipBranchAllfullstrong.c \
../lib/scip/interfaces/jni/src/JniScipBranchFullstrong.c \
../lib/scip/interfaces/jni/src/JniScipBranchInference.c \
../lib/scip/interfaces/jni/src/JniScipBranchLeastinf.c \
../lib/scip/interfaces/jni/src/JniScipBranchMostinf.c \
../lib/scip/interfaces/jni/src/JniScipBranchPscost.c \
../lib/scip/interfaces/jni/src/JniScipBranchRandom.c \
../lib/scip/interfaces/jni/src/JniScipBranchRelpscost.c \
../lib/scip/interfaces/jni/src/JniScipConflict.c \
../lib/scip/interfaces/jni/src/JniScipCons.c \
../lib/scip/interfaces/jni/src/JniScipConsAbspower.c \
../lib/scip/interfaces/jni/src/JniScipConsAnd.c \
../lib/scip/interfaces/jni/src/JniScipConsBivariate.c \
../lib/scip/interfaces/jni/src/JniScipConsBounddisjunction.c \
../lib/scip/interfaces/jni/src/JniScipConsConjunction.c \
../lib/scip/interfaces/jni/src/JniScipConsCountsols.c \
../lib/scip/interfaces/jni/src/JniScipConsCumulative.c \
../lib/scip/interfaces/jni/src/JniScipConsDisjunction.c \
../lib/scip/interfaces/jni/src/JniScipConsIndicator.c \
../lib/scip/interfaces/jni/src/JniScipConsIntegral.c \
../lib/scip/interfaces/jni/src/JniScipConsKnapsack.c \
../lib/scip/interfaces/jni/src/JniScipConsLinear.c \
../lib/scip/interfaces/jni/src/JniScipConsLinking.c \
../lib/scip/interfaces/jni/src/JniScipConsLogicor.c \
../lib/scip/interfaces/jni/src/JniScipConsNonlinear.c \
../lib/scip/interfaces/jni/src/JniScipConsOr.c \
../lib/scip/interfaces/jni/src/JniScipConsOrbitope.c \
../lib/scip/interfaces/jni/src/JniScipConsPseudoboolean.c \
../lib/scip/interfaces/jni/src/JniScipConsQuadratic.c \
../lib/scip/interfaces/jni/src/JniScipConsSetppc.c \
../lib/scip/interfaces/jni/src/JniScipConsSoc.c \
../lib/scip/interfaces/jni/src/JniScipConsSos1.c \
../lib/scip/interfaces/jni/src/JniScipConsSos2.c \
../lib/scip/interfaces/jni/src/JniScipConsSuperindicator.c \
../lib/scip/interfaces/jni/src/JniScipConsVarbound.c \
../lib/scip/interfaces/jni/src/JniScipConsXor.c \
../lib/scip/interfaces/jni/src/JniScipCutpool.c \
../lib/scip/interfaces/jni/src/JniScipDialog.c \
../lib/scip/interfaces/jni/src/JniScipDialogDefault.c \
../lib/scip/interfaces/jni/src/JniScipDisp.c \
../lib/scip/interfaces/jni/src/JniScipDispDefault.c \
../lib/scip/interfaces/jni/src/JniScipErrorCode.c \
../lib/scip/interfaces/jni/src/JniScipEvent.c \
../lib/scip/interfaces/jni/src/JniScipEventTypes.c \
../lib/scip/interfaces/jni/src/JniScipEventhdlr.c \
../lib/scip/interfaces/jni/src/JniScipException.c \
../lib/scip/interfaces/jni/src/JniScipExpr.c \
../lib/scip/interfaces/jni/src/JniScipHeur.c \
../lib/scip/interfaces/jni/src/JniScipHeurActconsdiving.c \
../lib/scip/interfaces/jni/src/JniScipHeurClique.c \
../lib/scip/interfaces/jni/src/JniScipHeurCoefdiving.c \
../lib/scip/interfaces/jni/src/JniScipHeurCrossover.c \
../lib/scip/interfaces/jni/src/JniScipHeurDins.c \
../lib/scip/interfaces/jni/src/JniScipHeurFeaspump.c \
../lib/scip/interfaces/jni/src/JniScipHeurFixandinfer.c \
../lib/scip/interfaces/jni/src/JniScipHeurFracdiving.c \
../lib/scip/interfaces/jni/src/JniScipHeurGuideddiving.c \
../lib/scip/interfaces/jni/src/JniScipHeurIntdiving.c \
../lib/scip/interfaces/jni/src/JniScipHeurIntshifting.c \
../lib/scip/interfaces/jni/src/JniScipHeurLinesearchdiving.c \
../lib/scip/interfaces/jni/src/JniScipHeurLocalbranching.c \
../lib/scip/interfaces/jni/src/JniScipHeurMutation.c \
../lib/scip/interfaces/jni/src/JniScipHeurNlpdiving.c \
../lib/scip/interfaces/jni/src/JniScipHeurObjpscostdiving.c \
../lib/scip/interfaces/jni/src/JniScipHeurOctane.c \
../lib/scip/interfaces/jni/src/JniScipHeurOneopt.c \
../lib/scip/interfaces/jni/src/JniScipHeurPscostdiving.c \
../lib/scip/interfaces/jni/src/JniScipHeurRens.c \
../lib/scip/interfaces/jni/src/JniScipHeurRins.c \
../lib/scip/interfaces/jni/src/JniScipHeurRootsoldiving.c \
../lib/scip/interfaces/jni/src/JniScipHeurRounding.c \
../lib/scip/interfaces/jni/src/JniScipHeurShiftandpropagate.c \
../lib/scip/interfaces/jni/src/JniScipHeurShifting.c \
../lib/scip/interfaces/jni/src/JniScipHeurSimplerounding.c \
../lib/scip/interfaces/jni/src/JniScipHeurSubnlp.c \
../lib/scip/interfaces/jni/src/JniScipHeurTrivial.c \
../lib/scip/interfaces/jni/src/JniScipHeurTrysol.c \
../lib/scip/interfaces/jni/src/JniScipHeurTwoopt.c \
../lib/scip/interfaces/jni/src/JniScipHeurUndercover.c \
../lib/scip/interfaces/jni/src/JniScipHeurVbounds.c \
../lib/scip/interfaces/jni/src/JniScipHeurVeclendiving.c \
../lib/scip/interfaces/jni/src/JniScipHeurZeroobj.c \
../lib/scip/interfaces/jni/src/JniScipHeurZirounding.c \
../lib/scip/interfaces/jni/src/JniScipImplics.c \
../lib/scip/interfaces/jni/src/JniScipLibraryLoader.c \
../lib/scip/interfaces/jni/src/JniScipLp.c \
../lib/scip/interfaces/jni/src/JniScipMessage.c \
../lib/scip/interfaces/jni/src/JniScipMessageDefault.c \
../lib/scip/interfaces/jni/src/JniScipMisc.c \
../lib/scip/interfaces/jni/src/JniScipNlp.c \
../lib/scip/interfaces/jni/src/JniScipNlpiIpopt.c \
../lib/scip/interfaces/jni/src/JniScipNodesel.c \
../lib/scip/interfaces/jni/src/JniScipNodeselBfs.c \
../lib/scip/interfaces/jni/src/JniScipNodeselDfs.c \
../lib/scip/interfaces/jni/src/JniScipNodeselEstimate.c \
../lib/scip/interfaces/jni/src/JniScipNodeselHybridestim.c \
../lib/scip/interfaces/jni/src/JniScipNodeselRestartdfs.c \
../lib/scip/interfaces/jni/src/JniScipParamset.c \
../lib/scip/interfaces/jni/src/JniScipPresol.c \
../lib/scip/interfaces/jni/src/JniScipPresolBoundshift.c \
../lib/scip/interfaces/jni/src/JniScipPresolComponents.c \
../lib/scip/interfaces/jni/src/JniScipPresolConvertinttobin.c \
../lib/scip/interfaces/jni/src/JniScipPresolDomcol.c \
../lib/scip/interfaces/jni/src/JniScipPresolDualfix.c \
../lib/scip/interfaces/jni/src/JniScipPresolGateextraction.c \
../lib/scip/interfaces/jni/src/JniScipPresolImplics.c \
../lib/scip/interfaces/jni/src/JniScipPresolInttobinary.c \
../lib/scip/interfaces/jni/src/JniScipPresolTrivial.c \
../lib/scip/interfaces/jni/src/JniScipPricer.c \
../lib/scip/interfaces/jni/src/JniScipProp.c \
../lib/scip/interfaces/jni/src/JniScipPropDualfix.c \
../lib/scip/interfaces/jni/src/JniScipPropGenvbounds.c \
../lib/scip/interfaces/jni/src/JniScipPropObbt.c \
../lib/scip/interfaces/jni/src/JniScipPropProbing.c \
../lib/scip/interfaces/jni/src/JniScipPropPseudoobj.c \
../lib/scip/interfaces/jni/src/JniScipPropRedcost.c \
../lib/scip/interfaces/jni/src/JniScipPropRootredcost.c \
../lib/scip/interfaces/jni/src/JniScipPropVbounds.c \
../lib/scip/interfaces/jni/src/JniScipReader.c \
../lib/scip/interfaces/jni/src/JniScipReaderBnd.c \
../lib/scip/interfaces/jni/src/JniScipReaderCcg.c \
../lib/scip/interfaces/jni/src/JniScipReaderCip.c \
../lib/scip/interfaces/jni/src/JniScipReaderCnf.c \
../lib/scip/interfaces/jni/src/JniScipReaderFix.c \
../lib/scip/interfaces/jni/src/JniScipReaderFzn.c \
../lib/scip/interfaces/jni/src/JniScipReaderGms.c \
../lib/scip/interfaces/jni/src/JniScipReaderLp.c \
../lib/scip/interfaces/jni/src/JniScipReaderMps.c \
../lib/scip/interfaces/jni/src/JniScipReaderOpb.c \
../lib/scip/interfaces/jni/src/JniScipReaderOsil.c \
../lib/scip/interfaces/jni/src/JniScipReaderPip.c \
../lib/scip/interfaces/jni/src/JniScipReaderPpm.c \
../lib/scip/interfaces/jni/src/JniScipReaderRlp.c \
../lib/scip/interfaces/jni/src/JniScipReaderSol.c \
../lib/scip/interfaces/jni/src/JniScipReaderWbo.c \
../lib/scip/interfaces/jni/src/JniScipReaderZpl.c \
../lib/scip/interfaces/jni/src/JniScipRelax.c \
../lib/scip/interfaces/jni/src/JniScipSepa.c \
../lib/scip/interfaces/jni/src/JniScipSepaCgmip.c \
../lib/scip/interfaces/jni/src/JniScipSepaClique.c \
../lib/scip/interfaces/jni/src/JniScipSepaClosecuts.c \
../lib/scip/interfaces/jni/src/JniScipSepaCmir.c \
../lib/scip/interfaces/jni/src/JniScipSepaFlowcover.c \
../lib/scip/interfaces/jni/src/JniScipSepaGomory.c \
../lib/scip/interfaces/jni/src/JniScipSepaImpliedbounds.c \
../lib/scip/interfaces/jni/src/JniScipSepaIntobj.c \
../lib/scip/interfaces/jni/src/JniScipSepaMcf.c \
../lib/scip/interfaces/jni/src/JniScipSepaOddcycle.c \
../lib/scip/interfaces/jni/src/JniScipSepaRapidlearning.c \
../lib/scip/interfaces/jni/src/JniScipSepaStrongcg.c \
../lib/scip/interfaces/jni/src/JniScipSepaZerohalf.c \
../lib/scip/interfaces/jni/src/JniScipSol.c \
../lib/scip/interfaces/jni/src/JniScipTcliqueColoring.c \
../lib/scip/interfaces/jni/src/JniScipTcliqueDef.c \
../lib/scip/interfaces/jni/src/JniScipTree.c \
../lib/scip/interfaces/jni/src/JniScipVar.c 

OBJS += \
./lib/scip/interfaces/jni/src/JniScip.o \
./lib/scip/interfaces/jni/src/JniScipBranch.o \
./lib/scip/interfaces/jni/src/JniScipBranchAllfullstrong.o \
./lib/scip/interfaces/jni/src/JniScipBranchFullstrong.o \
./lib/scip/interfaces/jni/src/JniScipBranchInference.o \
./lib/scip/interfaces/jni/src/JniScipBranchLeastinf.o \
./lib/scip/interfaces/jni/src/JniScipBranchMostinf.o \
./lib/scip/interfaces/jni/src/JniScipBranchPscost.o \
./lib/scip/interfaces/jni/src/JniScipBranchRandom.o \
./lib/scip/interfaces/jni/src/JniScipBranchRelpscost.o \
./lib/scip/interfaces/jni/src/JniScipConflict.o \
./lib/scip/interfaces/jni/src/JniScipCons.o \
./lib/scip/interfaces/jni/src/JniScipConsAbspower.o \
./lib/scip/interfaces/jni/src/JniScipConsAnd.o \
./lib/scip/interfaces/jni/src/JniScipConsBivariate.o \
./lib/scip/interfaces/jni/src/JniScipConsBounddisjunction.o \
./lib/scip/interfaces/jni/src/JniScipConsConjunction.o \
./lib/scip/interfaces/jni/src/JniScipConsCountsols.o \
./lib/scip/interfaces/jni/src/JniScipConsCumulative.o \
./lib/scip/interfaces/jni/src/JniScipConsDisjunction.o \
./lib/scip/interfaces/jni/src/JniScipConsIndicator.o \
./lib/scip/interfaces/jni/src/JniScipConsIntegral.o \
./lib/scip/interfaces/jni/src/JniScipConsKnapsack.o \
./lib/scip/interfaces/jni/src/JniScipConsLinear.o \
./lib/scip/interfaces/jni/src/JniScipConsLinking.o \
./lib/scip/interfaces/jni/src/JniScipConsLogicor.o \
./lib/scip/interfaces/jni/src/JniScipConsNonlinear.o \
./lib/scip/interfaces/jni/src/JniScipConsOr.o \
./lib/scip/interfaces/jni/src/JniScipConsOrbitope.o \
./lib/scip/interfaces/jni/src/JniScipConsPseudoboolean.o \
./lib/scip/interfaces/jni/src/JniScipConsQuadratic.o \
./lib/scip/interfaces/jni/src/JniScipConsSetppc.o \
./lib/scip/interfaces/jni/src/JniScipConsSoc.o \
./lib/scip/interfaces/jni/src/JniScipConsSos1.o \
./lib/scip/interfaces/jni/src/JniScipConsSos2.o \
./lib/scip/interfaces/jni/src/JniScipConsSuperindicator.o \
./lib/scip/interfaces/jni/src/JniScipConsVarbound.o \
./lib/scip/interfaces/jni/src/JniScipConsXor.o \
./lib/scip/interfaces/jni/src/JniScipCutpool.o \
./lib/scip/interfaces/jni/src/JniScipDialog.o \
./lib/scip/interfaces/jni/src/JniScipDialogDefault.o \
./lib/scip/interfaces/jni/src/JniScipDisp.o \
./lib/scip/interfaces/jni/src/JniScipDispDefault.o \
./lib/scip/interfaces/jni/src/JniScipErrorCode.o \
./lib/scip/interfaces/jni/src/JniScipEvent.o \
./lib/scip/interfaces/jni/src/JniScipEventTypes.o \
./lib/scip/interfaces/jni/src/JniScipEventhdlr.o \
./lib/scip/interfaces/jni/src/JniScipException.o \
./lib/scip/interfaces/jni/src/JniScipExpr.o \
./lib/scip/interfaces/jni/src/JniScipHeur.o \
./lib/scip/interfaces/jni/src/JniScipHeurActconsdiving.o \
./lib/scip/interfaces/jni/src/JniScipHeurClique.o \
./lib/scip/interfaces/jni/src/JniScipHeurCoefdiving.o \
./lib/scip/interfaces/jni/src/JniScipHeurCrossover.o \
./lib/scip/interfaces/jni/src/JniScipHeurDins.o \
./lib/scip/interfaces/jni/src/JniScipHeurFeaspump.o \
./lib/scip/interfaces/jni/src/JniScipHeurFixandinfer.o \
./lib/scip/interfaces/jni/src/JniScipHeurFracdiving.o \
./lib/scip/interfaces/jni/src/JniScipHeurGuideddiving.o \
./lib/scip/interfaces/jni/src/JniScipHeurIntdiving.o \
./lib/scip/interfaces/jni/src/JniScipHeurIntshifting.o \
./lib/scip/interfaces/jni/src/JniScipHeurLinesearchdiving.o \
./lib/scip/interfaces/jni/src/JniScipHeurLocalbranching.o \
./lib/scip/interfaces/jni/src/JniScipHeurMutation.o \
./lib/scip/interfaces/jni/src/JniScipHeurNlpdiving.o \
./lib/scip/interfaces/jni/src/JniScipHeurObjpscostdiving.o \
./lib/scip/interfaces/jni/src/JniScipHeurOctane.o \
./lib/scip/interfaces/jni/src/JniScipHeurOneopt.o \
./lib/scip/interfaces/jni/src/JniScipHeurPscostdiving.o \
./lib/scip/interfaces/jni/src/JniScipHeurRens.o \
./lib/scip/interfaces/jni/src/JniScipHeurRins.o \
./lib/scip/interfaces/jni/src/JniScipHeurRootsoldiving.o \
./lib/scip/interfaces/jni/src/JniScipHeurRounding.o \
./lib/scip/interfaces/jni/src/JniScipHeurShiftandpropagate.o \
./lib/scip/interfaces/jni/src/JniScipHeurShifting.o \
./lib/scip/interfaces/jni/src/JniScipHeurSimplerounding.o \
./lib/scip/interfaces/jni/src/JniScipHeurSubnlp.o \
./lib/scip/interfaces/jni/src/JniScipHeurTrivial.o \
./lib/scip/interfaces/jni/src/JniScipHeurTrysol.o \
./lib/scip/interfaces/jni/src/JniScipHeurTwoopt.o \
./lib/scip/interfaces/jni/src/JniScipHeurUndercover.o \
./lib/scip/interfaces/jni/src/JniScipHeurVbounds.o \
./lib/scip/interfaces/jni/src/JniScipHeurVeclendiving.o \
./lib/scip/interfaces/jni/src/JniScipHeurZeroobj.o \
./lib/scip/interfaces/jni/src/JniScipHeurZirounding.o \
./lib/scip/interfaces/jni/src/JniScipImplics.o \
./lib/scip/interfaces/jni/src/JniScipLibraryLoader.o \
./lib/scip/interfaces/jni/src/JniScipLp.o \
./lib/scip/interfaces/jni/src/JniScipMessage.o \
./lib/scip/interfaces/jni/src/JniScipMessageDefault.o \
./lib/scip/interfaces/jni/src/JniScipMisc.o \
./lib/scip/interfaces/jni/src/JniScipNlp.o \
./lib/scip/interfaces/jni/src/JniScipNlpiIpopt.o \
./lib/scip/interfaces/jni/src/JniScipNodesel.o \
./lib/scip/interfaces/jni/src/JniScipNodeselBfs.o \
./lib/scip/interfaces/jni/src/JniScipNodeselDfs.o \
./lib/scip/interfaces/jni/src/JniScipNodeselEstimate.o \
./lib/scip/interfaces/jni/src/JniScipNodeselHybridestim.o \
./lib/scip/interfaces/jni/src/JniScipNodeselRestartdfs.o \
./lib/scip/interfaces/jni/src/JniScipParamset.o \
./lib/scip/interfaces/jni/src/JniScipPresol.o \
./lib/scip/interfaces/jni/src/JniScipPresolBoundshift.o \
./lib/scip/interfaces/jni/src/JniScipPresolComponents.o \
./lib/scip/interfaces/jni/src/JniScipPresolConvertinttobin.o \
./lib/scip/interfaces/jni/src/JniScipPresolDomcol.o \
./lib/scip/interfaces/jni/src/JniScipPresolDualfix.o \
./lib/scip/interfaces/jni/src/JniScipPresolGateextraction.o \
./lib/scip/interfaces/jni/src/JniScipPresolImplics.o \
./lib/scip/interfaces/jni/src/JniScipPresolInttobinary.o \
./lib/scip/interfaces/jni/src/JniScipPresolTrivial.o \
./lib/scip/interfaces/jni/src/JniScipPricer.o \
./lib/scip/interfaces/jni/src/JniScipProp.o \
./lib/scip/interfaces/jni/src/JniScipPropDualfix.o \
./lib/scip/interfaces/jni/src/JniScipPropGenvbounds.o \
./lib/scip/interfaces/jni/src/JniScipPropObbt.o \
./lib/scip/interfaces/jni/src/JniScipPropProbing.o \
./lib/scip/interfaces/jni/src/JniScipPropPseudoobj.o \
./lib/scip/interfaces/jni/src/JniScipPropRedcost.o \
./lib/scip/interfaces/jni/src/JniScipPropRootredcost.o \
./lib/scip/interfaces/jni/src/JniScipPropVbounds.o \
./lib/scip/interfaces/jni/src/JniScipReader.o \
./lib/scip/interfaces/jni/src/JniScipReaderBnd.o \
./lib/scip/interfaces/jni/src/JniScipReaderCcg.o \
./lib/scip/interfaces/jni/src/JniScipReaderCip.o \
./lib/scip/interfaces/jni/src/JniScipReaderCnf.o \
./lib/scip/interfaces/jni/src/JniScipReaderFix.o \
./lib/scip/interfaces/jni/src/JniScipReaderFzn.o \
./lib/scip/interfaces/jni/src/JniScipReaderGms.o \
./lib/scip/interfaces/jni/src/JniScipReaderLp.o \
./lib/scip/interfaces/jni/src/JniScipReaderMps.o \
./lib/scip/interfaces/jni/src/JniScipReaderOpb.o \
./lib/scip/interfaces/jni/src/JniScipReaderOsil.o \
./lib/scip/interfaces/jni/src/JniScipReaderPip.o \
./lib/scip/interfaces/jni/src/JniScipReaderPpm.o \
./lib/scip/interfaces/jni/src/JniScipReaderRlp.o \
./lib/scip/interfaces/jni/src/JniScipReaderSol.o \
./lib/scip/interfaces/jni/src/JniScipReaderWbo.o \
./lib/scip/interfaces/jni/src/JniScipReaderZpl.o \
./lib/scip/interfaces/jni/src/JniScipRelax.o \
./lib/scip/interfaces/jni/src/JniScipSepa.o \
./lib/scip/interfaces/jni/src/JniScipSepaCgmip.o \
./lib/scip/interfaces/jni/src/JniScipSepaClique.o \
./lib/scip/interfaces/jni/src/JniScipSepaClosecuts.o \
./lib/scip/interfaces/jni/src/JniScipSepaCmir.o \
./lib/scip/interfaces/jni/src/JniScipSepaFlowcover.o \
./lib/scip/interfaces/jni/src/JniScipSepaGomory.o \
./lib/scip/interfaces/jni/src/JniScipSepaImpliedbounds.o \
./lib/scip/interfaces/jni/src/JniScipSepaIntobj.o \
./lib/scip/interfaces/jni/src/JniScipSepaMcf.o \
./lib/scip/interfaces/jni/src/JniScipSepaOddcycle.o \
./lib/scip/interfaces/jni/src/JniScipSepaRapidlearning.o \
./lib/scip/interfaces/jni/src/JniScipSepaStrongcg.o \
./lib/scip/interfaces/jni/src/JniScipSepaZerohalf.o \
./lib/scip/interfaces/jni/src/JniScipSol.o \
./lib/scip/interfaces/jni/src/JniScipTcliqueColoring.o \
./lib/scip/interfaces/jni/src/JniScipTcliqueDef.o \
./lib/scip/interfaces/jni/src/JniScipTree.o \
./lib/scip/interfaces/jni/src/JniScipVar.o 

C_DEPS += \
./lib/scip/interfaces/jni/src/JniScip.d \
./lib/scip/interfaces/jni/src/JniScipBranch.d \
./lib/scip/interfaces/jni/src/JniScipBranchAllfullstrong.d \
./lib/scip/interfaces/jni/src/JniScipBranchFullstrong.d \
./lib/scip/interfaces/jni/src/JniScipBranchInference.d \
./lib/scip/interfaces/jni/src/JniScipBranchLeastinf.d \
./lib/scip/interfaces/jni/src/JniScipBranchMostinf.d \
./lib/scip/interfaces/jni/src/JniScipBranchPscost.d \
./lib/scip/interfaces/jni/src/JniScipBranchRandom.d \
./lib/scip/interfaces/jni/src/JniScipBranchRelpscost.d \
./lib/scip/interfaces/jni/src/JniScipConflict.d \
./lib/scip/interfaces/jni/src/JniScipCons.d \
./lib/scip/interfaces/jni/src/JniScipConsAbspower.d \
./lib/scip/interfaces/jni/src/JniScipConsAnd.d \
./lib/scip/interfaces/jni/src/JniScipConsBivariate.d \
./lib/scip/interfaces/jni/src/JniScipConsBounddisjunction.d \
./lib/scip/interfaces/jni/src/JniScipConsConjunction.d \
./lib/scip/interfaces/jni/src/JniScipConsCountsols.d \
./lib/scip/interfaces/jni/src/JniScipConsCumulative.d \
./lib/scip/interfaces/jni/src/JniScipConsDisjunction.d \
./lib/scip/interfaces/jni/src/JniScipConsIndicator.d \
./lib/scip/interfaces/jni/src/JniScipConsIntegral.d \
./lib/scip/interfaces/jni/src/JniScipConsKnapsack.d \
./lib/scip/interfaces/jni/src/JniScipConsLinear.d \
./lib/scip/interfaces/jni/src/JniScipConsLinking.d \
./lib/scip/interfaces/jni/src/JniScipConsLogicor.d \
./lib/scip/interfaces/jni/src/JniScipConsNonlinear.d \
./lib/scip/interfaces/jni/src/JniScipConsOr.d \
./lib/scip/interfaces/jni/src/JniScipConsOrbitope.d \
./lib/scip/interfaces/jni/src/JniScipConsPseudoboolean.d \
./lib/scip/interfaces/jni/src/JniScipConsQuadratic.d \
./lib/scip/interfaces/jni/src/JniScipConsSetppc.d \
./lib/scip/interfaces/jni/src/JniScipConsSoc.d \
./lib/scip/interfaces/jni/src/JniScipConsSos1.d \
./lib/scip/interfaces/jni/src/JniScipConsSos2.d \
./lib/scip/interfaces/jni/src/JniScipConsSuperindicator.d \
./lib/scip/interfaces/jni/src/JniScipConsVarbound.d \
./lib/scip/interfaces/jni/src/JniScipConsXor.d \
./lib/scip/interfaces/jni/src/JniScipCutpool.d \
./lib/scip/interfaces/jni/src/JniScipDialog.d \
./lib/scip/interfaces/jni/src/JniScipDialogDefault.d \
./lib/scip/interfaces/jni/src/JniScipDisp.d \
./lib/scip/interfaces/jni/src/JniScipDispDefault.d \
./lib/scip/interfaces/jni/src/JniScipErrorCode.d \
./lib/scip/interfaces/jni/src/JniScipEvent.d \
./lib/scip/interfaces/jni/src/JniScipEventTypes.d \
./lib/scip/interfaces/jni/src/JniScipEventhdlr.d \
./lib/scip/interfaces/jni/src/JniScipException.d \
./lib/scip/interfaces/jni/src/JniScipExpr.d \
./lib/scip/interfaces/jni/src/JniScipHeur.d \
./lib/scip/interfaces/jni/src/JniScipHeurActconsdiving.d \
./lib/scip/interfaces/jni/src/JniScipHeurClique.d \
./lib/scip/interfaces/jni/src/JniScipHeurCoefdiving.d \
./lib/scip/interfaces/jni/src/JniScipHeurCrossover.d \
./lib/scip/interfaces/jni/src/JniScipHeurDins.d \
./lib/scip/interfaces/jni/src/JniScipHeurFeaspump.d \
./lib/scip/interfaces/jni/src/JniScipHeurFixandinfer.d \
./lib/scip/interfaces/jni/src/JniScipHeurFracdiving.d \
./lib/scip/interfaces/jni/src/JniScipHeurGuideddiving.d \
./lib/scip/interfaces/jni/src/JniScipHeurIntdiving.d \
./lib/scip/interfaces/jni/src/JniScipHeurIntshifting.d \
./lib/scip/interfaces/jni/src/JniScipHeurLinesearchdiving.d \
./lib/scip/interfaces/jni/src/JniScipHeurLocalbranching.d \
./lib/scip/interfaces/jni/src/JniScipHeurMutation.d \
./lib/scip/interfaces/jni/src/JniScipHeurNlpdiving.d \
./lib/scip/interfaces/jni/src/JniScipHeurObjpscostdiving.d \
./lib/scip/interfaces/jni/src/JniScipHeurOctane.d \
./lib/scip/interfaces/jni/src/JniScipHeurOneopt.d \
./lib/scip/interfaces/jni/src/JniScipHeurPscostdiving.d \
./lib/scip/interfaces/jni/src/JniScipHeurRens.d \
./lib/scip/interfaces/jni/src/JniScipHeurRins.d \
./lib/scip/interfaces/jni/src/JniScipHeurRootsoldiving.d \
./lib/scip/interfaces/jni/src/JniScipHeurRounding.d \
./lib/scip/interfaces/jni/src/JniScipHeurShiftandpropagate.d \
./lib/scip/interfaces/jni/src/JniScipHeurShifting.d \
./lib/scip/interfaces/jni/src/JniScipHeurSimplerounding.d \
./lib/scip/interfaces/jni/src/JniScipHeurSubnlp.d \
./lib/scip/interfaces/jni/src/JniScipHeurTrivial.d \
./lib/scip/interfaces/jni/src/JniScipHeurTrysol.d \
./lib/scip/interfaces/jni/src/JniScipHeurTwoopt.d \
./lib/scip/interfaces/jni/src/JniScipHeurUndercover.d \
./lib/scip/interfaces/jni/src/JniScipHeurVbounds.d \
./lib/scip/interfaces/jni/src/JniScipHeurVeclendiving.d \
./lib/scip/interfaces/jni/src/JniScipHeurZeroobj.d \
./lib/scip/interfaces/jni/src/JniScipHeurZirounding.d \
./lib/scip/interfaces/jni/src/JniScipImplics.d \
./lib/scip/interfaces/jni/src/JniScipLibraryLoader.d \
./lib/scip/interfaces/jni/src/JniScipLp.d \
./lib/scip/interfaces/jni/src/JniScipMessage.d \
./lib/scip/interfaces/jni/src/JniScipMessageDefault.d \
./lib/scip/interfaces/jni/src/JniScipMisc.d \
./lib/scip/interfaces/jni/src/JniScipNlp.d \
./lib/scip/interfaces/jni/src/JniScipNlpiIpopt.d \
./lib/scip/interfaces/jni/src/JniScipNodesel.d \
./lib/scip/interfaces/jni/src/JniScipNodeselBfs.d \
./lib/scip/interfaces/jni/src/JniScipNodeselDfs.d \
./lib/scip/interfaces/jni/src/JniScipNodeselEstimate.d \
./lib/scip/interfaces/jni/src/JniScipNodeselHybridestim.d \
./lib/scip/interfaces/jni/src/JniScipNodeselRestartdfs.d \
./lib/scip/interfaces/jni/src/JniScipParamset.d \
./lib/scip/interfaces/jni/src/JniScipPresol.d \
./lib/scip/interfaces/jni/src/JniScipPresolBoundshift.d \
./lib/scip/interfaces/jni/src/JniScipPresolComponents.d \
./lib/scip/interfaces/jni/src/JniScipPresolConvertinttobin.d \
./lib/scip/interfaces/jni/src/JniScipPresolDomcol.d \
./lib/scip/interfaces/jni/src/JniScipPresolDualfix.d \
./lib/scip/interfaces/jni/src/JniScipPresolGateextraction.d \
./lib/scip/interfaces/jni/src/JniScipPresolImplics.d \
./lib/scip/interfaces/jni/src/JniScipPresolInttobinary.d \
./lib/scip/interfaces/jni/src/JniScipPresolTrivial.d \
./lib/scip/interfaces/jni/src/JniScipPricer.d \
./lib/scip/interfaces/jni/src/JniScipProp.d \
./lib/scip/interfaces/jni/src/JniScipPropDualfix.d \
./lib/scip/interfaces/jni/src/JniScipPropGenvbounds.d \
./lib/scip/interfaces/jni/src/JniScipPropObbt.d \
./lib/scip/interfaces/jni/src/JniScipPropProbing.d \
./lib/scip/interfaces/jni/src/JniScipPropPseudoobj.d \
./lib/scip/interfaces/jni/src/JniScipPropRedcost.d \
./lib/scip/interfaces/jni/src/JniScipPropRootredcost.d \
./lib/scip/interfaces/jni/src/JniScipPropVbounds.d \
./lib/scip/interfaces/jni/src/JniScipReader.d \
./lib/scip/interfaces/jni/src/JniScipReaderBnd.d \
./lib/scip/interfaces/jni/src/JniScipReaderCcg.d \
./lib/scip/interfaces/jni/src/JniScipReaderCip.d \
./lib/scip/interfaces/jni/src/JniScipReaderCnf.d \
./lib/scip/interfaces/jni/src/JniScipReaderFix.d \
./lib/scip/interfaces/jni/src/JniScipReaderFzn.d \
./lib/scip/interfaces/jni/src/JniScipReaderGms.d \
./lib/scip/interfaces/jni/src/JniScipReaderLp.d \
./lib/scip/interfaces/jni/src/JniScipReaderMps.d \
./lib/scip/interfaces/jni/src/JniScipReaderOpb.d \
./lib/scip/interfaces/jni/src/JniScipReaderOsil.d \
./lib/scip/interfaces/jni/src/JniScipReaderPip.d \
./lib/scip/interfaces/jni/src/JniScipReaderPpm.d \
./lib/scip/interfaces/jni/src/JniScipReaderRlp.d \
./lib/scip/interfaces/jni/src/JniScipReaderSol.d \
./lib/scip/interfaces/jni/src/JniScipReaderWbo.d \
./lib/scip/interfaces/jni/src/JniScipReaderZpl.d \
./lib/scip/interfaces/jni/src/JniScipRelax.d \
./lib/scip/interfaces/jni/src/JniScipSepa.d \
./lib/scip/interfaces/jni/src/JniScipSepaCgmip.d \
./lib/scip/interfaces/jni/src/JniScipSepaClique.d \
./lib/scip/interfaces/jni/src/JniScipSepaClosecuts.d \
./lib/scip/interfaces/jni/src/JniScipSepaCmir.d \
./lib/scip/interfaces/jni/src/JniScipSepaFlowcover.d \
./lib/scip/interfaces/jni/src/JniScipSepaGomory.d \
./lib/scip/interfaces/jni/src/JniScipSepaImpliedbounds.d \
./lib/scip/interfaces/jni/src/JniScipSepaIntobj.d \
./lib/scip/interfaces/jni/src/JniScipSepaMcf.d \
./lib/scip/interfaces/jni/src/JniScipSepaOddcycle.d \
./lib/scip/interfaces/jni/src/JniScipSepaRapidlearning.d \
./lib/scip/interfaces/jni/src/JniScipSepaStrongcg.d \
./lib/scip/interfaces/jni/src/JniScipSepaZerohalf.d \
./lib/scip/interfaces/jni/src/JniScipSol.d \
./lib/scip/interfaces/jni/src/JniScipTcliqueColoring.d \
./lib/scip/interfaces/jni/src/JniScipTcliqueDef.d \
./lib/scip/interfaces/jni/src/JniScipTree.d \
./lib/scip/interfaces/jni/src/JniScipVar.d 


# Each subdirectory must supply rules for building sources it contributes
lib/scip/interfaces/jni/src/%.o: ../lib/scip/interfaces/jni/src/%.c
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C Compiler'
	gcc -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


