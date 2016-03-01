<p>
   SCIP-SDP is a plugin for SCIP to solve mixed integer semidefinite programs (MISDPs). It combines the branch-and-bound framework of SCIP with interior-point SDP-solvers. It provides the data handling, some presolving and propagation as well as a reader for a <a href="http://www.opt.tu-darmstadt.de/scipsdp/downloads/data_format.txt">modified sparse SDPA-format with additional lines for integrality constraints</a>. It is possible to solve the resulting SDP-relaxations using a linear approximation procedure, but for full functionality one of the following SDP-solvers needs to be installed:
</p>
<ul>
   <li><a href="http://www.mcs.anl.gov/hs/software/DSDP/">DSDP</a></li>
   <li><a href="http://sdpa.sourceforge.net/">SDPA</a></li>
</ul>
<p>
<p>
   Please note that the interface to SDPA is still in beta state. It works well for some instances and is faster than DSDP for those, but currently fails for others because of numerical problems, as some parameters need further tuning.
</p>
</p>
<p>The solution process of interior-point methods for SDPs is highly dependent on the Slater condition. One of the main purposes of the code is ensuring that the slater condition is not harmed by fixing variables in the branch-and-bound process. In some cases however the combination of variable fixings and specific linear or semidefinite constraints might still lead to relaxations for which the Slater condition no longer holds. In this case the SDP-solvers may be unable to solve the relaxations or even return wrong results, which cannot be compensated. For this purpose there is the possibility to check the Slater condition for the dual problem (which still does not guarantee it for the primal) before the solution of each SDP by setting a SCIP parameter, for details see the documentation.</p>
