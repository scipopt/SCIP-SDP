<p>
   SCIP-SDP is a plugin for SCIP to solve mixed integer semidefinite programs (MISDPs) of the form
</p>
<p style= "text-indent:1.5em;">
   <img src="http://www.opt.tu-darmstadt.de/scipsdp/downloads/misdp.png">
</p>
<p>
   It combines the branch-and-bound framework of SCIP with interior-point SDP-solvers to solve MISDPs using either a
   nonlinear branch-and-bound approach or an outer-approximation-based cutting-plane approach. In addition to providing
   a constraint handler for SDP-constraints and a relaxator to solve continuous SDP-relaxations using interior-point
   solvers, SCIPSDP adds several heuristics and propagators to SCIP. The MISDPs can be read in using either an 
   <a href="http://www.opt.tu-darmstadt.de/scipsdp/downloads/data_format.txt">extended SDPA-format</a> or the 
   <a href="http://cblib.zib.de/doc/format2.pdf">CBF-format</a>. To use the nonlinear 
   branch-and-bound approach one of the following SDP-solvers needs to be installed:
</p>
<ul>
   <li><a href="http://www.mcs.anl.gov/hs/software/DSDP/">DSDP</a></li>
   <li><a href="http://sdpa.sourceforge.net/">SDPA</a></li>
   <li><a href="https://www.mosek.com/">MOSEK</a></li>
</ul>
<p> The solution process of interior-point methods for SDPs is highly dependent on the Slater condition. One of the main
    purposes of the code is handling cases where the Slater condition does not hold using a penalty approach. However, in
    some cases the SDP-solvers may still fail because of numerical difficulties or even return wrong results, which cannot
    be compensated. For this purpose there is the possibility to check the Slater condition for the primal and dual problem
    before the solution of each SDP by setting a SCIP parameter, for details see the documentation.</p>
