<p>
   SCIP-SDP is a plugin for <a href="http://scip.zib.de">SCIP</a> to solve mixed integer semidefinite programs (MISDPs) of the form
   \[
   \begin{aligned}
	 \inf \quad \thinspace & b^\top y \\
	 \text{s.t.} \quad & \sum_{i=1}^m A_i\, y_i - A_0 \succeq 0, \\
	 & y_i \in \mathbb{Z} && \forall\, i \in \mathcal{I}.
   \end{aligned}
   \]
   SCIP-SDP allows to solve MISDPs using a nonlinear branch-and-bound
   approach or a linear programming cutting-plane approach. In the
   first case (the default), the semidefinite programming (SDP)
   relaxations are solve using interior-point SDP-solvers. In the
   second case, cutting planes based on eigenvector are generated.
   SCIP-SDP is based on <a href="scipopt.org">SCIP</a>.
</p>

<p>
   In addition to providing a constraint handler for SDP-constraints
   and a relaxator to solve continuous SDP-relaxations using
   interior-point solvers, SCIPSDP adds several heuristics and
   propagators to SCIP. The MISDPs can be read in using either an <a
   href="http://www.opt.tu-darmstadt.de/scipsdp/downloads/data_format.txt">extended
   SDPA-format</a> or the <a
   href="http://cblib.zib.de/doc/format2.pdf">CBF-format</a>. To use
   the nonlinear branch-and-bound approach one of the following
   SDP-solvers needs to be installed:
</p>
<ul>
   <li><a href="http://www.mcs.anl.gov/hs/software/DSDP/">DSDP</a></li>
   <li><a href="http://sdpa.sourceforge.net/">SDPA</a></li>
   <li><a href="https://www.mosek.com/">MOSEK</a></li>
</ul>

<p>
   Mixed-integer semidefinite programs are sometimes numerically
   challenging to solve. One reason is that the Slater Condition for
   the SDP-relaxations may not hold in some of the nodes. SCIP-SDP
   implements several methods to recover from a failure to accurately
   solve the relaxation.
</p>
