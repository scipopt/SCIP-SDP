SCIP-SDP - A framework for solving mixed-integer semidefinite programs
======================================================================

SCIP-SDP is a plugin for SCIP to solve mixed integer semidefinite
programs (MISDPs), i.e., semidefinite programs (SDPs) in which some
variables are required to be integral.

It combines the branch-and-bound framework of SCIP with interior-point
SDP-solvers to solve MISDPs using either a nonlinear branch-and-bound
approach or an outer-approximation-based cutting-plane approach using
linear programs (LPs). In addition to providing a constraint handler
for SDP-constraints and a relaxator to solve continuous
SDP-relaxations using interior-point solvers, SCIP-SDP adds several
heuristics and propagators to SCIP.

The MISDPs can be read in using either the CBF-format or an extended
SDPA-format with support for integrality as well as rank-1
constraints. For a description of the extended SDPA-format see the
file "sdpa_format.txt". The CBF-format is supported up to version 2,
see https://cblib.zib.de, and has also been extended to support rank-1
constraints.

To use the nonlinear branch-and-bound approach, one of the following
SDP-solvers needs to be installed: DSDP, SDPA, or MOSEK. For more
information about the installation of SCIP-SDP see the INSTALL file.

Features
--------

- SCIP-SDP can read MISDPs in CBF or extended SDPA format.

- SCIP-SDP contains presolving and propagation methods as well as primal
  heuristics.

- SCIP-SDP supports rank-1 constraints by adding quadratic constraints
  for each 2 by 2 minor. Such problems are usually (very) hard to
  solve.

- One can switch between solving SDPs or using an outer approximation
  via LP-relaxations.

- SCIP-SDP extends SCIP and can be incorporated into other codes
  similar to how this can be done for SCIP.


Interesting Parameters
----------------------

SCIP-SDP features many parameters to determine its behavior and it
inherits all SCIP parameters. Here, we highlight some:

- "misc/solvesdps": Determines whether SDPs (1) or LPs (0) are
  solved. It depends on the instance which approach is faster.

- "relaxing/SDP/sdpsolverthreads" sets the number of threads used for
  solving SDPs. By default it is set to 1, since this seems to be the
  fastest for most instances. For larger SDPs it might help to
  increase this number, where a value of "-1" corresponds to an
  automatic choice, if this is supported by the SDP solver.

- Tolerances: "numerics/feastol" and "numerics/dualfeastol" (default
  1e-5) are SCIP parameters that determine the feasibility
  tolerances. Note that the default is a bit weaker than the default
  SCIP parameter values (1e-6 and 1e-7, respectively).
  "relaxing/SDP/sdpsolverfeastol" and "relaxing/SDP/sdpsolvergaptol"
  (default 1e-5) determine the tolerances used for the SDP-solvers.
  For Mosek the feasibility tolerance is always tightened by 0.1,
  because this generated more reliable results. Depending on the
  instance, changing these parameters might have a dramatic effect on
  performance and correctness of the results (due to numerical
  issues).

- Slater condition: The solution process of interior-point methods for
  SDPs depends on the Slater condition. One can determine the Slater
  condition for the primal and dual problem by changing
  "relaxing/SDP/slatercheck" (0: no, 1: yes but only for statistics, 2:
  yes and print warning for every problem not satisfying primal and
  dual Slater condition).

Documentation
-------------

- A doxygen documentation is given at
  https://wwwopt.mathematik.tu-darmstadt.de/scipsdp/

- Reports for the SCIP Optimization Suites listed at
  https://www.scipopt.org/index.php#cite each contain a part on the
  developments in SCIP-SDP.

- Many methods used in SCIP-SDP are described in the following dissertations:

  Sonja Mars (2013) "Mixed-Integer Semidefinite Programming with an
  Application to Truss Topology Design"

  Tristan Gally (2019) "Computational Mixed-Integer Semidefinite
  Programming"

  Frederic Matter (2022) "Sparse Recovery Under Side Constraints Using
  Null Space Properties"
