/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the                                 */
/*      SDP-Package for SCIP: a solving framework for                        */
/*                            mixed-integer semidefinite programms           */
/*                                                                           */
/* Copyright (C) 2011-2012 Discrete Optimization, TU Darmstadt               */
/*                         EDOM, FAU Erlangen-Nürnberg                       */
/*                                                                           */
/* This program is free software; you can redistribute it and/or             */
/* modify it under the terms of the GNU Lesser General Public License        */
/* as published by the Free Software Foundation; either version 3            */
/* of the License, or (at your option) any later version.                    */
/*                                                                           */
/* This program is distributed in the hope that it will be useful,           */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of            */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             */
/* GNU Lesser General Public License for more details.                       */
/*                                                                           */
/* You should have received a copy of the GNU Lesser General Public License  */
/* along with this program; if not, write to the Free Software               */
/* Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301, USA.*/
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   objconshdlr_sdp.h
 * @brief  constrainthandler for SDPs
 * @author Sonja Mars, Lars Schewe
 */

#ifndef __SCIP_OBJCONSHDLR_SDP_H__
#define __SCIP_OBJCONSHDLR_SDP_H__


#include "objscip/objconshdlr.h"        // for ObjConshdlr
#include "scip/def.h"                   // for SCIP_Bool, FALSE, TRUE
#include "scip/type_cons.h"             // for SCIP_CONS, SCIP_CONSHDLR, etc
#include "scip/type_result.h"           // for SCIP_RESULT
#include "scip/type_retcode.h"          // for SCIP_RETCODE
#include "scip/type_scip.h"             // for SCIP
#include "scip/type_sol.h"              // for SCIP_SOL
class SdpCone;  // lines 20-20

/** C++ wrapper object for sdp-constraint handlers */
class ObjConshdlrSdp : public scip::ObjConshdlr
{

 public:
   /** default constructor */
 ObjConshdlrSdp(SCIP* scip)
    : scip::ObjConshdlr(scip, "SDP", "SDPconic-constraint",
       1000000, -2000000, -2000000, 1, -1, 1, -1,
       FALSE, FALSE, FALSE, TRUE, 1)
   {}

   /** destructor */
   virtual ~ObjConshdlrSdp() {}

   virtual SCIP_RETCODE scip_init(
      SCIP*              scip,               /**< SCIP data structure */
      SCIP_CONSHDLR*     conshdlr,           /**< the constraint handler itself */
      SCIP_CONS**        conss,              /**< array of constraints in transformed problem */
      int                nconss              /**< number of constraints in transformed problem */
      );

   /** variable rounding lock method of constraint handler
    *
    *  This method is called, after a constraint is added or removed from the transformed problem.
    *  It should update the rounding locks of all associated variables with calls to SCIPaddVarLocks(),
    *  depending on the way, the variable is involved in the constraint:
    *  - If the constraint may get violated by decreasing the value of a variable, it should call
    *    SCIPaddVarLocks(scip, var, nlockspos, nlocksneg), saying that rounding down is potentially rendering the
    *    (positive) constraint infeasible and rounding up is potentially rendering the negation of the constraint
    *    infeasible.
    *  - If the constraint may get violated by increasing the value of a variable, it should call
    *    SCIPaddVarLocks(scip, var, nlocksneg, nlockspos), saying that rounding up is potentially rendering the
    *    constraint's negation infeasible and rounding up is potentially rendering the constraint itself
    *    infeasible.
    *  - If the constraint may get violated by changing the variable in any direction, it should call
    *    SCIPaddVarLocks(scip, var, nlockspos + nlocksneg, nlockspos + nlocksneg).
    *
    *  Consider the linear constraint "3x -5y +2z <= 7" as an example. The variable rounding lock method of the
    *  linear constraint handler should call SCIPaddVarLocks(scip, x, nlocksneg, nlockspos),
    *  SCIPaddVarLocks(scip, y, nlockspos, nlocksneg) and SCIPaddVarLocks(scip, z, nlocksneg, nlockspos) to tell SCIP,
    *  that rounding up of x and z and rounding down of y can destroy the feasibility of the constraint, while rounding
    *  down of x and z and rounding up of y can destroy the feasibility of the constraint's negation "3x -5y +2z > 7".
    *  A linear constraint "2 <= 3x -5y +2z <= 7" should call
    *  SCIPaddVarLocks(scip, ..., nlockspos + nlocksneg, nlockspos + nlocksneg) on all variables, since rounding in both
    *  directions of each variable can destroy both the feasibility of the constraint and it's negation
    *  "3x -5y +2z < 2  or  3x -5y +2z > 7".
    *
    *  If the constraint itself contains other constraints as sub constraints (e.g. the "or" constraint concatenation
    *  "c(x) or d(x)"), the rounding lock methods of these constraints should be called in a proper way.
    *  - If the constraint may get violated by the violation of the sub constraint c, it should call
    *    SCIPaddConsLocks(scip, c, nlockspos, nlocksneg), saying that infeasibility of c may lead to infeasibility of
    *    the (positive) constraint, and infeasibility of c's negation (i.e. feasibility of c) may lead to infeasibility
    *    of the constraint's negation (i.e. feasibility of the constraint).
    *  - If the constraint may get violated by the feasibility of the sub constraint c, it should call
    *    SCIPaddConsLocks(scip, c, nlocksneg, nlockspos), saying that infeasibility of c may lead to infeasibility of
    *    the constraint's negation (i.e. feasibility of the constraint), and infeasibility of c's negation (i.e. feasibility
    *    of c) may lead to infeasibility of the (positive) constraint.
    *  - If the constraint may get violated by any change in the feasibility of the sub constraint c, it should call
    *    SCIPaddConsLocks(scip, c, nlockspos + nlocksneg, nlockspos + nlocksneg).
    *
    *  Consider the or concatenation "c(x) or d(x)". The variable rounding lock method of the or constraint handler
    *  should call SCIPaddConsLocks(scip, c, nlockspos, nlocksneg) and SCIPaddConsLocks(scip, d, nlockspos, nlocksneg)
    *  to tell SCIP, that infeasibility of c and d can lead to infeasibility of "c(x) or d(x)".
    *
    *  As a second example, consider the equivalence constraint "y <-> c(x)" with variable y and constraint c. The
    *  constraint demands, that y == 1 if and only if c(x) is satisfied. The variable lock method of the corresponding
    *  constraint handler should call SCIPaddVarLocks(scip, y, nlockspos + nlocksneg, nlockspos + nlocksneg) and
    *  SCIPaddConsLocks(scip, c, nlockspos + nlocksneg, nlockspos + nlocksneg), because any modification to the
    *  value of y or to the feasibility of c can alter the feasibility of the equivalence constraint.
    */
   virtual SCIP_RETCODE scip_lock(
      SCIP*              scip,               /**< SCIP data structure */
      SCIP_CONSHDLR*     conshdlr,           /**< the constraint handler itself */
      SCIP_CONS*         cons,               /**< the constraint that should lock rounding of its variables, or NULL if the
                                              *   constraint handler does not need constraints */
      int                nlockspos,          /**< no. of times, the roundings should be locked for the constraint */
      int                nlocksneg           /**< no. of times, the roundings should be locked for the constraint's negation */
      ) ;


   /** presolving initialization method of constraint handler (called when presolving is about to begin)
    *
    *  This method is called when the presolving process is about to begin, even if presolving is turned off.
    *  The constraint handler may use this call to initialize its presolving data, or to modify its constraints
    *  before the presolving process begins.
    *  Necessary constraint modifications that have to be performed even if presolving is turned off should be done here
    *  or in the presolving deinitialization call.
    *
    *  possible return values for *result:
    *  - SCIP_UNBOUNDED  : at least one variable is not bounded by any constraint in obj. direction -> problem is unbounded
    *  - SCIP_CUTOFF     : at least one constraint is infeasible in the variable's bounds -> problem is infeasible
    *  - SCIP_FEASIBLE   : no infeasibility nor unboundness could be found
    */
   virtual SCIP_RETCODE scip_initpre(
      SCIP*              scip,               /**< SCIP data structure */
      SCIP_CONSHDLR*     conshdlr,           /**< the constraint handler itself */
      SCIP_CONS**        conss,              /**< array of constraints in transformed problem */
      int                nconss             /**< number of constraints in transformed problem */
      );


   /** presolving method of constraint handler
    *
    *  The presolver should go through the variables and constraints and tighten the domains or
    *  constraints. Each tightening should increase the given total number of changes.
    *
    *  @note the counters state the changes since the last call including the changes of this presolving method during
    *        its call
    *
    *  possible return values for *result:
    *  - SCIP_UNBOUNDED  : at least one variable is not bounded by any constraint in obj. direction -> problem is unbounded
    *  - SCIP_CUTOFF     : at least one constraint is infeasible in the variable's bounds -> problem is infeasible
    *  - SCIP_SUCCESS    : the presolving method found a reduction
    *  - SCIP_DIDNOTFIND : the presolving method searched, but did not find a presolving change
    *  - SCIP_DIDNOTRUN  : the presolving method was skipped
    *  - SCIP_DELAYED    : the presolving method was skipped, but should be called again
    */
   virtual SCIP_RETCODE scip_presol(
      SCIP*              scip,               /**< SCIP data structure */
      SCIP_CONSHDLR*     conshdlr,           /**< the constraint handler itself */
      SCIP_CONS**        conss,              /**< array of constraints to process */
      int                nconss,             /**< no. of constraints to process */
      int                nrounds,            /**< no. of presolving rounds already done */
      int                nnewfixedvars,      /**< no. of variables fixed since last call to presolving method */
      int                nnewaggrvars,       /**< no. of variables aggregated since last call to presolving method */
      int                nnewchgvartypes,    /**< no. of variable type changes since last call to presolving method */
      int                nnewchgbds,         /**< no. of variable bounds tightend since last call to presolving method */
      int                nnewholes,          /**< no. of domain holes added since last call to presolving method */
      int                nnewdelconss,       /**< no. of deleted constraints since last call to presolving method */
      int                nnewaddconss,       /**< no. of added constraints since last call to presolving method */
      int                nnewupgdconss,      /**< no. of upgraded constraints since last call to presolving method */
      int                nnewchgcoefs,       /**< no. of changed coefficients since last call to presolving method */
      int                nnewchgsides,       /**< no. of changed left or right hand sides since last call to presolving method */
      int*               nfixedvars,         /**< pointer to count total number of variables fixed of all presolvers */
      int*               naggrvars,          /**< pointer to count total number of variables aggregated of all presolvers */
      int*               nchgvartypes,       /**< pointer to count total number of variable type changes of all presolvers */
      int*               nchgbds,            /**< pointer to count total number of variable bounds tightend of all presolvers */
      int*               naddholes,          /**< pointer to count total number of domain holes added of all presolvers */
      int*               ndelconss,          /**< pointer to count total number of deleted constraints of all presolvers */
      int*               naddconss,          /**< pointer to count total number of added constraints of all presolvers */
      int*               nupgdconss,         /**< pointer to count total number of upgraded constraints of all presolvers */
      int*               nchgcoefs,          /**< pointer to count total number of changed coefficients of all presolvers */
      int*               nchgsides,          /**< pointer to count total number of changed sides of all presolvers */
      SCIP_RESULT*       result              /**< pointer to store the result of the presolving call */
      );

   /** transforms constraint data into data belonging to the transformed problem */
   virtual SCIP_RETCODE scip_trans(
      SCIP*              scip,               /**< SCIP data structure */
      SCIP_CONSHDLR*     conshdlr,           /**< the constraint handler itself */
      SCIP_CONS*         sourcecons,         /**< source constraint to transform */
      SCIP_CONS**        targetcons          /**< pointer to store created target constraint */
      );


   /** feasibility check method of constraint handler for primal solutions
    *
    *  The given solution has to be checked for feasibility.
    *
    *  The check methods of the active constraint handlers are called in decreasing order of their check
    *  priorities until the first constraint handler returned with the result SCIP_INFEASIBLE.
    *  The integrality constraint handler has a check priority of zero. A constraint handler which can
    *  (or wants) to check its constraints only for integral solutions should have a negative check priority
    *  (e.g. the alldiff-constraint can only operate on integral solutions).
    *  A constraint handler which wants to check feasibility even on non-integral solutions must have a
    *  check priority greater than zero (e.g. if the check is much faster than testing all variables for
    *  integrality).
    *
    *  In some cases, integrality conditions or rows of the current LP do not have to be checked, because their
    *  feasibility is already checked or implicitly given. In these cases, 'checkintegrality' or
    *  'checklprows' is FALSE.
    *
    *  possible return values for *result:
    *  - SCIP_INFEASIBLE : at least one constraint of the handler is infeasible
    *  - SCIP_FEASIBLE   : all constraints of the handler are feasible
    */
   virtual SCIP_RETCODE scip_check(
      SCIP*              scip,               /**< SCIP data structure */
      SCIP_CONSHDLR*     conshdlr,           /**< the constraint handler itself */
      SCIP_CONS**        conss,              /**< array of constraints to process */
      int                nconss,             /**< number of constraints to process */
      SCIP_SOL*          sol,                /**< the solution to check feasibility for */
      SCIP_Bool          checkintegrality,   /**< has integrality to be checked? */
      SCIP_Bool          checklprows,        /**< have current LP rows to be checked? */
      SCIP_Bool          printreason,        /**< should the reason for the violation be printed? */
      SCIP_RESULT*       result              /**< pointer to store the result of the feasibility checking call */
      );

   /** constraint enforcing method of constraint handler for pseudo solutions
    *
    *  The method is called at the end of the node processing loop for a node where the LP was not solved.
    *  The pseudo solution has to be checked for feasibility. If possible, an infeasibility should be resolved by
    *  branching, reducing a variable's domain to exclude the solution or adding an additional constraint.
    *  Separation is not possible, since the LP is not processed at the current node. All LP informations like
    *  LP solution, slack values, or reduced costs are invalid and must not be accessed.
    *
    *  Like in the enforcing method for LP solutions, the enforcing methods of the active constraint handlers are
    *  called in decreasing order of their enforcing priorities until the first constraint handler returned with
    *  the value SCIP_CUTOFF, SCIP_REDUCEDDOM, SCIP_CONSADDED, SCIP_BRANCHED, or SCIP_SOLVELP.
    *
    *  The first nusefulconss constraints are the ones, that are identified to likely be violated. The enforcing
    *  method should process the useful constraints first. The other nconss - nusefulconss constraints should only
    *  be enforced, if no violation was found in the useful constraints.
    *
    *  If the pseudo solution's objective value is lower than the lower bound of the node, it cannot be feasible
    *  and the enforcing method may skip it's check and set *result to SCIP_DIDNOTRUN. However, it can also process
    *  its constraints and return any other possible result code.
    *
    *  possible return values for *result (if more than one applies, the first in the list should be used):
    *  - SCIP_CUTOFF     : the node is infeasible in the variable's bounds and can be cut off
    *  - SCIP_CONSADDED  : an additional constraint was generated
    *  - SCIP_REDUCEDDOM : a variable's domain was reduced
    *  - SCIP_BRANCHED   : no changes were made to the problem, but a branching was applied to resolve an infeasibility
    *  - SCIP_SOLVELP    : at least one constraint is infeasible, and this can only be resolved by solving the SCIP_LP
    *  - SCIP_INFEASIBLE : at least one constraint is infeasible, but it was not resolved
    *  - SCIP_FEASIBLE   : all constraints of the handler are feasible
    *  - SCIP_DIDNOTRUN  : the enforcement was skipped (only possible, if objinfeasible is TRUE)
    */
   virtual SCIP_RETCODE scip_enfops(
      SCIP*              scip,               /**< SCIP data structure */
      SCIP_CONSHDLR*     conshdlr,           /**< the constraint handler itself */
      SCIP_CONS**        conss,              /**< array of constraints to process */
      int                nconss,             /**< number of constraints to process */
      int                nusefulconss,       /**< number of useful (non-obsolete) constraints to process */
      SCIP_Bool          solinfeasible,      /**< was the solution already declared infeasible by a constraint handler? */
      SCIP_Bool          objinfeasible,      /**< is the solution infeasible anyway due to violating lower objective bound? */
      SCIP_RESULT*       result              /**< pointer to store the result of the enforcing call */
      );


   /** constraint enforcing method of constraint handler for LP solutions
    *
    *  The method is called at the end of the node processing loop for a node where the LP was solved.
    *  The LP solution has to be checked for feasibility. If possible, an infeasibility should be resolved by
    *  branching, reducing a variable's domain to exclude the solution or separating the solution with a valid
    *  cutting plane.
    *
    *  The enforcing methods of the active constraint handlers are called in decreasing order of their enforcing
    *  priorities until the first constraint handler returned with the value SCIP_CUTOFF, SCIP_SEPARATED,
    *  SCIP_REDUCEDDOM, SCIP_CONSADDED, or SCIP_BRANCHED.
    *  The integrality constraint handler has an enforcing priority of zero. A constraint handler which can
    *  (or wants) to enforce its constraints only for integral solutions should have a negative enforcing priority
    *  (e.g. the alldiff-constraint can only operate on integral solutions).
    *  A constraint handler which wants to incorporate its own branching strategy even on non-integral
    *  solutions must have an enforcing priority greater than zero (e.g. the SOS-constraint incorporates
    *  SOS-branching on non-integral solutions).
    *
    *  The first nusefulconss constraints are the ones, that are identified to likely be violated. The enforcing
    *  method should process the useful constraints first. The other nconss - nusefulconss constraints should only
    *  be enforced, if no violation was found in the useful constraints.
    *
    *  possible return values for *result (if more than one applies, the first in the list should be used):
    *  - SCIP_CUTOFF     : the node is infeasible in the variable's bounds and can be cut off
    *  - SCIP_CONSADDED  : an additional constraint was generated
    *  - SCIP_REDUCEDDOM : a variable's domain was reduced
    *  - SCIP_SEPARATED  : a cutting plane was generated
    *  - SCIP_BRANCHED   : no changes were made to the problem, but a branching was applied to resolve an infeasibility
    *  - SCIP_INFEASIBLE : at least one constraint is infeasible, but it was not resolved
    *  - SCIP_FEASIBLE   : all constraints of the handler are feasible
    */

   virtual SCIP_RETCODE scip_enfolp(
      SCIP*              scip,               /**< SCIP data structure */
      SCIP_CONSHDLR*     conshdlr,           /**< the constraint handler itself */
      SCIP_CONS**        conss,              /**< array of constraints to process */
      int                nconss,             /**< number of constraints to process */
      int                nusefulconss,       /**< number of useful (non-obsolete) constraints to process */
      SCIP_Bool          solinfeasible,      /**< was the solution already declared infeasible by a constraint handler? */
      SCIP_RESULT*       result              /**< pointer to store the result of the enforcing call */
      );


   /** separation method of constraint handler for arbitrary primal solution
    *
    *  Separates all constraints of the constraint handler. The method is called outside the LP solution loop (e.g., by
    *  a relaxator or a primal heuristic), which means that there is no valid LP solution.
    *  Instead, the method should produce cuts that separate the given solution.
    *
    *  The first nusefulconss constraints are the ones, that are identified to likely be violated. The separation
    *  method should process only the useful constraints in most runs, and only occasionally the remaining
    *  nconss - nusefulconss constraints.
    *
    *  possible return values for *result (if more than one applies, the first in the list should be used):
    *  - SCIP_CUTOFF     : the node is infeasible in the variable's bounds and can be cut off
    *  - SCIP_CONSADDED  : an additional constraint was generated
    *  - SCIP_REDUCEDDOM : a variable's domain was reduced
    *  - SCIP_SEPARATED  : a cutting plane was generated
    *  - SCIP_DIDNOTFIND : the separator searched, but did not find domain reductions, cutting planes, or cut constraints
    *  - SCIP_DIDNOTRUN  : the separator was skipped
    *  - SCIP_DELAYED    : the separator was skipped, but should be called again
    */
   virtual SCIP_RETCODE scip_sepasol(
      SCIP*              scip,               /**< SCIP data structure */
      SCIP_CONSHDLR*     conshdlr,           /**< the constraint handler itself */
      SCIP_CONS**        conss,              /**< array of constraints to process */
      int                nconss,             /**< number of constraints to process */
      int                nusefulconss,       /**< number of useful (non-obsolete) constraints to process */
      SCIP_SOL*          sol,                /**< primal solution that should be separated */
      SCIP_RESULT*       result              /**< pointer to store the result of the separation call */
      );


   /** separation method of constraint handler for LP solution
    *
    *  Separates all constraints of the constraint handler. The method is called in the LP solution loop,
    *  which means that a valid LP solution exists.
    *
    *  The first nusefulconss constraints are the ones that are identified to likely be violated. The separation
    *  method should process only the useful constraints in most runs, and only occasionally the remaining
    *  nconss - nusefulconss constraints.
    *
    *  possible return values for *result (if more than one applies, the first in the list should be used):
    *  - SCIP_CUTOFF     : the node is infeasible in the variable's bounds and can be cut off
    *  - SCIP_CONSADDED  : an additional constraint was generated
    *  - SCIP_REDUCEDDOM : a variable's domain was reduced
    *  - SCIP_SEPARATED  : a cutting plane was generated
    *  - SCIP_DIDNOTFIND : the separator searched, but did not find domain reductions, cutting planes, or cut constraints
    *  - SCIP_DIDNOTRUN  : the separator was skipped
    *  - SCIP_DELAYED    : the separator was skipped, but should be called again
    */
   virtual SCIP_RETCODE scip_sepalp(
      SCIP*              scip,               /**< SCIP data structure */
      SCIP_CONSHDLR*     conshdlr,           /**< the constraint handler itself */
      SCIP_CONS**        conss,              /**< array of constraints to process */
      int                nconss,             /**< number of constraints to process */
      int                nusefulconss,       /**< number of useful (non-obsolete) constraints to process */
      SCIP_RESULT*       result              /**< pointer to store the result of the separation call */
      );

   /** initialization method of constraint handler (called after problem has been transformed) */

   /**delete method of constraint handler*/
   virtual SCIP_RETCODE scip_delete(
      SCIP*              scip,               /**< SCIP data structure */
      SCIP_CONSHDLR*     conshdlr,           /**< the constraint handler itself */
      SCIP_CONS*         cons,               /**< the constraint belonging to the constraint data */
      SCIP_CONSDATA**    consdata            /**< pointer to the constraint data to free */
      );

};

/**gets (a pointer to) the sdpcone of a sdp-constraint*/
SCIP_RETCODE getSdpCone(
   SCIP* scip,             /**<SCIP data structure*/
   SCIP_CONS* conss,       /**<constraint to get data for*/
   SdpCone** sdpcone       /**<output: pointer to the cone*/);



/**Method for creating sdp constraints*/
SCIP_RETCODE SCIPcreateConsSdp(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           cons,               /**< pointer to hold the created constraint */
   const char*           name,               /**< name of constraint */
   const SdpCone&        sdpcone             /**< the sdpcone*/

   );

#endif
