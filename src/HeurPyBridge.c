/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the                                 */
/*      SDP-Package for SCIP: a solving framework for                        */
/*                            mixed-integer semidefinite programms           */
/*                                                                           */
/* Copyright (C) 2011-2014 Discrete Optimization, TU Darmstadt               */
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

/**@file   heur_xyz.c
 * @brief  xyz primal heuristic
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/


#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION

#include <assert.h>
#include <python2.7/Python.h>
#include <numpy/arrayobject.h>
#include "HeurPyBridge.h"


#define HEUR_NAME             "PyBridge"
#define HEUR_DESC             "interface to python heuristics"
#define HEUR_DISPCHAR         'P'
#define HEUR_PRIORITY         0
#define HEUR_FREQ             8
#define HEUR_FREQOFS          0
#define HEUR_MAXDEPTH         -1
#define HEUR_TIMING           SCIP_HEURTIMING_AFTERNODE
#define HEUR_USESSUBSCIP      FALSE  /**< does the heuristic use a secondary SCIP instance? */


SCIP_RETCODE callPyHeuristic(SCIP *scip, SCIP_HEURDATA *heurdata, SCIP_RESULT* result, SCIP_HEUR *heur);
PyObject *createPyEnvironment(char *filename);
void destroyPyEnvironment(PyObject * pHeurObj);
int8_t integerVarBounds(SCIP *scip, SCIP_VAR*);

/** primal heuristic data */
struct SCIP_HeurData
{
   SCIP_VAR **orig_vars;
   int num_orig_vars;
   PyObject *pyHeurObj;
   int num_bars;
   int num_bar_areas;
   int num_scenarios;
   int num_actors;
   int end_scenario_obj;
   int end_bars;
   int end_actor_str;
   int end_actor_pos;
   PyArrayObject *bar_mask;
   PyArrayObject *actor_mask;
   PyArrayObject *actor_bounds;
   PyArrayObject *scip_sol;
   PyArrayObject *relaxation;
   int runcounter;
};


/*
 * Local methods
 */


/* put your local methods here, and declare them static */
PyObject *createPyEnvironment(char *filename) {
   // start python and create the corresponding heuristic-object
   PyObject *pModuleName, *pModule, *pArgs, *pHeurObj, *pFunc;

   // initialize python
   Py_Initialize();
//   import_array();

//   PyRun_SimpleString(
//   "import sys\n"
//   "sys.path.insert(0, '/home/dassmann/projects/opti_hiwi/scipdsdp/Scipdsdp/src')\n");
   // PyRun_SimpleString("print sys.path\n");

   // get the names of the module / class to call
   pModuleName = PyString_FromString("truss_heuristics.scip_plugin.heurpybridge");

   // import module
   pModule = PyImport_Import(pModuleName);
   Py_DECREF(pModuleName);
   if (pModule == NULL){
      PyErr_Print();
      fprintf(stderr, "Failed to load heurpybridge.py\n" );
      SCIPABORT();
      return NULL;
   }

   // assume that module is loaded correctly, proceed with class
   pFunc = PyObject_GetAttrString(pModule, "Heuristic");

   if (!(pFunc && PyCallable_Check(pFunc))) { // error: nothing usable returned
      if (PyErr_Occurred())
         PyErr_Print();
      fprintf(stderr, "Cannot find class 'Heuristic'");
      SCIPABORT();
      return NULL;
   }

   // create the object
   pArgs = Py_BuildValue("(s)", filename);
   pHeurObj = PyObject_CallObject(pFunc, pArgs);
   Py_XDECREF(pArgs);
   if (!pHeurObj) {
      if (PyErr_Occurred())
         PyErr_Print();
      Py_DECREF(pModule);
      fprintf(stderr, "Could not instantiate heuristic object\n");
      SCIPABORT();
      return NULL;
   }

   Py_XDECREF(pFunc);
   Py_DECREF(pModule);

   return pHeurObj;
}


int8_t integerVarBounds(SCIP *scip, SCIP_VAR *var){
   // check lb & ub of var.
   int8_t value;
   // 0: lb=ub=0
   // 1: lb=ub=1
   // -1: lb=0, ub=1
   SCIP_Real lb = SCIPvarGetLbLocal(var);
   SCIP_Real ub = SCIPvarGetUbLocal(var);
   if(SCIPisEQ(scip, lb, 0) && SCIPisEQ(scip, ub, 0)) {
      value = 0;
   } else if(SCIPisEQ(scip, lb, 1) && SCIPisEQ(scip, ub, 1)) {
      value = 1;
   } else if(SCIPisEQ(scip, lb, 0) && SCIPisEQ(scip, ub, 1)) {
      value = -1;
   } else {
      SCIPABORT(); // "bounds not as excepted.. abort"
      value = -10;
   }
   return value;
}


SCIP_RETCODE callPyHeuristic(SCIP *scip, SCIP_HEURDATA* heurdata, SCIP_RESULT* result, SCIP_HEUR *heur){
   int i;
   SCIP_Bool valid_relaxation;
   int bar_mask_counter, actor_pos_counter, actor_str_counter;
   SCIP_VAR* var;
   SCIP_SOL* scip_sol;
   PyObject *ret, *has_valid_relaxation;
   PyArrayObject *bar_mask, *actor_mask, *actor_bounds, *relaxation;
   *result = SCIP_DIDNOTFIND;
//   PyObject * pyHeurObj = heurdata->pyHeurObj;

   bar_mask = heurdata->bar_mask;
   actor_mask = heurdata->actor_mask;
   actor_bounds = heurdata->actor_bounds;
   relaxation = heurdata->relaxation;
   valid_relaxation = SCIPisRelaxSolValid(scip);

   assert(PyArray_ISBEHAVED(bar_mask));
   assert(PyArray_ISBEHAVED(relaxation));
   if(heurdata->num_actors > 0){
      assert(PyArray_ISBEHAVED(actor_mask));
      assert(PyArray_ISBEHAVED(actor_bounds));
   }

   bar_mask_counter = heurdata->num_bars * heurdata->num_bar_areas;
   if(heurdata->num_actors > 0){
      actor_pos_counter = heurdata->num_bars;
      actor_str_counter = heurdata->num_bars * heurdata->num_scenarios;
   } else {
      actor_pos_counter = actor_str_counter = 0;
   }


   // loop over objective
   for (i = 0; i < 1; i++){
      var = heurdata->orig_vars[i];
      if(valid_relaxation){
         *((double *) PyArray_GETPTR1(relaxation, i)) = (double) SCIPgetRelaxSolVal(scip, var);
      }
   }

   // loop over scenarios
   for (i = 1; i < heurdata->end_scenario_obj; i++){
      var = heurdata->orig_vars[i];
      if(valid_relaxation){
         *((double *) PyArray_GETPTR1(relaxation, i)) = (double) SCIPgetRelaxSolVal(scip, var);
      }
   }

   // loop over bars
   for (i = heurdata->end_scenario_obj; i < heurdata->end_bars; i++){
      int bar, bar_area, value;
      var = heurdata->orig_vars[i];
      if(valid_relaxation){
         *((double *) PyArray_GETPTR1(relaxation, i)) = (double) SCIPgetRelaxSolVal(scip, var);
      }

      bar_area = (i - heurdata->end_scenario_obj) % heurdata->num_bar_areas;
      bar = ((i - heurdata->end_scenario_obj) - bar_area) / heurdata->num_bar_areas;
      value = integerVarBounds(scip, var);
      *((int8_t *) PyArray_GETPTR2(bar_mask, bar, bar_area)) = value;
      --bar_mask_counter;
   }

   // loop over actor strengths
   for (i = heurdata->end_bars; i < heurdata->end_actor_str; i++){
      int bar, scenario;
      double lb, ub;
      var = heurdata->orig_vars[i];
      if(valid_relaxation){
         *((double *) PyArray_GETPTR1(relaxation, i)) = (double) SCIPgetRelaxSolVal(scip, var);
      }

      scenario = (i - heurdata->end_bars) % heurdata->num_scenarios;
      bar = ((i - heurdata->end_bars) - scenario) / heurdata->num_scenarios;
      lb = (double) SCIPvarGetLbLocal(var);
      ub = (double) SCIPvarGetUbLocal(var);

      *((double *) PyArray_GETPTR3(actor_bounds, bar, scenario, 0)) = (double) lb;
      *((double *) PyArray_GETPTR3(actor_bounds, bar, scenario, 1)) = (double) ub;
      --actor_str_counter;
   }

   // loop over actor pos
   for (i = heurdata->end_actor_str; i < heurdata->end_actor_pos; i++){
      int value;
      int bar = i - heurdata->end_actor_str;
      var = heurdata->orig_vars[i];
      if(valid_relaxation){
         *((double *) PyArray_GETPTR1(relaxation, i)) = (double) SCIPgetRelaxSolVal(scip, var);
      }

      value = integerVarBounds(scip, var);
      *((int8_t *) PyArray_GETPTR1(actor_mask, bar)) = value;
      --actor_pos_counter;
   }

   assert(bar_mask_counter == 0);
   assert(actor_pos_counter == 0);
   assert(actor_str_counter == 0);

   SCIPdebugMessage("calling python heuristic\n");
   if(valid_relaxation){
      has_valid_relaxation = Py_True;
   } else {
      has_valid_relaxation = Py_False;
   }

   ret = PyObject_CallMethod(heurdata->pyHeurObj, (char *) "run_heuristic", (char *) "(i, O)", SCIPgetDepth(scip), has_valid_relaxation, NULL);
   if(ret == NULL){
      if (PyErr_Occurred())
         PyErr_Print();
      fprintf(stderr, "error occured, nothing returned from python heuristic\n");
   }
   if(ret == Py_True){
      SCIP_Bool stored;
      SCIP_CALL( SCIPcreateSol(scip, &scip_sol, heur) );
      SCIPdebugMessage("reading python heuristic solution\n");

      for(i = 0; i < heurdata->num_orig_vars; i++){
          SCIP_Real v = *((SCIP_Real *) PyArray_GETPTR1(heurdata->scip_sol, i));
          // SCIPdebugMessage("var %d: %f\n", i, v);
          SCIP_CALL( SCIPsetSolVal(scip, scip_sol, heurdata->orig_vars[i], v) );
      }
      //SCIP_CALL( SCIPprintSol(scip, scip_sol, NULL, FALSE) );

      SCIP_CALL(SCIPtrySol(scip, scip_sol, TRUE, FALSE, FALSE, TRUE, &stored));
      if(stored){
         SCIPdebugMessage("found feasible solution\n");
         // SCIP_CALL( SCIPprintSol(scip, scip_sol, NULL, FALSE) );
         *result = SCIP_FOUNDSOL;
      }
      SCIP_CALL( SCIPfreeSol(scip, &scip_sol) );

   }
   Py_XDECREF(ret);
   return SCIP_OKAY;
}


void destroyPyEnvironment(PyObject * pHeurObj){
   // dereference stuff
   Py_DECREF(pHeurObj);
   // end python
   Py_Finalize();
}


/*
 * Callback methods of primal heuristic
 */

/* TODO: Implement all necessary primal heuristic methods. The methods with an #if 0 ... #else #define ... are optional */

/** copy method for primal heuristic plugins (called when SCIP copies plugins) */
#if 0
static
SCIP_DECL_HEURCOPY(heurCopyPyBridge)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of xyz primal heuristic not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define heurCopyPyBridge NULL
#endif

/** destructor of primal heuristic to free user data (called when SCIP is exiting) */
#if 1
static
SCIP_DECL_HEURFREE(heurFreePyBridge)
{  /*lint --e{715}*/
   SCIP_HEURDATA* heurdata;
   assert(heur != NULL);
   assert(strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0);
   assert(scip != NULL);

   heurdata = SCIPheurGetData(heur);

   /* free heuristic data */
   assert(heurdata != NULL);

   /* stop python */
   destroyPyEnvironment(heurdata->pyHeurObj);
   SCIPfreeMemoryArray(scip, &(heurdata->orig_vars));
   SCIPfreeMemory(scip, &heurdata);
   SCIPheurSetData(heur, NULL);
   return SCIP_OKAY;
}
#else
#define heurFreePyBridge NULL
#endif


/** initialization method of primal heuristic (called after problem was transformed) */
#if 1
static
SCIP_DECL_HEURINIT(heurInitPyBridge)
{  /*lint --e{715}*/
   PyObject *ret;
   PyArrayObject *bar_mask, *actor_mask, *actor_bounds, *scip_sol, *relaxation;
   int i, num_orig_vars;
   SCIP_VAR **scip_orig_vars;
   char * probname = (char *) SCIPgetProbName(scip);
   SCIP_HEURDATA* heurdata = SCIPheurGetData(heur);
   heurdata->pyHeurObj = createPyEnvironment(probname);

   ret = PyObject_CallMethod(heurdata->pyHeurObj, (char *) "get_var_ranges", (char *) "()", NULL);
   if(ret == NULL || ret == Py_None){
      fprintf(stderr, "error while getting var_ranges");
   }
   // returns: num_bars, num_bar_areas, num_scenarios, num_actors, obj_end, bars_end, actors_end
   PyArg_ParseTuple(ret, "iiiiiiii", &(heurdata->num_bars),
                                    &(heurdata->num_bar_areas),
                                    &(heurdata->num_scenarios),
                                    &(heurdata->num_actors),
                                    &(heurdata->end_scenario_obj),
                                    &(heurdata->end_bars),
                                    &(heurdata->end_actor_str),
                                    &(heurdata->end_actor_pos));

   // maybe this could be done better, but we need the original variables in their original order
   SCIP_CALL(SCIPgetOrigVarsData(scip, &scip_orig_vars, &num_orig_vars, NULL, NULL, NULL, NULL));
   SCIP_CALL(SCIPallocMemoryArray(scip, &(heurdata->orig_vars), num_orig_vars));
   assert(heurdata->orig_vars != NULL);
   heurdata->runcounter = 0;
   for (i = 0; i < num_orig_vars; ++i){
      int n = SCIPvarGetIndex(scip_orig_vars[i]);
      heurdata->orig_vars[n] = scip_orig_vars[i];
   }
   heurdata->num_orig_vars = num_orig_vars;

   bar_mask = (PyArrayObject *) PyObject_GetAttrString(heurdata->pyHeurObj, "bar_mask");
   if(bar_mask == NULL || (PyObject *) bar_mask == Py_None){
      fprintf(stderr, "bar_mask == None or no bar_mask array found.");
      SCIPABORT();
   }
   heurdata->bar_mask = bar_mask;


   actor_mask = (PyArrayObject *) PyObject_GetAttrString(heurdata->pyHeurObj, "actor_mask");
   if(actor_mask == NULL){
      fprintf(stderr, "no actor_mask array found.");
      SCIPABORT();
   }
   heurdata->actor_mask = actor_mask;


   actor_bounds = (PyArrayObject *) PyObject_GetAttrString(heurdata->pyHeurObj, "actor_bounds");
   if(actor_bounds == NULL){
      fprintf(stderr, "no actor_bounds array found.");
      SCIPABORT();
   }
   heurdata->actor_bounds = actor_bounds;

   scip_sol = (PyArrayObject *) PyObject_GetAttrString(heurdata->pyHeurObj, "solution");
   if(scip_sol == NULL){
      fprintf(stderr, "no solution array found.");
      SCIPABORT();
   }
   heurdata->scip_sol = scip_sol;

   relaxation = (PyArrayObject *) PyObject_GetAttrString(heurdata->pyHeurObj, "relaxation");
   if(relaxation == NULL){
      fprintf(stderr, "no solution array found.");
      SCIPABORT();
   }
   heurdata->relaxation = relaxation;

   return SCIP_OKAY;
}
#else
#define heurInitPyBridge NULL
#endif


/** deinitialization method of primal heuristic (called before transformed problem is freed) */
#if 0
static
SCIP_DECL_HEUREXIT(heurExitPyBridge)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of xyz primal heuristic not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define heurExitPyBridge NULL
#endif


/** solving process initialization method of primal heuristic (called when branch and bound process is about to begin) */
#if 0
static
SCIP_DECL_HEURINITSOL(heurInitsolPyBridge)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of xyz primal heuristic not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define heurInitsolPyBridge NULL
#endif


/** solving process deinitialization method of primal heuristic (called before branch and bound process data is freed) */
#if 0
static
SCIP_DECL_HEUREXITSOL(heurExitsolPyBridge)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of xyz primal heuristic not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define heurExitsolPyBridge NULL
#endif


/** execution method of primal heuristic */
static
SCIP_DECL_HEUREXEC(heurExecPyBridge)
{  /*lint --e{715}*/
   SCIP_HEURDATA* heurdata;
   *result = SCIP_DIDNOTFIND;

   heurdata = SCIPheurGetData(heur);
   if(heurdata->runcounter < 0){
      return SCIP_OKAY;
   }
   ++heurdata->runcounter;
   SCIP_CALL(callPyHeuristic(scip, heurdata, result, heur));

   return SCIP_OKAY;
}


/*
 * primal heuristic specific interface methods
 */

/** creates the xyz primal heuristic and includes it in SCIP */
SCIP_RETCODE SCIPincludeHeurPyBridge(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_HEURDATA* heurdata;
   SCIP_HEUR* heur;

   /* create xyz primal heuristic data */
   SCIP_CALL( SCIPallocMemory(scip, &heurdata) );

   heur = NULL;

   /* include primal heuristic */
#if 0
   /* use SCIPincludeHeur() if you want to set all callbacks explicitly and realize (by getting compiler errors) when
    * new callbacks are added in future SCIP versions
    */
   SCIP_CALL( SCIPincludeHeur(scip, HEUR_NAME, HEUR_DESC, HEUR_DISPCHAR, HEUR_PRIORITY, HEUR_FREQ, HEUR_FREQOFS,
         HEUR_MAXDEPTH, HEUR_TIMING, HEUR_USESSUBSCIP,
         heurCopyPyBridge, heurFreePyBridge, heurInitPyBridge, heurExitPyBridge, heurInitsolPyBridge, heurExitsolPyBridge, heurExecPyBridge,
         heurdata) );
#else
   /* use SCIPincludeHeurBasic() plus setter functions if you want to set callbacks one-by-one and your code should
    * compile independent of new callbacks being added in future SCIP versions
    */
   SCIP_CALL( SCIPincludeHeurBasic(scip, &heur,
         HEUR_NAME, HEUR_DESC, HEUR_DISPCHAR, HEUR_PRIORITY, HEUR_FREQ, HEUR_FREQOFS,
         HEUR_MAXDEPTH, HEUR_TIMING, HEUR_USESSUBSCIP, heurExecPyBridge, heurdata) );

   assert(heur != NULL);

   /* set non fundamental callbacks via setter functions */
   SCIP_CALL( SCIPsetHeurCopy(scip, heur, heurCopyPyBridge) );
   SCIP_CALL( SCIPsetHeurFree(scip, heur, heurFreePyBridge) );
   SCIP_CALL( SCIPsetHeurInit(scip, heur, heurInitPyBridge) );
   SCIP_CALL( SCIPsetHeurExit(scip, heur, heurExitPyBridge) );
   SCIP_CALL( SCIPsetHeurInitsol(scip, heur, heurInitsolPyBridge) );
   SCIP_CALL( SCIPsetHeurExitsol(scip, heur, heurExitsolPyBridge) );
#endif

   /* add xyz primal heuristic parameters */
   /* TODO: (optional) add primal heuristic specific parameters with SCIPaddTypeParam() here */

   return SCIP_OKAY;
}
