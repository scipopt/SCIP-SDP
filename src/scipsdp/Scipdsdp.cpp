
/**@file   Scipdsdp.cpp
 * @brief  driver-file for solving MISDPs
 * @author Sonja Mars
 */
// standard library includes
#include <string>

#include "objscip/objscip.h"
#include "objscip/objscipdefplugins.h"

#include "objconshdlr_sdp.h"
#include "objrelax_sdp.h"
#include "objreader_sdpa.h"
//#include "objheur_sdp.h"
//#include "objpresol_binsdp.h"

using namespace scip;

/**run scip and set some parameters*/
static
SCIP_RETCODE runSCIP(int argc,
   char** argv)
{

   printf("Starting solver and do something.\n");
   SCIP* scip = NULL;


   printf("\n");

   SCIP_CALL( SCIPcreate(&scip) );
   SCIPprintVersion(scip, NULL);
   
   //int relax_freq = 0;

   SCIP_CALL( SCIPincludeObjReader(scip, new ObjReaderSDPA(scip), TRUE) );

   SCIP_CALL( SCIPincludeObjConshdlr(scip, new ObjConshdlrSdp(scip), TRUE) );

   //const char* nameint = "relax_freq";
   //const char* descint = "frequency how often the sdprelaxator is called";
   //int freq = 0;
   //SCIP_CALL( SCIPaddIntParam	(scip, nameint, descint, NULL, FALSE, freq, -1, 1000000, NULL, NULL));

   SCIP_CALL( SCIPincludeObjRelax(scip, new ObjRelaxSdp(scip), TRUE) );

   
   const char* name = "sdpsolver";
   const char * 	desc = "which sdpsolver should be called";
   //char** 	valueptr;
   //SCIP_Bool 	isadvanced;
  // const char* 	defaultvalue = "pensdp";
   SCIP_PARAMDATA* 	paramdata = NULL;

   SCIP_CALL( SCIPaddStringParam	(scip, name, desc, NULL, FALSE, "dsdp" , NULL, paramdata)	);	

   /* include default SCIP plugins */
   SCIP_CALL( SCIPincludeDefaultPlugins(scip) );
   
   SCIP_CALL(SCIPsetIntParam(scip,"relaxing/SDPRelax/freq",1));
   

   /**********************************
    * Process command line arguments *
    **********************************/
   // SCIP_CALL( SCIPprocessShellArguments(scip, argc, argv, "sciptsp.set") );
   SCIP_CALL( SCIPsetIntParam(scip, "display/verblevel", 5) );

   //Turn off lp relaxations
   SCIP_CALL( SCIPsetIntParam(scip, "lp/solvefreq", -1) );
   //relaxierer frequency als parameter einstellen

   //Do some stuff to be numerically stable
   SCIP_CALL( SCIPsetRealParam(scip, "numerics/epsilon", 1e-6) );//normal 10-9
   //SCIP_CALL( SCIPsetRealParam(scip, "numerics/sumepsilon", 1e-4) );//normal 10-6
   SCIP_CALL( SCIPsetRealParam(scip, "numerics/feastol", 1e-4) );//normal 10-6
   
   SCIP_CALL( SCIPsetStringParam(scip, "sdpsolver", "dsdp") );

   //SCIP_CALL( SCIPsetIntParam(scip, "heuristics/trivial/freq", -1));

   //SCIP_CALL( SCIPsetIntParam(scip, "propagating/maxrounds", 0));
   //SCIP_CALL( SCIPsetIntParam(scip, "presolving/maxrestarts", 0));

   //SCIP_CALL( SCIPsetLongintParam(scip, "limits/nodes",20));
   //SCIP_CALL( SCIPsetRealParam(scip, "limits/time",7200) );
   
   //SCIP_CALL( SCIPsetRealParam(scip, "separating/minefficacy ",0.001));kennt er nicht
   //SCIP_CALL( SCIPsetRealParam(scip, "separating/minefficacyroot",0.01));
   
   SCIP_CALL( SCIPsetBoolParam(scip, "lp/cleanuprows", FALSE) );
   SCIP_CALL( SCIPsetBoolParam(scip, "lp/cleanuprowsroot", FALSE) );
   SCIP_CALL( SCIPsetIntParam(scip, "lp/rowagelimit", 10) );
   
   //# maximum age a cut can reach before it is deleted from the global cut pool, or -1 to keep all cuts
   //# [type: int, range: [-1,2147483647], default: 100]
   SCIP_CALL( SCIPsetIntParam(scip, "separating/cutagelimit", 10));
   
   
   SCIP_CALL( SCIPsetIntParam(scip, "separating/maxrounds", 20));
   
   //# maximal number of consecutive separation rounds without objective or integrality improvement (-1: no additional restriction)
   //# [type: int, range: [-1,2147483647], default: 5]
   //SCIP_CALL( SCIPsetIntParam(scip, "separating/maxstallrounds", 15));
   

   
   //# separation frequency for the global cut pool (-1: disable global cut pool, 0: only separate pool at the root)
   //# [type: int, range: [-1,2147483647], default: 0]
//SCIP_CALL( SCIPsetIntParam(scip, "separating/poolfreq", 1));
   
   //# frequency for separating cuts (-1: never, 0: only in root node)
   //# [type: int, range: [-1,2147483647], default: 0]
   //SCIP_CALL( SCIPsetIntParam(scip, "constraints/linear/sepafreq", 0));
   //SCIP_CALL( SCIPsetIntParam(scip, "constraints/knapsack/sepafreq", 0));
   
   //# frequency for using all instead of only the useful constraints in separation, propagation and enforcement (-1: never, 0: only in first evaluation)
   //# [type: int, range: [-1,2147483647], default: 100]
   //SCIP_CALL( SCIPsetIntParam(scip, "constraints/linear/eagerfreq", 5));
   
   //# maximal number of separation rounds per node (-1: unlimited)
   //# [type: int, range: [-1,2147483647], default: 5]
   //SCIP_CALL( SCIPsetIntParam(scip, "constraints/linear/maxrounds",15));
   
   //SCIP_CALL( SCIPsetIntParam(scip, "nodeselection/dfs/stdpriority", 1000000));
   SCIP_CALL( SCIPsetIntParam(scip, "nodeselection/hybridestim/stdpriority", 1000000));
   SCIP_CALL( SCIPsetIntParam(scip, "nodeselection/hybridestim/maxplungedepth", 0));
   SCIP_CALL( SCIPsetRealParam(scip, "nodeselection/hybridestim/estimweight",0));
       
   SCIP_CALL( SCIPsetIntParam(scip, "branching/pscost/priority",-2000000));
   SCIP_CALL( SCIPsetIntParam(scip, "branching/relpscost/priority",-2000000));
   
   //# should presolving try to simplify inequalities
   //# [type: bool, range: {TRUE,FALSE}, default: FALSE]
   //SCIP_CALL( SCIPsetBoolParam(scip, "constraints/linear/simplifyinequalities", TRUE));
   
   //# frequency for calling separator <cmir> (-1: never, 0: only in root node)
   //# [type: int, range: [-1,2147483647], default: 0]
   //SCIP_CALL( SCIPsetIntParam(scip, "separating/cmir/freq", 0));
   //SCIP_CALL( SCIPsetIntParam(scip, "separating/flowcover/freq", 0));
   SCIP_CALL( SCIPsetIntParam(scip, "separating/intobj/freq", -1));
   //SCIP_CALL( SCIPsetIntParam(scip, "separating/closecuts/freq", 1));
   //SCIP_CALL( SCIPsetIntParam(scip, "separating/rapidlearning/freq", 0));
   //SCIP_CALL( SCIPsetBoolParam(scip, "separating/closecuts/separelint", TRUE));
   //SCIP_CALL( SCIPsetBoolParam(scip, "separating/closecuts/inclobjcutoff", TRUE));
   //SCIP_CALL( SCIPsetIntParam(scip, "separating/closecuts/maxunsuccessful", 10));
   //SCIP_CALL( SCIPsetRealParam(scip, "separating/closecuts/sepacombvalue", 0.5));
    //SCIP_CALL( SCIPsetBoolParam(scip, "separating/closecuts/recomputerelint", TRUE));

//   SCIP_CALL( SCIPsetIntParam(scip, "separating/mcf/freq", -1));
//   SCIP_CALL( SCIPsetIntParam(scip, "separating/zerohalf/freq", -1));
//   SCIP_CALL( SCIPsetIntParam(scip, "separating/clique/freq", -1));
//   SCIP_CALL( SCIPsetIntParam(scip, "separating/impliedbounds/freq", -1));
//   SCIP_CALL( SCIPsetIntParam(scip, "separating/cgmip/freq", -1));
//   SCIP_CALL( SCIPsetIntParam(scip, "separating/gomory/freq", -1));
//   SCIP_CALL( SCIPsetIntParam(scip, "separating/strongcg/freq", -1));
//   SCIP_CALL( SCIPsetIntParam(scip, "separating/oddcycle/freq", -1));


   //# priority of heuristic <simplerounding>
   //# [type: int, range: [-536870912,536870911], default: 0]
   //SCIP_CALL( SCIPsetIntParam(scip, "heuristics/simplerounding/priority",12000));
   //
   //# frequency for calling primal heuristic <simplerounding> (-1: never, 0: only at depth freqofs)
   //# [type: int, range: [-1,2147483647], default: 1]
   //SCIP_CALL( SCIPsetIntParam(scip, "heuristics/simplerounding/freq",1));
   //
   //# frequency offset for calling primal heuristic <simplerounding>
   //# [type: int, range: [0,2147483647], default: 0]
   //SCIP_CALL( SCIPsetIntParam(scip, "heuristics/simplerounding/freqofs",0));
   if (argc > 2 )
   {
      SCIP_CALL( SCIPreadParams(scip,argv[2]));
   }

 
   

   printf("\n read problem\n");
   printf("============\n");

   if (argc > 1)
   {
      SCIP_CALL( SCIPreadProb(scip, argv[1], NULL) );
   }
   else
   {
      SCIP_CALL( SCIPreadProb(scip, "/Users/sonja/EDOM/Scipdsdp/scipdsdp/Scipdsdp/data_for_testing/achtziger_stolpe06-6.1.dat-s", NULL));
   }
  // if (argc > 2 )
   //{
     // SCIP_CALL(SCIPreadProb(scip, argv[2], NULL));
   //}


   /* solve problem */
   printf("\nsolve problem\n");
   printf("=============\n");
   SCIP_CALL( SCIPsolve(scip) );

   SCIP_CALL( SCIPprintStatistics(scip, NULL));
   printf("\nprimal solution:\n" );
   printf("================\n\n");
   SCIP_CALL( SCIPprintBestSol(scip, NULL, FALSE) );

   FILE* solfile;
   std::string solfilename = "Solution.sol";
   if(argc > 1)
   {
      solfilename = std::string(argv[1]) + "_misdp.sol";
   }


   solfile = fopen(solfilename.c_str(), "w");
   if( solfile == NULL )
   {
      SCIPerrorMessage("error creating file <%s>\n", solfilename.c_str());
   }
   else
   {
      SCIP_CALL( SCIPprintBestSol(scip, solfile, FALSE) );

   }
   fclose(solfile);

   /********************
    * Deinitialization *
    ********************/
   SCIP_CALL( SCIPfree(&scip) );

   BMScheckEmptyMemory();

   return SCIP_OKAY;
}

/**main function */
int main (int argc,
   char** argv)
{
   SCIP_RETCODE retcode;

   retcode = runSCIP(argc, argv);
   if( retcode != SCIP_OKAY )
   {
    //  SCIPprintError(retcode);
      return -1;
   }

   return 0;
}
