using MathProgBase
using GZip
using Pajarito
using Mosek
using CPLEX
using ConicBenchmarkUtilities


function loadcbf(solver, filename::String; all_bin=false)
    dat = readcbfdata(filename) # .cbf.gz extension also accepted
    c, A, b, con_cones, var_cones, vartypes, sense, objoffset = cbftompb(dat)
    model = MathProgBase.ConicModel(solver)
    MathProgBase.loadproblem!(model, c, A, b, con_cones, var_cones)
    MathProgBase.setvartype!(model, vartypes)

    return model
end



Filename="/home/gally/MATLAB_libraries/SCIPSDPpaper_CBF_drives_Pajarito060_Mosek81054_CPLEX1261_tolerances.test"
f = open(Filename)
instances = readlines(f)
close(f)

for instance in instances
        instancename = normalize_string(instance, stripcc=true)
        print("$instancename\n")
        #model = loadcbf(PajaritoSolver(mip_solver=CplexSolver(CPX_PARAM_THREADS=1), cont_solver=MosekSolver(LOG=0), mip_solver_drives=true, dualize_relax=true, dualize_subp=true, timeout=3600), strip(instancename))
        model = loadcbf(PajaritoSolver(mip_solver=CplexSolver(CPX_PARAM_THREADS=1, CPX_PARAM_EPINT=1e-6, CPX_PARAM_EPRHS=1e-6), cont_solver=MosekSolver(LOG=0,MSK_DPAR_INTPNT_CO_TOL_DFEAS=1e-7, MSK_DPAR_INTPNT_CO_TOL_INFEAS=1e-7, MSK_DPAR_INTPNT_CO_TOL_PFEAS=1e-5, MSK_DPAR_INTPNT_CO_TOL_MU_RED=1e-5, MSK_DPAR_INTPNT_CO_TOL_REL_GAP=1e-5, MSK_IPAR_NUM_THREADS=1, MSK_IPAR_INTPNT_MULTI_THREAD=0), mip_solver_drives=true, dualize_relax=true, dualize_subp=true, timeout=3600, rel_gap=1e-5,prim_cut_feas_tol=1e-6), strip(instancename))
        MathProgBase.optimize!(model)
        status = MathProgBase.status(model)
        time =  MathProgBase.getsolvetime(model)
        obj = MathProgBase.getobjval(model)
        nodes = MathProgBase.getnodecount(model)
        conicsolves = model.logs[:n_conic]
        status = MathProgBase.status(model)
        print("$instance     $status     $obj     $time     $status     $nodes     $conicsolves\n")
        open("$Filename.pajaritoresults", "a") do g
                write(g, "$instancename     $status     $obj     $time     $status     $nodes     $conicsolves \n")
                close(g)
        end
end
