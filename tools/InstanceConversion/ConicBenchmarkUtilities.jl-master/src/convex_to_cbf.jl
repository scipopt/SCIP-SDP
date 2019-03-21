

function convex_to_cbf(problem, name, filename)
    # work-around for optional dependency on Convex
    c, A, b, cones, var_to_ranges, vartypes, conic_constraints = eval(Main,:(Convex.conic_problem($problem)))
    var_cones = fill((:Free, 1:size(A, 2)),1)

    dat = mpbtocbf(name, vec(full(c)), A, vec(full(b)), cones, var_cones, vartypes, :Min)

    writecbfdata(filename, dat)
end
