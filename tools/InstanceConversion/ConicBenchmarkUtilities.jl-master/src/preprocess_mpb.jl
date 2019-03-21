# Some MPB conic solvers don't support all equivalent representatons,
# although they should.

# A few utilities to normalize the MPB representation

# Remove :Zero cones from the variable cones
function remove_zero_varcones(c, A, b, con_cones, var_cones, vartypes)
    old_to_new_idx = zeros(Int,length(c))
    last_idx = 1
    new_varcones = Vector{Tuple{Symbol,Vector{Int}}}(0) 
    for (cname, cidx) in var_cones
        if cname != :Zero
            for i in cidx
                old_to_new_idx[i] = last_idx
                last_idx += 1
            end
            push!(new_varcones,(cname, old_to_new_idx[cidx]))
        else
            # dropped
        end
    end
    keep_indices = find(old_to_new_idx)
    c = c[keep_indices]
    A = A[:,keep_indices]
    vartypes = vartypes[keep_indices]
    var_cones = new_varcones

    return c, A, b, con_cones, var_cones, vartypes
end

# Introduce extra variables so that integer-constrained variables don't appear in nonlinear cones
function remove_ints_in_nonlinear_cones(c, A, b, con_cones, var_cones, vartypes)
    c = copy(c)
    b = copy(b)
    new_var_cones = Vector{Tuple{Symbol,Vector{Int}}}(0)
    con_cones = map((a) -> (a[1],vec(collect(a[2]))), con_cones)
    vartypes = copy(vartypes)

    nslack = 0
    dropped_idx = Int[]
    I = Int[]
    J = Int[]
    V = Float64[]
    for k in 1:length(var_cones)
        cname, cidx = var_cones[k]
        if !(cname ∈ (:Zero, :NonPos, :NonNeg, :Free))
            new_cidx = Int[]
            for i in cidx
                if vartypes[i] != :Cont
                    nslack += 1
                    push!(I,nslack)
                    push!(J,i)
                    push!(V,1.0)
                    push!(I,nslack)
                    push!(J,nslack+length(c))
                    push!(V,-1.0)
                    push!(new_cidx,nslack+length(c))
                    push!(dropped_idx, i)
                else
                    push!(new_cidx,i)
                end
            end
            push!(new_var_cones, (cname,new_cidx))
        else
            push!(new_var_cones, (cname, collect(cidx)))
        end
    end
    if length(dropped_idx) > 0
        push!(new_var_cones, (:Free, dropped_idx))
        append!(c, zeros(nslack))
        append!(b, zeros(nslack))
        append!(vartypes, fill(:Cont, nslack))
        push!(con_cones, (:Zero, collect((size(A,1)+1):(size(A,1)+nslack))))

        A = [A spzeros(size(A,1),nslack)]
        A = vcat(A, sparse(I,J,V,nslack,size(A,2)))
    end

    return c, A, b, con_cones, new_var_cones, vartypes

end

# Translate all SOCRotated cones to SOCs
function socrotated_to_soc(c, A, b, con_cones, var_cones, vartypes)

    con_cones = copy(con_cones)
    var_cones = copy(var_cones)
    vartypes = copy(vartypes)
    c = copy(c)
    b = copy(b)
    I, J, V = findnz(A)
    nslack = 0
    # introduce slack variables and put them into SOCRotated cones
    for i in 1:length(con_cones)
        cname, cidx = con_cones[i]
        if cname == :SOCRotated
            for j in cidx
                nslack += 1
                push!(I, j)
                push!(J, nslack+length(c))
                push!(V, 1.0)
                push!(vartypes,:Cont)
            end
            con_cones[i] = (:Zero, cidx)
            push!(var_cones, (:SOCRotated, (length(c)+nslack-length(cidx)+1):(length(c)+nslack)))
        end
    end
    append!(c, zeros(nslack))
    A = sparse(I,J,V,size(A,1),size(A,2)+nslack)

    # new rows to add to constraint matrix
    I = Int[]
    J = Int[]
    V = Float64[]
    rowidx = 1
    for i in 1:length(var_cones)
        cname, cidx = var_cones[i]
        if cname == :SOCRotated
            var_cones[i] = (:Free,cidx)
            # (y,z,x) in RSOC <=> (y+z,y-z,sqrt(2)*x) in SOC
            push!(I, rowidx)
            push!(J, cidx[1])
            push!(V, -1.0)
            push!(I, rowidx)
            push!(J, cidx[2])
            push!(V, -1.0)
            
            push!(I, rowidx+1)
            push!(J, cidx[1])
            push!(V, -1.0)
            push!(I, rowidx+1)
            push!(J, cidx[2])
            push!(V, 1.0)
            
            append!(I, (rowidx+2):(rowidx+length(cidx)-1))
            append!(J, cidx[3:end])
            append!(V, fill(-sqrt(2), length(cidx)-2))
            push!(con_cones, (:SOC, (size(A,1)+rowidx):(size(A,1)+rowidx+length(cidx)-1)))
            rowidx += length(cidx)
            append!(b, zeros(length(cidx)))
        end
    end

    A = vcat(A, sparse(I,J,V, rowidx-1, size(A,2)))

    return c, A, b, con_cones, var_cones, vartypes
end


