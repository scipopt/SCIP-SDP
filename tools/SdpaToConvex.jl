using Convex
using Pajarito
using Mosek

function SkipComments!(lines, pos)
        while lines[pos][1] == '*' || lines[pos][1] == '"'
                pos = pos + 1
        end
        return pos
end

function SdpaToConvex(Filename)
        f = open(Filename)
        lines = readlines(f)
        pos = 1

        pos = SkipComments!(lines, pos)
        nvars = parse(Int64, split(lines[pos])[1])
        pos = pos + 1
        #print("nvars: $nvars \n")

        pos = SkipComments!(lines, pos)
        nblocks = parse(Int64, split(lines[pos])[1])
        pos = pos + 1
        #print("nblocks: $nblocks \n")

        blocksizes = Array{Int64}(nblocks)
        pos = SkipComments!(lines, pos)
        maxblocksize = 0
        nsdpblocks = 0
        #entry i of sdpblockstart is the start index of the i-th sdp block in the large SDP matrix, it is equal to the sum of blocksizes of SDP-blocks 1 to i-1; the last entry gives the total size of the SDP matrix
        sdpblockstart = Array{Int64}(nblocks +1 )
        sdpblockstart[1] = 1
        for b in 1:nblocks
                blocksizes[b] = parse(Int64, split(lines[pos])[b])
                if abs(blocksizes[b]) > maxblocksize
                        maxblocksize = abs(blocksizes[b])
                end
                if blocksizes[b] > 0
                        nsdpblocks = nsdpblocks + 1
                        sdpblockstart[nsdpblocks + 1] = sdpblockstart[nsdpblocks] + blocksizes[b]
			#print("blockstart[$(b+1) ] = $(sdpblockstart[nsdpblocks + 1]) \n")
                end
        end
        #complete the totalsdpblocksizes-array
        if nsdpblocks < nblocks
                for b = nsdpblocks+1:nblocks
                        sdpblockstart[b + 1] = sdpblockstart[b]
                end
        end
        pos = pos + 1
        #print("blocksizes: $blocksizes \n")

        pos = SkipComments!(lines, pos)
        obj = zeros(nvars)
        for v in 1:nvars
                obj[v] = parse(Float64, split(lines[pos])[v])
        end
        pos = pos + 1
        #print("obj: $obj \n")

        #create A-array of size(nblocks,blocksize,blocksize,nvars) such that for block b: X[i,j] = sum A[b,i,j,k] * y[k] - const
        A = zeros(nblocks, maxblocksize, maxblocksize, nvars)
        cnst = zeros(nblocks, maxblocksize, maxblocksize)

        pos = SkipComments!(lines, pos)
        for line in lines[pos:end]
                # stop if we reached the integrality part
                if contains(line, "*INTEGER")
                        pos = pos + 1
                        break
                end
                # skip all comment lines
                if lines[pos][1] == '*' || lines[pos][1] == '"'
                        pos = pos + 1
                        continue
                end
                # check whether the entry belongs to the constant part
                if parse(Int64, split(line)[1]) > 0
                        # we add identical entries for lower and upper triangular part
                        A[parse(Int64, split(line)[2]), parse(Int64, split(line)[3]), parse(Int64, split(line)[4]), parse(Int64, split(line)[1])] = parse(Float64, split(line)[5])
                        A[parse(Int64, split(line)[2]), parse(Int64, split(line)[4]), parse(Int64, split(line)[3]), parse(Int64, split(line)[1])] = parse(Float64, split(line)[5])
                        #print("A[ $(parse(Int64, split(line)[2])), $(parse(Int64, split(line)[3])), $(parse(Int64, split(line)[4])), $(parse(Int64, split(line)[1])) ] = $(parse(Float64, split(line)[5])) \n")
                else
                        cnst[parse(Int64, split(line)[2]), parse(Int64, split(line)[3]), parse(Int64, split(line)[4])] = parse(Float64, split(line)[5])
                        cnst[parse(Int64, split(line)[2]), parse(Int64, split(line)[4]), parse(Int64, split(line)[3])] = parse(Float64, split(line)[5])
                        #print("cnst[ $(parse(Int64, split(line)[2])), $(parse(Int64, split(line)[3])), $(parse(Int64, split(line)[4])) ] = $(parse(Float64, split(line)[5])) \n")
                end
                pos = pos + 1
        end

        #create varmapper: varmapper[i] = j means y[i] =^= x[j], where x are the continuous variables, varmapper[i] = -j => y[i] =^= z[j], where z are the integer variables
        varmapper = Array{Int64}(nvars)

        #since we do not want to assume that the integrality constraints are ordered, we first generate a boolean array to decide which variables should be integer
        intvar = Array{Bool}(nvars)

        #initialize as false
        for v=1:nvars
                intvar[v] = false
        end

        #set values of intvar array
        for line in lines[pos:end]
                intvar[parse(Int64, line[2:end])] = true
                pos = pos + 1     
        end

        #compute varmapper
        nintvars = 0
        ncontvars = 0
        for v=1:nvars
                if intvar[v]
                        nintvars = nintvars + 1
                        varmapper[v] = -nintvars
                        #print("varmapper[$v] = $(-nintvars)\n")
                else
                        ncontvars = ncontvars + 1
                        varmapper[v] = ncontvars
                        #print("varmapper[$v] = $ncontvars\n")
                end
        end
        
        close(f)

        #start creating problem
        x = Variable(ncontvars)
        z = Variable(nintvars, :Int)
        X = Variable(sdpblockstart[nsdpblocks+1], sdpblockstart[nsdpblocks+1])
        
        # reorder objective vector according to varmapper
        ccont = zeros(ncontvars)
        cint = zeros(nintvars)
        
        for v=1:nvars
                pos = varmapper[v]
                if pos > 0
                        ccont[pos] = obj[v]
                else
                        cint[-pos] = obj[v]
                end
        end
        p = minimize(dot(ccont,x) + dot(cint,z))

	#add SDP-constraints
	for b=1:nsdpblocks
		p.constraints += (X[sdpblockstart[b]:sdpblockstart[b+1], sdpblockstart[b]:sdpblockstart[b+1]] in :SDP)
		#print("adding SDP-cone for indices $(sdpblockstart[b]) : $(sdpblockstart[b+1]) \n")
	end

        #iterate over all blocks
        for b = 1:nblocks
		#print("block $b \n")
                if blocksizes[b] > 0
                        #sdpblock
                        #for each matrix entry (i,j) create a vector of all entries of (A_k)_ij with correct indices according to varmapper
                        for i=1:blocksizes[b]
                               for j=1:blocksizes[b] 
                                        acont = zeros(ncontvars)
                                        aint = zeros(nintvars)
                                        for v=1:nvars
                                                pos = varmapper[v]
                                                if pos > 0
                                                        acont[pos] = A[b,i,j,v]
                                                else
                                                        aint[-pos] = A[b,i,j,v]
                                                end
					end
					#print("block $b : $acont , $aint \n")
					#print("X[$(sdpblockstart[b]+i-1), $(sdpblockstart[b]+j-1)]\n")
                                        p.constraints += dot(acont, x) + dot(aint,z) - cnst[b,i,j] == X[sdpblockstart[b]+i-1, sdpblockstart[b]+j-1] # -1 because both blockstart and i/j start from 1
                                end
                        end
                else
                        #lpblock
                        #for each diagonal entry create a vector of all lhs-coefficients with correct indices according to varmapper
                        for i=1:(-blocksizes[b])
                                acont = zeros(ncontvars)
                                aint = zeros(nintvars)
                                for v=1:nvars
                                        pos = varmapper[v]
                                        if pos > 0
                                                acont[pos] = A[b,i,i,v]
                                        else
                                                aint[-pos] = A[b,i,i,v]
                                        end
                                end
                                p.constraints += dot(acont, x) + dot(aint, z) >= cnst[b,i,i]
                        end
                end
        end
	solve!(p, PajaritoSolver(cont_solver=MosekSolver()))
end



SdpaToConvex("/home/gally/gally_SCIPSDP/instances/2x3_3bars.dat-s")
Convex.clearmemory()
SdpaToConvex("/home/gally/gally_SCIPSDP/instances/random_32_2_a.dat-s")
Convex.clearmemory()
SdpaToConvex("/home/gally/instances/cardLeastSquares/random_96_2_a.dat-s")
