type CBFData
    name::String
    sense::Symbol
    var::Vector{Tuple{String,Int}}
    psdvar::Vector{Int}
    con::Vector{Tuple{String,Int}}
    psdcon::Vector{Int}
    objacoord::Vector{Tuple{Int,Float64}}
    objfcoord::Vector{Tuple{Int,Int,Int,Float64}}
    objoffset::Float64
    fcoord::Vector{Tuple{Int,Int,Int,Int,Float64}}
    acoord::Vector{Tuple{Int,Int,Float64}} # linear coefficients
    bcoord::Vector{Tuple{Int,Float64}} # linear offsets
    hcoord::Vector{Tuple{Int,Int,Int,Int,Float64}}
    dcoord::Vector{Tuple{Int,Int,Int,Float64}}
    intlist::Vector{Int}
    nvar::Int
    nconstr::Int
end

CBFData() = CBFData("",:xxx,[],[],[],[],[],[],0.0,[],[],[],[],[],[],0,0)

function parse_matblock(fd,outputmat,num_indices)
    nextline = readline(fd)
    nnz = parse(Int,strip(nextline))
    for k in 1:nnz
        nextline = readline(fd)
        tup = split(strip(nextline))
        push!(outputmat, (map(s->parse(Int,s)+1,tup[1:num_indices])...,parse(Float64,tup[end])))
    end
end


function readcbfdata(filename)

    if endswith(filename,"cbf.gz")
        fd = gzopen(filename,"r")
    else
        @assert endswith(filename, "cbf")
        fd = open(filename,"r")
    end

    dat = CBFData()
    dat.name = split(basename(filename),".")[1]

    while !eof(fd)
        line = readline(fd)
        startswith(line,"#") && continue # comments
        length(line) == 1 && continue # blank lines
        
        # new block
        
        if startswith(line,"VER")
            nextline = readline(fd)
            @assert startswith(nextline,"1") || startswith(nextline,"2") 
            continue
        end

        if startswith(line,"OBJSENSE")
            nextline = readline(fd)
            if strip(nextline) == "MIN"
                dat.sense = :Min
            else
                dat.sense = :Max
            end
            continue
        end

        if startswith(line,"VAR")
            nextline = readline(fd)
            totalvars, lines = split(nextline)
            totalvars = parse(Int,strip(totalvars))
            lines = parse(Int,strip(lines))
            varcnt = 0

            for k in 1:lines
                nextline = readline(fd)
                cone, sz = split(nextline)
                sz = parse(Int,strip(sz))
                push!(dat.var, (cone, sz))
                varcnt += sz
            end
            @assert totalvars == varcnt
            dat.nvar = varcnt
            continue
        end

        if startswith(line, "INT")
            nextline = readline(fd)
            intvar = parse(Int,strip(nextline))
            for k in 1:intvar
                nextline = readline(fd)
                idx = parse(Int,strip(nextline))
                push!(dat.intlist,idx+1)
            end
            continue
        end

        if startswith(line,"CON")
            nextline = readline(fd)
            totalconstr, lines = split(nextline)
            totalconstr = parse(Int,strip(totalconstr))
            lines = parse(Int,strip(lines))
            constrcnt = 0

            for k in 1:lines
                nextline = readline(fd)
                cone, sz = split(nextline)
                sz = parse(Int,strip(sz))
                push!(dat.con, (cone, sz))
                constrcnt += sz
            end
            @assert totalconstr == constrcnt
            dat.nconstr = constrcnt
            continue
        end

        if startswith(line,"PSDVAR")
            nextline = readline(fd)
            lines = parse(Int,strip(nextline))

            for k in 1:lines
                nextline = readline(fd)
                sz = parse(Int,strip(nextline))
                push!(dat.psdvar, sz)
            end
            continue
        end

        if startswith(line,"PSDCON")
            nextline = readline(fd)
            lines = parse(Int,strip(nextline))

            for k in 1:lines
                nextline = readline(fd)
                sz = parse(Int,strip(nextline))
                push!(dat.psdcon, sz)
            end
            continue
        end

        if startswith(line,"OBJACOORD")
            parse_matblock(fd,dat.objacoord,1)
        end

        if startswith(line,"OBJBCOORD")
            nextline = readline(fd)
            dat.objoffset = float(strip(nextline))
            warn("Instance has objective offset")
        end

        if startswith(line,"BCOORD")
            parse_matblock(fd,dat.bcoord,1)
        end

        if startswith(line,"ACOORD")
            parse_matblock(fd,dat.acoord,2)
        end

        if startswith(line,"OBJFCOORD")
            parse_matblock(fd,dat.objfcoord,3)
        end

        if startswith(line,"FCOORD")
            parse_matblock(fd,dat.fcoord,4)
        end

        if startswith(line,"HCOORD")
            parse_matblock(fd,dat.hcoord,4)
        end

        if startswith(line,"DCOORD")
            parse_matblock(fd,dat.dcoord,3)
        end


    end

    return dat
end


