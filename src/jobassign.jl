# jobassign.jl
"""
job assignment related functions
"""

using StatsBase

# assigntable = nothing

function setassigntable(atbl,itbl,gtbl)
    global assigntable = atbl
    global intervaltable = itbl
    global grouptable = gtbl # basically not used
    return
end

function isassignable(x, lst)
    if x in lst
        return true
    else
        return false
    end
end

function findassignnum(x,atbl = assigntable)
    k = []
    for i in 1:length(atbl)
        if isassignable(x,atbl[i])
            push!(k,i)
        end
    end
    k
end

# make correct assignment 
function assignjob!(jtbl::Array{Array{Int,1},1},atbl = assigntable)
    for i in 1:length(jtbl)
        for j in 1:length(jtbl[i])
            x = jtbl[i][j]
            if !isassignable(x,atbl[i])
                rw = findassignnum(x,atbl)
                k2 = 0
                for r in rw
                    k2 = findfirst(x->x==0, jtbl[r])
                    if k2 > 0
                        jtbl[i][j], jtbl[r][k2] = jtbl[r][k2], jtbl[i][j]
                        break
                    end
                end
            end
        end # j
    end # i
    jtbl
end

# randomly switch column,row of elements if they are switchable.
function switchjob!(jtbl::Array{Array{Int,1},1}, atbl = assigntable)
    rw = length(jtbl)
    col = validlength(jtbl)
    
    r1 = rand(1:rw)
    c1 = rand(1:col)
    r2 = sample(findassignnum(r1,atbl),1)[1]
    c2 = sample(1:col,1)[1]
    jtbl[r1][c1],jtbl[r2][c2] = jtbl[r2][c2], jtbl[r1][c1]
    jtbl
end

