# jobassign.jl
"""
job assignment related functions
"""

using StatsBase

# assigntable = nothing

# function setassigntable(atbl,itbl,gtbl)
#     global assigntable = atbl
#     global intervaltable = itbl
#     global grouptable = gtbl # basically not used
#     return
# end

function isassignable(x, lst)
    if x in lst
        return true
    else
        return false
    end
end

function findassignnum(x,atbl)
    k = []
    for i in 1:length(atbl)
        if isassignable(x,atbl[i])
            push!(k,i)
        end
    end
    k
end

# make correct assignment 
function assignjob!(jtbl::DVector{Int},atbl)
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

# penalty list for bad assignment
function listbadassign(jlst::Vector{Int}, alst::Vector{Int})
    st = vcat(alst,0)
    [(x in st) ? 0 : 1 for x in jlst]
end
function listbadassign(tbl::DVector{Int}, atbl::DVector{Int})
    [listbadassign(t,a) for (t,a) in zip(tbl,atbl)]
end

