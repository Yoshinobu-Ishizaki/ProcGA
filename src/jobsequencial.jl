# jobsequencial.jl
"""
sequencial job table handling
"""

using LinearAlgebra, Statistics
using Random
using DataFrames

# create list of jobnumbers using proc table
# example
# lst = [2,3,3] ==> [1,1,2,2,2,3,3,3]
function joblist(lst)
    a = Array{Int,1}()
    for i in 1:length(lst)
        append!(a,repeat([i],lst[i]))
    end
    a
end

# make job time table base (ignoring sequence)
function jobtablebase(ptable::Array{Int,2})
    sz = size(ptable)
    lp = sum(ptable)
    
    tb = DVector{Int}()
    for i in 1:sz[1]
        jl = joblist(ptable[i,:])
        lst = vcat(jl,zeros(Int,lp - length(jl)))
        push!(tb,lst)
    end
    tb
end


# return set of grouping number for given list
function getgroups(lst)
    gs = Int[]
    for gp in grouptable
        g = gp[:grp]
        for x in lst
            if x in gp[:id]
                push!(gs,g)
            end
        end
    end
    Set(gs)
end

function evolutecross!(pptbl::Array{DVector{Int},1}, survival=0.8, elite=0.2, mutant=0.05)
    s0 = length(pptbl)
    sortpopulation!(pptbl)
    survive!(pptbl,survival)
    s = length(pptbl)
    en = Int(round(elite*s))

    k = (s0-s)
    while k > 0
        if rand() < mutant
            i1 = rand(1:en)
            i2 = rand(1:en)
        else
            i1 = rand(1:s)
            i2 = rand(1:s)
        end

        c1,c2 = crossover(pptbl[i1],pptbl[i2])
        if rand() < mutant
            swapjobrow!(c1)
            swapjobrow!(c2)
        end
        push!(pptbl,c1)
        push!(pptbl,c2)
        k -= 2
    end
    v = penalty.(pptbl)
    [minimum(v),Int(round(median(v))),maximum(v)]
end

function seqevolution!(pptbl::Array{DVector{Int},1}, n::Int, nr=10, survival=0.8, elite=0.2, mutant=0.05; stopat = 0)
    rep = DVector{Int}()

    for i in 1:n
        # sortjob!.(pptbl)
        # sortjobcol!.(pptbl)
        dt = evolutecross!(pptbl,survival,elite,mutant)
        push!(rep,dt)

        if i % nr == 0
            println("i:$i => $dt")
        end
        if dt[1] <= stopat
            break
        end
    end
    rep
end

# give exclusive group duplication penalty 
# function grdupenalty(jtbl::Array{Int,2})
#     p = 0
#     for j in 1:validlength(jtbl)
#         lst = jtbl[:,j]
#         x = length(getgroups(lst))
#         p += ( x>1 ? (x-1) : 0 )
#     end
#     p
# end

# give if duplicate occurs at not groupable item
# function listduplicate(jtbl::Array{Array{Int,1},1}, glst::Array{Int,1})

#     p = zeros(Int,length(jtbl[1]))
#     for i in 1:validlength(jtbl)
#         col = coltake(jtbl,i)
#         cc = []
#         for c in col
#             if !(c in glist)
#                 push!(cc,c)
#             end
#         end

#         cs = Set(cc)
#         p[i] = (length(cc) - length(cs))
#     end
#     [p for x in jtbl]
# end

