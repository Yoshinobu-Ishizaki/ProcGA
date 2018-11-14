# procga_base.jl
"""
Core routine for ProcGA
    gene must be givin as Array{Int,2} rowwise.
"""

using LinearAlgebra, Statistics
using Random
using DataFrames

# sort job table in order while keeping position of 0 entity
function orderjob!(jtbl::Array{Int,2})
    for i in 1:size(jtbl)[1]
        lst = jtbl[i,:]
        st = sort(lst[lst.>0]) # sorted value
        pos = range(1,stop=length(lst))[lst.>0] # index of non zero
    
        for k in 1:length(pos)
            lst[pos[k]] = st[k]
        end
        jtbl[i,:] = lst
    end
    jtbl
end

# validation 

# get valid length of job time table
function validlength(jtb::Array{Int,2})
    l0 = size(jtb)[2]
    l = l0
    for i in l0:-1:1
        if sum(jtb[:,i]) == 0
            l-=1
        else
            break
        end
    end
    return(l)
end

# this must be overriden by user
function penalty()
    return(0)
end

# shrink job table by skipping 0 column
function shrinkjob1!(jtbl::Array{Int,2})
    flg = false
    cols = validlength(jtbl)
    for j in 1:cols
        if sum(jtbl[:,j]) == 0
            jtbl[:,j:end] = circshift(jtbl[:,j:end],(0,-1))
            flg = true
        end
    end
    return(flg)
end

function shrinkjob!(jtbl::Array{Int,2})
    while true
        if !shrinkjob1!(jtbl)
            break
        end
    end
    jtbl
end    

function clipjob(jtbl::Array{Int,2})
    m = validlength(jtbl)
    jtbl[:,1:m]
end

# crossover and mutation
# crossover exchange rows of two genes at middle
function crossover(jtbl1::Array{Int,2}, jtbl2::Array{Int,2})
    if jtbl1 != jtbl2
        rows = size(jtbl1)[1]
        m = rows ÷ 2
        c1 = vcat(jtbl1[1:m,:],jtbl2[m+1:end,:])
        c2 = vcat(jtbl2[1:m,:],jtbl1[m+1:end,:])
    else
        c1 = copy(jtbl1)
        c2 = copy(jtbl1)
    end
    return(c1,c2)
end

# swap columns of each row independently
function mutatejob!(jtbl::Array{Int,2})
    c1 = validlength(jtbl)
    rw,c2 = size(jtbl)
    for i in 1:rw
        m1 = rand(1:c1)
        m2 = rand(1:c2) # swap valid cell with all including zero
        jtbl[i,m1], jtbl[i,m2] = jtbl[i,m2], jtbl[i,m1] 
    end
    return(jtbl)
end

function shufflejob!(jtbl::Array{Int,2})
    for i in 1:size(jtbl)[1]
        jtbl[i,:] = shuffle(jtbl[i,:])
    end
    jtbl
end

# population 

# OBSOLETE!
# create initial population of size n
# function initpopulation(n)
#     # jt = jobtablebase()
#     popu = [] 
#     for i in 1:n
#         jt = jobtableshuffle() # init with shuffled data
#         # mutatejob!(jt,1) # make variation of base jobtable
#         shrinkjob!(jt)
#         push!(popu,jt)
#     end
#     popu
# end

function populateshuffle(jtbl::Array{Int,2},n::Int)
    popu = [copy(jtbl)]
    for i in 2:n
        jt = copy(jtbl)
        shufflejob!(jt)
        push!(popu,jt)
    end
    popu
end

    # create initial population of size n from source data
function populatefrom(jtbl::Array{Int,2},n::Int)
    popu = [copy(jtbl)]
    for i in 2:n
        jt = copy(jtbl)
        mutatejob!(jt)
        push!(popu, jt)
    end    
    popu
end

# sort population table by its valid length
function sortpopulation!(pptbl::Array{Array{Int,2},1})
    sort!(pptbl, by = x->penalty(x))
end

# survival
# inferior genes cannot live long
function survive!(pptbl::Array{Array{Int,2},1}, survival = 0.8)
    ll = Int(floor(length(pptbl) * survival))
    splice!(pptbl,ll:length(pptbl))
    pptbl
end

# fill up population using elite children and mutant genes
function fillgeneration!(pptbl::Array{Array{Int,2},1},n::Int,elite = 0.2, mutant = 0.05, jobswitch = false)
    s = length(pptbl)
    ss = n - s
    
    en = Int(floor(s*elite))
    while ss > 0
        if rand() < elite
            i1 = rand(1:en)
            i2 = rand(1:en)
        else
            i1 = rand(1:s)
            i2 = rand(1:s)
        end
        c1,c2 = crossover(pptbl[i1],pptbl[i2])
        if rand() < mutant
            if jobswitch
                # defined in jobassign.jl
                switchjob!(c1)
                switchjob!(c2)
            else
                mutatejob!(c1)
                mutatejob!(c2)
            end
        end
        push!(pptbl,c1,c2)
        ss -= 2
    end
    while length(pptbl) > n
        pop!(pptbl)
    end
end

# evolution
function evolution!(pptbl::Array{Array{Int,2},1}, n::Int, nr=10, survival=0.8, elite=0.2, mutant=0.05, jobswitch = false)
    rep = Array{Array{Int,1},1}()
    s = length(pptbl)
    for i in 1:n
        sortpopulation!(pptbl)
        survive!(pptbl,survival)
        fillgeneration!(pptbl,s,elite, mutant)
        shrinkjob!.(pptbl)
        orderjob!.(pptbl)

        v = penalty.(pptbl)
        dt = [minimum(v),Int(round(median(v))),maximum(v)]
        
        if i % nr == 0
            println("i:$i => $dt")
        end
        push!(rep,dt)
#        pptbl = clipjoblength(pptbl,maximum(v))
    end
    rep
end

