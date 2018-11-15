# procga_base.jl
"""
Core routine for ProcGA
    gene must be givin as Array{Int,2} rowwise.
"""

using LinearAlgebra, Statistics
using Random
using DataFrames

# this must be overriden by user
function listpenalty()
    []
end

# taking column  
function coltake(lst::Array{Int,1},idx)
    lst[idx]
end

function coltake(tbl::Array{Array{Int,1},1},idx)
    coltake.(tbl,idx)
end

# return list of sum of each column
function colsumlist(tbl::Array{Array{Int,1},1})
    [sum(coltake(tbl,i)) for i in 1:length(tbl[1])]
end

# validation 
# get valid length of job time table
function validlength(lst::Array{Int,1})
    x = findlast(x->(x>0),lst)
    if x == nothing
        return 0
    else
        return x
    end
end

function validlength(tbl::Array{Array{Int,1},1})
    v = validlength.(tbl)
    maximum(v)
end

# calculate all penalty adding its validlength 
function penalty(tbl::Array{Array{Int,1},1})
    ptbl = listpenalty(tbl)
    sum(sum(ptbl)) + validlength(ptbl)
end

# ============================== penalty =================================

# consequtive use limit 
# if id continues more than lmt give penalty
function listcontinuity(lst::Array{Int,1},id,lmt)
    ln = length(lst)
    p = zeros(Int,ln)
    chk = [(x == id ? 1 : 0) for x in lst]
    
    cnt = 0
    x0 = 0
    for i in 1:ln
        x = chk[i]
        if x > 0
            cnt += 1
            p[i] = cnt
        else
            cnt = 0
        end
    end
    
    [ x > lmt ? (x-lmt) : 0 for x in p]
end

function listcontinuity(tbl::Array{Array{Int,1},1},id,lmt)
    [listcontinuity(lst,id,lmt) for lst in tbl] 
end   

function listcontinuity(tbl::Array{Array{Int,1},1}, lmt::Array{Int,1})
    sum([listcontinuity(tbl,i,lmt[i]) for i in 1:length(lmt)])
end

function listcoldup(tbl::Array{Array{Int,1},1})
    a = Array{Int,1}()
    for i in 1:length(tbl[1])
        v = coltake(tbl,i)
        v2 = v[v.>0]
        l = sum(v)>0 ? length(Set(v2)) : 0
        s = length(v2) - l
        push!(a,s)
    end
    [x .* a for x in tbl]
end


# ==================== job table editing =============================

# swap job list based on penaly list
function swapjob!(jlst::Array{Int,1},plst::Array{Int,1})
    pos = findmax(plst)[2]
    if pos > 0
        p2 = rand(1:length(jlst))
        jlst[pos],jlst[p2] = jlst[p2],jlst[pos]
    end
end
swapjob!(tbl::Array{Array{Int,1},1},ptbl::Array{Array{Int,1},1}) = swapjob!.(tbl,ptbl)

function swapjob!(tbl::Array{Array{Int,1},1})
    ptbl = listpenalty(tbl)
    swapjob!(tbl,ptbl)
end

# sort job table in order while keeping position of 0 entity
function orderjob!(jlst::Array{Int,1})
    st = sort(jlst[jlst.>0]) # sorted value
    pos = range(1,stop=length(jlst))[jlst.>0] # index of non zero

    for k in 1:length(pos)
        jlst[pos[k]] = st[k]
    end
    jlst
end
orderjob!(tbl::Array{Array{Int,1},1}) = orderjob!.(tbl)

# shrink job table by skipping 0 column
function shrinkjob!(tbl::Array{Array{Int,1},1})
    cols = validlength(tbl)
    for j in 1:cols
        v = coltake(tbl,j)
        if sum(v) == 0
            for i in length(tbl)
                tbl[i][j:end] = circshift(tbl[i][j:end],-1)
            end
        end
    end
end

# this has same effect as shrinkjob
function colsortjob!(tbl::Array{Array{Int,1},1})
    col = length(tbl[1])
    vs = [sum(coltake(tbl,i)) for i in 1:col]
    ky = sortperm(vs, by = x->(x>0 ? x : Inf))

    for i in 1:length(tbl)
        tbl[i] = tbl[i][ky]
    end
    tbl
    # [vs,ky]
end

function clipjob(jlst::Array{Int,1})
    jlst[1:validlength(jlst)]
end
clipjob!(tbl::Array{Array{Int,1},1}) = clipjob!.(tbl)


# ============================== crossover and mutation ==============================

# crossover exchange rows of two genes at middle
function crossover(jtbl1::Array{Array{Int,1},1}, jtbl2::Array{Array{Int,1},1})
    c1 = deepcopy(jtbl1)
    c2 = deepcopy(jtbl2)
    if c1 != c2
        rows = length(c1)
        m = rows ÷ 2
        c1 = vcat(c1[1:m],c2[m+1:end])
        c2 = vcat(c2[1:m],c1[m+1:end])
    end
    return(c1,c2)
end

# swap elements between valid index span and all
function mutatejob!(jlst::Array{Int,1})
    c1 = validlength(jlst)
    c2 = length(jlst)

    m1 = rand(1:c1)
    m2 = rand(1:c2) # swap valid cell with all including zero
    jlst[m1], jlst[m2] = jlst[m2], jlst[m1] 
    return(jlst)
end
mutatejob!(tbl::Array{Array{Int,1},1}) = mutatejob!.(tbl)

# ============================== population ==============================

function populateshuffle(tbl::Array{Array{Int,1},1},n::Int)
    popu = [deepcopy(tbl)]
    for i in 2:n
        jt = deepcopy(tbl)
        shuffle!.(jt)
        push!(popu,jt)
    end
    popu
end

# create initial population of size n from source data
function populatefrom(jlst::Array{Array{Int,1},1},n::Int)
    popu = [deepcopy(jlst)]
    for i in 2:n
        jt = deepcopy(jlst)
        mutatejob!.(jt)
        push!(popu, jt)
    end    
    popu
end

# sort population table by its valid length
function sortpopulation!(pptbl::Array{Array{Array{Int,1},1},1})
    sort!(pptbl, by = x->penalty(x))
end

# survival
# inferior genes cannot live long
function survive!(pptbl::Array{Array{Array{Int,1},1},1}, survival = 0.8)
    ll = Int(floor(length(pptbl) * survival))
    splice!(pptbl,ll:length(pptbl))
    pptbl
end

# fill up population using elite children and mutant genes
function fillgeneration!(pptbl::Array{Array{Array{Int,1},1},1},n::Int,elite = 0.2, mutant = 0.05, jobswitch = false)
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
            end
            # mutatejob!.(c1)
            # mutatejob!.(c2)
            swapjob!(c1)
            swapjob!(c2)
        end
        push!(pptbl,c1,c2)
        ss -= 2
    end
    while length(pptbl) > n
        pop!(pptbl)
    end
end

# evolution
function evolution!(pptbl::Array{Array{Array{Int,1},1},1}, n::Int, nr=10, survival=0.8, elite=0.2, mutant=0.05, jobswitch = false)
    rep = Array{Array{Int,1},1}()
    s = length(pptbl)
    for i in 1:n
        sortpopulation!(pptbl)
        survive!(pptbl,survival)
        fillgeneration!(pptbl,s,elite, mutant)
        shrinkjob!.(pptbl)
        # orderjob!.(pptbl)
        # colsortjob!.(pptbl)

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

