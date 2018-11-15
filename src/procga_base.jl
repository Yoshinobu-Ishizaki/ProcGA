# procga_base.jl
"""
Core routine for ProcGA
    gene must be givin as Array{Int,2} rowwise.
"""

using LinearAlgebra, Statistics
using Random
using DataFrames

# abstract function ==========================
# these must be overriden by user
function listpenalty()
    []
end
# calculate all penalty adding its validlength 
function penalty()
    return 0
end
# ============================================

# taking column  
function coltake(lst::Vector{Int},idx)
    lst[idx]
end

function coltake(tbl::DVector{Int},idx)
    [coltake(v,idx) for v in tbl]
end

# return list of sum of each column
function colsumlist(tbl::DVector{Int})
    # [sum(coltake(tbl,i)) for i in 1:length(tbl[1])]
    sum(tbl) # this returns list 
end

# validation 
# get valid length of job time table
function validlength(lst::Vector{Int})
    x = findlast(x->(x>0),lst)
    if x == nothing
        return 0
    else
        return x
    end
end

function validlength(tbl::DVector{Int})
    v = validlength.(tbl)
    maximum(v)
end


# ============================== penalty =================================

# consequtive use limit 
# if id continues more than lmt, give penalty
# function listcontinuity(lst::Vector{Int},id,lmt)
#     ln = length(lst)
#     p = zeros(Int,ln)
#     chk = [(x == id ? 1 : 0) for x in lst]
    
#     cnt = 0
#     x0 = 0
#     for i in 1:ln
#         x = chk[i]
#         if x > 0
#             cnt += 1
#             p[i] = cnt
#         else
#             cnt = 0
#         end
#     end
    
#     [ x > lmt ? (x-lmt) : 0 for x in p]
# end

# function listcontinuity(tbl::DVector{Int},id,lmt)
#     [listcontinuity(lst,id,lmt) for lst in tbl] 
# end   

# function listcontinuity(tbl::DVector{Int}, lmt::Vector{Int})
#     sum([listcontinuity(tbl,i,lmt[i]) for i in 1:length(lmt)])
# end

function listcoldup(tbl::DVector{Int})
    n = length(tbl)
    pp = [ zeros(Int,length(tbl[1])) for x in tbl]

    for i in 1:length(tbl[1])
        p = zeros(Int,length(tbl))
        v = coltake(tbl,i)
        v2 = v[v.>0]
        for x in unique(v2)
            ck = v[v.==x]
            if length(ck) > 1 # duplicate
                p .+= [ c == x ? 1 : 0 for c in v ]
            end
        end
        for j in 1:n
            pp[j][i] = p[j]
        end
    end
    pp
end

# add penalty if each column is not same as preceeding one
function listhomogenious(tbl::DVector{Int})
    v0 = coltake(tbl,1)
    p = zeros(Int,length(tbl[1]))
    for i in 2:length(tbl[1])
        v1 = coltake(tbl,i)
        if (v1 != v0) & (Set(v1) == Set(v0))
            p[i] = 1
        end
        v0 = v1
    end
    [p for x in tbl] # return as table form
end

function listzero(lst::Vector{Int})
    l = validlength(lst)
    p = zeros(Int,length(lst))
    for i in 1:l
        if lst[i] == 0
            p[i] = 1
        end
    end
    p
end
listzero(tbl::DVector{Int}) = listzero.(tbl)

function listusing(lst,id)
    [ x == id ? 1 : 0 for x in lst]
end

# check overuse of same job overall
function listoveruse(tbl::DVector{Int}, id::Int, lmt::Int)
    vn = listusing.(tbl,id)
    v = colsumlist(vn)
    pl = zeros(Int,length(tbl[1]))

    used = false
    cnt = 0
    for i in 1:validlength(v)
        if v[i] > 0
            if used 
                cnt += 1
            else
                used = true
                cnt = 1
            end
            pl[i] = cnt > lmt ? (cnt-lmt) : 0
        else
            if used
                used = false
                cnt = 0
            end
        end
    end
    [pl .* x for x in vn]
end

function listoveruse(tbl::DVector{Int}, lmt::Vector{Int})
    pl = listoveruse(tbl,1,lmt[1])
    for i in 2:length(lmt)
        pl .+= listoveruse(tbl,i,lmt[i])
    end
    pl
end

# check short interval between same work 
function listshortinterval(tbl::DVector{Int}, id::Int, lmt::Int)
    vn = listusing.(tbl,id)
    v = colsumlist(vn)
    pl = zeros(Int,length(tbl[1]))

    used = false
    cnt = 0
    for i in 1:validlength(v)
        if v[i] > 0
            if !used 
                pl[i] = cnt > 0 ? cnt : 0
            end
            used = true
        else
            if !used
                cnt -=  1
            else
                cnt = lmt
            end
            used = false
        end
    end
    [pl for x in tbl]
end

# uses default value
function listshortinterval(tbl::DVector{Int}, lmt::Vector{Int})
    pl = listshortinterval(tbl,1,lmt[1])
    for i in 2:length(lmt)
        pl .+= listshortinterval(tbl,i,lmt[i])
    end
    pl
end

# ==================== job table editing =============================

# randomly switch column,row of elements based on penalty list
function switchjob!(jtbl::DVector{Int},plst::DVector{Int})
    rw = length(jtbl)
    col = validlength(jtbl)
    
    c1 = maximum((x->x[2]).(findmax.(plst)))
    r1 = rand(1:rw)

    c2 = rand(1:col)
    r2 = rand(1:rw)
    
    jtbl[r1][c1],jtbl[r2][c2] = jtbl[r2][c2], jtbl[r1][c1]
    jtbl
end

# swap job list based on penaly list
function swapjob!(jlst::Vector{Int},plst::Vector{Int})
    pos = findmax(plst)[2]
    if pos > 0
        p2 = rand(findall(x->x==0, plst)) # swap with zero penalty element
        jlst[pos],jlst[p2] = jlst[p2],jlst[pos]
    end
end
swapjob!(tbl::DVector{Int},ptbl::DVector{Int}) = swapjob!.(tbl,ptbl)

function swapjob!(tbl::DVector{Int})
    ptbl = listpenalty(tbl)
    swapjob!(tbl,ptbl)
end

# sort job table in order while keeping position of 0 entity
function orderjob!(jlst::Vector{Int})
    st = sort(jlst[jlst.>0]) # sorted value
    pos = range(1,stop=length(jlst))[jlst.>0] # index of non zero

    for k in 1:length(pos)
        jlst[pos[k]] = st[k]
    end
    jlst
end
orderjob!(tbl::DVector{Int}) = orderjob!.(tbl)

# shift job table by skipping 0 column
function shiftjob!(tbl::DVector{Int})
    cols = validlength(tbl)
    v = sum(tbl)
    for j in 1:cols
        if sum(v[j]) == 0
            for i in length(tbl)
                tbl[i][j:end] = circshift(tbl[i][j:end],-1)
            end
        end
    end
end

# this has same effect as shiftjob
function colsortjob!(tbl::DVector{Int})
    col = length(tbl[1])
    vs = sum(tbl)
    ky = sortperm(vs, by = x->(x>0 ? x : Inf))

    for i in 1:length(tbl)
        tbl[i] = tbl[i][ky]
    end
    tbl
    # [vs,ky]
end

# condense all by skipping 0 element
function condensejob!(lst::Vector{Int})
    v = lst[lst.>0]
    vm = length(v)

    for i in 1:length(lst)
        if i > vm
            lst[i] = 0
        else
            lst[i] = v[i]
        end
    end
    lst
end
condensejob!(tbl::DVector{Int}) = condensejob!.(tbl)

# remove trailing zeros
function clipjob(jlst::Vector{Int})
    jlst[1:validlength(jlst)]
end
clipjob!(tbl::DVector{Int}) = clipjob!.(tbl)

# ============================== crossover and mutation ==============================

# crossover exchange rows of two genes at middle
# function crossover(jtbl1::DVector{Int}, jtbl2::DVector{Int})
#     c1 = deepcopy(jtbl1)
#     c2 = deepcopy(jtbl2)
#     if c1 != c2
#         rows = length(c1)
#         m = rows ÷ 2
#         c1 = vcat(c1[1:m],c2[m+1:end])
#         c2 = vcat(c2[1:m],c1[m+1:end])
#     end
#     return(c1,c2)
# end

# swap elements between valid index span and all
# obsolete: use swapjob
# function mutatejob!(jlst::Vector{Int})
#     c1 = validlength(jlst)
#     c2 = length(jlst)

#     m1 = rand(1:c1)
#     m2 = rand(1:c2) # swap valid cell with all including zero
#     jlst[m1], jlst[m2] = jlst[m2], jlst[m1] 
#     return(jlst)
# end
# mutatejob!(tbl::DVector{Int}) = mutatejob!.(tbl)

# ============================== population ==============================

function populateshuffle(tbl::DVector{Int},n::Int)
    popu = [deepcopy(tbl)]
    for i in 2:n
        jt = deepcopy(tbl)
        shuffle!.(jt)
        push!(popu,jt)
    end
    popu
end

# create initial population of size n from source data
function populatefrom(jlst::DVector{Int},n::Int)
    popu = [deepcopy(jlst)]
    for i in 2:n
        jt = deepcopy(jlst)
        mutatejob!.(jt)
        push!(popu, jt)
    end    
    popu
end

# sort population table by its valid length
function sortpopulation!(pptbl::Array{DVector{Int},1})
    sort!(pptbl, by = x->penalty(x))
end

# survival
# inferior genes cannot live long
function survive!(pptbl::Array{DVector{Int},1}, survival = 0.8)
    ll = Int(floor(length(pptbl) * survival))
    splice!(pptbl,ll:length(pptbl))
    pptbl
end

# new generation from parent (uni-sex generation)
function childjob(tbl::DVector{Int}, mutant = 0.05, jobswitch = false)
    # make a copy of selected gene
    c1 = deepcopy(tbl)
    if rand() < mutant
        if jobswitch
            switchjob!(c1)
        else
            swapjob!(c1) # change job for eachrow
        end
    end
    c1
end

# evolution
function evolution!(pptbl::Array{DVector{Int},1}, n::Int, nr=10, survival=0.8, elite=0.2, mutant=0.05, jobswitch = false)
    rep = DVector{Int}()
    s0 = length(pptbl)

    for i in 1:n
        sortpopulation!(pptbl)
        survive!(pptbl,survival)
        
        s = length(pptbl)
        en = Int(round(elite*s))
        for i in 1:(s0-s)
            ic = rand(1:en)
            c = childjob(pptbl[ic],elite, mutant)
            shiftjob!(c)
            push!(pptbl,c)
        end

        v = penalty.(pptbl)
        dt = [minimum(v),Int(round(median(v))),maximum(v)]
        
        if i % nr == 0
            println("i:$i => $dt")
        end
        push!(rep,dt)
    end
    rep
end

