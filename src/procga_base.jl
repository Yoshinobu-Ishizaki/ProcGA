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
# calcutblate all penalty adding its validlength 
function penalty()
    return 0
end
# ============================================

function nzminimum(lst::Vector{Int})::Int
    l = lst[lst.>0]
    if length(l) > 0
        minimum(l)
    else
        zero(Int)
    end
end

# taking column  
function coltake(lst::Vector{Int},idx)
    lst[idx]
end

function coltake(tbl::DVector{Int},idx)
    [coltake(v,idx) for v in tbl]
end

# return list of sum of each column
# function colsumlist(tbl::DVector{Int})
#     # [sum(coltake(tbl,i)) for i in 1:length(tbl[1])]
#     sum(tbl) # this returns list 
# end

function zerotbl(tbl::DVector{Int})
    [ zeros(Int,length(tbl[1])) for x in tbl ]
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

function listusing(lst,id)
    [ x == id ? 1 : 0 for x in lst]
end

# ============================== penalty =================================

# job switching between assignment
function listjobswitch(tbl::DVector{Int},id::Int)
    utbl = listusing.(tbl,id)
    p = zerotbl(utbl)
    
    v1 = coltake(utbl,1)
    for i in 2:length(utbl[1])
        v = coltake(utbl,i)
        if (sum(v)>0) & (sum(v1)>0) & (v != v1)
            cc = findall(x->x>0, v)
            for x in cc
                p[x][i] = 1
            end
        end
        v1 = v
    end
    p
end
function listjobswitch(tbl::DVector{Int})
    # take unique non-zero elements
    xx = Int[]
    [append!(xx,unique(x)) for x in tbl]
    sort!(unique!(xx))
    el = xx[xx.>0]

    p = zerotbl(tbl)
    for e in el
        p .+= listjobswitch(tbl,e)
    end
    p
end

function listcoldup(tbl::DVector{Int})
    n = length(tbl)
    pp = zerotbl(tbl)

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
function listinhomogeneous(tbl::DVector{Int})
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

function listzero(lst::Vector{Int},n::Int)
    p = zeros(Int,length(lst))
    if n < 0 
        l = length(lst)
    else
        l = n
    end 
    for i in 1:l
        if lst[i] == 0
            p[i] = 1
        end
    end
    p
end
function listzero(lst::Vector{Int})
    l = validlength(lst)
    listzero(lst,l)
end
function listzero(tbl::DVector{Int}, n::Int)
    [ listzero(lst,n) for lst in tbl]
end
listzero(tbl::DVector{Int}) = listzero.(tbl)

# check overuse of same job overall
function listoveruse(tbl::DVector{Int}, id::Int, lmt::Int)
    vn = listusing.(tbl,id)
    v = sum(vn)
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
    v = sum(vn)
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
                cnt = lmt-1
            end
            used = false
        end
    end
    [pl .* x for x in vn]
end

# uses defautblt value
function listshortinterval(tbl::DVector{Int}, lmt::Vector{Int})
    pl = listshortinterval(tbl,1,lmt[1])
    for i in 2:length(lmt)
        pl .+= listshortinterval(tbl,i,lmt[i])
    end
    pl
end

# returns 1 if item is in grp
function listgroupable(lst::Vector{Int}, grp::Vector{Int})
    p = [ (x in grp) ? 1 : 0 for x in lst]   
end
function listgroupable(tbl::DVector{Int}, grp::Vector{Int})
    p = [ listgroupable(lst,grp) for lst in tbl ]
end
function listgroupable(tbl::DVector{Int}, gtbl::DVector{Int})
    p = zerotbl(tbl)
    for g in gtbl
        p .+= listgroupable(tbl,g)
    end
    p
end

# filter groupable items only
function listgroupableitems(tbl::DVector{Int}, grp::DVector{Int})
    gg = listgroupable(tbl,grp)
    [ t .* g for (t,g) in zip(tbl,gg )]
end

# filter out groupable
function listnongroupableitems(tbl::DVector{Int},grp::DVector{Int})
    gg = listgroupable(tbl,grp)
    gn = listzero(gg,-1)
    [ t .* g for (t,g) in zip(tbl,gn )]
end

# sequencial order
function listbadsequence(lst::Vector{Int})
    p = zeros(Int,length(lst))
    for i in 1:(length(lst)-1)
        m = nzminimum(lst[(i+1):end])
        if (lst[i] > 0) & (lst[i] > m)
            p[i] = m > 0 ? lst[i] - m : 0
        end
    end
    p
end
function listbadsequence(tbl::DVector{Int})
    listbadsequence.(tbl)
end

# penalty for bad grouping.
# if element is appear less than n (for each group) or more than n give penalty as its difference.
function listbadgroup(tbl::DVector{Int},grp::Vector{Int}, n::Int)
    nn = min(n,length(tbl))

    gg = listgroupable(tbl,grp)
    gl = sum(gg) # column sum
    pl = [ x == 0 ? 0 : abs(x-nn) for x in gl ]
    [ g .* pl for g in gg]
end
function listbadgroup(tbl::DVector{Int},grp::DVector{Int}, n::Vector{Int})
    p = zerotbl(tbl)
    for i in 1:length(n)
        p .+= listbadgroup(tbl,grp[i],n[i])
    end
    p
end

# if jobs are located sparce it is not good
function listlaziness(lst::Vector{Int})
    p = zeros(Int,length(lst))
    
    # simply give 1 to zero element between first and last non-zero elemant
    i1 = findfirst(x->x>0, lst)
    i2 = findlast(x->x>0, lst)
    for i in i1:i2
        p[i] = (lst[i] > 0) ? 0 : 1
    end
    
    # also give penalty to such that 7,0,7,7 
    xl0 = unique(lst)
    xl = xl0[xl0.>0] 
    for x in xl
        id = (1:length(lst))[lst.==x]
        sort!(id)
        i0 = id[1]
        for i in 2:length(id)
            i1 = id[i]
            p[i1] = (i1 - i0) > 1 ? (i1-i0-1) : 0
            i1 = i0
        end
    end
    p
end
function listlaziness(tbl::DVector{Int})
    listlaziness.(tbl)
end

# ==================== job table editing =============================

# enlarge job length and rotate last item
function addjobcol!(tbl::DVector{Int}, n::Int)
    t = [append!(x, zeros(Int,n)) for x in tbl]   
    t[end] = circshift(t[end],n)
    t
end

# randomly switch column,row of elements based on penalty list
function swapjob!(jtbl::DVector{Int})
    rw = length(jtbl)
    col = validlength(jtbl)
    plst = listpenalty(jtbl)
    
    c1,r1 = findmax((x->x[2]).(findmax.(plst))) # col, row of maximum

    zeropos = findall.(x->x==0, plst)
    r2 = rand(1:rw)
    c2 = rand(zeropos[r2])
    
    jtbl[r1][c1],jtbl[r2][c2] = jtbl[r2][c2], jtbl[r1][c1]
    jtbl
end

# swap job list based on penaly list
function swapjobrow!(jlst::Vector{Int},plst::Vector{Int})
    pos = findmax(plst)[2]
    if pos > 0
        pos2 = findall(x->x==0, plst)
        if length(pos2) > 0
            p2 = rand(pos2) # swap with zero penalty element
        else
            p2 = rand(1:length(jlst))
        end
        jlst[pos],jlst[p2] = jlst[p2],jlst[pos]
    end
    jlst
end
function swapjobrow!(tbl::DVector{Int})
    ptbl = listpenalty(tbl)
    swapjobrow!.(tbl,ptbl)
end

# sort job table in order while keeping position of 0 entity
function sortjob!(jlst::Vector{Int})
    st = sort(jlst[jlst.>0]) # sorted value
    pos = range(1,stop=length(jlst))[jlst.>0] # index of non zero

    for k in 1:length(pos)
        jlst[pos[k]] = st[k]
    end
    jlst
end
sortjob!(tbl::DVector{Int}) = sortjob!.(tbl)

# move zero column to end
function sortjobcol!(tbl::DVector{Int})
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
function clipjob(jlst::Vector{Int},n)
    jlst[1:n]
end
function clipjob(jlst::Vector{Int})
    clipjob(jlst,validlength(jlst))
end
clipjob(tbl::DVector{Int},n) = clipjob.(tbl, n)
clipjob(tbl::DVector{Int}) = clipjob.(tbl, validlength(tbl))
function clipjob!(tbl::DVector{Int},n)
    tbl = clipjob(tbl,n)
end

# ============================== crossover and mutation ==============================

# crossover exchange rows of two genes at middle
function crossover(jtbl1::DVector{Int}, jtbl2::DVector{Int})
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

# ============================== poputblation ==============================

function populateshuffle(tbl::DVector{Int},n::Int)
    popu = [deepcopy(tbl)]
    for i in 2:n
        jt = deepcopy(tbl)
        shuffle!.(jt)
        push!(popu,jt)
    end
    popu
end

# create initial poputblation of size n from source data
function populatefrom(jlst::DVector{Int},n::Int)
    jt = deepcopy(jlst)
    popu = [jt]
    for i in 2:n
        jt = deepcopy(jt)
        swapjobrow!(jt)
        push!(popu, jt)
    end
    popu
end

# sort poputblation table by its valid length
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
            swapjob!(c1)
        else
            swapjobrow!(c1) # change job for eachrow
        end
    end
    c1
end

# evolution

function evoluteone!(pptbl::Array{DVector{Int},1}, survival=0.8, elite=0.2, mutant=0.05, jobswitch = false; stopat = 0)
    s0 = length(pptbl)
    sortpopulation!(pptbl)
    survive!(pptbl,survival)
    s = length(pptbl)
    en = Int(round(elite*s))
    for i in 1:(s0-s)
        ic = rand(1:en)
        c = childjob(pptbl[ic],mutant,jobswitch)
        sortjobcol!(c)
        push!(pptbl,c)
    end
    v = penalty.(pptbl)
    [minimum(v),Int(round(median(v))),maximum(v)]
end

function evolution!(pptbl::Array{DVector{Int},1}, n::Int, nr=10, survival=0.8, elite=0.2, mutant=0.05, jobswitch = false; stopat = 0)
    rep = DVector{Int}()

    for i in 1:n
        dt = evoluteone!(pptbl,survival,elite,mutant, jobswitch)
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

function sortevolution!(pptbl::Array{DVector{Int},1}, n::Int, nr=10, survival=0.8, elite=0.2, mutant=0.05, jobswitch = false; stopat = 0)
    rep = DVector{Int}()

    for i in 1:n
        sortjob!.(pptbl)

        dt = evoluteone!(pptbl,survival,elite,mutant, jobswitch)
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