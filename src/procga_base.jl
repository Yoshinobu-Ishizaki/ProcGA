# procga_base.jl
"""
Core routine for ProcGA
"""

using LinearAlgebra, Statistics
using Random
using DataFrames

# proctable has columns corresponding each process
# rows correspond to each material objects

# job table's columns are time units, therefore it is a schedule matrix
# rows correspond to each material

# create list of jobnumbers using proc table
function joblist(lst)
    a = []
    for i in 1:length(lst)
        append!(a,repeat([i],lst[i]))
    end
    return(a)
end

# make job time table base (ignoring sequence)
function jobtablebase(ptable)
    lp = sum(ptable)
    sz = size(ptable)
    
    tb = zeros(Int,sz[1],lp)
    for k in 1:sz[1]
        jl = joblist(ptable[k,:])
        tb[k,1:length(jl)] = jl
    end
    tb
end

# make shuffled job table
function jobtable(ptable)
    jtable=jobtablebase(ptable)
    for i in 1:size(jtable)[1]
        jtable[i,:] = shuffle(jtable[i,:])
    end
    jtable
end

# sort job table in order while keeping position of 0 entity
function possort!(jtbl)
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

# create initial population of size n
function initpopulation(n,ptable)
    jt = jobtable(ptable)
    popu = []
    for i in 1:n
        p = jobtable(ptable)
        possort!(p)
        push!(popu,p)
    end
    popu
end

# validation 

# basically each job number appearing in same column must be exclusive except groupable items until amount of count is below limit.
# sample: groupable = [Dict(:id =>[2,6], :cnt => 150)] # job 2,6 can be done simultaneously
# this table must be given manually

function isgroupable(x, gp)
    for i in 1:length(gp)
        if x in gp[i][:id]
            return(true)
        end
    end
    return(false)
end

# return capacity of amount of max quantity in a group
function groupcapacity(x,gp)
    for i in 1:length(gp)
        if x in gp[i][:id]
            return gp[i][:cnt]
        end
    end
    return(0)
end

# check if each column is valid combination of process
function checkvalidity(jtbl,gp,mt)
    flg = true
    k,ll = size(jtbl)
    for i in 1:ll
        q = mt[1,2]
        for j in 1:k
            a = jtbl[j,i]
            for m in (j+1):k
                if jtbl[m,i] == a
                    if isgroupable(a,gp)
                        q += mt[m,2]
                    else
                        flg = false
                    end
                end
            end
            if q > groupcapacity(a,gp)
                flg = false
            end
            if flg == false
                return(flg)
            end
        end
    end
    return(flg)
end

function validatejob1(jtbl,gp,mt)
    k,ll = size(jtbl)
    for i in 1:ll
        q = mt[1,2]
        for j in 1:k
            a = jtbl[j,i]
            for m in (j+1):k
                if jtbl[m,i] == a
                    if isgroupable(a,gp)
                        q += mt[m,2]
                        if q > groupcapacity(a,gp)
                            jtbl[m,i:end] = circshift(jtbl[m,i:end],1)
                            break
                        end
                    else
                        # circular shift latter data
                        jtbl[m,i:end] = circshift(jtbl[m,i:end],1)
                        break
                    end
                end
            end
        end
    end
    return(possort!(jtbl))
end

# 
function validatejob(jtbl,gp,mt,maxcount = 100)
    i = maxcount
    while !checkvalidity(jtbl,gp,mt) | i > 0
        jtbl = validatejob1(jtbl,gp,mt)
        i -= 1
    end
    if i == 0
        println("maxcount reached")
    end
    return(jtbl)
end

# get valid length of job time table
function validlength(jtb)
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

# shrink job table by skipping 0 column
function shrinkjob1!(jtbl)
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

function shrinkjob(jtbl)
    while true
        if !shrinkjob1!(jtbl)
            break
        end
    end
    jtbl
end

# sort population table by its valid length
function sortpopulation!(pptbl)
    sort!(pptbl, by = x->validlength(x))
end

# survival
# inferior genes cannot live long
function survive!(pptbl, survival = 0.8)
    ll = Int(floor(length(pptbl) * survival))
    splice!(pptbl,ll:length(pptbl))
    pptbl
end

# crossover and mutation
# crossover exchange rows of two genes at middle
function crossover(jtbl1, jtbl2)
    rows = size(jtbl1)[1]
    m = rows ÷ 2
    c1 = vcat(jtbl1[1:m,:],jtbl2[m+1:end,:])
    c2 = vcat(jtbl2[1:m,:],jtbl1[m+1:end,:])
    return(c1,c2)
end

function mutatejob!(jtbl,rate = 0.05)
    if rand() < rate
        m = ProcGA.validlength(jtbl)
        m1 = rand(1:m)
        m2 = rand(1:m)
        jtbl[:,m1], jtbl[:,m2] = jtbl[:,m2], jtbl[:,m1] 
    end
    return(jtbl)
end

# fill up population using elite children and mutant genes
function fillgeneration!(pptbl,n,elite = 0.2, mutant = 0.05)
    s = length(pptbl)
    ss = n - s
    
    while ss > 0
        en = Int(floor(s*elite))
        i1 = rand() < 0.2 ? rand(1:en) : rand(1:s)
        i2 = rand() < 0.2 ? rand(1:en) : rand(1:s)
        c1,c2 = crossover(pptbl[i1],pptbl[i2])
        mutatejob!(c1,mutant)
        mutatejob!(c2,mutant)
        push!(pptbl,c1,c2)
        ss -= 2
    end
    if length(pptbl) > n
        pop!(pptbl)
    end
    sortpopulation!(pptbl)
end

# evolution
function evolution!(pptbl, n, survival = 0.8, elite = 0.2, mutant = 0.05)
    rep = []
    s = length(pptbl)
    for i in 1:n
        survive!(pptbl,survival)
        fillgeneration!(pptbl,s,elite, mutant)
        v = validlength.(pptbl)
        dt = (minimum(v),median(v),maximum(v))
        #println(dt)
        push!(rep,dt)
    end
    rep
end
