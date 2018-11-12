# procga_base.jl
"""
Core routine for ProcGA
"""

using LinearAlgebra, Statistics
using Random
using DataFrames

# globals
# proctable has columns corresponding each process
# rows correspond to each material objects
# set common proc_table, group_table, mat_table
function settable(ptbl, gtbl, mtbl)
    global proctable = ptbl
    global grouptable = gtbl
    global mattable = mtbl
    return()
end

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
    sz = size(ptable)
    lp = sum(ptable)
    
    tb = zeros(Int,sz[1],lp)
    for i in 1:sz[1]
        #jl = findall(ptable[i,:].>0)
        jl = joblist(ptable[i,:])
        tb[i,1:length(jl)] = jl
    end
    tb
end
jobtablebase() = jobtablebase(proctable)

# make shuffled job table
function jobtableshuffle(ptable)
    jtable=jobtablebase(ptable)
    for i in 1:size(jtable)[1]
        jtable[i,:] = shuffle(jtable[i,:])
    end
    jtable
end
jobtableshuffle() = jobtableshuffle(proctable)

# sort job table in order while keeping position of 0 entity
function orderjob!(jtbl)
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
# basically each job number appearing in same column must be exclusive except groupable items until amount of count is below limit.
# sample: groupable = [Dict(:id =>[2,6], :cnt => 150)] # job 2,6 can be done simultaneously
# this table must be given manually

function isgroupable(x)
    gp = grouptable
    for i in 1:length(gp)
        if x in gp[i][:id]
            return(true)
        end
    end
    return(false)
end

# return capacity of amount of max quantity in a group
function groupcapacity(x)
    gp = grouptable
    for i in 1:length(gp)
        if x in gp[i][:id]
            return gp[i][:cnt]
        end
    end
    return(0)
end

# check if each column is valid combination of process
function checkvalidity(jtbl)
    mt = mattable
    flg = true
    k,ll = size(jtbl)
    for i in 1:ll
        if sum(jtbl[:,i]) > 0
            q = mt[1,2]
            for j in 1:k
                a = jtbl[j,i]
                if a > 0
                    for m in (j+1):k
                        if jtbl[m,i] == a
                            if isgroupable(a)
                                q += mt[m,2]
                            else
                                flg = false
                            end
                        end
                    end
                    if isgroupable(a) && q > groupcapacity(a)
                        flg = false
                    end
                end
                if flg == false
                    return(flg)
                end
            end
        end
    end
    return(flg)
end

function validatejob1!(jtbl)
    mt = mattable
    k,ll = size(jtbl)
    for i in 1:ll
        if sum(jtbl[:,i]) > 0
            q = mt[1,2]
            for j in 1:k
                a = jtbl[j,i]
                for m in (j+1):k
                    if jtbl[m,i] == a
                        if isgroupable(a)
                            q += mt[m,2]
                            if q > groupcapacity(a)
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
    end
end

# 
function validatejob!(jtbl,maxcount = 100)
    i = maxcount
    orderjob!(jtbl)
    while !checkvalidity(jtbl) && i > 0
        validatejob1!(jtbl)
        orderjob!(jtbl)
        i -= 1

    end
    if i == 0
        # println("maxcount reached")
        orderjob!(jtbl) # at least let it in order
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

# penalty to non-consecutive numbers!

# check incontinuity of lst corresponding to process list plst
function penalty0(lst,plst)
    p = 0
    e0 = 0 # element 
    cm = 0
    idx = findlast(x->(x>0),lst)
    
    for i in 1:idx
        e1 = lst[i]
        if e1 != e0
            if e1 > 0 
                e0 = e1
                cm = plst[e0] - 1 # max time units for e0-th process
            else
                # add penalty because 0 occured between numbers, cm can be zero
                p += cm
            end
        else 
            cm -= 1
        end
    end # for
    p 
end

# calculate group penalty 
function grpenalty(jtbl)
    p = 0
    for gp in grouptable
        cc = gp[:timespan]
        g1 = (x->(x in gp[:id]) ? x : 0).(jtbl)

        v0 = [0,0,0]
        for i in 1:validlength(g1)
            v1 = g1[:,i]
            if sum(v1) > 0
                if v1 != v0
                    v0 = v1
                    # add penalty 
                    p += cc-1
                else
                    p -= 1
                end

                # check count overflow
                m = (x->(x>0) ? 1 : 0).(v1)
                q = sum(m.*mattable[:,2])
                if q > gp[:cnt]
                    p += 10
                end
                # println("$i=>$p")
            else
                if v1 != v0
                    v0 = v1
                end
            end
        end # for i
    end # for k
    p
end

function penalty(jtbl,ptable = proctable)
    p = validlength(jtbl)
    for i in 1:size(jtbl)[1]
        p += penalty0(jtbl[i,:],ptable[i,:])
    end

    p += grpenalty(jtbl)
    p
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

function shrinkjob!(jtbl)
    while true
        if !shrinkjob1!(jtbl)
            break
        end
    end
    jtbl
end    

function clipjob(jtbl)
    m = validlength(jtbl)
    jtbl[:,1:m]
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

# swap columns of each row independently
function mutatejob!(jtbl,rate = 0.05)
    # m = validlength(jtbl)
    rw,m = size(jtbl)
    for i in 1:rw
        if rand() < rate
            m1 = rand(1:m)
            m2 = rand(1:m)
           jtbl[i,m1], jtbl[i,m2] = jtbl[i,m2], jtbl[i,m1] 
        end
    end
    return(jtbl)
end

# population 

# create initial population of size n
function initpopulation(n)
    # jt = jobtablebase()
    popu = [] 
    for i in 1:n
        jt = jobtableshuffle() # init with shuffled data
        # mutatejob!(jt,1) # make variation of base jobtable
        shrinkjob!(jt)
        # validatejob!(jt) # validate it after mutation
        push!(popu,jt)
    end
    popu
end

# create initial population of size n from source data
function initpopulationfrom(jtbl,n, mutant=1.0)
    popu = [copy(jtbl)]
    
    for i in 2:n
        jt = copy(jtbl)
        mutatejob!(jt,mutant)
        # validatejob!(jt)
        push!(popu, jt)
    end    
    popu
end

# sort population table by its valid length
function sortpopulation!(pptbl)
    sort!(pptbl, by = x->penalty(x))
end

# survival
# inferior genes cannot live long
function survive!(pptbl, survival = 0.8)
    ll = Int(floor(length(pptbl) * survival))
    splice!(pptbl,ll:length(pptbl))
    pptbl
end

# fill up population using elite children and mutant genes
function fillgeneration!(pptbl,n,elite = 0.2, mutant = 0.05)
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
        mutatejob!(c1,mutant)
        mutatejob!(c2,mutant)
        push!(pptbl,c1,c2)
        ss -= 2
    end
    if length(pptbl) > n
        pop!(pptbl)
    end
end

# evolution
# validatejob must be done before this
function evolution!(pptbl, n, nr=10, survival=0.8, elite=0.2, mutant=0.05)
    rep = []
    s = length(pptbl)
    for i in 1:n
        sortpopulation!(pptbl)
        survive!(pptbl,survival)
        fillgeneration!(pptbl,s,elite, mutant)
        shrinkjob!.(pptbl)
        orderjob!.(pptbl)
        # validatejob!.(pptbl)

        v = penalty.(pptbl)
        dt = (minimum(v),median(v),maximum(v))
        
        if i % nr == 0
            println("i:$i => $dt")
        end
        push!(rep,dt)
#        pptbl = clipjoblength(pptbl,maximum(v))
    end
    rep
end

