# jobsequencial.jl
"""
sequencial job table handling
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
# function jobtableshuffle(ptable)
#     jtable=jobtablebase(ptable)
#     for i in 1:size(jtable)[1]
#         jtable[i,:] = shuffle(jtable[i,:])
#     end
#     jtable
# end
# jobtableshuffle() = jobtableshuffle(proctable)

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
function checkvalidity(jtbl::Array{Int,2})
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

# penalty of non-consecutive occurence
function serpenalty(jtbl::Array{Int,2})
    p = 0
    for k in 1:size(jtbl)[1]
        lst = jtbl[k,:]
        plst = proctable[k,:]

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
        end # i
    end # k
    p
end

# calculate group penalty 
function grpenalty(jtbl::Array{Int,2})
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

# give exclusive group duplication penalty 
function grdupenalty(jtbl::Array{Int,2})
    p = 0
    for j in 1:validlength(jtbl)
        lst = jtbl[:,j]
        x = length(getgroups(lst))
        p += ( x>1 ? (x-1) : 0 )
    end
    p
end

# give if duplicate occurs at not groupable item
function dupenalty(jtbl::Array{Int,2})
    p = 0
    for i in 1:validlength(jtbl)
        col = jtbl[:,i]
        col = col[col .> 0]
        cc = col[@. ~isgroupable(col)]
        cs = Set(cc)
        p += (length(cc) - length(cs))
    end
    p
end

# You can give your own penalty function by overriding this.
function seqpenalty(jtbl::Array{Int,2})
    p = validlength(jtbl)
    
    p += serpenalty(jtbl)
    p += grpenalty(jtbl)
    p += grdupenalty(jtbl)
    p += 5*dupenalty(jtbl) # weight

    p
end
