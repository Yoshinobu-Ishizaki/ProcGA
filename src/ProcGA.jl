__precompile__()
"""
ProcGA

Process scheduling program with GA written in Julia.
"""
module ProcGA

# type definition for easy reading
DVector{T} = Vector{Vector{T}} where T

include("jobassign.jl") # assignment job table
include("jobsequencial.jl") # sequencial job table 
include("procga_base.jl") # basics for population

end