__precompile__()
"""
ProcGA

Process scheduling program with GA written in Julia.
"""
module ProcGA

include("jobassign.jl") # assignment job table
include("jobsequencial.jl") # sequencial job table 
include("procga_base.jl") # basics for population

end