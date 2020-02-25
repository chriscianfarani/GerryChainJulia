# using GerryChain
using Test

include("../src/graph.jl")
include("../src/partition.jl")

const testdir = dirname(@__FILE__)
filepath = "./test_grid_4x4.json"

tests = [
    "graph",
    "partition"
]

@testset "GerryChainJulia" begin
    for t in tests
        tp = joinpath(testdir, "$(t).jl")
        include(tp)
    end
end
