

"""
    kruskal_mst(graph::BaseGraph,
                edges::Array{Int, 1},
                nodes::Array{Int, 1},
                weights::Array{Float64, 1})::BitSet

Generates and returns a minimum spanning tree from the subgraph induced
by `edges` and `nodes`, using Kruskal's MST algorithm. The `edges` are weighted
by `weights`.

## Note:
The `graph` represents the entire graph of the plan, where as `edges` and
`nodes` represent only the sub-graph on which we want to draw the MST.

*Arguments:*
- graph: Underlying Graph object
- edges: Array of edges of the sub-graph
- nodes: Set of nodes of the sub-graph
- weights: Array of weights of `length(edges)` where `weights[i]` is the
           weight of `edges[i]`

*Returns* a BitSet of edges that form a mst.
"""
function kruskal_mst(
    graph::BaseGraph,
    edges::Array{Int,1},
    nodes::Array{Int,1},
    weights::Array{Float64,1},
)::BitSet
    num_nodes = length(nodes)

    # sort the edges arr by their weights
    sorted_indices = sortperm(weights)
    sorted_edges = edges[sorted_indices]

    mst = BitSet()
    connected_vs = DisjointSets{Int}(nodes)

    for edge in sorted_edges
        if !in_same_set(connected_vs, graph.edge_src[edge], graph.edge_dst[edge])
            union!(connected_vs, graph.edge_src[edge], graph.edge_dst[edge])
            push!(mst, edge)
            (length(mst) >= num_nodes - 1) && break
        end
    end
    return mst
end

"""
    random_kruskal_mst(graph::BaseGraph,
                       edges::Array{Int, 1},
                       nodes::Array{Int, 1},
                       rng::AbstractRNG=Random.default_rng())

Generates and returns a random minimum spanning tree from the subgraph induced
by `edges` and `nodes`, using Kruskal's MST algorithm.

## Note:
The `graph` represents the entire graph of the plan, where as `edges` and
`nodes` represent only the sub-graph on which we want to draw the MST.

*Arguments:*
- graph: Underlying Graph object
- edges: Array of edges of the sub-graph
- nodes: Set of nodes of the sub-graph
- rng: A random number generator that implements the [AbstractRNG type](https://docs.julialang.org/en/v1/stdlib/Random/#Random.AbstractRNG) (e.g. `Random.default_rng()` or `MersenneTwister(1234)`)

*Returns* a BitSet of edges that form a mst.
"""
function random_kruskal_mst(
    graph::BaseGraph,
    edges::Array{Int,1},
    nodes::Array{Int,1},
    rng::AbstractRNG = Random.default_rng(),
)::BitSet
    weights = rand(rng, length(edges))
    return kruskal_mst(graph, edges, nodes, weights)
end

"""
    random_region_weighted_mst(graph::BaseGraph,
                       edges::Array{Int, 1},
                       nodes::Array{Int, 1},
                       rng::AbstractRNG=Random.default_rng())

Generates and returns a random minimum spanning tree from the subgraph induced
by `edges` and `nodes`, using Kruskal's MST algorithm. The weight associated
with each edge reflects the number of regional boundaries that are crossed
by that edge

## Note:
The `graph` represents the entire graph of the plan, where as `edges` and
`nodes` represent only the sub-graph on which we want to draw the MST.

*Arguments:*
- graph: Underlying Graph object
- edges: Array of edges of the sub-graph
- nodes: Set of nodes of the sub-graph
- rng: A random number generator that implements the [AbstractRNG type](https://docs.julialang.org/en/v1/stdlib/Random/#Random.AbstractRNG) (e.g. `Random.default_rng()` or `MersenneTwister(1234)`)
- region_weights: Array of 2 weights: 1st is the weight for counties, second is the weight for places

*Returns* a BitSet of edges that form a mst.
"""
function random_region_weighted_mst(
    graph::BaseGraph,
    edges::Array{Int,1},
    nodes::Array{Int,1},
    rng::AbstractRNG = Random.default_rng(),
    region_weights::Array{Tuple{String, Float64},1} = [("COUNTYFP10", 1.0), ("PLACE", 1.0)],
)::Tuple{BitSet, Vector{Float64}}
    weights = rand(rng,length(edges))
    for edge ∈ 1:length(edges)
        for (region, weight) ∈ region_weights
            if graph.attributes[graph.edge_src[edges[edge]]][region] != graph.attributes[graph.edge_dst[edges[edge]]][region]
                weights[edge] += weight
            end
        end
    end
    return (kruskal_mst(graph, edges, nodes, weights), weights)
end
