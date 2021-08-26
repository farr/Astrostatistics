module Astrostatistics688

using MCMCChains

export append_generated_quantities

function _names(k, d)
    string(k)
end
function _names(k, d::Vector)
    ["$(k)[$(i)]" for i in 1:length(d)]
end

function _3Dify(d::Array{Float64, 2})
    ns, nc = size(d)
    reshape(d, (ns, 1, nc))
end
function _3Dify(d::Array{Array{Float64, 1}, 2})
    ns, nc = size(d)
    np, = size(d[1,1])

    out = zeros(ns, np, nc)
    for k in 1:nc
        for j in 1:np
            for i in 1:ns
                out[i,j,k] = d[i,k][j]
            end
        end
    end
    out
end

"""
    append_generated_quantities(trace, genq)

Given a trace (`MCMCChains` object) and a 2D array of named tuples corresponding
to `generated_quantities` output, concatenate the two together into a single
`MCMCChains` object.
"""
function append_generated_quantities(trace, genq)
    ks = keys(genq[1,1])
    nms = []
    ds = []
    for k in ks
        d = map(x -> getindex(x, k), genq)
        push!(nms, _names(k, d[1,1]))
        push!(ds, _3Dify(d))
    end
    nms = vcat(nms...)
    ds = cat(ds...; dims=2)

    nmap = trace.name_map
    atrace = cat(Array(trace, keys(nmap), append_chains=false)..., dims=3)
    Chains(cat(atrace, ds, dims=2), vcat(map(string, names(trace)), nms), nmap)
end

end
