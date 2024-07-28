using WormQMC
using LinearAlgebra, Statistics, Dates, DelimitedFiles, Logging, Accessors
# Logging.disable_logging(Logging.Info)

# Example 8×8 hard-core Hubbard (roughly in SF regime)
# U ≫ 1 and nmax ≫ 1 can also be set. Here U is useless.
H = BH_Trimer(nmax=1, Lx=30, Ly=30, U=0.0, J1=1.0, J2=1.0, V=10.0, μ=-0.8)
β = 20.0 # T = 1/β

# Update constants. Can be fine tuned.
update_const = UpdateConsts(0.5, 1.0, 1.0)
cycle_prob = CycleProb(1, 1, 1, 1)

# Thermalization and simulation time. All need to be in second.
time_ther = 60 |> Second
time_simu = 300 |> Second

# initialize the world line config and measurement
# green_lmax for imaginary-time green's function precision
x = Wsheet(β, H)
m = WormMeasure(x, H, update_const; green_lmax=1)
ψsnaps = Array{i8,3}[]

#! Do the simulation (should finish in time_ther+time_simu)
onesimu!(x, H, m, ψsnaps, update_const, cycle_prob, time_ther, time_simu)

ρmap = density_histogram(m.Sfact)
using Plots
heatmap(ρmap[2])
# calculate the density matrix with Gfunc buffer
G0 = normalize_density_matrix(m.Gfunc)[:, :, 1, 1]
# calculate the real-space correlation Cij = ⟨ninj⟩

ψ = [l[end].n_L for l ∈ x.wl]

function parse_bonds(ψ)
    ψA = ψ[:,:,1]
    ψB = ψ[:,:,2]
    xb = Union{Float64,Missing}[]
    yb = Union{Float64,Missing}[]
    wb = Union{Int64}[]
    cb = Union{Int64}[]
    Lx, Ly = size(ψA)
    for r ∈ CartesianIndices(ψA)
        x, y = Tuple(r)
        push!(xb, x - 1, x + 1, missing)
        push!(yb, y, y, missing)
        push!(wb, ψA[r], ψA[r], ψA[r])
        push!(cb, 1, 1, 1)

        push!(xb, x, x, missing)
        push!(yb, y - 1, y + 1, missing)
        push!(wb, ψB[r], ψB[r], ψB[r])
        push!(cb, 2, 2, 2)
    end
    pop!(wb)
    pop!(cb)
    return (xb, yb, wb, cb)
end

anidata = [parse_bonds(ψ) for ψ ∈ ψsnaps]

# ψsnaps = datas[100].ψsnaps
# anidata = [parse_bonds(ψ[:, :, 1], ψ[:, :, 2]) for ψ ∈ ψsnaps]
# anidata[1][end]

using Plots
function plotsegs(anidata, t)
    xb, yb, wb, cb = anidata[t]
    L = (H.Lx, H.Ly)
    plt = plot()
    vline!(plt, [1:L[1]], lw=1, lc=:grey, la=0.2, labels=false)
    hline!(plt, [1:L[2]], lw=1, lc=:grey, la=0.2, labels=false)
    plot!(plt, xb, yb, lw=2wb, lc=cb,
        # yflip=true,
        xlabel="x", ylabel="y", xlims=(0, L[1] + 1), ylims=(0, L[2] + 1), labels=false,
        framestyle=:box, aspect_ratio=:equal,
        # size=(500, 500),
        grid=false
    )
    # scatter!(holes, mc=4, marker=(2, stroke(0)), labels=false)
    annotate!((0, 1.05), text("t=$(t)", :top))
end
plotsegs(anidata, 1)

ani = @animate for i ∈ eachindex(anidata)[1:5:end]
    plotsegs(anidata, i)
end
gif_ani = gif(ani, fps=12)
gif_ani |> display

heatmap(0:17, 0:17, Cij',
    yflip=true,
    ylabel="x", xlabel="y", xlims=(0, 18), ylims=(0, 18),
    labels=false,
    framestyle=:box, aspect_ratio=:equal, size=(500, 500), grid=false
)
Cij = Cab(m.Sfact, 1, 1)
heatmap(Cij)

G0 = normalize_density_matrix(m.Gfunc)[:, :, 2, 2]
heatmap(0:17, 0:17, G0',
    yflip=true,
    ylabel="x", xlabel="y", xlims=(0, 18), ylims=(0, 18),
    labels=false,
    framestyle=:box, aspect_ratio=:equal, size=(500, 500), grid=false
)
plot(G0[1:end,1])
length(-1.14:0.04:-0.9)
# using Plots
# heatmap(G0)


# # Now display the result
# begin
#     @info "QMC result"
#     m.simple |> display
#     m.winding |> display
#     println("┌ DensMat")
#     writedlm(stdout, G0)
#     println("┌ DensCor")
#     writedlm(stdout, Cij)
# end

# mean(m.simple.N) / (H.Lx * H.Ly)
# m.winding


# # Plot the imaginary green's function
# using Plots,LaTeXStrings
# τgrid = LinRange(0.0, β, 100)
# G0τ = cal_Gτ(m.Gfunc, τgrid)[1]
# let G00 = G0τ[1] # normalization with density matrix
#     G0τ .*= (G0[1,1] ./ G00)
# end
# plot(τgrid,G0τ,
#     xlims=(0.,β),ylims=(0,0.8),
#     xlabel=L"τ",ylabel=L"G(τ,0)",
#     framestyle=:box, label = false
# )

# # check that for hard-core bosons, G(0⁺) + G(0⁻) == 1
# @assert ≈(G0τ[1] + G0τ[end], 1, atol=0.05)


# using BenchmarkTools
# @btime get_nbs($H, 1)
