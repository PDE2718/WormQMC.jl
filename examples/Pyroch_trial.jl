using WormQMC
using LinearAlgebra, Statistics, Dates, DelimitedFiles, Logging, Accessors
# Logging.disable_logging(Logging.Info)

# Example 8×8 hard-core Hubbard (roughly in SF regime)
# U ≫ 1 and nmax ≫ 1 can also be set. Here U is useless.
H = BH_Pyroch(nmax=1, Lx=32, Ly=32, U=0.0, J1=1.0, J2=0.125, V=10.0, μ=2.0)
β = 4.0 # T = 1/β

# Update constants. Can be fine tuned.
update_const = UpdateConsts(0.5, 4.0, 1.0)
cycle_prob = CycleProb(1, 1, 1, 1)

# Thermalization and simulation time. All need to be in second.
time_ther = 5 |> Minute |> Second
time_simu = 5 |> Minute |> Second

# initialize the world line config and measurement
# green_lmax for imaginary-time green's function precision
x = Wsheet(β, H)
m = WormMeasure(x, update_const; green_lmax=1)

#! Do the simulation (should finish in time_ther+time_simu)
onesimu!(x, H, m, update_const, cycle_prob, time_ther, time_simu)

# calculate the density matrix with Gfunc buffer
G0 = normalize_density_matrix(m.Gfunc)[:, :, 1, 1]
# calculate the real-space correlation Cij = ⟨ninj⟩
Cij = Cab(m.Sfact, 1, 1)

# Now display the result
begin
    @info "QMC result"
    m.simple |> display
    m.winding |> display
    println("┌ DensMat")
    writedlm(stdout, G0)
    println("┌ DensCor")
    writedlm(stdout, Cij)
end

mean(m.simple.N) / (H.Lx * H.Ly)
m.winding


# Plot the imaginary green's function
using Plots,LaTeXStrings
τgrid = LinRange(0.0, β, 100)
G0τ = cal_Gτ(m.Gfunc, τgrid)[1]
let G00 = G0τ[1] # normalization with density matrix
    G0τ .*= (G0[1,1] ./ G00)
end
plot(τgrid,G0τ,
    xlims=(0.,β),ylims=(0,0.8),
    xlabel=L"τ",ylabel=L"G(τ,0)",
    framestyle=:box, label = false
)

# check that for hard-core bosons, G(0⁺) + G(0⁻) == 1
@assert ≈(G0τ[1] + G0τ[end], 1, atol=0.05)


using BenchmarkTools
@btime get_nbs($H, 1)

