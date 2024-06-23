# WormQMC.jl

A Julia implementation of the worm algorithm quantum Monte-Carlo (QMC). Largely inspired by [the original program by Nicolas Sadoune & Lode Pollet](https://github.com/LodePollet/worm).

## Overview
> Warning⚠️: The code is still under development.

`WormQMC.jl` can be installed with the Julia package manager. From the Julia REPL, type `]` to enter the Pkg REPL mode and run:
```bash
pkg> add https://github.com/PDE2718/WormQMC.jl
```

## Minimum example

Here is a minimum example: A Bose-Hubbard model on 2D square lattice with $L_x=Ly=β=8$ and PBC. Other parameters are $t=1.0, U=40.0, J=1.0, V=0.25, μ=0.0$. Note that in the program we replace the symbol $t$ with $J$ since the latter is more distinct against the label for imaginary time.

```julia
using WormQMC
using LinearAlgebra, Statistics, Dates, DelimitedFiles, Logging
# Logging.disable_logging(Logging.Info)

# Example 8×8 hard-core Hubbard (roughly in SF regime)
# U ≫ 1 and nmax ≫ 1 can also be set. Here U is useless.
H = BH_Square(nmax=1, Lx=8, Ly=8, U=40.0, J=1.0, V=0.25, μ=0.0)
β = 8.0 # T = 1/β

# Update constants. Can be fine tuned.
update_const = UpdateConsts(0.5, 1.0, 1.0)
cycle_prob = CycleProb(1, 1, 1, 1)

# Thermalization and simulation time. All need to be in second.
time_ther = 1 |> Minute |> Second
time_simu = 10 |> Minute |> Second

# initialize the world line config and measurement
# green_lmax for imaginary-time green's function precision
x = Wsheet(β, H)
m = WormMeasure(x, update_const; green_lmax=100)

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
@assert ≈(G0τ[1] + G0τ[end-1], 1, atol=0.03)
```




