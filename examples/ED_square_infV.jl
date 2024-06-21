pwd()
@static if false include("./src/WormQMC.jl") end
using Revise, WormQMC
using EDTools, LinearAlgebra, Statistics, SparseArrays, Dates, DelimitedFiles, Graphs
using Logging
# Logging.disable_logging(Logging.Info)

# Example 3×3 hard-core Hubbard
H = BH_Square(nmax=1, Lx=5, Ly=4, J=1.0, V=1e6, μ=0.2)
β = 24.0

begin
    function nextPXPconfig!(conf, n, G)
        nbs = neighbors(G, n)
        nbs = nbs[nbs.<n]
        N = length(conf)
        for i ∈ 1:N
            if ~any(conf[i][nbs])
                push!(conf, vcat(conf[i], true))
            end
            push!(conf[i], false)
        end
    end
    function PXPconfig(G)
        nmax = nv(G)
        conf = [[false] .== true, [true] .== true]
        for n = 2:nmax
            nextPXPconfig!(conf, n, G)
        end
        return conf
    end
end

function get_kets(H)
    L = (H.Lx, H.Ly)
    n_site = prod(L)
    Vac = zeros(Bool, L)
    lattice = CartesianIndices(Vac)
    Lids = LinearIndices(Vac)
    G = SimpleGraph(length(Vac))
    for r ∈ lattice
        rU = pbcshift(r, (0, 1), L)
        rR = pbcshift(r, (1, 0), L)
        add_edge!(G, Lids[r], Lids[rU])
        add_edge!(G, Lids[r], Lids[rR])
    end
    kets = sort(PXPconfig(G); by=sum)
    return kets
end

# L = Int.((H.Lx, H.Ly))
# ϕ = FBbasis(prod(L), 0, :hcboson, false)
function BH_Square_ED(H::BH_Square, β::f64)
    L = Int.((H.Lx, H.Ly))
    ϕ = FBbasis(get_kets(H), :hcboson, false)
    # ϕ = FBbasis(prod(L), 0, :hcboson, false)
    Lids = LinearIndices(L)
    ni = [densities(i; ϕ=ϕ) for i ∈ Lids];
    bi = [annihilation(i; ϕ=ϕ) for i ∈ Lids];
    b̂i = [creation(i; ϕ=ϕ) for i ∈ Lids];
    x_bonds = Pair.(circshift(LinearIndices(L), (0, 0)), circshift(LinearIndices(L), (-1, 0)))
    y_bonds = Pair.(circshift(LinearIndices(L), (0, 0)), circshift(LinearIndices(L), (0, -1)))
    bonds = vcat(x_bonds..., y_bonds...)
    bib̂jx = [hopping(i, j; ϕ=ϕ) for (i, j) ∈ x_bonds]
    bib̂jy = [hopping(i, j; ϕ=ϕ) for (i, j) ∈ y_bonds]
    Jx = let u = sum(bib̂jx)
        (u - u')
    end
    Jy = let u = sum(bib̂jy)
        (u - u')
    end
    T_term = let u = sum(bib̂jx) + sum(bib̂jy)
        u + u'
    end
    V_term = sum(ni[i] * ni[j] for (i, j) ∈ bonds) |> dropzeros
    @assert V_term.nzval |> isempty
    μ_term = sum(ni)
    ℋ = (- H.J * T_term - H.μ * μ_term) |> Matrix |> Hermitian
    G1j = [1*hopping(1, j; ϕ=ϕ) for j ∈ Lids]
    C1j = [1*densities(1, j; ϕ=ϕ) for j ∈ Lids]
    res = eigen(ℋ)
    Λ = res.values
    P = res.vectors
    @assert eltype(P) <: Real
    w_rel = @. exp.(-β * Λ)
    Z = sum(w_rel)
    normalize!(w_rel, 1)
    ρmat = Diagonal(w_rel)
    for i ∈ eachindex(bi, b̂i)
        matrix!(bi[i])
        matrix!(b̂i[i])
    end
    # function to_Ebasis(x::OpTerm)
    #     matrix!(x)
    #     return P' * (x.matrix * P)
    # end
    # to_Ebasis(x::AbstractMatrix) = P' * x * P
    # bi_ebasis = [b |> to_Ebasis for b ∈ bi]
    # b̂i_ebasis = [b̂ |> to_Ebasis for b̂ ∈ b̂i]
    # ni_op = [n |> to_Ebasis for n ∈ ni]
    # N_op = sum(ni) |> to_Ebasis
    # T_op = T_term |> to_Ebasis
    # μ_op = μ_term |> to_Ebasis
    Jx2 = Jx^2
    @assert ishermitian(Jx2)
    Jy2 = Jy^2
    @assert ishermitian(Jy2)

    # WxSquare = inv(abs2(H.Lx)).*Jx2 |> to_Ebasis |> hermitianpart
    # WySquare = inv(abs2(H.Ly)).*Jy2 |> to_Ebasis |> hermitianpart
    expval(o) = tr(ρmat * P' * (1*o) * P)
    res_dict = (
        N = μ_term |> expval,
        N2= (μ_term^2) |> expval,
        E = sum(Λ .* w_rel),
        K = expval(-H.J * T_term),
        U = 0.,
        μ = expval(- H.μ * μ_term),
        Wx2 = expval(Jx2 .* inv(abs2(H.Lx))),
        Wy2 = expval(Jy2 .* inv(abs2(H.Ly))),
        DensMat = [expval(g) for g ∈ G1j],
        DensCor = [expval(c) for c ∈ C1j],
    )
    return res_dict
end

update_const = UpdateConsts(1.0, 2.0, 1.0)
cycle_prob = CycleProb(1, 1, 1, 1)
# time_ther = Second(Minute(20))
# time_simu = Second(Minute(60))
time_ther = Second(20)
time_simu = Second(60)
x = Wsheet(β, H)
m = WormMeasure(x, update_const; green_lmax=1)
onesimu!(x, H, m, update_const, cycle_prob, time_ther, time_simu)

# using BenchmarkTools
# @btime worm_cycle!($x, $H, $update_const, $cycle_prob)
G0 = normalize_density_matrix(m.Gfunc)[:,:,1,1]
Cij = Cab(m.Sfact,1,1)
# using Plots
# heatmap(Cij)

# Benchmark with ED
res_ED = BH_Square_ED(H, β)

begin
@info "QMC result"
m.simple |> display
m.winding |> display
println("┌ DensMat")
writedlm(stdout, G0)
println("┌ DensCor")
writedlm(stdout, Cij)
end

begin
@info "ED result"
for k ∈ keys(res_ED)
    v = getfield(res_ED, k)
    println(k => getfield(res_ED, k))
end
end
println(
(mean(m.winding.Wy2) + mean(m.winding.Wx2)) / β
)
# G0
# mean(G0)
# 0.026864143447406823 / 0.02262594575520018

# Jx2_qmc = mean(m.winding.Wx2) * abs(H.Lx * H.Ly) * inv(2 * abs(β))
# Jy2_qmc = mean(m.winding.Wy2) * abs(H.Ly * H.Ly) * inv(2 * abs(β))
# res_ED.Wx2
# res_ED.Wy2
# (m.winding.Wy2 |> mean) / (2 * H.Ly^2 * β^2) * (H.Ly * H.Lx)

# res_dict["JxSquare"] / res_dict["JySquare"]
# res_dict["JxSquare"]
# #  / 4 * 3
# res_dict["JySquare"]
# #  / 4 * 3

# ΔΛ = [exp(Ei - Ej) for Ei ∈ Λ, Ej ∈ Λ]
# Gτ11(τ) = tr(ρmat * ((ΔΛ .^ τ) .* b̂i_ebasis[1]) * bi_ebasis[1])
# Gτ12(τ) = tr(ρmat * ((ΔΛ .^ τ) .* b̂i_ebasis[2]) * bi_ebasis[1])
# Ḡτ11(τ) = tr(ρmat * ((ΔΛ .^ τ) .* bi_ebasis[1]) * b̂i_ebasis[1])

# using Plots
# τgrid = LinRange(0.0, β, 100)
# Gτs = Gτ11.(τgrid)
# Ḡτs = Ḡτ11.(τgrid)
# Gτs[1]
# Gτs[end]
# Gτs[1] + Gτs[end]
# G12s = Gτ12.(τgrid)

# using LaTeXStrings
# scalefontsizes();
# scalefontsizes(1.5);
# plot(τgrid, Gτs, xlims=(0.0, β),
#     ylims=(0.0, 1.0),
#     xlabel=L"\tau", label=L"G_{ii}(\tau)", lw=2, xticks=([0.0, β / 2, β], [L"0", L"β/2", L"β"])
# )
# plot!(τgrid, reverse(Gτs), lw=2, label=L"G_{ii}(\beta - \tau)")
# plot!(τgrid, Ḡτs, label=L"\bar{G}_{ii}(\tau)")
# plot!(τgrid, G12s, label=L"G_{12}(\tau)")
# ylims!(0.1, 0.4)