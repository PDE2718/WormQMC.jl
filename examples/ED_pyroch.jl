pwd()
@static if false include("./src/WormQMC.jl") end
using Revise, WormQMC
using EDTools, LinearAlgebra, Statistics, SparseArrays, Dates, DelimitedFiles
using Logging
# Logging.disable_logging(Logging.Info)

# Example 3×3 hard-core Hubbard
H = BH_Pyroch(nmax=1, Lx=2, Ly=3,
    J1=1.0, J2=1.0, V=0.25, μ=1.0
)
β = 4.0

function BH_ED(H::BH_Pyroch, β::f64)
    L = Int.((H.Lx, H.Ly, 2))
    ϕ = FBbasis(prod(L), 0, :hcboson, false)
    Lids = LinearIndices(L)
    ni = [densities(i; ϕ=ϕ) for i ∈ Lids]
    bi = [annihilation(i; ϕ=ϕ) for i ∈ Lids]
    b̂i = [creation(i; ϕ=ϕ) for i ∈ Lids]
    subA = collect(LinearIndices(L))[:, :, 1]
    subB = collect(LinearIndices(L))[:, :, 2]
    bonds_J1 = [
        Pair.(subA, circshift(subA, (+0, +1)))...,
        Pair.(subB, circshift(subB, (-1, +0)))...
    ]
    bonds_J2 = [
        Pair.(subA, circshift(subB, (+0, +0)))...
        Pair.(subA, circshift(subB, (+1, +0)))...
        Pair.(subB, circshift(subA, (-1, +1)))...
        Pair.(subB, circshift(subA, (+0, +1)))...
    ]
    bib̂j_J1 = [hopping(i, j; ϕ=ϕ) for (i, j) ∈ bonds_J1]
    bib̂j_J2 = [hopping(i, j; ϕ=ϕ) for (i, j) ∈ bonds_J2]
    T1 = sum(t for t ∈ bib̂j_J1)
    T1 = T1 + T1'
    T2 = sum(t for t ∈ bib̂j_J2)
    T2 = T2 + T2'
    @assert ishermitian(T1) && ishermitian(T2)
    V_term = sum(ni[i] * ni[j] for (i, j) ∈ vcat(bonds_J1, bonds_J2))
    μ_term = sum(ni)
    ℋ = (-H.J1 * T1 - H.J2 * T2 + H.V * V_term - H.μ * μ_term) |> Matrix |> Hermitian
    
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
    function to_Ebasis(x::OpTerm)
        matrix!(x)
        return P' * (x.matrix * P)
    end
    to_Ebasis(x::AbstractMatrix) = P' * x * P
    bi_ebasis = [b |> to_Ebasis for b ∈ bi]
    b̂i_ebasis = [b̂ |> to_Ebasis for b̂ ∈ b̂i]
    ni_op = [n |> to_Ebasis for n ∈ ni]
    N_op = sum(ni) |> to_Ebasis
    T1_op = T1 |> to_Ebasis
    T2_op = T2 |> to_Ebasis
    V_op = V_term |> to_Ebasis
    μ_op = μ_term |> to_Ebasis

    # Jx2 = Jx^2
    # @assert ishermitian(Jx2)
    # Jy2 = Jy^2
    # @assert ishermitian(Jy2)
    # WxSquare = inv(abs2(H.Lx)).*Jx2 |> to_Ebasis |> hermitianpart
    # WySquare = inv(abs2(H.Ly)).*Jy2 |> to_Ebasis |> hermitianpart
    
    D0 = zeros(H.Lx, H.Ly, 2, 2)
    for aa ∈ 1:2, bb ∈ 1:2
        for (i, j) ∈ Tuple.(CartesianIndices((H.Lx, H.Ly)))
            D0[i, j, aa, bb] = tr(
                ρmat * b̂i_ebasis[i, j, bb] * bi_ebasis[1, 1, aa]
            )
        end
    end
    C0 = zeros(H.Lx, H.Ly, 2, 2)
    for aa ∈ 1:2, bb ∈ 1:2
        for (i, j) ∈ Tuple.(CartesianIndices((H.Lx, H.Ly)))
            C0[i, j, aa, bb] = tr(
                ρmat * ni_op[i, j, bb] * ni_op[1, 1, aa]
            )
        end
    end

    res_dict = (
        N = tr(ρmat * N_op),
        E = sum(Λ .* w_rel),
        K = - tr(ρmat * (H.J1 * T1_op + H.J2 * T2_op)),
        U = 0.,
        μ = - H.μ * tr(ρmat * μ_op),
        V = + H.V * tr(ρmat * V_op),
        DensMat = D0,
        DensCor = C0,
    )
    return res_dict
end

update_const = UpdateConsts(0.5, 2.0, 1.0)
cycle_prob = CycleProb(1, 1, 1, 1)
time_ther = Second(20)
time_simu = Second(120)
x = Wsheet(β, H)
m = WormMeasure(x, update_const)
onesimu!(x, H, m, update_const, cycle_prob, time_ther, time_simu)
G0 = normalize_density_matrix(m.Gfunc)
# Benchmark with ED
res_ED = BH_ED(H, β)

begin
@info "QMC result"
m.simple |> display
m.winding |> display
println("┌ DensMat")
for aa ∈ 1:2, bb ∈ 1:2
    println("G_{$aa,$bb}")
    writedlm(stdout, G0[:,:,aa,bb])
end
println("┌ DensCor")
for aa ∈ 1:2, bb ∈ 1:2
    println("C_{$aa,$bb}")
    writedlm(stdout, Cab(m.Sfact, aa, bb))
end

@info "ED result"
for k ∈ keys(res_ED)
    v = getfield(res_ED, k)
    if v isa AbstractArray
        # println(k,":")
        # display(v)
    else
        println(k => v)
    end
end
for aa ∈ 1:2, bb ∈ 1:2
    println("G_{$aa,$bb}")
    writedlm(stdout, res_ED.DensMat[:, :, aa, bb])
end
println("┌ DensCor")
for aa ∈ 1:2, bb ∈ 1:2
    println("C_{$aa,$bb}")
    writedlm(stdout, res_ED.DensCor[:, :, aa, bb])
end

end
