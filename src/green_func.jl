# Ref: https://ar5iv.labs.arxiv.org/html/1104.3215
# The idea : Legendre Polynomial Expansion
# Gτ = β⁻¹ ∑ₗ P̄ₗ(x) Gₗ
# Gₗ = ∫₀ᵝ dt P̄ₗ(x) Gτ
# P̄(x) = √(2l+1) P(x) : 
#   Legendre Polynomials with schmidt normalization, P(x) is the standard
# here: x(τ) = 2τ/β - 1   ∈ [-1,1]

# accumulator of Matsubara Green Functions & density matrix.
# if not is_full, then measure only the on site one
# only for systems with translational symmetry
# 
mutable struct GreenFuncBin
    const is_full::Bool
    const β::f64
    const Cw::f64
    const _Pl::Vector{f64} # Pₗ buffer for speed
    const Gl::AbstractArray{Vector{f64}}
    const G0::AbstractArray{Complex{Int}}
    insertion_trial::Int
end
import Base:empty!
function empty!(G::GreenFuncBin)::Nothing
    for Gl ∈ G.Gl
        fill!(Gl, 0.)
    end
    fill!(G.G0, complex(0))
    G.insertion_trial = 0
    return nothing
end
function GreenFuncBin(x::Wsheet, lmax::Int, Cw::f64; is_full=false)
    @assert 2 ≤ ndims(x.wl) ≤ 3
    _Pl = zeros(lmax+1)
    Nsub::Int = ndims(x.wl) == 3 ? size(x.wl)[end] : 1
    G_lattice = CartesianIndices((size(x.wl, 1), size(x.wl, 2), Nsub, Nsub))
    Gl = [zeros(lmax + 1) for _ ∈ (is_full ? G_lattice : CartesianIndices((1,1,Nsub,Nsub)))]
    G0 = [complex(0) for _ ∈ G_lattice]
    return GreenFuncBin(is_full, x.β, Cw, _Pl, Gl, G0, 0)
end
function Base.show(io::IO, m::GreenFuncBin)
    println(io, "┌ GreenFuncBin:")
    println(io, "│ is_full = $(m.is_full), β = $(m.β)")
    println(io, "│ Cw = $(m.Cw), lmax = $(length(m._Pl)-1)")
    println(io, "│ G0 → $(summary(m.G0))")
    println(io, "└ Gl → $(summary(m.Gl))")
end
# HH = BH_Square(Lx=3,Ly=3)
# xx = Wsheet(100., HH)
# GG = GreenFuncBin(xx, 10; is_full = false)
# GG.Gl
# GG.G0

@inline function cal_Δτ(τ::f64, τ̂::f64, β::f64)::f64
    return mod(τ̂ - τ, β)
end
@inline function xτmap(Δτ::f64, β::f64)::f64
    return (2Δτ/β) - 1.
end
Base.@propagate_inbounds function accum_green!(G::GreenFuncBin, w::Worm, H::BH_Parameters)::Nothing
    b::Element = b̂::Element = w.head::Element
    if w.head.op == b_
        b̂ = w.tail
    else
        b = w.tail
    end
    Gid = site_diff(H, b.i, b̂.i)
    #measure density matrix
    @inbounds begin
    if w.loc == _at_green
        G.G0[Gid] += complex(1,0)
    elseif w.loc == _at_stop
        G.G0[Gid] += (b̂.t > b.t) ? complex(1, 0) : complex(0, 1)
    elseif w.loc == _at_free && (G.is_full || Gid[1]==Gid[2]==1)
        Δτ = cal_Δτ(b.t, b̂.t, G.β)
        x = xτmap(Δτ, G.β)
        collectPl!(G._Pl, x, norm=Val(:schmidt)) # norm with √(2l+1)
        G.Gl[Gid] .+= G._Pl
    end
    end
    return nothing
end

function cal_Gτ(G::GreenFuncBin, τgrid, lmax::Int=-1)
    Gτs = [zeros(length(τgrid)) for i ∈ CartesianIndices(G.Gl)]
    Plx::Vector{f64} = G._Pl
    lmax::Int = lmax < 0 ? (length(Plx) - 1) : clamp(lmax, 0, length(Plx) - 1)
    for (iτ, τ) ∈ enumerate(τgrid)
        collectPl!(Plx, xτmap(τ, G.β), norm=Val(:schmidt)) # now they are order 0...lmax
        for (Gτ, Gl) ∈ zip(Gτs, G.Gl)
            @inbounds for lp1 ∈ 1:(lmax+1)
                Gτ[iτ] += Gl[lp1] * Plx[lp1]
            end
        end
    end
    Γ = inv(4 * G.Cw * G.insertion_trial)
    for Gτ ∈ Gτs
        Gτ .*= Γ
    end
    return Gτs
end

# two ways of normalizing greens function, as a cross check
function normalize_density_matrix(G::GreenFuncBin)
    @assert size(G.G0, 3) == size(G.G0, 4)
    Nsub = size(G.G0, 3)
    Γ = Nsub * inv(4 * G.Cw * G.insertion_trial)
    g = @. Γ*(real(G.G0))
    return g
end
function normalize_density_matrix(Dmat::Array{Complex{Int},4}, ishardcore=true)
    @assert size(Dmat,3) == size(Dmat,4)
    Nsub = size(Dmat, 3)
    d = [Dmat[1,1,i,i] for i ∈ 1:Nsub]
    D_NormCoeff::f64 = mean(imag, d) + (ishardcore ? +1 : -1) * mean(real, d)
    Dmat_n = let D = Dmat ./ D_NormCoeff
        for i ∈ 1:Nsub
            D[1, 1, i, i] = real(D[1, 1, i, i])
        end
        Dp = real.(D) .+ imag.(D)
        Dp
    end
    return Dmat_n
end






















# function normalize_density_matrix(Dmat, ishardcore=true)
#     d = Dmat[diagind(Dmat)]
#     D_NormCoeff::f64 = mean(imag, d) + (ishardcore ? +1 : -1) * mean(real, d)
#     Dmat_n = let D = Dmat ./ D_NormCoeff
#         Dp = real.(D)
#         for i ∈ CartesianIndices(D)
#             if i[1] ≠ i[2]
#                 Dp[i] += D[i].im
#             end
#         end
#         Dp
#     end
#     return Dmat_n
# end
# function translation_average(H::BH_Parameters,D::AbstractMatrix)
#     @assert size(D) |> allequal
#     Lx, Ly = H.Lx, H.Ly
#     N = Lx * Ly
#     L = (Lx, Ly)
#     Nsub, Rsub = divrem(size(D,1), N)
#     @assert Rsub == 0 "the dimension of H and D not match"
#     subids = Tuple.(CartesianIndices((Nsub, Nsub)))
#     Dab = [D[(a*N-N+1):(a*N),(b*N-N+1):(b*N)] for (a, b) ∈ subids]
#     display(Dab)
#     Cab = [zeros(L) for _ ∈ subids]
#     lattice = CartesianIndices(L)
#     lid = LinearIndices(L)
#     for (c,d) ∈ zip(Cab, Dab)
#         for (ic, Δr) ∈ enumerate(lattice)
#             for (iR, R) ∈ enumerate(lattice)
#                 c[ic] += d[iR, lid[mod1(R[1] + Δr[1] - 1, Lx), mod1(R[2] + Δr[2] - 1, Ly)]]
#             end
#         end
#         c ./= N
#     end
#     return Dict(zip(subids, Cab))
# end














# function translation_average(H::BH_Square, D)
#     Lx, Ly = H.Lx, H.Ly
#     @assert size(D, 1) == size(D, 2) == Lx * Ly
#     lattice = CartesianIndices((Lx, Ly))
#     lattice = ((Lx, Ly))

#     A = zeros(Lx, Ly)
#     for (ir, r) ∈ enumerate(CartesianIndices(A))
#         for (i, Ri) ∈ enumerate(CartesianIndices(A))
#             A[ir] += D[i, lattice[mod1(Ri[1] + r[1] - 1, Lx), mod1(Ri[2] + r[2] - 1, Ly)]]
#         end
#     end
#     return A ./ size(D,1)
# end
# function translation_average(H::BH_Pyroch, D)
#     Lx, Ly = H.Lx, H.Ly
#     N::Int = Lx * Ly
#     @assert size(D) = (2N,2N)
#     Dab = [
#         D[1:N, 1:N]
#         D[1:N, (N+1):2N]
#         D[(N+1):2N, 1:N]
#         D[(N+1):2N, (N+1):2N]
#     ]
#     A = [zeros(Lx, Ly) for i ∈ 1:4]
#     lattice = CartesianIndices((Lx, Ly))
#     lattice_lid = LinearIndices((Lx, Ly))
#     for (a,d) ∈ zip(A,Dab)
#         for (ir, r) ∈ enumerate(lattice)
#             for (i, Ri) ∈ enumerate(lattice)
#                 a[ir] += d[i, lattice_lid[mod1(Ri[1] + r[1]-1, Lx), mod1(Ri[2] + r[2], Ly)]]
#             end
#         end
#     end

#     return circshift(A, (1, 1)) ./ size(D, 1)
# end

# Dmat_n = Dmat |> normalize_density_matrix
# symmetric_sum(H, Dmat_n) ./ 9