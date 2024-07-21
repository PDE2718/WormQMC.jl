import Base: push!, empty!, merge
# function simple_measure(x::Wsheet{N}, H::Union{BH_Square, BH_Pyroch, BH_Trimer}, bond_buffer::Wline) where {N}
#     Ū::f64 = μ̄::f64 = V̄::f64 = K̄::f64 = 0.
#     Nkink::Int = Npar::Int = 0
#     for i ∈ eachindex(x.wl)
#         li::Wline = x[i]
#         Ū, μ̄ = (Ū, μ̄) .+ measure_site(li, H.U, H.μ)
#         for nb ∈ get_nbs(H, i)
#             if nb < i
#                 V̄ += measure_bond(li, x[nb], H.V, bond_buffer)
#             end
#         end
#         Nkink += (length(li)-1)
#         Npar += li[end].n_R
#     end
#     @assert Nkink |> iseven
#     K̄ = - (Nkink÷2) / x.β
#     Ē = Ū + μ̄ + V̄ + K̄
#     return (f64(Npar), f64(abs2(Npar)), Ē, K̄, Ū, μ̄, V̄)
# end
# @kwdef mutable struct SimpleMeasure
#     N::Accum{f64} = Accum(0.)
#     N2::Accum{f64} = Accum(0.)
#     E::Accum{f64} = Accum(0.)
#     K::Accum{f64} = Accum(0.)
#     U::Accum{f64} = Accum(0.)
#     μ::Accum{f64} = Accum(0.)
#     V::Accum{f64} = Accum(0.)
# end

# function push!(M::SimpleMeasure, one_measure::NTuple{7, f64})
#     @inbounds for (i, m) ∈ enumerate(one_measure)
#         push!(getfield(M, i)::Accum{f64}, m)
#     end
#     return nothing
# end
# function empty!(M::SimpleMeasure)
#     for i ∈ 1:fieldcount(SimpleMeasure)
#         empty!(getfield(M, i)::Accum{f64})
#     end
#     return nothing
# end
# function merge(M1::SimpleMeasure, M2::SimpleMeasure)
#     return SimpleMeasure(
#         Tuple(merge(getfield(M1,i), getfield(M2,i)) for i ∈ 1:fieldcount(SimpleMeasure))...)
# end

### Measurement for the structure factor / Density Correlation

function SimpleMeasure(H::Ham) where {Ham<:BH_Parameters}
    names = simple_measure_names(H)
    return SimpleMeasure((x -> 0.0).(names), names, 0)
end

@kwdef mutable struct StructureFactor2D{Nsub,NSk}
    const ψs::NTuple{Nsub, Matrix{Float64}}
    const ψks::NTuple{Nsub, Matrix{ComplexF64}}
    const abinds::NTuple{NSk, Pair{Int, Int}}
    const Sk::NTuple{NSk, Matrix{ComplexF64}}
    const Sk_::Matrix{ComplexF64}
    n_measure::Int
end
function empty!(S::StructureFactor2D)
    for Ski ∈ S.Sk
        fill!(Ski, complex(0.0))
    end
    S.n_measure = 0
    return nothing
end

function rfft_buffer(L::NTuple{N,Int}) where {N}
    L = Tuple(i == 1 ? (1 + l >> 1) : l for (i, l) ∈ enumerate(L))
    return zeros(ComplexF64, L)
end
rfftplan(S::StructureFactor2D) = plan_rfft(S.ψs |> first)

function StructureFactor2D(Lx::Integer, Ly::Integer, Nsub::Integer)
    @assert Nsub > 0
    sz = Int.((Lx,Ly))
    ψs = Tuple(zeros(Float64, sz) for i ∈ 1:Nsub)
    ψks = Tuple(rfft_buffer(sz) for i ∈ 1:Nsub)
    NSk::Int = (Nsub*(Nsub+1))÷2
    abinds = Tuple(unique(Pair(minmax(a, b)...) for a ∈ 1:Nsub, b ∈ 1:Nsub))
    Sk = Tuple(rfft_buffer(sz) for i ∈ 1:NSk)
    Sk_ = rfft_buffer(sz)
    return StructureFactor2D{Nsub, NSk}(ψs, ψks, abinds, Sk, Sk_, 0)
end
StructureFactor2D(x::Wsheet{2}) = StructureFactor2D(size(x.wl)...,1)
StructureFactor2D(x::Wsheet{3}) = StructureFactor2D(size(x.wl)...)

function cal_Sk!(S::StructureFactor2D{Nsub,NSk}, P::FFTW.rFFTWPlan) where {Nsub,NSk}
    for (ψk, ψ) ∈ zip(S.ψks, S.ψs)
        mul!(ψk, P, ψ)
    end
    for ((a, b), Skab) ∈ zip(S.abinds, S.Sk)
        map!(dot, S.Sk_, S.ψks[a], S.ψks[b])
        Skab .+= S.Sk_
    end
    S.n_measure += 1
    return nothing
end

function get_slice!(S::StructureFactor2D{1,1}, x::Wsheet{2})
    map!(x -> f64(last(x).n_L), S.ψs[1], x.wl)
    return nothing
end
function get_slice!(S::StructureFactor2D, x::Wsheet{3})
    for i ∈ axes(x.wl, 3)
        map!(x -> f64(last(x).n_L), S.ψs[i], view(x.wl, :, :, i))
    end
    return nothing
end
function measure_Sk2D!(S::StructureFactor2D, x::Wsheet, P::FFTW.rFFTWPlan)
    get_slice!(S, x)
    cal_Sk!(S, P)
    return nothing
end
function Cab(S::StructureFactor2D, a::Int, b::Int)
    P = rfftplan(S)
    N = length(S.ψs[1])
    a_b = Pair(minmax(a,b)...)
    ab = findfirst(x->x==a_b, S.abinds)
    if a ≤ b
        return inv(N * S.n_measure)*(P \ S.Sk[ab])
    else
        return inv(N * S.n_measure)*(P \ conj.(S.Sk[ab]))
    end
end

# @kwdef struct WindingMeasure
#     Wx2::Accum{Int} = Accum(0)
#     Wy2::Accum{Int} = Accum(0)
#     Wxy::Accum{Int} = Accum(0)
# end
# function push!(M::WindingMeasure, WxWy::NTuple{2,Int})
#     Wx, Wy = WxWy
#     push!(M.Wx2, abs2(Wx))
#     push!(M.Wy2, abs2(Wy))
#     push!(M.Wxy, abs(Wx*Wy))
#     return nothing
# end
# function empty!(M::WindingMeasure)
#     for i ∈ 1:fieldcount(WindingMeasure)
#         empty!(getfield(M, i)::Accum{f64})
#     end
#     return nothing
# end
# function merge(M1::WindingMeasure, M2::WindingMeasure)
#     return WindingMeasure(Tuple(merge(getfield(M1, i), getfield(M2, i)) for i ∈ 1:fieldcount(WindingMeasure))...)
# end

struct WormMeasure
    simple::SimpleMeasure
    Gfunc::GreenFuncBin
    Sfact::StructureFactor2D
end
function WormMeasure(x::Wsheet, H::BH_Parameters, update_const::UpdateConsts;
    green_lmax::Int=100, is_full::Bool= false)::WormMeasure
    return WormMeasure(
        SimpleMeasure(H),
        GreenFuncBin(x, green_lmax, update_const.Cw; is_full = is_full),
        StructureFactor2D(x),
    )
end
function measure!(m::WormMeasure, x::Wsheet{Nwl}, H::Ham,
    bond_buffer::Wline, Sk_plan::FFTW.rFFTWPlan) where {Nwl,Ham<:BH_Parameters}
    tic::Int = time_ns()
    simple_measure!(m.simple, x, H, bond_buffer)
    measure_Sk2D!(m.Sfact, x, Sk_plan)
    return UInt64(time_ns()-tic)
end

function Base.show(io::IO, m::StructureFactor2D{Nsub, NSk}) where {Nsub, NSk}
    println(io, "┌ StructureFactor2D:")
    println(io, "│ ψs (Nsub = $(Nsub)): $(summary(m.ψs[1]))")
    println(io, "│ Sk (NSk  = $(NSk)): $(summary(m.Sk_))")
    println(io, "└ indices : $(m.abinds)")
end
function Base.show(io::IO, m::WormMeasure)
    T = WormMeasure
    println(io, T, ":")
    for (i, x) ∈ enumerate(fieldnames(T))
        show(io, getfield(m, x))
    end
end

function merge(s1::SimpleMeasure{Np}, s2::SimpleMeasure{Np}
    )::SimpleMeasure{Np} where {Np}
    @assert s1.names == s2.names
    return SimpleMeasure(
        s1.props .+ s2.props,
        s1.names,
        s1.n_measure + s2.n_measure
    )
end

function merge(s1::StructureFactor2D{Nsub,NSk}, s2::StructureFactor2D{Nsub,NSk}
    )::StructureFactor2D{Nsub,NSk} where {Nsub,NSk}

    @assert s1.abinds == s2.abinds
    @assert size(first(s1.ψs)) == size(first(s2.ψs))
    abinds = deepcopy(s1.abinds)
    ψs = ntuple(i -> similar(s1.ψs[i]), Nsub)
    ψks = ntuple(i -> similar(s1.ψks[i]), Nsub)
    Sk_ = similar(s1.Sk_)
    Sk = ntuple(n -> s1.Sk[n] + s2.Sk[n], NSk)
    num = s1.n_measure + s2.n_measure
    return StructureFactor2D(ψs, ψks, abinds, Sk, Sk_, num)
end
function merge(g1::GreenFuncBin, g2::GreenFuncBin)::GreenFuncBin
    @assert g1.is_full == g2.is_full
    @assert g1.β == g2.β
    @assert g1.Cw == g2.Cw
    @assert length(g1._Pl) == length(g2._Pl)
    @assert size(g1.Gl) == size(g2.Gl)
    @assert size(g1.G0) == size(g2.G0)
    _Pl = similar(g1._Pl)
    Gl = similar(g1.Gl)
    for i ∈ eachindex(Gl, g1.Gl, g2.Gl)
        Gl[i] = g1.Gl[i] + g2.Gl[i]
    end
    G0 = g1.G0 + g2.G0
    num = g1.insertion_trial + g2.insertion_trial
    return GreenFuncBin(g1.is_full, g1.β, g1.Cw, _Pl, Gl, G0, num)
end

function merge(m1::WormMeasure, m2::WormMeasure)::WormMeasure
    return WormMeasure(
        merge(m1.simple, m2.simple),
        merge(m1.Gfunc, m2.Gfunc),
        merge(m1.Sfact, m2.Sfact)
    )
end