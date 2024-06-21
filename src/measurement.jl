import Base: push!, empty!, merge
# function get_slice!(ψ::Matrix{T}, x::Wsheet, τ::f64 = 0.) where {T<:Real}
#     if iszero(τ)
#         for i ∈ eachindex(ψ, x.wl)
#             ψ[i] = x[i][end].n_L
#         end
#     else
#         for i ∈ eachindex(ψ, x.wl)
#             ψ[i] = element_around(x[i], τ, +1).n_L
#         end
#     end
#     return nothing
# end
# get_slice(li::Wline)::StateType = li[end].n_L
# get_slice(li::Wline, τ::f64)::StateType = element_around(li, τ, +1).n_L

function measure_site(l::Wline, U::f64, μ::f64)::NTuple{2,f64}
    @assert issorted(l)
    t::f64 = 0.
    n::StateType = l[end].n_R
    # β::f64 = l[end].t
    Δt::f64 = 0.
    U_val::f64 = 0.
    μ_val::f64 = 0.
    for e ∈ l
        Δt = e.t - t
        U_val += Δt * 0.5 * n * (n-1) * U
        μ_val -= Δt * n * μ
        t = e.t
        n = e.n_R
    end
    return (U_val, μ_val)
end
function merge_sorted!(C, A, B)::Nothing
    empty!(C)
    # @assert issorted(A) && issorted(B)
    iA::Int = 1
    lA::Int = lastindex(A)
    iB::Int = 1
    lB::Int = lastindex(B)
    @inbounds while true
        if iA ≤ lA && iB ≤ lB
            if A[iA] < B[iB]
                push!(C, A[iA])
                iA += 1
            else
                push!(C, B[iB])
                iB += 1
            end
        elseif iA ≤ lA
            push!(C, A[iA])
            iA += 1
        elseif iB ≤ lB
            push!(C, B[iB])
            iB += 1
        else
            return nothing
        end
    end
end
function measure_bond(li::Wline, lj::Wline, Vij::f64, bond_buffer::Wline)::f64
    t::f64 = 0.
    ni::StateType = li[end].n_R
    nj::StateType = lj[end].n_R
    i::IndexType = li[end].i
    j::IndexType = lj[end].i
    Δt::f64 = 0.0
    V_val::f64 = 0.
    merge_sorted!(bond_buffer, li, lj)
    for e ∈ bond_buffer
        if e.t > t
            Δt = e.t - t
            V_val += Δt * Vij * ni * nj
            t = e.t
        end
        if t == 1.0
            break
        elseif e.i == i
            ni = e.n_R
        elseif e.i == j
            nj = e.n_R
        else
            error("how can e.i ≠ i or j")
        end
    end
    return V_val
end
function simple_measure(x::Wsheet{N}, H::Union{BH_Square, BH_Pyroch}, bond_buffer::Wline) where {N}
    Ū::f64 = μ̄::f64 = V̄::f64 = K̄::f64 = 0.
    Nkink::Int = Npar::Int = 0
    # ρ::f64 = U::f64 = μ::f64 = V::f64 = K::f64 = 0.
    for i ∈ eachindex(x.wl)
        li::Wline = x[i]
        Ū, μ̄ = (Ū, μ̄) .+ measure_site(li, H.U, H.μ)
        for nb ∈ get_half_nbs(H, i)
            V̄ += measure_bond(li, x[nb], H.V, bond_buffer)
        end
        Nkink += (length(li)-1)
        Npar += li[end].n_R
    end
    @assert Nkink |> iseven
    K̄ = - (Nkink÷2) / x.β
    Ē = Ū + μ̄ + V̄ + K̄
    return (f64(Npar), f64(abs2(Npar)), Ē, K̄, Ū, μ̄, V̄)
end
@kwdef struct SimpleMeasure
    N::Accum{f64} = Accum(0.)
    N2::Accum{f64} = Accum(0.)
    E::Accum{f64} = Accum(0.)
    K::Accum{f64} = Accum(0.)
    U::Accum{f64} = Accum(0.)
    μ::Accum{f64} = Accum(0.)
    V::Accum{f64} = Accum(0.)
end
function push!(M::SimpleMeasure, one_measure::NTuple{7, f64})
    @inbounds for (i, m) ∈ enumerate(one_measure)
        push!(getfield(M, i)::Accum{f64}, m)
    end
    return nothing
end
function empty!(M::SimpleMeasure)
    for i ∈ 1:fieldcount(SimpleMeasure)
        empty!(getfield(M, i)::Accum{f64})
    end
    return nothing
end
function merge(M1::SimpleMeasure, M2::SimpleMeasure)
    return SimpleMeasure(
        Tuple(merge(getfield(M1,i), getfield(M2,i)) for i ∈ 1:fieldcount(SimpleMeasure))...)
end

### Measurement for the structure factor / Density Correlation
function rfft_buffer(L::NTuple{N,Int}) where {N}
    L = Tuple(i == 1 ? (1 + l >> 1) : l for (i, l) ∈ enumerate(L))
    return zeros(ComplexF64, L)
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

## Measurement for the winding_number
function winding_number(x::Wsheet{2}, H::BH_Square)::NTuple{2, Int}
    Wx::Int = Wy::Int = 0
    for i ∈ eachindex(x)
        for e::Element ∈ x[i]::Wline
            if e.op == b_
                dir = findfirst(==(e.j),get_nbs(H,i))
                Wy = Wy - (dir==1) + (dir==2)
                Wx = Wx - (dir==3) + (dir==4)
            end
        end
    end
    @assert Wx % H.Lx == Wy % H.Ly == 0
    WxWy = (Wx÷H.Lx, Wy÷H.Ly)
    return WxWy
end
function winding_number(x::Wsheet{3}, H::BH_Pyroch)::NTuple{2, Int}
    Wx::Int = Wy::Int = 0
    N::Int = H.Lx * H.Ly
    for i ∈ eachindex(x)
        for e::Element ∈ x[i]::Wline
            if e.op == b_
                dir = findfirst(==(e.j), get_nbs(H, i))
                # @assert dir isa Int
                sid =  i ≤ N ? 1 : 2
                if dir == 1
                    Wx += (-1)
                    Wy += (-1)
                elseif dir == 2
                    Wx += (+1)
                    Wy += (-1)
                elseif dir == 3
                    Wx += (-1)
                    Wy += (+1)
                elseif dir == 4
                    Wx += (+1)
                    Wy += (+1)
                elseif dir == 5
                    if sid == 1
                        Wy -= 2
                    else
                        Wx -= 2
                    end
                elseif dir == 6
                    if sid == 1
                        Wy += 2
                    else
                        Wx += 2
                    end
                else
                    error("dir illegal")
                end
            end
        end
    end
    @assert Wx % (2H.Lx) == Wy % (2H.Ly) == 0
    WxWy = (Wx ÷ (2H.Lx), Wy ÷ (2H.Ly))
    return WxWy
end

@kwdef struct WindingMeasure
    Wx2::Accum{Int} = Accum(0)
    Wy2::Accum{Int} = Accum(0)
    Wxy::Accum{Int} = Accum(0)
end
function push!(M::WindingMeasure, WxWy::NTuple{2,Int})
    Wx, Wy = WxWy
    push!(M.Wx2, abs2(Wx))
    push!(M.Wy2, abs2(Wy))
    push!(M.Wxy, abs(Wx*Wy))
    return nothing
end
function empty!(M::WindingMeasure)
    for i ∈ 1:fieldcount(WindingMeasure)
        empty!(getfield(M, i)::Accum{f64})
    end
    return nothing
end
function merge(M1::WindingMeasure, M2::WindingMeasure)
    return WindingMeasure(Tuple(merge(getfield(M1, i), getfield(M2, i)) for i ∈ 1:fieldcount(WindingMeasure))...)
end

struct WormMeasure
    simple::SimpleMeasure
    winding::WindingMeasure
    Gfunc::GreenFuncBin
    Sfact::StructureFactor2D
end
function WormMeasure(x::Wsheet, update_const::UpdateConsts;
    green_lmax::Int=100, is_full::Bool= false)::WormMeasure

    return WormMeasure(
        SimpleMeasure(),
        WindingMeasure(),
        GreenFuncBin(x, green_lmax, update_const.Cw; is_full = is_full),
        StructureFactor2D(x),
    )
end
function measure!(m::WormMeasure, x::Wsheet, H::BH_Parameters, bond_buffer, Sk_plan)
    tic::Int = time_ns()
    push!(m.simple, simple_measure(x, H, bond_buffer))
    push!(m.winding, winding_number(x, H))
    measure_Sk2D!(m.Sfact, x, Sk_plan)
    return UInt64(time_ns()-tic)
end

function Base.show(io::IO, m::T) where {T<:Union{SimpleMeasure, WindingMeasure}}
    println(io, "┌ $(T):")
    nterm = fieldcount(T)
    for (i,x) ∈ enumerate(fieldnames(T))
        println(io, i < nterm ? "│ " : "└ ", x, " = ", mean(getfield(m, x)))
    end
end
function Base.show(io::IO, m::StructureFactor2D{Nsub, NSk}) where {Nsub, NSk}
    println(io, "┌ StructureFactor2D:")
    println(io, "│ ψs (Nsub = $(Nsub)): $(summary(m.Sk_))")
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