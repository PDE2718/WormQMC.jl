mutable struct SimpleMeasure{Np}
    props::NTuple{Np,f64}
    const names::NTuple{Np,Symbol}
    n_measure::Int
end
import Base:empty!
function empty!(m::SimpleMeasure)
    m.props = 0.0 .* m.props
    m.n_measure = 0
    return m
end
import Base:show
function Base.show(io::IO, m::SimpleMeasure)
    println(io, "SimpleMeasure: ", NamedTuple(zip(m.names, m.props ./ m.n_measure)))
end

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


