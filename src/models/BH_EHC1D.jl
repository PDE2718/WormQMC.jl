abstract type BH_Parameters end
@kwdef struct BH_EHC1D <: BH_Parameters
    nmax::StateType = 1 # nmax must be implemented
    Lx::IndexType = 0 # the size should also be implemented
    Ly::IndexType = 1 # So far only 2D is supported
    J::f64 = 1.0 # hopping in the parallel direction
    V::f64 = 0.0 # off site repulsion, set to a large value for shaped hard core
    μ::f64 = 0.0 # chemical potential A
end
N_sublatt(::Type{BH_EHC1D}) = 1
N_wldim(::Type{BH_EHC1D})::Int = 2 # can only be 2 or 3
N_nbs(::Type{BH_EHC1D})::Int = 2
N_hops(::Type{BH_EHC1D})::Int = 2

function get_nbs(H::BH_EHC1D, i::Integer)::NTuple{2,Int}
    return mod1.((i - 1, i + 1), H.Lx)
end
function get_hops(H::BH_EHC1D, i::Integer)::NTuple{2,Int}
    return mod1.((i - 1, i + 1), H.Lx)
end
function diagE(H::BH_EHC1D, i::Integer, ni::StateType, njs::NTuple{2,StateType})::f64
    return ni * (H.V * sum(njs) - H.μ)
end
function bond_weight(H::BH_EHC1D, i::Integer, j::Integer)::f64
    return H.J
end
function site_diff(H::BH_EHC1D, i::Integer, j::Integer)::CartesianIndex{4}
    dx = mod1(j - i + 1, H.Lx)
    return CartesianIndex(dx, 1, 1, 1)
end
function Wsheet(β::f64, H::BH_EHC1D)
    @assert H.Ly == 1
    Wsheet(β, zeros(StateType, H.Lx, H.Ly))
end

function simple_measure!(m::SimpleMeasure, x::Wsheet{2}, H::BH_EHC1D,
    bond_buffer::Wline)::Nothing
    μ = V = 0.0
    N = 0
    Wx = NK_J = 0
    lattice = CartesianIndices(x.wl)
    for (i, c) ∈ enumerate(lattice)
        li::Wline = x[i]
        μ += measure_site(li, 0.0, H.μ)[2]
        for nb ∈ get_nbs(H, i)
            if 0 < nb < i
                V += measure_bond(li, x[nb], H.V, bond_buffer)
            end
        end
        ni = li[end].n_R
        N += ni
        for e ∈ li
            if e.op == b_
                j = Int(e.j)
                hops = get_hops(H, i)
                NK_J += 1
                if j == hops[1]
                    Wx -= 1 
                elseif j == hops[2]
                    Wx += 1
                else
                    error("something wrong with nb index")
                end
            end
        end
    end
    @assert Wx % H.Lx == 0
    Wx ÷= H.Lx
    Wx² = abs2(Wx)
    N² = abs2(N)
    K = -NK_J / x.β
    E = μ + V + K
    m.props = m.props .+ (E, N, N², K, μ, V, Wx²)
    m.n_measure += 1
    return nothing
end

function simple_measure_names(::BH_EHC1D)
    return (:E, :N, :N², :K,  :μ, :V, :Wx²)
end