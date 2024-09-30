abstract type BH_Parameters end
@kwdef struct BH_Trim1D <: BH_Parameters
    nmax::StateType = 1 # nmax must be implemented
    Lx::IndexType = 0 # the size should also be implemented
    Ly::IndexType = 1 # So far only 2D is supported
    J::f64 = 1.0 # hopping in the parallel direction
    R::f64 = 1.0 # strength of rotation
    V::f64 = 10.0 # off site repulsion, set to a large value for shaped hard core
    μ1::f64 = 0.0 # chemical potential A
    μ2::f64 = 0.0 # chemical potential B
end
N_sublatt(::Type{BH_Trim1D}) = 2
N_wldim(::Type{BH_Trim1D})::Int = 3
N_nbs(::Type{BH_Trim1D})::Int = 7
N_hops(::Type{BH_Trim1D})::Int = 3

function get_nbs(H::BH_Trim1D, i)::NTuple{7,Int}
    L = Int(H.Lx)
    x, y, s = CartesianIndices((L, 1, 2))[i] |> Tuple
    if s == 1
        return (
            mod1(x - 2, L), mod1(x - 1, L), mod1(x + 1, L), mod1(x + 2, L),
            mod1(x - 1, L) + L, x + L, mod1(x + 1, L) + L
        )
    else # s0 == 2
        return (mod1(x - 1, L), x, mod1(x + 1, L), 0, 0, 0, 0)
    end
end
function get_hops(H::BH_Trim1D, i)::NTuple{3,Int}
    L = Int(H.Lx)
    x, y, s = CartesianIndices((L, 1, 2))[i] |> Tuple
    if s == 1
        return (mod1(x - 1, L), mod1(x + 1, L), x + L) # J,J,R
    else # s0 == 2
        return (0, 0, x)
    end
end
function diagE(H::BH_Trim1D, i, ni::StateType, njs::NTuple{7,StateType})::f64
    μ = i > H.Lx ? H.μ2 : H.μ1
    return ni * (H.V * sum(njs) - μ)
end

function bond_weight(H::BH_Trim1D, i, j)::f64
    C = CartesianIndices((H.Lx, H.Ly, 2))
    (C[i][3] == C[j][3] == 1) ? H.J : H.R
end

function site_diff(H::BH_Trim1D, i, j)::CartesianIndex{4}
    lattice = CartesianIndices((H.Lx, H.Ly, 2))
    xi, yi, si = Tuple(lattice[i])
    xj, yj, sj = Tuple(lattice[j])
    dx = mod1(xj - xi + 1, H.Lx)
    dy = mod1(yj - yi + 1, H.Ly)
    return CartesianIndex(dx, dy, si, sj)
end

function Wsheet(β::f64, H::BH_Trim1D)::Wsheet{3}
    @assert H.Ly == 1
    Wsheet(β, zeros(StateType, H.Lx, H.Ly, 2))
end

function simple_measure!(m::SimpleMeasure, x::Wsheet{3}, H::BH_Trim1D,
    bond_buffer::Wline)::Nothing
    μ = V = 0.0
    N₁ = N₂ = 0
    Ns = H.Lx * H.Ly
    Wx = NK_J = NK_R = 0
    lattice = CartesianIndices(x.wl)
    for (i, c) ∈ enumerate(lattice)
        li::Wline = x[i]
        μc::f64 = c[3]==1 ? H.μ1 : H.μ2
        μ += measure_site(li, 0.0, μc)[2]
        for nb ∈ get_nbs(H, i)
            if 0 < nb < i
                V += measure_bond(li, x[nb], H.V, bond_buffer)
            end
        end
        ni = li[end].n_R
        if c[3] == 1
            N₁ += ni
        else #if c[3] == 2
            N₂ += ni
        end

        for e ∈ li
            if e.op == b_
                j = Int(e.j)
                if c[3] == 1
                    hops = get_hops(H, i)
                    if j == hops[1]
                        Wx -= 1
                        NK_J += 1
                    elseif j == hops[2]
                        Wx += 1
                        NK_J += 1
                    else
                        @assert j == hops[3]
                        NK_R += 1
                    end
                else # c[3] == 2
                    @assert 0 < j ≤ Ns
                    NK_R += 1
                end
            end
        end
    end

    @assert Wx % H.Lx == 0
    Wx ÷= H.Lx
    J = -NK_J / x.β
    R = -NK_R / x.β
    N = N₁ + N₂
    K = J + R
    E = μ + V + J + R
    Wx² = abs2(Wx)
    m.props = m.props .+ (E,
        N, abs2(N), N₁, abs2(N₁), N₂, abs2(N₂),
        K, J, R, μ, V, abs2(Wx)
    )
    m.n_measure += 1
    return nothing
end

function simple_measure_names(::BH_Trim1D)
    return (:E,
        :N, :N², :N₁, :N₁², :N₂, :N₂², 
        :K, :J, :R, :μ, :V, :Wx²
    )
end