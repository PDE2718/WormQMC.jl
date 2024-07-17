# Example: all the functions need to be implement for a 2D lattice:
# parameters for 2D Square lattice with only J1
@kwdef struct BH_Square <: BH_Parameters
    nmax::StateType = 1 # nmax must be implemented
    Lx::IndexType = 0 # the size should also be implemented
    Ly::IndexType = 0 # 
    J::f64 = 1.0 # other parameters
    U::f64 = 0.0
    V::f64 = 1.0
    μ::f64 = 0.0
end
N_sublatt(::Type{BH_Square}) = 1
N_wldim(::Type{BH_Square})::Int = 2
N_nbs(::Type{BH_Square})::Int = 4
N_hops(::Type{BH_Square})::Int = 4

function get_nbs(H::BH_Square, i::Integer)::NTuple{4,Int}
    @inbounds begin
        Lx = Int(H.Lx)
        Ly = Int(H.Ly)
        x0, y0 = CartesianIndices((Lx, Ly))[i] |> Tuple
        xm, xp = mod1.((x0 - 1, x0 + 1), Lx)
        ym, yp = mod1.((y0 - 1, y0 + 1), Ly)
        idmap = LinearIndices((Lx, Ly))
        return (
            idmap[x0, ym],
            idmap[x0, yp],
            idmap[xm, y0],
            idmap[xp, y0],
        )
    end
end
get_hops(H::BH_Square, i::Integer)::NTuple{4,Int} = get_nbs(H, i)

function diagE(H::BH_Square, ni::StateType, njs::NTuple{4,StateType})::f64
    return 0.5 * H.U * ni * (ni - 1) - H.μ * ni + H.V * ni * sum(njs)
end
bond_weight(H::BH_Square, i::Integer, j::Integer)::f64 = H.J
Wsheet(β::f64, H::BH_Square) = Wsheet(β, zeros(StateType, H.Lx, H.Ly))
function site_diff(H::BH_Square, i, j)::CartesianIndex{4}
    lattice = CartesianIndices((H.Lx, H.Ly))
    xi, yi = Tuple(lattice[i])
    xj, yj = Tuple(lattice[j])
    dx = mod1(xj - xi + 1, H.Lx)
    dy = mod1(yj - yi + 1, H.Ly)
    return CartesianIndex(dx, dy, 1, 1)
end

## Measurement for the winding_number
function winding_number(x::Wsheet{2}, H::BH_Square)::NTuple{2,Int}
    Wx::Int = Wy::Int = 0
    for i ∈ eachindex(x)
        for e::Element ∈ x[i]::Wline
            if e.op == b_
                hops = get_hops(H, i)
                j = Int(e.j)
                if j == hops[1]
                    Wy -= 1
                elseif j == hops[2]
                    Wy += 1
                elseif j == hops[3]
                    Wx -= 1
                elseif j == hops[4]
                    Wx += 1
                else
                    error("illegal operator")
                end
                # Wy = Wy - (dir == 1) + (dir == 2)
                # Wx = Wx - (dir == 3) + (dir == 4)
            end
        end
    end
    @assert Wx % H.Lx == Wy % H.Ly == 0
    WxWy = (Wx ÷ H.Lx, Wy ÷ H.Ly)
    return WxWy
end
# function get_half_nbs(H::BH_Pyroch, i::Integer)::NTuple{3,Int}
#     nbs = get_nbs(H, i)
#     return (nbs[1], nbs[2], nbs[5])
# end
