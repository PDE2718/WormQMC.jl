
# function accum_green!(G::GreenFuncBin, w::Worm, H::T_Ham)::Nothing where {T_Ham<:BH_Parameters}
#     b::Element = b̂::Element = w.head::Element
#     if w.head.op == b_
#         b̂ = w.tail
#     else
#         b = w.tail
#     end
#     Gid = site_diff(H, b.i, b̂.i)
#     #measure density matrix
#     @inbounds begin
#         if w.loc == _at_green
#             G.G0[Gid] += complex(1, 0)
#         elseif w.loc == _at_stop
#             G.G0[Gid] += (b̂.t > b.t) ? complex(1, 0) : complex(0, 1)
#         elseif w.loc == _at_free && (G.is_full || Gid[1] == Gid[2] == 1)
#             Δτ = cal_Δτ(b.t, b̂.t, G.β)
#             x = xτmap(Δτ, G.β)
#             collectPl!(G._Pl, x, norm=Val(:schmidt)) # norm with √(2l+1)
#             G.Gl[Gid] .+= G._Pl
#         end
#     end
#     return nothing
# end
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