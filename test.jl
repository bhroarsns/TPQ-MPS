using ITensors
using Printf

include("./tpq_mps.jl")

# physical feature
J = 1.0 # Heisenberg model parameter
# J = 4.0
# Γ = 2.0

# mTPQ parameters
l = 1.0
kmax = 500
# cTPQ parameters
tempstep = 0.05
numtemps = 80

# MPS parameters
n = 4 # system size
χ = 1 # auxiliary site dimension
ξ = 20 # bond dimension

# # prepare random MPS
# ψ, sites = randomTPQMPS(n, ξ, χ)

# # prepare MPO
# # l_h = l_heisenberg(sites, n, J, l)
# # h = heisenberg(sites, n, J)
# l_h = l_transverseising(sites, n, J, Γ, l)
# h = transverseising(sites, n, J, Γ)

# canonicalsummation(ψ, l_h, h, ξ, l; χ = χ, kmax=kmax, filename="ising/bigl")

for i in 1:50
    # prepare random MPS
    local ψ, sites = randomTPQMPS(n, ξ, χ)
    
    # prepare MPO
    local l_h = l_heisenberg(sites, n, J, l)
    local h = heisenberg(sites, n, J)
    # l_h = l_transverseising(sites, n, J, Γ, l)
    # h = transverseising(sites, n, J, Γ)

    canonicalsummation(ψ, l_h, h, ξ, l; χ = χ, kmax=kmax, filename="heisen4-$χ/run$i")
end

# χ = 20 # auxiliary site dimension

# for i in 1:10
#     # prepare random MPS
#     local ψ, sites = randomTPQMPS(n, ξ, χ)
    
#     # prepare MPO
#     local l_h = l_heisenberg(sites, n, J, l)
#     local h = heisenberg(sites, n, J)
#     # l_h = l_transverseising(sites, n, J, Γ, l)
#     # h = transverseising(sites, n, J, Γ)

#     canonicalsummation(ψ, l_h, h, ξ, l; χ = χ, kmax=kmax, filename="heisenberg-chi$χ/run$i")
# end

# χ = 40 # auxiliary site dimension

# for i in 1:10
#     # prepare random MPS
#     local ψ, sites = randomTPQMPS(n, ξ, χ)
    
#     # prepare MPO
#     local l_h = l_heisenberg(sites, n, J, l)
#     local h = heisenberg(sites, n, J)
#     # l_h = l_transverseising(sites, n, J, Γ, l)
#     # h = transverseising(sites, n, J, Γ)

#     canonicalsummation(ψ, l_h, h, ξ, l; χ = χ, kmax=kmax, filename="heisenberg-chi$χ/run$i")
# end