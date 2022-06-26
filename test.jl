using ITensors
using Printf

include("./tpq_mps.jl")

# physical feature
n = 4 # system size
J = 1.0 # Heisenberg model parameter
# J = 4.0
# Γ = 2.0

# MPS parameters
χ = 20 # auxiliary site dimension (for future updates)
ξ = 20 # bond dimension
# mTPQ parameters
l = 1.0
kmax = 500
# cTPQ parameters
tempstep = 0.05
numtemps = 80

# prepare random MPS
ψ, sites = randomTPQMPS(n, χ, ξ)

# prepare MPO
l_h = l_heisenberg(sites, n, J, l)
h = heisenberg(sites, n, J)
# l_h = l_transverseising(sites, n, J, Γ, l)
# h = transverseising(sites, n, J, Γ)
m = magnetization(sites, n)

canonicalsummation(ψ, h, ξ, l)

filename = "output"
println("start calculating canonical summation")

n = length(ψ) - 2

energy = zeros(Complex{Float64}, numtemps, 2kmax + 2)
energy2 = zeros(Complex{Float64}, numtemps, 2kmax + 2)
magnet = zeros(Complex{Float64}, numtemps, 2kmax + 2)
magnet2 = zeros(Complex{Float64}, numtemps, 2kmax + 2)
beta_norm = zeros(Complex{Float64}, numtemps, 2kmax + 2)
factors = zeros(Complex{Float64}, numtemps, kmax + 2)

# temperatures to calculate cTPQ
temps = (numtemps * tempstep):-tempstep:tempstep
println("temperature range: from $(tempstep) to $(numtemps * tempstep)")

# ψ = deepcopy(init)
h2 = h' * h
m2 = m' * m
# h2 = contract(h', h, cutoff=-Inf)
# m2 = contract(m', m, cutoff=-Inf)

mtpq_data = zeros(Complex{Float64}, kmax+1, 2)

for k = 0:kmax
    @printf "\rstep %03i" k
    
    khk = log(inner(ψ', h, ψ))
    kh2k = log(inner(ψ'', h2, ψ))
    kmk = log(inner(ψ', m, ψ))
    km2k = log(inner(ψ'', m2, ψ))

    mtpq_data[k+1, 1] = exp(khk)
    mtpq_data[k+1, 2] = exp(kh2k)
    
    global ϕ = apply(l_h, ψ)
    canonicalform!(ϕ, ξ, χ)

    khk1 = log(inner(ψ', h, ϕ))
    kh2k1 = log(inner(ψ'', h2, ϕ))
    kmk1 = log(inner(ψ', m, ϕ))
    km2k1 = log(inner(ψ'', m2, ϕ))
    kk1 = logdot(ψ, ϕ)

    global ψ = deepcopy(ϕ)
    kk = norm(ψ)
    normalize!(ψ)

    for (i, t) in enumerate(temps)
        factor = factors[i, k+1]
        new_factor = factor + log(n) - log(t) - log(2k + 1)
        energy[i, 2k+1] = khk + factor
        energy[i, 2k+2] = khk1 + new_factor
        energy2[i, 2k+1] = kh2k + factor
        energy2[i, 2k+2] = kh2k1 + new_factor
        magnet[i, 2k+1] = kmk + factor
        magnet[i, 2k+2] = kmk1 + new_factor
        magnet2[i, 2k+1] = km2k + factor
        magnet2[i, 2k+2] = km2k1 + new_factor
        beta_norm[i, 2k+1] = factor
        beta_norm[i, 2k+2] = kk1 + new_factor
        factors[i, k+2] = new_factor + log(n) - log(t) - log(2k + 2) + 2.0 * log(kk)
    end
end

open(string(filename, "_mtpq.dat"), "w") do io
    for k in 1:kmax+1
        temp = real(n*(l - mtpq_data[k, 1])/2(k-1))
        ene = real(mtpq_data[k, 1])
        ene2 = real(mtpq_data[k, 2])

        dev_ene = ene2 - ene^2

        write(io, "$temp\t$(ene/n)\t$(dev_ene  / temp^2 / n)\n")
    end
end

open(string(filename, "_ctpq.dat"), "w") do io
    for (i, t) in enumerate(temps)
        energyt = sort(energy[i, :]; by = x -> real(x))
        energy2t = sort(energy2[i, :]; by = x -> real(x))
        magnett = sort(magnet[i, :]; by = x -> real(x))
        magnet2t = sort(magnet2[i, :]; by = x -> real(x))
        beta_normt = sort(beta_norm[i, :]; by = x -> real(x))

        divider = real(last(beta_normt))
    
        ene = 0.0
        for e in energyt
            ene += exp(e - divider)
        end

        ene2 = 0.0
        for e in energy2t
            ene2 += exp(e - divider)
        end

        mag = 0.0
        for e in magnett
            mag += exp(e - divider)
        end

        mag2 = 0.0
        for e in magnet2t
            mag2 += exp(e - divider)
        end
    
        bet_nor = 0.0
        for b in beta_normt
            bet_nor += exp(b - divider)
        end

        ene = ene / bet_nor
        mag = mag / bet_nor

        dev_ene = ene2 / bet_nor - ene^2
        dev_mag = mag2 / bet_nor - mag^2
        
        write(io, "$t\t$(real(ene) / n)\t$(real(mag) / n)\t$(real(dev_ene) / t^2 / n)\t$(real(dev_mag) / t / n)\n")
    end
end
