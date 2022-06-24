using ITensors
using Printf

include("./tpq_mps.jl")

# physical feature
n = 8 # system size
J = 1.0 # Heisenberg model parameter
# J = 4.0
# Γ = 2.0

# MPS parameters
χ = 20 # auxiliary site dimension (for future updates)
ξ = χ # bond dimension
# mTPQ parameters
l = 1.0
kmax = 500
# cTPQ parameters
tempstep = 0.05
numtemps = 80

# prepare random MPS
ψ, sites = randomTPQMPS(n, χ, ξ)

# prepare MPO
h = heisenberg(sites, n, J)
m = magnetization(sites, n)

canonicalsummation(ψ, h, ξ, χ, l, kmax, tempstep, numtemps, m)

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
normalize!(ψ)

h2 = h' * h
m2 = m' * m

mtpq_data = zeros(Complex{Float64}, kmax+1, 2)

for k = 0:kmax
    if k % 50 == 0
        @printf "\rstep %03i" k
    end

    khk = inner(ψ', h, ψ)
    kh2k = inner(ψ'', h2, ψ)
    kmk = log(inner(ψ', m, ψ))
    km2k = log(inner(ψ'', m2, ψ))

    mtpq_data[k+1, 1] = khk
    mtpq_data[k+1, 2] = kh2k

    ϕ = l * ψ - apply(h, ψ)
    # @printf "\tmTPQ state(k = %03i) generated" k+1
    canonicalform!(ϕ, ξ, χ)
    # print("\tcanonicalization finished")

    khk1 = log(l * khk - kh2k) # ⟨k|h|k+1⟩ = ⟨k|h(l - h)|k⟩ = l⟨k|h|k⟩ - ⟨k|h²|k⟩
    kh2k = log(kh2k)
    kh2k1 = log(inner(ψ'', h2, ϕ))
    kmk1 = log(inner(ψ', m, ϕ))
    km2k1 = log(inner(m, ψ, m, ϕ))
    kk1 = log(l - khk) # ⟨k|k+1⟩ = ⟨k|(l - h)|k⟩ = l⟨k|k⟩ - ⟨k|h|k⟩ = l⟨k|k⟩ - ⟨k|h|k⟩
    khk = log(khk)

    ψ = deepcopy(ϕ)
    norm_k = []
    normalize!(ψ; (lognorm!)=norm_k)

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
        factors[i, k+2] = new_factor + log(n) - log(t) - log(2k + 2) + 2 * norm_k[1]
    end
end

open(string(filename, "_mtpq.dat"), "w") do io
    for k in 1:kmax+1
        write(io, "$(real(n*(l - mtpq_data[k, 1])/2(k-1)))\t$(real(mtpq_data[k, 1]))\t$(real(mtpq_data[k, 2] - mtpq_data[k, 1]^2))\n")
    end
end

println("\ncanonical summation finished")

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
        
        write(io, "$t\t$(real(ene))\t$(real(mag))\t$(n*real(dev_ene)/t^2)\t$(n*real(dev_mag)/t)\n")
    end
end