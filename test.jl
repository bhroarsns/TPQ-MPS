using ITensors
using Printf

include("./mps_tpq.jl")

# physical feature
n = 2 # system size
J = 1.0 # Heisenberg model parameter
# J = 4.0
# Γ = 2.0

# MPS parameters
χ = 2 # auxiliary site dimension (for future updates)
ξ = χ # bond dimension
# mTPQ parameters
l = 1.0
k_max = 2
# cTPQ parameters
temp_delta = 0.05
num_temps = 1

# prepare random MPS
ψ, sites = randomTPQMPS(n, χ, ξ)

# prepare MPO
h = heisenberg(sites, n, J)
h2 = h' * h
m = magnetization(sites, n)
m2 = m' * m

println("start calculating canonical summation")
# temperatures to calculate cTPQ
temps = (num_temps * temp_delta):-temp_delta:temp_delta

print("step 001")
khk = log(inner(ψ', h, ψ))
kh2k = log(inner(ψ'', h2, ψ))
kmk = log(inner(ψ', m, ψ))
km2k = log(inner(ψ'', m2, ψ))
kk = lognorm(ψ) * 2

ϕ = l * ψ - apply(h, ψ)
print("\tmTPQ state(k = 002) generated")
canonicalform!(ϕ, χ)
khk1 = log(inner(ψ', h, ϕ))
kh2k1 = log(inner(ψ'', h2, ϕ))
kmk1 = log(inner(ψ', m, ϕ))
km2k1 = log(inner(ψ'', m2, ϕ))
kk1 = logdot(ψ, ϕ)

energy = collect(Iterators.map(t -> [khk, khk1 + log(n) - log(t)], temps))
energy_2 = collect(Iterators.map(t -> [kh2k, kh2k1 + log(n) - log(t)], temps))
magnet = collect(Iterators.map(t -> [kmk, kmk1 + log(n) - log(t)], temps))
magnet_2 = collect(Iterators.map(t -> [km2k, km2k1 + log(n) - log(t)], temps))
beta_norm = collect(Iterators.map(t -> [kk, kk1 + log(n) - log(t)], temps))
factors = fill(0.0, length(temps))

for k = 1:k_max-1
    @printf "\rstep %03i\t" k+1
    global ψ = deepcopy(ϕ)
    local km2k = log(inner(ψ'', m2, ψ))
    local kmk = log(inner(ψ', m, ψ))
    local kh2k = log(inner(ψ'', h2, ψ))
    local khk = log(inner(ψ', h, ψ))
    local kk = lognorm(ψ) * 2
    # println(inner(ψ', h, ψ) / norm(ψ) / norm(ψ))

    global ϕ = l * ψ - apply(h, ψ)
    @printf "mTPQ state(k = %03i) generated" k+2
    canonicalform!(ϕ, χ)
    local km2k1 = log(inner(ψ'', m2, ϕ))
    local kmk1 = log(inner(ψ', m, ϕ))
    local kh2k1 = log(inner(ψ'', h2, ϕ))
    local khk1 = log(inner(ψ', h, ϕ))
    local kk1 = logdot(ψ, ϕ)

    for (i, t) in Iterators.enumerate(temps)
        factors[i] += 2*log(n) - 2*log(t) - log(2k-1) - log(2k)
        magnet_2[i] = cat(magnet_2[i], [km2k + factors[i], km2k1 + factors[i] + log(n) - log(t) - log(2k + 1)], dims=1)
        magnet[i] = cat(magnet[i], [kmk + factors[i], kmk1 + factors[i] + log(n) - log(t) - log(2k + 1)], dims=1)
        energy_2[i] = cat(energy_2[i], [kh2k + factors[i], kh2k1 + factors[i] + log(n) - log(t) - log(2k + 1)], dims=1)
        energy[i] = cat(energy[i], [khk + factors[i], khk1 + factors[i] + log(n) - log(t) - log(2k + 1)], dims=1)
        beta_norm[i] = cat(beta_norm[i], [kk + factors[i], kk1 + factors[i] + log(n) - log(t) - log(2k + 1)], dims=1)
    end
end
print("\ncanonical summation finished\n")

open("output.dat", "w") do io
    for (i, t) in Iterators.enumerate(temps)
        sort!(magnet_2[i]; by = x -> real(x), rev = true)
        sort!(magnet[i]; by = x -> real(x), rev = true)
        sort!(energy_2[i]; by = x -> real(x), rev = true)
        sort!(energy[i]; by = x -> real(x), rev = true)
        sort!(beta_norm[i]; by = x -> real(x), rev = true)
    
        divider = real(first(beta_norm[i]))
    
        mag2 = 0.0
        for e in reverse(magnet_2[i])
            mag2 += exp(e - divider)
        end

        mag = 0.0
        for e in reverse(magnet[i])
            mag += exp(e - divider)
        end

        ene2 = 0.0
        for e in reverse(energy_2[i])
            ene2 += exp(e - divider)
        end

        ene = 0.0
        for e in reverse(energy[i])
            ene += exp(e - divider)
        end
    
        bet_nor = 0.0
        for b in reverse(beta_norm[i])
            bet_nor += exp(b - divider)
        end

        dev_ene = ene2 / bet_nor - (ene / bet_nor)^2
        dev_mag = mag2 / bet_nor - (mag / bet_nor)^2
        
        write(io, "$t\t$(real(ene/bet_nor))\t$(real(mag/bet_nor))\t$(real(dev_ene)/t^2)\t$(n*real(dev_mag)/t)\n")
    end
end