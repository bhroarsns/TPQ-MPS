using ITensors
using Printf

#calculation parameters
n = 16 # system size
χ = 20 # bond dimension
J = 1.0
l = 1.0
k_max = 100
temps = 4.0:-0.05:0.05

# prepare random MPS
sites = siteinds("S=1/2", n)
sizehint!(sites, n + 2)
pushfirst!(sites, Index(χ))
push!(sites, Index(χ))
ψ = randomMPS(Complex{Float64}, sites, χ)
ψ[1] = delta(Complex{Float64}, sites[1], linkinds(ψ, 1))
ψ[n + 2] = delta(Complex{Float64}, linkinds(ψ, n + 1), sites[n + 2])
print("initial state successfully generated: norm = $(norm(ψ)) equals to sqrt(χ) = $(sqrt(χ))\n")

# prepare MPO
h = let
    JJ = J / n
    ampo = OpSum()
    for j = 2:n
        ampo += JJ / 2.0, "S+", j, "S-", j+1
        ampo += JJ / 2.0, "S-", j, "S+ ", j+1
        ampo += JJ, "Sz", j, "Sz", j+1
    end
    MPO(ampo, sites)
end

h2 = h' * h

m = let
    ampo = OpSum()
    for j = 2:n+1
        ampo += 1/n,"Sz", j
    end
    MPO(ampo, sites)
end

m2 = m' * m

print("matrix product operators successfully generated\n")

function canonical_form!(A::MPS, max)
    len = length(A)
    sites = siteinds(A)

    if len != 1
        U1, s1, V1 = LinearAlgebra.svd(A[1] * A[2], sites[1])
        prev_link, _ = inds(s1)
        A[1] = U1
        A[2] = s1 * V1

        for j = 2:len-1
            Uj, sj, Vj = LinearAlgebra.svd(A[j] * A[j + 1], sites[j], prev_link)
            prev_link, _ = inds(sj)
            A[j] = Uj
            A[j + 1] = sj * Vj
        end

        Vn, sn, Un = LinearAlgebra.svd(A[len - 1] * A[len], sites[len]; lefttags = "Link,l=$(len-1)", maxdim = max)
        prev_link, _ = inds(sn)
        A[len] = Vn
        A[len - 1] = Un * sn

        for j = len-1:-1:3
            Vj, sj, Uj = LinearAlgebra.svd(A[j - 1] * A[j], sites[j], prev_link; lefttags = "Link,l=$(j-1)", maxdim = max)
            prev_link, _ = inds(sj)
            A[j] = Vj
            A[j - 1] = Uj * sj
        end
    end

    # print("\tcanonicalization finished")
end

print("start calculating canonical summation\n")

print("step 001")
km2k = log(inner(ψ'', m2, ψ))
kmk = log(inner(ψ', m, ψ))
kh2k = log(inner(ψ'', h2, ψ))
khk = log(inner(ψ', h, ψ))
kk = lognorm(ψ) * 2

ϕ = l * ψ - apply(h, ψ)
print("\tmTPQ state(k = 002) generated")
canonical_form!(ϕ, χ)
km2k1 = log(inner(ψ'', m2, ϕ))
kmk1 = log(inner(ψ', m, ϕ))
kh2k1 = log(inner(ψ'', h2, ϕ))
khk1 = log(inner(ψ', h, ϕ))
kk1 = logdot(ψ, ϕ)

energy_2 = collect(Iterators.map(t -> [kh2k, kh2k1 + log(n) - log(t)], temps))
energy = collect(Iterators.map(t -> [khk, khk1 + log(n) - log(t)], temps))
magnet_2 = collect(Iterators.map(t -> [km2k, km2k1 + log(n) - log(t)], temps))
magnet = collect(Iterators.map(t -> [kmk, kmk1 + log(n) - log(t)], temps))
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
    canonical_form!(ϕ, χ)
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