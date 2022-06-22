using ITensors

#calculation parameters
n = 16 # system size
χ = 20 # bond dimension
J = 1.0
l = 1.0
k_max = 500
temps = 4.0:-0.05:3.95

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
print("energy density operator successfully generated\n")

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

        for j = len-1:-1:2
            Vj, sj, Uj = LinearAlgebra.svd(A[j - 1] * A[j], sites[j], prev_link; lefttags = "Link,l=$(j-1)", maxdim = max)
            prev_link, _ = inds(sj)
            A[j] = Vj
            A[j - 1] = Uj * sj
        end
    end

    print("\tcanonicalization finished")
end

print("start calculating canonical summation\n")
print("step 1")
khk = log(inner(ψ', h, ψ))
kk = lognorm(ψ)

ϕ = l * ψ - apply(h, ψ)
print("\tmTPQ state(k = 2) generated")
canonical_form!(ϕ, χ)
khk1 = log(inner(ψ', h, ϕ))
kk1 = logdot(ψ, ϕ)

energy_density = Iterators.map(t -> [khk, khk1 + log(n) - log(t)], temps)
beta_norm = Iterators.map(t -> [kk, kk1 + log(n) - log(t)], temps)
factors = fill(0.0, length(temps))

for k = 1:k_max-1
    print("\rstep $(k+1)\t")
    global ψ = deepcopy(ϕ)
    local khk = log(inner(ψ', h, ψ))
    local kk = lognorm(ψ)

    global ϕ = l * ψ - apply(h, ψ)
    print("mTPQ state(k = $(k+2)) generated")
    canonical_form!(ϕ, χ)
    local khk1 = log(inner(ψ', h, ϕ))
    local kk1 = logdot(ψ, ϕ)

    for (i, t) in Iterators.enumerate(temps)
        factors[i] += 2*log(n) - 2*log(t) - log(2k-1) - log(2k)
        energy_density[i] = cat(energy_density[i], [khk + factors[i], khk1 + factors[i] + log(n) - log(t) - log(2k + 1)], dims=1)
        beta_norm[i] = cat(beta_norm[i], [kk + factors[i], kk1 + factors[i] + log(n) - log(t) - log(2k + 1)], dims=1)
    end
end
print("\ncanonical summation finished\n")

open("output.dat", "w") do io
    for (i, t) in Iterators.enumerate(temps)
        sort!(energy_density[i]; by = x -> real(x), rev = true)
        sort!(beta_norm[i]; by = x -> real(x), rev = true)
    
        divider = min(real(first(energy_density[i])), real(first(beta_norm[i])))
    
        ene_den = 0.0
        for e in reverse(energy_density[i])
            ene_den += exp(e - divider)
        end
    
        bet_nor = 0.0
        for b in reverse(beta_norm[i])
            bet_nor += exp(b - divider)
        end
        
        write(io, "$t\t$(ene_den / bet_nor)\n")
    end
end