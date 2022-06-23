export heisenberg, transverseising, magnetization, randomTPQMPS, canonicalform!, canonicalsummation

using ITensors

function heisenberg(sites, n, J)
    JJ = J / n
    mpo = OpSum()
    for j = 2:n
        mpo += JJ / 2.0, "S+", j, "S-", j+1
        mpo += JJ / 2.0, "S-", j, "S+ ", j+1
        mpo += JJ, "Sz", j, "Sz", j+1
    end
    println("MPO of the energy density operator of the Heisenberg model generated")
    return MPO(mpo, sites)
end

function transverseising(sites, n, J, Γ)
    JJ = J / n
    ΓΓ = Γ / n
    mpo = OpSum()
    for j = 2:n
        mpo += JJ, "Sz", j, "Sz", j+1
        mpo += ΓΓ, "Sx", j
    end
    mpo += ΓΓ, "Sx", n + 1
    println("MPO of the energy density operator of the transverse Ising model generated")
    return MPO(mpo, sites)
end

function magnetization(sites, n)
    mpo = OpSum()
    for j = 2:n+1
        mpo += 1/n,"Sz", j
    end
    println("MPO of magnetization generated")
    return MPO(mpo, sites)
end

function randomTPQMPS(system_size, χ, ξ=χ, sitetype="S=1/2")
    sites = siteinds(sitetype, system_size)
    sizehint!(sites, system_size + 2)
    pushfirst!(sites, Index(χ, tags="Site,aux-L"))
    push!(sites, Index(χ, tags="Site,aux-R"))
    ψ = randomMPS(Complex{Float64}, sites, ξ)
    ψ[1] = delta(Complex{Float64}, sites[1], linkinds(ψ, 1))
    ψ[system_size + 2] = delta(Complex{Float64}, linkinds(ψ, system_size + 1), sites[system_size + 2])
    println("initial state successfully generated: norm = $(norm(ψ)) vs sqrt(χ) = $(sqrt(χ))")

    return (ψ, sites)
end

function canonicalform!(A::MPS, max)
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
end

function canonicalsummation(init::MPS, h::MPO, χ, l, k_max, temp_delta, num_temps, m)
    println("start calculating canonical summation")

    # temperatures to calculate cTPQ
    temps = (num_temps * temp_delta):-temp_delta:temp_delta
    println("temperature range: $(temp_delta) to $(num_temps * temp_delta)")

    print("step 000")
    ψ = deepcopy(init)
    khk = inner(ψ', h, ψ)
    kh2k = inner(h, ψ, h, ψ)
    kmk = inner(ψ', m, ψ)
    km2k = inner(m, ψ, m, ψ)
    kk = 1

    ϕ = l * ψ - apply(h, ψ)
    print("\tmTPQ state(k = 001) generated")
    canonicalform!(ϕ, χ)

    khk1 = l * khk - kh2k # ⟨k|h|k+1⟩ = ⟨k|h(l - h)|k⟩ = l⟨k|h|k⟩ - ⟨k|h²|k⟩
    kh2k1 = inner(h, ψ, h, ϕ)
    kmk1 = inner(ψ', m, ϕ)
    km2k1 = inner(m, ψ, m, ϕ)
    kk1 = l - khk # ⟨k|k+1⟩ = ⟨k|(l - h)|k⟩ = l⟨k|k⟩ - ⟨k|h|k⟩

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
end