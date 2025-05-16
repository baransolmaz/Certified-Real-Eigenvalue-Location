using Plots
using LinearAlgebra
gr()

function gershgorin_disks(A::AbstractMatrix{<:Number})
    n, m = size(A)
    disks = []
    for i in 1:n
        center = A[i, i]
        radius1 = sum(abs(A[i, j]) for j in 1:n if j != i)
        radius2 = sum(abs(A[k, i]) for k in 1:n if k != i)
        push!(disks, (center=center, radius=min(radius1,radius2))) # minimumu al
    end
    return disks
end

function in_gershgorin_disk(x::Number, center::Number, radius::Real)
    return abs(complex(x) - complex(center)) <= radius
end

function g(center,radius , x , y)
    return (x-center)^2 +(y^2) - (radius^2) 
end

function draw_gershgorin_disks(disks)
    plt = plot(aspect_ratio=1, legend=false, title="Gershgorin Diskleri")

    for d in disks
        c = d.center
        r = d.radius
        θ = range(0, 2π, length=200)
        x = real(c) .+ r .* cos.(θ)
        y = imag(c) .+ r .* sin.(θ)
        plot!(x, y, label="", fill=(0.2, :blue), linecolor=:blue)
        if r == 0
            scatter!([real(c)], [imag(c)], color=:black, markersize=4)
        else
            # Draw the center of the disk
            scatter!([real(c)], [imag(c)], color=:black, markersize=4)
        end
    end

    savefig(plt, "images/gershgorin.png")
    println("Görsel kaydedildi: images/gershgorin.png")
end

function eigenvalue_norm_based_bounds(A::AbstractMatrix)
    # Norm-based bounds
    bounds = Dict{String,Float64}()
    bounds["Spectral norm (2-norm)"] = opnorm(A, 2)
    bounds["Frobenius norm"] = norm(vec(A))
    bounds["1-norm (max column sum)"] = opnorm(A, 1)
    bounds["∞-norm (max row sum)"] = opnorm(A, Inf)

    return bounds
end

function check_matrix_requirements(A::AbstractMatrix{<:Number})
    n,m =size(A)
    if n != m
        return false
    end
    return true
end

function signature(M::Matrix{<:Number})
    # Ensure the matrix is symmetric
   # @assert M == transpose(M) "Matrix must be symmetric."

    # Compute characteristic polynomial coefficients (descending order: x^k to x^0)
    coeff = charpoly_faddeev_leverrier(M)

    # Function to count sign variations, ignoring zeros
    function count_sign_variations(c)
        filtered = filter(!iszero, c)
        isempty(filtered) && return 0
        count = 0
        for i in 1:length(filtered)-1
            if sign(filtered[i]) != sign(filtered[i+1])
                count += 1
            end
        end
        return count
    end

    # Sign variations for p(x)
    v1 = count_sign_variations(coeff)

    # Coefficients for p(-x)
    n = length(coeff)
    coeff_neg = [coeff[i] * (-1)^(n - i) for i in 1:n]

    # Sign variations for p(-x)
    v2 = count_sign_variations(coeff_neg)

    # Signature: difference in sign variations
    σ = v1 - v2

    return σ
end

function newton_girard_power_sums(coeffs; max_k::Int=(length(coeffs) * 2) - 1)
    # Check leading coefficient and convert coefficients to fractions
    a₀ = coeffs[1]
    a₀ == 0 && error("Leading coefficient must be non-zero.")
    T = Rational{BigInt}  # Use exact fractions (adjust Int to BigInt for large numbers)
    coeffs_frac = T.(coeffs)

    # Degree of polynomial and σ_k calculations
    n = length(coeffs_frac) - 1
    sigma = [(-1)^k * coeffs_frac[k+1] // coeffs_frac[1] for k in 1:n]
    # Initialize power sums (S₀ = degree, stored as a fraction)
    S = Vector{T}(undef, max_k + 1)
    S[1] = T(n)

    # Compute S₁ to S_max_k using Newton-Girard recursion
    for k in 1:max_k
        if k == 1
            S[2] = isempty(sigma) ? zero(T) : sigma[1]
        else
            sum_term = zero(T)
            for i in 1:min(k - 1, n)
                sum_term += (-1)^(i - 1) * sigma[i] * S[k-i+1]
            end
            sigma_k = k <= n ? sigma[k] : zero(T)
            S[k+1] = sum_term + (-1)^(k - 1) * k * sigma_k
        end
    end

    return S
end

function setHermite_1(degree, power_sums)
    mat = Matrix{eltype(power_sums)}(undef, degree, degree)  # Preallocate m×n matrix

    for i in 1:degree
        for j in 1:degree
            mat[i, j] = power_sums[i+j-1]
        end
    end

    return mat
end
function g(M)
    Mr = Rational.(M)  # Convert all entries to Rational
    Ir = Matrix{Rational}(I, size(M, 1), size(M, 2))  # Rational identity matrix
    return (Mr - (1 // 10) * Ir) * (Mr + (1 // 10) * Ir)
end

function getMx(h1, k)
    return inv(h1) * (h1^k)
end

function getHg(h1, mx)
    return h1 * g(mx)
end


function companion_matrix(coeffs::Vector{T}) where {T<:Number}
    n = length(coeffs) - 1  # degree of polynomial
    n < 1 && throw(ArgumentError("Polynomial must have degree at least 1"))

    # Normalize to monic form (leading coefficient = 1)
    normalized_coeffs = coeffs ./ coeffs[1]

    # Create companion matrix
    C = zeros(Rational, n, n)

    # Fill subdiagonal with 1s
    for i in 2:n
        C[i, i-1] = one(T)
    end

    # Fill last column with negated coefficients
    for i in 1:n
        C[i, n] = -normalized_coeffs[n+2-i]
    end

    return C
end

function charpoly_faddeev_leverrier(A::Matrix{<:Rational})
    n = size(A, 1)
    T = eltype(A)
    M = Matrix{T}(I, n, n)  # Initialize with identity matrix
    c = Vector{T}(undef, n) # Coefficients [c_{n-1}, c_{n-2}, ..., c₀]

    for k in 1:n
        AM = A * M           # Matrix multiplication
        trace_AM = tr(AM)    # Trace of the product
        c_k = -trace_AM // k # Compute coefficient
        c[k] = c_k

        if k < n
            M = AM + c_k * I  # Update matrix for next iteration
        end
    end

    # Return monic polynomial coefficients: [1, c_{n-1}, ..., c₀]
    return vcat(one(T), c)
end

# Main fonksiyonu
function main()
    A = Matrix{Rational{Int}}(
            [0 -1 0 ;
            1 0 0;
            0 0 1]
    )

    pa = charpoly_faddeev_leverrier(A)
    println("Pa")
    display(pa)
    power_sums = newton_girard_power_sums(pa)
    h1 = setHermite_1(length(pa) - 1, power_sums)
    println("H1")
    display(h1)
    cp= companion_matrix(pa)
    println("\nCp")
    display(cp)
    signH1 = signature(h1)
    println("Sign H1 : $(signH1)")

    if check_matrix_requirements(A)
        disks = gershgorin_disks(A)
        for d in disks
            println("Center: $(d.center) | Radius: $(d.radius)")
            #gx = in_gershgorin_disk(1.2, d.center, d.radius)
            #println(gx)
            
            if d.radius == 0
                println("EigenValue at $(d.center)") 
            else
                g(x) = (d.radius * I)^2 - (x - d.center * I)^2
        
                hg = h1 * g(cp)
                println("\nHg")
                display(hg)
    
                signHg = signature(hg)#!!
                println("Sign Hg : $(signHg)")
                
                if signH1 != signHg
                    println("Circle contains at least one real eigen")
                else
    
                    println("Circle does not contain any real eigen")
                end
    
    
    
            end
        end
        draw_gershgorin_disks(disks)

        #bounds = eigenvalue_norm_based_bounds(A)

        #for (key, value) in bounds
        #    println(rpad(key, 40), ": ", round(value, digits=4))
        #end






    else
        println("---------- Matris Uygun Değil! ----------")
    end

    



end












#Runner
main()