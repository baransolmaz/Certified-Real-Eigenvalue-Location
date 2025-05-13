using LinearAlgebra
using Polynomials
using AbstractAlgebra

function newton_girard_power_sums(coeffs; max_k::Int=(length(coeffs)*2) - 1)
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
            for i in 1:min(k-1, n)
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

function getMx(h1,k)
    return inv(h1)* (h1^k)
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

function signature(M::Matrix{<:Number})
    # Ensure the matrix is symmetric
    @assert M == transpose(M) "Matrix must be symmetric."

    # Compute characteristic polynomial coefficients (descending order: x^k to x^0)
    coeff = charpoly_faddeev_leverrier(M)
    display(coeff)

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



coeffs = [16.0, 0.0, -10.0 , 0.0, 1]
power_sums = newton_girard_power_sums(coeffs)
println("Power Sums: ", power_sums)

h_1 = setHermite_1(length(coeffs)-1, power_sums)
display(h_1)

h_1_bigfloat = AbstractMatrix{Float64}(h_1)
rank_h1 = rank(h_1_bigfloat)
println("Rank of h_1: ", rank_h1)

#m_x = getMx(h_1, rank_h1)
m_x = companion_matrix(coeffs)
display(m_x)
hg = getHg(h_1,m_x)
display(hg)

sig = signature(hg)
display(sig)
sig1 = signature(h_1)
display(sig1)
