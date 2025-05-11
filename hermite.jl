using LinearAlgebra
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
    return (M-(5*I))*(M-(2*I))
end

function getMx(h1,k)
    return inv(h1)* (h1^k)
end

function getHg(h1, mx)
    return h1 * g(mx)
end

coeffs = [16.0, 0.0, -10.0 , 0.0, 1]
power_sums = newton_girard_power_sums(coeffs)
println("Power Sums: ", power_sums)

h_1 = setHermite_1(length(coeffs), power_sums)
display(h_1)

h_1_bigfloat = AbstractMatrix{Float64}(h_1)
rank_h1 = rank(h_1_bigfloat)
println("Rank of h_1: ", rank_h1)

m_x = getMx(h_1, rank_h1)
display(m_x)
hg = getHg(h_1,m_x)
display(hg)


