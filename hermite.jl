function newton_girard_power_sums(coeffs; max_k::Int=(length(coeffs)*2) - 1)
    # Normalize polynomial to be monic and get degree
    a0 = coeffs[1]
    a0 == 0 && error("Leading coefficient must be non-zero.")
    normalized_coeffs = coeffs ./ a0
    n = length(normalized_coeffs) - 1  # Degree of polynomial
    
    # Extract elementary symmetric functions σ_k (with sign adjustment)
    σ = [(-1)^k * normalized_coeffs[k+1] for k in 1:n]
    
    # Initialize power sums array to include S₀
    S = zeros(typeof(normalized_coeffs[1]), max_k + 1)
    S[1] = n  # S₀ = number of roots (degree)
    
    # Compute S₁ to S_max_k recursively
    for k in 1:max_k
        if k == 1
            # S₁ = σ₁ directly (no recursion needed)
            S[k+1] = isempty(σ) ? 0.0 : σ[1]
        else
            sum_term = zero(eltype(S))
            for i in 1:min(k-1, n)
                sum_term += (-1)^(i-1) * σ[i] * S[k - i + 1]
            end
            σ_k = k <= n ? σ[k] : zero(eltype(S))
            S[k+1] = sum_term + (-1)^(k-1) * k * σ_k
        end
    end
    
    return S[1:max_k+1]
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


coeffs = [16.0, 0.0, -10.0 , 0.0, 1]
power_sums = newton_girard_power_sums(coeffs)
println("Power Sums: ", power_sums)

h_1 = setHermite_1(length(coeffs), power_sums)
display(h_1)