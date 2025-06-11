using Plots
using LinearAlgebra
gr()

# Polynomial utilities --------------------------------------------------------
"""
    charpoly_faddeev_leverrier(A::Matrix{<:Rational})

Compute characteristic polynomial coefficients using Faddeev-Leverrier method.
Returns monic polynomial coefficients: [1, c_{n-1}, ..., c₀]
"""
function charpoly_faddeev_leverrier(A::Matrix{<:Rational})
    n = size(A, 1)
    T = eltype(A)
    M = Matrix{T}(I, n, n)  # Identity matrix
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

    return vcat(one(T), c)  # [1, c_{n-1}, ..., c₀]
end

"""
    newton_girard_power_sums(coeffs; max_k::Int=(length(coeffs) * 2) - 1)

Compute power sums of polynomial roots using Newton-Girard formulas.
"""
function newton_girard_power_sums(coeffs; max_k::Int=(length(coeffs) * 2) - 1)
    # Check leading coefficient and convert to fractions
    a₀ = coeffs[1]
    a₀ == 0 && error("Leading coefficient must be non-zero.")
    T = Rational{BigInt}
    coeffs_frac = T.(coeffs)

    n = length(coeffs_frac) - 1
    sigma = [(-1)^k * coeffs_frac[k+1] // coeffs_frac[1] for k in 1:n]
    S = Vector{T}(undef, max_k + 1)
    S[1] = T(n)  # S₀ = degree

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

"""
    setHermite_1(degree, power_sums)

Construct the Hermite matrix from power sums.
"""
function setHermite_1(degree, power_sums)
    mat = Matrix{eltype(power_sums)}(undef, degree, degree)
    for i in 1:degree
        for j in 1:degree
            mat[i, j] = power_sums[i+j-1]
        end
    end
    return mat
end

"""
    companion_matrix(coeffs::Vector{T}) where {T<:Number}

Create companion matrix for a monic polynomial with given coefficients.
"""
function companion_matrix(coeffs::Vector{T}) where {T<:Number}
    n = length(coeffs) - 1
    n < 1 && throw(ArgumentError("Polynomial must have degree at least 1"))

    # Normalize to monic form
    normalized_coeffs = coeffs ./ coeffs[1]
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

# Signature and matrix analysis -----------------------------------------------
"""
    signature(M::Matrix{<:Number})

Compute the signature of a symmetric matrix using sign variations.
"""
function signature(M::Matrix{<:Number})
    coeff = charpoly_faddeev_leverrier(M)
    
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
    v2 = count_sign_variations(coeff_neg)

    return v1 - v2  # Signature
end

"""
    gershgorin_disks(A::AbstractMatrix{<:Number})

Compute Gershgorin disks for a square matrix.
"""
function gershgorin_disks(A::AbstractMatrix{<:Number})
    n, m = size(A)
    @assert n == m "Matrix must be square"
    disks = []
    for i in 1:n
        center = A[i, i]
        radius = sum(abs(A[i, j]) for j in 1:n if j != i)
        push!(disks, (center=center, radius=radius)) 
    end
    return disks
end

# Visualization functions -----------------------------------------------------
"""
    plot_gershgorin_disks(disks; title="Gershgorin Disks")

Plot Gershgorin disks with optional filling.
"""
function plot_gershgorin_disks(disks; title="Gershgorin Disks", filled=false, filepath=nothing)
    plt = plot(aspect_ratio=1, title=title, xlabel="Re", ylabel="Im", legend=false)
    
    for d in disks
        c = d.center
        r = d.radius
        θ = range(0, 2π, length=200)
        x = real(c) .+ r .* cos.(θ)
        y = imag(c) .+ r .* sin.(θ)
        
        if filled && get(d, :fill, false)
            plot!(x, y, fill=(0.2, :blue), linecolor=:blue)
        else
            plot!(x, y, linecolor=:blue)
        end
        scatter!([real(c)], [imag(c)], color=:red, markersize=4)
    end
    
    filepath !== nothing && savefig(plt, filepath)
    return plt
end

"""
    plot_scanned_disks(disks; title="Scanned Disks")

Plot Gershgorin disks with scanning lines.
"""
function plot_scanned_disks(disks; title="Scanned Disks", filepath=nothing)
    plt = plot(aspect_ratio=1, title=title, legend=false)
    
    for d in disks
        c = d.center
        r = d.radius
        cx = real(c)
        cy = imag(c)
        θ = range(0, 2π, length=200)
        x = cx .+ r .* cos.(θ)
        y = cy .+ r .* sin.(θ)
        plot!(x, y, linecolor=:blue)
        scatter!([cx], [cy], color=:red, markersize=4)

        if r > 0
            n_lines = 20
            a = cx - r
            b = cx + r
            x_vals = range(a, b, length=n_lines)
            line_color = get(d, :fill, false) ? :green : :red
            
            if !get(d, :fill, false)
                x_vals = reverse(x_vals)
            end

            for x_val in x_vals
                dx = x_val - cx
                offset = sqrt(r^2 - dx^2)
                y1 = cy - offset
                y2 = cy + offset
                plot!([x_val, x_val], [y1, y2], 
                      color=line_color, linewidth=1.5, alpha=0.6)
            end
        end
    end
    
    filepath !== nothing && savefig(plt, filepath)
    return plt
end

"""
    plot_intervals(intervals, plt; title="Eigenvalue Intervals")

Plot eigenvalue intervals on an existing plot.
"""
function plot_intervals(intervals, plt; title="Eigenvalue Intervals", filepath=nothing)
    # Get current plot boundaries
    ymin, ymax = ylims()
    xmin, xmax = xlims()
    plot_height = ymax - ymin

    for interval in intervals
        a = interval.startP
        b = interval.endP

        # Draw vertical boundary lines for all intervals
        vline!(plt,[a], line=(:solid, 2, :black), alpha=0.8)
        vline!([b], line=(:solid, 2, :black), alpha=0.8)

        # Add diagonal scans for eigenvalue-containing intervals
        if interval.isExist
            width = b - a

            # Calculate adaptive line density (more lines for narrower intervals)
            base_density = 40
            min_density = 30
            max_density = 150
            density_factor = min(max_density, max(min_density, base_density * ((xmax - xmin) / width)))
            n_lines = round(Int, density_factor)

            # Calculate step size for diagonal lines
            step = width / n_lines

            # Draw dotted diagonal lines at 45 degrees
            for i in 0:n_lines
                x_start = a + i * step

                # Line starting from bottom boundary
                y_end = min(ymax, ymin + (b - x_start))
                x_end = x_start + (y_end - ymin)
                plot!([x_start, x_end], [ymin, y_end],
                    line=(:dot, 1.2, :brown),
                    alpha=0.7, label="")

                # Line starting from left boundary (for tall plots)
                if i > 0
                    y_start = ymin + i * step
                    if y_start < ymax
                        x_end2 = min(b, a + (ymax - y_start))
                        y_end2 = y_start + (x_end2 - a)
                        plot!([a, x_end2], [y_start, y_end2],
                            line=(:dot, 1.2, :brown),
                            alpha=0.7, label="")
                    end
                end
            end
        end
    end

    filepath !== nothing && savefig(plt, filepath)
    return plt
end

# Main analysis functions ----------------------------------------------------
"""
    analyze_disks(disks, h1, signH1, cp)

Analyze Gershgorin disks for real eigenvalue containment.
"""
function analyze_disks(disks, h1, signH1, cp)
    contained_disks = []
    candidate_points = []
    
    for d in disks
        if d.radius == 0
            # Exact eigenvalue case
            push!(candidate_points, d.center)
            push!(contained_disks, (center=d.center, radius=d.radius, fill=true))
            println("Eigenvalue at $(d.center)")
        else
            # Test using quadratic form
            g(x) = (x - d.center * I)^2 - (d.radius * I)^2
            hg = h1 * g(cp)
            signHg = signature(hg)
            
            push!(candidate_points, d.center - d.radius)
            push!(candidate_points, d.center + d.radius)
            if signH1 != signHg
                push!(contained_disks, (center=d.center, radius=d.radius, fill=true))
                println("Disk contains real eigenvalue: center=$(d.center), radius=$(d.radius)")
            else
                push!(contained_disks, (center=d.center, radius=d.radius, fill=false))
                println("Disk contains no real eigenvalue: center=$(d.center), radius=$(d.radius)")
            end
        end
    end
    
    # Sort and deduplicate points
    unique!(candidate_points)
    sort!(candidate_points)
    return contained_disks, candidate_points
end

"""
    analyze_intervals(points, h1, signH1, cp)

Analyze intervals between critical points for eigenvalue containment.
"""
function analyze_intervals(points, h1, signH1, cp)
    intervals = []
    sort!(points)
    
    for i in 1:(length(points)-1)
        a = points[i]
        b = points[i+1]
        g(x) = (x - a * I) * (x - b * I)
        hg = h1 * g(cp)
        signHg = signature(hg)
        
        contains_eigen = (signH1 != signHg)
        push!(intervals, (startP=a, endP=b, isExist=contains_eigen))
        
        println("Interval [$a, $b]: $(contains_eigen ? "Contains" : "No") real eigenvalue")
    end
    
    return intervals
end

# Main application -----------------------------------------------------------
function main()
    # Define input matrix
    B3 = Matrix{Rational{Int}}([
        5//4  1     3//4  1//2  1//4;
        1     0     0     0     0;
        -1    1     0     0     0;
        0     0     1     3     0;
        0     0     0     1//2  5
    ])
    
    # Characteristic polynomial and power sums
    pa = charpoly_faddeev_leverrier(B3)
    power_sums = newton_girard_power_sums(pa)
    
    # Hermite matrix and signature
    n = length(pa) - 1
    h1 = setHermite_1(n, power_sums)
    signH1 = signature(h1)
    println("Signature H1: $signH1")
    
    # Companion matrix for polynomial
    cp = companion_matrix(pa)
    
    # Gershgorin analysis
    disks = gershgorin_disks(B3)
    plot_gershgorin_disks(disks, filepath="images/all_disks.png")
    
    # Disk analysis
    contained_disks, candidate_points = analyze_disks(disks, h1, signH1, cp)
    scanned_plot = plot_scanned_disks(contained_disks, filepath="images/scanned.png")
    plot_gershgorin_disks(contained_disks, filled=true, filepath="images/remain_disks.png")
    
    # Interval analysis
    intervals = analyze_intervals(candidate_points, h1, signH1, cp)
    result_plot = plot_intervals(intervals, scanned_plot, filepath="images/intervals.png")
end

# Run application
main()