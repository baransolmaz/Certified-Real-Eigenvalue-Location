using Plots
using LinearAlgebra
import GenericSchur: hessenberg

gr()

# Polynomial utilities --------------------------------------------------------
"""
    la_budde_general(A::Matrix{<:BigFloat})

Compute characteristic polynomial coefficients using La Budde's method.
Returns monic polynomial coefficients: [1, c_{n-1}, ..., c₀]
"""
function la_budde_general(A::Matrix{<:BigFloat})
    n = size(A, 1)

    # Reduce matrix to upper Hessenberg form
    F = hessenberg(A)
    H = F.H

    # Extract subdiagonal elements
    β = [0; diag(H, -1)]

    # Initialize coefficient matrix
    C = zeros(n, n)
    C[1, 1] = -H[1, 1]

    for i in 2:n
        # j = 1 case
        C[i, 1] = C[i-1, 1] - H[i, i]

        # j from 2 to i-1
        for j in 2:i-1
            s = 0.0
            # Sum for m=1 to j-2
            for m in 1:j-2
                prod = 1.0
                for t in 0:m-1
                    prod *= β[i-t]
                end
                s += H[i-m, i] * prod * C[i-m-1, j-m-1]
            end

            # Term for m = j-1
            prod = 1.0
            for t in 0:j-2
                prod *= β[i-t]
            end
            s += H[i-j+1, i] * prod

            C[i, j] = C[i-1, j] - H[i, i] * C[i-1, j-1] - s
        end

        # j = i case
        s = 0.0
        for m in 1:i-2
            prod = 1.0
            for t in 0:m-1
                prod *= β[i-t]
            end
            s += H[i-m, i] * prod * C[i-m-1, i-m-1]
        end

        # Term for m = i-1
        prod = 1.0
        for t in 0:i-2
            prod *= β[i-t]
        end
        s += H[1, i] * prod

        C[i, i] = -H[i, i] * C[i-1, i-1] - s
    end

    return [1; C[end, :]]  # Return coefficients including leading 1
end


"""
    newton_girard_power_sums(coeffs; max_k::Int=(length(coeffs) * 2) - 1)

Compute power sums of polynomial roots using Newton-Girard formulas.
"""
function newton_girard_power_sums(coeffs; max_k::Int=(length(coeffs) * 2) - 1)
    # Check leading coefficient and convert to fractions
    a₀ = coeffs[1]
    a₀ == 0 && error("Leading coefficient must be non-zero.")
    T = BigFloat
    coeffs_frac = T.(coeffs)

    n = length(coeffs_frac) - 1
    sigma = [(-1)^k * coeffs_frac[k+1] / coeffs_frac[1] for k in 1:n]
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
    C = zeros(BigFloat, n, n)

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
function signature(M::Matrix{<:BigFloat})
    coeff = la_budde_general(M)
    #display(coeff)
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
            plot!(x, y, fill=(0.2, :brown), linecolor=:brown)
        else
            plot!(x, y, linecolor=:brown)
        end
        scatter!([real(c)], [imag(c)], color=:red, markersize=3)
    end
    
    filepath !== nothing && savefig(plt, filepath)
    return plt
end

"""
    plot_scanned_disks(disks; title="Scanned Disks")

Plot Gershgorin disks with scanning lines.
"""
function plot_scanned_disks(disks; title="Scanned Disks", filepath=nothing)
    plt = plot(aspect_ratio=1, title=title,xlabel="Re", ylabel="Im", legend=false)

    for d in disks
        c = d.center
        r = d.radius
        cx = real(c)
        cy = imag(c)
        θ = range(0, 2π, length=200)
        x = cx .+ r .* cos.(θ)
        y = cy .+ r .* sin.(θ)
        plot!(x, y, linecolor=:brown)
        scatter!([cx], [cy], color=:red, markersize=3)

        if r > 0
            n_lines = 20
            t_vals = range(-r, r, length=n_lines)
            angle = get(d, :fill, false) ? π / 4 : 3π / 4  # 45° or 135°
            line_color = get(d, :fill, false) ? :green : :red

            cos_a = cos(angle)
            sin_a = sin(angle)

            for t in t_vals
                # Parametric line centered at (cx, cy)
                dx = cos_a * t
                dy = sin_a * t

                # Extend line segment symmetrically within the circle
                length_max = sqrt(r^2 - t^2)
                x1 = cx + dx - length_max * sin_a
                y1 = cy + dy + length_max * cos_a
                x2 = cx + dx + length_max * sin_a
                y2 = cy + dy - length_max * cos_a

                plot!([x1, x2], [y1, y2], color=line_color, linewidth=1.5, alpha=0.6)
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

    for interval in intervals
        a = interval.startP
        b = interval.endP

        # Draw vertical boundary lines for all intervals
        vline!(plt, [a], line=(:solid, 1, :black), alpha=0.8)
        vline!(plt, [b], line=(:solid, 1, :black), alpha=0.8)
        

        # Draw semi-transparent rectangle if eigenvalue is known to exist
        if interval.isExist
            plot!(plt, [a, b], [0, 0], line=(:solid, 2, :black))

            plot!(
                [a, b, b, a, a],                 # x-coordinates of the rectangle
                [ymin, ymin, ymax, ymax, ymin],  # y-coordinates
                fill=(0, :blue),
                alpha=0.25,
                line=:transparent,
                label=""
            )
        end
    end

    title!(plt, title)
    filepath !== nothing && savefig(plt, filepath)
    return plt
end


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
            push!(candidate_points, d.center) # TODO
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
    M = Matrix((BigFloat)[
        1.25  1     0.75  0.5  0.25;
        1     0     0     0     0;
        -1    1     0     0     0;
        0     0     1     3     0;
        0     0     0     0.5  5
    ])

    inputMatrix = M
    # Characteristic polynomial and power sums
    pa = la_budde_general(inputMatrix)
    #display(pa)
    power_sums = newton_girard_power_sums(pa)
    
    # Hermite matrix and signature
    n = length(pa) - 1
    h1 = setHermite_1(n, power_sums)
    #display(h1)
    signH1 = signature(h1)
    #println("Signature H1: $signH1")
    
    # Companion matrix for polynomial
    cp = companion_matrix(pa)
    #display(cp)
    
    # Gershgorin analysis
    disks = gershgorin_disks(inputMatrix)
    
    # Disk analysis
    contained_disks, candidate_points = analyze_disks(disks, h1, signH1, cp)
    scanned_plot = plot_scanned_disks(contained_disks, filepath="images/scanned.png")
    plot_gershgorin_disks(contained_disks, filled=true, filepath="images/all_disks.png")
    
    # Interval analysis
    intervals = analyze_intervals(candidate_points, h1, signH1, cp)
    benchmark_intervals(candidate_points, h1, signH1, cp, 0.0000001)
    plot_intervals(intervals, scanned_plot, filepath="images/intervals.png")

end

"""
    benchmark_intervals(points, h1, signH1, cp; tol=1e-6)

Apply binary interval refinement on each [points[i], points[i+1]].
"""
function benchmark_intervals(points, h1, signH1, cp, tol)
    sort!(points)
    all_intervals = []

    for i in 1:(length(points)-1)
        a, b = points[i], points[i+1]
        refined = benchmark_interval(a, b, h1, signH1, cp, tol)
        append!(all_intervals, refined)
    end
    display(all_intervals)
    return all_intervals
end

"""
    benchmark_interval(a, b, h1, signH1, cp; tol=1e-6)

Divide & conquer search:
- If interval [a,b] may contain eigenvalue, split until width ≤ tol.
- Return only the smallest intervals with isExist=true.
"""
function benchmark_interval(a, b, h1, signH1, cp, tol)
    # If interval is smaller than tolerance → final check
    #print("a:\t$a\nb:\t$b\t")
    if abs(b - a) <= tol
        contains_eigen = benchmark_check_interval(a, b, h1, signH1, cp)
        #println(contains_eigen)
        return contains_eigen ? [(startP=a, endP=b, isExist=true)] : []
    end

    # Check if this interval may contain eigenvalue
    contains_eigen = benchmark_check_interval(a, b, h1, signH1, cp)
    #println(contains_eigen)
    if contains_eigen
        # Divide & conquer: split into halves
        mid = (a + b) / 2
        left = benchmark_interval(a, mid, h1, signH1, cp, tol)
        right = benchmark_interval(mid, b, h1, signH1, cp, tol)
        return vcat(left, right)
    else
        # If no eigenvalue, discard
        return []
    end
end

function benchmark_check_interval(a, b, h1, signH1, cp)
    g(x) = (x - a * I) * (x - b * I)
    hg = h1 * g(cp)
    signHg = signature(hg)
    return (signH1 != signHg)
end

# Run application
t1 = time()
main()
elapsed_time = time() - t1;
println("Elapsed time: ", elapsed_time, " seconds");