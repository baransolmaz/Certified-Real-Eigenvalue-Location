using Plots
using LinearAlgebra
using GenericLinearAlgebra
gr()
"""
    charpoly_faddeev_leverrier(A::AbstractMatrix{T}) where {T}

Compute characteristic polynomial coefficients using Faddeev-Leverrier method.
Returns monic polynomial coefficients: [1, c_{n-1}, ..., c₀]
"""
# 1. Update charpoly_faddeev_leverrier for complex numbers
function charpoly_faddeev_leverrier(A::AbstractMatrix{T}) where {T}
    n = size(A, 1)
    # Initialize matrix M and coefficients vector
    M = zeros(T, n, n)
    coeffs = T[1]  # Start with coefficient for λⁿ (always 1)
    
    for k in 1:n
        # Update M: M = A*M + cₖ*I
        M = A * M + coeffs[k] * I
        
        # Compute trace of A*M efficiently
        tr_val = zero(T)
        for i in 1:n
            for j in 1:n
                tr_val += A[i, j] * M[j, i]  # Equivalent to tr(A*M)
            end
        end

        # Calculate next coefficient: cₖ₊₁ = -tr(A*M)/k
        push!(coeffs, -tr_val // k)
    end
    
    return coeffs  # [a₀, a₁, ..., aₙ] for a₀λⁿ + a₁λⁿ⁻¹ + ... + aₙ]
end

"""
    newton_girard_power_sums(coeffs::Vector{T}, m::Int) where {T}

Compute power sums of polynomial roots using Newton-Girard formulas.
"""
function newton_girard_power_sums(coeffs::Vector{T}, m::Int) where {T}
    # Determine appropriate type for arbitrary-precision arithmetic
    U = if T <: Complex
        Complex{Rational{BigInt}}
    else
        Rational{BigInt}
    end
    
    # Convert coefficients to the appropriate arbitrary-precision type
    coeffs_big = map(coeffs) do c
        if c isa Complex
            complex(Rational{BigInt}(real(c)), Rational{BigInt}(imag(c)))
        else
            Rational{BigInt}(c)
        end
    end
    
    n_degree = length(coeffs_big) - 1
    s = Vector{U}(undef, m + 1)
    s[1] = n_degree  # s₀ = number of roots
    
    n_degree == 0 && return s
    
    c = coeffs_big[2:end]  # Non-leading coefficients
    
    for k in 1:m
        jmax = k ≤ n_degree ? k - 1 : n_degree
        term = zero(U)
        
        # Sum: cⱼ * sₖ₋ⱼ for j=1 to jmax
        for j in 1:jmax
            term += c[j] * s[k - j + 1]
        end
        
        # Extra term for k ≤ degree
        if k ≤ n_degree
            term += k * c[k]
        end
        
        s[k + 1] = -term
    end
    
    return s  # [s₀, s₁, ..., sₘ]
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
# 2. Update companion_matrix for complex numbers
function companion_matrix(coeffs::Vector{T}) where {T}
    n = length(coeffs) - 1
    n < 1 && throw(ArgumentError("Polynomial must have degree at least 1"))
    C = zeros(T, n, n)  

    for i in 2:n
        C[i, i-1] = one(eltype(C))
    end

    for i in 1:n
        C[i, n] = -coeffs[n+2-i] // coeffs[1]
    end

    return C
end

# Signature and matrix analysis -----------------------------------------------
"""
    signature(M::Matrix{<:Number})

Compute the signature of a symmetric matrix using sign variations.
"""
function signature(M::AbstractMatrix{T}) where {T}
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
# 3. Update gershgorin_disks for complex centers
function gershgorin_disks(A::AbstractMatrix{T}) where {T}
    n = size(A, 1)
    disks = []
    for i in 1:n
        center = A[i, i]
        sums = sum(abs, A[i, j] for j in 1:n if j != i)
        radius = ceil(sums, digits=6)
        #display(sums)
        #display(radius)
        #display(rationalize(radius))
        push!(disks, (center=center, radius=rationalize(radius)))
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
        c = real(d.center)
        r = real(d.radius)
        θ = range(0, 2π, length=200)
        x = real(c) .+ r .* cos.(θ)
        y = imag(c) .+ r .* sin.(θ)

        if filled
            plot!(x, y, fill=(0.2, :blue), linecolor=:blue, alpha=0.5)
        else
            plot!(x, y, linecolor=:blue)
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
function plot_scanned_disks(disks; title="Gershgorin Disks", filepath=nothing)
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
function plot_intervals(intervals, plt; title="", filepath=nothing)
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

    title!(plt, "Eigenvalue Intervals of $title")
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
        radius_real = real(d.radius)
        #display(d.center)
        #display(d.radius)
        if radius_real == 0
            push!(candidate_points, d.center)
            push!(contained_disks, (center=d.center, radius=d.radius, fill=true))
        else
            push!(candidate_points, d.center)
            push!(candidate_points, d.center - d.radius)
            push!(candidate_points, d.center + d.radius)
            push!(contained_disks, (center=d.center, radius=d.radius, fill=true))
           #if signH1 != signHg
            #    push!(contained_disks, (center=d.center, radius=d.radius, fill=true))
           # else
            #    push!(contained_disks, (center=d.center, radius=d.radius, fill=false))
            #end
        end
    end
    #display(candidate_points)
    unique!(candidate_points)
    sort!(candidate_points; by= x -> (real(x), imag(x)))

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
        #display("$a $b")
        g(x) = (x - a * I) * (x - b * I)
        #display(g(cp))
        hg = h1 * g(cp)
       # display(hg)
        signHg = signature(hg)
        #display(signHg)
        
        contains_eigen = (signH1 != signHg)
        push!(intervals, (startP=a, endP=b, isExist=contains_eigen))
        
        println("Interval [$a, $b]: $(contains_eigen ? "Contains" : "No") real eigenvalue")
    end
    
    return intervals
end

function compute_s(A::AbstractMatrix)
    n, m = size(A)
    return sqrt(((tr(A * A) - ((tr(A)^2)) / n) / n))
end

function mean_m(A::AbstractMatrix)
    n, m = size(A)
    return tr(A) / n
end

function all_eigenvalue_bounds(A::AbstractMatrix)
    n = size(A, 1)
    m = mean_m(A)
    s = compute_s(A)
    bounds = Vector{NamedTuple}(undef, n)

    for k in 1:n
        if k == 1
            #1.10 2.3
            # Bounds for largest eigenvalue (λ₁)
            lower = m + s / sqrt(n - 1) # min_bound
            upper = m + s * sqrt(n - 1) # max_bound
        elseif k == n #2.2
            # Bounds for smallest eigenvalue (λₙ)
            lower = m - s * sqrt(n - 1)
            upper = m - s / sqrt(n - 1)
        else
            # Bounds for intermediate eigenvalues (Theorem 2.2)
            lower = m - s * sqrt((k - 1) / (n - k + 1))
            upper = m + s * sqrt((n - k) / k)
        end
        bounds[k] = (lambda="lambda_$k", lower=real(lower), upper=real(upper))
    end

    return bounds
end

# Main application -----------------------------------------------------------
function main()
    # Test with complex matrix
    M_complex = Matrix([
        7+3im  -4-6im   -4;
        -1-6im   7    -2-6im;
          2     4-6im   13-3im
    ])
    #check_hermittian(M_complex, "ComplexMatrix")


    # Original real matrix
    M_real = Matrix{Rational{BigInt}}([
        5//4 1 3//4 1//2 1//4;
        1 0 0 0 0;
        -1 1 0 0 0;
        0 0 1 3 0;
        0 0 0 1//2 5
    ])
    check_hermittian(M_real, "RealMatrix")

    M_complex2 = Matrix([
        8 1 1 1 1 1;
        0 3-3im 1 1 1 1;
        5 0 3+3im 1 1 1;
        3 2 0 5 1 1;
        0 0 0 0 3-3im 1;
        0 1im 1 0 0 3+3im
        ])
    #check_hermittian(M_complex2, "ComplexMatrix_2") #Sadece Complex

    M_real2 = Matrix{Rational{BigInt}}([
        3 1 1 1 1 1;
        0 5 1 1 1 1;
        8 0 3 1 1 1;
        4 0 0 2 1 1;
        0 0 1.5 0 2.5 1;
        3.75 0 0 0 0 5
    ])
    #check_hermittian(M_real2, "RealMatrix_2") #Sadece Real
    M_real3 = Matrix([
        1 0 2;
        -1 1 3;
        0 0 2
    ])
    #check_hermittian(M_real3, "RealMatrix_3")
    BOTH = Matrix([
        3+4im 1 1 1 1 1;
        0 6 1 1 1 1;
        -2im 0 3-4im 1 1 1;
        2 0 0 3-4im 1 1;
        0 1 0 0 6 1;
        0 2im 0 0 0 3+4im
    ])
    #check_hermittian(BOTH, "BOTH") #Hem  Real hem de complex
end

function get_candidate_points(disks)
    candidate_points = []
    for d in disks
        center_real = real(d.center)
        radius_real = real(d.radius)
        push!(candidate_points, center_real)
        push!(candidate_points, center_real - radius_real)
        push!(candidate_points, center_real + radius_real)
    end
    unique!(candidate_points)
    sort!(candidate_points)
    return candidate_points
end
function check_eigenvalue_locations(inputMatrix::AbstractMatrix{T}, name::AbstractString) where {T}
    ensure_dir("images/$name")
    display(name)
    # Characteristic polynomial and power sums
    #display(inputMatrix)
    pa = charpoly_faddeev_leverrier(inputMatrix)
    #display(pa)
    power_sum = newton_girard_power_sums(pa, (length(pa) * 2) - 1)
    #display(power_sum)


    # Hermite matrix and signature
    n = length(pa) - 1
    h1 = setHermite_1(n, power_sum)
    #display(h1)
    signH1 = signature(h1)
    #println("Signature H1: $signH1")

    # Companion matrix for polynomial
    cp = companion_matrix(pa)
    #display(cp)

    # Gershgorin analysis
    row_disks = gershgorin_disks(inputMatrix)
    #display(row_disks)
    plot_gershgorin_disks(row_disks, filepath="images/$name/all_disks.png")

    # Disk analysis
    contained_disks, candidate_points = analyze_disks(row_disks, h1, signH1, cp)
    #display(candidate_points)
    scanned_plot = plot_scanned_disks(contained_disks, filepath="images/$name/scanned.png")
    candidate_points = get_candidate_points(row_disks)

    plot_gershgorin_disks(contained_disks, filled=true, filepath="images/$name/remain_disks.png")

    # Interval analysis
    intervals = analyze_intervals(candidate_points, h1, signH1, cp)
    result_plot = plot_intervals(intervals, scanned_plot,title=name, filepath="images/$name/intervals.png")

    bounds = all_eigenvalue_bounds(inputMatrix)
    draw_intervals(result_plot, bounds, filepath="images/$name/intervals_with_bounds.png")
    draw_only_intervals(bounds, filepath="images/$name/bounds.png")

    for b in bounds
      # println("$(b.lambda)    $(b.lower)  $(b.upper)")
    end


    #=for interval in intervals
        if interval.isExist
            #println("Interval [$(interval.startP), $(interval.endP)]: ")
            for b in bounds
                i = interval
                case1 = b.lower <= i.startP && i.endP <= b.upper
                case2 = i.startP < b.upper && b.upper < i.endP
                case3 = b.lower > i.startP && i.endP > b.lower
                case4 = b.lower >= i.startP && i.endP >= b.upper

                if case1 || case2 || case3 || case4
                #    println("\t $(b.lambda)")
                end
            end
        end

    end=#
    return row_disks, filter!(x -> x.isExist, intervals)
end

function draw_intervals(plt, V::Vector{NamedTuple}; filepath="images/bounds.png")
    # Get current plot boundaries
    ymin, ymax = ylims()
    xmin, xmax = xlims()

    # Calculate offset for bounds (place at bottom of plot)
    bounds_y = ymin + 0.13 * (ymax - ymin)
    bounds_height = 0.03 * (ymax - ymin)

    colors = palette(:tab10)

    for (i, b) in enumerate(V)
        color = colors[(i-1)%length(colors)+1]

        # Draw horizontal line for bound
        plot!(plt, [b.lower, b.upper], [bounds_y, bounds_y],
            linewidth=2,
            color=color,
            label=b.lambda)

        # Draw vertical markers at boundaries
        vline!(plt, [b.lower], line=(:dash, 1, color), label="")
        vline!(plt, [b.upper], line=(:dash, 1, color), label="")

        # Add text label
        mid = (b.lower + b.upper) / 2
        annotate!(plt, mid, bounds_y - bounds_height, text(b.lambda, 8, :center, color))

        # Move down for next bound
        bounds_y -= bounds_height * 2
    end
    
    # Reset y-axis limits to include bounds
    new_ymin = bounds_y - bounds_height
    ylims!(new_ymin, ymax)

    if filepath !== nothing
        savefig(plt, filepath)
        println("Bounds plot saved to $filepath")
    end

    return plt
end

function draw_only_intervals(V::Vector{NamedTuple}; filepath="images/bounds.png")
    colors = palette(:tab10)
    plt = plot(title="Eigenvalue Bounds", xlabel="Value", ylabel="Eigenvalue", legend=:right)

    for (i, b) in enumerate(V)
        color = colors[(i-1)%length(colors)+1]
        y = length(V) - i + 1
        plot!(plt, [b.lower, b.upper], [y,y],
            linewidth=2,
            color=color,
            label=b.lambda)

        # Draw vertical markers at boundaries
        vline!(plt, [b.lower], line=(:dash, 1, color), label="")
        vline!(plt, [b.upper], line=(:dash, 1, color), label="")
    end

    yticks!(reverse(1:length(V)), [lambda for (lambda, _, _) in V])
    if filepath !== nothing
        savefig(plt, filepath)
        println("Plot saved to $filepath")
    end
end

function plot_intersections(disk_A, interval_B, interval_C,eigenvalues, name; filepath="images/intersections_$name.png")
    # Determine plot limits from disks
    all_x = Float64[]
    all_y = Float64[]
    for d in disk_A
        cx = real(d.center)
        cy = imag(d.center)
        r = d.radius
        push!(all_x, cx - r, cx + r)
        push!(all_y, cy - r, cy + r)
    end
    xmin, xmax = extrema(all_x)
    ymin, ymax = extrema(all_y)
    padding = 0.1
    x_pad = padding * (xmax - xmin)
    y_pad = padding * (ymax - ymin)
    xlims = (xmin - x_pad, xmax + x_pad)
    ylims = (ymin - y_pad, ymax + y_pad)

    plt = plot(xlim=xlims, ylim=ylims, aspect_ratio=1,
        title="Eigenvalue Locations: $name",
        xlabel="Re", ylabel="Im",legend=false)

    # Extract constraint regions
    real_constraints = [intv for intv in interval_B if intv.isExist]
    imag_constraints = [intv for intv in interval_C if intv.isExist]

    # Plot real constraints as YELLOW vertical bands
    for intv in real_constraints
        a, b = intv.startP, intv.endP
        plot!(plt, [a, a, b, b, a], [ylims[1], ylims[2], ylims[2], ylims[1], ylims[1]],
            seriestype=:shape, fill=true, color=:yellow, alpha=0.15, label="Real Constraint")
    end

    # Plot imaginary constraints as BLUE horizontal bands
    for intv in imag_constraints
        c, d = intv.startP, intv.endP
        plot!(plt, [xlims[1], xlims[1], xlims[2], xlims[2], xlims[1]], [c, d, d, c, c],
            seriestype=:shape, fill=true, color=:blue, alpha=0.15, label="Imag Constraint")
    end

    # Plot Gershgorin disks (filled with light purple)
    for d in disk_A
        c = d.center
        r = d.radius
        cx = real(c)
        cy = imag(c)
        θ = range(0, 2π, length=200)
        x = cx .+ r .* cos.(θ)
        y = cy .+ r .* sin.(θ)
        plot!(plt, x, y, seriestype=:shape, fill=true, color=:purple, alpha=0.2, label="Gershgorin Disk")
    end

    # Highlight intersection areas with contour
    resolution = 200
    x_grid = range(xlims[1], xlims[2], length=resolution)
    y_grid = range(ylims[1], ylims[2], length=resolution)

    # Create mask for real constraints (YELLOW)
    real_mask = zeros(length(x_grid), length(y_grid))
    for intv in real_constraints
        a, b = intv.startP, intv.endP
        for (i, x) in enumerate(x_grid)
            if a ≤ x ≤ b
                real_mask[i, :] .= 1
            end
        end
    end

    # Create mask for imaginary constraints (BLUE)
    imag_mask = zeros(length(x_grid), length(y_grid))
    for intv in imag_constraints
        c, d = intv.startP, intv.endP
        for (j, y) in enumerate(y_grid)
            if c ≤ y ≤ d
                imag_mask[:, j] .= 1
            end
        end
    end

    # Create mask for disks (PURPLE)
    disk_mask = zeros(length(x_grid), length(y_grid))
    for d in disk_A
        cx = real(d.center)
        cy = imag(d.center)
        r = d.radius
        for (i, x) in enumerate(x_grid)
            for (j, y) in enumerate(y_grid)
                if (x - cx)^2 + (y - cy)^2 ≤ r^2
                    disk_mask[i, j] = 1
                end
            end
        end
    end

    # Combined mask (intersection)
    combined_mask = real_mask .* imag_mask .* disk_mask

    # Plot intersection contour in RED
    contour!(plt, x_grid, y_grid, combined_mask',
        levels=[0.5], color=:red, linewidth=3,
        label="Intersection Region")

    # Add center points
    for d in disk_A
        cx = real(d.center)
        cy = imag(d.center)
        scatter!([cx], [cy], color=:white, markersize=2, label="")
    end
    for d in eigenvalues
        cx = real(d)
        cy = imag(d)
        scatter!([cx], [cy], color=:black, markersize=3, label="")
    end

    savefig(plt, filepath)
    println("Saved intersection plot to $filepath")
    return plt
end

function ensure_dir(path::AbstractString)
    if !isdir(path)
        mkpath(path)
        println("Created directory: $path")
    end
end

function check_hermittian(matrix::AbstractMatrix{T}, name) where {T}
    A = matrix
    B = (A + A') // 2
    #display(B)
    C = (A - A') // (2im)
    #display(C)

    # Analyze original complex matrix
    disk_A, interval_A = check_eigenvalue_locations(A, "$name-A")

    # Analyze symmetric part (real eigenvalues)
    _, interval_B = check_eigenvalue_locations(B, "$name-B")

   # display(interval_intersection(interval_A,interval_B))
    # Analyze skew-symmetric part (imaginary eigenvalues)
    _, interval_C = check_eigenvalue_locations(C, "$name-C")

    # Generate combined plot
    eigenvalues = eigvals(A)
    display("Eigen Values: ")
    display(eigenvalues)
    #plot_intersections(disk_A, interval_B, interval_C,eigenvalues, name)
    #real_b , imag_b = eigen_bounds(A)

    for realP in interval_B
        for imgP in interval_C 
            min = realP.startP + (imgP.startP*1im)
            max = realP.endP + (imgP.endP * 1im)
            println("$(min)  $(max)")
            
        end 
        
    end

    #draw_bounds(interval_B,interval_C,
    #    name="SampleMatrix",
    #    filepath="custom_intersection.png")


end

function interval_intersection(v1, v2)
    # Filter out intervals where isExist is false
    filtered_v1 = filter(x -> x.isExist, v1)
    filtered_v2 = filter(x -> x.isExist, v2)

    # Sort intervals by startP
    sort!(filtered_v1, by=x -> x.startP)
    sort!(filtered_v2, by=x -> x.startP)

    # Initialize pointers and result vector
    i, j = 1, 1
    result = []
    n1, n2 = length(filtered_v1), length(filtered_v2)

    # Traverse both vectors
    while i ≤ n1 && j ≤ n2
        a = filtered_v1[i]
        b = filtered_v2[j]

        # Compute potential overlap
        start_ov = max(a.startP, b.startP)
        end_ov = min(a.endP, b.endP)

        # If overlap exists, add to results
        if start_ov <= end_ov
            push!(result, (startP=start_ov, endP=end_ov, isExist=true))
        end

        # Move pointer based on end points
        if a.endP < b.endP
            i += 1
        elseif b.endP < a.endP
            j += 1
        else
            i += 1
            j += 1
        end
    end

    return result
end

function draw_bounds(realB, imagB; name="Matrix", filepath="images/bounds_intersection_$name.png")
    # Determine plot limits from both constraint sets
    all_x = [lo for (lo, hi) in realB]
    append!(all_x, [hi for (lo, hi) in realB])
    all_y = [lo for (lo, hi) in imagB]
    append!(all_y, [hi for (lo, hi) in imagB])

    xmin, xmax = extrema(all_x)
    ymin, ymax = extrema(all_y)
    padding = 0.1
    x_pad = padding * (xmax - xmin)
    y_pad = padding * (ymax - ymin)
    xlims = (xmin - x_pad, xmax + x_pad)
    ylims = (ymin - y_pad, ymax + y_pad)

    plt = plot(xlim=xlims, ylim=ylims, aspect_ratio=1,
        title="Eigenvalue Constraints: $name",
        xlabel="Re", ylabel="Im", legend=:topright)

    # Plot real constraints as YELLOW vertical bands
    for (lo, hi) in realB
        plot!(plt, [lo, lo, hi, hi, lo], [ylims[1], ylims[2], ylims[2], ylims[1], ylims[1]],
            seriestype=:shape, fill=true, color=:yellow, alpha=0.2,
            label=(lo == realB[1][1] ? "Real Constraint" : ""))
    end

    # Plot imaginary constraints as BLUE horizontal bands
    for (lo, hi) in imagB
        plot!(plt, [xlims[1], xlims[1], xlims[2], xlims[2], xlims[1]], [lo, hi, hi, lo, lo],
            seriestype=:shape, fill=true, color=:blue, alpha=0.2,
            label=(lo == imagB[1][1] ? "Imag Constraint" : ""))
    end

    # Highlight intersection areas
    resolution = 200
    x_grid = range(xlims[1], xlims[2], length=resolution)
    y_grid = range(ylims[1], ylims[2], length=resolution)

    # Create mask for real constraints
    real_mask = zeros(length(x_grid), length(y_grid))
    for (lo, hi) in realB
        for (i, x) in enumerate(x_grid)
            if lo ≤ x ≤ hi
                real_mask[i, :] .= 1
            end
        end
    end

    # Create mask for imaginary constraints
    imag_mask = zeros(length(x_grid), length(y_grid))
    for (lo, hi) in imagB
        for (j, y) in enumerate(y_grid)
            if lo ≤ y ≤ hi
                imag_mask[:, j] .= 1
            end
        end
    end

    # Create combined mask for intersections
    combined_mask = real_mask .* imag_mask

    # Plot intersection regions in RED
    contourf!(plt, x_grid, y_grid, combined_mask',
        levels=[0.5, 1.0],
        color=:red, alpha=0.4,
        label="Intersection Region")

    # Add boundary lines for clarity
    for (lo, hi) in realB
        vline!([lo, hi], line=(:dash, 1, :yellow), label="")
    end

    for (lo, hi) in imagB
        hline!([lo, hi], line=(:dash, 1, :blue), label="")
    end

    # Save and show result
    savefig(plt, filepath)
    println("Saved bounds intersection plot to $filepath")
    return plt
end
function eigen_bounds(A::AbstractMatrix{T}) where {T<:Number}
    n,i = size(A)
    B = (A + A') / 2
    C = (A - A') / (2im)

    # Real part bounds
    trB = real(tr(B))
    trB² = real(tr(B^2))
    m_b = trB / n
    s_b = sqrt((trB² - trB^2 / n) / n)

    # Imaginary part bounds
    trC = real(tr(C))
    trC² = real(tr(C^2))
    m_c = trC / n
    s_c = sqrt((trC² - trC^2 / n) / n)

    # Modulus bounds
    trAA = real(tr(A' * A))
    m_a = abs(tr(A) / n)
    s_a = sqrt((trAA / n) - m_a^2)

    # Precompute coefficients
    k_vals = 1:n
    sqrt_terms_real = [sqrt((k - 1) / (n - k + 1)) for k in k_vals]
    sqrt_terms_imag = [sqrt((n - k) / k) for k in k_vals]

    # Generate bounds
    real_bounds = [(m_b - s_b * sqrt_terms_real[k], m_b + s_b * sqrt_terms_imag[k]) for k in k_vals]
    imag_bounds = [(m_c - s_c * sqrt_terms_real[k], m_c + s_c * sqrt_terms_imag[k]) for k in k_vals]
    mod_bounds = [(max(m_a - s_a * sqrt(n - 1), 0), sqrt(trAA / n) + s_a * sqrt((n - 1) / n)) for _ in k_vals]

    # Sort bounds descendingly
    sort!(real_bounds, rev=true)
    sort!(imag_bounds, rev=true)
    sort!(mod_bounds, rev=true)

    return real_bounds, imag_bounds#, mod_bounds
end
# Run application
main()