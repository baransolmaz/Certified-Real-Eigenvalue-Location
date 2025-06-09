using LinearAlgebra
using Printf
using Plots
using ColorSchemes

function compute_s(A::AbstractMatrix)
    n, m = size(A)
    x = ((tr(A * A) - ((tr(A)^2)) / n) / n)
    if x <  0.0
        return 0
    end
    return sqrt(x)
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
        bounds[k] = (lambda="lambda_$k", lower=lower, upper=upper)
    end

    return bounds
end

function gershgorin_disks(A::AbstractMatrix{<:Number})
    n = size(A, 1)
    disks = []
    for i in 1:n
        center = A[i, i]
        radius1 = sum(j -> j == i ? 0.0 : abs(A[i, j]), 1:n)  # Row sum
        radius2 = sum(k -> k == i ? 0.0 : abs(A[k, i]), 1:n)  # Column sum
        push!(disks, (center=center, radius=min(radius1, radius2)))
    end
    return disks
end
function match_disks_to_eigenvalues(A)
    n = size(A, 1)

    # Get statistical eigenvalue bounds (ordered λ₁ largest to λₙ smallest)
    stat_bounds = all_eigenvalue_bounds(A)

    # Get Gershgorin disks
    disks = gershgorin_disks(A)

    # Compute disk intervals and midpoints
    disk_intervals = [(d.center - d.radius, d.center + d.radius) for d in disks]
    disk_midpoints = [d.center for d in disks]

    # Compute statistical midpoints (for eigenvalue ordering)
    stat_midpoints = [(b.lower + b.upper) / 2 for b in stat_bounds]

    # Sort eigenvalues by their statistical midpoints (largest first)
    eig_order = sort(stat_midpoints,by = real,rev=true)

    # Sort disks by their centers (largest first)
    disk_order = sortperm(disk_midpoints, rev=true)

    # Match largest eigenvalue to largest disk center, etc.
    matches = Vector{Dict}(undef, n)
    for (idx, k) in enumerate(eig_order)
        disk_idx = disk_order[idx]
        disk = disks[disk_idx]
        stat_bound = stat_bounds[idx]

        # Create intersection of bounds
        lower_intersect = max(stat_bound.lower, disk.center - disk.radius)
        upper_intersect = min(stat_bound.upper, disk.center + disk.radius)

        # Ensure valid interval
        lower_intersect = min(lower_intersect, upper_intersect)
        upper_intersect = max(lower_intersect, upper_intersect)

        matches[idx] = Dict(
            :eigenvalue => stat_bound.lambda,
            :stat_lower => stat_bound.lower,
            :stat_upper => stat_bound.upper,
            :disk_center => disk.center,
            :disk_radius => disk.radius,
            :intersect_lower => lower_intersect,
            :intersect_upper => upper_intersect
        )
    end

    # Print all matches with type information
    display(A)
    for k in 1:n
        m = matches[k]
        println("\nEigenvalue: ", m[:eigenvalue])
        println("Statistical bounds: [", m[:stat_lower], ", ", m[:stat_upper], "]")
        println("Matched disk: center=", m[:disk_center], " radius=", m[:disk_radius])
        println("Intersection: [", m[:intersect_lower], ", ", m[:intersect_upper], "]")
    end

    return matches
end

function main()
    A = [
        -10.0 8.0 8.0 8.0;
         8.0 2.0 0.0 0.0;
         8.0 0.0 3.0 0.0;
         8.0 0.0 0.0 4.0
        ]

    all_bounds = all_eigenvalue_bounds(A)
    disks = gershgorin_disks(A)

    println(A)
    for b in all_bounds
        @printf("%s ∈ (%5.3f, %5.3f)\n", b.lambda, b.lower, b.upper)
    end
    for d in disks
        @printf("Center : %s ,Radius: %s\n", d.center, d.radius)
    end
    draw_gershgorin_disks_and_bounds("name",disks,all_bounds)
    
    matches = match_disks_to_eigenvalues(A)

    # Access results for λ₁ (largest eigenvalue)
    #= λ1_match = matches[1]
    println("Matched disk center: ", λ1_match[:disk_center])
    println("Intersection bounds: [", λ1_match[:intersect_lower], ", ", λ1_match[:intersect_upper], "]")

    # Print all matches
    for k in 1:size(A, 1)
        m = matches[k]
        println("\nEigenvalue: ", m[:eigenvalue])
        println("Statistical bounds: [", m[:stat_lower], ", ", m[:stat_upper], "]")
        println("Matched disk: center=", m[:disk_center], " radius=", m[:disk_radius])
        println("Intersection: [", m[:intersect_lower], ", ", m[:intersect_upper], "]")
    end =#


end
function drawBounds(name,match)
    colors = resample(ColorSchemes.turbo, length(match)*6)
    plt = plot(aspect_ratio=1, title="Bounds & Discs & Zipped Bounds", xlabel="X", ylabel="Y", legend=false)

    for (i, m) in enumerate(match)
        c = m[:disk_center]
        r = m[:disk_radius]
        θ = range(0, 2π, length=200)
        x = real(c) .+ r .* cos.(θ)
        y = imag(c) .+ r .* sin.(θ)
        plot!(x, y, label="", fill=(0.2, colors[((i-1)*6)+1]), linecolor=colors[((i-1)*6)+1])
        scatter!([real(c)], [imag(c)], color=colors[((i-1)*6)+1], markersize=4)
    
        #Static
        y = length(match) - (2*i) + 1

        plot!([m[:stat_lower], m[:stat_upper]], [y, y], label="$(m[:eigenvalue])-s", color=colors[((i-1)*6)+3], linewidth=3)
        vline!([m[:stat_lower]], line=(:dot, 2, colors[((i-1)*6)+3]), label="", alpha=0.8)
        vline!([m[:stat_upper]], line=(:dot, 2, colors[((i-1)*6)+3]), label="", alpha=0.8)

        #Intersection
        plot!([m[:intersect_lower], m[:intersect_upper]], [y-1, y-1], label="$(m[:eigenvalue])-i", color=colors[((i-1)*6)+5], linewidth=3)
        vline!([m[:intersect_lower]], line=(:dot, 2, colors[((i-1)*6)+5]), label="", alpha=0.8)
        vline!([m[:intersect_upper]], line=(:dot, 2, colors[((i-1)*6)+5]), label="", alpha=0.8)

    end
    file_name = "images/$name.png"
    savefig(plt, file_name)
    println("Görsel kaydedildi: $file_name")
end

function draw_gershgorin_disks_and_bounds(name,disks, bounds)
    colors = palette(:tab10)
    plt = plot(aspect_ratio=1, title="Eigenvalue Bounds & Discs", xlabel="X", ylabel="Y", legend=false)

    for d in disks
        c = d.center
        r = d.radius
        θ = range(0, 2π, length=200)
        x = real(c) .+ r .* cos.(θ)
        y = imag(c) .+ r .* sin.(θ)
        plot!(x, y, label="", fill=(0.2, :blue), linecolor=:blue)
        scatter!([real(c)], [imag(c)], color=:red, markersize=4)
    end

    for (i, (lambda, lower, upper)) in enumerate(bounds)
        color = colors[(i-1)%length(colors)+1]
        y = length(bounds) - i + 1
        #plot!([lower, upper], [y, y], label=lambda, color=color, linewidth=3)
        # Add vertical dotted lines at the interval boundaries
        vline!([lower], line=(:dot, 2, color), label="", alpha=0.8)
        vline!([upper], line=(:dot, 2, color), label="", alpha=0.8)
    end
    file_name="images / "+name+".png"
    savefig(plt, file_name)
    println("Görsel kaydedildi: $file_name")
end

function test()
    # Symmetric Positive Definite (2×2)
    # Eigenvalues: 4.0, 2.0
    A1 = [
        3.0 1.0;
        1.0 3.0
    ]

    # Non-Symmetric Complex Eigenvalues (2×2)
    # Eigenvalues: ±1.4142im
    A2 = [
        0.0 -2.0;
        1.0 0.0
    ]

    # Diagonal Matrix (3×3)
    # Eigenvalues: 4.0, 1.0, 7.0
    A3 = [
        4.0 0.0 0.0;
        0.0 1.0 0.0;
        0.0 0.0 7.0
    ]

    # Symmetric Indefinite (3×3)
    # Eigenvalues: ≈8.41, 1.53, -1.94
    A4 = [
        2.0 3.0 0.0;
        3.0 1.0 4.0;
        0.0 4.0 5.0
    ]

    # Non-Symmetric Defective Eigenvalue (3×3)
    A5 = [
        5.0 1.0 0.0;
        0.0 5.0 1.0;
        0.0 0.0 5.0
    ]

    # Weakly Diagonally Dominant (4×4)
    A6 = [
        8.0 1.0 0.0 1.0;
        0.0 6.0 2.0 0.0;
        1.0 0.0 9.0 1.0;
        -1.0 0.0 1.0 7.0
    ]

    # Non-Symmetric Complex Eigenvalues (4×4)
    A7 = [
        1.0 0.0 0.0 -5.0;
        0.0 2.0 0.0 0.0;
        0.0 0.0 3.0 0.0;
        1.0 0.0 0.0 4.0
    ]#Radius 0 olması durumu eşleştirmeyi bozuyor

    # Symmetric Tridiagonal Positive Definite (5×5)
    A8 = [
        4.0 1.0 0.0 0.0 0.0;
        1.0 5.0 1.0 0.0 0.0;
        0.0 1.0 6.0 1.0 0.0;
        0.0 0.0 1.0 7.0 1.0;
        0.0 0.0 0.0 1.0 8.0
    ]

    # Triangular Matrix (5×5)
    A9 = [
        10.0 1.0 0.0 0.0 0.0;
        0.0 9.0 2.0 0.0 0.0;
        0.0 0.0 8.0 3.0 0.0;
        0.0 0.0 0.0 7.0 4.0;
        0.0 0.0 0.0 0.0 6.0
    ]

    # Near-Singular (5×5)
    A10 = [
        1.0 0.0 0.0 0.0 100.0;
        0.0 2.0 0.0 0.0 0.0;
        0.0 0.0 3.0 0.0 0.0;
        0.0 0.0 0.0 4.0 0.0;
        0.0 0.0 0.0 0.0 0.01
    ]
    matrices = [A1, A2, A3, A4, A5, A6, A7, A8, A9, A10]

    for (i, A) in enumerate(matrices)
        println("\nMatrix A$i results:")
        matches = match_disks_to_eigenvalues(A)  # Convert to Float64 matrix
        drawBounds("A$i",matches)
    end
end
#main()
test()