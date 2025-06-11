using Plots
using LinearAlgebra
gr()

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

# Main fonksiyonu
function main()
    B3 = Matrix{Rational{Int}}([5/4 1 3/4 1/2 1/4;
        1 0 0 0 0;
        -1 1 0 0 0;
        0 0 1 3 0;
        0 0 0 1/2 5])

    input = B3
    pa = charpoly_faddeev_leverrier(input)
    power_sums = newton_girard_power_sums(pa)
    h1 = setHermite_1(length(pa) - 1, power_sums)
    println("H1")
    cp = companion_matrix(pa)
    signH1 = signature(h1)
    println("Sign H1 : $(signH1)")

    if check_matrix_requirements(input)
        disks = gershgorin_disks(input)
        draw_gershgorin_disks_and_bounds("all_disks",disks)
        containedDisks =[]
        points = []
        for d in disks
            println("Center: $(d.center) | Radius: $(d.radius)")

            if d.radius == 0
                println("EigenValue at $(d.center)")
            else
                g(x) = (x - d.center * I)^2 - (d.radius * I)^2

                hg = h1 * g(cp)
                signHg = signature(hg)#!!
                println("Sign Hg : $(signHg)")

                if signH1 != signHg

                    push!(containedDisks, (center=d.center, radius=d.radius,fill = true))
                    push!(points, d.center - d.radius)
                    push!(points, d.center + d.radius)
                    println("Circle contains at least one real eigen")
                else
                    push!(containedDisks, (center=d.center, radius=d.radius, fill=false))
                    println("Circle does not contain any real eigen")
                end



            end
        end
        plt = scan_gershgorin_disks("scanned",containedDisks)
        draw_gershgorin_disks_and_bounds("remain_disks", containedDisks)
        sort!(points)
        intervals = []
        for i in 1:(length(points)-1)
            d1 = points[i]
            d2 = points[i+1]
            println("Interval $d1 - $d2")
            g(x) = (x - d1 * I) * (x - d2 * I)

            hg = h1 * g(cp)
            signHg = signature(hg)#!!
            println("Sign Hg : $(signHg)")

            if signH1 != signHg
                push!(intervals,(startP = d1,endP=d2,isExist = true))
                println("Interval contains at least one real eigen")
            else
                push!(intervals, (startP=d1, endP=d2, isExist=false))
                println("Interval does not contain any real eigen")
            end

        end

        draw_intervals("intervals", plt ,intervals)




    else
        println("---------- Matris Uygun Değil! ----------")
    end





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

function check_matrix_requirements(A::AbstractMatrix{<:Number})
    n, m = size(A)
    return n == m
end

function gershgorin_disks(A::AbstractMatrix{<:Number})
    n, m = size(A)
    disks = []
    for i in 1:n
        center = A[i, i]
        radius1 = sum(abs(A[i, j]) for j in 1:n if j != i)
        radius2 = sum(abs(A[k, i]) for k in 1:n if k != i)
        push!(disks, (center=center, radius=radius1)) 
    end
    return disks
end

function draw_gershgorin_disks_and_bounds(name, disks)
    colors = palette(:tab10)
    plt = plot(aspect_ratio=1, title="Discs", xlabel="X", ylabel="Y", legend=false)

    
    for d in disks
        c = d.center
        r = d.radius
        θ = range(0, 2π, length=200)
        x = real(c) .+ r .* cos.(θ)
        y = imag(c) .+ r .* sin.(θ)
        plot!(x, y, label="", fill=(0.2, :blue), linecolor=:blue)
        scatter!([real(c)], [imag(c)], color=:red, markersize=4)
    end
    #=for (i, (lambda, lower, upper)) in enumerate(bounds)
        color = colors[(i-1)%length(colors)+1]
        y = length(bounds) - i + 1
        #plot!([lower, upper], [y, y], label=lambda, color=color, linewidth=3)
        # Add vertical dotted lines at the interval boundaries
        vline!([lower], line=(:dot, 2, color), label="", alpha=0.8)
        vline!([upper], line=(:dot, 2, color), label="", alpha=0.8)
    end=#
    file_name = "images/$name.png"
    savefig(plt, file_name)
    println("Görsel kaydedildi: $file_name")
end

function draw_intervals(name, plt, intervals)
    colors = palette(:tab10)
    ymin, ymax = ylims()  # Get current y-axis limits

    for i in intervals
        vline!([i.startP], line=(:dot, 2, :black), label="", alpha=0.8)
        vline!([i.endP], line=(:dot, 2, :black), label="", alpha=0.8)

        if i.isExist
            n_lines = 20
            width = i.endP - i.startP
            height = ymax - ymin

            # Calculate diagonal scanning parameters
            diagonal_length = √(width^2 + height^2)
            spacing = diagonal_length / (n_lines + 1)

            # Draw diagonal lines at 45° angle
            for j in 1:n_lines
                offset = j * spacing

                # Calculate starting point (bottom-left)
                start_x = i.startP + offset
                start_y = ymin

                # Calculate ending point (top-right)
                end_x = i.startP
                end_y = ymin + offset

                # Draw the diagonal segment within the interval
                if offset <= width
                    # First segment: bottom-left to top boundary
                    plot!([start_x, i.startP], [ymin, end_y],
                        color=:blue, linewidth=1.5, alpha=0.6, label="")
                else
                    # Second segment: left boundary to top-right
                    plot!([i.startP, i.startP + (offset - width)],
                        [ymin + (offset - width), end_y],
                        color=:blue, linewidth=1.5, alpha=0.6, label="")
                end
            end
        end
    end

    file_name = "images/$name.png"
    savefig(plt, file_name)
    println("Görsel kaydedildi: $file_name")
end

function scan_gershgorin_disks(name, disks)
    colors = palette(:tab10)
    plt = plot(aspect_ratio=1, title="Discs", xlabel="X", ylabel="Y", legend=false)

    for d in disks
        c = d.center
        r = d.radius
        cx = real(c)
        cy = imag(c)
        θ = range(0, 2π, length=200)
        x = cx .+ r .* cos.(θ)
        y = cy .+ r .* sin.(θ)
        plot!(x, y, label="", fill=(0.2, :blue), linecolor=:blue)
        scatter!([cx], [cy], color=:red, markersize=4)

        n_lines = 20
        # Get x-range from disk properties or default to full disk
        a = hasproperty(d, :a) ? d.a : cx - r
        b = hasproperty(d, :b) ? d.b : cx + r

        # Ensure a < b
        a, b = min(a, b), max(a, b)
        # Clamp to disk's x-range
        a = max(a, cx - r)
        b = min(b, cx + r)

        x_vals = range(a, b, length=n_lines)
        line_color = d.fill ? :green : :red

        # Reverse scanning direction for fill=false
        if !d.fill
            x_vals = reverse(x_vals)
        end

        for x in x_vals
            dx = x - cx
            offset_sq = r^2 - dx^2
            offset_sq < 0 && continue  # Skip if outside disk

            offset = sqrt(offset_sq)
            y1 = cy - offset
            y2 = cy + offset

            plot!([x, x], [y1, y2],
                color=line_color,
                linewidth=1.5,
                alpha=0.6,
                label="")
        end
    end

    file_name = "images/$name.png"
    savefig(plt, file_name)
    println("Görsel kaydedildi: $file_name")
end

#Runner
main()

function circleShape(h, k, r) 
    theta= LinRange(0, 2*π, 500)
    h .+r*sin.(theta), k.+r*cos.(theta)
end
plot(circleShape(0, 0, 1), seriestype=[:shape,], lw=0.5,
    c=:blue, linecolor=:black,
    legend=false, fillalpha=0.2, aspect_ratio=1)