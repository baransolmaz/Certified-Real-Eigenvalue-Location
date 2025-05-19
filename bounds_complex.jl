using LinearAlgebra
using Printf
using Plots

#=
Bounds for Elgenvalues Using Traces*
    Hemy WoIkowicz
    George P. H. Styan
=#
#=1.1
Ref[9,  p. 144]
M.  Marcus  and  H. Mint, A Survey of Matrix Theory and Matrix Inequalities, 
Prindle,  Weber &  Schmidt,  Boston, 1964.
=#
function compute_R_C(A::AbstractMatrix)
    R = maximum(sum(abs.(A), dims=2))
    C = maximum(sum(abs.(A), dims=1))
    return R, C
end

#=1.2
Ref  [9,  p.  145]
M.  Marcus  and  H. Mint, A Survey of Matrix Theory and Matrix Inequalities, 
Prindle,  Weber &  Schmidt,  Boston, 1964.
=#
function brauer_inequality(lambda::Number, A::AbstractMatrix)
    return abs.(lambda) <= minimum(compute_R_C(A))
end

#1.3
function eq1_3(lambda::Number, A::AbstractMatrix)
    #TODO
end

#=1.4
Ref  [1,  p.  133]
A. R  Amir-Moez and A. L.  Fass, Elements of Linear  Spaces, Pergamon, Oxford, 1962
=#
function matrix_B(A::AbstractMatrix)
    return (1 / 2) * (A + adjoint(A))
end
#=1.4
Ref  [1,  p.  133]
A. R  Amir-Moez and A. L.  Fass, Elements of Linear  Spaces, Pergamon, Oxford, 1962
=#
function matrix_C(A::AbstractMatrix)
    return (1 / 2) * (A - adjoint(A)) / 1im
end

#1.5
function eq1_5(A::AbstractMatrix)
    #TODO
end

#1.6
function eq1_6(A::AbstractMatrix)
    #TODO
end

#1.7
function eq1_7(A::AbstractMatrix)
    #TODO
end
#=1.8 S= Standart Deviation 
[4,p. 2271
F. A. Graybill, Zntroductionto Matrices with Applicutions in Statistics, Wadsworth, Belmont, Cal&, 196
=#
function compute_s(A::AbstractMatrix)
    n, m = size(A)
    return sqrt(((tr(A * A) - ((tr(A)^2)) / n) / n))
end

#1.10
#=
H. Wolkowicz and G. P. H. Styan, Extensions of Samuelson’s inequality, Aw.
stat. 33:143-144 (1979)
=#
function min_bound(A::AbstractMatrix)
    n, m = size(A)
    m = mean_m(A)
    s = compute_s(A)
    return m + (s / sqrt(n - 1))
end

#1.10
#=
H. Wolkowicz and G. P. H. Styan, Extensions of Samuelson’s inequality, Aw.
stat. 33:143-144 (1979)
=#
function max_bound(A::AbstractMatrix)
    n, m = size(A)
    m = mean_m(A)
    s = compute_s(A)
    return m + (s * sqrt(n - 1))
end

#1.11
#=
H. Wolkowicz and G. P. H. Styan, Extensions of Samuelson’s inequality, Aw.
stat. 33:143-144 (1979)
=#
function mean_m(A::AbstractMatrix)
    n, m = size(A)
    return tr(A) / n
end

struct EigenBounds
    real::Vector{Tuple{Float64, Float64}}
    imag::Vector{Tuple{Float64, Float64}}
    modulus::Vector{Tuple{Float64, Float64}}
end

function eigen_bounds(A::AbstractMatrix{T}) where {T<:Number}
    n = checksquare(A)
    B = (A + A')/2
    C = (A - A')/(2im)
    
    # Real part bounds
    trB = real(tr(B))
    trB² = real(tr(B^2))
    m_b = trB / n
    s_b = sqrt((trB² - trB^2/n)/n)
    
    # Imaginary part bounds
    trC = real(tr(C))
    trC² = real(tr(C^2))
    m_c = trC / n
    s_c = sqrt((trC² - trC^2/n)/n)
    
    # Modulus bounds
    trAA = real(tr(A'*A))
    m_a = abs(tr(A)/n)
    s_a = sqrt((trAA/n) - m_a^2)
    
    # Precompute coefficients
    k_vals = 1:n
    sqrt_terms_real = [sqrt((k-1)/(n-k+1)) for k in k_vals]
    sqrt_terms_imag = [sqrt((n-k)/k) for k in k_vals]
    
    # Generate bounds
    real_bounds = [(m_b - s_b*sqrt_terms_real[k], m_b + s_b*sqrt_terms_imag[k]) for k in k_vals]
    imag_bounds = [(m_c - s_c*sqrt_terms_real[k], m_c + s_c*sqrt_terms_imag[k]) for k in k_vals]
    mod_bounds = [(max(m_a - s_a*sqrt(n-1), 0), sqrt(trAA/n) + s_a*sqrt((n-1)/n)) for _ in k_vals]
    
    # Sort bounds descendingly
    sort!(real_bounds, rev=true)
    sort!(imag_bounds, rev=true)
    sort!(mod_bounds, rev=true)
    
    EigenBounds(real_bounds, imag_bounds, mod_bounds)
end

# Helper functions
function checksquare(A)
    m, n = size(A)
    m == n || throw(DimensionMismatch("Matrix must be square"))
    n
end

function print_bounds(b::EigenBounds)
    println("Real Part Bounds:")
    for (i, (lo, hi)) in enumerate(b.real)
        @printf("λ_%d ∈ (%6.3f, %6.3f)\n", i, lo, hi)
    end
    
    println("\nImaginary Part Bounds:")
    for (i, (lo, hi)) in enumerate(b.imag)
        @printf("λ_%d ∈ (%6.3f, %6.3f)\n", i, lo, hi)
    end
    
    println("\nModulus Bounds:")
    for (i, (lo, hi)) in enumerate(b.modulus)
        @printf("|λ_%d| ∈ (%6.3f, %6.3f)\n", i, lo, hi)
    end
end


function draw_intervals(V::Vector{NamedTuple})
    colors = palette(:tab10)
    plt=plot(title="Eigenvalue Bounds", xlabel="Value", ylabel="Eigenvalue", legend=:right)

    for (i, (lambda, lower, upper)) in enumerate(V)
        color = colors[(i-1)%length(colors)+1]
        y = length(V) - i + 1
        plot!([lower, upper], [y, y], label=lambda, color=color, linewidth=3)
        # Add vertical dotted lines at the interval boundaries
        vline!([lower], line=(:dot, 1, color), label="", alpha=0.6)
        vline!([upper], line=(:dot, 1, color), label="", alpha=0.6)
    end

    yticks!(reverse(1:length(V)), [lambda for (lambda, _, _) in V])
    savefig(plt, "images/bounds.png")
    println("Görsel kaydedildi: images/bounds.png")
end
function print_bounds(b::EigenBounds)
    println("Real Part Bounds:")
    for (i, (lo, hi)) in enumerate(b.real)
        @printf("λ_%d ∈ (%6.3f, %6.3f)\n", i, lo, hi)
    end

    println("\nImaginary Part Bounds:")
    for (i, (lo, hi)) in enumerate(b.imag)
        @printf("λ_%d ∈ (%6.3f, %6.3f)\n", i, lo, hi)
    end

    println("\nModulus Bounds:")
    for (i, (lo, hi)) in enumerate(b.modulus)
        @printf("|λ_%d| ∈ (%6.3f, %6.3f)\n", i, lo, hi)
    end
end

function main()
    A = [1.0 3.0 5.0 9.0;
        2.0 4.0 6.0 8.0;
        1.0 2.0 3.0 4.0;
        5.0 6.0 7.0 8.0]

    B = [0.0 -1.0 0.0;
        1.0 0.0 0.0;
        0.0 0.0 1.0]

    C = [1.0 -1.0;
        1.0 1.0 ]


    bounds = eigen_bounds(B)
    print_bounds(bounds)

    #==println(C)
    for b in all_bounds
        display(b)
    end
    ==#

    #draw_intervals(all_bounds)
    #draw_gershgorin_disks_and_bounds(disks,all_bounds)

end

function gershgorin_disks(A::AbstractMatrix{<:Number})
    n, m = size(A)
    disks = []
    for i in 1:n
        center = A[i, i]
        radius1 = sum(abs(A[i, j]) for j in 1:n if j != i)
        radius2 = sum(abs(A[k, i]) for k in 1:n if k != i)
        push!(disks, (center=center, radius=min(radius1, radius2))) # minimumu al
    end
    return disks
end

function draw_gershgorin_disks_and_bounds(disks,bounds)
    colors = palette(:tab10)
    plt = plot(aspect_ratio=1,title="Eigenvalue Bounds & Discs", xlabel="X", ylabel="Y",legend=false)

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
    savefig(plt, "images/discs&bounds.png")
    println("Görsel kaydedildi: images/discs&bounds.png")
end


#Runner
main()