using LinearAlgebra

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
function s_square(A::AbstractMatrix)
    n,m= size(A) 
    return ((tr(A*A) - ((tr(A)^2))/n)/n)
end

#1.10
#=
H. Wolkowicz and G. P. H. Styan, Extensions of Samuelson’s inequality, Aw.
stat. 33:143-144 (1979)
=#
function min_bound(A::AbstractMatrix)
    n, m = size(A)
    m = mean_m(A)
    s = sqrt(s_square(A))
    return m + (s /sqrt(n-1))
end

#1.10
#=
H. Wolkowicz and G. P. H. Styan, Extensions of Samuelson’s inequality, Aw.
stat. 33:143-144 (1979)
=#
function max_bound(A::AbstractMatrix)
    n, m = size(A)
    m = mean_m(A)
    s = sqrt(s_square(A))
    return m + (s * sqrt(n - 1))
end

#1.11
#=
H. Wolkowicz and G. P. H. Styan, Extensions of Samuelson’s inequality, Aw.
stat. 33:143-144 (1979)
=#
function mean_m(A::AbstractMatrix)
    n, m = size(A)
    return tr(A)/n
end



function main()
    A = [1.0 2.0 3.0;
        1.0 2.0 3.0;
        1.0 1.0 3.0]


    print(matrix_B( A) + (1im * matrix_C(A)))




end





#Runner
main()