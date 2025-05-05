using Plots
gr()  # GR backend ile çalış, otomatik pencere açılır


function gershgorin_disks(A::AbstractMatrix{<:Number})
    n, m = size(A)
    if n != m
        error("Matris kare değil!")
    end

    disks = []
    for i in 1:n
        center = A[i, i]
        radius = sum(abs(A[i, j]) for j in 1:n if j != i)
        push!(disks, (center=center, radius=radius))
    end
    return disks
end

function ciz_gershgorin_disks(disks)
    plt = plot(aspect_ratio=1, legend=false, title="Gershgorin Diskleri")

    for d in disks
        c = d.center
        r = d.radius
        θ = range(0, 2π, length=200)
        x = real(c) .+ r .* cos.(θ)
        y = imag(c) .+ r .* sin.(θ)
        plot!(x, y, label="", fill=(0.2, :blue), linecolor=:blue)
        scatter!([real(c)], [imag(c)], color=:red, markersize=4)
    end

    savefig(plt, "gershgorin.png")
    println("Görsel kaydedildi: gershgorin.png")
end

# Main fonksiyonu
function main()
    A = [4.0 2.0 0.0;
        1.0 3.0 3.0;
        0.0 1.0 2.0]

    disks = gershgorin_disks(A)
    for d in disks
        println("Center: $(d.center) | Radius: $(d.radius)")
    end
    ciz_gershgorin_disks(disks)
end

main()