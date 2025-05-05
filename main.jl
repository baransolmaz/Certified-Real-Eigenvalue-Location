using Plots
gr()


function gershgorin_disks(A::AbstractMatrix{<:Number})
    n, m = size(A)
    if n != m
        error("Matris kare değil!")
    end

    disks = []
    for i in 1:n
        center = A[i, i]
        radius1 = sum(abs(A[i, j]) for j in 1:n if j != i)
        radius2 = sum(abs(A[k, i]) for k in 1:n if k != i)
        push!(disks, (center=center, radius=min(radius1,radius2))) # minimumu al
    end
    return disks
end

function draw_gershgorin_disks(disks)
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

    savefig(plt, "disk_images/gershgorin.png")
    println("Görsel kaydedildi: disk_images/gershgorin.png")
end

# Main fonksiyonu
function main()
    A = [1.0 2.0 3.0;
        1.0 2.0 3.0;
        1.0 1.0 3.0]

    disks = gershgorin_disks(A)
    for d in disks
        println("Center: $(d.center) | Radius: $(d.radius)")
    end
    draw_gershgorin_disks(disks)






end












#Runner
main()