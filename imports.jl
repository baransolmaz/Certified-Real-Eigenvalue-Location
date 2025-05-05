import Pkg
Pkg.add("Plots")
Pkg.add("LinearAlgebra")
Pkg.add("GR")   # Otomatik görüntüleme için
ENV["GKSwstype"] = "100"  # Pencere açılmasını sağlar (X11 penceresi gibi düşün)

# GR'yi yeniden derle:
Pkg.build("GR")