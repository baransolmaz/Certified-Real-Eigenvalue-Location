# Certified Real Eigenvalue Location

## ℹ️ About This Project

This project focuses on the certified computation of real eigenvalues for matrices. It provides numerical and symbolic tools to locate and verify eigenvalues with high precision. The project is designed to support research in numerical linear algebra and symbolic computation.

## 📁 Files

- `imports.jl` — contains all package imports.
- `num_main.jl` — the numerical main program that includes `imports.jl` before executing numerical computations.
- `sym_main.jl` — the symbolic main program that includes `imports.jl` before executing symbolic computations.

## ▶️ Running the Project

1. Make sure [Julia](https://julialang.org/downloads/) is installed.
2. Open a terminal and navigate to the project directory.
3. Run the numerical program with:
    ```sh
    julia imports.jl
    julia num_main.jl
    ```
4. Run the symbolic program with:
    ```sh
    julia imports.jl
    julia sym_main.jl
    ```