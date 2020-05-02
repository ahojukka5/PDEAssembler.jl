# This file is a part of JuliaFEM.
# License is MIT: see https://github.com/JuliaFEM/PDEAssembles.jl/blob/master/LICENSE

if VERSION < v"1.0.0"
    push!(LOAD_PATH, joinpath(homedir(), ".julia", "dev"))
    using Base.Test
else
    using Test
end

using PDEAssembler
using PDEAssembler: get_unit_square, Poisson

field_elements, boundary_elements = get_unit_square()

field_problem = Problem(Poisson, "Poisson problem in 1x1 square", 1)
update!(field_elements, "source", 10.0)
update!(field_elements, "density", 6.0)
add_elements!(field_problem, field_elements)

boundary_problem = Problem(Poisson, "Poisson boundary", 1)
update!(boundary_elements, "fixed u", 0.0)
update!(boundary_elements, "density", 0.0)
add_elements!(boundary_problem, boundary_elements)

problems = (field_problem, boundary_problem)
M, C, K, f = get_global_matrices(problems, 0.0)
u = cholfact(Symmetric(K)) \ full(f)
X = first(problems)("geometry", 0.0)

N = length(u)
x = [X[i][1] for i=1:N]
y = [X[i][2] for i=1:N]

#using Plots
#surface(x, y, u)

@test isapprox(maximum(u), 0.7381696594268394)
