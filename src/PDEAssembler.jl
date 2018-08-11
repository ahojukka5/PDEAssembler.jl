# This file is a part of JuliaFEM.
# License is MIT: see https://github.com/JuliaFEM/PDEAssembler.jl/blob/master/LICENSE

module PDEAssembler

using Reexport
@reexport using FEMBase
using TimerOutputs

include("poisson.jl")

"""
    get_global_matrices(problems, time)

Assemble a set of problems and return assembles global matrices M for mass,
C for damping, K for stiffness, and right hand side vector f. Semi-discrete
equation of motion is then

    Mu'' + Cu' + Ku = f

Problems can describe a PDE in whole domain or in part of it. It can also describe
a boundary of the domain.
"""

function get_global_matrices(problems, time)

    @timeit "assemble problems" for problem in problems
        empty!(problem.assembly)
        assemble!(problem, time)
        assemble_mass_matrix!(problem, time)
    end

    @timeit "combine field problems" begin
        M = SparseMatrixCOO()
        K = SparseMatrixCOO()
        f = SparseMatrixCOO()
        C1 = SparseMatrixCOO()
        C2 = SparseMatrixCOO()
        D = SparseMatrixCOO()
        g = SparseMatrixCOO()
        for problem in problems
            A = problem.assembly
            append!(M, A.M)
            append!(K, A.K)
            append!(K, A.Kg)
            append!(f, A.f)
            append!(C1, A.C1)
            append!(C2, A.C2)
            append!(D, A.D)
            append!(g, A.g)
        end
    end

    @timeit "create LinearSystem(dim)" begin

        # Based on the nodal/dof mapping of PDEs, there might be empty rows in
        # stiffness matrix giving singular solution. For that reason we have
        # to find nonzero rows and add 1 to diagonal of those to indicate that
        # this particular dof is "disabled"

        dim = size(K, 1)
        nonzero_rows = zeros(dim)
        for i in K.I
            nonzero_rows[i] = 1.0
        end
        zero_indices = find(nonzero_rows .== 0.0)
        append!(K.I, zero_indices)
        append!(K.J, zero_indices)
        append!(K.V, ones(zero_indices))

        ls = LinearSystem(dim)
        ls.M = sparse(M, dim, dim)
        ls.K = sparse(K, dim, dim)
        ls.f = sparse(f, dim, 1)
        ls.C1 = sparse(C1, dim, dim)
        ls.C2 = sparse(C2, dim, dim)
        ls.D = sparse(D, dim, dim)
        ls.g = sparse(g, dim, 1)
    end

    # Now we have a linear system of equations to solve problem
    #
    # M*u'' + C*u' + K*u + C1*la = f
    # C2*u + D*la = g
    #
    # for u and la. In a case C1 == C2 and D = 0, we have a saddle
    # point problem where it is possible to either eliminate Lagrange
    # multipliers or solve constrained optimization problem using penalty
    # method. In simples kind of boundary conditions, constraint matrix is
    # a diagonal matrix having 1.0 in those diagonal entries that are
    # constrained.

    if ls.C1 == ls.C2' && nnz(ls.D) == 0
        ls.K += 1e36*ls.C1'*ls.C1
        ls.f += 1e36*ls.C1'*ls.g
    else
        error("Cannot handle boundary conditions using penalty method.")
    end

    # TODO: we are missing damping matrix ?!
    return ls.M, spzeros(dim,dim), ls.K, ls.f
end

export get_global_matrices

end
