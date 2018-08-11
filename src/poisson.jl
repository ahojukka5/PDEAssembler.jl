# This file is a part of JuliaFEM.
# License is MIT: see https://github.com/JuliaFEM/PDEAssembler.jl/blob/master/LICENSE

# Poisson's problem as example

struct Poisson <: FieldProblem end

FEMBase.get_unknown_field_name(::Problem{Poisson}) = "u"

function FEMBase.assemble_elements!(problem::Problem{Poisson},
                                    assembly::Assembly,
                                    elements::Vector{Element{E}},
                                    time::Float64) where E

    bi = BasisInfo(E)
    ndofs = length(bi)
    Ke = zeros(ndofs, ndofs)
    fe = zeros(ndofs)

    for element in elements
        fill!(Ke, 0.0)
        fill!(fe, 0.0)
        for ip in get_integration_points(element)
            J, detJ, N, dN = element_info!(bi, element, ip, time)
            k = 1.0
            if haskey(element, "coefficient")
                k = element("coefficient", ip, time)
            end
            Ke += ip.weight * k*dN'*dN * detJ
            if haskey(element, "source")
                f = element("source", ip, time)
                fe += ip.weight * N'*f * detJ
            end
        end
        gdofs = get_gdofs(problem, element)
        add!(assembly.K, gdofs, gdofs, Ke)
        add!(assembly.f, gdofs, fe)
    end

    return nothing

end

function FEMBase.assemble_elements!(problem::Problem{Poisson},
                                    assembly::Assembly,
                                    elements::Vector{Element{E}},
                                    time::Float64) where E<:Union{Seg2,Seg3}

    bi = BasisInfo(E)
    ndofs = length(bi)
    Ce = zeros(ndofs, ndofs)
    fe = zeros(ndofs)
    ge = zeros(ndofs)

    for element in elements
        fill!(Ce, 0.0)
        fill!(fe, 0.0)
        fill!(ge, 0.0)
        for ip in get_integration_points(element)
            J, detJ, N, dN = element_info!(bi, element, ip, time)
            if haskey(element, "flux")
                f = element("flux", ip, time)
                fe += ip.weight * N'*f * detJ
            end
            if haskey(element, "fixed u")
                f = element("fixed u", ip, time)
                Ce += ip.weight * N'*N * detJ
                ge += ip.weight * N'*f * detJ
            end
        end
        gdofs = get_gdofs(problem, element)
        add!(assembly.C1, gdofs, gdofs, Ce')
        add!(assembly.C2, gdofs, gdofs, Ce)
        add!(assembly.f, gdofs, fe)
        add!(assembly.g, gdofs, ge)
    end

    return nothing

end

"""
    get_unit_square(nel_x=10, nel_y=10; lx=1.0, ly=1.0)

Return elements and boundary elements for a standard test problem, unit square.
"""
function get_unit_square(nel_x=20, nel_y=20; lx=1.0, ly=1.0)

    nnodes_x = nel_x+1
    nnodes_y = nel_y+1
    nnode = nnodes_x*nnodes_y
    nodemap = reshape(1:nnode, nnodes_x, nnodes_y)
    nodes_1 = vec(nodemap[1:nnodes_x-1, 1:nnodes_y-1])
    nodes_2 = vec(nodemap[2:nnodes_x, 1:nnodes_y-1])
    nodes_3 = vec(nodemap[2:nnodes_x, 2:nnodes_y])
    nodes_4 = vec(nodemap[1:nnodes_x-1, 2:nnodes_y])

    X = Dict{Int64, Vector{Float64}}()

    # Create nodes
    nid = 1
    for y in linspace(0, ly, nnodes_y)
        for x in linspace(0, lx, nnodes_x)
            X[nid] = [x, y]
            nid += 1
        end
    end

    # Create elements for volume
    field_elements = []
    for c in zip(nodes_1, nodes_2, nodes_3, nodes_4)
        push!(field_elements, Element(Quad4, collect(c)))
    end

    # add boundary elements to the left side of domain
    boundary_elements = []
    nids = 1:nel_x+1:(nel_x+1)*(nel_y+1) # LEFT and RIGHT
    for c in zip(nids[1:end-1], nids[2:end])
        connectivity = collect(c)
        push!(boundary_elements, Element(Seg2, connectivity))
        push!(boundary_elements, Element(Seg2, connectivity+nel_x))
    end
    nids = 1:nel_x+1 # BOTTOM and TOP
    for c in zip(nids[1:end-1], nids[2:end])
        connectivity = collect(c)
        push!(boundary_elements, Element(Seg2, connectivity))
        push!(boundary_elements, Element(Seg2, connectivity+(nel_x+1)*nel_y))
    end

    update!(field_elements, "geometry", X)
    update!(boundary_elements, "geometry", X)

    return field_elements, boundary_elements
end
