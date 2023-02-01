 
using LinearAlgebra
using ExaPF
using DelimitedFiles

datafile = "/Users/leonard/opf/case9.m"
DEST_DIR = "/Users/leonard/opf"
RESOLUTION = 1000

"""
Plot feasibility space w.r.t. control u1 (given index iu) and contorl u2 (given by index ju).
"""
function compute_feasibility(datafile, iu, ju; resolution=10, load_factor=1.0)
    # Instantiate model
    polar = PolarForm(datafile)
    stack = ExaPF.NetworkStack(polar)
    # Power flow equations
    pflow = ExaPF.PowerFlowBalance(polar) ∘ ExaPF.PolarBasis(polar)
    # Operational constraints
    opcons = [
        ExaPF.VoltageMagnitudeBounds(polar) ∘ ExaPF.PolarBasis(polar),
        ExaPF.PowerGenerationBounds(polar) ∘ ExaPF.PolarBasis(polar),
        ExaPF.LineFlows(polar) ∘ ExaPF.PolarBasis(polar),
    ]
    
    nx = ExaPF.number(polar, State())
    nu = ExaPF.number(polar, Control())
    @assert (iu <= nu) && (ju <= nu)
    mapu = ExaPF.mapping(polar, Control())
    bounds = [ExaPF.bounds(polar, c) for c in opcons]
    values = [zeros(length(c)) for c in opcons]

    feasibility = zeros(Int, resolution+1, resolution+1)

    algo = NewtonRaphson(tol=1e-10)
    buffer = ExaPF.NLBuffer{Vector{Float64}}(nx)

    jac = ExaPF.Jacobian(polar, pflow, State())

    stack.pload .*= load_factor
    stack.qload .*= load_factor
    x0 = copy(stack.input)
    # We set voltage magnitude found at solution of OPF
    x0[mapu[2]] = 1.1
    x0[mapu[3]] = 1.1
    ExaPF.set_params!(jac, stack)

    i_, j_ = mapu[iu], mapu[ju]
    lbs, ubs = ExaPF.bounds(polar, stack)
    u1 = range(lbs[i_], ubs[i_]; length=resolution+1)
    u2 = range(lbs[j_], ubs[j_]; length=resolution+1)

    for (i, p1) in enumerate(reverse(u1)), (j, p2) in enumerate(u2)
        copyto!(stack.input, x0)
        stack.input[i_] = p1
        stack.input[j_] = p2

        # Solve power flow equations.
        conv = ExaPF.nlsolve!(algo, jac, stack; nl_buffer=buffer)
        feasibility[i, j] += !conv.has_converged
        # Evaluate operational constraints
        for (k, cons) in enumerate(opcons)
            l, u = bounds[k]
            cons(values[k], stack)
            is_feas = all(l .<= values[k] .<= u)
            feasibility[i, j] += 2^k * !is_feas
        end
    end

    return feasibility
end

pg2 = 4
pg3 = 5

feas = compute_feasibility(datafile, pg2, pg3; resolution=RESOLUTION, load_factor=1.0)

writedlm(joinpath(DEST_DIR, "feasibility.txt"), feas)
