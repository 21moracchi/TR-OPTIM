
using ExaPF
using JuMP
using Ipopt
using LinearAlgebra
const PS = ExaPF.PowerSystem


constraint_index(cons::Vector{NonlinearConstraintRef{ScalarShape}}) = getfield.(JuMP.index.(cons), :value)

function build_opf_model(polar, buffer, solver; line_constraints=true, w)
    nbus = ExaPF.get(polar, PS.NumberOfBuses())
    ngen = ExaPF.get(polar, PS.NumberOfGenerators())
    nlines = ExaPF.get(polar, PS.NumberOfLines())

    pf = polar.network
    baseMVA = pf.baseMVA

    # Bounds
    pg_min, pg_max = PS.bounds(pf, PS.Generators(), PS.ActivePower())
    qg_min, qg_max = PS.bounds(pf, PS.Generators(), PS.ReactivePower())
    vm_min, vm_max = PS.bounds(pf, PS.Buses(), PS.VoltageMagnitude())

    # Max flow through lines
    flow_min, flow_max = PS.bounds(pf, PS.Lines(), PS.ActivePower())
    flow_max = min.(1e5, flow_max)

    # Initial values
    vm0 = buffer.vmag
    va0 = buffer.vang
    pg0 = buffer.pgen

    # Loads
    pload = w[1:nlines]#buffer.pload
    qload = w[nlines+1:2nlines]#buffer.qload

    # Bus-admittance matrix
    Ybus = pf.Ybus
    rows = Ybus.rowval
    yvals = Ybus.nzval
    g_ij = real.(yvals)
    b_ij = imag.(yvals)

    # Line characteristics (Matpower manual, Section 3.2, page 25)
    yff_re = real.(pf.lines.Yff)
    yff_im = imag.(pf.lines.Yff)
    yft_re = real.(pf.lines.Yft)
    yft_im = imag.(pf.lines.Yft)
    ytf_re = real.(pf.lines.Ytf)
    ytf_im = imag.(pf.lines.Ytf)
    ytt_re = real.(pf.lines.Ytt)
    ytt_im = imag.(pf.lines.Ytt)

    # Objective
    cost_coefs = PS.get_costs_coefficients(pf)

    bus2gen = PS.get_bus_generators(pf.buses, pf.generators, pf.bus_to_indexes)

    #=
        Build model
    =#

    opfmodel = Model(solver)

    # VARIABLES
    @variable(opfmodel, pg_min[i] <= Pg[i=1:ngen] <= pg_max[i], start=pg0[i])
    @variable(opfmodel, qg_min[i] <= Qg[i=1:ngen] <= qg_max[i])
    @variable(opfmodel, vm_min[i] <= Vm[i=1:nbus] <= vm_max[i], start=vm0[i])
    @variable(opfmodel, Va[i=1:nbus], start=va0[i])

    @variable(opfmodel, Pd[i=1:nbus])
    @variable(opfmodel, Qd[i=1:nbus])
    JuMP.fix.(Pd, pload)
    JuMP.fix.(Qd, qload)

    # Power-flow constraints
    ## active
    opfmodel.ext[:active_pf] = @NLconstraint(
        opfmodel, [b=1:nbus],
        Vm[b] * sum(
            Vm[rows[c]] * (g_ij[c] * cos(Va[b] - Va[rows[c]]) + b_ij[c] * sin(Va[b] - Va[rows[c]]))
            for c in (Ybus.colptr[b]):(Ybus.colptr[b+1]-1)
        ) == (sum(Pg[g] for g in get(bus2gen, b, Int[])) - Pd[b])
    )
    ## reactive
    opfmodel.ext[:reactive_pf] = @NLconstraint(
        opfmodel, [b=1:nbus],
        Vm[b] * sum(
            Vm[rows[c]] * (g_ij[c] * sin(Va[b] - Va[rows[c]]) - b_ij[c] * cos(Va[b] - Va[rows[c]]))
            for c in (Ybus.colptr[b]):(Ybus.colptr[b+1]-1)) == (sum(Qg[g] for g in get(bus2gen, b, Int[])) - Qd[b])
    )

    # Line constraints
    f = pf.lines.from_buses
    t = pf.lines.to_buses

    ## from lines
    yff_abs = yff_re.^2 .+ yff_im.^2
    yft_abs = yft_re.^2 .+ yft_im.^2
    yre_fr =   yff_re .* yft_re .+ yff_im .* yft_im
    yim_fr = - yff_re .* yft_im .+ yff_im .* yft_re

    opfmodel.ext[:line_fr] = @NLconstraint(
        opfmodel, [???=1:nlines],
        Vm[f[???]]^2 * (
            yff_abs[???] * Vm[f[???]]^2 + yft_abs[???] * Vm[t[???]]^2 +
            2.0 * Vm[f[???]] * Vm[t[???]] * (yre_fr[???] * cos(Va[f[???]]-Va[t[???]]) - yim_fr[???] * sin(Va[f[???]]-Va[t[???]]))
        ) <= flow_max[???]
    )

    ## to lines
    ytf_abs = ytf_re.^2 .+ ytf_im.^2
    ytt_abs = ytt_re.^2 .+ ytt_im.^2
    yre_to =   ytf_re .* ytt_re .+ ytf_im .* ytt_im
    yim_to = - ytf_re .* ytt_im .+ ytf_im .* ytt_re

    opfmodel.ext[:line_to] = @NLconstraint(
        opfmodel, [???=1:nlines],
        Vm[t[???]]^2 * (
            ytf_abs[???] * Vm[f[???]]^2 + ytt_abs[???] * Vm[t[???]]^2 +
            2.0 * Vm[f[???]] * Vm[t[???]] * (yre_to[???] * cos(Va[f[???]]-Va[t[???]]) - yim_to[???] * sin(Va[f[???]]-Va[t[???]]))
        ) <= flow_max[???]
    )

    # Objective
    @objective(
        opfmodel,
        Min,
        sum(
            cost_coefs[g, 4] * Pg[g]^2 + cost_coefs[g, 3] * Pg[g] + cost_coefs[g, 2]
            for g in 1:ngen
        )
    )

    return opfmodel
end

function main(datafile::String,w)
    polar = ExaPF.PolarForm(datafile)
    stack = ExaPF.NetworkStack(polar)
    m = build_opf_model(polar, stack, Ipopt.Optimizer; line_constraints=true,w)
    JuMP.optimize!(m)
    return m
end
w_0 = zeros(18)
w_0[1] = 1
model = main("/Users/leonard/opf/case9.m",w_0)

# Get primal solution
pg = JuMP.value.(model[:Pg])
print(pg)
qg = JuMP.value.(model[:Qg])

# Get dual solution
## Power flow equations
yp = JuMP.dual.(model.ext[:active_pf])
yq = JuMP.dual.(model.ext[:reactive_pf])
