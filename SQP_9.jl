
using ExaPF
using JuMP
using Ipopt
using LinearAlgebra
using ForwardDiff
using FiniteDiff
using SparseArrays
using JuMP
using Plots
using HiGHS
using Polyhedra
using CDDLib
using LaTeXStrings
const PS = ExaPF.PowerSystem



constraint_index(cons::Vector{NonlinearConstraintRef{ScalarShape}}) = getfield.(JuMP.index.(cons), :value)
#file = "/Users/leonard/opf/case9.m"
#file = "/Users/leonard/opf/pglib_opf_case14_ieee__api.m"
#file = "/Users/leonard/opf/pglib_opf_case14_ieee.m"
#file = "/Users/leonard/opf/case14.m"

file = "/Users/leonard/opf/case5_critical_opf.m"
polar = ExaPF.PolarForm(file)
buffer = ExaPF.NetworkStack(polar)
w_0 = vcat([buffer.pload, buffer.qload]...)
load_factor = 1.139332 # load_factor critique pour le case5
w_1 = w_0 * load_factor
NBUS = ExaPF.get(polar, PS.NumberOfBuses())
NGEN = ExaPF.get(polar, PS.NumberOfGenerators())
NLINES = ExaPF.get(polar, PS.NumberOfLines())
N_C = 2NBUS
N_X = 2NBUS + 2NGEN
N_W = 2NBUS

x_c = [2.615304394048055, 9.854626876142665e-09, 1.8000000068306052, 1.5089721452039435, 1.5107513979100616e-09, 1.13840668720261, 2.4585048970061872e-08, 1.500000022913465, 0.5209280611395186, 3.852929057151792e-09, 1.0500003343970223, 0.9499996286190487, 1.049656489901429, 1.050000101176089, 0.9962415371635135, -1.169383592929008, -1.3094485371191, -1.2342555657720744, -1.1778664571949722, -1.2513778588833584, 0.92285892, 1.62924476, 1.07097208, 0.99121884, 1.20769192, 0.42155284, 0.854499, 0.58105932, 0.42155284, 0.569666]

y_c = [-543060.8796857579, -11815476986.996782, -4264872167.3605523, -915383.2870773544, -1689979463.1183052, -3.4281061260832822e-06, -29733120430.007057, -14478093003.212645, -6.9543424289791724e-06, -4553078798.200016, -2.4707197651038695e-06, -2.0859999287756096e-06, -1.8482376189047735e-06, -1.8072002350862042e-06, -1.8021910100100372e-06, -2.2608864531466414e-06, -2.3338951360472733e-06, -2.140430167460519e-06, -1.8484164712018849e-06, -1.8070563582008217e-06, -1.7975511410092602e-06, -2.3169431234289117e-06, -0.000873593136274457, -11815477956.910263, -4264112167.3578258, -0.0, -1689985455.9990644, -2.895163730804003e-06, -969.9134817346932, -5.363986783935588e-06, -8.778721722862993e-06, -5992.880759091349, -7.72304883266667e-06, -29733120821.141308, -14478093003.212648, -1.0816809586791887e-05, -4553081240.655457, -3.500256798366805e-06, -391.1342527717352, -3.3293710936689615e-06, -4.778152498352236e-06, -2442.455440701451, -9.655143952090981e-05, -237089573285.09555, -9.51703922041367e-05, -9.655019920743077e-05, -0.00020767934952760074, -212782254491.90118, -9.655107154099313e-05, -0.17238132623176805, -64349713300.20704, -0.00018043116485910295]


function get_opf(; w=w_0, data=file)
    polar = ExaPF.PolarForm(data)
    buffer = ExaPF.NetworkStack(polar)

    nbus = ExaPF.get(polar, PS.NumberOfBuses())
    ngen = ExaPF.get(polar, PS.NumberOfGenerators())
    nlines = ExaPF.get(polar, PS.NumberOfLines())
    solver = Ipopt.Optimizer
    line_constraints = true
    pf = polar.network
    baseMVA = pf.baseMVA

    # Initial consuption (parameter)
    pload = w[1:nbus]
    qload = w[nbus+1:2nbus]

    # DATA from .m file
    ## Bounds
    pg_min, pg_max = PS.bounds(pf, PS.Generators(), PS.ActivePower())
    qg_min, qg_max = PS.bounds(pf, PS.Generators(), PS.ReactivePower())
    vm_min, vm_max = PS.bounds(pf, PS.Buses(), PS.VoltageMagnitude())

    ## Max flow through lines
    flow_min, flow_max = PS.bounds(pf, PS.Lines(), PS.ActivePower())
    flow_max = min.(1e5, flow_max)

    ## Initial values
    vm0 = buffer.vmag
    va0 = buffer.vang
    pg0 = buffer.pgen

    ## Bus-admittance matrix
    Ybus = pf.Ybus
    rows = Ybus.rowval
    yvals = Ybus.nzval
    g_ij = real.(yvals)
    b_ij = imag.(yvals)

    ## Line characteristics (Matpower manual, Section 3.2, page 25)
    yff_re = real.(pf.lines.Yff)
    yff_im = imag.(pf.lines.Yff)
    yft_re = real.(pf.lines.Yft)
    yft_im = imag.(pf.lines.Yft)
    ytf_re = real.(pf.lines.Ytf)
    ytf_im = imag.(pf.lines.Ytf)
    ytt_re = real.(pf.lines.Ytt)
    ytt_im = imag.(pf.lines.Ytt)

    ## Line constraints
    f = pf.lines.from_buses
    t = pf.lines.to_buses

    ### from lines
    yff_abs = yff_re .^ 2 .+ yff_im .^ 2
    yft_abs = yft_re .^ 2 .+ yft_im .^ 2
    yre_fr = yff_re .* yft_re .+ yff_im .* yft_im
    yim_fr = -yff_re .* yft_im .+ yff_im .* yft_re

    ### to lines
    ytf_abs = ytf_re .^ 2 .+ ytf_im .^ 2
    ytt_abs = ytt_re .^ 2 .+ ytt_im .^ 2
    yre_to = ytf_re .* ytt_re .+ ytf_im .* ytt_im
    yim_to = -ytf_re .* ytt_im .+ ytf_im .* ytt_re

    ## Objective
    cost_coefs = PS.get_costs_coefficients(pf)
    bus2gen = PS.get_bus_generators(pf.buses, pf.generators, pf.bus_to_indexes)

    #=
        Build model
    =#

    opfmodel = Model(solver)
    # VARIABLES
    @variable(opfmodel, Pg[i=1:ngen], start = pg0[i])
    @variable(opfmodel, Qg[i=1:ngen])
    @variable(opfmodel, Vm[i=1:nbus], start = vm0[i])
    @variable(opfmodel, Va[i=1:nbus], start = va0[i])
    @variable(opfmodel, Pd[i=1:nbus])
    @variable(opfmodel, Qd[i=1:nbus])
    JuMP.fix.(Pd, pload)
    JuMP.fix.(Qd, qload)


    # Power-flow constraints

    ## active
    opfmodel.ext[:active_pf] = @NLconstraint(
        opfmodel, [b = 1:nbus],
        Vm[b] * sum(
            Vm[rows[c]] * (g_ij[c] * cos(Va[b] - Va[rows[c]]) + b_ij[c] * sin(Va[b] - Va[rows[c]]))
            for c in (Ybus.colptr[b]):(Ybus.colptr[b+1]-1)
        ) - (sum(Pg[g] for g in get(bus2gen, b, Int[])) - Pd[b]) == 0)
    ## reactive
    opfmodel.ext[:reactive_pf] = @NLconstraint(
        opfmodel, [b = 1:nbus],
        Vm[b] * sum(
            Vm[rows[c]] * (g_ij[c] * sin(Va[b] - Va[rows[c]]) - b_ij[c] * cos(Va[b] - Va[rows[c]]))
            for c in (Ybus.colptr[b]):(Ybus.colptr[b+1]-1)) - (sum(Qg[g] for g in get(bus2gen, b, Int[])) - Qd[b]) == 0
    )

    opfmodel.ext[:line_fr] = @NLconstraint(
        opfmodel, [ℓ = 1:nlines],
        Vm[f[ℓ]]^2 * (
            yff_abs[ℓ] * Vm[f[ℓ]]^2 + yft_abs[ℓ] * Vm[t[ℓ]]^2 +
            2.0 * Vm[f[ℓ]] * Vm[t[ℓ]] * (yre_fr[ℓ] * cos(Va[f[ℓ]] - Va[t[ℓ]]) - yim_fr[ℓ] * sin(Va[f[ℓ]] - Va[t[ℓ]]))
        ) - flow_max[ℓ] <= 0
    )

    opfmodel.ext[:line_to] = @NLconstraint(
        opfmodel, [ℓ = 1:nlines],
        Vm[t[ℓ]]^2 * (
            ytf_abs[ℓ] * Vm[f[ℓ]]^2 + ytt_abs[ℓ] * Vm[t[ℓ]]^2 +
            2.0 * Vm[f[ℓ]] * Vm[t[ℓ]] * (yre_to[ℓ] * cos(Va[f[ℓ]] - Va[t[ℓ]]) - yim_to[ℓ] * sin(Va[f[ℓ]] - Va[t[ℓ]]))
        ) - flow_max[ℓ] <= 0
    )

    ###power and voltage limitation
    constraints_limit = vcat([
        @NLconstraint(opfmodel, [b = 1:ngen], Pg[b] <= pg_max[b]),
        @NLconstraint(opfmodel, [b = 1:ngen], -Pg[b] + pg_min[b] <= 0),
        @NLconstraint(opfmodel, [b = 1:ngen], Qg[b] - qg_max[b] <= 0),
        @NLconstraint(opfmodel, [b = 1:ngen], -Qg[b] + qg_min[b] <= 0),
        @NLconstraint(opfmodel, [b = 1:nbus], -Vm[b] + vm_min[b] <= 0),
        @NLconstraint(opfmodel, [b = 1:nbus], Vm[b] - vm_max[b] <= 0),
    ]...)

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


function solve_opf(; w=w_0, data=file, get_sensitivity_jacobian=true, get_termination_status=false,
    get_constraints_jacobian=false)

    #get the opf model
    opfmodel = get_opf(w=w, data=data)
    JuMP.set_optimizer_attribute(opfmodel, "tol", 1e-9)

    polar = ExaPF.PolarForm(data)
    buffer = ExaPF.NetworkStack(polar)

    nbus = ExaPF.get(polar, PS.NumberOfBuses())
    ngen = ExaPF.get(polar, PS.NumberOfGenerators())
    nlines = ExaPF.get(polar, PS.NumberOfLines())
    #Optimize
    set_silent(opfmodel)
    JuMP.set_optimizer_attribute(opfmodel, "print_level", 0)
    JuMP.set_optimizer_attribute(opfmodel, "nlp_scaling_method", "none")
    JuMP.optimize!(opfmodel)


    # Get primal-dual solution
    pg = JuMP.value.(opfmodel[:Pg])
    qg = JuMP.value.(opfmodel[:Qg])
    vm = JuMP.value.(opfmodel[:Vm])
    va = JuMP.value.(opfmodel[:Va])


    nvar = JuMP.num_variables(opfmodel)
    ncon = JuMP.num_constraints(opfmodel; count_variable_in_set_constraints=false)


    # Get MOI model sent to Ipopt
    moi_model = JuMP.backend(opfmodel)
    ipopt_model = moi_model.optimizer.model

    # Current primal solution
    x = Float64[MOI.get(moi_model, MOI.VariablePrimal(), MOI.VariableIndex(vi)) for vi in 1:nvar]

    # Current dual solution
    y = MOI.get(ipopt_model, MOI.NLPBlockDual())
    if get_termination_status == true
        return termination_status(opfmodel)
    end

    if get_sensitivity_jacobian == false

        return [x; y]

    else
        nx = 2ngen + 2nbus #dimension of variables
        nw = 2nbus #dimension of parameters

        dim_c = 2nbus
        dim_g = ncon - dim_c

        # Constraints
        cons = zeros(ncon)
        MOI.eval_constraint(ipopt_model, cons, x)

        # Gradient
        grad = zeros(nvar)
        MOI.eval_objective_gradient(ipopt_model, grad, x)

        # Jacobian
        jac_struct = MOI.jacobian_structure(ipopt_model)
        i_jac = [j[1] for j in jac_struct]
        j_jac = [j[2] for j in jac_struct]
        nnzj = length(i_jac)
        v_jac = zeros(nnzj)
        MOI.eval_constraint_jacobian(ipopt_model, v_jac, x)
        J = Matrix(sparse(i_jac, j_jac, v_jac))

        ∇L = grad .- J' * y
        #@assert norm(∇L[1:nx], Inf) <= 1e-8

        if get_constraints_jacobian == true
            return J
        end
        #Partial jacobians
        c_x = J[1:dim_c, 1:nx]
        c_x = convert(Matrix{Float64}, c_x)
        c_w = J[1:dim_c, nx+1:nx+nw]
        c_w = convert(Matrix{Float64}, c_w)
        g_x = J[dim_c+1:dim_c+dim_g, 1:nx]
        g_x = convert(Matrix{Float64}, g_x)
        g_w = J[dim_c+1:dim_c+dim_g, nx+1:nx+nw]
        g_w = convert(Matrix{Float64}, g_w)

        g_k = cons[dim_c+1:dim_g+dim_c]
        #Hessian
        hess_struct = MOI.hessian_lagrangian_structure(ipopt_model)
        i_hess = [j[1] for j in hess_struct]
        j_hess = [j[2] for j in hess_struct]
        nnzh = length(i_hess)
        v_hess = zeros(nnzh)
        σ = 1.0 # scaling of objective
        MOI.eval_hessian_lagrangian(ipopt_model, v_hess, x, σ, -y)

        H = Matrix(sparse(i_hess, j_hess, v_hess))
        for i in 1:nx
            for j in i+1:nx
                H[i, j] = H[j, i]
            end
        end

        #QP

        function sensitivity(variation) #returns the variation of the primal-dual solution in response to a variation of the parameters.
            sqp = Model(Ipopt.Optimizer)
            JuMP.set_optimizer_attribute(sqp, "tol", 1e-9)

            @variable(sqp, dx[i=1:nx])

            @objective(sqp, Min, 1 / 2 * transpose(dx) * H * dx) #pas de L_xw car rien ne dépend des deux.

            constraint_eq = @constraint(sqp, c_x * dx + c_w * variation .== 0)
            for i in 1:dim_g
                gi = cons[i+dim_c]

                is_active = gi > -1e-7
                is_mult_null = abs(y[i+dim_c]) < 1e-7
                if is_active & !is_mult_null
                    @constraint(sqp, dot(g_x[i, :], dx) + dot(g_w[i, :], variation) == 0.0)
                else
                    
                    # TODO: commented right now, as it returns the wrong result.
                    # Ipopt is not converging if we include the inequality constraints.
                    # We should check whether the polyhedra is nonempty
                    @constraint(sqp, gi + dot(g_x[i, :], dx) + dot(g_w[i, :], variation) <= 0.0)
                end
            end

            set_silent(sqp)
            JuMP.optimize!(sqp)
            print(termination_status(sqp), "\n")
            ds = [JuMP.value.(dx)[1:2NGEN+NBUS]; JuMP.dual.(constraint_eq)]
            #return -JuMP.dual.(constraint_eq)
            return ds
        end

        function jacobian_computation() #returns the jacobian with respect to the parameters
            dwi = zeros(nw)
            dwi[1] = 1e-5
            jacobian = sensitivity(dwi)
            for i in 2:nw
                var_coeff_i = zeros(nw)
                var_coeff_i[i] = 1e-5
                jacobian = hcat(jacobian, sensitivity(var_coeff_i))

            end
            return 1e5* jacobian
        end

        return jacobian_computation()
    end


end
function inf_point(; starting_point=w_0, direction=w_0, data=file)
    alpha = 0
    termination_status = MOI.LOCALLY_SOLVED
    while termination_status == MOI.LOCALLY_SOLVED
        alpha += 0.0001
        termination_status = solve_opf(; w=starting_point + alpha * direction, get_sensitivity_jacobian=false, get_termination_status=true,
            get_constraints_jacobian=false)
    end
    return alpha
end


#w_limit = [0.0, 0.0, 0.0, 0.0, 1.9632599999998979, 0.0, 2.1813999999998863, 0.0, 2.726749999999858, 0.0, 0.0, 0.0, 0.0, 0.6544199999999658, 0.0, 0.7634899999999601, 0.0, 1.0906999999999432]
#case14_api : 
#w_limit = [0.92285892, 1.62924476, 1.07097208, 0.99121884, 1.20769192, 0.42155284, 0.854499, 0.58105932, 0.42155284, 0.569666]
w_limit = w_1
function get_active_constraints(; w=w_0, data=file)
    opf_model = get_opf(w=w, data=data)
    set_silent(opf_model)
    optimize!(opf_model)
    #print(termination_status(opf_model))
    constraint_index = 1
    active_set = []
    active_constraints_dict = primal_feasibility_report(opf_model)
    for constraint in all_nonlinear_constraints(opf_model)
        if haskey(active_constraints_dict, constraint)
            push!(active_set, constraint_index)

        end
        constraint_index += 1
    end
    return active_set

end

function rank_jacobian(; w=w_0, data=file)
    J = solve_opf(w=w, data=data, get_constraints_jacobian=true)
    J = J[:, 1:2NGEN+2NBUS] #uniquement le gradient par rapport à x et pas w
    active_set = get_active_constraints(w=w, data=data)
    J_active = reshape([], 0, 2NGEN + 2NBUS)
    for i in active_set
        J_active = vcat(J_active, transpose(J[i, :]))
    end
    J_active = convert(Matrix{Float64}, J_active)

    return rank(J_active, atol=10e-4) / size(active_set)[1]
end

function plot_price(; w=w_0, data=file)
    T = range(0.95, 1.03, length=100)
    W = []
    for t in T
        push!(W, t * w_limit + (1 - t) * w)
    end
    L = reshape([], 0, NBUS)
    for w in W
        solution = -solve_opf(w=w, get_sensitivity_jacobian=false)

        L = vcat(L, solution[4NBUS+2NGEN+1:5NBUS+2NGEN]') # on récupère les lambda des équations powerflow actives

    end
    L = convert(Matrix{Float64}, L)
    plot(T, abs.(L), yaxis=:log)
    
end

function plot_jacobian(; w=w_0, data=file)
    T = range(0.6, 1.05, length=40)
    W = []
    for t in T
        push!(W, t * w_limit + (1 - t) * w)
    end

    R = []
    for w in W
        push!(R, rank_jacobian(w=w, data=data))
    end
    plot(T, R, label = L"rank(tp_{limit})",lw = 1, dpi=300)
    xlabel!(L"t")
    xticks!(0.95:0.01:1.1)
    savefig("plot_jacobian.png")
    
end


function SVD_jacobian(; w=w_0, data=file)
    J = solve_opf(w=w, data=data, get_constraints_jacobian=true)
    J = J[:, 1:2NGEN+2NBUS] #uniquement le gradient par rapport à x et pas w
    active_set = get_active_constraints(w=w, data=data)
    J_active = reshape([], 0, 2NGEN + 2NBUS)
    for i in active_set
        J_active = vcat(J_active, transpose(J[i, :]))
    end
    J_active = convert(Matrix{Float64}, J_active)
    s_values = svd(J_active).S
    return s_values
end

function plot_min_SVD(; w=w_0, data=file)
    T = range(0.95, 1.05, length=10)
    W = []
    for t in T
        push!(W, t * w_limit + (1 - t) * w)
    end
    L = []
    for w in W
        push!(L, minimum(abs.(SVD_jacobian(w=w, data=data))))
    end

    plot(T, L, yaxis=:log,label = "Minimal value of the SVD",dpi=300)
    xlabel!(L"t")
    
    xticks!(0.95:0.01:1.05)

end

function double_plot(; w=w_0,data= file)
    T = range(0.95, 1.05, length=100)
    W = []
    for t in T
        push!(W, t * w_limit + (1 - t) * w)
    end

    R = []
    for w in W
        push!(R, rank_jacobian(w=w, data=data))
    end

    L = []
    for w in W
        push!(L, minimum(abs.(SVD_jacobian(w=w, data=data))))
    end
    
    
    plot(T,R, label = "Normalized rank", color = :red, legend = :topleft, ylims = (0.92,1.02),dpi = 300 )
    vline!([1,1], label = "Feasibility border",color = :black)
    p = twinx()
    plot!(p,T,L,yaxis =:log, label = "Minimal value of the SVD", color = :blue, legend = :topright,xlabel = L"t", dpi=300)
    
    savefig("double_plot.png")
end
#J = solve_opf(w=w_0)

function check_MFCQ(; w=w_0, data=file)
    grad = solve_opf(w=w, data=data, get_constraints_jacobian=true)
    opfmodel = get_opf(w=w, data=data)
    ncon = JuMP.num_constraints(opfmodel; count_variable_in_set_constraints=false)
    J = grad[:, 1:N_X]

    active_set = get_active_constraints(w=w, data=data)
    J_active = reshape([], 0, 2NGEN + 2NBUS)
    for i in active_set
        J_active = vcat(J_active, transpose(J[i, :]))
    end
    n_g_active = size(active_set)[1] - N_C
    J_active = convert(Matrix{Float64}, J_active)
    J_C = J_active[1:N_C, :]
    J_G_active = J_active[N_C+1:N_C+n_g_active, :]

    model_MFCQ = Model(HiGHS.Optimizer)
    @variable(model_MFCQ, u[1:N_X])
    for i in 1:N_C
        @constraint(model_MFCQ, transpose(J_C[i, :]) * u == 0)
    end
    for i in 1:n_g_active
        @constraint(model_MFCQ, transpose(J_G_active[i, :]) * u <= -1)
    end
    @objective(model_MFCQ, Min, 0)
    optimize!(model_MFCQ)
    if termination_status(model_MFCQ) == MOI.OPTIMAL
        return true
    else
        return false
    end
end


function check_SMFCQ(; w=w_0, data=file)
    grad = solve_opf(w=w, data=data, get_constraints_jacobian=true)
    opfmodel = get_opf(w=w, data=data)

    ncon = JuMP.num_constraints(opfmodel; count_variable_in_set_constraints=false)
    J = grad[:, 1:N_X]

    active_set = get_active_constraints(w=w, data=data)

    prim_dual = solve_opf(w=w, data=data, get_sensitivity_jacobian=false)
    z = prim_dual[N_X+N_W+1:N_X+N_W+ncon] #multiplicateurs associés aux contraintes égalité puis inégalité
    print(z)
    model_SMFCQ = Model(HiGHS.Optimizer)
    @variable(model_SMFCQ, u[1:N_X])

    for i in active_set
        if i <= N_C || abs(z[i]) > 1e-5 #si c'est une contrainte égalité ou une contrainte fortement active
            @constraint(model_SMFCQ, transpose(J[i, :]) * u == 0)

        else
            print("ok")
            @constraint(model_SMFCQ, transpose(J[i, :]) * u <= -1)
        end
    end

    @objective(model_SMFCQ, Min, 0)
    optimize!(model_SMFCQ)
    return JuMP.value.(u)
end


function get_primal_dual_solution(moi_model, ipopt_model, nvar, ncon)
    x = Float64[MOI.get(moi_model, MOI.VariablePrimal(), MOI.VariableIndex(vi)) for vi in 1:nvar]
    y =  MOI.get(ipopt_model, MOI.NLPBlockDual())
    y = -y #our convention
    return x,y
end

function get_constraint_jacobian(x,ipopt_model)
    jac_struct = MOI.jacobian_structure(ipopt_model)
    i_jac = [j[1] for j in jac_struct]
    j_jac = [j[2] for j in jac_struct]
    nnzj = length(i_jac)
    v_jac = zeros(nnzj)
    MOI.eval_constraint_jacobian(ipopt_model, v_jac, x)
    J = Matrix(sparse(i_jac, j_jac, v_jac))
    return J
end

function get_gradient(x,ipopt_model,nvar)
    grad = zeros(nvar)
    MOI.eval_objective_gradient(ipopt_model, grad, x)
    return grad
end

function get_constraints_value(x,ncon,ipopt_model)
    cons = zeros(ncon)
    MOI.eval_constraint(ipopt_model, cons, x)
    return cons
end

function polyhedra(; w=w_0, data=file)

    opfmodel = get_opf(w=w, data=data)
    set_silent(opfmodel)
    optimize!(opfmodel)
    
    nvar = JuMP.num_variables(opfmodel)
    ncon = JuMP.num_constraints(opfmodel; count_variable_in_set_constraints=false)

    # Get MOI model sent to Ipopt
    moi_model = JuMP.backend(opfmodel)
    ipopt_model = moi_model.optimizer.model
    #primal_dual_solution
    x,y = get_primal_dual_solution(moi_model, ipopt_model, nvar, ncon)
    #Gradient
    grad = get_gradient(x,ipopt_model,nvar)
    #Constraints
    cons = get_constraints_value(x,ncon,ipopt_model)
    n_g = ncon-N_C
    #Jacobian
    J = get_constraint_jacobian(x,ipopt_model)
    #Partial jacobian
    J_x = J[:, 1:N_X] 

    
    z = y[N_C+1:N_C+n_g] #distinguishing lagrangian multipliers
    y = y[1:N_C]
    
    #Find multipliers 
    lin_model = Model(HiGHS.Optimizer)
    @variable(lin_model, y[1:N_C])
    @variable(lin_model,z[1:n_g] >= 0)
    @constraint(lin_model, -1e-5.<= sum(J_x[i,:]*y[i] for i in 1:N_C)
     + sum(J_x[N_C+i,:]*z[i] for i in 1:n_g) + grad[1:N_X].<=1e-5)
    
    for i in 1:n_g
        @constraint(lin_model,-1e-6 <=cons[N_C+i]*z[i] <= 1e-6)
    end
    
    @objective(lin_model,Min,0)
    optimize!(lin_model)

    poly = polyhedron(lin_model, CDDLib.Library(:exact))
    vrep(poly)
    
end

function max_multiplier(;dw=1e-5*w_0/norm(w_0),w=w_0,data=file)
    opfmodel = get_opf(w=w, data=data)
    set_silent(opfmodel)
    optimize!(opfmodel)

    nvar = JuMP.num_variables(opfmodel)
    ncon = JuMP.num_constraints(opfmodel; count_variable_in_set_constraints=false)

    # Get MOI model sent to Ipopt
    moi_model = JuMP.backend(opfmodel)
    ipopt_model = moi_model.optimizer.model
    #primal_dual_solution (thanks to global variable)
    if w == w_1
        x,y = x_c,y_c
    else 
        x,y = get_primal_dual_solution(moi_model, ipopt_model, nvar, ncon)
    end
    #Gradient
    grad = get_gradient(x,ipopt_model,nvar)
    #Constraints
    cons = get_constraints_value(x,ncon,ipopt_model)
    n_g = ncon-N_C
    #Jacobian
    J = get_constraint_jacobian(x,ipopt_model)
    #Partial jacobian
    J_x = J[:, 1:N_X] 
    J_w = J[:,N_X+1:N_X+N_W]

    max_model = Model(HiGHS.Optimizer)
    @variable(max_model, y[1:N_C])
    @variable(max_model, z[1:n_g]>=0)
    @objective(max_model,Max, y'*J_w[1:N_C,:]*dw + z'*J_w[N_C+1:n_g+N_C,:]*dw)
    @constraint(max_model,  sum(J_x[i,:]*y[i] for i in 1:N_C)
     + sum(J_x[N_C+i,:]*z[i] for i in 1:n_g) + grad[1:N_X] .== 0)
    
    for i in 1:n_g
        @constraint(max_model,-1e-3 <=cons[N_C+i]*z[i] <= 1e-3)

    end
    optimize!(max_model)
    @assert(termination_status(max_model) == MOI.OPTIMAL)
    return JuMP.value.(-1e5*[y;z])
end



