
using ExaPF
using JuMP
using Ipopt
using LinearAlgebra
using ForwardDiff
using FiniteDiff
using SparseArrays
using JuMP
using Plots
const PS = ExaPF.PowerSystem


constraint_index(cons::Vector{NonlinearConstraintRef{ScalarShape}}) = getfield.(JuMP.index.(cons), :value)
polar = ExaPF.PolarForm("/Users/leonard/opf/case9.m")
buffer = ExaPF.NetworkStack(polar)
w_0 = vcat([buffer.pload,buffer.qload]...)
NBUS = ExaPF.get(polar, PS.NumberOfBuses())
NGEN = ExaPF.get(polar, PS.NumberOfGenerators())
NLINES = ExaPF.get(polar, PS.NumberOfLines())

function get_opf(;w = w_0,data = "/Users/leonard/opf/case9.m")
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
    yff_abs = yff_re.^2 .+ yff_im.^2
    yft_abs = yft_re.^2 .+ yft_im.^2
    yre_fr =   yff_re .* yft_re .+ yff_im .* yft_im
    yim_fr = - yff_re .* yft_im .+ yff_im .* yft_re

    ### to lines
    ytf_abs = ytf_re.^2 .+ ytf_im.^2
    ytt_abs = ytt_re.^2 .+ ytt_im.^2
    yre_to =   ytf_re .* ytt_re .+ ytf_im .* ytt_im
    yim_to = - ytf_re .* ytt_im .+ ytf_im .* ytt_re

    ## Objective
    cost_coefs = PS.get_costs_coefficients(pf)
    bus2gen = PS.get_bus_generators(pf.buses, pf.generators, pf.bus_to_indexes)
        
    #=
        Build model
    =#

    opfmodel = Model(solver)

    # VARIABLES
    @variable(opfmodel,  Pg[i=1:ngen] , start=pg0[i])
    @variable(opfmodel,  Qg[i=1:ngen] )
    @variable(opfmodel, Vm[i=1:nbus] , start=vm0[i])
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
        ) - (sum(Pg[g] for g in get(bus2gen, b, Int[])) - Pd[b]) == 0 )
    ## reactive
    opfmodel.ext[:reactive_pf] = @NLconstraint(
        opfmodel, [b=1:nbus],
        Vm[b] * sum(
            Vm[rows[c]] * (g_ij[c] * sin(Va[b] - Va[rows[c]]) - b_ij[c] * cos(Va[b] - Va[rows[c]]))
            for c in (Ybus.colptr[b]):(Ybus.colptr[b+1]-1)) - (sum(Qg[g] for g in get(bus2gen, b, Int[])) - Qd[b]) ==0
    )

    opfmodel.ext[:line_fr] = @NLconstraint(
        opfmodel, [ℓ=1:nlines],
        Vm[f[ℓ]]^2 * (
            yff_abs[ℓ] * Vm[f[ℓ]]^2 + yft_abs[ℓ] * Vm[t[ℓ]]^2 +
            2.0 * Vm[f[ℓ]] * Vm[t[ℓ]] * (yre_fr[ℓ] * cos(Va[f[ℓ]]-Va[t[ℓ]]) - yim_fr[ℓ] * sin(Va[f[ℓ]]-Va[t[ℓ]]))
        ) - flow_max[ℓ] <= 0
    )

    opfmodel.ext[:line_to] = @NLconstraint(
        opfmodel, [ℓ=1:nlines],
        Vm[t[ℓ]]^2 * (
            ytf_abs[ℓ] * Vm[f[ℓ]]^2 + ytt_abs[ℓ] * Vm[t[ℓ]]^2 +
            2.0 * Vm[f[ℓ]] * Vm[t[ℓ]] * (yre_to[ℓ] * cos(Va[f[ℓ]]-Va[t[ℓ]]) - yim_to[ℓ] * sin(Va[f[ℓ]]-Va[t[ℓ]]))
        ) - flow_max[ℓ] <= 0
    )

    ##power and voltage limitation
    constraints_limit = vcat([@NLconstraint(opfmodel, [b=1:ngen],Pg[b]<=pg_max[b]),
    @NLconstraint(opfmodel, [b=1:ngen],-Pg[b] + pg_min[b]<= 0),
    @NLconstraint(opfmodel, [b=1:ngen],Qg[b]-qg_max[b]<= 0),
    @NLconstraint(opfmodel, [b=1:ngen],-Qg[b]+qg_min[b]<= 0),
    @NLconstraint(opfmodel, [b=1:nbus],Vm[b]-vm_max[b]<= 0),
    @NLconstraint(opfmodel, [b=1:nbus],-Vm[b]+ vm_min[b]<= 0) ]...)

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


function solve_opf(;w = w_0,data = "/Users/leonard/opf/case9.m",get_sensitivity_jacobian = true, get_termination_status = false,
    get_constraints_jacobian = false)

    #get the opf model
    opfmodel = get_opf( w= w, data = data)
    polar = ExaPF.PolarForm(data)
    buffer = ExaPF.NetworkStack(polar)

    nbus = ExaPF.get(polar, PS.NumberOfBuses())
    ngen = ExaPF.get(polar, PS.NumberOfGenerators())
    nlines = ExaPF.get(polar, PS.NumberOfLines())
    #Optimize
    set_silent(opfmodel)
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

        return [x;y]

    else 
        nx = 2ngen + 2nbus #dimension of variables
        nw = 2nbus #dimension of parameters

        dim_c = 2nbus
        dim_g = 2nlines + 2nbus + 4ngen

    
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
        if get_constraints_jacobian == true
            return J
        end
        #Partial jacobians
        c_x = J[1:dim_c,1:nx]
        c_x = convert(Matrix{Float64},c_x)
        c_w = J[1:dim_c , nx + 1 : nx + nw]
        c_w = convert(Matrix{Float64},c_w)
        g_x = J[dim_c+1:dim_c+dim_g,1:nx]
        g_x = convert(Matrix{Float64},g_x)
        g_w = J[dim_c+1: dim_c + dim_g , nx + 1 : nx + nw]
        g_w = convert(Matrix{Float64},g_w)

        #Hessian
        hess_struct = MOI.hessian_lagrangian_structure(ipopt_model)
        i_hess = [j[1] for j in hess_struct]
        j_hess = [j[2] for j in hess_struct]
        nnzh = length(i_hess)
        v_hess = zeros(nnzh)
        σ = 1.0 # scaling of objective
        MOI.eval_hessian_lagrangian(ipopt_model, v_hess, x, σ, -y)

        H = Matrix(sparse(i_hess, j_hess, v_hess))
        for i in 1:2nbus + 2ngen
            for j in i+1:2nbus + 2ngen
                H[i,j] = H[j,i]
            end
        end
        
        #=
            SQP Computation
        =# 

        

        #QP 

        function sensitivity(variation) #returns the variation of the primal-dual solution in response to a variation of the parameters. 
            sqp = Model(Ipopt.Optimizer)

            @variable(sqp, dx[i=1:nx])

            @objective(sqp, Min,  1/2 * transpose(dx) * H * dx ) #pas de L_xw car rien ne dépend des deux. 

            constraint_eq = @constraint(sqp, c_x * dx + c_w * variation  .== 0 )
            constraint_in = @constraint(sqp, cons[dim_c+1:dim_g+dim_c] + (g_x * dx) + (g_w * variation)  .<= 0 )

            set_silent(sqp)
            JuMP.optimize!(sqp)
            
            ds = [JuMP.value.(dx);JuMP.dual.(constraint_eq)]
            return -JuMP.dual.(constraint_eq)
            return ds
        end


        function jacobian_computation() #returns the jacobian with respect to the parameters
            dwi = zeros(nw)
            dwi[1] = 0.01
            jacobian = sensitivity(dwi)
            for i in 2:nw
                var_coeff_i = zeros(nw)
                var_coeff_i[i] = 0.01
                jacobian = hcat(jacobian,sensitivity(var_coeff_i))

            end
            return jacobian
        end
            
        return 100*jacobian_computation()
    end


end
function inf_point(;starting_point = w_0, direction = w_0, data = "/Users/leonard/opf/case9.m")
    alpha = 0
    termination_status = MOI.LOCALLY_SOLVED
    while termination_status == MOI.LOCALLY_SOLVED 
        alpha += 0.001
        termination_status = solve_opf(;w = starting_point + alpha * direction,get_sensitivity_jacobian = false, get_termination_status = true,
        get_constraints_jacobian = false)
    end 
    return starting_point + alpha_limit * direction
end


w_limit = [0.0, 0.0, 0.0, 0.0, 1.9632599999998979, 0.0, 2.1813999999998863, 
0.0, 2.726749999999858, 0.0, 0.0, 0.0, 0.0, 0.6544199999999658, 0.0, 0.7634899999999601, 0.0, 1.0906999999999432]

function get_active_constraints(;w = w_0, data = "/Users/leonard/opf/case9.m")
    opf_model = get_opf(w = w, data = data )
    set_silent(opf_model)
    optimize!(opf_model)
    #print(termination_status(opf_model))
    constraint_index = 0
    active_set = []
    active_constraints_dict = primal_feasibility_report(opf_model)
    for constraint in all_nonlinear_constraints(opf_model)
        if haskey(active_constraints_dict, constraint)
            push!(active_set, constraint_index)
        end
        constraint_index += 1
    end
    return active_set .+ 1
    
end

function print_active_constraints(;w = w_0, data = "/Users/leonard/opf/case9.m")
    opf_model = get_opf(w = w, data = data )
    set_silent(opf_model)
    optimize!(opf_model)
    return primal_feasibility_report(opf_model)
end

function rank_jacobian(;w = w_0, data = "/Users/leonard/opf/case9.m")
    J = solve_opf(w = w, data = data, get_constraints_jacobian = true)
    J = J[:,1:2NGEN + 2NBUS] #uniquement le gradient par rapport à x et pas w
    active_set = get_active_constraints(w = w, data = data)
    J_active = reshape([],0,24)
    for i in active_set
        J_active = vcat(J_active, transpose(J[i,:]))
    end
    J_active = convert(Matrix{Float64},J_active)
    
    return rank(J_active,atol = 0.0000001)/size(active_set)[1]
end

function plot_price(;w=w_0,data = "/Users/leonard/opf/case9.m")
    T = range(0.97,1.01,length = 1000)
    W = []
    for t in T
        push!(W, t*w_limit+(1-t)*w)
    end
    L = reshape([],0,NBUS)
    for w in W
        solution = -solve_opf(w=w, get_sensitivity_jacobian = false)

        L = vcat(L, solution[4NBUS + 2NGEN + 1: 5NBUS + 2NGEN]') # on récupère les lambda des équations powerflow actives
        
    end
    L = convert(Matrix{Float64},L)
    plot(T,L)
end

function plot_jacobian(;w=w_0,data = "/Users/leonard/opf/case9.m")
    T = range(0.95,1.1,length = 100)
    W = []
    for t in T
        push!(W, t*w_limit+(1-t)*w)
    end
 
    R = []
    for w in W 
        push!(R, rank_jacobian(w=w))
    end
    plot(T,R)
end


function SVD_jacobian(;w=w_0,data = "/Users/leonard/opf/case9.m")
    J = solve_opf(w = w, data = data, get_constraints_jacobian = true)
    J = J[:,1:2NGEN + 2NBUS] #uniquement le gradient par rapport à x et pas w
    active_set = get_active_constraints(w = w, data = data)
    J_active = reshape([],0,24)
    for i in active_set
        J_active = vcat(J_active, transpose(J[i,:]))
    end
    J_active = convert(Matrix{Float64},J_active)
    s_values = svd(J_active).S
    return s_values
end



function error(A,B)
    return (A.-B)./max.(1,A)
end

