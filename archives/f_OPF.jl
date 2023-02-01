#=
    Copyright (c) Léonard Moracchini 2023
=#

using ExaPF
using JuMP
using Ipopt
using LinearAlgebra
using ForwardDiff
using FiniteDiff
const PS = ExaPF.PowerSystem

nbus = 2
ngen = 2
nlines = 1

pload = [2,1]
qload = [0.1,0.1]
w_0 = [pload;qload]

function f_OPF(;w = w_0)
    alpha = 0.5

    opf = Model(Ipopt.Optimizer)
    pg0 = [alpha, -alpha]
    @variable(opf, Pg[i=1:ngen], start=pg0[i])
    @variable(opf, Qg[i=1:ngen])
    @variable(opf, Vm2)
    @variable(opf, Va2)
    var = [Pg; Qg; [Vm2, Va2]]

    Pd, Qd = w[1:2], w[3:4]

    constraints = [
                   @NLconstraint(opf, Pg[1] - Pd[1]+ Vm2*sin(Va2) == 0 ),
                   @NLconstraint(opf, Pg[2] - Pd[2] - Vm2*sin(Va2) == 0 ),
                   @NLconstraint(opf, Qg[1] - Qd[1 ]+ Vm2*cos(Va2)-1 == 0),
                   @NLconstraint(opf, Qg[2] - Qd[2] + Vm2*cos(Va2)-Vm2^2 == 0),
                   @NLconstraint(opf, alpha * Pg[2] ==  Qg[2]),
                   @NLconstraint(opf, Vm2 <= sqrt(alpha^2 +1)),
                   @NLconstraint(opf,Va2<= pi),
                   @NLconstraint(opf,-pi<= Va2)
                  ]
    @objective(opf, Min, 1/2 * Pg[1]^2 + alpha*Pg[2])
    set_silent(opf)
    JuMP.optimize!(opf)
    x_0 = JuMP.value.(var)

    lambda_0 = vcat(JuMP.dual.(constraints)...)

    return vcat([x_0,lambda_0]...)
end



