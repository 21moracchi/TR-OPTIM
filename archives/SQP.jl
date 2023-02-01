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
nvar = 6
npar = 4
ncon = 8
pload = [2,1]
qload = [0.1,0.1]
w_0 = [pload; qload]
alpha = 0.5
#print("initial consuption : ",w_0, "\n")
opf = Model(Ipopt.Optimizer)
pg0 = [alpha, -alpha]

@variable(opf, Pg[i=1:ngen], start=pg0[i]),
@variable(opf, Qg[i=1:ngen]),
@variable(opf, Vm2),
@variable(opf, Va2)
var = [Pg; Qg; [Vm2, Va2]]

Pd, Qd = pload, qload
Vm1, Va1 = 1.0, 0.0

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
# set_silent(opf)
JuMP.optimize!(opf)
x_0 = JuMP.value.(var)

lambda_0 = -vcat(JuMP.dual.(constraints)...)

t_0 = vcat([x_0,w_0]...)

function f(t)
    return 1/2 * t[1]^2 + alpha * t[2]
end
function c(t)
    return [t[1]-t[7] + t[5]*sin(t[6]),t[2]-t[8]-t[5]*sin(t[6]),t[3]-t[9]+t[5]*cos(t[6])-1,
    t[4]-t[10]+t[5]*cos(t[6])-t[5]^2,-t[4]+alpha*t[2]]
end

function g(t)
    return [t[5]-sqrt(alpha^2+1),t[6]-pi,-pi-t[6]]
end

Lag(t) = f(t) +  transpose(lambda_0)*vcat([c(t),g(t)]...)


grad_c = FiniteDiff.finite_difference_jacobian(c,t_0)

L = FiniteDiff.finite_difference_hessian(Lag,t_0)

∇L = FiniteDiff.finite_difference_jacobian(Lag,t_0)

#print(
L_xx = L[1:nvar,1:nvar]

L_xw = L[nvar+1:nvar+npar,1:nvar]

c_x = grad_c[:,1:6]
c_w = grad_c[:,7:10]

grad_g = FiniteDiff.finite_difference_jacobian(g,t_0)

g_x = grad_g[:,1:6]
g_w = grad_g[:,7:10]


dw = Float64[0,1,0,0]
#print("variation of the comsuption : ", dw, "\n\n")

sqp = Model(Ipopt.Optimizer)
@variable(sqp, dx[i=1:6], start = 0)
@objective(sqp, Min, 1/2 * transpose(dx) * L_xx * dx + transpose(dw) * L_xw * dx)

constraint_eq = @constraint(sqp, c_x * dx + c_w * dw .== 0 )
# TODO: a corriger
constraint_in = @constraint(sqp, g(t_0) + (g_x) * dx + (g_w) * dw .<= 0)

JuMP.optimize!(sqp)

Δs = [JuMP.value.(dx); JuMP.dual.(constraint_eq)]
println(Δs)

