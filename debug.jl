include("/Users/leonard/opf/SQP_9.jl")

include("/Users/leonard/opf/f_OPF_9.jl")

#file = "/Users/leonard/opf/case9.m"
file = "/Users/leonard/opf/case14.m"
#file = "/Users/leonard/opf/pglib_opf_case14_ieee__api.m"
#file = "/Users/leonard/opf/pglib_opf_case14_ieee.m"
function get_error(;w=w_0)
    J = solve_opf(w=w)
    J_ref =FiniteDiff.finite_difference_jacobian(f_gen, w)
    return (norm((J-J_ref)./max.(1,J_ref),Inf))
end

w_limit = [0.0, 0.42706926, 1.8537368699999999, 0.4785258, 0.14956434000000002, 0.22044222, 0.0, 0.0, 0.58053789, 0.17709459, 0.06887568, 0.12003189, 0.26569193999999996, 0.29322219, 0.0, 0.1271397, 0.190209, -0.0390429, 0.0160176, 0.0750825, 0.0, 0.0, 0.1661826, 0.0580638, 0.018019800000000002, 0.0160176, 0.0580638, 0.050055]


function plot_error(; w=w_0, data=file)
    T = range(0.50, 0.99, length=50)
    W = []
    for t in T
        push!(W, t * w_limit)
    end
    L = []
    for w in W
        push!(L,get_error(w=w))
    end
    
    plot(T, L, yaxis=:log)
end
