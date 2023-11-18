#=
    Use the branch and bound algorithm to solve mixed-integer programing problems (small or medium size)

    The standard format is defined as follows:
    
        \min_{x} c'*x
    s.t. A*x <= b
         lb <= x <= ub
         x_{i} \in J

=#

include("MixedIntegerLinearProgramming.jl")






lb = [0.0, 0.0, 0.0]
ub = [1.0, 1.0, 1.0]
cobj = [1.0, 1.0, 2.0]
vtype = [GRB_BINARY, GRB_BINARY, GRB_BINARY]
A = [1.0 2.0 3.0; -1.0 -1.0 0.0]
b = [4.0; -1.0]





result = @time mixed_integer_linear_programming(cobj, A, b, lb, ub, vtype, "max")
print(result)