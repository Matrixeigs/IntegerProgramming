#=
    Use the branch and bound algorithm to solve mixed-integer programing problems (small or medium size)

    The standard format is defined as follows:
    
        \max_{x} c'*x
    s.t. A*x <= b
         lb <= x <= ub
         x_{i} \in J

=#

include("LinearProgramming.jl")


function branch_and_bound(cobj::Vector, A::Matrix, b::Vector, lb::Vector, ub::Vector, vtype::Vector, model_sense::String)
    # Branch and bound function for the mixed integer programming problem
    # 1. solve the relaxed problem
    integerity_tol = Float16(1e-3)
    result0 = linear_programming(cobj, A, b, lb, ub, model_sense)

    root_node = Dict("bound" => result0["objval"], 
    "x" => result0["x"], 
    "depth"=>0,
    "feasibility" => sum(integerity_feasibility(result0["x"], vtype, integerity_tol)))
    # This algorithm 
    if root_node["feasibility"] <= integerity_tol
        x_incumbent = root_node["x"]
        f_opt = root_node["bound"]
    else
        x_incumbent = zeros(nx)
        f_opt = -Inf
    end
    # initial the queue
    queue = [root_node]
    while length(queue)>0 # The queue is not empty
        if f_opt > -Inf # A feasible solution is derived.
            (~,index) = queue_management(queue, "bound")
        else 
            (~,index) = queue_management(queue, "feasibility")
        end
    end
    print(root_node)


end

function integerity_feasibility(x::Vector, vtype::Vector, integerity_tol::Float16)
    # 
    nx = length(x)
    xfrac = zeros(nx)
    for i in 1:nx
        if vtype[i] != "C"
            # check the lower boundary
            x_floor = floor(x[i])
            if x_floor - integerity_tol <= x[i] <= x_floor + integerity_tol
                xfrac[i] = 0.0
            else
                xfrac[i] = max(x[i] - x_floor, xfrac[i])
            end
            # check the upper boundary
            x_ceil = ceil(x[i])
            if x_ceil - integerity_tol <= x[i] <= x_ceil + integerity_tol
                xfrac[i] = 0.0
            else
                xfrac[i] = max(x_floor - x[i], xfrac[i])
            end
        end
    end

    return xfrac
end

function queue_management(queue::Vector, key_to_check::String, sense::String)
    if sense == "max"
        max_index = 0
        max_value = nothing
    
        # Iterate through the vector to find the maximum index according to the key
        for (index, dict) in enumerate(queue)
            if haskey(dict, key_to_check)
                current_value = dict[key_to_check]
                if max_value === nothing || current_value >= max_value
                    max_value = current_value
                    max_index = index
                end
            end
        end

        return max_value, max_index
    else
        min_index = 0
        min_value = nothing
    
        # Iterate through the vector to find the maximum index according to the key
        for (index, dict) in enumerate(queue)
            if haskey(dict, key_to_check)
                current_value = dict[key_to_check]
                if min_value === nothing || current_value <= min_value
                    min_value = current_value
                    min_index = index
                end
            end
        end

        return min_value, min_index
    end


end



lb = [0.0, 0.0, 0.0]
ub = [1.0, 1.0, 1.0]
cobj = [1.0, 1.0, 2.0]
vtype = ['B', 'B', 'B']
A = [1.0 2.0 3.0; -1.0 -1.0 0.0]
b = [4.0; -1.0]

integerity_tol = 1e-3;
optimality_gap = 1e-2;





result = branch_and_bound(cobj, A, b, lb, ub, vtype, "max")
print(result)