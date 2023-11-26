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
    problem = Dict("cojb" => cobj, "A"=>A,"b"=>b, "lb"=>lb, "ub"=>ub, "model_sense"=>model_sense)
    root_node = Dict("bound" => result0["objval"], 
    "x" => result0["x"], 
    "depth"=>0,
    "feasibility" => sum(integerity_feasibility(result0["x"], vtype, integerity_tol)),
    "problem"=>problem)
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
    iter = 0
    while length(queue)>0 # The queue is not empty
        #  Step 3.0) Pop the top node from the queue
        if f_opt > -Inf # A feasible solution is derived.
            (~,index) = queue_management(queue, "bound", "max")
        else # An infeasbile solution is derived.
            (~,index) = queue_management(queue, "feasibility", "min")
        end
        node = queue[index]
        splice!(queue, index)# Pop out this node 
        if node["depth"] > 0 && node["bound"] < f_opt # 3.1) Check if the node can be pruned
            continue
        else
            # 3.2) Update incumbent solution
            x = node["x"]
            problem = node["problem"]
            if node["feasibility"] <= integerity_tol
                obj = node["bound"]
                if obj>f_opt
                    f_opt = obj
                    x_incumbent = x
                end
            else
                # 3.3) branch on given variables
                [~, branch_idx] = max(abs(x_frac.*problem["cobj"]))
                A_x = zeros(1, nx)
                A_x(1,branch_idx) = 1
                # 3.3.1) Create two subproblems
                # (1) x(branch_idx)<= fix(x(branch_idx))
                problem_left = problem
                problem_left["A"] = [problem_left["A"]; A_x]
                problem_left["b"] = [problem_left["b"]; floor(x[branch_idx])]
                # (2) x(branch_idx)>=ceil(x(branch_idx))
                problem_right = problem
                problem_right["A"] = [problem_left["A"]; -A_x]
                problem_right["b"] = [problem_right["b"]; -ceil(x[branch_idx])]
                # 3.3.2) Solve these two problems
                solution_left = gurobi(problem_left, options);
                solution_right = gurobi(problem_right, options);
                if strcmp(solution_left.status, "OPTIMAL") # 3.4) Prune by infeasibility
                    node_left.bound = solution_left.objval;
                    node_left.x = solution_left.x;
                    node_left.problem = problem_left;
                    node_left.depth = node.depth + 1;
                    node_left.feasibility = sum(integrity_check(node_left.x, vtype) >= integerity_feasibility);
                    if node_left.bound > f_opt # Prune by optimality
                        queue = [queue, node_left];
                    end
                end
                if strcmp(solution_right.status, "OPTIMAL") # 3.4) Prune by infeasibility
                    node_right.bound = solution_right.objval;
                    node_right.x = solution_right.x;
                    node_right.problem = problem_right;
                    node_right.depth = node.depth + 1;
                    node_right.feasibility = sum(integrity_check(node_right.x, vtype) >= integerity_feasibility);
                    if node_right.bound > f_opt # Prune by optimality
                        queue = [queue, node_right];
                    end
                end
            end
        end
        
        # 3.5) Update the gap
        if strcmp(problem.modelsense, "max")
            f_best = max([queue.bound]);
        else
            f_best = min([queue.bound]);
        end
        gap = (f_best-f_opt)/max([abs(f_best), abs(f_opt), 0]);
        disp(gap);
        if gap<optimality_gap
            queue([queue.feasibility]>0) = [ ]; # Delete infeasible solutions
            queue([queue.bound] < f_opt) = [ ]; # Delete suboptimal solutions
            break
        end
        disp(node.feasibility);
        disp(node.bound);
        iter = iter + 1;
    
    end
    
    return x_incumbent, f_opt

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