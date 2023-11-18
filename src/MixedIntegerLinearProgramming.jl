using Gurobi

function mixed_integer_linear_programming(cobj::Vector, A::Matrix, b::Vector, lb::Vector, ub::Vector, vtype::Vector, model_sense::String)
    # 0: initialize model
    env_p = Ref{Ptr{Cvoid}}()
    error = GRBloadenv(env_p, "lp.log")
    env = env_p[]

    GRBsetparam(env, "OutputFlag", "0") # Update environment parameters
    model_p = Ref{Ptr{Cvoid}}()
    error = GRBnewmodel(env, model_p, "milp", 0, C_NULL, C_NULL, C_NULL, C_NULL, C_NULL)
    model = model_p[]

    # 1:  add variables
    NumVars = Ref{Cint}()
    NumConstrs = Ref{Cint}()
    (NumConstrs, NumVars) = size(A)
    vtype = map(x -> Int8(x), vtype) # Change type value
    error = GRBaddvars(
        model, # model
        NumVars,      # : numvars
        0,      # : numnz
        C_NULL, # : *vbeg
        C_NULL, # : *vind
        C_NULL, # : *vval
        cobj,   # : *obj
        lb,    # : *lb
        ub,    # : *ub
        vtype, # : *vtype
        C_NULL   # : **varnames
        )
    # 2:  add constraints
    for i in 1 : NumConstrs
        nonzero_indices = findall(!iszero, A[i,:])
        numnz = length(nonzero_indices)
        val = zeros(numnz)
        for j in 1:numnz
            val[j] = A[i, nonzero_indices[j]]
        end
        error = GRBaddconstr(
            model,   # : *model
            numnz,       # : numnz
            nonzero_indices,   # : *cind
            val,   # : *cval
            GRB_LESS_EQUAL, # : sense
            b[i],   # : rhs
            C_NULL,    # : *constrname
            )
    end
    # 3: update model parameters
    if cmp(model_sense, "min") == 0
        error = GRBsetintattr(model, "ModelSense", GRB_MINIMIZE)
    else
        error = GRBsetintattr(model, "ModelSense", GRB_MAXIMIZE)
    end
    # error = GRBsetintparam(model, "OutputFlag", 0)
    # error = GRBsetintparam(model, "OutputFlag", 0)
    error = GRBoptimize(model)
    # error = GRBwrite(model, "lp.lp");

    pinfeas = Ref{Cdouble}()
    dinfeas = Ref{Cdouble}()
    relgap = Ref{Cdouble}()

    # 4: obtain results
    IterCount = Ref{Cint}() # simplex iters
    BarIterCount = Ref{Cint}() # barrier iters
    optimstatus = Ref{Cint}() # barrier iters
    objval = Ref{Cdouble}() # barrier iters
    sol = ones(NumVars)

    GRBgetdblattr(model, "ConstrVio", pinfeas) # maximum (primal) constraint violation
    GRBgetdblattr(model, "MaxVio", dinfeas) # sum of (dual) constraint violations
    GRBgetdblattr(model, "ComplVio", relgap) # complementarity violation
    GRBgetintattr(model, "NumVars", NumVars) # sum of (dual) constraint violations
    GRBgetintattr(model, "NumConstrs", NumConstrs) # sum of (dual) constraint violations
    GRBgetintattr(model, "BarIterCount", BarIterCount) # sum of (dual) constraint violations

    error = GRBgetintattr(model, GRB_INT_ATTR_STATUS, optimstatus);
    error = GRBgetdblattr(model, GRB_DBL_ATTR_OBJVAL, objval);
    error = GRBgetdblattrarray(model, GRB_DBL_ATTR_X, 0, NumVars, sol);

    GRBfreemodel(model)
    GRBfreeenv(env)
    sol = Dict("x" => sol, "objval" => objval[])
    return sol
end

# lb = [0.0, 0.0]
# ub = [Inf, Inf]
# cobj = [1.0, 1.0]
# vtype = [GRB_INTEGER, GRB_INTEGER]
# A = [2.0 -2.0; -8.0 10.0]
# b = [-1.0; 13.0]

lb = [0.0, 0.0, 0.0]
ub = [1.0, 1.0, 1.0]
cobj = [1.0, 1.0, 2.0]
vtype = [GRB_BINARY, GRB_BINARY, GRB_BINARY]
A = [1.0 2.0 3.0; -1.0 -1.0 0.0]
b = [4.0; -1.0]

result = mixed_integer_linear_programming(cobj, A, b, lb, ub, vtype, "max")
print(result)