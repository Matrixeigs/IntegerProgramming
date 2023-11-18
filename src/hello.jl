using Gurobi

lb = [0.0, 0.0]
ub = [Inf, Inf]
cobj = Cdouble[-1.0, -1.0]
vtype = Cchar[GRB_INTEGER, GRB_INTEGER]
A = [2.0 -2.0; -8.0 10.0]
b = [-1.0; 13.0]


# initialize model
env_p = Ref{Ptr{Cvoid}}()
error = GRBloadenv(env_p, "lp.log")
env = env_p[]
model_p = Ref{Ptr{Cvoid}}()
error = GRBnewmodel(env, model_p, "lp", 0, C_NULL, C_NULL, C_NULL, C_NULL, C_NULL)
model = model_p[]

error = GRBaddvars(
  model, # model
  2,      # : numvars
  0,      # : numnz
  C_NULL, # : *vbeg
  C_NULL, # : *vind
  C_NULL, # : *vval
  cobj,   # : *obj
  lb,    # : *lb
  ub,    # : *ub
  vtype, # : *vtype
  # "xyz" # : **varnames
  ["xy"]   # : **varnames
  # Cuchar["xyz"]   # : **varnames
  # Base.unsafe_convert(Cstring,"xyz")
  # Base.cconvert(Cstring,"xyz")
)

for i in 1 : 2
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
# (4) Solve problem
error = GRBoptimize(model)
error = GRBwrite(model, "lp.lp");

pinfeas = Ref{Cdouble}()
dinfeas = Ref{Cdouble}()
relgap = Ref{Cdouble}()
NumVars = Ref{Cint}()
NumConstrs = Ref{Cint}()
IterCount = Ref{Cint}() # simplex iters
BarIterCount = Ref{Cint}() # barrier iters
optimstatus = Ref{Cint}() # barrier iters
objval = Ref{Cdouble}() # barrier iters
sol = ones(2)

GRBgetdblattr(model, "ConstrVio", pinfeas) # maximum (primal) constraint violation
GRBgetdblattr(model, "MaxVio", dinfeas) # sum of (dual) constraint violations
GRBgetdblattr(model, "ComplVio", relgap) # complementarity violation
GRBgetintattr(model, "NumVars", NumVars) # sum of (dual) constraint violations
GRBgetintattr(model, "NumConstrs", NumConstrs) # sum of (dual) constraint violations
GRBgetintattr(model, "BarIterCount", BarIterCount) # sum of (dual) constraint violations

error = GRBgetintattr(model, GRB_INT_ATTR_STATUS, optimstatus);
error = GRBgetdblattr(model, GRB_DBL_ATTR_OBJVAL, objval);
error = GRBgetdblattrarray(model, GRB_DBL_ATTR_X, 0, 2, sol);

print(sol)