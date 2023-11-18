#=
minimize    x^2 + x*y + y^2 + y*z + z^2 + 2 x
subject to  x + 2 y + 3 z >= 4
            x +   y       >= 1
            x, y, z non-negative
=#

function linear_programming(c::Vector, A::Matrix, b::Vector, lb::Vector, ub::Vector)
    # A linear programming wrapper to call Gurobi
    env_p = Ref{Ptr{Cvoid}}()
    error = GRBloadenv(env_p, "lp.log")
    env = env_p[]
    model_p = Ref{Ptr{Cvoid}}()
    error = GRBnewmodel(env, model_p, "lp", 0, C_NULL, C_NULL, C_NULL, C_NULL, C_NULL)
    model = model_p[]
    
  end
  
  
  
  using Gurobi
  ENV["JULIA_NUM_THREADS"]
  # initialize model
  env_p = Ref{Ptr{Cvoid}}()
  error = GRBloadenv(env_p, "qp.log")
  env = env_p[]
  model_p = Ref{Ptr{Cvoid}}()
  error = GRBnewmodel(env, model_p, "qp", 0, C_NULL, C_NULL, C_NULL, C_NULL, C_NULL)
  model = model_p[]
  
  # variables bounds
  vlb = Cdouble[-Inf, -Inf, -Inf]
  vub = Cdouble[+Inf, +Inf, +Inf]
  
  # linear objective coeffs
  cobj = Cdouble[2.0, 0.0, 0.0]
  
  # quadratic objective coeffs
  qrow = Cint[0, 0, 1, 1, 2]
  qcol = Cint[0, 1, 1, 2, 2]
  qval = Cdouble[1, 1, 1, 1, 1]
  
  # linear constraint1 coeffs
  c1ind = Cint[0, 1, 2]
  c1val = Cdouble[1.0, 2.0, 3.0]
  c1sense = GRB_GREATER_EQUAL
  c1rhs = 4.0
  
  # linear constraint2 coeffs
  c2ind = Cint[0, 1]
  c2val = Cdouble[1.0, 1.0]
  c2sense = GRB_GREATER_EQUAL
  c2rhs = 1.0
  
  # linear constraint3 coeffs: nonnegativity
  c3xind = Cint[0]
  c3xval = Cdouble[1.0]
  c3xsense = GRB_GREATER_EQUAL
  c3xrhs = 0.0
  c3yind = Cint[1]
  c3yval = Cdouble[1.0]
  c3ysense = GRB_GREATER_EQUAL
  c3yrhs = 0.0
  c3zind = Cint[2]
  c3zval = Cdouble[1.0]
  c3zsense = GRB_GREATER_EQUAL
  c3zrhs = 0.0
  c3ind = [c3xind, c3yind, c3zind]
  c3val = [c3xval, c3yval, c3zval]
  c3sense = [c3xsense, c3ysense, c3zsense]
  c3rhs = [c3xrhs, c3yrhs, c3zrhs]
  
  # variables and objective
  vname = []
  vname = append!(vname,["xyz"])
  error = GRBaddvars(
    model, # model
    3,      # : numvars
    0,      # : numnz
    C_NULL, # : *vbeg
    C_NULL, # : *vind
    C_NULL, # : *vval
    cobj,   # : *obj
    vlb,    # : *lb
    vub,    # : *ub
    C_NULL, # : *vtype
    # "xyz" # : **varnames
    ["xyz"]   # : **varnames
    # Cuchar["xyz"]   # : **varnames
    # Base.unsafe_convert(Cstring,"xyz")
    # Base.cconvert(Cstring,"xyz")
  )
  error = GRBaddqpterms(model, 5, qrow, qcol, qval)
  
  # constraint1
  error = GRBaddconstr(
    model,   # : *model
    3,       # : numnz
    c1ind,   # : *cind
    c1val,   # : *cval
    c1sense, # : sense
    c1rhs,   # : rhs
    "c1",    # : *constrname
  )
  
  # constraint2
  error = GRBaddconstr(
    model,   # : *model
    2,       # : numnz
    c2ind,   # : *cind
    c2val,   # : *cval
    c2sense, # : sense
    c2rhs,   # : rhs
    "c2",    # : *constrname
  )
  
  # constraint3: nonnegativity
  for i in 1:3
    error = GRBaddconstr(
      model,   # : *model
      1,       # : numnz
      c3ind[i],   # : *cind
      c3val[i],   # : *cval
      c3sense[i], # : sense
      c3rhs[i],   # : rhs
      "c3_$i",    # : *constrname
    )
  end
  error = GRBoptimize(model)
  error = GRBwrite(model, "qp.lp");
  
  pinfeas = Ref{Cdouble}()
  dinfeas = Ref{Cdouble}()
  relgap = Ref{Cdouble}()
  NumVars = Ref{Cint}()
  NumConstrs = Ref{Cint}()
  IterCount = Ref{Cint}() # simplex iters
  BarIterCount = Ref{Cint}() # barrier iters
  optimstatus = Ref{Cint}() # barrier iters
  objval = Ref{Cdouble}() # barrier iters
  sol = ones(3)
  
  GRBgetdblattr(model, "ConstrVio", pinfeas) # maximum (primal) constraint violation
  GRBgetdblattr(model, "MaxVio", dinfeas) # sum of (dual) constraint violations
  GRBgetdblattr(model, "ComplVio", relgap) # complementarity violation
  GRBgetintattr(model, "NumVars", NumVars) # sum of (dual) constraint violations
  GRBgetintattr(model, "NumConstrs", NumConstrs) # sum of (dual) constraint violations
  GRBgetintattr(model, "BarIterCount", BarIterCount) # sum of (dual) constraint violations
  
  error = GRBgetintattr(model, GRB_INT_ATTR_STATUS, optimstatus);
  error = GRBgetdblattr(model, GRB_DBL_ATTR_OBJVAL, objval);
  error = GRBgetdblattrarray(model, GRB_DBL_ATTR_X, 0, 3, sol);