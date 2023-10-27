# Integer programming models
Solving integer programs is a difficult task in general. This document give some basic definations on the mixed-integer programming problems.

The standard format of binary programming problem is defined as follows:
```
$$ \max_{\mathbf{x}} \mathbf{c}^{\text{T}}\mathbf{x} $$
$$ \text{s.t.} \mathbf{A}\mathbf{x} \leq \mathbf{b}\\
$$ \mathbf{x} \in \{0, 1\}^{n},$$
```

The mixed-integer linear programming model, including both nonnegative discrete vector x and continuous vector y, with shared constraints, as follows:
