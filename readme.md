# Levenberg-Marquardt optimziation

Note: This implementation uses an inexact line-search.

### Step 1:

Assume the nonlinear energy to be

```
EQ(x) = ||A(x) + b||^2 + g(A(x))
```

where ```A(x)``` is a nonlinear vector of ```x``` and ```g``` is a scalar function of ```A(x)```.

### Step 2:
We can then linearize the vector ```A(x)``` around an operating point ```x0``` as

```
A(d) ~ A(x0) + sum_i dAdi(x0) d_i
```

where

```
d = (d_1 d_2 ...)
```

is a step from the current operating point ```x0``` as in ```d := x - x0```.

### Step 3:
Replacing the nonlinear ```A(x)``` by the linearized vector ```a(d)``` into the energy yields the linearized energy ```EQLIN(x0,d)``` which we rewrite as

```
EQLIN(xb,d) = 0.5 d^T H(xb)^T d  + f(xb)^T d + 0.5 c
```

Minimizing ```EQLIN``` essentially comes down to solving a linear system. This is then done many times around subsequent, hopefully improving operating points.