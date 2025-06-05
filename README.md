# P3P.jl

A minimal package for solving P3P problems in Julia. To use in Julia, run:

```Julia
pkg> add https://github.com/lopenguin/P3P.jl
```

## How to use
The main interface is the function `p3p(y, b, camK)`. The arguments of this function are:
- `y`: pixel detection coordinates `[u; v; 1]` arranged as columns of a 3x3 matrix
- `b`: corresponding 3D coordinates `[x; y; z]` arranged as columns of a 3x3 matrix
- `camK`: camera calibration matrix [3 x 3]

This function returns a list of rotations and translations that solve the P3P equations. It may be empty.

## Solvers
Currently only the [Nanako solver](https://github.com/g9nkn/p3p_problem/tree/main) is implemented. It is a direct method and ported from the [MATLAB code](https://github.com/g9nkn/p3p_problem/blob/main/p3p_nakano_bmvc2019.m) released by the author.

Contributions welcome! The Nanako solver is enough for my needs, and I hope releasing my implementation is useful for you.

# References
1. Gaku Nakano, "A Simple Direct Solution to the Perspective-Three-Point Problem," BMVC2019.  
<https://bmvc2019.org/wp-content/uploads/papers/0533-paper.pdf>
