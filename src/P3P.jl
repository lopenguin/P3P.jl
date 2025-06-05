module P3P

using LinearAlgebra
import Polynomials, FastPolynomialRoots

"""
    p3p(y, b, camK)

Unproject pixel coordinates and run p3p.

# Inputs:
- `y`: homogenized pixel coordinates for each column `[u; v; 1]` [3 x 3]
- `b`: corresponding 3D points in each column [3 x 3]
- `camK`: camera calibration matrix [3 x 3]
"""
function p3p(y, b, camK)

    # unproject pixel coordinates
    y_normed = camK \ y
    # normalize columns for Nakano
    y_normed = reduce(hcat,normalize.(eachcol(y_normed)))

    # run p3p
    Rs, ts = p3p_nakano(y_normed, b)

    return Rs, ts
end

"""
Runs Nakano p3p with normalized coordinates.

# Inputs
- `y`: unprojected pixels in each column, unit column norm [3 x 3]
- `b`: cooresponding 3D points in each column [3 x 3]
- `polishing`: number of polishing iterates to run
"""
function p3p_nakano(y, b; polishing=1)
    # permute y and b so the longest distance is between 1 and 2
    d = norm.(eachcol([b[:,1]-b[:,2]  b[:,2]-b[:,3]  b[:,1]-b[:,3]]))
    if argmax(d) == 2
        y = y[:,[2,3,1]]
        b = b[:,[2,3,1]]
    elseif argmax(d) == 3
        y = y[:,[1,3,2]]
        b = b[:,[1,3,2]]
    end

    # rigid transform so all points on z=0 plane
    b21 = b[:,2] - b[:,1]
    b31 = b[:,3] - b[:,1]
    nx = normalize(b21)
    nz = normalize(cross(nx, b31))
    n = [nx  cross(nz, nx)  nz]

    # polynomial cooefficients
    pa = n[:,1]'*b21
    pb = n[:,1]'*b31
    pc = n[:,2]'*b31

    y12 = y[:,1]'*y[:,2]
    y13 = y[:,1]'*y[:,3]
    y23 = y[:,2]'*y[:,3]
    p = pb/pa
    q = (pb^2+pc^2)/(pa^2)

    f = [p; -y23; 0; -y12*(2*p-1); y13; p-1]
    g = [q; 0; -1; -2*y12*q; 2*y13; q-1]

    h = [-f[1]^2 + g[1]*f[2]^2;
         f[2]^2*g[4] - 2*f[1]*f[4] - 2*f[1]*f[2]*f[5] + 2*f[2]*f[5]*g[1];
         f[5]^2*g[1] - 2*f[1]*f[5]^2 - 2*f[1]*f[6] + f[2]^2*g[6] - f[4]^2 - 2*f[2]*f[4]*f[5] + 2*f[2]*f[5]*g[4];
         f[5]^2*g[4] - 2*f[4]*f[5]^2 - 2*f[4]*f[6] - 2*f[2]*f[5]*f[6] + 2*f[2]*f[5]*g[6];
         -2*f[5]^2*f[6] + g[6]*f[5]^2 - f[6]^2
        ]

    # root finding
    # h(1)*x^4 + h(2)*x^3 + h(3)*x^2 + h(4)*x + h(5) = 0
    x = Polynomials.roots(Polynomials.Polynomial(reverse(h)))
    x = real(x[real(x) .> 0 .&& abs.(imag(x)) .< 1e-8])
    z = -((f[1]*x .+ f[4]).*x .+ f[6]) ./ (f[5] .+ f[2]*x)

    if polishing > 0
        x, z = rootpolishing(f, g, x, z, polishing)
    end

    # recover motion
    nsols = length(x)
    A = y .* [-1 1 0]
    B = y .* [-1 0 1]
    C = B - p*A

    R = zeros(3,3,nsols)
    t = zeros(3,nsols)
    for i = 1:nsols
        lambda = [1; x[i]; z[i]]
        s = norm(A*lambda) / pa
        d = lambda / s

        r1 = (A*d) / pa         
        r2 = (C*d) / pc

        Rc = [r1 r2 cross(r1,r2)]
        tc = d[1]*y[:,1]

        # if isdefined(Main, :Infiltrator)
        #     Main.infiltrate(@__MODULE__, Base.@locals, @__FILE__, @__LINE__)
        # end

        R[:,:,i] = Rc*n'
        t[:,i]   = tc - Rc*n'*b[:,1]
    end

    return R, t
end

"""
    rootpolishing(f, g, x, y, maxitr)

Improve root quality via Gauss-Newton.

From Nakano p3p code.
"""
function rootpolishing(f, g, x, y, maxitr=1)

    for _ = 1:maxitr
        x2 = x.^2
        xy = x.*y
        y2 = y.^2

        fv = f[1]*x2 + f[2]*xy + f[4]*x + f[5]*y .+ f[6]
        gv = g[1]*x2 - y2 + g[4]*x + g[5]*y .+ g[6]

        ii = abs.(fv).<1e-15 .&& abs.(gv).<1e-15

        dfdx = 2*f[1]*x + f[2]*y .+ f[4]
        dfdy =   f[2]*x          .+ f[5]
        dgdx = 2*g[1]*x          .+ g[4]
        dgdy =              -2*y .+ g[5]

        inv_detJ  = 1 ./ (dfdx.*dgdy - dfdy.*dgdx);

        dx = ( dgdy.*fv - dfdy.*gv) .* inv_detJ;
        dy = (-dgdx.*fv + dfdx.*gv) .* inv_detJ;

        dx[ii] .= 0
        dy[ii] .= 0
        
        x = x - dx
        y = y - dy
    end

    return x, y
end

export p3p

end # module P3P
