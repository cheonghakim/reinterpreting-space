"""
mp_model.jl

Exact algebraic model for the Majumdar-Papapetrou lapse potential
    U(x,y;M1,M2,M3) = 1 + M1/|x-P1| + M2/|x-P2| + lambda*...
for three sources at the vertices of a circumradius-1 equilateral triangle,
and the symmetry-adapted (r,s=t^2) reduction used for the symmetric family
(M1,M2,M3) = (1,1,lambda).

Every function here is generic in the element type T (Float64, Interval{Float64},
Interval{BigFloat}, Arb, ...) so the SAME code path is used for plain float
sanity checks and for outward-rounded / ball-arithmetic certification. This
is deliberate: it removes an entire class of "the interval code computes a
different formula than the float code" bugs.

Cartesian source positions (exact):
    P1 = ( 1,        0,        0)
    P2 = (-1/2,  sqrt(3)/2,    0)
    P3 = (-1/2, -sqrt(3)/2,    0)

This file is `include`d into the `ValidatedMP` module (src/ValidatedMP.jl);
it is not itself a standalone module, so that functions here, in
derivatives.jl, and in interval_utils.jl all share one namespace without
fragile relative-import paths.
"""

# ---- Cartesian model (used for the full 2D completeness / degree sections) ----

"Exact source positions as a function of the scalar type T (so BigFloat / Interval
constants are constructed with correct precision, not truncated Float64 literals)."
function cartesian_sources(::Type{T}) where {T}
    s3 = sqrt(T(3))
    P1 = (T(1), T(0))
    P2 = (T(-1)/2, s3/2)
    P3 = (T(-1)/2, -s3/2)
    return (P1, P2, P3)
end

"U(x,y) for sources M1,M2,M3 at the standard triangle, generic element type."
function U_cartesian(x, y, M1, M2, M3)
    T = typeof(x)
    P1, P2, P3 = cartesian_sources(T)
    r1 = sqrt((x - P1[1])^2 + (y - P1[2])^2)
    r2 = sqrt((x - P2[1])^2 + (y - P2[2])^2)
    r3 = sqrt((x - P3[1])^2 + (y - P3[2])^2)
    return 1 + M1 / r1 + M2 / r2 + M3 / r3
end

"""
U_symmetric(x,y,lambda) = U_cartesian(x,y,lambda,1,1): the symmetric family
(M1,M2,M3)=(1,1,lambda), with lambda attached to P1=(1,0,0) per the
verified (r,t)=(x,y) correspondence documented above `R1sq`.
"""
U_symmetric(x, y, lambda) = U_cartesian(x, y, lambda, one(lambda), one(lambda))

"""
U_asymmetric(x,y,eps,lambda) = U_cartesian(x,y,lambda,1+eps,1-eps): the
asymmetric family (M1,M2,M3)=(1+eps,1-eps,lambda) used in §12, with the
SAME lambda-at-P1 convention, and the unequal masses 1+-eps at P2,P3 (i.e.
attached to the two sources that carried equal mass 1 in the symmetric
family — equivalent up to an M1<->M2 relabeling of the (P1,P2,P3) vertices,
which is immaterial: it is the same triangle up to a relabeling of which
vertex is "P1").
"""
function U_asymmetric(x, y, eps, lambda)
    return U_cartesian(x, y, lambda, one(lambda) + eps, one(lambda) - eps)
end

# ---- Symmetry-adapted (r, t) axis coordinates ----
#
# VERIFIED CORRESPONDENCE (checked by direct algebraic expansion, not assumed):
# (r,t) IS the same plane as Cartesian (x,y) above, i.e. x=r, y=t directly,
# under the labeling where the lambda-weighted source sits at P1=(1,0,0) and
# the two mass-1 sources sit at P2=(-1/2,sqrt(3)/2), P3=(-1/2,-sqrt(3)/2):
#   dist((r,t),P2)^2 = (r+1/2)^2+(t-sqrt(3)/2)^2 = r^2+r+1+t^2-sqrt(3)t = R1^2
#   dist((r,t),P3)^2 = (r+1/2)^2+(t+sqrt(3)/2)^2 = r^2+r+1+t^2+sqrt(3)t = R2^2
#   dist((r,t),P1)^2 = (r-1)^2+t^2                                      = R3^2
# So V(r,s,lambda) below equals U_cartesian(r,t,1,1,lambda) exactly with s=t^2,
# lambda attached to the source at (1,0,0). This is used directly in
# cross-validation between the axis-reduced and full-2D code paths.

"R1^2(r,t) exact."
R1sq(r, t) = r^2 + r + t^2 - sqrt(oftype(r, 3)) * t + 1

"R2^2(r,t) exact."
R2sq(r, t) = r^2 + r + t^2 + sqrt(oftype(r, 3)) * t + 1

"R3^2(r,t) exact."
R3sq(r, t) = (1 - r)^2 + t^2

"""
V(r,s,lambda) using the EXACT identity for the symmetric pair of sources,
s standing for t^2 (s >= 0 required; D = a^2 - 3s must be > 0, i.e. the
point must not coincide with either off-axis source).

    a = r^2 + r + 1 + s
    D = a^2 - 3s
    1/sqrt(R1^2) + 1/sqrt(R2^2) = sqrt( 2a/D + 2/sqrt(D) )
    V = 1 + sqrt(2a/D + 2/sqrt(D)) + lambda/sqrt((1-r)^2+s)

This is proven algebraically identical to U(r,t,lambda) with t^2=s (see
README.md, "Symmetric-axis reduction"), not a numerical approximation.
"""
function V_rs(r, s, lambda)
    a = r^2 + r + 1 + s
    D = a^2 - 3s
    sumterm = sqrt(2a / D + 2 / sqrt(D))
    tail = lambda / sqrt((1 - r)^2 + s)
    return 1 + sumterm + tail
end

"""
lambda_onaxis(r)

On-axis (s=0, i.e. t=0) stationarity condition. At s=0, V_rs(r,0,lambda)
reduces exactly to 1 + 2/sqrt(r^2+r+1) + lambda/(1-r) (checked: a=r^2+r+1,
D=a^2, so sqrt(2a/D+2/sqrt(D)) = sqrt(4/a) = 2/sqrt(a)). Setting
d/dr V_rs(r,0,lambda) = 0 for r<1 gives

    lambda(r) = (2r+1)*(1-r)^2 / (r^2+r+1)^(3/2)

which is what this function computes, written as a*sqrt(a) rather than
a^1.5 so it works uniformly for Interval, Arb, and plain float element
types without relying on fractional-power support.
"""
function lambda_onaxis(r)
    a = r^2 + r + 1
    return (2r + 1) * (1 - r)^2 / (a * sqrt(a))
end

"""
onaxis_r_derivs(r, lambda)

Analytic closed-form r-derivatives of g(r,lambda) := V_rs(r,0,lambda)
= 1 + 2/sqrt(a) + lambda/(1-r), a = r^2+r+1, derived by hand (not AD) as an
INDEPENDENT cross-check against the ForwardDiff nested-derivative path in
derivatives.jl:

    U_r    =  -(2r+1)/a^(3/2)          + lambda/(1-r)^2
    U_rr   =  -2/a^(3/2) + (3/2)(2r+1)^2/a^(5/2)      + 2*lambda/(1-r)^3
    U_rrr  =  9(2r+1)/a^(5/2) - (15/4)(2r+1)^3/a^(7/2) + 6*lambda/(1-r)^4
    U_rlambda = 1/(1-r)^2   (= d/dlambda of the lambda/(1-r)^2 term in U_r;
                              independent of lambda itself)

All half-integer powers of a are written as a^k * sqrt(a) (integer k) so the
formula uses only +,-,*,/,^Int,sqrt and works for Float64, Interval, and Arb
element types alike. Only valid for r<1 (needed since |1-r| was replaced by
1-r); true throughout this project's r_+, r_- boxes.
"""
function onaxis_r_derivs(r, lambda)
    a = r^2 + r + 1
    sqa = sqrt(a)
    b = 2r + 1
    q = 1 - r
    Ur = -b / (a * sqa) + lambda / q^2
    Urr = -2 / (a * sqa) + (3 // 2) * b^2 / (a^2 * sqa) + 2 * lambda / q^3
    Urrr = 9 * b / (a^2 * sqa) - (15 // 4) * b^3 / (a^3 * sqa) + 6 * lambda / q^4
    Urlambda = 1 / q^2
    return (Ur = Ur, Urr = Urr, Urrr = Urrr, Urlambda = Urlambda)
end
