"""
derivatives.jl

Rigorous derivatives of V(r,s,lambda) and U(x,y) via forward-mode automatic
differentiation (ForwardDiff.jl) run DIRECTLY on interval / Arb element
types. This is sound (not merely "probably fine"): IntervalArithmetic.jl's
Interval{T} implements +,-,*,/,sqrt,^ with correctly outward-rounded
interval results, so ForwardDiff's dual-number chain rule, evaluated over
Interval{T}, produces a rigorous outward-rounded ENCLOSURE of the true
derivative over the whole box — not just the derivative at one point. The
same holds for Arb ball arithmetic in Arblib.jl.

We deliberately do NOT hand-code separate symbolic derivative formulas here:
using the same V_rs / U_cartesian definitions from mp_model.jl for both the
value and (via AD) the derivatives removes an entire class of transcription
bugs where the analytic derivative formula silently drifts from the
function actually being certified.

`include`d into ValidatedMP after mp_model.jl (see that file's docstring).
"""

"""
V_second_order(r,s,lambda) -> NamedTuple(Vr, Vs, Vrr, Vrs, Vss)

First and second partial derivatives of V(r,s,lambda) w.r.t. (r,s), generic
over the element type T of r,s,lambda (Float64, Interval{Float64},
Interval{BigFloat}, ...).
"""
function V_second_order(r::T, s::T, lambda::T) where {T}
    f(v) = V_rs(v[1], v[2], lambda)
    v0 = [r, s]
    g = ForwardDiff.gradient(f, v0)
    H = ForwardDiff.hessian(f, v0)
    return (Vr = g[1], Vs = g[2], Vrr = H[1, 1], Vrs = H[1, 2], Vss = H[2, 2])
end

"""
U_grad_hess(x,y,M1,M2,M3) -> NamedTuple(Ux, Uy, Uxx, Uxy, Uyy)

Full 2D Cartesian gradient and Hessian of U at (x,y), generic over element
type. Used for the completeness / Brouwer-degree / phase-slab sections
(§9-11), which are not restricted to the symmetric axis.
"""
function U_grad_hess(x::T, y::T, M1::T, M2::T, M3::T) where {T}
    f(v) = U_cartesian(v[1], v[2], M1, M2, M3)
    v0 = [x, y]
    g = ForwardDiff.gradient(f, v0)
    H = ForwardDiff.hessian(f, v0)
    return (Ux = g[1], Uy = g[2], Uxx = H[1, 1], Uxy = H[1, 2], Uyy = H[2, 2])
end

"""
V_full_second_order(r,s,lambda) -> NamedTuple with the full gradient/Hessian
of V_rs w.r.t. ALL THREE variables (r,s,lambda), from a SINGLE (non-nested)
ForwardDiff.hessian call on the 3-vector [r,s,lambda]. Used in §7 to get the
mixed partial V_slambda needed for d/dlambda[Q(0,lambda)] without stacking
multiple levels of Dual numbers (nesting works too, as onaxis_r_derivs_ad
shows, but a single flat 3-variable AD call is simpler and more robust, so
it is preferred here).
"""
function V_full_second_order(r::T, s::T, lambda::T) where {T}
    f(v) = V_rs(v[1], v[2], v[3])
    v0 = [r, s, lambda]
    g = ForwardDiff.gradient(f, v0)
    H = ForwardDiff.hessian(f, v0)
    return (Vr = g[1], Vs = g[2], Vlambda = g[3],
            Vrr = H[1, 1], Vrs = H[1, 2], Vrlambda = H[1, 3],
            Vss = H[2, 2], Vslambda = H[2, 3], Vlambdalambda = H[3, 3])
end

"""
onaxis_r_derivs_ad(r, lambda) -> NamedTuple(Ur, Urr, Urrr, Urlambda)

AD cross-check of onaxis_r_derivs (mp_model.jl's hand-derived closed form),
computed instead by NESTED ForwardDiff.derivative on
g(r,lambda) = V_rs(r, zero(r), lambda). Independent derivation path (AD vs.
analytic algebra) for the same quantities used in the fold nondegeneracy
certificate (§8) — used to cross-validate onaxis_r_derivs, not to replace
it.
"""
function onaxis_r_derivs_ad(r::T, lambda::T) where {T}
    g(rr, ll) = V_rs(rr, zero(rr), ll)

    d1(rr, ll) = ForwardDiff.derivative(x -> g(x, ll), rr)
    d2(rr, ll) = ForwardDiff.derivative(x -> d1(x, ll), rr)
    d3(rr, ll) = ForwardDiff.derivative(x -> d2(x, ll), rr)

    Ur = d1(r, lambda)
    Urr = d2(r, lambda)
    Urrr = d3(r, lambda)
    Urlambda = ForwardDiff.derivative(ll -> d1(r, ll), lambda)
    return (Ur = Ur, Urr = Urr, Urrr = Urrr, Urlambda = Urlambda)
end
