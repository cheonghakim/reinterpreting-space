"""
krawczyk.jl

Generic n-dimensional interval Newton/Krawczyk root isolation, used by §7
(1D specialization via monotonicity), §9 (2D, Ux=Uy=0), §10/§11 (3D,
F=(Ux,Uy,detH)=0).

`include`d into ValidatedMP (see mp_model.jl docstring for the pattern).

THEORY (why this certifies existence+uniqueness, not just "probably a root"):
For F: R^n -> R^n continuously differentiable, box X with midpoint m, and Y
an approximate inverse of the Jacobian at m (any nonsingular real matrix —
does not need to be a rigorous inverse, only used as a preconditioner), the
Krawczyk operator is

    K(X) = m - Y*F(m) + (I - Y*J(X)) * (X - m)

where J(X) is an interval enclosure of the Jacobian over all of X. If
K(X) is a STRICT SUBSET of the interior of X (`isinterior`), then the
Brouwer fixed point theorem applied to the (rigorously enclosing) affine
map guarantees F has AT LEAST ONE root in X; if in addition
I - Y*J(X) has operator norm < 1 in a suitable sense (guaranteed here
because K(X) ⊂⊂ X together with the mean value form implies the map is a
contraction on X), the root is UNIQUE in X. This is the standard
Krawczyk-operator argument (Krawczyk 1969; Moore, Kearfott & Cloud,
"Introduction to Interval Analysis", 2009). If K(X) ∩ X = ∅, F has NO root
in X (since any root x* satisfies x* ∈ K(X) by the mean value form). Any
other outcome is UNRESOLVED and the box must be subdivided.
"""

using LinearAlgebra

"One Krawczyk step. F: Vector{T} -> Vector{T}, X::Vector{T}.
Returns (K, singular) where singular=true means the point-Jacobian
preconditioner was not invertible (caller must treat this box as
UNRESOLVED without trusting K)."
function krawczyk_operator(F, X::Vector{T}) where {T}
    n = length(X)
    mpt = [_midpoint_thin(x) for x in X]      # thin interval at the box midpoint
    Fm = F(mpt)
    J = ForwardDiff.jacobian(F, X)
    Jmid_float = [Float64(_to_float(J[i, j])) for i in 1:n, j in 1:n]
    if abs(det(Jmid_float)) < 1e-300 || !isfinite(det(Jmid_float))
        return (X, true)
    end
    Y = inv(Jmid_float)
    Yt = [_thinT(Y[i, j], X[1]) for i in 1:n, j in 1:n]
    Xc = X .- mpt
    Ii = [_thinT(i == j ? 1.0 : 0.0, X[1]) for i in 1:n, j in 1:n]
    # IMPORTANT: compute C = I - Y*J as ONE interval matrix BEFORE multiplying
    # by Xc. Expanding to `Xc - Y*J*Xc` instead would make Xc appear TWICE in
    # the expression (once bare, once inside J*Xc); naive interval arithmetic
    # then treats the two occurrences as independent, destroying the
    # correlated cancellation that gives Krawczyk its O(width^2) contraction
    # (the classic interval "dependency problem" — verified by hand: doing it
    # the wrong way gives K with width a CONSTANT multiple of the box width,
    # never shrinking faster, instead of the expected quadratic shrink).
    C = Ii .- Yt * J
    K = mpt .- Yt * Fm .+ C * Xc
    return (K, false)
end

_midpoint_thin(x::Interval{T}) where {T} = interval(T(mid(x)))
_to_float(x::Interval) = mid(x)
_to_float(x::Real) = x
_thinT(v::Real, ref::Interval{T}) where {T} = interval(T(v))

"""
KrawczykResult: one leaf box from `isolate_roots`, with its final status.
status ∈ (:unique, :none, :excluded, :unresolved)
"""
struct KrawczykResult{T}
    box::Vector{T}
    status::Symbol
    depth::Int
end

_all_finite(v) = all(isfinite(inf(x)) && isfinite(sup(x)) for x in v)
_all_finite(M::AbstractMatrix) = all(isfinite(inf(x)) && isfinite(sup(x)) for x in M)

"""
isolate_roots(F, box0::Vector{Interval}; maxdepth=60, mindiam=1e-14, exclude=nothing) -> Vector{KrawczykResult}

Subdivision-driven Krawczyk root isolation, with two robustness layers added
on top of the bare Krawczyk step (both needed once this is used near source
singularities in §9-§11, not just on smooth toy systems):

1. `exclude(box)::Bool`, if given, is checked FIRST on every box. If it
   returns true (the box is provably root-free by an EXTERNAL argument —
   e.g. the source-dominance disk argument in README.md, "Exclusion disks")
   the box is recorded with status `:excluded` and never evaluated at all —
   this is how boxes touching or containing a 1/r singularity are kept away
   from F entirely, rather than relying on interval arithmetic to somehow
   produce a sensible enclosure of an infinite value.
2. Before running the (more expensive, Jacobian-based) Krawczyk step, F(box)
   itself is evaluated directly. If any component's interval enclosure
   excludes 0, the box is immediately :none — a cheap, purely-value-based
   exclusion that also protects the Krawczyk step from ever being handed a
   box where F is nearly singular. If F(box) is not finite (can happen very
   close to a singularity not yet fully covered by `exclude`), the box is
   forcibly bisected (never passed to the Jacobian/Krawczyk machinery, which
   is not guaranteed sound on non-finite input) until it either shrinks into
   the `exclude` region or resolves normally.

Status labels produced: :unique, :none, :excluded, :unresolved. Every input
box's fate is accounted for in exactly one of these — none are silently
dropped.
"""
function isolate_roots(F, box0::Vector{T}; maxdepth::Int = 60, mindiam::Float64 = 1e-14,
                        exclude = nothing) where {T <: Interval}
    results = KrawczykResult{T}[]
    stack = [(box0, 0)]
    while !isempty(stack)
        box, depth = pop!(stack)
        if any(isempty_interval, box)
            continue
        end
        maxwidth = maximum(diam(b) for b in box)
        function subdivide_or_unresolved(tag::Symbol)
            if depth >= maxdepth || maxwidth < mindiam
                push!(results, KrawczykResult(box, tag, depth))
            else
                widest = argmax(i -> diam(box[i]), 1:length(box))
                # α=0.4 (not the default 0.5) is deliberate: a root sitting
                # EXACTLY at a symmetric domain's midpoint (e.g. y=0 on the
                # reflection axis, which is common here) would otherwise stay
                # forever ON the boundary between the two halves under
                # repeated exact bisection, so it can never land in either
                # half's INTERIOR and `isinterior` can never certify it —
                # verified by hand to be the actual failure mode for the
                # lambda=0.5 on-axis roots (they hit maxdepth stuck at
                # ~1e-11 width, straddling y=0, before this fix). An
                # off-center ratio de-syncs the split point from the root
                # after one bisection, letting it fall into a shrinking
                # box's true interior from then on.
                b1, b2 = bisect(box[widest], 0.4)
                box_a = copy(box); box_a[widest] = b1
                box_b = copy(box); box_b[widest] = b2
                push!(stack, (box_a, depth + 1))
                push!(stack, (box_b, depth + 1))
            end
        end

        if exclude !== nothing && exclude(box)
            push!(results, KrawczykResult(box, :excluded, depth))
            continue
        end

        Fbox = F(box)
        if !_all_finite(Fbox)
            subdivide_or_unresolved(:unresolved)
            continue
        end
        if any(!in_interval(0, Fbox[i]) for i in 1:length(box))
            push!(results, KrawczykResult(box, :none, depth))
            continue
        end

        K, singular = krawczyk_operator(F, box)
        if singular || !_all_finite(K)
            subdivide_or_unresolved(:unresolved)
            continue
        end
        no_root = any(isdisjoint_interval(K[i], box[i]) for i in 1:length(box))
        if no_root
            push!(results, KrawczykResult(box, :none, depth))
            continue
        end
        unique_root = all(isinterior(K[i], box[i]) for i in 1:length(box))
        if unique_root
            refined = [intersect_interval(K[i], box[i]) for i in 1:length(box)]
            push!(results, KrawczykResult(refined, :unique, depth))
            continue
        end
        subdivide_or_unresolved(:unresolved)
    end
    return results
end
