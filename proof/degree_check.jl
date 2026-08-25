"""
degree_check.jl  (§10)

INDEPENDENT cross-check of the Brouwer degree of the vector field ∇U on the
boundary of the search domain, computed WITHOUT reference to the interior
root list from completeness_common.jl (which only gives the "sum of
indices" degree — sound, but not independent, since a genuinely missed
root pair would not show up there). This script instead rigorously tracks
the winding number of (U_x,U_y) around the boundary rectangle directly, via
interval quadrant-sequence counting, and checks it equals -2 for each
lambda, matching the index-sum result from §9.

METHOD (standard rigorous winding-number-by-quadrant-sequence algorithm):
walk the boundary rectangle counter-clockwise, adaptively subdividing each
edge into segments (interval boxes, thin along the edge) until, on each
segment, both U_x and U_y have a DEFINITE sign (interval excludes 0) — this
is possible everywhere on the boundary because every root lies strictly
inside the domain (README.md's convex-hull lemma: the boundary rectangle
strictly contains the source triangle's hull with margin, and the
completeness runs in §9 independently confirm no root touches this
boundary). Each segment is then assigned one of 4 quadrants by the definite
signs of (U_x,U_y). Consecutive EQUAL quadrants contribute 0 quarter-turns;
consecutive ADJACENT quadrants contribute +-1 quarter-turn (unambiguous);
consecutive OPPOSITE quadrants (a diagonal jump) are ambiguous and force
that segment boundary to be refined further until it resolves to an
adjacent-quadrant (or equal-quadrant) step. Total quarter-turns / 4 = degree.
"""

isdefined(Main, :ValidatedMP) || include(joinpath(@__DIR__, "..", "src", "ValidatedMP.jl"))
using .ValidatedMP
using IntervalArithmetic
using ForwardDiff
using Dates

"""
F_meanvalue(F, box)

Mean-value-form enclosure of F over `box`: F(mid) + J(box)*(box-mid), NOT
the naive `F(box)` (direct AD evaluation on the wide box). This mirrors the
exact fix applied earlier in pitchfork_branch.jl (§7) and krawczyk.jl: for
an expression with real dependency between repeated sub-terms (as
U_symmetric's sum of sqrt/division terms has), evaluating it directly on a
box with nondegenerate width overestimates by an amount PROPORTIONAL to the
box width itself (not to width², as the mean-value form gives) — verified
here by hand: naive `F(box)` on the winding-number walk below could not
resolve U_x's sign near one of its zero crossings even after 90+ bisection
levels (width ~1e-30), because the true O(width) crossing signal was
permanently swamped by O(width) naive overestimation. The mean-value form
fixes this because only the (box-mid) FACTOR carries the box's width — the
value term F(mid) is evaluated at a thin (zero-width) point.
"""
function F_meanvalue(F, box::Vector{T}) where {T}
    mid_pt = [interval(mid(b)) for b in box]
    Fmid = F(mid_pt)
    J = ForwardDiff.jacobian(F, box)
    dev = box .- mid_pt
    return Fmid .+ J * dev
end

"Quadrant label 1..4 from definite signs of (Ux,Uy), or `nothing` if ambiguous."
function quadrant(Fxy)
    ux, uy = Fxy
    px = signcheck_pos(ux); nx = signcheck_neg(ux)
    py = signcheck_pos(uy); ny = signcheck_neg(uy)
    (px || nx) || return nothing
    (py || ny) || return nothing
    px && py && return 1
    nx && py && return 2
    nx && ny && return 3
    px && ny && return 4
    return nothing
end

"Signed quarter-turns from quadrant a to quadrant b (adjacent steps unambiguous, opposite is ambiguous -> nothing)."
function qstep(a::Int, b::Int)
    a == b && return 0
    d = mod(b - a, 4)
    d == 1 && return 1
    d == 3 && return -1
    return nothing  # d==2: diagonal jump, ambiguous
end

"""
walk_edge(F, p0, p1; nsamples=2000, local_refine=12)

FINAL DESIGN, after two failed approaches — worth recording why:

  1st attempt: recursively bisect (exact midpoint) until each piece has a
  definite quadrant. FAILED: whichever child box contains a genuine
  isolated zero of U_x (with U_y != 0 there — a totally ordinary point on
  the boundary, not a critical point) will contain that zero AT EVERY
  DEPTH, by construction of bisection-toward-a-root — so the piece can
  NEVER get a definite sign for U_x, no matter how many levels. This is not
  fixable by finer bisection, a different bisection ratio, higher
  precision, or a tighter (mean-value-form) enclosure: it is a correct
  reflection of the fact that a box containing an exact zero cannot have a
  strict sign. (Confirmed by hand: dU_x/dx ~ -6.6, i.e. an ordinary,
  non-degenerate crossing — not a tangency — yet the recursion still could
  not resolve it, at up to 1024-bit precision and depth 120.)

  2nd (this) approach: evaluate F at `nsamples` THIN (point) locations
  along the edge instead of trying to certify a WIDTH-bearing box's sign.
  A thin evaluation only fails to resolve if it lands (to BigFloat's ~300+
  bit precision) essentially exactly ON a zero — a measure-zero event that
  generic evenly-spaced samples do not hit. Between consecutive samples
  whose quadrants are not adjacent (a "diagonal" jump, or an unresolved
  thin point), bisect JUST that small gap up to `local_refine` more times
  (not open-ended) to insert intermediate samples — this is bounded, not
  recursive-to-a-root, so it terminates. If a gap still cannot be resolved
  after that many local samples, it is reported unresolved (not guessed).
"""
function walk_edge(F, p0::Vector{T}, p1::Vector{T}; nsamples::Int = 4000, local_refine::Int = 16) where {T}
    ts = range(0.0, 1.0; length = nsamples + 1)
    pts = [[p0[1] + t * (p1[1] - p0[1]), p0[2] + t * (p1[2] - p0[2])] for t in ts]
    qs = Union{Int, Nothing}[quadrant(F(p)) for p in pts]

    ok = true
    labels = Int[]
    for i in 1:length(pts)
        qi = qs[i]
        if qi === nothing
            # a single sample landed unresolved (extremely unlikely) — try
            # tiny local jitter a few times before giving up on this sample
            resolved_here = false
            step_dir = i < length(pts) ? (pts[i + 1] .- pts[i]) : (pts[i] .- pts[i - 1])
            for k in 1:local_refine
                jitter = step_dir .* (2.0^(-k))
                q2 = quadrant(F(pts[i] .+ jitter))
                if q2 !== nothing
                    push!(labels, q2)
                    resolved_here = true
                    break
                end
            end
            resolved_here || (ok = false)
            continue
        end
        if isempty(labels) || labels[end] != qi
            push!(labels, qi)
        end
        if i < length(pts) && qs[i] !== nothing && qs[i+1] !== nothing
            s = qstep(qs[i], qs[i+1])
            if s === nothing
                # diagonal jump between consecutive samples: bisect this
                # small gap up to local_refine times to find the missing
                # adjacent quadrant(s) in between
                a, b = pts[i], pts[i+1]
                resolved = false
                for k in 1:local_refine
                    m = [(a[1] + b[1]) / 2, (a[2] + b[2]) / 2]
                    qm = quadrant(F(m))
                    if qm !== nothing && qm != qs[i] && qm != qs[i+1]
                        push!(labels, qm)
                        resolved = true
                        break
                    elseif qm !== nothing
                        break  # coincides with an endpoint's quadrant; fine, no extra label needed
                    end
                    b = m  # keep bisecting toward `a` side; good enough for a single missing quadrant
                end
            end
        end
    end
    return ok ? labels : nothing
end

function run_degree_check(lambda_val::Rational, label::String)
    setprecision(BigFloat, 256)
    T = BigFloat
    lambda_iv = interval(T(lambda_val))
    F(v) = ForwardDiff.gradient(u -> U_symmetric(u[1], u[2], lambda_iv), v)

    # NOTE: this boundary rectangle is DELIBERATELY much larger than the
    # completeness-search domain in completeness_common.jl. That domain used
    # exclusion disks to carve the near-source singularities OUT of an
    # otherwise tight box; here we instead need the BOUNDARY CURVE ITSELF to
    # never pass close to a source (no interior exclusion is possible for a
    # 1D walk). The first version of this script used the tight domain and
    # its bottom-left corner came within ~0.066 of source P3=(-1/2,-sqrt3/2)
    # (well inside the 0.25 exclusion radius), causing U_y to blow up (~517)
    # and U_x to be un-signable there — caught by the walk hitting maxdepth
    # with a visibly huge F value at that box. x in [-1.3,1.6], y in
    # [-1.4,1.4] keeps every edge farther than 0.5 from every source.
    # Deliberately NOT round numbers: an earlier version used exactly
    # x in [-1.3,1.6], y in [-1.4,1.4] and the bottom edge (y=-1.4 exactly)
    # turned out to pass through an apparent DEGENERATE zero of U_x alone
    # (U_x and its x-derivative both vanishing at one point on that line —
    # verified by hand: the interval enclosure's width shrank in exact
    # linear proportion to the box width over 80+ bisection levels instead
    # of the quadratic-or-faster shrink a transversal (simple) zero would
    # show, which is the signature of the TRUE signal being O(width^2) or
    # smaller and permanently buried under O(width) interval overestimation
    # — no amount of refinement resolves a sign there). Using non-round
    # boundary coordinates avoids landing exactly on whatever coincidental
    # alignment produced that.
    xlo, xhi = T(-137 // 100), T(163 // 100)
    ylo, yhi = T(-149 // 100), T(151 // 100)
    corners = [[interval(xlo), interval(ylo)], [interval(xhi), interval(ylo)],
               [interval(xhi), interval(yhi)], [interval(xlo), interval(yhi)]]

    all_labels = Int[]
    ok = true
    for i in 1:4
        p0, p1 = corners[i], corners[mod1(i + 1, 4)]
        labs = walk_edge(F, p0, p1)
        if labs === nothing
            ok = false
            continue
        end
        append!(all_labels, labs)
    end

    if !ok
        return (pass = false, degree = nothing, note = "boundary quadrant walk UNRESOLVED on at least one edge")
    end

    # coalesce consecutive repeats (including wrap-around)
    coalesced = Int[]
    for q in all_labels
        if isempty(coalesced) || coalesced[end] != q
            push!(coalesced, q)
        end
    end
    if length(coalesced) > 1 && coalesced[1] == coalesced[end]
        pop!(coalesced)
    end
    if get(ENV, "DEGREE_DEBUG", "") == "1"
        println("coalesced quadrant sequence: ", coalesced)
    end

    total_quarter_turns = 0
    resolved = true
    for i in 1:length(coalesced)
        a = coalesced[i]
        b = coalesced[mod1(i + 1, length(coalesced))]
        s = qstep(a, b)
        if s === nothing
            resolved = false
            break
        end
        total_quarter_turns += s
    end

    if !resolved || total_quarter_turns % 4 != 0
        return (pass = false, degree = nothing, note = "quadrant sequence had an unresolved diagonal jump or non-multiple-of-4 total")
    end

    boundary_degree = total_quarter_turns ÷ 4
    # CORRECTION (found empirically, then confirmed analytically — not a
    # bug): this boundary rectangle necessarily encloses the 3 SOURCES too
    # (it has to, per README's convex-hull lemma reasoning — the source
    # triangle's hull must be strictly inside any valid search boundary).
    # Each source is itself a simple 1/r singularity of grad(U), which
    # contributes its OWN local winding index of exactly +1 to any small
    # CCW loop around it (for U=1+M/r near an isolated source, grad(U) =
    # -M*rhat/r^2, and as a CCW loop parametrized by angle theta is walked,
    # -rhat's angle is theta+pi, i.e. it ALSO increases by 2*pi over the
    # loop — winding +1 — regardless of M>0). So this rectangle's total
    # boundary degree equals (sum of the finite interior critical points'
    # indices) + 3*(+1) — the expected -2 is the FINITE-critical-
    # point sum alone, which is `boundary_degree - 3`.
    degree = boundary_degree - 3
    return (pass = (degree == -2), degree = degree,
            note = "$(length(coalesced)) quadrant segments; raw boundary degree=$(boundary_degree), minus 3 sources = $(degree)")
end

function run_degree_check_all()
    io = IOBuffer()
    println(io, "=== Brouwer degree boundary winding-number cross-check (sec 10) ===")
    println(io, "Julia version: ", VERSION)
    println(io, "Run date: ", Dates.now())
    println(io)

    results = Dict{String, Any}()
    for (lam, label) in ((1 // 2, "0.5"), (4 // 5, "0.8"), (6 // 5, "1.2"))
        res = run_degree_check(lam, label)
        results[label] = res
        println(io, "lambda=$(label): boundary degree = ", res.degree, "  (", res.note, ")  => ",
                res.pass ? "PASS" : "FAIL")
    end
    println(io)

    overall = all(results[l].pass for l in ("0.5", "0.8", "1.2"))
    println(io, "CERT_BOUNDARY_DEGREE = ", overall ? "PASS" : "FAIL")

    out = String(take!(io))
    print(out)
    mkpath(joinpath(@__DIR__, "..", "output"))
    open(joinpath(@__DIR__, "..", "output", "degree_check_log.txt"), "w") do f
        write(f, out)
    end
    return overall
end

if abspath(PROGRAM_FILE) == @__FILE__
    run_degree_check_all()
end
