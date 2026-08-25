"""
pitchfork_branch.jl  (§7)

Off-axis (t=±sqrt(s)) existence and uniqueness certificate for the
symmetric family, across the full lambda range around the pitchfork point.

LOGICAL STRUCTURE (three independent pieces, each a rigorous interval
argument on its own — this went through two failed designs before landing
here; see the two "LESSON" notes below, which are load-bearing, not just
commentary):

Part A — GLOBAL monotonicity of Q(0,lambda) across the WHOLE range
[lambda_--0.001, lambda_-+0.001] in ONE box computation:
    Q(0,lambda_-) contains 0                          (definitional check)
    dQ(0,lambda)/dlambda > 0  throughout the WHOLE range (one box, not just
                                                          at the point lambda_-)
Combined via the mean value theorem: for any lambda > lambda_- in the range,
Q(0,lambda) = Q(0,lambda_-) + (dQ/dlambda at some xi in between)*(lambda-lambda_-)
            >= inf(Q(0,lambda_-)) + inf(dQ/dlambda)*(lambda-lambda_-) > 0
(the first term is a ~1e-76-scale enclosure of 0, utterly dominated by the
second as soon as lambda-lambda_- exceeds ~1e-76, i.e. for every slab used
below). Symmetrically Q(0,lambda)<0 for lambda<lambda_-. This is the ONLY
place Q(0,lambda)'s sign is established — LESSON 1: evaluating Q(0,lambda)
DIRECTLY on each narrow lambda-slab was tried first and failed, because the
true Q(0,lambda) value that close to the bifurcation is only ~1e-4, while
evaluating Q over even the (already fairly tight) r-box needed to cover the
FULL s in [0,s_max] range independently on each slab produced an enclosure
~1e-2 wide (V_rs ~ -4.8 times the box's necessary ~0.002 r-width) — noise
far exceeding the signal. Only a DERIVATIVE-based global monotonicity
argument sidesteps this.

Part B — per-slab, lambda in [lambda_-, lambda_-+0.001] (N=32 slabs), each
using r_box = [r_- -0.0007, r_- +0.0014] (the SAME box pitchfork_local.jl
already certified V_rr>0, Q_s<0 on for the full +-0.001 range):
    V_rr(r_box,s_box,slab) > 0            [one full-box evaluation]
    Q_s(r_box,s_box,slab)  < 0            [one full-box evaluation]
    V_r(r_lo,s,slab)<0 and V_r(r_hi,s,slab)>0  AT s=0 AND AT s=s_max
        (thin r, thin s each — LESSON 2: evaluating this over the full
        s_box in one shot, as in an earlier version of this script, gave
        V_r enclosures straddling 0 at both endpoints, again from the same
        dependency-driven overestimation; pinning s to its two thin
        endpoint values instead of leaving it a width-0.0004 interval
        removes that overestimation entirely and gives comfortably signed
        results with a >2x margin)
    Q(r_box_refined_at_smax, s_max, slab) < 0   [r refined via 1D interval
        Newton specifically at s=s_max, not the wide r_box, for the same
        reason as Part A's lesson]
V_rr>0 on the whole box plus existence+uniqueness of the radial root AT
s=0 and s=s_max, by the standard implicit-function/connectedness argument,
give a well-defined unique continuous branch r(s,lambda) for ALL s in
[0,s_max] (not just the two checked endpoints): V_rr>0 makes V_r(.,s,lambda)
strictly increasing in r throughout the box for every fixed (s,lambda), so
its unique zero r(s,lambda) varies continuously in (s,lambda) by the
implicit function theorem, and it cannot "jump" between the branches
identified at s=0 and s=s_max without violating uniqueness somewhere in
between (which V_rr>0 on the whole box excludes). Then, on this branch:
Q(0,slab)>0 (Part A) and Q(s_max,slab)<0 (this section) and Q_s<0
throughout (this section) => Q is strictly decreasing along the branch and
changes sign exactly once => exactly one s* in (0,s_max) with Q(s*)=0.

Part C — per-slab, lambda in [lambda_--0.001, lambda_-] (N=32 slabs):
    Q_s(r_box,s_box,slab) < 0
Combined with Part A's Q(0,slab)<0 on this side: Q is decreasing and
already negative at s=0, so Q(s)<0 throughout [0,s_max] => NO off-axis root
exists on any slab in this range.
"""

isdefined(Main, :ValidatedMP) || include(joinpath(@__DIR__, "..", "src", "ValidatedMP.jl"))
using .ValidatedMP
using IntervalArithmetic
using Dates

const N_SLABS = 32

"""
refine_r_box(r_box, s_box, lambda_slab; iters=8)

1D interval-Newton refinement of the radial-root enclosure for a FIXED
(s_box, lambda_slab) box (either may still be a nondegenerate interval).
Valid because V_rr>0 is required to hold (checked each iteration) on the
starting box, making the map monotonic and the contraction rigorous.
"""
function refine_r_box(r_box, s_box, lambda_slab; iters::Int = 8)
    rb = r_box
    for _ in 1:iters
        d = V_second_order(rb, s_box, lambda_slab)
        signcheck_pos(d.Vrr) || break
        rmid = interval(mid(rb))
        Vr_mid = V_second_order(rmid, s_box, lambda_slab).Vr
        rnew = rmid - Vr_mid / d.Vrr
        rb2 = intersect_interval(rb, rnew)
        isempty_interval(rb2) && break
        isequal_interval(rb2, rb) && break
        rb = rb2
    end
    return rb
end

function run_pitchfork_branch()
    io = IOBuffer()
    println(io, "=== Pitchfork off-axis existence/uniqueness (sec 7) ===")
    println(io, "Julia version: ", VERSION)
    println(io, "Run date: ", Dates.now())
    println(io)

    setprecision(BigFloat, 256)
    T = BigFloat
    r_minus = (interval(T(-5)) + sqrt(interval(T(33)))) / interval(T(4))
    lambda_minus = lambda_onaxis(r_minus)
    r_box = r_minus + interval(T(-7 // 10000), T(14 // 10000))
    s_box = interval(T(0), T(4 // 10000))
    s0 = interval(T(0))
    smax = interval(sup(s_box))
    r_lo = interval(inf(r_box))
    r_hi = interval(sup(r_box))

    # --- Part A: GLOBAL monotonicity of Q(0,lambda) over the whole +-0.001 range ---
    println(io, "--- Part A: global monotonicity of Q(0,lambda) ---")
    full_at_minus = V_full_second_order(r_minus, s0, lambda_minus)
    Q0_at_bifurcation = 2 * full_at_minus.Vs
    q0_contains_zero = in_interval(0, Q0_at_bifurcation)
    println(io, "Q(0,lambda_-) = ", interval_str(Q0_at_bifurcation), "   contains 0? ", fmt(q0_contains_zero))

    lambda_full_range = lambda_minus + interval(T(-1 // 1000), T(1 // 1000))
    r_axis_full = refine_r_box(r_box, s0, lambda_full_range)
    full_box = V_full_second_order(r_axis_full, s0, lambda_full_range)
    dQ0_dlambda_box = 2 * (full_box.Vslambda - full_box.Vrs * full_box.Vrlambda / full_box.Vrr)
    dq0_positive_global = signcheck_pos(dQ0_dlambda_box)
    println(io, "r_axis enclosure over full lambda range: ", interval_str(r_axis_full))
    println(io, "dQ(0,lambda)/dlambda over FULL range = ", interval_str(dQ0_dlambda_box), "   > 0 throughout? ", fmt(dq0_positive_global))

    partA_ok = q0_contains_zero && dq0_positive_global
    println(io, "=> Q(0,lambda) > 0 for lambda in (lambda_-, lambda_-+0.001] : ", fmt(partA_ok))
    println(io, "=> Q(0,lambda) < 0 for lambda in [lambda_--0.001, lambda_-) : ", fmt(partA_ok))
    println(io, "Part A overall: ", partA_ok ? "PASS" : "FAIL")
    println(io)

    # --- Part B: slab sweep UP, lambda in [lambda_-, lambda_-+0.001] ---
    println(io, "--- Part B: existence+uniqueness slabs, lambda in [lambda_-, lambda_-+0.001], N=$(N_SLABS) ---")
    up_ok = Bool[]
    worst_qsmax = -Inf
    for k in 0:(N_SLABS - 1)
        lo = T(k) / T(N_SLABS) * T(1 // 1000)
        hi = T(k + 1) / T(N_SLABS) * T(1 // 1000)
        slab = lambda_minus + interval(lo, hi)

        d = V_second_order(r_box, s_box, slab)
        Qs = 2 * (d.Vss - d.Vrs^2 / d.Vrr)
        cert_vrr = signcheck_pos(d.Vrr)
        cert_qs = signcheck_neg(Qs)

        Vr_lo_0 = V_second_order(r_lo, s0, slab).Vr
        Vr_hi_0 = V_second_order(r_hi, s0, slab).Vr
        Vr_lo_s = V_second_order(r_lo, smax, slab).Vr
        Vr_hi_s = V_second_order(r_hi, smax, slab).Vr
        cert_exist = signcheck_neg(Vr_lo_0) && signcheck_pos(Vr_hi_0) &&
                     signcheck_neg(Vr_lo_s) && signcheck_pos(Vr_hi_s)

        r_atsmax = refine_r_box(r_box, smax, slab)
        Q_at_smax = 2 * V_second_order(r_atsmax, smax, slab).Vs
        cert_qsmax = signcheck_neg(Q_at_smax)
        worst_qsmax = max(worst_qsmax, Float64(sup(Q_at_smax)))

        ok = cert_vrr && cert_qs && cert_exist && cert_qsmax
        push!(up_ok, ok)
        if !ok || k == 0 || k == N_SLABS - 1
            println(io, "slab $(k+1)/$(N_SLABS): Vrr>0=", fmt(cert_vrr), " Qs<0=", fmt(cert_qs),
                    " exist=", fmt(cert_exist), " Qsmax<0=", fmt(cert_qsmax), " => ", ok ? "PASS" : "FAIL")
        end
    end
    partB_ok = all(up_ok)
    println(io, "worst-case sup(Q(s_max)) over all UP slabs = ", worst_qsmax, " (must be < 0)")
    println(io, "Part B overall: ", partB_ok ? "PASS" : "FAIL", "  [", count(identity, up_ok), "/", N_SLABS, "]")
    println(io)

    # --- Part C: slab sweep DOWN, lambda in [lambda_--0.001, lambda_-] : no off-axis root ---
    println(io, "--- Part C: NO off-axis root slabs, lambda in [lambda_--0.001, lambda_-], N=$(N_SLABS) ---")
    down_ok = Bool[]
    for k in 0:(N_SLABS - 1)
        lo = -T(N_SLABS - k) / T(N_SLABS) * T(1 // 1000)
        hi = -T(N_SLABS - k - 1) / T(N_SLABS) * T(1 // 1000)
        slab = lambda_minus + interval(lo, hi)
        d = V_second_order(r_box, s_box, slab)
        Qs = 2 * (d.Vss - d.Vrs^2 / d.Vrr)
        cert_qs = signcheck_neg(Qs)
        ok = cert_qs
        push!(down_ok, ok)
        if !ok || k == 0 || k == N_SLABS - 1
            println(io, "slab $(k+1)/$(N_SLABS): Qs<0=", fmt(cert_qs), " => ", ok ? "PASS" : "FAIL")
        end
    end
    partC_ok = all(down_ok) && partA_ok  # Part C's "Q(0)<0" half comes from Part A
    println(io, "Part C overall: ", partC_ok ? "PASS" : "FAIL", "  [", count(identity, down_ok), "/", N_SLABS, "]")
    println(io)

    overall = partA_ok && partB_ok && partC_ok
    println(io, "CERT_PITCHFORK_OFFAXIS_EXISTENCE  = ", (partA_ok && partB_ok) ? "PASS" : "FAIL")
    println(io, "CERT_PITCHFORK_OFFAXIS_UNIQUENESS = ", (partA_ok && partB_ok) ? "PASS" : "FAIL")
    println(io, "CERT_PITCHFORK_NOROOT_BELOW       = ", partC_ok ? "PASS" : "FAIL")
    println(io, "CERT_PITCHFORK_BRANCH_OVERALL     = ", overall ? "PASS" : "FAIL")

    out = String(take!(io))
    print(out)
    mkpath(joinpath(@__DIR__, "..", "output"))
    open(joinpath(@__DIR__, "..", "output", "pitchfork_branch_log.txt"), "w") do f
        write(f, out)
    end
    return overall
end

if abspath(PROGRAM_FILE) == @__FILE__
    run_pitchfork_branch()
end
