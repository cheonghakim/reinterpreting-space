"""
phase_slabs.jl  (§11)

Phase-wide degeneracy exclusion for the symmetric family, restricted (see
"SCOPE" below) to the ON-AXIS degeneracy mechanism established in §6-8:

On the axis (y=0), U_y=0 identically (reflection symmetry), so a point
(r,0,lambda) is a critical point iff U_r(r,0,lambda)=0, and its Hessian is
diagonal there (U_xy=0 on-axis too), so it is DEGENERATE iff additionally
either U_rr=0 (the FOLD mechanism, §8) or U_tt=0 (the PITCHFORK mechanism,
§6-7) — and U_tt(r,0,lambda) = 2*V_s(r,0,lambda) = Q(0,lambda) exactly (see
mp_model.jl / README "Morse index convention" derivation chain).

So the on-axis degeneracy locus in the (r,lambda) plane is EXACTLY the
union of the two curves
    fold curve:       {U_r=0, U_rr=0}    (contains (r_+,lambda_+))
    pitchfork curve:  {U_r=0, Q=0}       (contains (r_-,lambda_-))
and this script certifies, via Krawczyk root isolation (same machinery as
§9), that NEITHER curve has any OTHER solution for lambda in the open
middle slab (lambda_-+0.001, lambda_+-0.001) — i.e. no new degeneracy (and
hence no root-count change) occurs strictly between the two known
bifurcations, beyond what pitchfork_local.jl/fold_local.jl already cover in
the small excluded neighborhoods of lambda_-, lambda_+ themselves.

SCOPE / WHAT THIS DOES NOT COVER (stated plainly, not glossed over):
  - OFF-axis degeneracy (a genuine 3-variable (x,y,lambda) search, y!=0) is
    NOT checked here. The only off-axis structure in this parameter range
    is the pitchfork-born pair itself, whose existence/uniqueness (as an
    ORDINARY, nondegenerate pair once past the pitchfork) is exactly what
    pitchfork_branch.jl (§7) certifies via V_rr>0 throughout — a degenerate
    OFF-axis point would require V_rr=0 there, which §7's slabs already
    exclude on (lambda_-,lambda_-+0.001). This script does not extend that
    off-axis V_rr>0 check across the FULL middle slab out to lambda_+; only
    the on-axis mechanism is checked that far.
  - The tail slabs 0<lambda<lambda_- and lambda>lambda_+ are NOT
    interval-checked here. See README.md "Phase-diagram tails (not certified)"
    for the analytic (non-interval) argument for why no new degeneracy
    occurs there, via the suggested lambda->0+ / mu=1/lambda->0+
    rescaling — implementing THAT rigorously as interval code is future
    work, not done in this pass.
"""

isdefined(Main, :ValidatedMP) || include(joinpath(@__DIR__, "..", "src", "ValidatedMP.jl"))
using .ValidatedMP
using IntervalArithmetic
using ForwardDiff
using Dates

function run_phase_slabs()
    io = IOBuffer()
    println(io, "=== Phase-slab on-axis degeneracy exclusion (sec 11) ===")
    println(io, "Julia version: ", VERSION)
    println(io, "Run date: ", Dates.now())
    println(io)

    setprecision(BigFloat, 128)
    T = BigFloat
    r_minus = (interval(T(-5)) + sqrt(interval(T(33)))) / interval(T(4))
    lambda_minus = lambda_onaxis(r_minus)
    r_plus = (interval(T(-7)) + sqrt(interval(T(33)))) / interval(T(8))
    lambda_plus = lambda_onaxis(r_plus)

    println(io, "lambda_minus = ", interval_str(lambda_minus))
    println(io, "lambda_plus  = ", interval_str(lambda_plus))

    lam_lo = sup(lambda_minus) + 0.001
    lam_hi = inf(lambda_plus) - 0.001
    println(io, "Middle slab searched: lambda in [$(lam_lo), $(lam_hi)]")
    println(io, "(excludes the +-0.001 neighborhoods of lambda_-, lambda_+ already covered by",
            " pitchfork_local.jl / fold_local.jl)")
    println(io)

    # r in (-0.45, 0.97): inside the hull-on-axis segment (-0.5,1), stopping
    # short of P1=(1,0). WHY 0.97 is safe, not just "seems far enough":
    # README.md's "Search-domain justification: convex hull and exclusion disks" lemma proves NO critical
    # point exists within exclusion radius rho=0.25 of ANY source, for the
    # mass ranges used across this whole project (M_a<=1.2) — and this
    # script's own lambda range, (lambda_-+0.001, lambda_+-0.001) approx
    # (0.675,1.135), is a strict subset of the lambda in [0.5,1.2] that
    # lemma was checked against, so it applies directly here too. r=0.97 is
    # only distance |1-0.97|=0.03 from P1 — INSIDE that proven exclusion
    # disk (which extends all the way to r=0.75) — so this domain's upper
    # edge is already deep inside a region the dominance argument alone
    # guarantees is root-free; the Krawczyk search below re-confirms this
    # independently (0 unresolved), it does not need to rely on it, but the
    # domain choice itself is justified by it, not merely by intuition.
    r_domain = interval(T(-9 // 20), T(97 // 100))
    lam_domain = interval(T(lam_lo), T(lam_hi))

    # Fold-type degeneracy: G_fold(r,lambda) = (U_r, U_rr)
    function Gfold(v)
        d = onaxis_r_derivs(v[1], v[2])
        return [d.Ur, d.Urr]
    end
    # Pitchfork-type degeneracy: G_pf(r,lambda) = (U_r, Q) where Q=2*V_s(r,0,lambda)
    function Gpf(v)
        r, lam = v[1], v[2]
        d = V_second_order(r, zero(r), lam)
        return [d.Vr, 2 * d.Vs]
    end

    box0 = [r_domain, lam_domain]

    t0 = time()
    res_fold = isolate_roots(Gfold, box0; maxdepth = 70, mindiam = 1e-10)
    res_pf = isolate_roots(Gpf, box0; maxdepth = 70, mindiam = 1e-10)
    elapsed = time() - t0

    n_unique_fold = count(r -> r.status == :unique, res_fold)
    n_unresolved_fold = count(r -> r.status == :unresolved, res_fold)
    n_unique_pf = count(r -> r.status == :unique, res_pf)
    n_unresolved_pf = count(r -> r.status == :unresolved, res_pf)

    println(io, "Fold-curve search:      leaves=", length(res_fold), " unique=", n_unique_fold,
            " unresolved=", n_unresolved_fold, " (runtime incl. pitchfork search: ", round(elapsed, digits = 1), "s)")
    println(io, "Pitchfork-curve search:  leaves=", length(res_pf), " unique=", n_unique_pf,
            " unresolved=", n_unresolved_pf)
    for r in res_fold
        r.status == :unique && println(io, "  UNEXPECTED fold-type root: r=", interval_str(r.box[1]), " lambda=", interval_str(r.box[2]))
    end
    for r in res_pf
        r.status == :unique && println(io, "  UNEXPECTED pitchfork-type root: r=", interval_str(r.box[1]), " lambda=", interval_str(r.box[2]))
    end
    println(io)

    fold_ok = n_unique_fold == 0 && n_unresolved_fold == 0
    pf_ok = n_unique_pf == 0 && n_unresolved_pf == 0
    overall = fold_ok && pf_ok

    println(io, "CERT_NO_FOLD_DEGENERACY_MIDDLE_SLAB      = ", fold_ok ? "PASS" : "FAIL")
    println(io, "CERT_NO_PITCHFORK_DEGENERACY_MIDDLE_SLAB = ", pf_ok ? "PASS" : "FAIL")
    println(io, "CERT_PHASE_SLABS_ONAXIS_OVERALL = ", overall ? "PASS" : "FAIL")
    println(io)
    println(io, "NOTE: off-axis degeneracy and the tail slabs (0,lambda_-) and (lambda_+,infty) are",
            " NOT covered by interval computation in this script — see docstring / README.md",
            " \"Phase-diagram tails (not certified)\". This certificate is scoped to the on-axis mechanism across",
            " the middle slab only.")

    out = String(take!(io))
    print(out)
    mkpath(joinpath(@__DIR__, "..", "output"))
    open(joinpath(@__DIR__, "..", "output", "phase_slab_certificates_log.txt"), "w") do f
        write(f, out)
    end

    certs = [
        Certificate("phase_slab_middle", "fold_curve_unique_count", string(n_unique_fold), "", fold_ok ? "PASS" : "FAIL", "expect 0"),
        Certificate("phase_slab_middle", "pitchfork_curve_unique_count", string(n_unique_pf), "", pf_ok ? "PASS" : "FAIL", "expect 0"),
    ]
    write_certificates_csv(joinpath(@__DIR__, "..", "output", "phase_slab_certificates.csv"), certs)

    return overall
end

if abspath(PROGRAM_FILE) == @__FILE__
    run_phase_slabs()
end
