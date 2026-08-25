"""
fold_local.jl  (§8)

Fold nondegeneracy certificate at (r_+, 0, lambda_+):

    U_r  = 0       (should hold, by construction of lambda_+ = lambda_onaxis(r_+))
    U_rr = 0       (should hold, since r_+ is a critical point of lambda_onaxis(r);
                     see derivation note below)
    CERT_URLAMBDA_NONZERO:  0 not in U_rlambda enclosure
    CERT_URRR_NONZERO:      0 not in U_rrr enclosure

r_+ and lambda_+ are recomputed independently here (closed-form radical +
on-axis stationarity formula), not copied from constants_arb.jl's output.

Two INDEPENDENT derivations of U_rlambda, U_rrr are cross-checked against
each other before any sign certificate is trusted:
  (a) onaxis_r_derivs      — hand-derived analytic closed form (mp_model.jl)
  (b) onaxis_r_derivs_ad   — nested ForwardDiff automatic differentiation
If these two enclosures do not overlap, that indicates a bug in one of the
derivations and this script reports FAIL rather than picking one.

Run in BOTH IntervalArithmetic.jl (256-bit BigFloat) and Arblib.jl
(256-bit Arb ball) — two independent arithmetic libraries, since
fold-nondegeneracy values are reported as Arb balls, not point floats.

Why U_rr=0 at the fold (not just U_r=0): on the on-axis curve lambda(r) =
lambda_onaxis(r), implicit differentiation of U_r(r,lambda(r)) = 0 gives
U_rr + U_rlambda * lambda'(r) = 0, i.e. lambda'(r) = -U_rr/U_rlambda. Since
U_rlambda = 1/(1-r)^2 is never 0 for r<1, lambda'(r_+) = 0 (r_+ is by
definition a critical point / fold of the r -> lambda_onaxis(r) curve, the
quadratic 4r^2+7r-1=0 that r_+ solves) forces U_rr(r_+,lambda_+) = 0.
"""

isdefined(Main, :ValidatedMP) || include(joinpath(@__DIR__, "..", "src", "ValidatedMP.jl"))
using .ValidatedMP
using IntervalArithmetic
using Arblib
using Dates

# ---- Arb-side sign certification helper (ball -> strict sign verdict) ----
"Strict: whole Arb ball lies in (0,infty)."
arb_pos(x::Arb) = (BigFloat(Arblib.midref(x)) - BigFloat(Arblib.radius(Arb, x))) > 0
"Strict: whole Arb ball lies in (-infty,0)."
arb_neg(x::Arb) = (BigFloat(Arblib.midref(x)) + BigFloat(Arblib.radius(Arb, x))) < 0
arb_nonzero(x::Arb) = arb_pos(x) || arb_neg(x)
arb_str(x::Arb) = "[$(BigFloat(Arblib.midref(x)) - BigFloat(Arblib.radius(Arb,x))), $(BigFloat(Arblib.midref(x)) + BigFloat(Arblib.radius(Arb,x)))]"

function run_interval_path(io)
    setprecision(BigFloat, 256)
    T = BigFloat
    r_plus = (interval(T(-7)) + sqrt(interval(T(33)))) / interval(T(8))
    lambda_plus = lambda_onaxis(r_plus)

    an = onaxis_r_derivs(r_plus, lambda_plus)
    ad = onaxis_r_derivs_ad(r_plus, lambda_plus)

    println(io, "--- IntervalArithmetic.jl path (256-bit BigFloat) ---")
    println(io, "r_plus      = ", interval_str(r_plus))
    println(io, "lambda_plus = ", interval_str(lambda_plus))
    println(io)
    println(io, "U_r    (analytic) = ", interval_str(an.Ur), "   (AD) = ", interval_str(ad.Ur))
    println(io, "U_rr   (analytic) = ", interval_str(an.Urr), "   (AD) = ", interval_str(ad.Urr))
    println(io, "U_rlambda (analytic) = ", interval_str(an.Urlambda), "   (AD) = ", interval_str(ad.Urlambda))
    println(io, "U_rrr  (analytic) = ", interval_str(an.Urrr), "   (AD) = ", interval_str(ad.Urrr))
    println(io)

    ov_ur = overlaps_iv(an.Ur, ad.Ur)
    ov_urr = overlaps_iv(an.Urr, ad.Urr)
    ov_urlam = overlaps_iv(an.Urlambda, ad.Urlambda)
    ov_urrr = overlaps_iv(an.Urrr, ad.Urrr)
    println(io, "analytic/AD overlap: U_r=$(fmt(ov_ur)) U_rr=$(fmt(ov_urr)) U_rlambda=$(fmt(ov_urlam)) U_rrr=$(fmt(ov_urrr))")

    ur_zero = in_interval(0, an.Ur) && in_interval(0, ad.Ur)
    urr_zero = in_interval(0, an.Urr) && in_interval(0, ad.Urr)
    println(io, "0 in_interval U_r (both derivations)?  ", fmt(ur_zero))
    println(io, "0 in_interval U_rr (both derivations)? ", fmt(urr_zero))

    cert_urlambda = signcheck_pos(an.Urlambda) && signcheck_pos(ad.Urlambda)
    cert_urrr = signcheck_nonzero(an.Urrr) && signcheck_nonzero(ad.Urrr)
    println(io, "CERT_URLAMBDA_NONZERO = ", fmt(cert_urlambda))
    println(io, "CERT_URRR_NONZERO     = ", fmt(cert_urrr))
    println(io)

    overall = ov_ur && ov_urr && ov_urlam && ov_urrr && ur_zero && urr_zero && cert_urlambda && cert_urrr
    return (pass = overall, r_plus = r_plus, lambda_plus = lambda_plus, an = an, ad = ad)
end

function run_arb_path(io)
    setprecision(Arb, 256)
    r_plus = (-7 + sqrt(Arb(33))) / 8
    lambda_plus = begin
        a = r_plus^2 + r_plus + 1
        (2r_plus + 1) * (1 - r_plus)^2 / (a * sqrt(a))
    end

    a = r_plus^2 + r_plus + 1
    sqa = sqrt(a)
    b = 2r_plus + 1
    q = 1 - r_plus
    Ur = -b / (a * sqa) + lambda_plus / q^2
    Urr = -2 / (a * sqa) + Arb(3 // 2) * b^2 / (a^2 * sqa) + 2 * lambda_plus / q^3
    Urlambda = 1 / q^2
    Urrr = 9 * b / (a^2 * sqa) - Arb(15 // 4) * b^3 / (a^3 * sqa) + 6 * lambda_plus / q^4

    println(io, "--- Arblib.jl path (256-bit Arb ball, analytic closed form) ---")
    println(io, "r_plus      = ", r_plus)
    println(io, "lambda_plus = ", lambda_plus)
    println(io, "U_r    = ", Ur, "   as [lo,hi] = ", arb_str(Ur))
    println(io, "U_rr   = ", Urr, "   as [lo,hi] = ", arb_str(Urr))
    println(io, "U_rlambda = ", Urlambda, "   as [lo,hi] = ", arb_str(Urlambda))
    println(io, "U_rrr  = ", Urrr, "   as [lo,hi] = ", arb_str(Urrr))

    ur_zero = Arblib.contains(Ur, Arb(0))
    urr_zero = Arblib.contains(Urr, Arb(0))
    cert_urlambda = arb_pos(Urlambda)
    cert_urrr = arb_nonzero(Urrr)
    println(io, "0 contained in U_r ball?  ", fmt(ur_zero))
    println(io, "0 contained in U_rr ball? ", fmt(urr_zero))
    println(io, "CERT_URLAMBDA_NONZERO (Arb) = ", fmt(cert_urlambda))
    println(io, "CERT_URRR_NONZERO (Arb)     = ", fmt(cert_urrr))
    println(io)

    overall = ur_zero && urr_zero && cert_urlambda && cert_urrr
    return (pass = overall, r_plus = r_plus, lambda_plus = lambda_plus,
            Ur = Ur, Urr = Urr, Urlambda = Urlambda, Urrr = Urrr)
end

function run_fold_local()
    io = IOBuffer()
    println(io, "=== Fold nondegeneracy certificate (sec 8) ===")
    println(io, "Julia version: ", VERSION)
    println(io, "Run date: ", Dates.now())
    println(io)

    iv_res = run_interval_path(io)
    arb_res = run_arb_path(io)

    # cross-check IntervalArithmetic vs Arblib on the headline nonzero quantities
    iv_urlam_lo, iv_urlam_hi = inf(iv_res.an.Urlambda), sup(iv_res.an.Urlambda)
    arb_urlam_lo = BigFloat(Arblib.midref(arb_res.Urlambda)) - BigFloat(Arblib.radius(Arb, arb_res.Urlambda))
    arb_urlam_hi = BigFloat(Arblib.midref(arb_res.Urlambda)) + BigFloat(Arblib.radius(Arb, arb_res.Urlambda))
    cross_urlam = iv_urlam_lo <= arb_urlam_hi && arb_urlam_lo <= iv_urlam_hi

    iv_urrr_lo, iv_urrr_hi = inf(iv_res.an.Urrr), sup(iv_res.an.Urrr)
    arb_urrr_lo = BigFloat(Arblib.midref(arb_res.Urrr)) - BigFloat(Arblib.radius(Arb, arb_res.Urrr))
    arb_urrr_hi = BigFloat(Arblib.midref(arb_res.Urrr)) + BigFloat(Arblib.radius(Arb, arb_res.Urrr))
    cross_urrr = iv_urrr_lo <= arb_urrr_hi && arb_urrr_lo <= iv_urrr_hi

    println(io, "=== Cross-library consistency (IntervalArithmetic vs Arblib) ===")
    println(io, "U_rlambda overlap: ", fmt(cross_urlam))
    println(io, "U_rrr overlap:     ", fmt(cross_urrr))

    # reference-value sanity (reference values ~0.74711263784 and ~10.6199070453 — not
    # hardcoded into the certificate logic, only reported here for comparison).
    # NOTE: these references have only 11-12 significant digits, i.e. an implied
    # tolerance (~1e-11/1e-12) many orders wider than our ~1e-76 enclosure, so we
    # check agreement to within that implied tolerance (full ULP of the last
    # printed digit) rather than strict interval containment, which would
    # spuriously "fail" on a correct, far tighter enclosure (see constants_arb.jl
    # for the same issue encountered with r_-, r_+, lambda_-, lambda_+).
    function agrees_to_printed_digits(mid::BigFloat, ref_str::String)
        dot = findfirst('.', ref_str)
        ndig = dot === nothing ? 0 : length(ref_str) - dot
        ulp = big(10.0)^(-ndig)
        return abs(mid - parse(BigFloat, ref_str)) <= ulp
    end
    ref_urlam = "0.74711263784"
    ref_urrr = "10.6199070453"
    mid_urlam = (inf(iv_res.an.Urlambda) + sup(iv_res.an.Urlambda)) / 2
    mid_urrr = (inf(iv_res.an.Urrr) + sup(iv_res.an.Urrr)) / 2
    println(io)
    println(io, "reference U_rlambda ~ $(ref_urlam) : agrees to printed precision? ",
            fmt(agrees_to_printed_digits(mid_urlam, ref_urlam)))
    println(io, "reference U_rrr ~ $(ref_urrr) : agrees to printed precision? ",
            fmt(agrees_to_printed_digits(mid_urrr, ref_urrr)))

    overall = iv_res.pass && arb_res.pass && cross_urlam && cross_urrr
    println(io)
    println(io, "CERT_FOLD_NONDEGENERACY = ", overall ? "PASS" : "FAIL")

    out = String(take!(io))
    print(out)

    certs = [
        Certificate("fold_local_interval256", "U_rlambda", string(inf(iv_res.an.Urlambda)), string(sup(iv_res.an.Urlambda)),
                    iv_res.pass ? "PASS" : "FAIL", "nonzero required"),
        Certificate("fold_local_interval256", "U_rrr", string(inf(iv_res.an.Urrr)), string(sup(iv_res.an.Urrr)),
                    iv_res.pass ? "PASS" : "FAIL", "nonzero required"),
        Certificate("fold_local_arb256", "U_rlambda", arb_str(arb_res.Urlambda), "", arb_res.pass ? "PASS" : "FAIL", "nonzero required"),
        Certificate("fold_local_arb256", "U_rrr", arb_str(arb_res.Urrr), "", arb_res.pass ? "PASS" : "FAIL", "nonzero required"),
    ]
    mkpath(joinpath(@__DIR__, "..", "output"))
    write_certificates_csv(joinpath(@__DIR__, "..", "output", "fold_certificate.csv"), certs)
    open(joinpath(@__DIR__, "..", "output", "fold_local_log.txt"), "w") do f
        write(f, out)
    end

    return overall
end

if abspath(PROGRAM_FILE) == @__FILE__
    run_fold_local()
end
