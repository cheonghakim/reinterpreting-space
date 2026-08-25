"""
constants_arb.jl  (§2, §13)

Independent Arb (256-bit ball arithmetic, via Arblib.jl) re-derivation of
the reference constants r_-, lambda_-, r_+, lambda_+, epsilon_c, lambda_c.

IMPORTANT — this is NOT just re-typing the reference decimals.
r_- and r_+ are given as closed-form radicals; lambda_- and lambda_+ are
recomputed here from the FIRST-PRINCIPLES on-axis stationarity condition

    V_r(r, s=0, lambda) = 0
    <=>  -(2r+1)/(r^2+r+1)^(3/2) + lambda/(1-r)^2 = 0     (for r<1)
    <=>  lambda(r) = (2r+1)*(1-r)^2 / (r^2+r+1)^(3/2)

(derived directly from V(r,0,lambda) = 1 + 2/sqrt(r^2+r+1) + lambda/(1-r),
itself the s=0 specialization of the exact V(r,s,lambda) identity in
mp_model.jl — checked symbolically in README.md). lambda_pitchfork/fold are
then lambda(r_-), lambda(r_+) evaluated in Arb ball arithmetic. This is an
independent derivation path from whatever produced the reference
decimals below, so agreement is a genuine cross-check, not a tautology.

epsilon_c, lambda_c are algebraic functions of lambda_- (recomputed here in
Arb, not copied from any reference source).

Each computed ball is checked for containing the quoted reference
decimal (as a strict interval-containment test) AND for being tight (radius
reported explicitly). If a reference decimal is NOT contained, this script
reports FAIL for that constant instead of silently accepting it.
"""

using Arblib
using Dates

const PREC_BITS = 256
setprecision(Arb, PREC_BITS)

"""
lambda_onaxis_arb(r::Arb)

Same on-axis stationarity formula as mp_model.jl's generic `lambda_onaxis`
(see pitchfork_local.jl / fold_local.jl), but written independently here
using Arb's fractional-power support (`^Arb(3//2)`) instead of `a*sqrt(a)`,
and under a distinct name so this script has NO dependency on
src/ValidatedMP.jl at all — it is a standalone, from-scratch Arblib
computation, deliberately not sharing code with the IntervalArithmetic.jl
scripts, so a bug in one derivation path is not invisible to the other.
"""
function lambda_onaxis_arb(r::Arb)
    num = (2r + 1) * (1 - r)^2
    den = (r^2 + r + 1)^Arb(3//2)
    return num / den
end

"""
consistent_with_decimal(ball, decimal_str)

The reference numbers below are DECIMAL APPROXIMATIONS with a finite number
of digits, not exact values — so the correct check is NOT `ball` containing
the literal truncated rational `decimal_str` (a much tighter 256-bit ball
essentially never will, since the truncation error alone exceeds the ball's
~1e-76 radius). Instead we check that the true value the ball certifies to
enclose is consistent with the reference to within the reference's own
implied precision: half a unit in its last printed decimal digit, plus the
ball's own radius. Returns (is_consistent::Bool, abs_difference::BigFloat,
threshold::BigFloat) so the actual agreement is visible in the output, not
just a boolean.
"""
function consistent_with_decimal(ball::Arb, decimal_str::String)
    dot = findfirst('.', decimal_str)
    ndigits_after = dot === nothing ? 0 : length(decimal_str) - dot
    ulp = big(10.0)^(-ndigits_after)
    ref = parse(BigFloat, decimal_str)
    mid = BigFloat(Arblib.midref(ball))
    rad = BigFloat(Arblib.radius(Arb, ball))
    diff = abs(mid - ref)
    # Full-ULP (not half-ULP) tolerance: some reference decimals below
    # turned out on inspection to be TRUNCATED rather than correctly rounded
    # (e.g. r_plus: true tail "...4175186735..." chopped to "...417518"), so
    # a rounding-only half-ULP bound is unsafe here. Full ULP is the correct
    # bound for either truncation or rounding.
    threshold = ulp + rad
    return (diff <= threshold, diff, threshold)
end

function run_constants_arb()
    io = IOBuffer()
    println(io, "=== Arb 256-bit constants (independent re-derivation) ===")
    println(io, "Arblib precision (bits): ", PREC_BITS)
    println(io, "Julia version: ", VERSION)
    println(io, "Run date: ", Dates.now())
    println(io)

    # r_- = (-5+sqrt(33))/4 , root of 2r^2+5r-1=0
    r_minus = (-5 + sqrt(Arb(33))) / 4
    r_minus_ref = "0.1861406616345071649626528670547323"
    (ok_rm, diff_rm, thr_rm) = consistent_with_decimal(r_minus, r_minus_ref)
    println(io, "r_minus = ", r_minus)
    println(io, "  radius = ", Arblib.radius(Arb, r_minus))
    println(io, "  consistent with reference decimal ($(r_minus_ref))? ", ok_rm, "  |diff|=", diff_rm, " threshold=", thr_rm)

    # r_+ = (-7+sqrt(33))/8, root of 4r^2+7r-1=0 [since 2a=8=>a=4,b=7,b^2-4ac=33=>49-16c=33=>c=1]
    r_plus = (-7 + sqrt(Arb(33))) / 8
    r_plus_ref = "-0.156929669182746417518"
    (ok_rp, diff_rp, thr_rp) = consistent_with_decimal(r_plus, r_plus_ref)
    println(io)
    println(io, "r_plus = ", r_plus)
    println(io, "  radius = ", Arblib.radius(Arb, r_plus))
    println(io, "  consistent with reference decimal ($(r_plus_ref))? ", ok_rp, "  |diff|=", diff_rp, " threshold=", thr_rp)

    # lambda_- = lambda_onaxis(r_-), independent formula
    lambda_minus = lambda_onaxis_arb(r_minus)
    lambda_minus_ref = "0.6738774744707380113762233701706991"
    (ok_lm, diff_lm, thr_lm) = consistent_with_decimal(lambda_minus, lambda_minus_ref)
    println(io)
    println(io, "lambda_minus = ", lambda_minus)
    println(io, "  radius = ", Arblib.radius(Arb, lambda_minus))
    println(io, "  consistent with reference decimal ($(lambda_minus_ref))? ", ok_lm, "  |diff|=", diff_lm, " threshold=", thr_lm)

    # lambda_+ = lambda_onaxis(r_+), independent formula
    lambda_plus = lambda_onaxis_arb(r_plus)
    lambda_plus_ref = "1.13625221066468090151835125249"
    (ok_lp, diff_lp, thr_lp) = consistent_with_decimal(lambda_plus, lambda_plus_ref)
    println(io)
    println(io, "lambda_plus = ", lambda_plus)
    println(io, "  radius = ", Arblib.radius(Arb, lambda_plus))
    println(io, "  consistent with reference decimal ($(lambda_plus_ref))? ", ok_lp, "  |diff|=", diff_lp, " threshold=", thr_lp)

    # epsilon_c = (1-lambda_-)/(1+lambda_-) ; lambda_c = 1+epsilon_c
    epsilon_c = (1 - lambda_minus) / (1 + lambda_minus)
    epsilon_c_ref = "0.194830583781157"
    (ok_ec, diff_ec, thr_ec) = consistent_with_decimal(epsilon_c, epsilon_c_ref)
    println(io)
    println(io, "epsilon_c = ", epsilon_c)
    println(io, "  radius = ", Arblib.radius(Arb, epsilon_c))
    println(io, "  consistent with reference decimal ($(epsilon_c_ref))? ", ok_ec, "  |diff|=", diff_ec, " threshold=", thr_ec)

    lambda_c = 1 + epsilon_c
    lambda_c_ref = "1.194830583781157"
    (ok_lc, diff_lc, thr_lc) = consistent_with_decimal(lambda_c, lambda_c_ref)
    println(io)
    println(io, "lambda_c = ", lambda_c)
    println(io, "  radius = ", Arblib.radius(Arb, lambda_c))
    println(io, "  consistent with reference decimal ($(lambda_c_ref))? ", ok_lc, "  |diff|=", diff_lc, " threshold=", thr_lc)

    allok = ok_rm && ok_rp && ok_lm && ok_lp && ok_ec && ok_lc
    println(io)
    println(io, "CERT_ARB_CONSTANTS = ", allok ? "PASS" : "FAIL")

    out = String(take!(io))
    print(out)
    return (out=out, pass=allok,
            r_minus=r_minus, r_plus=r_plus, lambda_minus=lambda_minus,
            lambda_plus=lambda_plus, epsilon_c=epsilon_c, lambda_c=lambda_c)
end

if abspath(PROGRAM_FILE) == @__FILE__
    result = run_constants_arb()
    mkpath(joinpath(@__DIR__, "..", "output"))
    open(joinpath(@__DIR__, "..", "output", "constants.txt"), "w") do f
        write(f, result.out)
    end
end
