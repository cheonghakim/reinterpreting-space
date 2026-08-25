"""
pitchfork_local.jl  (§6)

Local pitchfork certificate: on the explicit box

    r in [r_- - 0.0007, r_- + 0.0014]
    s in [0, 0.0004]
    lambda in [lambda_- - 0.001, lambda_- + 0.001]

certify

    CERT_VRR_POSITIVE:  inf(V_rr) > 0   over the WHOLE box
    CERT_QS_NEGATIVE:   sup(Q_s)  < 0   over the WHOLE box

where Q_s = 2*(V_ss - V_sr^2/V_rr), evaluated directly by outward-rounded
interval arithmetic (via ForwardDiff automatic differentiation run ON
Interval{BigFloat} arguments — see src/derivatives.jl for why this is sound).

r_- and lambda_- are NOT hardcoded: they are rebuilt in interval arithmetic
here from the closed-form radical / on-axis-stationarity formulas, using
IntervalArithmetic.jl — independently of the Arblib computation in
constants_arb.jl. Box literal offsets (0.0007, 0.0014, 0.0004, 0.001) are
built from exact Rationals converted to BigFloat, so a single BigFloat
rounding (not a Float64-literal-then-BigFloat double rounding) sets the box
endpoints.

Run at BOTH 128-bit and 256-bit BigFloat precision; both must independently
certify, and the certified sign must agree between them (a real consistency
check, not just running twice for show).
"""

isdefined(Main, :ValidatedMP) || include(joinpath(@__DIR__, "..", "src", "ValidatedMP.jl"))
using .ValidatedMP
using IntervalArithmetic
using Dates

function run_at_precision(prec_bits::Int)
    setprecision(BigFloat, prec_bits)
    T = BigFloat

    r_minus = (interval(T(-5)) + sqrt(interval(T(33)))) / interval(T(4))
    lambda_minus = lambda_onaxis(r_minus)

    r_box = r_minus + interval(T(-7 // 10000), T(14 // 10000))
    s_box = interval(T(0), T(4 // 10000))
    lambda_box = lambda_minus + interval(T(-1 // 1000), T(1 // 1000))

    d = V_second_order(r_box, s_box, lambda_box)
    Vrr = d.Vrr
    Vrs = d.Vrs
    Vss = d.Vss
    Qs = 2 * (Vss - Vrs^2 / Vrr)

    cert_vrr = signcheck_pos(Vrr)
    cert_qs = signcheck_neg(Qs)

    return (prec_bits = prec_bits,
            r_minus = r_minus, lambda_minus = lambda_minus,
            r_box = r_box, s_box = s_box, lambda_box = lambda_box,
            Vrr = Vrr, Vrs = Vrs, Vss = Vss, Qs = Qs,
            cert_vrr = cert_vrr, cert_qs = cert_qs)
end

function run_pitchfork_local()
    io = IOBuffer()
    println(io, "=== Pitchfork local certificate (sec 6) ===")
    println(io, "Julia version: ", VERSION)
    println(io, "Run date: ", Dates.now())
    println(io)

    results = Dict{Int, Any}()
    for prec in (128, 256)
        res = run_at_precision(prec)
        results[prec] = res
        println(io, "--- precision = $(prec) bits (BigFloat) ---")
        println(io, "r_minus     = ", interval_str(res.r_minus))
        println(io, "lambda_minus= ", interval_str(res.lambda_minus))
        println(io, "r_box       = ", interval_str(res.r_box))
        println(io, "s_box       = ", interval_str(res.s_box))
        println(io, "lambda_box  = ", interval_str(res.lambda_box))
        println(io, "V_rr = ", interval_str(res.Vrr))
        println(io, "V_rs = ", interval_str(res.Vrs))
        println(io, "V_ss = ", interval_str(res.Vss))
        println(io, "Q_s  = ", interval_str(res.Qs))
        println(io, "CERT_VRR_POSITIVE = ", fmt(res.cert_vrr))
        println(io, "CERT_QS_NEGATIVE  = ", fmt(res.cert_qs))
        println(io)
    end

    r128, r256 = results[128], results[256]
    consistent = (r128.cert_vrr == r256.cert_vrr) && (r128.cert_qs == r256.cert_qs)
    # a genuine tightness check: the 256-bit enclosure of V_rr, Q_s should be
    # no wider than the 128-bit one (both are outward-rounded enclosures of
    # the same true range over the same mathematical box, so higher working
    # precision must not produce a LARGER interval)
    tighter = (diam(r256.Vrr) <= diam(r128.Vrr) + 1e-6) && (diam(r256.Qs) <= diam(r128.Qs) + 1e-6)

    overall = r128.cert_vrr && r128.cert_qs && r256.cert_vrr && r256.cert_qs && consistent
    println(io, "128-bit vs 256-bit verdict consistent: ", fmt(consistent))
    println(io, "256-bit enclosure at least as tight as 128-bit: ", fmt(tighter))
    println(io, "CERT_PITCHFORK_LOCAL = ", overall ? "PASS" : "FAIL")

    out = String(take!(io))
    print(out)

    certs = [
        Certificate("pitchfork_local_128", "V_rr", string(inf(r128.Vrr)), string(sup(r128.Vrr)),
                    r128.cert_vrr ? "PASS" : "FAIL", "inf>0 required"),
        Certificate("pitchfork_local_128", "Q_s", string(inf(r128.Qs)), string(sup(r128.Qs)),
                    r128.cert_qs ? "PASS" : "FAIL", "sup<0 required"),
        Certificate("pitchfork_local_256", "V_rr", string(inf(r256.Vrr)), string(sup(r256.Vrr)),
                    r256.cert_vrr ? "PASS" : "FAIL", "inf>0 required"),
        Certificate("pitchfork_local_256", "Q_s", string(inf(r256.Qs)), string(sup(r256.Qs)),
                    r256.cert_qs ? "PASS" : "FAIL", "sup<0 required"),
    ]
    mkpath(joinpath(@__DIR__, "..", "output"))
    write_certificates_csv(joinpath(@__DIR__, "..", "output", "pitchfork_certificate.csv"), certs)
    open(joinpath(@__DIR__, "..", "output", "pitchfork_local_log.txt"), "w") do f
        write(f, out)
    end

    return overall
end

if abspath(PROGRAM_FILE) == @__FILE__
    run_pitchfork_local()
end
