"""
completeness_common.jl  (§9 shared driver, called by roots_lambda_{05,08,12}.jl)

For a given lambda, certifies:
  - completeness: box-subdivision of a padded bounding box of the source
    triangle, with exclusion disks (radius 0.25, see README.md
    "Search-domain justification: convex hull and exclusion disks") around each singular source, Krawczyk
    root isolation (src/krawczyk.jl) on the rest => every leaf box is
    :unique, :none, or :excluded. UNRESOLVED must be 0.
  - each :unique root's Hessian determinant sign (2D, in-plane block of the
    full 3D U — see README.md "Morse index convention", derived from
    U_zz<0 always at z=0, NOT >0 — an easy sign slip this project's own
    first attempt made and had to correct against the lambda=0.8 data):
    det<0 => S1, det>0 & trace>0 => S2, det>0 & trace<0 => S0 (not expected
    for any lambda used here; flagged loudly if it occurs, never silently
    mislabeled).
  - root count and Morse multiset match the expected values.

Search domain: x in [-0.55,1.05], y in [-0.91,0.91] — a small pad around
the bounding box of P1=(1,0), P2=(-1/2,sqrt(3)/2), P3=(-1/2,-sqrt(3)/2).
Per README's convex-hull lemma, every critical point is inside the
triangle's hull, itself inside this box, so nothing outside the physically
relevant region is missed; extra area outside the hull but inside this box
is harmless (must and does resolve to :none).
"""

isdefined(Main, :ValidatedMP) || include(joinpath(@__DIR__, "..", "src", "ValidatedMP.jl"))
using .ValidatedMP
using IntervalArithmetic
using ForwardDiff
using Dates

const RHO_EXCL2 = 1 // 16  # (0.25)^2, exact rational

function run_completeness(lambda_val::Rational, label::String;
                           expected_count::Int, expected_morse::Dict{Symbol, Int},
                           maxdepth::Int = 70, mindiam::Float64 = 1e-11)
    io = IOBuffer()
    println(io, "=== lambda=$(label) completeness (sec 9) ===")
    println(io, "Julia version: ", VERSION)
    println(io, "Run date: ", Dates.now())
    println(io)

    setprecision(BigFloat, 128)
    T = BigFloat
    lambda_iv = interval(T(lambda_val))
    P1, P2, P3 = cartesian_sources(T)
    sources = (P1, P2, P3)
    rho2 = interval(T(RHO_EXCL2))

    F(v) = ForwardDiff.gradient(u -> U_symmetric(u[1], u[2], lambda_iv), v)

    domain = [interval(T(-11 // 20), T(21 // 20)), interval(T(-91 // 100), T(91 // 100))]

    function exclude(box)
        x, y = box
        for P in sources
            d2 = (x - interval(P[1]))^2 + (y - interval(P[2]))^2
            if sup(d2) < inf(rho2)
                return true
            end
        end
        return false
    end

    t0 = time()
    results = isolate_roots(F, domain; maxdepth = maxdepth, mindiam = mindiam, exclude = exclude)
    elapsed = time() - t0

    n_unique = count(r -> r.status == :unique, results)
    n_none = count(r -> r.status == :none, results)
    n_excluded = count(r -> r.status == :excluded, results)
    n_unresolved = count(r -> r.status == :unresolved, results)

    println(io, "Total leaf boxes: ", length(results), "  (runtime ", round(elapsed, digits = 1), " s)")
    println(io, "UNIQUE=", n_unique, "  NONE=", n_none, "  EXCLUDED=", n_excluded, "  UNRESOLVED=", n_unresolved)
    println(io)

    morse_counts = Dict{Symbol, Int}(:S0 => 0, :S1 => 0, :S2 => 0)
    root_rows = Any[]
    all_classified = true
    for r in results
        r.status == :unique || continue
        x, y = r.box
        d = U_grad_hess(x, y, lambda_iv, interval(T(1)), interval(T(1)))
        detH = d.Uxx * d.Uyy - d.Uxy^2
        trace = d.Uxx + d.Uyy
        if signcheck_neg(detH)
            morse_counts[:S1] += 1
            label_s = :S1
        elseif signcheck_pos(detH) && signcheck_pos(trace)
            morse_counts[:S2] += 1
            label_s = :S2
        elseif signcheck_pos(detH) && signcheck_neg(trace)
            morse_counts[:S0] += 1
            label_s = :S0
        else
            all_classified = false
            label_s = :UNRESOLVED_MORSE
        end
        push!(root_rows, (x = x, y = y, detH = detH, trace = trace, morse = label_s))
        println(io, "root at x=", interval_str(x), " y=", interval_str(y),
                "  detH=", interval_str(detH), " trace=", interval_str(trace), " => ", label_s)
    end
    println(io)
    println(io, "Morse counts: S0=", morse_counts[:S0], " S1=", morse_counts[:S1], " S2=", morse_counts[:S2])
    println(io, "Expected: ", expected_morse)

    count_ok = n_unique == expected_count
    unresolved_ok = n_unresolved == 0
    morse_ok = all_classified && morse_counts[:S1] == get(expected_morse, :S1, 0) &&
               morse_counts[:S2] == get(expected_morse, :S2, 0) &&
               morse_counts[:S0] == get(expected_morse, :S0, 0)
    degree = sum(signcheck_neg(rr.detH) ? -1 : 1 for rr in root_rows; init = 0)

    println(io, "Brouwer degree (sum of sgn det H_U over roots) = ", degree, "  (expected -2)")
    println(io)
    println(io, "CERT_COUNT_$(label) = ", count_ok ? "PASS" : "FAIL", " ($(n_unique) vs expected $(expected_count))")
    println(io, "CERT_UNRESOLVED_$(label) = ", unresolved_ok ? "PASS" : "FAIL")
    println(io, "CERT_MORSE_$(label) = ", morse_ok ? "PASS" : "FAIL")
    println(io, "CERT_DEGREE_$(label) = ", degree == -2 ? "PASS" : "FAIL")

    overall = count_ok && unresolved_ok && morse_ok && degree == -2
    println(io, "CERT_COMPLETENESS_$(label)_OVERALL = ", overall ? "PASS" : "FAIL")

    out = String(take!(io))
    print(out)
    mkpath(joinpath(@__DIR__, "..", "output"))
    open(joinpath(@__DIR__, "..", "output", "completeness_$(label)_log.txt"), "w") do f
        write(f, out)
    end

    certs = Certificate[]
    for (i, rr) in enumerate(root_rows)
        push!(certs, Certificate("lambda_$(label)_root_$(i)", "x", string(inf(rr.x)), string(sup(rr.x)), "PASS", string(rr.morse)))
        push!(certs, Certificate("lambda_$(label)_root_$(i)", "y", string(inf(rr.y)), string(sup(rr.y)), "PASS", string(rr.morse)))
    end
    write_certificates_csv(joinpath(@__DIR__, "..", "output", "root_certificates_$(label).csv"), certs)

    return (pass = overall, n_unique = n_unique, n_unresolved = n_unresolved, degree = degree,
            morse_counts = morse_counts, root_rows = root_rows)
end
