"""
asymmetric_cusp.jl  (§12)

Asymmetric family (M1,M2,M3) = (lambda, 1+eps, 1-eps) [lambda at P1, per the
U_asymmetric convention in mp_model.jl]. For each small eps, locate and
VALIDATE the fold point of the degeneracy system

    F(x,y,lambda) = (U_x, U_y, det Hess(U)) = 0

(3 equations, 3 unknowns), then compute

    R(eps) = (lambda_fold(eps) - lambda_-) / |eps|^(2/3)

and check it is consistent with the predicted cusp-unfolding constant
C ~ 1.58535311676 (a known reference value; NOT re-derived here — this script
tests numerical consistency with it, it does not prove the asymptotic
constant analytically).

TWO-PHASE METHOD (float pre-solve, then interval certify — standard
practice, and necessary here): the fold point moves by O(eps^(2/3)) in
lambda as eps shrinks from 0, which for eps=1e-5 is already a shift of
~4.6e-4 — far too large to guess a tight box directly, and a plain Newton
solve from the eps=0 point (r_-,0,lambda_-) DIVERGES for eps>=1e-3 (its
basin of attraction is too narrow — checked by hand). So:
  1. Float64 Newton, via a fine CONTINUATION LADDER in eps from 1e-6 up to
     the target eps (each step's solution seeds the next), locates an
     approximate fold point with residual ~1e-15 to 1e-16.
  2. A small interval box (relative half-width given below) is built
     around that Float64 point, converted to BigFloat, and `isolate_roots`
     (src/krawczyk.jl, the SAME machinery as §9/§11, here used for the
     first time on a genuinely 3-variable system) attempts to certify
     existence+uniqueness there.
"""

isdefined(Main, :ValidatedMP) || include(joinpath(@__DIR__, "..", "src", "ValidatedMP.jl"))
using .ValidatedMP
using IntervalArithmetic
using ForwardDiff
using LinearAlgebra
using Dates

function Fdeg_float(v, eps)
    x, y, lam = v
    g = ForwardDiff.gradient(u -> U_asymmetric(u[1], u[2], eps, lam), [x, y])
    H = ForwardDiff.hessian(u -> U_asymmetric(u[1], u[2], eps, lam), [x, y])
    detH = H[1, 1] * H[2, 2] - H[1, 2]^2
    return [g[1], g[2], detH]
end

function newton_solve(v0::Vector{Float64}, eps::Float64; iters::Int = 60)
    v = copy(v0)
    for _ in 1:iters
        Fv = Fdeg_float(v, eps)
        J = ForwardDiff.jacobian(u -> Fdeg_float(u, eps), v)
        dv = J \ Fv
        v = v .- dv
        norm(dv) < 1e-15 && break
    end
    return v
end

"Float64 pre-solve via a fine continuation ladder from eps=1e-6 up to the target eps."
function float_presolve(r_minus_f::Float64, lambda_minus_f::Float64, eps_target::Float64)
    v = [r_minus_f, 0.001, lambda_minus_f + 1.5 * abs(1e-6)^(2 / 3)]
    v = newton_solve(v, 1e-6)
    ladder = exp10.(range(-6, log10(eps_target), length = 30))
    for e in ladder
        v = newton_solve(v, e)
    end
    return v
end

"3-variable degeneracy system in interval/BigFloat arithmetic, generic element type."
function Fdeg_iv(v, eps)
    x, y, lam = v[1], v[2], v[3]
    f(u) = U_asymmetric(u[1], u[2], eps, lam)
    g = ForwardDiff.gradient(f, [x, y])
    H = ForwardDiff.hessian(f, [x, y])
    detH = H[1, 1] * H[2, 2] - H[1, 2]^2
    return [g[1], g[2], detH]
end

function certify_fold(eps_rational::Rational, r_minus_f::Float64, lambda_minus_f::Float64;
                       relhalf::Float64 = 1e-6)
    eps_f = Float64(eps_rational)
    vfloat = float_presolve(r_minus_f, lambda_minus_f, eps_f)
    resid = norm(Fdeg_float(vfloat, eps_f))

    setprecision(BigFloat, 256)
    T = BigFloat
    eps_iv = interval(T(eps_rational))
    # additive half-width scaled to each coordinate's own magnitude (plus a
    # tiny floor term) rather than a naive multiplicative `lo*(1-r),hi*(1+r)`
    # form, which inverts for negative coordinates (y can be negative here)
    box = [interval(T(vfloat[i]) - abs(T(vfloat[i])) * relhalf - T(relhalf) * T(1 // 1000),
                     T(vfloat[i]) + abs(T(vfloat[i])) * relhalf + T(relhalf) * T(1 // 1000)) for i in 1:3]

    F(v) = Fdeg_iv(v, eps_iv)
    results = isolate_roots(F, box; maxdepth = 60, mindiam = 1e-16)
    uniq = filter(r -> r.status == :unique, results)

    return (eps = eps_f, vfloat = vfloat, resid = resid, box0 = box,
            n_unique = length(uniq), n_total = length(results),
            n_unresolved = count(r -> r.status == :unresolved, results),
            root = isempty(uniq) ? nothing : uniq[1])
end

function run_asymmetric_cusp()
    io = IOBuffer()
    println(io, "=== Asymmetric family cusp scaling (sec 12) ===")
    println(io, "Julia version: ", VERSION)
    println(io, "Run date: ", Dates.now())
    println(io)

    setprecision(BigFloat, 256)
    T = BigFloat
    r_minus = (interval(T(-5)) + sqrt(interval(T(33)))) / interval(T(4))
    lambda_minus = lambda_onaxis(r_minus)
    r_minus_f = Float64(mid(r_minus))
    lambda_minus_f = Float64(mid(lambda_minus))

    C_ref = 1.58535311676
    eps_list = [1 // 100000, 1 // 10000, 1 // 1000, 1 // 100]

    println(io, "lambda_minus = ", interval_str(lambda_minus))
    println(io, "reference C  = ", C_ref)
    println(io)

    all_ok = true
    Rvals = Float64[]
    for eps_r in eps_list
        res = certify_fold(eps_r, r_minus_f, lambda_minus_f)
        println(io, "--- eps = $(res.eps) ---")
        println(io, "float presolve: x=$(res.vfloat[1]) y=$(res.vfloat[2]) lambda=$(res.vfloat[3])  residual=$(res.resid)")
        println(io, "Krawczyk search: total_leaves=$(res.n_total) unique=$(res.n_unique) unresolved=$(res.n_unresolved)")
        if res.root === nothing
            println(io, "CERTIFICATION FAILED for eps=$(res.eps): no unique root isolated in the search box")
            all_ok = false
            continue
        end
        lam_fold = res.root.box[3]
        dlam = lam_fold - lambda_minus
        R = dlam / interval(T(abs(res.eps)))^interval(T(2 // 3))
        push!(Rvals, Float64(mid(R)))
        println(io, "lambda_fold(eps) = ", interval_str(lam_fold))
        println(io, "R(eps) = (lambda_fold-lambda_-)/|eps|^(2/3) = ", interval_str(R))
        println(io, "  midpoint = ", Float64(mid(R)), "   vs C_ref = ", C_ref,
                "   |diff| = ", abs(Float64(mid(R)) - C_ref))
        println(io)
    end

    println(io, "R(eps) sequence, in INCREASING eps order [1e-5,1e-4,1e-3,1e-2] (|R-C_ref| should",
            " therefore be INCREASING along this list, i.e. R approaches C_ref as eps shrinks",
            " towards the start of the list): ", Rvals)
    monotone_approaching = length(Rvals) == length(eps_list) &&
                            all(abs(Rvals[i] - C_ref) < abs(Rvals[i + 1] - C_ref) for i in 1:(length(Rvals) - 1))
    println(io, "Strictly monotone approach to C_ref as eps shrinks? ", fmt(monotone_approaching))

    overall = all_ok && monotone_approaching
    println(io, "CERT_ASYMMETRIC_CUSP_SCALING = ", overall ? "PASS" : "FAIL")

    out = String(take!(io))
    print(out)
    mkpath(joinpath(@__DIR__, "..", "output"))
    open(joinpath(@__DIR__, "..", "output", "asymmetric_cusp_log.txt"), "w") do f
        write(f, out)
    end
    return overall
end

if abspath(PROGRAM_FILE) == @__FILE__
    run_asymmetric_cusp()
end
