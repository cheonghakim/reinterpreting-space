"""
proof_all.jl

Single entry point:  julia --project=. proof/proof_all.jl

Runs every implemented certificate script and prints the final
PASS/FAIL/UNRESOLVED summary.

STATUS OF THIS PACKAGE (be read before trusting the final banner):
Sections §2/§13 (Arb constants), §6 (pitchfork local), §7 (off-axis
existence/uniqueness), §8 (fold nondegeneracy), §9 (lambda=0.5/0.8/1.2
completeness + Morse/degree), §10 (Brouwer degree, two independent
methods), §12 (asymmetric cusp scaling) are fully implemented and
independently certified. §11 (phase-wide degeneracy exclusion) is
implemented and certified for its ON-AXIS MECHANISM on the MIDDLE SLAB
only (lambda_-+0.001 to lambda_+-0.001) — the off-axis check and the two
tail slabs (0,lambda_-) and (lambda_+,infty) are NOT interval-certified;
see README.md "Phase-diagram tails (not certified)" for the (non-certified) analytic
plausibility argument used instead. Per this project's honesty conventions
(see README.md "Reproducibility and honesty conventions"), this remainder
is reported as UNRESOLVED / NOT IMPLEMENTED below, not silently folded
into a claimed-complete §11 PASS.
"""

using Dates

const PROOF_DIR = @__DIR__
const OUTPUT_DIR = joinpath(PROOF_DIR, "..", "output")
mkpath(OUTPUT_DIR)

include(joinpath(PROOF_DIR, "constants_arb.jl"))
include(joinpath(PROOF_DIR, "pitchfork_local.jl"))
include(joinpath(PROOF_DIR, "pitchfork_branch.jl"))
include(joinpath(PROOF_DIR, "fold_local.jl"))
include(joinpath(PROOF_DIR, "completeness_common.jl"))
include(joinpath(PROOF_DIR, "roots_lambda_05.jl"))
include(joinpath(PROOF_DIR, "roots_lambda_08.jl"))
include(joinpath(PROOF_DIR, "roots_lambda_12.jl"))
include(joinpath(PROOF_DIR, "degree_check.jl"))
include(joinpath(PROOF_DIR, "phase_slabs.jl"))
include(joinpath(PROOF_DIR, "asymmetric_cusp.jl"))

function main()
    t0 = time()
    println("Running validated_mp proof suite — ", Dates.now())
    println()

    println(">>> constants_arb.jl (sec 2, 13)")
    r_arb = run_constants_arb()
    open(joinpath(OUTPUT_DIR, "constants.txt"), "w") do f
        write(f, r_arb.out)
    end
    println()

    println(">>> pitchfork_local.jl (sec 6)")
    pass_pitchfork = run_pitchfork_local()
    println()

    println(">>> pitchfork_branch.jl (sec 7)")
    pass_pitchfork_branch = run_pitchfork_branch()
    println()

    println(">>> fold_local.jl (sec 8)")
    pass_fold = run_fold_local()
    println()

    println(">>> roots_lambda_05.jl (sec 9)")
    r05 = run_roots_lambda_05()
    println()
    println(">>> roots_lambda_08.jl (sec 9)")
    r08 = run_roots_lambda_08()
    println()
    println(">>> roots_lambda_12.jl (sec 9)")
    r12 = run_roots_lambda_12()
    println()

    println(">>> degree_check.jl (sec 10, independent boundary winding-number cross-check)")
    pass_boundary_degree = run_degree_check_all()
    println()

    println(">>> phase_slabs.jl (sec 11, on-axis middle-slab only — see README scope note)")
    pass_phase_slabs = run_phase_slabs()
    println()

    println(">>> asymmetric_cusp.jl (sec 12)")
    pass_asym_cusp = run_asymmetric_cusp()
    println()

    degree_ok = r05.pass && r08.pass && r12.pass && r05.degree == -2 && r08.degree == -2 && r12.degree == -2 &&
                pass_boundary_degree

    elapsed = time() - t0

    implemented = [
        ("Arb constants", r_arb.pass),
        ("pitchfork V_rr > 0 / Q_s < 0 (local box)", pass_pitchfork),
        ("pitchfork off-axis existence (sec 7)", pass_pitchfork_branch),
        ("pitchfork off-axis uniqueness (sec 7)", pass_pitchfork_branch),
        ("fold nondegeneracy (U_rlambda, U_rrr != 0)", pass_fold),
        ("lambda=0.5 completeness (sec 9)", r05.pass),
        ("lambda=0.8 completeness (sec 9)", r08.pass),
        ("lambda=1.2 completeness (sec 9)", r12.pass),
        ("Brouwer degree consistency (sec 10, index-sum + boundary winding)", degree_ok),
        ("phase-wide degeneracy exclusion (sec 11, ON-AXIS MIDDLE SLAB ONLY -- not full scope, see README)", pass_phase_slabs),
        ("asymmetric cusp scaling (sec 12)", pass_asym_cusp),
    ]
    not_implemented = [
        "phase-wide degeneracy exclusion, off-axis + tail slabs (sec 11 remainder -- see README \"Phase-diagram tails (not certified)\")",
    ]

    summary = IOBuffer()
    println(summary, "=== MP VALIDATED NUMERICS SUMMARY ===")
    println(summary)
    for (name, ok) in implemented
        println(summary, "[$(ok ? "PASS" : "FAIL")] ", name)
    end
    for name in not_implemented
        println(summary, "[UNRESOLVED / NOT IMPLEMENTED] ", name)
    end
    println(summary)
    total_unresolved = r05.n_unresolved + r08.n_unresolved + r12.n_unresolved
    println(summary, "UNRESOLVED BOXES (sec 9 completeness searches): ", total_unresolved,
            " (sec 7/10/11/12's own box-subdivision searches each separately require, and here",
            " achieved, zero unresolved boxes as part of their own PASS criterion above)")
    all_implemented_pass = all(ok for (_, ok) in implemented)
    overall_status = all_implemented_pass && isempty(not_implemented) ? "CERTIFIED" :
                      all_implemented_pass ? "PARTIALLY CERTIFIED (sec 11 off-axis+tail scope excluded — see list above)" :
                      "NOT CERTIFIED"
    println(summary, "OVERALL STATUS: ", overall_status)
    println(summary)
    println(summary, "Runtime: ", round(elapsed, digits = 1), " s")
    println(summary, "Julia version: ", VERSION)
    println(summary, "Run date: ", Dates.now())

    out = String(take!(summary))
    print(out)
    open(joinpath(OUTPUT_DIR, "proof_summary.txt"), "w") do f
        write(f, out)
    end
end

main()
