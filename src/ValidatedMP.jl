"""
ValidatedMP

Umbrella module for the validated_mp certificate package. Pulls together
mp_model.jl (exact algebraic model), derivatives.jl (rigorous AD-based
derivatives over interval/Arb types), and interval_utils.jl (sign-check and
certificate I/O helpers) into one namespace, so proof/*.jl scripts do

    include(joinpath(@__DIR__, "..", "src", "ValidatedMP.jl"))
    using .ValidatedMP

and get everything without fragile relative module paths.
"""
module ValidatedMP

using IntervalArithmetic
using ForwardDiff
using Logging

# IntervalArithmetic's generic BigFloat matrix multiplication (used throughout
# the Krawczyk machinery below) logs an @info "Fast multiplication is only
# supported for Float16/32/64..." on EVERY matrix multiply. This is expected
# (BigFloat/Interval always takes the generic path) and floods the console
# during any box-subdivision search; suppress Info-level messages for this
# module's computations while leaving Warn/Error visible.
global_logger(ConsoleLogger(stderr, Logging.Warn))

export cartesian_sources, U_cartesian, U_symmetric, U_asymmetric, V_rs, R1sq, R2sq, R3sq, lambda_onaxis, onaxis_r_derivs,
       V_second_order, U_grad_hess, onaxis_r_derivs_ad, V_full_second_order,
       signcheck_pos, signcheck_neg, signcheck_nonzero, overlaps_iv,
       Certificate, write_certificates_csv, interval_str, fmt,
       krawczyk_operator, isolate_roots, KrawczykResult

include("mp_model.jl")
include("derivatives.jl")
include("interval_utils.jl")
include("krawczyk.jl")

end # module
