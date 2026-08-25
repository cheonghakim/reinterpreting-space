"""
roots_lambda_12.jl (§9): completeness certificate for lambda=1.2. Expected: 2 roots, both S1.
"""

isdefined(Main, :run_completeness) || include(joinpath(@__DIR__, "completeness_common.jl"))

function run_roots_lambda_12()
    run_completeness(6 // 5, "1.2"; expected_count = 2, expected_morse = Dict(:S1 => 2, :S2 => 0, :S0 => 0))
end

if abspath(PROGRAM_FILE) == @__FILE__
    run_roots_lambda_12()
end
