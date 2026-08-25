"""
roots_lambda_08.jl (§9): completeness certificate for lambda=0.8. Expected: 4 roots, 3xS1 + 1xS2.
"""

isdefined(Main, :run_completeness) || include(joinpath(@__DIR__, "completeness_common.jl"))

function run_roots_lambda_08()
    run_completeness(4 // 5, "0.8"; expected_count = 4, expected_morse = Dict(:S1 => 3, :S2 => 1, :S0 => 0))
end

if abspath(PROGRAM_FILE) == @__FILE__
    run_roots_lambda_08()
end
