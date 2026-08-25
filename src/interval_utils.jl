"""
interval_utils.jl

Small utilities shared by all certificate scripts:
- exact-rational-safe interval construction (`@interval`, checked),
- a `Certificate` record type with strict PASS/FAIL/UNRESOLVED semantics,
- CSV/text writers that always dump the full enclosure, never just a verdict.

STRICT RULE ENFORCED HERE: `signcheck_pos`/`signcheck_neg` below are the ONLY
sanctioned way to turn an interval into a sign verdict in this codebase. They
require the interval to be a *proper, finite* interval (not empty, not
containing NaN) and return true only when the WHOLE interval has that sign
— never silently guess a sign when 0 is in the interval.

`include`d into ValidatedMP (see mp_model.jl docstring for the pattern).
"""

"Strict: returns true iff the WHOLE interval is guaranteed > 0."
signcheck_pos(x) = !isempty_interval(x) && inf(x) > 0

"Strict: returns true iff the WHOLE interval is guaranteed < 0."
signcheck_neg(x) = !isempty_interval(x) && sup(x) < 0

"Strict: returns true iff the WHOLE interval is guaranteed to exclude 0."
signcheck_nonzero(x) = signcheck_pos(x) || signcheck_neg(x)

"""
overlaps_iv(a,b)

True iff the two intervals share at least one point. Used to cross-check two
INDEPENDENTLY derived enclosures of the same true quantity (e.g. an
AD-derived enclosure vs. a hand-derived analytic-formula enclosure): if they
provably do not overlap, at least one derivation is wrong and this must be
investigated, not papered over.
"""
overlaps_iv(a, b) = inf(a) <= sup(b) && inf(b) <= sup(a)

"One row of certificate output: never store just a verdict, always store the enclosure."
struct Certificate
    name::String
    quantity::String         # e.g. "V_rr", "Q_s"
    lo::String
    hi::String
    verdict::String          # "PASS", "FAIL", "UNRESOLVED"
    note::String
end

function interval_str(x)
    return "[$(inf(x)), $(sup(x))]"
end

fmt(x::Bool) = x ? "true" : "false"

function write_certificates_csv(path::String, certs::Vector{Certificate})
    open(path, "w") do io
        println(io, "name,quantity,lo,hi,verdict,note")
        for c in certs
            println(io, "$(c.name),$(c.quantity),$(c.lo),$(c.hi),$(c.verdict),\"$(c.note)\"")
        end
    end
end
