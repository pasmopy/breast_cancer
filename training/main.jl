using BioMASS

const model = load_model("erbb_network_jl")

if abspath(PROGRAM_FILE) == @__FILE__
    optimize(
        model,
        parse(Int64, ARGS[1]),
        max_generation=50,
        allowable_error=6.0,
        popsize=3,
        local_search_method="DE",
        maxiter=200,
    )
end