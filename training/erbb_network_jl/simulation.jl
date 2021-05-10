module Sim
include("./name2idx/parameters.jl")
include("./name2idx/species.jl")
include("./set_model.jl")
include("./observable.jl")
include("./load_csv.jl")

using .C
using .V

using Sundials
using SteadyStateDiffEq

# Options for ODE solver
const ABSTOL = 1e-8
const RELTOL = 1e-8

normalization = Dict{String,Dict{}}()

const dt = 1.0
const t = collect(0.0:dt:120.0)  # 0, 1, 2, ..., 5400 [sec.]

const conditions = [
    "MCF7_EGF", "MCF7_HRG",
    "BT474_EGF", "BT474_HRG",
    "MDAMB231_EGF", "MDAMB231_HRG",
    "SKBR3_EGF", "SKBR3_HRG"
]

simulations = Array{Float64,3}(
    undef, length(observables), length(t), length(conditions)
)


function solveode(
        f::Function,
        u0::Vector{Float64},
        t::Vector{Float64},
        p::Vector{Float64})::Union{ODESolution{},Nothing}
    local sol::ODESolution{}, is_successful::Bool
    try
        prob = ODEProblem(f, u0, (t[1], t[end]), p)
        sol = solve(
            prob,CVODE_BDF(),
            abstol=ABSTOL,reltol=RELTOL,saveat=dt,dtmin=eps(),verbose=false
        )
        is_successful = ifelse(sol.retcode === :Success, true, false)
    catch
        is_successful = false
    finally
        if !is_successful
            GC.gc()
        end
    end
    return is_successful ? sol : nothing
end


function get_steady_state(
        f::Function,
        u0::Vector{Float64},
        p::Vector{Float64})::Vector{Float64}
    local sol::SteadyStateSolution{}, is_successful::Bool
    try
        prob = ODEProblem(diffeq, u0, (0.0, Inf), p)
        prob = SteadyStateProblem(prob)
        sol = solve(
            prob,
            DynamicSS(
                CVODE_BDF();abstol=ABSTOL,reltol=RELTOL
            ),
            dt=dt,verbose=false
        )
        is_successful = ifelse(sol.retcode === :Success, true, false)
    catch
        is_successful = false
    finally
        if !is_successful
            GC.gc()
        end
    end
    return is_successful ? sol.u : []
end


function simulate!(p::Vector{Float64}, u0::Vector{Float64})::Union{Bool,Nothing}
    p_in_mcf7::Vector{Float64} = copy(p)
    u0_in_mcf7::Vector{Float64} = copy(u0)
    # add ligand
    @inbounds for (i, condition) in enumerate(conditions)
        if condition == "MCF7_EGF"
            p .= p_in_mcf7
            u0 .= u0_in_mcf7
            # u0 = get_steady_state(diffeq,u0,p)
            if isempty(u0)
                return false
            else
                u0[V.EGF] = 10.0
                u0[V.HRG] = 0.0
            end
        elseif condition == "MCF7_HRG"
            p .= p_in_mcf7
            u0 .= u0_in_mcf7
            # u0 = get_steady_state(diffeq,u0,p)
            if isempty(u0)
                return false
            else
                u0[V.EGF] = 0.0
                u0[V.HRG] = 10.0
            end
        elseif condition == "BT474_EGF"
            p .= p_in_mcf7
            u0 .= u0_in_mcf7
            (p, u0) = scale_using_mcf7!(p, u0, "BT474_BREAST")
            # u0 = get_steady_state(diffeq,u0,p)
            if isempty(u0)
                return false
            else
                u0[V.EGF] = 10.0
                u0[V.HRG] = 0.0
            end
        elseif condition == "BT474_HRG"
            p .= p_in_mcf7
            u0 .= u0_in_mcf7
            (p, u0) = scale_using_mcf7!(p, u0, "BT474_BREAST")
            # u0 = get_steady_state(diffeq,u0,p)
            if isempty(u0)
                return false
            else
                u0[V.EGF] = 0.0
                u0[V.HRG] = 10.0
            end
        elseif condition == "MDAMB231_EGF"
            p .= p_in_mcf7
            u0 .= u0_in_mcf7
            (p, u0) = scale_using_mcf7!(p, u0, "MDAMB231_BREAST")
            # u0 = get_steady_state(diffeq,u0,p)
            if isempty(u0)
                return false
            else
                u0[V.EGF] = 10.0
                u0[V.HRG] = 0.0
            end
        elseif condition == "MDAMB231_HRG"
            p .= p_in_mcf7
            u0 .= u0_in_mcf7
            (p, u0) = scale_using_mcf7!(p, u0, "MDAMB231_BREAST")
            # u0 = get_steady_state(diffeq,u0,p)
            if isempty(u0)
                return false
            else
                u0[V.EGF] = 0.0
                u0[V.HRG] = 10.0
            end
        elseif condition == "SKBR3_EGF"
            p .= p_in_mcf7
            u0 .= u0_in_mcf7
            (p, u0) = scale_using_mcf7!(p, u0, "SKBR3_BREAST")
            # u0 = get_steady_state(diffeq,u0,p)
            if isempty(u0)
                return false
            else
                u0[V.EGF] = 10.0
                u0[V.HRG] = 0.0
            end
        elseif condition == "SKBR3_HRG"
            p .= p_in_mcf7
            u0 .= u0_in_mcf7
            (p, u0) = scale_using_mcf7!(p, u0, "SKBR3_BREAST")
            # u0 = get_steady_state(diffeq,u0,p)
            if isempty(u0)
                return false
            else
                u0[V.EGF] = 0.0
                u0[V.HRG] = 10.0
            end
        end
        sol = solveode(diffeq, u0, t, p)
        if sol === nothing
            return false
        else
            @inbounds @simd for j in eachindex(t)
                simulations[observables_index("Phosphorylated_Akt"),j,i] = (
                    sol.u[j][V.AktP]
                )
                simulations[observables_index("Phosphorylated_ERK"),j,i] = (
                    sol.u[j][V.ERKP] + sol.u[j][V.ERKPP]
                )
                simulations[observables_index("Phosphorylated_c-Myc"),j,i] = (
                    sol.u[j][V.cMycP] + sol.u[j][V.cMycPP]
                )
            end
        end
    end
end
end # module