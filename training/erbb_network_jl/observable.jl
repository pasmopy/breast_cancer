const observables = [
    "Phosphorylated_Akt",
    "Phosphorylated_ERK",
    "Phosphorylated_c-Myc",
]

function observables_index(observable_name::String)::Union{Int,Nothing}
    if !(observable_name in observables)
        error("$observable_name is not defined in observables.")
    end
    return findfirst(isequal(observable_name), observables)
end