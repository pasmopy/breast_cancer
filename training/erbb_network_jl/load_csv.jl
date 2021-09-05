using DataFrames, CSV


const postComBat = CSV.read("../transcriptomic_data/TPM_RLE_postComBat.csv", DataFrame);


function get_tpm(gene::String, cell_line::String)::Float64
    return postComBat[postComBat[!,Symbol("Description")] .== gene,Symbol(cell_line)][1]
end


function incorporate_ErbB1(p::Vector{Float64}, cell_line::String)::Float64
    return p[C.w_EGFR] * get_tpm("EGFR", cell_line)
end


function incorporate_ErbB2(p::Vector{Float64}, cell_line::String)::Float64
    return p[C.w_ERBB2] * get_tpm("ERBB2", cell_line)
end


function incorporate_ErbB3(p::Vector{Float64}, cell_line::String)::Float64
    return p[C.w_ERBB3] * get_tpm("ERBB3", cell_line)
end


function incorporate_Erb4(p::Vector{Float64}, cell_line::String)::Float64
    return p[C.w_ERBB4] * get_tpm("ERBB4", cell_line)
end


function incorporate_Grb2(p::Vector{Float64}, cell_line::String)::Float64
    return p[C.w_GRB2] * get_tpm("GRB2", cell_line)
end


function incorporate_Shc(p::Vector{Float64}, cell_line::String)::Float64
    return (p[C.w_SHC1] * get_tpm("SHC1", cell_line) +
            p[C.w_SHC2] * get_tpm("SHC2", cell_line) +
            p[C.w_SHC3] * get_tpm("SHC3", cell_line) +
            p[C.w_SHC4] * get_tpm("SHC4", cell_line))
end


function incorporate_RasGAP(p::Vector{Float64}, cell_line::String)::Float64
    return (p[C.w_RASA1] * get_tpm("RASA1", cell_line) +
            p[C.w_RASA2] * get_tpm("RASA2", cell_line) +
            p[C.w_RASA3] * get_tpm("RASA3", cell_line))
end


function incorporate_PI3K(p::Vector{Float64}, cell_line::String)::Float64
    return (p[C.w_PIK3CA] * get_tpm("PIK3CA", cell_line) +
            p[C.w_PIK3CB] * get_tpm("PIK3CB", cell_line) +
            p[C.w_PIK3CD] * get_tpm("PIK3CD", cell_line) +
            p[C.w_PIK3CG] * get_tpm("PIK3CG", cell_line))
end

function incorporate_PTEN(p::Vector{Float64}, cell_line::String)::Float64
    return p[C.w_PTEN] * get_tpm("PTEN", cell_line)
end


function incorporate_SOS(p::Vector{Float64}, cell_line::String)::Float64
    return (p[C.w_SOS1] * get_tpm("SOS1", cell_line) +
            p[C.w_SOS2] * get_tpm("SOS2", cell_line))
end


function incorporate_Gab1(p::Vector{Float64}, cell_line::String)::Float64
    return p[C.w_GAB1] * get_tpm("GAB1", cell_line)
end


function incorporate_RasGDP(p::Vector{Float64}, cell_line::String)::Float64
    return (p[C.w_HRAS] * get_tpm("HRAS", cell_line) +
            p[C.w_KRAS] * get_tpm("KRAS", cell_line) +
            p[C.w_NRAS] * get_tpm("NRAS", cell_line))
end


function incorporate_Raf(p::Vector{Float64}, cell_line::String)::Float64
    return (p[C.w_ARAF] * get_tpm("ARAF", cell_line) +
            p[C.w_BRAF] * get_tpm("BRAF", cell_line) +
            p[C.w_RAF1] * get_tpm("RAF1", cell_line))
end


function incorporate_MEK(p::Vector{Float64}, cell_line::String)::Float64
    return (p[C.w_MAP2K1] * get_tpm("MAP2K1", cell_line) +
            p[C.w_MAP2K2] * get_tpm("MAP2K2", cell_line))
end


function incorporate_ERK(p::Vector{Float64}, cell_line::String)::Float64
    return (p[C.w_MAPK1] * get_tpm("MAPK1", cell_line) +
            p[C.w_MAPK3] * get_tpm("MAPK3", cell_line))
end


function incorporate_Akt(p::Vector{Float64}, cell_line::String)::Float64
    return (p[C.w_AKT1] * get_tpm("AKT1", cell_line) +
            p[C.w_AKT2] * get_tpm("AKT2", cell_line))
end


function incorporate_PTP1B(p::Vector{Float64}, cell_line::String)::Float64
    return p[C.w_PTPN1] * get_tpm("PTPN1", cell_line)
end


function incorporate_GSK3b(p::Vector{Float64}, cell_line::String)::Float64
    return p[C.w_GSK3B] * get_tpm("GSK3B", cell_line)
end


function incorporate_DUSP(p::Vector{Float64}, cell_line::String)::Float64
    return (p[C.w_DUSP5] * get_tpm("DUSP5", cell_line) +
            p[C.w_DUSP6] * get_tpm("DUSP6", cell_line) +
            p[C.w_DUSP7] * get_tpm("DUSP7", cell_line))
end


function incorporate_cMyc(p::Vector{Float64}, cell_line::String)::Float64
    return p[C.w_MYC] * get_tpm("MYC", cell_line)
end


function individualize!(
        p::Vector{Float64},
        u0::Vector{Float64},
        cell_line::String)::Tuple{Vector{Float64},Vector{Float64}}
    # as initial concentration of model species
    u0[V.ErbB1] *= incorporate_ErbB1(p, cell_line)
    u0[V.ErbB3] *= incorporate_ErbB3(p, cell_line)
    u0[V.ErbB4] *= incorporate_Erb4(p, cell_line)
    u0[V.ErbB2] *= incorporate_ErbB2(p, cell_line)
    u0[V.Grb2] *= incorporate_Grb2(p, cell_line)
    u0[V.Shc] *= incorporate_Shc(p, cell_line)
    u0[V.RasGAP] *= incorporate_RasGAP(p, cell_line)
    u0[V.PI3K] *= incorporate_PI3K(p, cell_line)
    u0[V.PTP1B] *= incorporate_PTP1B(p, cell_line)
    u0[V.SOS] *= incorporate_SOS(p, cell_line)
    u0[V.Gab1] *= incorporate_Gab1(p, cell_line)
    u0[V.RasGDP] *= incorporate_RasGDP(p, cell_line)
    u0[V.Raf] *= incorporate_Raf(p, cell_line)
    u0[V.MEK] *= incorporate_MEK(p, cell_line)
    u0[V.ERK] *= incorporate_ERK(p, cell_line)
    u0[V.PTEN] *= incorporate_PTEN(p, cell_line)
    u0[V.Akt] *= incorporate_Akt(p, cell_line)
    u0[V.GSK3b] *= incorporate_GSK3b(p, cell_line)

    # as transcription rates
    p[C.V291] *= incorporate_DUSP(p, cell_line)
    p[C.V310] *= incorporate_cMyc(p, cell_line)

    return p, u0
end


function scale_using_mcf7!(
        p::Vector{Float64},
        u0::Vector{Float64},
        cell_line_2::String)::Tuple{Vector{Float64},Vector{Float64}}
    # as initial concentration of model species
    u0[V.ErbB1] *= incorporate_ErbB1(p, cell_line_2) / incorporate_ErbB1(p, "MCF7_BREAST")
    u0[V.ErbB3] *= incorporate_ErbB3(p, cell_line_2) / incorporate_ErbB3(p, "MCF7_BREAST")
    u0[V.ErbB4] *= incorporate_Erb4(p, cell_line_2) / incorporate_Erb4(p, "MCF7_BREAST")
    u0[V.ErbB2] *= incorporate_ErbB2(p, cell_line_2) / incorporate_ErbB2(p, "MCF7_BREAST")
    u0[V.Grb2] *= incorporate_Grb2(p, cell_line_2) / incorporate_Grb2(p, "MCF7_BREAST")
    u0[V.Shc] *= incorporate_Shc(p, cell_line_2) / incorporate_Shc(p, "MCF7_BREAST")
    u0[V.RasGAP] *= incorporate_RasGAP(p, cell_line_2) / incorporate_RasGAP(p, "MCF7_BREAST")
    u0[V.PI3K] *= incorporate_PI3K(p, cell_line_2) / incorporate_PI3K(p, "MCF7_BREAST")
    u0[V.PTP1B] *= incorporate_PTP1B(p, cell_line_2) / incorporate_PTP1B(p, "MCF7_BREAST")
    u0[V.SOS] *= incorporate_SOS(p, cell_line_2) / incorporate_SOS(p, "MCF7_BREAST")
    u0[V.Gab1] *= incorporate_Gab1(p, cell_line_2) / incorporate_Gab1(p, "MCF7_BREAST")
    u0[V.RasGDP] *= incorporate_RasGDP(p, cell_line_2) / incorporate_RasGDP(p, "MCF7_BREAST")
    u0[V.Raf] *= incorporate_Raf(p, cell_line_2) / incorporate_Raf(p, "MCF7_BREAST")
    u0[V.MEK] *= incorporate_MEK(p, cell_line_2) / incorporate_MEK(p, "MCF7_BREAST")
    u0[V.ERK] *= incorporate_ERK(p, cell_line_2) / incorporate_ERK(p, "MCF7_BREAST")
    u0[V.PTEN] *= incorporate_PTEN(p, cell_line_2) / incorporate_PTEN(p, "MCF7_BREAST")
    u0[V.Akt] *= incorporate_Akt(p, cell_line_2) / incorporate_Akt(p, "MCF7_BREAST")
    u0[V.GSK3b] *= incorporate_GSK3b(p, cell_line_2) / incorporate_GSK3b(p, "MCF7_BREAST")

    # as transcription rates
    p[C.V291] *= incorporate_DUSP(p, cell_line_2) / incorporate_DUSP(p, "MCF7_BREAST")
    p[C.V310] *= incorporate_cMyc(p, cell_line_2) / incorporate_cMyc(p, "MCF7_BREAST")

    return p, u0
end