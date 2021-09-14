include("./load_csv.jl")

# Specify model parameters and/or initial values to optimize
function get_search_index()::Tuple{Array{Int64,1},Array{Int64,1}}
    # parameters
    search_idx_params::Vector{Int} = [
        C.kf1,
        C.kr1,
        C.kf2,
        C.kr2,
        C.kf3,
        C.kr3,
        C.kf4,
        C.kr4,
        C.kf5,
        C.kr5,
        C.kf6,
        C.kr6,
        C.kf7,
        C.kr7,
        C.kf8,
        C.kr8,
        C.kf9,
        C.kr9,
        C.kf10,
        C.kr10,
        C.kf11,
        C.kr11,
        C.kf12,
        C.kr12,
        C.kf13,
        C.kr13,
        C.kf14,
        C.kr14,
        C.kf15,
        C.kr15,
        C.kf16,
        C.kr16,
        C.kf17,
        C.kr17,
        C.kf18,
        C.kr18,
        C.kf19,
        C.kr19,
        C.kf20,
        C.kr20,
        C.kf21,
        C.kr21,
        C.kf22,
        C.kr22,
        C.kf23,
        C.kr23,
        C.kf24,
        C.kr24,
        C.kf25,
        C.kr25,
        C.kf26,
        C.kr26,
        C.kf27,
        C.kr27,
        C.kf28,
        C.kr28,
        C.kf29,
        C.kr29,
        C.kf30,
        C.kr30,
        C.kf31,
        C.kr31,
        C.kf32,
        C.kr32,
        C.kf33,
        C.kr33,
        C.kf34,
        C.kr34,
        C.kf35,
        C.kr35,
        C.kf36,
        C.kr36,
        C.kf37,
        C.kr37,
        C.kf38,
        C.kr38,
        # C.kf39,
        # C.kr39,
        # C.kf40,
        # C.kr40,
        # C.kf41,
        # C.kr41,
        # C.kf42,
        # C.kr42,
        # C.kf43,
        # C.kr43,
        C.kf44,
        # C.kr44,
        # C.kf45,
        # C.kr45,
        # C.kf46,
        # C.kr46,
        # C.kf47,
        # C.kr47,
        # C.kf48,
        # C.kr48,
        # C.kf49,
        # C.kr49,
        C.kf50,
        C.kr50,
        # C.kf51,
        # C.kr51,
        # C.kf52,
        # C.kr52,
        # C.kf53,
        # C.kr53,
        # C.kf54,
        # C.kr54,
        # C.kf55,
        # C.kr55,
        C.kf56,
        # C.kr56,
        # C.kf57,
        # C.kr57,
        # C.kf58,
        # C.kr58,
        # C.kf59,
        # C.kr59,
        # C.kf60,
        # C.kr60,
        # C.kf61,
        # C.kr61,
        C.kf62,
        C.kr62,
        # C.kf63,
        # C.kr63,
        # C.kf64,
        # C.kr64,
        # C.kf65,
        # C.kr65,
        # C.kf66,
        # C.kr66,
        # C.kf67,
        # C.kr67,
        # C.kf68,
        # C.kr68,
        # C.kf69,
        # C.kr69,
        # C.kf70,
        # C.kr70,
        # C.kf71,
        # C.kr71,
        # C.kf72,
        # C.kr72,
        # C.kf73,
        # C.kr73,
        C.kf74,
        C.kr74,
        # C.kf75,
        # C.kr75,
        # C.kf76,
        # C.kr76,
        # C.kf77,
        # C.kr77,
        # C.kf78,
        # C.kr78,
        # C.kf79,
        # C.kr79,
        # C.kf80,
        # C.kr80,
        # C.kf81,
        # C.kr81,
        # C.kf82,
        # C.kr82,
        # C.kf83,
        # C.kr83,
        # C.kf84,
        # C.kr84,
        # C.kf85,
        # C.kr85,
        C.kf86,
        C.kr86,
        # C.kf87,
        # C.kr87,
        # C.kf88,
        # C.kr88,
        # C.kf89,
        # C.kr89,
        # C.kf90,
        # C.kr90,
        # C.kf91,
        # C.kr91,
        # C.kf92,
        # C.kr92,
        # C.kf93,
        # C.kr93,
        # C.kf94,
        # C.kr94,
        # C.kf95,
        # C.kr95,
        # C.kf96,
        # C.kr96,
        # C.kf97,
        # C.kr97,
        C.kf98,
        # C.kr98,
        # C.kf99,
        # C.kr99,
        # C.kf100,
        # C.kr100,
        # C.kf101,
        # C.kr101,
        # C.kf102,
        # C.kr102,
        # C.kf103,
        # C.kr103,
        # C.kf104,
        # C.kr104,
        # C.kf105,
        # C.kr105,
        # C.kf106,
        # C.kr106,
        # C.kf107,
        # C.kr107,
        # C.kf108,
        # C.kr108,
        # C.kf109,
        # C.kr109,
        C.V110,
        C.K110,
        C.V111,
        C.K111,
        C.V112,
        C.K112,
        C.V113,
        C.K113,
        C.V114,
        C.K114,
        C.V115,
        C.K115,
        C.V116,
        C.K116,
        C.V117,
        C.K117,
        C.V118,
        C.K118,
        C.V119,
        C.K119,
        C.V120,
        C.K120,
        C.V121,
        C.K121,
        C.kf122,
        # C.kr122,
        # C.kf123,
        # C.kr123,
        # C.kf124,
        # C.kr124,
        # C.kf125,
        # C.kr125,
        # C.kf126,
        # C.kr126,
        # C.kf127,
        # C.kr127,
        # C.kf128,
        # C.kr128,
        # C.kf129,
        # C.kr129,
        # C.kf130,
        # C.kr130,
        # C.kf131,
        # C.kr131,
        # C.kf132,
        # C.kr132,
        # C.kf133,
        # C.kr133,
        # C.kf134,
        # C.kr134,
        C.kf135,
        C.kr135,
        # C.kf136,
        # C.kr136,
        # C.kf137,
        # C.kr137,
        # C.kf138,
        # C.kr138,
        # C.kf139,
        # C.kr139,
        # C.kf140,
        # C.kr140,
        # C.kf141,
        # C.kr141,
        # C.kf142,
        # C.kr142,
        # C.kf143,
        # C.kr143,
        # C.kf144,
        # C.kr144,
        # C.kf145,
        # C.kr145,
        # C.kf146,
        # C.kr146,
        # C.kf147,
        # C.kr147,
        C.kf148,
        # C.kr148,
        # C.kf149,
        # C.kr149,
        # C.kf150,
        # C.kr150,
        # C.kf151,
        # C.kr151,
        # C.kf152,
        # C.kr152,
        # C.kf153,
        # C.kr153,
        # C.kf154,
        # C.kr154,
        # C.kf155,
        # C.kr155,
        # C.kf156,
        # C.kr156,
        # C.kf157,
        # C.kr157,
        # C.kf158,
        # C.kr158,
        # C.kf159,
        # C.kr159,
        # C.kf160,
        # C.kr160,
        C.kf161,
        C.kr161,
        # C.kf162,
        # C.kr162,
        # C.kf163,
        # C.kr163,
        # C.kf164,
        # C.kr164,
        # C.kf165,
        # C.kr165,
        # C.kf166,
        # C.kr166,
        # C.kf167,
        # C.kr167,
        # C.kf168,
        # C.kr168,
        # C.kf169,
        # C.kr169,
        # C.kf170,
        # C.kr170,
        # C.kf171,
        # C.kr171,
        # C.kf172,
        # C.kr172,
        # C.kf173,
        # C.kr173,
        C.kf174,
        C.kr174,
        # C.kf175,
        # C.kr175,
        # C.kf176,
        # C.kr176,
        # C.kf177,
        # C.kr177,
        # C.kf178,
        # C.kr178,
        # C.kf179,
        # C.kr179,
        # C.kf180,
        # C.kr180,
        # C.kf181,
        # C.kr181,
        # C.kf182,
        # C.kr182,
        # C.kf183,
        # C.kr183,
        # C.kf184,
        # C.kr184,
        # C.kf185,
        # C.kr185,
        # C.kf186,
        # C.kr186,
        # C.kf187,
        # C.kr187,
        # C.kf188,
        # C.kr188,
        # C.kf189,
        # C.kr189,
        # C.kf190,
        # C.kr190,
        # C.kf191,
        # C.kr191,
        # C.kf192,
        # C.kr192,
        C.kf193,
        # C.kr193,
        # C.kf194,
        # C.kr194,
        # C.kf195,
        # C.kr195,
        # C.kf196,
        # C.kr196,
        # C.kf197,
        # C.kr197,
        # C.kf198,
        # C.kr198,
        # C.kf199,
        # C.kr199,
        # C.kf200,
        # C.kr200,
        # C.kf201,
        # C.kr201,
        # C.kf202,
        # C.kr202,
        # C.kf203,
        # C.kr203,
        # C.kf204,
        # C.kr204,
        # C.kf205,
        # C.kr205,
        # C.kf206,
        # C.kr206,
        # C.kf207,
        # C.kr207,
        # C.kf208,
        # C.kr208,
        # C.kf209,
        # C.kr209,
        # C.kf210,
        # C.kr210,
        # C.kf211,
        # C.kr211,
        C.kf212,
        C.kr212,
        # C.kf213,
        # C.kr213,
        # C.kf214,
        # C.kr214,
        # C.kf215,
        # C.kr215,
        # C.kf216,
        # C.kr216,
        # C.kf217,
        # C.kr217,
        # C.kf218,
        # C.kr218,
        # C.kf219,
        # C.kr219,
        # C.kf220,
        # C.kr220,
        # C.kf221,
        # C.kr221,
        # C.kf222,
        # C.kr222,
        # C.kf223,
        # C.kr223,
        # C.kf224,
        # C.kr224,
        C.kf225,
        C.kr225,
        # C.kf226,
        # C.kr226,
        # C.kf227,
        # C.kr227,
        # C.kf228,
        # C.kr228,
        # C.kf229,
        # C.kr229,
        # C.kf230,
        # C.kr230,
        # C.kf231,
        # C.kr231,
        # C.kf232,
        # C.kr232,
        # C.kf233,
        # C.kr233,
        # C.kf234,
        # C.kr234,
        # C.kf235,
        # C.kr235,
        # C.kf236,
        # C.kr236,
        # C.kf237,
        # C.kr237,
        # C.kf238,
        # C.kr238,
        # C.kf239,
        # C.kr239,
        # C.kf240,
        # C.kr240,
        # C.kf241,
        # C.kr241,
        C.kf242,
        # C.kr242,
        # C.kf243,
        # C.kr243,
        # C.kf244,
        # C.kr244,
        # C.kf245,
        # C.kr245,
        # C.kf246,
        # C.kr246,
        # C.kf247,
        # C.kr247,
        # C.kf248,
        # C.kr248,
        # C.kf249,
        # C.kr249,
        # C.kf250,
        # C.kr250,
        # C.kf251,
        # C.kr251,
        # C.kf252,
        # C.kr252,
        # C.kf253,
        # C.kr253,
        # C.kf254,
        # C.kr254,
        # C.kf255,
        # C.kr255,
        # C.kf256,
        # C.kr256,
        # C.kf257,
        # C.kr257,
        # C.kf258,
        # C.kr258,
        C.kf259,
        C.kr259,
        C.V260,
        C.K260,
        C.V261,
        C.K261,
        C.V262,
        C.K262,
        C.V263,
        C.K263,
        C.V264,
        C.K264,
        C.V265,
        C.K265,
        C.V266,
        C.K266,
        C.kf267,
        C.kf268,
        # C.kf269,
        # C.kf270,
        # C.kf271,
        # C.kf272,
        # C.kf273,
        # C.kf274,
        # C.kf275,
        # C.kf276,
        # C.kf277,
        # C.kf278,
        # C.kf279,
        # C.kf280,
        # C.kf281,
        # C.kf282,
        # C.kf283,
        # C.kf284,
        # C.kf285,
        # C.kf286,
        # C.kf287,
        # C.kf288,
        # C.kf289,
        # C.kf290,
        C.V291,
        C.K291,
        # C.n291,
        C.kf292,
        # C.kr292,
        C.kf293,
        C.kf294,
        C.V295,
        C.K295,
        C.kf296,
        C.kf297,
        C.kf298,
        C.kr298,
        C.kf299,
        C.kr299,
        C.kf300,
        C.kr300,
        C.kf301,
        C.kr301,
        C.kf302,
        C.kr302,
        # C.kf303,
        # C.kr303,
        # C.kf304,
        # C.kr304,
        # C.kf305,
        # C.kr305,
        # C.kf306,
        # C.kr306,
        # C.kf307,
        # C.kr307,
        C.V308,
        C.K308,
        C.V309,
        C.K309,
        C.V310,
        C.K310,
        # C.n310,
        C.KF310,
        # C.nF310,
        C.kf311,
        # C.kr311,
        C.kf312,
        C.kf313,
        C.V314,
        C.K314,
        C.V315,
        C.K315,
        C.V316,
        C.K316,
        C.V317,
        C.K317,
        C.kf318,
        C.kf319,
        #
        C.w_EGFR,
        C.w_ERBB2,
        C.w_ERBB3,
        C.w_ERBB4,
        #
        C.w_GRB2,
        #
        C.w_SHC1,
        C.w_SHC2,
        C.w_SHC3,
        C.w_SHC4,
        #
        C.w_RASA1,
        C.w_RASA2,
        C.w_RASA3,
        #
        C.w_PIK3CA,
        C.w_PIK3CB,
        C.w_PIK3CD,
        C.w_PIK3CG,
        #
        C.w_PTEN,
        #
        C.w_SOS1,
        C.w_SOS2,
        #
        C.w_GAB1,
        #
        C.w_HRAS,
        C.w_KRAS,
        C.w_NRAS,
        #
        C.w_ARAF,
        C.w_BRAF,
        C.w_RAF1,
        #
        C.w_MAP2K1,
        C.w_MAP2K2,
        #
        C.w_MAPK1,
        C.w_MAPK3,
        #
        C.w_AKT1,
        C.w_AKT2,
        #
        C.w_PTPN1,
        #
        C.w_GSK3B,
        #
        C.w_DUSP5,
        C.w_DUSP6,
        C.w_DUSP7,
        #
        C.w_MYC,
    ]

    # initial values
    search_idx_initials::Vector{Int} = [
        V.PIP2
    ]

    return search_idx_params, search_idx_initials
end


function get_search_region()::Matrix{Float64}
    p::Vector{Float64} = param_values()
    u0::Vector{Float64} = initial_values()

    search_idx::Tuple{Array{Int64,1},Array{Int64,1}} = get_search_index()
    search_param::Vector{Float64} = init_search_param(search_idx, p, u0)

    search_rgn::Matrix{Float64} = zeros(2, length(p) + length(u0))

    # Default: 0.1 ~ 10x
    for (i, j) in enumerate(search_idx[1])
        search_rgn[1,j] = 1e-4  # lower bound
        search_rgn[2,j] = 1e+5  # upper bound
    end

    # Default: 0.5 ~ 2x
    for (i, j) in enumerate(search_idx[2])
        search_rgn[1,j + length(p)] = search_param[i + length(search_idx[1])] * 1  # lower bound
        search_rgn[2,j + length(p)] = search_param[i + length(search_idx[1])] * 1000  # upper bound
    end

    # search_rgn[:, C.param_name] = [lower_bound, upper_bound]
    # search_rgn[:, V.var_name+length(p)] = [lower_bound, upper_bound]

    for input_layer in [
        C.kf1,
        C.kr1,
        C.kf2,
        C.kr2,
        C.kf3,
        C.kr3,
        C.kf4,
        C.kr4,
        C.kf5,
        C.kr5,
        C.kf6,
        C.kr6,
        C.kf7,
        C.kr7,
        C.kf8,
        C.kr8,
        C.kf9,
        C.kr9,
    ]
        search_rgn[:, input_layer] = [1e-2, 1e+2]
    end

    for receptor_phosphorylation in [
        C.kf10,
        C.kr10,
        C.kf11,
        C.kr11,
        C.kf12,
        C.kr12,
        C.kf13,
        C.kr13,
        C.kf14,
        C.kr14,
        C.kf15,
        C.kr15,
    ]
        search_rgn[:, receptor_phosphorylation] = [1e-2, 1e+2]
    end

    for recepotor_adaptor in [
        C.kf16,
        C.kr16,
        C.kf17,
        C.kr17,
        C.kf18,
        C.kr18,
        C.kf19,
        C.kr19,
        C.kf20,
        C.kr20,
        C.kf21,
        C.kr21,
        C.kf22,
        C.kr22,
        C.kf23,
        C.kr23,
        C.kf24,
        C.kr24,
        C.kf25,
        C.kr25,
        C.kf26,
        C.kr26,
        C.kf27,
        C.kr27,
        C.kf28,
        C.kr28,
        C.kf29,
        C.kr29,
        C.kf30,
        C.kr30,
        C.kf31,
        C.kr31,
        C.kf32,
        C.kr32,
        C.kf33,
        C.kr33,
        C.kf34,
        C.kr34,
        C.kf35,
        C.kr35,
        C.kf36,
        C.kr36,
        C.kf37,
        C.kr37,
    ]
        search_rgn[:, recepotor_adaptor] = [1e-4, 1e+2]
    end

    for degradation in [
        C.kf267,
        C.kf268,
        C.kf294,
        C.kf296,
        C.kf297,
        C.kf312,
        C.kf318,
        C.kf319,
    ]
        search_rgn[:, degradation] = [1e-4, 1e+0]
    end

    for EC50 in [
        C.K110,
        C.K111,
        C.K112,
        C.K113,
        C.K114,
        C.K115,
        C.K116,
        C.K117,
        C.K118,
        C.K119,
        C.K120,
        C.K121,
        C.K260,
        C.K261,
        C.K262,
        C.K263,
        C.K264,
        C.K265,
        C.K266,
        C.K291,
        C.K295,
        C.K308,
        C.K309,
        C.K310,
        C.K314,
        C.K315,
        C.K316,
        C.K317,
    ]
        search_rgn[:, EC50] = [1e+0, 1e+4]
    end

    search_rgn[:, C.w_EGFR] = [0.01, 100]
    search_rgn[:, C.w_ERBB2] = [0.01, 100]
    search_rgn[:, C.w_ERBB3] = [0.01, 100]
    search_rgn[:, C.w_ERBB4] = [0.01, 100]
    #
    search_rgn[:, C.w_GRB2] = [0.01, 100]
    #
    search_rgn[:, C.w_SHC1] = [0.01, 100]
    search_rgn[:, C.w_SHC2] = [0.01, 100]
    search_rgn[:, C.w_SHC3] = [0.01, 100]
    search_rgn[:, C.w_SHC4] = [0.01, 100]
    #
    search_rgn[:, C.w_RASA1] = [0.01, 100]
    search_rgn[:, C.w_RASA2] = [0.01, 100]
    search_rgn[:, C.w_RASA3] = [0.01, 100]
    #
    search_rgn[:, C.w_PIK3CA] = [0.01, 100]
    search_rgn[:, C.w_PIK3CB] = [0.01, 100]
    search_rgn[:, C.w_PIK3CD] = [0.01, 100]
    search_rgn[:, C.w_PIK3CG] = [0.01, 100]
    #
    search_rgn[:, C.w_PTEN] = [0.01, 100]
    #
    search_rgn[:, C.w_SOS1] = [0.01, 100]
    search_rgn[:, C.w_SOS2] = [0.01, 100]
    #
    search_rgn[:, C.w_GAB1] = [0.01, 100]
    #
    search_rgn[:, C.w_HRAS] = [0.01, 100]
    search_rgn[:, C.w_KRAS] = [0.01, 100]
    search_rgn[:, C.w_NRAS] = [0.01, 100]
    #
    search_rgn[:, C.w_ARAF] = [0.01, 100]
    search_rgn[:, C.w_BRAF] = [0.01, 100]
    search_rgn[:, C.w_RAF1] = [0.01, 100]
    #
    search_rgn[:, C.w_MAP2K1] = [0.01, 100]
    search_rgn[:, C.w_MAP2K2] = [0.01, 100]
    #
    search_rgn[:, C.w_MAPK1] = [0.01, 100]
    search_rgn[:, C.w_MAPK3] = [0.01, 100]
    #
    search_rgn[:, C.w_AKT1] = [0.01, 100]
    search_rgn[:, C.w_AKT2] = [0.01, 100]
    #
    # search_rgn[:, C.w_RPS6KA1] = [0.01, 100]
    # search_rgn[:, C.w_RPS6KA2] = [0.01, 100]
    # search_rgn[:, C.w_RPS6KA3] = [0.01, 100]
    # search_rgn[:, C.w_RPS6KA6] = [0.01, 100]
    #
    # search_rgn[:, C.w_PPP2CA] = [0.01, 100]
    # search_rgn[:, C.w_PPP2CB] = [0.01, 100]
    #
    search_rgn[:, C.w_PTPN1] = [0.01, 100]
    #
    # search_rgn[:, C.w_CREB1] = [0.01, 100]
    #
    # search_rgn[:, C.w_ELK1] = [0.01, 100]
    #
    search_rgn[:, C.w_GSK3B] = [0.01, 100]
    #
    search_rgn[:, C.w_DUSP5] = [0.01, 100]
    search_rgn[:, C.w_DUSP6] = [0.01, 100]
    search_rgn[:, C.w_DUSP7] = [0.01, 100]
    #
    # search_rgn[:, C.w_FOS] = [0.01, 100]
    #
    search_rgn[:, C.w_MYC] = [0.01, 100]
    

    search_rgn = conv_lin2log!(search_rgn, search_idx)

    return search_rgn
end


function update_param(indiv::Vector{Float64})::Tuple{Array{Float64,1},Array{Float64,1}}
    p::Vector{Float64} = param_values()
    u0::Vector{Float64} = initial_values()

    search_idx::Tuple{Array{Int64,1},Array{Int64,1}} = get_search_index()

    for (i, j) in enumerate(search_idx[1])
        @inbounds p[j] = indiv[i]
    end
    for (i, j) in enumerate(search_idx[2])
        @inbounds u0[j] = indiv[i + length(search_idx[1])]
    end

    (p, u0) = individualize!(p, u0, "MCF7_BREAST")

    # constraints --------------------------------------------------------------
    p[C.kf39] = p[C.kf38]
    p[C.kr39] = p[C.kr38]
    p[C.kf40] = p[C.kf38]
    p[C.kr40] = p[C.kr38]
    p[C.kf41] = p[C.kf38]
    p[C.kr41] = p[C.kr38]
    p[C.kf42] = p[C.kf38]
    p[C.kr42] = p[C.kr38]
    p[C.kf43] = p[C.kf38]
    p[C.kr43] = p[C.kr38]
    p[C.kf45] = p[C.kf44]
    p[C.kr45] = p[C.kr44]
    p[C.kf46] = p[C.kf44]
    p[C.kr46] = p[C.kr44]
    p[C.kf47] = p[C.kf44]
    p[C.kr47] = p[C.kr44]
    p[C.kf48] = p[C.kf44]
    p[C.kr48] = p[C.kr44]
    p[C.kf49] = p[C.kf44]
    p[C.kr49] = p[C.kr44]
    p[C.kf51] = p[C.kf50]
    p[C.kr51] = p[C.kr50]
    p[C.kf52] = p[C.kf50]
    p[C.kr52] = p[C.kr50]
    p[C.kf53] = p[C.kf50]
    p[C.kr53] = p[C.kr50]
    p[C.kf54] = p[C.kf50]
    p[C.kr54] = p[C.kr50]
    p[C.kf55] = p[C.kf50]
    p[C.kr55] = p[C.kr50]
    p[C.kf57] = p[C.kf56]
    p[C.kr57] = p[C.kr56]
    p[C.kf58] = p[C.kf56]
    p[C.kr58] = p[C.kr56]
    p[C.kf59] = p[C.kf56]
    p[C.kr59] = p[C.kr56]
    p[C.kf60] = p[C.kf56]
    p[C.kr60] = p[C.kr56]
    p[C.kf61] = p[C.kf56]
    p[C.kr61] = p[C.kr56]
    p[C.kf63] = p[C.kf62]
    p[C.kr63] = p[C.kr62]
    p[C.kf64] = p[C.kf62]
    p[C.kr64] = p[C.kr62]
    p[C.kf65] = p[C.kf62]
    p[C.kr65] = p[C.kr62]
    p[C.kf66] = p[C.kf62]
    p[C.kr66] = p[C.kr62]
    p[C.kf67] = p[C.kf62]
    p[C.kr67] = p[C.kr62]
    p[C.kf68] = p[C.kf62]
    p[C.kr68] = p[C.kr62]
    p[C.kf69] = p[C.kf62]
    p[C.kr69] = p[C.kr62]
    p[C.kf70] = p[C.kf62]
    p[C.kr70] = p[C.kr62]
    p[C.kf71] = p[C.kf62]
    p[C.kr71] = p[C.kr62]
    p[C.kf72] = p[C.kf62]
    p[C.kr72] = p[C.kr62]
    p[C.kf73] = p[C.kf62]
    p[C.kr73] = p[C.kr62]
    p[C.kf75] = p[C.kf74]
    p[C.kr75] = p[C.kr74]
    p[C.kf76] = p[C.kf74]
    p[C.kr76] = p[C.kr74]
    p[C.kf77] = p[C.kf74]
    p[C.kr77] = p[C.kr74]
    p[C.kf78] = p[C.kf74]
    p[C.kr78] = p[C.kr74]
    p[C.kf79] = p[C.kf74]
    p[C.kr79] = p[C.kr74]
    p[C.kf80] = p[C.kf74]
    p[C.kr80] = p[C.kr74]
    p[C.kf81] = p[C.kf74]
    p[C.kr81] = p[C.kr74]
    p[C.kf82] = p[C.kf74]
    p[C.kr82] = p[C.kr74]
    p[C.kf83] = p[C.kf74]
    p[C.kr83] = p[C.kr74]
    p[C.kf84] = p[C.kf74]
    p[C.kr84] = p[C.kr74]
    p[C.kf85] = p[C.kf74]
    p[C.kr85] = p[C.kr74]
    p[C.kf87] = p[C.kf86]
    p[C.kr87] = p[C.kr86]
    p[C.kf88] = p[C.kf86]
    p[C.kr88] = p[C.kr86]
    p[C.kf89] = p[C.kf86]
    p[C.kr89] = p[C.kr86]
    p[C.kf90] = p[C.kf86]
    p[C.kr90] = p[C.kr86]
    p[C.kf91] = p[C.kf86]
    p[C.kr91] = p[C.kr86]
    p[C.kf92] = p[C.kf86]
    p[C.kr92] = p[C.kr86]
    p[C.kf93] = p[C.kf86]
    p[C.kr93] = p[C.kr86]
    p[C.kf94] = p[C.kf86]
    p[C.kr94] = p[C.kr86]
    p[C.kf95] = p[C.kf86]
    p[C.kr95] = p[C.kr86]
    p[C.kf96] = p[C.kf86]
    p[C.kr96] = p[C.kr86]
    p[C.kf97] = p[C.kf86]
    p[C.kr97] = p[C.kr86]
    p[C.kf99] = p[C.kf98]
    p[C.kr99] = p[C.kr98]
    p[C.kf100] = p[C.kf98]
    p[C.kr100] = p[C.kr98]
    p[C.kf101] = p[C.kf98]
    p[C.kr101] = p[C.kr98]
    p[C.kf102] = p[C.kf98]
    p[C.kr102] = p[C.kr98]
    p[C.kf103] = p[C.kf98]
    p[C.kr103] = p[C.kr98]
    p[C.kf104] = p[C.kf98]
    p[C.kr104] = p[C.kr98]
    p[C.kf105] = p[C.kf98]
    p[C.kr105] = p[C.kr98]
    p[C.kf106] = p[C.kf98]
    p[C.kr106] = p[C.kr98]
    p[C.kf107] = p[C.kf98]
    p[C.kr107] = p[C.kr98]
    p[C.kf108] = p[C.kf98]
    p[C.kr108] = p[C.kr98]
    p[C.kf109] = p[C.kf98]
    p[C.kr109] = p[C.kr98]
    p[C.kf123] = p[C.kf122]
    p[C.kr123] = p[C.kr122]
    p[C.kf124] = p[C.kf122]
    p[C.kr124] = p[C.kr122]
    p[C.kf125] = p[C.kf122]
    p[C.kr125] = p[C.kr122]
    p[C.kf126] = p[C.kf122]
    p[C.kr126] = p[C.kr122]
    p[C.kf127] = p[C.kf122]
    p[C.kr127] = p[C.kr122]
    p[C.kf128] = p[C.kf122]
    p[C.kr128] = p[C.kr122]
    p[C.kf129] = p[C.kf122]
    p[C.kr129] = p[C.kr122]
    p[C.kf130] = p[C.kf122]
    p[C.kr130] = p[C.kr122]
    p[C.kf131] = p[C.kf122]
    p[C.kr131] = p[C.kr122]
    p[C.kf132] = p[C.kf122]
    p[C.kr132] = p[C.kr122]
    p[C.kf133] = p[C.kf122]
    p[C.kr133] = p[C.kr122]
    p[C.kf134] = p[C.kf122]
    p[C.kr134] = p[C.kr122]
    p[C.kf136] = p[C.kf135]
    p[C.kr136] = p[C.kr135]
    p[C.kf137] = p[C.kf135]
    p[C.kr137] = p[C.kr135]
    p[C.kf138] = p[C.kf135]
    p[C.kr138] = p[C.kr135]
    p[C.kf139] = p[C.kf135]
    p[C.kr139] = p[C.kr135]
    p[C.kf140] = p[C.kf135]
    p[C.kr140] = p[C.kr135]
    p[C.kf141] = p[C.kf135]
    p[C.kr141] = p[C.kr135]
    p[C.kf142] = p[C.kf135]
    p[C.kr142] = p[C.kr135]
    p[C.kf143] = p[C.kf135]
    p[C.kr143] = p[C.kr135]
    p[C.kf144] = p[C.kf135]
    p[C.kr144] = p[C.kr135]
    p[C.kf145] = p[C.kf135]
    p[C.kr145] = p[C.kr135]
    p[C.kf146] = p[C.kf135]
    p[C.kr146] = p[C.kr135]
    p[C.kf147] = p[C.kf135]
    p[C.kr147] = p[C.kr135]
    p[C.kf149] = p[C.kf148]
    p[C.kr149] = p[C.kr148]
    p[C.kf150] = p[C.kf148]
    p[C.kr150] = p[C.kr148]
    p[C.kf151] = p[C.kf148]
    p[C.kr151] = p[C.kr148]
    p[C.kf152] = p[C.kf148]
    p[C.kr152] = p[C.kr148]
    p[C.kf153] = p[C.kf148]
    p[C.kr153] = p[C.kr148]
    p[C.kf154] = p[C.kf148]
    p[C.kr154] = p[C.kr148]
    p[C.kf155] = p[C.kf148]
    p[C.kr155] = p[C.kr148]
    p[C.kf156] = p[C.kf148]
    p[C.kr156] = p[C.kr148]
    p[C.kf157] = p[C.kf148]
    p[C.kr157] = p[C.kr148]
    p[C.kf158] = p[C.kf148]
    p[C.kr158] = p[C.kr148]
    p[C.kf159] = p[C.kf148]
    p[C.kr159] = p[C.kr148]
    p[C.kf160] = p[C.kf148]
    p[C.kr160] = p[C.kr148]
    p[C.kf162] = p[C.kf161]
    p[C.kr162] = p[C.kr161]
    p[C.kf163] = p[C.kf161]
    p[C.kr163] = p[C.kr161]
    p[C.kf164] = p[C.kf161]
    p[C.kr164] = p[C.kr161]
    p[C.kf165] = p[C.kf161]
    p[C.kr165] = p[C.kr161]
    p[C.kf166] = p[C.kf161]
    p[C.kr166] = p[C.kr161]
    p[C.kf167] = p[C.kf161]
    p[C.kr167] = p[C.kr161]
    p[C.kf168] = p[C.kf161]
    p[C.kr168] = p[C.kr161]
    p[C.kf169] = p[C.kf161]
    p[C.kr169] = p[C.kr161]
    p[C.kf170] = p[C.kf161]
    p[C.kr170] = p[C.kr161]
    p[C.kf171] = p[C.kf161]
    p[C.kr171] = p[C.kr161]
    p[C.kf172] = p[C.kf161]
    p[C.kr172] = p[C.kr161]
    p[C.kf173] = p[C.kf161]
    p[C.kr173] = p[C.kr161]
    p[C.kf175] = p[C.kf174]
    p[C.kr175] = p[C.kr174]
    p[C.kf176] = p[C.kf174]
    p[C.kr176] = p[C.kr174]
    p[C.kf177] = p[C.kf174]
    p[C.kr177] = p[C.kr174]
    p[C.kf178] = p[C.kf174]
    p[C.kr178] = p[C.kr174]
    p[C.kf179] = p[C.kf174]
    p[C.kr179] = p[C.kr174]
    p[C.kf180] = p[C.kf174]
    p[C.kr180] = p[C.kr174]
    p[C.kf181] = p[C.kf174]
    p[C.kr181] = p[C.kr174]
    p[C.kf182] = p[C.kf174]
    p[C.kr182] = p[C.kr174]
    p[C.kf183] = p[C.kf174]
    p[C.kr183] = p[C.kr174]
    p[C.kf184] = p[C.kf174]
    p[C.kr184] = p[C.kr174]
    p[C.kf185] = p[C.kf174]
    p[C.kr185] = p[C.kr174]
    p[C.kf186] = p[C.kf174]
    p[C.kr186] = p[C.kr174]
    p[C.kf187] = p[C.kf174]
    p[C.kr187] = p[C.kr174]
    p[C.kf188] = p[C.kf174]
    p[C.kr188] = p[C.kr174]
    p[C.kf189] = p[C.kf174]
    p[C.kr189] = p[C.kr174]
    p[C.kf190] = p[C.kf174]
    p[C.kr190] = p[C.kr174]
    p[C.kf191] = p[C.kf174]
    p[C.kr191] = p[C.kr174]
    p[C.kf192] = p[C.kf174]
    p[C.kr192] = p[C.kr174]
    p[C.kf194] = p[C.kf193]
    p[C.kr194] = p[C.kr193]
    p[C.kf195] = p[C.kf193]
    p[C.kr195] = p[C.kr193]
    p[C.kf196] = p[C.kf193]
    p[C.kr196] = p[C.kr193]
    p[C.kf197] = p[C.kf193]
    p[C.kr197] = p[C.kr193]
    p[C.kf198] = p[C.kf193]
    p[C.kr198] = p[C.kr193]
    p[C.kf199] = p[C.kf193]
    p[C.kr199] = p[C.kr193]
    p[C.kf200] = p[C.kf193]
    p[C.kr200] = p[C.kr193]
    p[C.kf201] = p[C.kf193]
    p[C.kr201] = p[C.kr193]
    p[C.kf202] = p[C.kf193]
    p[C.kr202] = p[C.kr193]
    p[C.kf203] = p[C.kf193]
    p[C.kr203] = p[C.kr193]
    p[C.kf204] = p[C.kf193]
    p[C.kr204] = p[C.kr193]
    p[C.kf205] = p[C.kf193]
    p[C.kr205] = p[C.kr193]
    p[C.kf206] = p[C.kf193]
    p[C.kr206] = p[C.kr193]
    p[C.kf207] = p[C.kf193]
    p[C.kr207] = p[C.kr193]
    p[C.kf208] = p[C.kf193]
    p[C.kr208] = p[C.kr193]
    p[C.kf209] = p[C.kf193]
    p[C.kr209] = p[C.kr193]
    p[C.kf210] = p[C.kf193]
    p[C.kr210] = p[C.kr193]
    p[C.kf211] = p[C.kf193]
    p[C.kr211] = p[C.kr193]
    p[C.kf213] = p[C.kf212]
    p[C.kr213] = p[C.kr212]
    p[C.kf214] = p[C.kf212]
    p[C.kr214] = p[C.kr212]
    p[C.kf215] = p[C.kf212]
    p[C.kr215] = p[C.kr212]
    p[C.kf216] = p[C.kf212]
    p[C.kr216] = p[C.kr212]
    p[C.kf217] = p[C.kf212]
    p[C.kr217] = p[C.kr212]
    p[C.kf218] = p[C.kf212]
    p[C.kr218] = p[C.kr212]
    p[C.kf219] = p[C.kf212]
    p[C.kr219] = p[C.kr212]
    p[C.kf220] = p[C.kf212]
    p[C.kr220] = p[C.kr212]
    p[C.kf221] = p[C.kf212]
    p[C.kr221] = p[C.kr212]
    p[C.kf222] = p[C.kf212]
    p[C.kr222] = p[C.kr212]
    p[C.kf223] = p[C.kf212]
    p[C.kr223] = p[C.kr212]
    p[C.kf224] = p[C.kf212]
    p[C.kr224] = p[C.kr212]
    p[C.kf226] = p[C.kf225]
    p[C.kr226] = p[C.kr225]
    p[C.kf227] = p[C.kf225]
    p[C.kr227] = p[C.kr225]
    p[C.kf228] = p[C.kf225]
    p[C.kr228] = p[C.kr225]
    p[C.kf229] = p[C.kf225]
    p[C.kr229] = p[C.kr225]
    p[C.kf230] = p[C.kf225]
    p[C.kr230] = p[C.kr225]
    p[C.kf231] = p[C.kf225]
    p[C.kr231] = p[C.kr225]
    p[C.kf232] = p[C.kf225]
    p[C.kr232] = p[C.kr225]
    p[C.kf233] = p[C.kf225]
    p[C.kr233] = p[C.kr225]
    p[C.kf234] = p[C.kf225]
    p[C.kr234] = p[C.kr225]
    p[C.kf235] = p[C.kf225]
    p[C.kr235] = p[C.kr225]
    p[C.kf236] = p[C.kf225]
    p[C.kr236] = p[C.kr225]
    p[C.kf237] = p[C.kf225]
    p[C.kr237] = p[C.kr225]
    p[C.kf238] = p[C.kf225]
    p[C.kr238] = p[C.kr225]
    p[C.kf239] = p[C.kf225]
    p[C.kr239] = p[C.kr225]
    p[C.kf240] = p[C.kf225]
    p[C.kr240] = p[C.kr225]
    p[C.kf241] = p[C.kf225]
    p[C.kr241] = p[C.kr225]
    p[C.kf243] = p[C.kf242]
    p[C.kr243] = p[C.kr242]
    p[C.kf244] = p[C.kf242]
    p[C.kr244] = p[C.kr242]
    p[C.kf245] = p[C.kf242]
    p[C.kr245] = p[C.kr242]
    p[C.kf246] = p[C.kf242]
    p[C.kr246] = p[C.kr242]
    p[C.kf247] = p[C.kf242]
    p[C.kr247] = p[C.kr242]
    p[C.kf248] = p[C.kf242]
    p[C.kr248] = p[C.kr242]
    p[C.kf249] = p[C.kf242]
    p[C.kr249] = p[C.kr242]
    p[C.kf250] = p[C.kf242]
    p[C.kr250] = p[C.kr242]
    p[C.kf251] = p[C.kf242]
    p[C.kr251] = p[C.kr242]
    p[C.kf252] = p[C.kf242]
    p[C.kr252] = p[C.kr242]
    p[C.kf253] = p[C.kf242]
    p[C.kr253] = p[C.kr242]
    p[C.kf254] = p[C.kf242]
    p[C.kr254] = p[C.kr242]
    p[C.kf255] = p[C.kf242]
    p[C.kr255] = p[C.kr242]
    p[C.kf256] = p[C.kf242]
    p[C.kr256] = p[C.kr242]
    p[C.kf257] = p[C.kf242]
    p[C.kr257] = p[C.kr242]
    p[C.kf258] = p[C.kf242]
    p[C.kr258] = p[C.kr242]
    p[C.kf269] = p[C.kf268]
    p[C.kf270] = p[C.kf268]
    p[C.kf271] = p[C.kf268]
    p[C.kf272] = p[C.kf268]
    p[C.kf273] = p[C.kf268]
    p[C.kf274] = p[C.kf268]
    p[C.kf275] = p[C.kf268]
    p[C.kf276] = p[C.kf268]
    p[C.kf277] = p[C.kf268]
    p[C.kf278] = p[C.kf268]
    p[C.kf279] = p[C.kf268]
    p[C.kf280] = p[C.kf268]
    p[C.kf281] = p[C.kf268]
    p[C.kf282] = p[C.kf268]
    p[C.kf283] = p[C.kf268]
    p[C.kf284] = p[C.kf268]
    p[C.kf285] = p[C.kf268]
    p[C.kf286] = p[C.kf268]
    p[C.kf287] = p[C.kf268]
    p[C.kf288] = p[C.kf268]
    p[C.kf289] = p[C.kf268]
    p[C.kf290] = p[C.kf268]
    p[C.kf303] = p[C.kf298]
    p[C.kr303] = p[C.kr298]
    p[C.kf304] = p[C.kf299]
    p[C.kr304] = p[C.kr299]
    p[C.kf305] = p[C.kf300]
    p[C.kr305] = p[C.kr300]
    p[C.kf306] = p[C.kf301]
    p[C.kr306] = p[C.kr301]
    p[C.kf307] = p[C.kf302]
    p[C.kr307] = p[C.kr302]
    # --------------------------------------------------------------------------

    return p, u0
end


function decode_gene2val(indiv_gene::Vector{Float64})::Vector{Float64}
    search_rgn::Matrix{Float64} = get_search_region()
    indiv::Vector{Float64} = zeros(length(indiv_gene))

    for i in eachindex(indiv_gene)
        indiv[i] = 10^(
            indiv_gene[i] * (
                search_rgn[2,i] - search_rgn[1,i]
            ) + search_rgn[1,i]
        )
    end

    return round.(indiv, sigdigits=7)
end


function encode_val2gene(indiv::Vector{Float64})
    search_rgn::Matrix{Float64} = get_search_region()
    indiv_gene::Vector{Float64} = zeros(length(indiv))

    for i in eachindex(indiv)
        indiv_gene[i] = (
            log10(indiv[i]) - search_rgn[1,i]
        ) / (
            search_rgn[2,i] - search_rgn[1,i]
        )
    end

    return indiv_gene
end


function encode_bestIndivVal2randGene(
        idx::Int64,
        best_indiv::Vector{Float64},
        p0_bounds::Vector{Float64})::Float64
    search_rgn::Matrix{Float64} = get_search_region()
    rand_gene::Float64 = (
        log10(
            best_indiv[idx] * 10^(
                rand() * log10(p0_bounds[2] / p0_bounds[1]) + log10(p0_bounds[1])
            )
        ) - search_rgn[1,idx]
    ) / (
        search_rgn[2,idx] - search_rgn[1,idx]
    )
    return rand_gene
end


function init_search_param(
        search_idx::Tuple{Array{Int64,1},Array{Int64,1}},
        p::Vector{Float64},
        u0::Vector{Float64})::Vector{Float64}
    duplicate::Vector{String} = []
    if length(search_idx[1]) != length(unique(search_idx[1]))
        for idx in findall([count(x -> x == i, search_idx[1])
                            for i in unique(search_idx[1])] .!= 1)
            push!(duplicate, C.NAMES[search_idx[1][idx]])
        end
        error(
            "Duplicate parameters (C.): $duplicate"
        )
    elseif length(search_idx[2]) != length(unique(search_idx[2]))
        for idx in findall([count(x -> x == i, search_idx[2])
                            for i in unique(search_idx[2])] .!= 1)
            push!(duplicate, V.NAMES[search_idx[2][idx]])
        end
        error(
            "Duplicate species (V.): $duplicate"
        )
    end
    search_param = zeros(
        length(search_idx[1]) + length(search_idx[2])
    )
    for (i, j) in enumerate(search_idx[1])
        @inbounds search_param[i] = p[j]
    end
    for (i, j) in enumerate(search_idx[2])
        @inbounds search_param[i + length(search_idx[1])] = u0[j]
    end

    if any(x -> x == 0.0, search_param)
        msg::String = "search_param must not contain zero."
        for idx in search_idx[1]
            if p[idx] == 0.0
                error(
                    @sprintf(
                        "`C.%s` in search_idx_params: ", C.NAMES[idx]
                    ) * msg
                )
            end
        end
        for idx in search_idx[2]
            if u0[idx] == 0.0
                error(
                    @sprintf(
                        "`V.%s` in search_idx_initials: ", V.NAMES[idx]
                    ) * msg
                )
            end
        end
    end

    return search_param
end


function conv_lin2log!(
        search_rgn::Matrix{Float64},
        search_idx::Tuple{Array{Int64,1},Array{Int64,1}})::Matrix{Float64}
    for i = 1:size(search_rgn, 2)
        if minimum(search_rgn[:,i]) < 0.0
            msg = "search_rgn[lower_bound,upper_bound] must be positive.\n"
            if i <= C.NUM
                error(@sprintf("`C.%s` ", C.NAMES[i]) * msg)
            else
                error(@sprintf("`V.%s` ", V.NAMES[i - C.NUM]) * msg)
            end
        elseif minimum(search_rgn[:,i]) == 0.0 && maximum(search_rgn[:,i]) != 0.0
            msg = "lower_bound must be larger than 0.\n"
            if i <= C.NUM
                error(@sprintf("`C.%s` ", C.NAMES[i]) * msg)
            else
                error(@sprintf("`V.%s` ", V.NAMES[i - C.NUM]) * msg)
            end
        elseif search_rgn[2,i] - search_rgn[1,i] < 0.0
            msg = "lower_bound must be smaller than upper_bound.\n"
            if i <= C.NUM
                error(@sprintf("`C.%s` ", C.NAMES[i]) * msg)
            else
                error(@sprintf("`V.%s` ", V.NAMES[i - C.NUM]) * msg)
            end
        end
    end

    nonzero_idx::Vector{Int} = []
    for i = 1:size(search_rgn, 2)
        if search_rgn[:,i] != [0.0,0.0]
            push!(nonzero_idx, i)
        end
    end
    difference::Vector{Int} = collect(
        symdiff(
            Set(nonzero_idx),
            Set(append!(search_idx[1], C.NUM .+ search_idx[2]))
        )
    )
    if length(difference) > 0
        for idx in difference
            if idx <= C.NUM
                println(@sprintf("`C.%s`", C.NAMES[Int(idx)]))
            else
                println(@sprintf("`V.%s`", V.NAMES[Int(idx) - C.NUM]))
            end
        end
        error(
            "Set these search_params in both search_idx and search_rgn."
        )
    end

    search_rgn = search_rgn[:,nonzero_idx]

    return log10.(search_rgn)
end