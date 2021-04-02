import os

import numpy as np
from dyaus_dev import Individualization

from . import __path__
from .name2idx import C, V
from .set_model import initial_values, param_values

incorporating_gene_expression_levels = Individualization(
    parameters=C.NAMES,
    species=V.NAMES,
    tpm_values="transcriptomic_data/TPM_RLE_postComBat.tar.xz",
    structure={
        "ErbB1": ["EGFR"],
        "ErbB2": ["ERBB2"],
        "ErbB3": ["ERBB3"],
        "ErbB4": ["ERBB4"],
        "Grb2": ["GRB2"],
        "Shc": ["SHC1", "SHC2", "SHC3", "SHC4"],
        "RasGAP": ["RASA1", "RASA2", "RASA3"],
        "PI3K": ["PIK3CA", "PIK3CB", "PIK3CD", "PIK3CG"],
        "PTEN": ["PTEN"],
        "SOS": ["SOS1", "SOS2"],
        "Gab1": ["GAB1"],
        "RasGDP": ["HRAS", "KRAS", "NRAS"],
        "Raf": ["ARAF", "BRAF", "RAF1"],
        "MEK": ["MAP2K1", "MAP2K2"],
        "ERK": ["MAPK1", "MAPK3"],
        "Akt": ["AKT1", "AKT2"],
        "PTP1B": ["PTPN1"],
        "GSK3b": ["GSK3B"],
        "DUSP": ["DUSP5", "DUSP6", "DUSP7"],
        "cMyc": ["MYC"],
    },
)


class SearchParam(object):
    """Specify model parameters and/or initial values to optimize"""

    # parameters
    idx_params = [
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
        # C.w_RPS6KA1,
        # C.w_RPS6KA2,
        # C.w_RPS6KA3,
        # C.w_RPS6KA6,
        #
        # C.w_PPP2CA,
        # C.w_PPP2CB,
        #
        C.w_PTPN1,
        #
        # C.w_CREB1,
        #
        # C.w_ELK1,
        #
        C.w_GSK3B,
        #
        C.w_DUSP5,
        C.w_DUSP6,
        C.w_DUSP7,
        #
        # C.w_FOS,
        #
        C.w_MYC,
    ]
    # initial values
    idx_initials = [
        V.PIP2,
    ]

    def get_region(self):
        x = param_values()
        y0 = initial_values()

        search_param = self._init_search_param(x, y0)

        search_rgn = np.zeros((2, len(x) + len(y0)))
        # Default: 0.1 ~ 10
        for i, j in enumerate(self.idx_params):
            search_rgn[0, j] = 1e-4  # lower bound
            search_rgn[1, j] = 1e5  # upper bound
        # Default: 0.5 ~ 2
        for i, j in enumerate(self.idx_initials):
            search_rgn[0, j + len(x)] = search_param[i + len(self.idx_params)] * 0.5  # lower bound
            search_rgn[1, j + len(x)] = search_param[i + len(self.idx_params)] * 2.0  # upper bound

        # search_rgn[:,C.parameter] = [lower_bound, upper_bound]
        # search_rgn[:,V.specie+len(x)] = [lower_bound, upper_bound]

        search_rgn = self._conv_lin2log(search_rgn)

        return search_rgn

    def update(self, indiv):
        x = param_values()
        y0 = initial_values()

        for i, j in enumerate(self.idx_params):
            x[j] = indiv[i]
        for i, j in enumerate(self.idx_initials):
            y0[j] = indiv[i + len(self.idx_params)]

        x[C.V291] = incorporating_gene_expression_levels.as_reaction_rate(
            __path__[0].split(os.sep)[-1].replace("_", "."), x, "V291", "DUSP"
        )
        x[C.V310] = incorporating_gene_expression_levels.as_reaction_rate(
            __path__[0].split(os.sep)[-1].replace("_", "."), x, "V310", "cMyc"
        )
        y0 = incorporating_gene_expression_levels.as_initial_conditions(
            __path__[0].split(os.sep)[-1].replace("_", "."), x, y0
        )

        # parameter constraints
        x[C.kf39] = x[C.kf38]
        x[C.kr39] = x[C.kr38]
        x[C.kf40] = x[C.kf38]
        x[C.kr40] = x[C.kr38]
        x[C.kf41] = x[C.kf38]
        x[C.kr41] = x[C.kr38]
        x[C.kf42] = x[C.kf38]
        x[C.kr42] = x[C.kr38]
        x[C.kf43] = x[C.kf38]
        x[C.kr43] = x[C.kr38]
        x[C.kf45] = x[C.kf44]
        x[C.kr45] = x[C.kr44]
        x[C.kf46] = x[C.kf44]
        x[C.kr46] = x[C.kr44]
        x[C.kf47] = x[C.kf44]
        x[C.kr47] = x[C.kr44]
        x[C.kf48] = x[C.kf44]
        x[C.kr48] = x[C.kr44]
        x[C.kf49] = x[C.kf44]
        x[C.kr49] = x[C.kr44]
        x[C.kf51] = x[C.kf50]
        x[C.kr51] = x[C.kr50]
        x[C.kf52] = x[C.kf50]
        x[C.kr52] = x[C.kr50]
        x[C.kf53] = x[C.kf50]
        x[C.kr53] = x[C.kr50]
        x[C.kf54] = x[C.kf50]
        x[C.kr54] = x[C.kr50]
        x[C.kf55] = x[C.kf50]
        x[C.kr55] = x[C.kr50]
        x[C.kf57] = x[C.kf56]
        x[C.kr57] = x[C.kr56]
        x[C.kf58] = x[C.kf56]
        x[C.kr58] = x[C.kr56]
        x[C.kf59] = x[C.kf56]
        x[C.kr59] = x[C.kr56]
        x[C.kf60] = x[C.kf56]
        x[C.kr60] = x[C.kr56]
        x[C.kf61] = x[C.kf56]
        x[C.kr61] = x[C.kr56]
        x[C.kf63] = x[C.kf62]
        x[C.kr63] = x[C.kr62]
        x[C.kf64] = x[C.kf62]
        x[C.kr64] = x[C.kr62]
        x[C.kf65] = x[C.kf62]
        x[C.kr65] = x[C.kr62]
        x[C.kf66] = x[C.kf62]
        x[C.kr66] = x[C.kr62]
        x[C.kf67] = x[C.kf62]
        x[C.kr67] = x[C.kr62]
        x[C.kf68] = x[C.kf62]
        x[C.kr68] = x[C.kr62]
        x[C.kf69] = x[C.kf62]
        x[C.kr69] = x[C.kr62]
        x[C.kf70] = x[C.kf62]
        x[C.kr70] = x[C.kr62]
        x[C.kf71] = x[C.kf62]
        x[C.kr71] = x[C.kr62]
        x[C.kf72] = x[C.kf62]
        x[C.kr72] = x[C.kr62]
        x[C.kf73] = x[C.kf62]
        x[C.kr73] = x[C.kr62]
        x[C.kf75] = x[C.kf74]
        x[C.kr75] = x[C.kr74]
        x[C.kf76] = x[C.kf74]
        x[C.kr76] = x[C.kr74]
        x[C.kf77] = x[C.kf74]
        x[C.kr77] = x[C.kr74]
        x[C.kf78] = x[C.kf74]
        x[C.kr78] = x[C.kr74]
        x[C.kf79] = x[C.kf74]
        x[C.kr79] = x[C.kr74]
        x[C.kf80] = x[C.kf74]
        x[C.kr80] = x[C.kr74]
        x[C.kf81] = x[C.kf74]
        x[C.kr81] = x[C.kr74]
        x[C.kf82] = x[C.kf74]
        x[C.kr82] = x[C.kr74]
        x[C.kf83] = x[C.kf74]
        x[C.kr83] = x[C.kr74]
        x[C.kf84] = x[C.kf74]
        x[C.kr84] = x[C.kr74]
        x[C.kf85] = x[C.kf74]
        x[C.kr85] = x[C.kr74]
        x[C.kf87] = x[C.kf86]
        x[C.kr87] = x[C.kr86]
        x[C.kf88] = x[C.kf86]
        x[C.kr88] = x[C.kr86]
        x[C.kf89] = x[C.kf86]
        x[C.kr89] = x[C.kr86]
        x[C.kf90] = x[C.kf86]
        x[C.kr90] = x[C.kr86]
        x[C.kf91] = x[C.kf86]
        x[C.kr91] = x[C.kr86]
        x[C.kf92] = x[C.kf86]
        x[C.kr92] = x[C.kr86]
        x[C.kf93] = x[C.kf86]
        x[C.kr93] = x[C.kr86]
        x[C.kf94] = x[C.kf86]
        x[C.kr94] = x[C.kr86]
        x[C.kf95] = x[C.kf86]
        x[C.kr95] = x[C.kr86]
        x[C.kf96] = x[C.kf86]
        x[C.kr96] = x[C.kr86]
        x[C.kf97] = x[C.kf86]
        x[C.kr97] = x[C.kr86]
        x[C.kf99] = x[C.kf98]
        x[C.kr99] = x[C.kr98]
        x[C.kf100] = x[C.kf98]
        x[C.kr100] = x[C.kr98]
        x[C.kf101] = x[C.kf98]
        x[C.kr101] = x[C.kr98]
        x[C.kf102] = x[C.kf98]
        x[C.kr102] = x[C.kr98]
        x[C.kf103] = x[C.kf98]
        x[C.kr103] = x[C.kr98]
        x[C.kf104] = x[C.kf98]
        x[C.kr104] = x[C.kr98]
        x[C.kf105] = x[C.kf98]
        x[C.kr105] = x[C.kr98]
        x[C.kf106] = x[C.kf98]
        x[C.kr106] = x[C.kr98]
        x[C.kf107] = x[C.kf98]
        x[C.kr107] = x[C.kr98]
        x[C.kf108] = x[C.kf98]
        x[C.kr108] = x[C.kr98]
        x[C.kf109] = x[C.kf98]
        x[C.kr109] = x[C.kr98]
        x[C.kf123] = x[C.kf122]
        x[C.kr123] = x[C.kr122]
        x[C.kf124] = x[C.kf122]
        x[C.kr124] = x[C.kr122]
        x[C.kf125] = x[C.kf122]
        x[C.kr125] = x[C.kr122]
        x[C.kf126] = x[C.kf122]
        x[C.kr126] = x[C.kr122]
        x[C.kf127] = x[C.kf122]
        x[C.kr127] = x[C.kr122]
        x[C.kf128] = x[C.kf122]
        x[C.kr128] = x[C.kr122]
        x[C.kf129] = x[C.kf122]
        x[C.kr129] = x[C.kr122]
        x[C.kf130] = x[C.kf122]
        x[C.kr130] = x[C.kr122]
        x[C.kf131] = x[C.kf122]
        x[C.kr131] = x[C.kr122]
        x[C.kf132] = x[C.kf122]
        x[C.kr132] = x[C.kr122]
        x[C.kf133] = x[C.kf122]
        x[C.kr133] = x[C.kr122]
        x[C.kf134] = x[C.kf122]
        x[C.kr134] = x[C.kr122]
        x[C.kf136] = x[C.kf135]
        x[C.kr136] = x[C.kr135]
        x[C.kf137] = x[C.kf135]
        x[C.kr137] = x[C.kr135]
        x[C.kf138] = x[C.kf135]
        x[C.kr138] = x[C.kr135]
        x[C.kf139] = x[C.kf135]
        x[C.kr139] = x[C.kr135]
        x[C.kf140] = x[C.kf135]
        x[C.kr140] = x[C.kr135]
        x[C.kf141] = x[C.kf135]
        x[C.kr141] = x[C.kr135]
        x[C.kf142] = x[C.kf135]
        x[C.kr142] = x[C.kr135]
        x[C.kf143] = x[C.kf135]
        x[C.kr143] = x[C.kr135]
        x[C.kf144] = x[C.kf135]
        x[C.kr144] = x[C.kr135]
        x[C.kf145] = x[C.kf135]
        x[C.kr145] = x[C.kr135]
        x[C.kf146] = x[C.kf135]
        x[C.kr146] = x[C.kr135]
        x[C.kf147] = x[C.kf135]
        x[C.kr147] = x[C.kr135]
        x[C.kf149] = x[C.kf148]
        x[C.kr149] = x[C.kr148]
        x[C.kf150] = x[C.kf148]
        x[C.kr150] = x[C.kr148]
        x[C.kf151] = x[C.kf148]
        x[C.kr151] = x[C.kr148]
        x[C.kf152] = x[C.kf148]
        x[C.kr152] = x[C.kr148]
        x[C.kf153] = x[C.kf148]
        x[C.kr153] = x[C.kr148]
        x[C.kf154] = x[C.kf148]
        x[C.kr154] = x[C.kr148]
        x[C.kf155] = x[C.kf148]
        x[C.kr155] = x[C.kr148]
        x[C.kf156] = x[C.kf148]
        x[C.kr156] = x[C.kr148]
        x[C.kf157] = x[C.kf148]
        x[C.kr157] = x[C.kr148]
        x[C.kf158] = x[C.kf148]
        x[C.kr158] = x[C.kr148]
        x[C.kf159] = x[C.kf148]
        x[C.kr159] = x[C.kr148]
        x[C.kf160] = x[C.kf148]
        x[C.kr160] = x[C.kr148]
        x[C.kf162] = x[C.kf161]
        x[C.kr162] = x[C.kr161]
        x[C.kf163] = x[C.kf161]
        x[C.kr163] = x[C.kr161]
        x[C.kf164] = x[C.kf161]
        x[C.kr164] = x[C.kr161]
        x[C.kf165] = x[C.kf161]
        x[C.kr165] = x[C.kr161]
        x[C.kf166] = x[C.kf161]
        x[C.kr166] = x[C.kr161]
        x[C.kf167] = x[C.kf161]
        x[C.kr167] = x[C.kr161]
        x[C.kf168] = x[C.kf161]
        x[C.kr168] = x[C.kr161]
        x[C.kf169] = x[C.kf161]
        x[C.kr169] = x[C.kr161]
        x[C.kf170] = x[C.kf161]
        x[C.kr170] = x[C.kr161]
        x[C.kf171] = x[C.kf161]
        x[C.kr171] = x[C.kr161]
        x[C.kf172] = x[C.kf161]
        x[C.kr172] = x[C.kr161]
        x[C.kf173] = x[C.kf161]
        x[C.kr173] = x[C.kr161]
        x[C.kf175] = x[C.kf174]
        x[C.kr175] = x[C.kr174]
        x[C.kf176] = x[C.kf174]
        x[C.kr176] = x[C.kr174]
        x[C.kf177] = x[C.kf174]
        x[C.kr177] = x[C.kr174]
        x[C.kf178] = x[C.kf174]
        x[C.kr178] = x[C.kr174]
        x[C.kf179] = x[C.kf174]
        x[C.kr179] = x[C.kr174]
        x[C.kf180] = x[C.kf174]
        x[C.kr180] = x[C.kr174]
        x[C.kf181] = x[C.kf174]
        x[C.kr181] = x[C.kr174]
        x[C.kf182] = x[C.kf174]
        x[C.kr182] = x[C.kr174]
        x[C.kf183] = x[C.kf174]
        x[C.kr183] = x[C.kr174]
        x[C.kf184] = x[C.kf174]
        x[C.kr184] = x[C.kr174]
        x[C.kf185] = x[C.kf174]
        x[C.kr185] = x[C.kr174]
        x[C.kf186] = x[C.kf174]
        x[C.kr186] = x[C.kr174]
        x[C.kf187] = x[C.kf174]
        x[C.kr187] = x[C.kr174]
        x[C.kf188] = x[C.kf174]
        x[C.kr188] = x[C.kr174]
        x[C.kf189] = x[C.kf174]
        x[C.kr189] = x[C.kr174]
        x[C.kf190] = x[C.kf174]
        x[C.kr190] = x[C.kr174]
        x[C.kf191] = x[C.kf174]
        x[C.kr191] = x[C.kr174]
        x[C.kf192] = x[C.kf174]
        x[C.kr192] = x[C.kr174]
        x[C.kf194] = x[C.kf193]
        x[C.kr194] = x[C.kr193]
        x[C.kf195] = x[C.kf193]
        x[C.kr195] = x[C.kr193]
        x[C.kf196] = x[C.kf193]
        x[C.kr196] = x[C.kr193]
        x[C.kf197] = x[C.kf193]
        x[C.kr197] = x[C.kr193]
        x[C.kf198] = x[C.kf193]
        x[C.kr198] = x[C.kr193]
        x[C.kf199] = x[C.kf193]
        x[C.kr199] = x[C.kr193]
        x[C.kf200] = x[C.kf193]
        x[C.kr200] = x[C.kr193]
        x[C.kf201] = x[C.kf193]
        x[C.kr201] = x[C.kr193]
        x[C.kf202] = x[C.kf193]
        x[C.kr202] = x[C.kr193]
        x[C.kf203] = x[C.kf193]
        x[C.kr203] = x[C.kr193]
        x[C.kf204] = x[C.kf193]
        x[C.kr204] = x[C.kr193]
        x[C.kf205] = x[C.kf193]
        x[C.kr205] = x[C.kr193]
        x[C.kf206] = x[C.kf193]
        x[C.kr206] = x[C.kr193]
        x[C.kf207] = x[C.kf193]
        x[C.kr207] = x[C.kr193]
        x[C.kf208] = x[C.kf193]
        x[C.kr208] = x[C.kr193]
        x[C.kf209] = x[C.kf193]
        x[C.kr209] = x[C.kr193]
        x[C.kf210] = x[C.kf193]
        x[C.kr210] = x[C.kr193]
        x[C.kf211] = x[C.kf193]
        x[C.kr211] = x[C.kr193]
        x[C.kf213] = x[C.kf212]
        x[C.kr213] = x[C.kr212]
        x[C.kf214] = x[C.kf212]
        x[C.kr214] = x[C.kr212]
        x[C.kf215] = x[C.kf212]
        x[C.kr215] = x[C.kr212]
        x[C.kf216] = x[C.kf212]
        x[C.kr216] = x[C.kr212]
        x[C.kf217] = x[C.kf212]
        x[C.kr217] = x[C.kr212]
        x[C.kf218] = x[C.kf212]
        x[C.kr218] = x[C.kr212]
        x[C.kf219] = x[C.kf212]
        x[C.kr219] = x[C.kr212]
        x[C.kf220] = x[C.kf212]
        x[C.kr220] = x[C.kr212]
        x[C.kf221] = x[C.kf212]
        x[C.kr221] = x[C.kr212]
        x[C.kf222] = x[C.kf212]
        x[C.kr222] = x[C.kr212]
        x[C.kf223] = x[C.kf212]
        x[C.kr223] = x[C.kr212]
        x[C.kf224] = x[C.kf212]
        x[C.kr224] = x[C.kr212]
        x[C.kf226] = x[C.kf225]
        x[C.kr226] = x[C.kr225]
        x[C.kf227] = x[C.kf225]
        x[C.kr227] = x[C.kr225]
        x[C.kf228] = x[C.kf225]
        x[C.kr228] = x[C.kr225]
        x[C.kf229] = x[C.kf225]
        x[C.kr229] = x[C.kr225]
        x[C.kf230] = x[C.kf225]
        x[C.kr230] = x[C.kr225]
        x[C.kf231] = x[C.kf225]
        x[C.kr231] = x[C.kr225]
        x[C.kf232] = x[C.kf225]
        x[C.kr232] = x[C.kr225]
        x[C.kf233] = x[C.kf225]
        x[C.kr233] = x[C.kr225]
        x[C.kf234] = x[C.kf225]
        x[C.kr234] = x[C.kr225]
        x[C.kf235] = x[C.kf225]
        x[C.kr235] = x[C.kr225]
        x[C.kf236] = x[C.kf225]
        x[C.kr236] = x[C.kr225]
        x[C.kf237] = x[C.kf225]
        x[C.kr237] = x[C.kr225]
        x[C.kf238] = x[C.kf225]
        x[C.kr238] = x[C.kr225]
        x[C.kf239] = x[C.kf225]
        x[C.kr239] = x[C.kr225]
        x[C.kf240] = x[C.kf225]
        x[C.kr240] = x[C.kr225]
        x[C.kf241] = x[C.kf225]
        x[C.kr241] = x[C.kr225]
        x[C.kf243] = x[C.kf242]
        x[C.kr243] = x[C.kr242]
        x[C.kf244] = x[C.kf242]
        x[C.kr244] = x[C.kr242]
        x[C.kf245] = x[C.kf242]
        x[C.kr245] = x[C.kr242]
        x[C.kf246] = x[C.kf242]
        x[C.kr246] = x[C.kr242]
        x[C.kf247] = x[C.kf242]
        x[C.kr247] = x[C.kr242]
        x[C.kf248] = x[C.kf242]
        x[C.kr248] = x[C.kr242]
        x[C.kf249] = x[C.kf242]
        x[C.kr249] = x[C.kr242]
        x[C.kf250] = x[C.kf242]
        x[C.kr250] = x[C.kr242]
        x[C.kf251] = x[C.kf242]
        x[C.kr251] = x[C.kr242]
        x[C.kf252] = x[C.kf242]
        x[C.kr252] = x[C.kr242]
        x[C.kf253] = x[C.kf242]
        x[C.kr253] = x[C.kr242]
        x[C.kf254] = x[C.kf242]
        x[C.kr254] = x[C.kr242]
        x[C.kf255] = x[C.kf242]
        x[C.kr255] = x[C.kr242]
        x[C.kf256] = x[C.kf242]
        x[C.kr256] = x[C.kr242]
        x[C.kf257] = x[C.kf242]
        x[C.kr257] = x[C.kr242]
        x[C.kf258] = x[C.kf242]
        x[C.kr258] = x[C.kr242]
        x[C.kf269] = x[C.kf268]
        x[C.kf270] = x[C.kf268]
        x[C.kf271] = x[C.kf268]
        x[C.kf272] = x[C.kf268]
        x[C.kf273] = x[C.kf268]
        x[C.kf274] = x[C.kf268]
        x[C.kf275] = x[C.kf268]
        x[C.kf276] = x[C.kf268]
        x[C.kf277] = x[C.kf268]
        x[C.kf278] = x[C.kf268]
        x[C.kf279] = x[C.kf268]
        x[C.kf280] = x[C.kf268]
        x[C.kf281] = x[C.kf268]
        x[C.kf282] = x[C.kf268]
        x[C.kf283] = x[C.kf268]
        x[C.kf284] = x[C.kf268]
        x[C.kf285] = x[C.kf268]
        x[C.kf286] = x[C.kf268]
        x[C.kf287] = x[C.kf268]
        x[C.kf288] = x[C.kf268]
        x[C.kf289] = x[C.kf268]
        x[C.kf290] = x[C.kf268]
        x[C.kf303] = x[C.kf298]
        x[C.kr303] = x[C.kr298]
        x[C.kf304] = x[C.kf299]
        x[C.kr304] = x[C.kr299]
        x[C.kf305] = x[C.kf300]
        x[C.kr305] = x[C.kr300]
        x[C.kf306] = x[C.kf301]
        x[C.kr306] = x[C.kr301]
        x[C.kf307] = x[C.kf302]
        x[C.kr307] = x[C.kr302]

        return x, y0

    def gene2val(self, indiv_gene):
        search_rgn = self.get_region()
        indiv = 10 ** (indiv_gene * (search_rgn[1, :] - search_rgn[0, :]) + search_rgn[0, :])

        return indiv

    def val2gene(self, indiv):
        search_rgn = self.get_region()
        indiv_gene = (np.log10(indiv) - search_rgn[0, :]) / (search_rgn[1, :] - search_rgn[0, :])

        return indiv_gene

    def _init_search_param(self, x, y0):
        """Initialize search_param"""
        if len(self.idx_params) != len(set(self.idx_params)):
            raise ValueError(
                "Duplicate parameters (C.): {}".format(
                    [
                        C.NAMES[idx]
                        for idx in [
                            name
                            for name in set(self.idx_params)
                            if self.idx_params.count(name) > 1
                        ]
                    ]
                )
            )
        elif len(self.idx_initials) != len(set(self.idx_initials)):
            raise ValueError(
                "Duplicate species (V.): {}".format(
                    [
                        V.NAMES[idx]
                        for idx in [
                            name
                            for name in set(self.idx_initials)
                            if self.idx_initials.count(name) > 1
                        ]
                    ]
                )
            )
        search_param = np.empty(len(self.idx_params) + len(self.idx_initials))
        for i, j in enumerate(self.idx_params):
            search_param[i] = x[j]
        for i, j in enumerate(self.idx_initials):
            search_param[i + len(self.idx_params)] = y0[j]

        if np.any(search_param == 0.0):
            message = "search_param must not contain zero."
            for idx in self.idx_params:
                if x[int(idx)] == 0.0:
                    raise ValueError('"C.{}" in idx_params: '.format(C.NAMES[int(idx)]) + message)
            for idx in self.idx_initials:
                if y0[int(idx)] == 0.0:
                    raise ValueError(
                        '"V.{}" in idx_initials: '.format(V.NAMES[int(idx)]) + message
                    )

        return search_param

    def _conv_lin2log(self, search_rgn):
        """Convert Linear scale to Logarithmic scale"""
        for i in range(search_rgn.shape[1]):
            if np.min(search_rgn[:, i]) < 0.0:
                msg = "search_rgn[lower_bound, upper_bound] must be positive."
                if i <= C.NUM:
                    raise ValueError('"C.{}": '.format(C.NAMES[i]) + msg)
                else:
                    raise ValueError('"V.{}": '.format(V.NAMES[i - C.NUM]) + msg)
            elif np.min(search_rgn[:, i]) == 0 and np.max(search_rgn[:, i]) > 0:
                msg = "lower_bound must be larger than 0."
                if i <= C.NUM:
                    raise ValueError('"C.{}" '.format(C.NAMES[i]) + msg)
                else:
                    raise ValueError('"V.{}" '.format(V.NAMES[i - C.NUM]) + msg)
            elif search_rgn[1, i] - search_rgn[0, i] < 0.0:
                msg = "lower_bound must be smaller than upper_bound."
                if i <= C.NUM:
                    raise ValueError('"C.{}" : '.format(C.NAMES[i]) + msg)
                else:
                    raise ValueError('"V.{}" : '.format(V.NAMES[i - C.NUM]) + msg)
        difference = list(
            set(np.where(np.any(search_rgn != 0.0, axis=0))[0])
            ^ set(np.append(self.idx_params, [C.NUM + idx for idx in self.idx_initials]))
        )
        if len(difference) > 0:
            msg = "in both search_idx and search_rgn"
            for idx in difference:
                if idx <= C.NUM:
                    raise ValueError('Set "C.{}" '.format(C.NAMES[int(idx)]) + msg)
                else:
                    raise ValueError('Set "V.{}" '.format(V.NAMES[int(idx - C.NUM)]) + msg)
        search_rgn = search_rgn[:, np.any(search_rgn != 0.0, axis=0)]

        return np.log10(search_rgn)
