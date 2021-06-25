from .name2idx import C, V


class DifferentialEquation(object):
    def __init__(self, perturbation):
        super(DifferentialEquation, self).__init__()
        self.perturbation = perturbation

    def diffeq(self, t, y, *x):
        """Kinetic equations"""
        # v : flux vector
        v = {}
        v[1] = x[C.kf1] * y[V.EGF] * y[V.ErbB1] - x[C.kr1] * y[V.EGF_ErbB1]
        v[2] = x[C.kf2] * y[V.HRG] * y[V.ErbB3] - x[C.kr2] * y[V.HRG_ErbB3]
        v[3] = x[C.kf3] * y[V.HRG] * y[V.ErbB4] - x[C.kr3] * y[V.HRG_ErbB4]
        v[4] = x[C.kf4] * y[V.EGF_ErbB1] * y[V.EGF_ErbB1] - x[C.kr4] * y[V.E11]
        v[5] = x[C.kf5] * y[V.EGF_ErbB1] * y[V.ErbB2] - x[C.kr5] * y[V.E12]
        v[6] = x[C.kf6] * y[V.HRG_ErbB3] * y[V.ErbB2] - x[C.kr6] * y[V.E23]
        v[7] = x[C.kf7] * y[V.HRG_ErbB3] * y[V.HRG_ErbB4] - x[C.kr7] * y[V.E34]
        v[8] = x[C.kf8] * y[V.HRG_ErbB4] * y[V.ErbB2] - x[C.kr8] * y[V.E24]
        v[9] = x[C.kf9] * y[V.HRG_ErbB4] * y[V.HRG_ErbB4] - x[C.kr9] * y[V.E44]
        v[10] = x[C.kf10] * y[V.E11] - x[C.kr10] * y[V.E11P]
        v[11] = x[C.kf11] * y[V.E12] - x[C.kr11] * y[V.E12P]
        v[12] = x[C.kf12] * y[V.E23] - x[C.kr12] * y[V.E23P]
        v[13] = x[C.kf13] * y[V.E24] - x[C.kr13] * y[V.E24P]
        v[14] = x[C.kf14] * y[V.E34] - x[C.kr14] * y[V.E34P]
        v[15] = x[C.kf15] * y[V.E44] - x[C.kr15] * y[V.E44P]
        v[16] = x[C.kf16] * y[V.E11P] * y[V.Grb2] - x[C.kr16] * y[V.E11P_Grb2]
        v[17] = x[C.kf17] * y[V.E11P] * y[V.Shc] - x[C.kr17] * y[V.E11P_Shc]
        v[18] = x[C.kf18] * y[V.E11P] * y[V.RasGAP] - x[C.kr18] * y[V.E11P_RasGAP]
        v[19] = x[C.kf19] * y[V.E12P] * y[V.Grb2] - x[C.kr19] * y[V.E12P_Grb2]
        v[20] = x[C.kf20] * y[V.E12P] * y[V.Shc] - x[C.kr20] * y[V.E12P_Shc]
        v[21] = x[C.kf21] * y[V.E12P] * y[V.RasGAP] - x[C.kr21] * y[V.E12P_RasGAP]
        v[22] = x[C.kf22] * y[V.E23P] * y[V.Grb2] - x[C.kr22] * y[V.E23P_Grb2]
        v[23] = x[C.kf23] * y[V.E23P] * y[V.Shc] - x[C.kr23] * y[V.E23P_Shc]
        v[24] = x[C.kf24] * y[V.E23P] * y[V.PI3K] - x[C.kr24] * y[V.E23P_PI3K]
        v[25] = x[C.kf25] * y[V.E23P] * y[V.RasGAP] - x[C.kr25] * y[V.E23P_RasGAP]
        v[26] = x[C.kf26] * y[V.E24P] * y[V.Grb2] - x[C.kr26] * y[V.E24P_Grb2]
        v[27] = x[C.kf27] * y[V.E24P] * y[V.Shc] - x[C.kr27] * y[V.E24P_Shc]
        v[28] = x[C.kf28] * y[V.E24P] * y[V.PI3K] - x[C.kr28] * y[V.E24P_PI3K]
        v[29] = x[C.kf29] * y[V.E24P] * y[V.RasGAP] - x[C.kr29] * y[V.E24P_RasGAP]
        v[30] = x[C.kf30] * y[V.E34P] * y[V.Grb2] - x[C.kr30] * y[V.E34P_Grb2]
        v[31] = x[C.kf31] * y[V.E34P] * y[V.Shc] - x[C.kr31] * y[V.E34P_Shc]
        v[32] = x[C.kf32] * y[V.E34P] * y[V.PI3K] - x[C.kr32] * y[V.E34P_PI3K]
        v[33] = x[C.kf33] * y[V.E34P] * y[V.RasGAP] - x[C.kr33] * y[V.E34P_RasGAP]
        v[34] = x[C.kf34] * y[V.E44P] * y[V.Grb2] - x[C.kr34] * y[V.E44P_Grb2]
        v[35] = x[C.kf35] * y[V.E44P] * y[V.Shc] - x[C.kr35] * y[V.E44P_Shc]
        v[36] = x[C.kf36] * y[V.E44P] * y[V.PI3K] - x[C.kr36] * y[V.E44P_PI3K]
        v[37] = x[C.kf37] * y[V.E44P] * y[V.RasGAP] - x[C.kr37] * y[V.E44P_RasGAP]
        v[38] = x[C.kf38] * y[V.E11P_Shc] * y[V.Grb2] - x[C.kr38] * y[V.E11P_Shc_Grb2]
        v[39] = x[C.kf39] * y[V.E12P_Shc] * y[V.Grb2] - x[C.kr39] * y[V.E12P_Shc_Grb2]
        v[40] = x[C.kf40] * y[V.E23P_Shc] * y[V.Grb2] - x[C.kr40] * y[V.E23P_Shc_Grb2]
        v[41] = x[C.kf41] * y[V.E24P_Shc] * y[V.Grb2] - x[C.kr41] * y[V.E24P_Shc_Grb2]
        v[42] = x[C.kf42] * y[V.E34P_Shc] * y[V.Grb2] - x[C.kr42] * y[V.E34P_Shc_Grb2]
        v[43] = x[C.kf43] * y[V.E44P_Shc] * y[V.Grb2] - x[C.kr43] * y[V.E44P_Shc_Grb2]
        v[44] = x[C.kf44] * y[V.E11P_Shc_Grb2] - x[C.kr44] * y[V.E11P_ShcP_Grb2]
        v[45] = x[C.kf45] * y[V.E12P_Shc_Grb2] - x[C.kr45] * y[V.E12P_ShcP_Grb2]
        v[46] = x[C.kf46] * y[V.E23P_Shc_Grb2] - x[C.kr46] * y[V.E23P_ShcP_Grb2]
        v[47] = x[C.kf47] * y[V.E24P_Shc_Grb2] - x[C.kr47] * y[V.E24P_ShcP_Grb2]
        v[48] = x[C.kf48] * y[V.E34P_Shc_Grb2] - x[C.kr48] * y[V.E34P_ShcP_Grb2]
        v[49] = x[C.kf49] * y[V.E44P_Shc_Grb2] - x[C.kr49] * y[V.E44P_ShcP_Grb2]
        v[50] = (
            x[C.kf50] * y[V.E11P_ShcP_Grb2] * y[V.PTP1B] - x[C.kr50] * y[V.E11P_ShcP_Grb2_PTP1B]
        )
        v[51] = (
            x[C.kf51] * y[V.E12P_ShcP_Grb2] * y[V.PTP1B] - x[C.kr51] * y[V.E12P_ShcP_Grb2_PTP1B]
        )
        v[52] = (
            x[C.kf52] * y[V.E23P_ShcP_Grb2] * y[V.PTP1B] - x[C.kr52] * y[V.E23P_ShcP_Grb2_PTP1B]
        )
        v[53] = (
            x[C.kf53] * y[V.E24P_ShcP_Grb2] * y[V.PTP1B] - x[C.kr53] * y[V.E24P_ShcP_Grb2_PTP1B]
        )
        v[54] = (
            x[C.kf54] * y[V.E34P_ShcP_Grb2] * y[V.PTP1B] - x[C.kr54] * y[V.E34P_ShcP_Grb2_PTP1B]
        )
        v[55] = (
            x[C.kf55] * y[V.E44P_ShcP_Grb2] * y[V.PTP1B] - x[C.kr55] * y[V.E44P_ShcP_Grb2_PTP1B]
        )
        v[56] = x[C.kf56] * y[V.E11P_ShcP_Grb2_PTP1B] - x[C.kr56] * y[V.E11P_Shc_Grb2] * y[V.PTP1B]
        v[57] = x[C.kf57] * y[V.E12P_ShcP_Grb2_PTP1B] - x[C.kr57] * y[V.E12P_Shc_Grb2] * y[V.PTP1B]
        v[58] = x[C.kf58] * y[V.E23P_ShcP_Grb2_PTP1B] - x[C.kr58] * y[V.E23P_Shc_Grb2] * y[V.PTP1B]
        v[59] = x[C.kf59] * y[V.E24P_ShcP_Grb2_PTP1B] - x[C.kr59] * y[V.E24P_Shc_Grb2] * y[V.PTP1B]
        v[60] = x[C.kf60] * y[V.E34P_ShcP_Grb2_PTP1B] - x[C.kr60] * y[V.E34P_Shc_Grb2] * y[V.PTP1B]
        v[61] = x[C.kf61] * y[V.E44P_ShcP_Grb2_PTP1B] - x[C.kr61] * y[V.E44P_Shc_Grb2] * y[V.PTP1B]
        v[62] = x[C.kf62] * y[V.E11P_Grb2] * y[V.SOS] - x[C.kr62] * y[V.E11P_Grb2_SOS]
        v[63] = x[C.kf63] * y[V.E12P_Grb2] * y[V.SOS] - x[C.kr63] * y[V.E12P_Grb2_SOS]
        v[64] = x[C.kf64] * y[V.E23P_Grb2] * y[V.SOS] - x[C.kr64] * y[V.E23P_Grb2_SOS]
        v[65] = x[C.kf65] * y[V.E24P_Grb2] * y[V.SOS] - x[C.kr65] * y[V.E24P_Grb2_SOS]
        v[66] = x[C.kf66] * y[V.E34P_Grb2] * y[V.SOS] - x[C.kr66] * y[V.E34P_Grb2_SOS]
        v[67] = x[C.kf67] * y[V.E44P_Grb2] * y[V.SOS] - x[C.kr67] * y[V.E44P_Grb2_SOS]
        v[68] = x[C.kf68] * y[V.E11P_ShcP_Grb2] * y[V.SOS] - x[C.kr68] * y[V.E11P_ShcP_Grb2_SOS]
        v[69] = x[C.kf69] * y[V.E12P_ShcP_Grb2] * y[V.SOS] - x[C.kr69] * y[V.E12P_ShcP_Grb2_SOS]
        v[70] = x[C.kf70] * y[V.E23P_ShcP_Grb2] * y[V.SOS] - x[C.kr70] * y[V.E23P_ShcP_Grb2_SOS]
        v[71] = x[C.kf71] * y[V.E24P_ShcP_Grb2] * y[V.SOS] - x[C.kr71] * y[V.E24P_ShcP_Grb2_SOS]
        v[72] = x[C.kf72] * y[V.E34P_ShcP_Grb2] * y[V.SOS] - x[C.kr72] * y[V.E34P_ShcP_Grb2_SOS]
        v[73] = x[C.kf73] * y[V.E44P_ShcP_Grb2] * y[V.SOS] - x[C.kr73] * y[V.E44P_ShcP_Grb2_SOS]
        v[74] = x[C.kf74] * y[V.E11P_Grb2] * y[V.Gab1] - x[C.kr74] * y[V.E11P_Grb2_Gab1]
        v[75] = x[C.kf75] * y[V.E12P_Grb2] * y[V.Gab1] - x[C.kr75] * y[V.E12P_Grb2_Gab1]
        v[76] = x[C.kf76] * y[V.E23P_Grb2] * y[V.Gab1] - x[C.kr76] * y[V.E23P_Grb2_Gab1]
        v[77] = x[C.kf77] * y[V.E24P_Grb2] * y[V.Gab1] - x[C.kr77] * y[V.E24P_Grb2_Gab1]
        v[78] = x[C.kf78] * y[V.E34P_Grb2] * y[V.Gab1] - x[C.kr78] * y[V.E34P_Grb2_Gab1]
        v[79] = x[C.kf79] * y[V.E44P_Grb2] * y[V.Gab1] - x[C.kr79] * y[V.E44P_Grb2_Gab1]
        v[80] = x[C.kf80] * y[V.E11P_ShcP_Grb2] * y[V.Gab1] - x[C.kr80] * y[V.E11P_ShcP_Grb2_Gab1]
        v[81] = x[C.kf81] * y[V.E12P_ShcP_Grb2] * y[V.Gab1] - x[C.kr81] * y[V.E12P_ShcP_Grb2_Gab1]
        v[82] = x[C.kf82] * y[V.E23P_ShcP_Grb2] * y[V.Gab1] - x[C.kr82] * y[V.E23P_ShcP_Grb2_Gab1]
        v[83] = x[C.kf83] * y[V.E24P_ShcP_Grb2] * y[V.Gab1] - x[C.kr83] * y[V.E24P_ShcP_Grb2_Gab1]
        v[84] = x[C.kf84] * y[V.E34P_ShcP_Grb2] * y[V.Gab1] - x[C.kr84] * y[V.E34P_ShcP_Grb2_Gab1]
        v[85] = x[C.kf85] * y[V.E44P_ShcP_Grb2] * y[V.Gab1] - x[C.kr85] * y[V.E44P_ShcP_Grb2_Gab1]
        v[86] = (
            x[C.kf86] * y[V.E11P_Grb2_SOS] * y[V.RasGDP] - x[C.kr86] * y[V.E11P_Grb2_SOS_RasGDP]
        )
        v[87] = (
            x[C.kf87] * y[V.E12P_Grb2_SOS] * y[V.RasGDP] - x[C.kr87] * y[V.E12P_Grb2_SOS_RasGDP]
        )
        v[88] = (
            x[C.kf88] * y[V.E23P_Grb2_SOS] * y[V.RasGDP] - x[C.kr88] * y[V.E23P_Grb2_SOS_RasGDP]
        )
        v[89] = (
            x[C.kf89] * y[V.E24P_Grb2_SOS] * y[V.RasGDP] - x[C.kr89] * y[V.E24P_Grb2_SOS_RasGDP]
        )
        v[90] = (
            x[C.kf90] * y[V.E34P_Grb2_SOS] * y[V.RasGDP] - x[C.kr90] * y[V.E34P_Grb2_SOS_RasGDP]
        )
        v[91] = (
            x[C.kf91] * y[V.E44P_Grb2_SOS] * y[V.RasGDP] - x[C.kr91] * y[V.E44P_Grb2_SOS_RasGDP]
        )
        v[92] = (
            x[C.kf92] * y[V.E11P_ShcP_Grb2_SOS] * y[V.RasGDP]
            - x[C.kr92] * y[V.E11P_ShcP_Grb2_SOS_RasGDP]
        )
        v[93] = (
            x[C.kf93] * y[V.E12P_ShcP_Grb2_SOS] * y[V.RasGDP]
            - x[C.kr93] * y[V.E12P_ShcP_Grb2_SOS_RasGDP]
        )
        v[94] = (
            x[C.kf94] * y[V.E23P_ShcP_Grb2_SOS] * y[V.RasGDP]
            - x[C.kr94] * y[V.E23P_ShcP_Grb2_SOS_RasGDP]
        )
        v[95] = (
            x[C.kf95] * y[V.E24P_ShcP_Grb2_SOS] * y[V.RasGDP]
            - x[C.kr95] * y[V.E24P_ShcP_Grb2_SOS_RasGDP]
        )
        v[96] = (
            x[C.kf96] * y[V.E34P_ShcP_Grb2_SOS] * y[V.RasGDP]
            - x[C.kr96] * y[V.E34P_ShcP_Grb2_SOS_RasGDP]
        )
        v[97] = (
            x[C.kf97] * y[V.E44P_ShcP_Grb2_SOS] * y[V.RasGDP]
            - x[C.kr97] * y[V.E44P_ShcP_Grb2_SOS_RasGDP]
        )
        v[98] = (
            x[C.kf98] * y[V.E11P_Grb2_SOS_RasGDP] - x[C.kr98] * y[V.E11P_Grb2_SOS] * y[V.RasGTP]
        )
        v[99] = (
            x[C.kf99] * y[V.E12P_Grb2_SOS_RasGDP] - x[C.kr99] * y[V.E12P_Grb2_SOS] * y[V.RasGTP]
        )
        v[100] = (
            x[C.kf100] * y[V.E23P_Grb2_SOS_RasGDP] - x[C.kr100] * y[V.E23P_Grb2_SOS] * y[V.RasGTP]
        )
        v[101] = (
            x[C.kf101] * y[V.E24P_Grb2_SOS_RasGDP] - x[C.kr101] * y[V.E24P_Grb2_SOS] * y[V.RasGTP]
        )
        v[102] = (
            x[C.kf102] * y[V.E34P_Grb2_SOS_RasGDP] - x[C.kr102] * y[V.E34P_Grb2_SOS] * y[V.RasGTP]
        )
        v[103] = (
            x[C.kf103] * y[V.E44P_Grb2_SOS_RasGDP] - x[C.kr103] * y[V.E44P_Grb2_SOS] * y[V.RasGTP]
        )
        v[104] = (
            x[C.kf104] * y[V.E11P_ShcP_Grb2_SOS_RasGDP]
            - x[C.kr104] * y[V.E11P_ShcP_Grb2_SOS] * y[V.RasGTP]
        )
        v[105] = (
            x[C.kf105] * y[V.E12P_ShcP_Grb2_SOS_RasGDP]
            - x[C.kr105] * y[V.E12P_ShcP_Grb2_SOS] * y[V.RasGTP]
        )
        v[106] = (
            x[C.kf106] * y[V.E23P_ShcP_Grb2_SOS_RasGDP]
            - x[C.kr106] * y[V.E23P_ShcP_Grb2_SOS] * y[V.RasGTP]
        )
        v[107] = (
            x[C.kf107] * y[V.E24P_ShcP_Grb2_SOS_RasGDP]
            - x[C.kr107] * y[V.E24P_ShcP_Grb2_SOS] * y[V.RasGTP]
        )
        v[108] = (
            x[C.kf108] * y[V.E34P_ShcP_Grb2_SOS_RasGDP]
            - x[C.kr108] * y[V.E34P_ShcP_Grb2_SOS] * y[V.RasGTP]
        )
        v[109] = (
            x[C.kf109] * y[V.E44P_ShcP_Grb2_SOS_RasGDP]
            - x[C.kr109] * y[V.E44P_ShcP_Grb2_SOS] * y[V.RasGTP]
        )
        v[110] = x[C.V110] * y[V.RasGTP] * y[V.Raf] / (x[C.K110] + y[V.Raf])
        v[111] = x[C.V111] * y[V.RafP] / (x[C.K111] + y[V.RafP])
        v[112] = x[C.V112] * y[V.RafP] * y[V.MEK] / (x[C.K112] + y[V.MEK])
        v[113] = x[C.V113] * y[V.MEKP] / (x[C.K113] + y[V.MEKP])
        v[114] = x[C.V114] * y[V.MEKP] * y[V.ERK] / (x[C.K114] + y[V.ERK])
        v[115] = x[C.V115] * y[V.MEKP] * y[V.ERKP] / (x[C.K115] + y[V.ERKP])
        v[116] = x[C.V116] * y[V.ERKP] / (x[C.K116] + y[V.ERKP])
        v[117] = x[C.V117] * y[V.ERKPP] / (x[C.K117] + y[V.ERKPP])
        v[118] = x[C.V118] * y[V.ERKPP] * y[V.SOS] / (x[C.K118] + y[V.SOS])
        v[119] = x[C.V119] * y[V.SOSP] / (x[C.K119] + y[V.SOSP])
        v[120] = x[C.V120] * y[V.ERKPP] * y[V.Gab1] / (x[C.K120] + y[V.Gab1])
        v[121] = x[C.V121] * y[V.Gab1P] / (x[C.K121] + y[V.Gab1P])
        v[122] = x[C.kf122] * y[V.E11P_Grb2_Gab1] - x[C.kr122] * y[V.E11P_Grb2_Gab1P]
        v[123] = x[C.kf123] * y[V.E12P_Grb2_Gab1] - x[C.kr123] * y[V.E12P_Grb2_Gab1P]
        v[124] = x[C.kf124] * y[V.E23P_Grb2_Gab1] - x[C.kr124] * y[V.E23P_Grb2_Gab1P]
        v[125] = x[C.kf125] * y[V.E24P_Grb2_Gab1] - x[C.kr125] * y[V.E24P_Grb2_Gab1P]
        v[126] = x[C.kf126] * y[V.E34P_Grb2_Gab1] - x[C.kr126] * y[V.E34P_Grb2_Gab1P]
        v[127] = x[C.kf127] * y[V.E44P_Grb2_Gab1] - x[C.kr127] * y[V.E44P_Grb2_Gab1P]
        v[128] = x[C.kf128] * y[V.E11P_ShcP_Grb2_Gab1] - x[C.kr128] * y[V.E11P_ShcP_Grb2_Gab1P]
        v[129] = x[C.kf129] * y[V.E12P_ShcP_Grb2_Gab1] - x[C.kr129] * y[V.E12P_ShcP_Grb2_Gab1P]
        v[130] = x[C.kf130] * y[V.E23P_ShcP_Grb2_Gab1] - x[C.kr130] * y[V.E23P_ShcP_Grb2_Gab1P]
        v[131] = x[C.kf131] * y[V.E24P_ShcP_Grb2_Gab1] - x[C.kr131] * y[V.E24P_ShcP_Grb2_Gab1P]
        v[132] = x[C.kf132] * y[V.E34P_ShcP_Grb2_Gab1] - x[C.kr132] * y[V.E34P_ShcP_Grb2_Gab1P]
        v[133] = x[C.kf133] * y[V.E44P_ShcP_Grb2_Gab1] - x[C.kr133] * y[V.E44P_ShcP_Grb2_Gab1P]
        v[134] = x[C.kf134] * y[V.PIP3_Gab1] - x[C.kr134] * y[V.PIP3_Gab1P]
        v[135] = (
            x[C.kf135] * y[V.E11P_Grb2_Gab1P] * y[V.PTP1B]
            - x[C.kr135] * y[V.E11P_Grb2_Gab1P_PTP1B]
        )
        v[136] = (
            x[C.kf136] * y[V.E12P_Grb2_Gab1P] * y[V.PTP1B]
            - x[C.kr136] * y[V.E12P_Grb2_Gab1P_PTP1B]
        )
        v[137] = (
            x[C.kf137] * y[V.E23P_Grb2_Gab1P] * y[V.PTP1B]
            - x[C.kr137] * y[V.E23P_Grb2_Gab1P_PTP1B]
        )
        v[138] = (
            x[C.kf138] * y[V.E24P_Grb2_Gab1P] * y[V.PTP1B]
            - x[C.kr138] * y[V.E24P_Grb2_Gab1P_PTP1B]
        )
        v[139] = (
            x[C.kf139] * y[V.E34P_Grb2_Gab1P] * y[V.PTP1B]
            - x[C.kr139] * y[V.E34P_Grb2_Gab1P_PTP1B]
        )
        v[140] = (
            x[C.kf140] * y[V.E44P_Grb2_Gab1P] * y[V.PTP1B]
            - x[C.kr140] * y[V.E44P_Grb2_Gab1P_PTP1B]
        )
        v[141] = (
            x[C.kf141] * y[V.E11P_ShcP_Grb2_Gab1P] * y[V.PTP1B]
            - x[C.kr141] * y[V.E11P_ShcP_Grb2_Gab1P_PTP1B]
        )
        v[142] = (
            x[C.kf142] * y[V.E12P_ShcP_Grb2_Gab1P] * y[V.PTP1B]
            - x[C.kr142] * y[V.E12P_ShcP_Grb2_Gab1P_PTP1B]
        )
        v[143] = (
            x[C.kf143] * y[V.E23P_ShcP_Grb2_Gab1P] * y[V.PTP1B]
            - x[C.kr143] * y[V.E23P_ShcP_Grb2_Gab1P_PTP1B]
        )
        v[144] = (
            x[C.kf144] * y[V.E24P_ShcP_Grb2_Gab1P] * y[V.PTP1B]
            - x[C.kr144] * y[V.E24P_ShcP_Grb2_Gab1P_PTP1B]
        )
        v[145] = (
            x[C.kf145] * y[V.E34P_ShcP_Grb2_Gab1P] * y[V.PTP1B]
            - x[C.kr145] * y[V.E34P_ShcP_Grb2_Gab1P_PTP1B]
        )
        v[146] = (
            x[C.kf146] * y[V.E44P_ShcP_Grb2_Gab1P] * y[V.PTP1B]
            - x[C.kr146] * y[V.E44P_ShcP_Grb2_Gab1P_PTP1B]
        )
        v[147] = x[C.kf147] * y[V.PIP3_Gab1P] * y[V.PTP1B] - x[C.kr147] * y[V.PIP3_Gab1P_PTP1B]
        v[148] = (
            x[C.kf148] * y[V.E11P_Grb2_Gab1P_PTP1B] - x[C.kr148] * y[V.E11P_Grb2_Gab1] * y[V.PTP1B]
        )
        v[149] = (
            x[C.kf149] * y[V.E12P_Grb2_Gab1P_PTP1B] - x[C.kr149] * y[V.E12P_Grb2_Gab1] * y[V.PTP1B]
        )
        v[150] = (
            x[C.kf150] * y[V.E23P_Grb2_Gab1P_PTP1B] - x[C.kr150] * y[V.E23P_Grb2_Gab1] * y[V.PTP1B]
        )
        v[151] = (
            x[C.kf151] * y[V.E24P_Grb2_Gab1P_PTP1B] - x[C.kr151] * y[V.E24P_Grb2_Gab1] * y[V.PTP1B]
        )
        v[152] = (
            x[C.kf152] * y[V.E34P_Grb2_Gab1P_PTP1B] - x[C.kr152] * y[V.E34P_Grb2_Gab1] * y[V.PTP1B]
        )
        v[153] = (
            x[C.kf153] * y[V.E44P_Grb2_Gab1P_PTP1B] - x[C.kr153] * y[V.E44P_Grb2_Gab1] * y[V.PTP1B]
        )
        v[154] = (
            x[C.kf154] * y[V.E11P_ShcP_Grb2_Gab1P_PTP1B]
            - x[C.kr154] * y[V.E11P_ShcP_Grb2_Gab1] * y[V.PTP1B]
        )
        v[155] = (
            x[C.kf155] * y[V.E12P_ShcP_Grb2_Gab1P_PTP1B]
            - x[C.kr155] * y[V.E12P_ShcP_Grb2_Gab1] * y[V.PTP1B]
        )
        v[156] = (
            x[C.kf156] * y[V.E23P_ShcP_Grb2_Gab1P_PTP1B]
            - x[C.kr156] * y[V.E23P_ShcP_Grb2_Gab1] * y[V.PTP1B]
        )
        v[157] = (
            x[C.kf157] * y[V.E24P_ShcP_Grb2_Gab1P_PTP1B]
            - x[C.kr157] * y[V.E24P_ShcP_Grb2_Gab1] * y[V.PTP1B]
        )
        v[158] = (
            x[C.kf158] * y[V.E34P_ShcP_Grb2_Gab1P_PTP1B]
            - x[C.kr158] * y[V.E34P_ShcP_Grb2_Gab1] * y[V.PTP1B]
        )
        v[159] = (
            x[C.kf159] * y[V.E44P_ShcP_Grb2_Gab1P_PTP1B]
            - x[C.kr159] * y[V.E44P_ShcP_Grb2_Gab1] * y[V.PTP1B]
        )
        v[160] = x[C.kf160] * y[V.PIP3_Gab1P_PTP1B] - x[C.kr160] * y[V.PIP3_Gab1] * y[V.PTP1B]
        v[161] = (
            x[C.kf161] * y[V.E11P_Grb2_Gab1P] * y[V.RasGAP]
            - x[C.kr161] * y[V.E11P_Grb2_Gab1P_RasGAP]
        )
        v[162] = (
            x[C.kf162] * y[V.E12P_Grb2_Gab1P] * y[V.RasGAP]
            - x[C.kr162] * y[V.E12P_Grb2_Gab1P_RasGAP]
        )
        v[163] = (
            x[C.kf163] * y[V.E23P_Grb2_Gab1P] * y[V.RasGAP]
            - x[C.kr163] * y[V.E23P_Grb2_Gab1P_RasGAP]
        )
        v[164] = (
            x[C.kf164] * y[V.E24P_Grb2_Gab1P] * y[V.RasGAP]
            - x[C.kr164] * y[V.E24P_Grb2_Gab1P_RasGAP]
        )
        v[165] = (
            x[C.kf165] * y[V.E34P_Grb2_Gab1P] * y[V.RasGAP]
            - x[C.kr165] * y[V.E34P_Grb2_Gab1P_RasGAP]
        )
        v[166] = (
            x[C.kf166] * y[V.E44P_Grb2_Gab1P] * y[V.RasGAP]
            - x[C.kr166] * y[V.E44P_Grb2_Gab1P_RasGAP]
        )
        v[167] = (
            x[C.kf167] * y[V.E11P_ShcP_Grb2_Gab1P] * y[V.RasGAP]
            - x[C.kr167] * y[V.E11P_ShcP_Grb2_Gab1P_RasGAP]
        )
        v[168] = (
            x[C.kf168] * y[V.E12P_ShcP_Grb2_Gab1P] * y[V.RasGAP]
            - x[C.kr168] * y[V.E12P_ShcP_Grb2_Gab1P_RasGAP]
        )
        v[169] = (
            x[C.kf169] * y[V.E23P_ShcP_Grb2_Gab1P] * y[V.RasGAP]
            - x[C.kr169] * y[V.E23P_ShcP_Grb2_Gab1P_RasGAP]
        )
        v[170] = (
            x[C.kf170] * y[V.E24P_ShcP_Grb2_Gab1P] * y[V.RasGAP]
            - x[C.kr170] * y[V.E24P_ShcP_Grb2_Gab1P_RasGAP]
        )
        v[171] = (
            x[C.kf171] * y[V.E34P_ShcP_Grb2_Gab1P] * y[V.RasGAP]
            - x[C.kr171] * y[V.E34P_ShcP_Grb2_Gab1P_RasGAP]
        )
        v[172] = (
            x[C.kf172] * y[V.E44P_ShcP_Grb2_Gab1P] * y[V.RasGAP]
            - x[C.kr172] * y[V.E44P_ShcP_Grb2_Gab1P_RasGAP]
        )
        v[173] = x[C.kf173] * y[V.PIP3_Gab1P] * y[V.RasGAP] - x[C.kr173] * y[V.PIP3_Gab1P_RasGAP]
        v[174] = x[C.kf174] * y[V.E11P_RasGAP] * y[V.RasGTP] - x[C.kr174] * y[V.E11P_RasGAP_RasGTP]
        v[175] = x[C.kf175] * y[V.E12P_RasGAP] * y[V.RasGTP] - x[C.kr175] * y[V.E12P_RasGAP_RasGTP]
        v[176] = x[C.kf176] * y[V.E23P_RasGAP] * y[V.RasGTP] - x[C.kr176] * y[V.E23P_RasGAP_RasGTP]
        v[177] = x[C.kf177] * y[V.E24P_RasGAP] * y[V.RasGTP] - x[C.kr177] * y[V.E24P_RasGAP_RasGTP]
        v[178] = x[C.kf178] * y[V.E34P_RasGAP] * y[V.RasGTP] - x[C.kr178] * y[V.E34P_RasGAP_RasGTP]
        v[179] = x[C.kf179] * y[V.E44P_RasGAP] * y[V.RasGTP] - x[C.kr179] * y[V.E44P_RasGAP_RasGTP]
        v[180] = (
            x[C.kf180] * y[V.E11P_Grb2_Gab1P_RasGAP] * y[V.RasGTP]
            - x[C.kr180] * y[V.E11P_Grb2_Gab1P_RasGAP_RasGTP]
        )
        v[181] = (
            x[C.kf181] * y[V.E12P_Grb2_Gab1P_RasGAP] * y[V.RasGTP]
            - x[C.kr181] * y[V.E12P_Grb2_Gab1P_RasGAP_RasGTP]
        )
        v[182] = (
            x[C.kf182] * y[V.E23P_Grb2_Gab1P_RasGAP] * y[V.RasGTP]
            - x[C.kr182] * y[V.E23P_Grb2_Gab1P_RasGAP_RasGTP]
        )
        v[183] = (
            x[C.kf183] * y[V.E24P_Grb2_Gab1P_RasGAP] * y[V.RasGTP]
            - x[C.kr183] * y[V.E24P_Grb2_Gab1P_RasGAP_RasGTP]
        )
        v[184] = (
            x[C.kf184] * y[V.E34P_Grb2_Gab1P_RasGAP] * y[V.RasGTP]
            - x[C.kr184] * y[V.E34P_Grb2_Gab1P_RasGAP_RasGTP]
        )
        v[185] = (
            x[C.kf185] * y[V.E44P_Grb2_Gab1P_RasGAP] * y[V.RasGTP]
            - x[C.kr185] * y[V.E44P_Grb2_Gab1P_RasGAP_RasGTP]
        )
        v[186] = (
            x[C.kf186] * y[V.E11P_ShcP_Grb2_Gab1P_RasGAP] * y[V.RasGTP]
            - x[C.kr186] * y[V.E11P_ShcP_Grb2_Gab1P_RasGAP_RasGTP]
        )
        v[187] = (
            x[C.kf187] * y[V.E12P_ShcP_Grb2_Gab1P_RasGAP] * y[V.RasGTP]
            - x[C.kr187] * y[V.E12P_ShcP_Grb2_Gab1P_RasGAP_RasGTP]
        )
        v[188] = (
            x[C.kf188] * y[V.E23P_ShcP_Grb2_Gab1P_RasGAP] * y[V.RasGTP]
            - x[C.kr188] * y[V.E23P_ShcP_Grb2_Gab1P_RasGAP_RasGTP]
        )
        v[189] = (
            x[C.kf189] * y[V.E24P_ShcP_Grb2_Gab1P_RasGAP] * y[V.RasGTP]
            - x[C.kr189] * y[V.E24P_ShcP_Grb2_Gab1P_RasGAP_RasGTP]
        )
        v[190] = (
            x[C.kf190] * y[V.E34P_ShcP_Grb2_Gab1P_RasGAP] * y[V.RasGTP]
            - x[C.kr190] * y[V.E34P_ShcP_Grb2_Gab1P_RasGAP_RasGTP]
        )
        v[191] = (
            x[C.kf191] * y[V.E44P_ShcP_Grb2_Gab1P_RasGAP] * y[V.RasGTP]
            - x[C.kr191] * y[V.E44P_ShcP_Grb2_Gab1P_RasGAP_RasGTP]
        )
        v[192] = (
            x[C.kf192] * y[V.PIP3_Gab1P_RasGAP] * y[V.RasGTP]
            - x[C.kr192] * y[V.PIP3_Gab1P_RasGAP_RasGTP]
        )
        v[193] = x[C.kf193] * y[V.E11P_RasGAP_RasGTP] - x[C.kr193] * y[V.E11P_RasGAP] * y[V.RasGDP]
        v[194] = x[C.kf194] * y[V.E12P_RasGAP_RasGTP] - x[C.kr194] * y[V.E12P_RasGAP] * y[V.RasGDP]
        v[195] = x[C.kf195] * y[V.E23P_RasGAP_RasGTP] - x[C.kr195] * y[V.E23P_RasGAP] * y[V.RasGDP]
        v[196] = x[C.kf196] * y[V.E24P_RasGAP_RasGTP] - x[C.kr196] * y[V.E24P_RasGAP] * y[V.RasGDP]
        v[197] = x[C.kf197] * y[V.E34P_RasGAP_RasGTP] - x[C.kr197] * y[V.E34P_RasGAP] * y[V.RasGDP]
        v[198] = x[C.kf198] * y[V.E44P_RasGAP_RasGTP] - x[C.kr198] * y[V.E44P_RasGAP] * y[V.RasGDP]
        v[199] = (
            x[C.kf199] * y[V.E11P_Grb2_Gab1P_RasGAP_RasGTP]
            - x[C.kr199] * y[V.E11P_Grb2_Gab1P_RasGAP] * y[V.RasGDP]
        )
        v[200] = (
            x[C.kf200] * y[V.E12P_Grb2_Gab1P_RasGAP_RasGTP]
            - x[C.kr200] * y[V.E12P_Grb2_Gab1P_RasGAP] * y[V.RasGDP]
        )
        v[201] = (
            x[C.kf201] * y[V.E23P_Grb2_Gab1P_RasGAP_RasGTP]
            - x[C.kr201] * y[V.E23P_Grb2_Gab1P_RasGAP] * y[V.RasGDP]
        )
        v[202] = (
            x[C.kf202] * y[V.E24P_Grb2_Gab1P_RasGAP_RasGTP]
            - x[C.kr202] * y[V.E24P_Grb2_Gab1P_RasGAP] * y[V.RasGDP]
        )
        v[203] = (
            x[C.kf203] * y[V.E34P_Grb2_Gab1P_RasGAP_RasGTP]
            - x[C.kr203] * y[V.E34P_Grb2_Gab1P_RasGAP] * y[V.RasGDP]
        )
        v[204] = (
            x[C.kf204] * y[V.E44P_Grb2_Gab1P_RasGAP_RasGTP]
            - x[C.kr204] * y[V.E44P_Grb2_Gab1P_RasGAP] * y[V.RasGDP]
        )
        v[205] = (
            x[C.kf205] * y[V.E11P_ShcP_Grb2_Gab1P_RasGAP_RasGTP]
            - x[C.kr205] * y[V.E11P_ShcP_Grb2_Gab1P_RasGAP] * y[V.RasGDP]
        )
        v[206] = (
            x[C.kf206] * y[V.E12P_ShcP_Grb2_Gab1P_RasGAP_RasGTP]
            - x[C.kr206] * y[V.E12P_ShcP_Grb2_Gab1P_RasGAP] * y[V.RasGDP]
        )
        v[207] = (
            x[C.kf207] * y[V.E23P_ShcP_Grb2_Gab1P_RasGAP_RasGTP]
            - x[C.kr207] * y[V.E23P_ShcP_Grb2_Gab1P_RasGAP] * y[V.RasGDP]
        )
        v[208] = (
            x[C.kf208] * y[V.E24P_ShcP_Grb2_Gab1P_RasGAP_RasGTP]
            - x[C.kr208] * y[V.E24P_ShcP_Grb2_Gab1P_RasGAP] * y[V.RasGDP]
        )
        v[209] = (
            x[C.kf209] * y[V.E34P_ShcP_Grb2_Gab1P_RasGAP_RasGTP]
            - x[C.kr209] * y[V.E34P_ShcP_Grb2_Gab1P_RasGAP] * y[V.RasGDP]
        )
        v[210] = (
            x[C.kf210] * y[V.E44P_ShcP_Grb2_Gab1P_RasGAP_RasGTP]
            - x[C.kr210] * y[V.E44P_ShcP_Grb2_Gab1P_RasGAP] * y[V.RasGDP]
        )
        v[211] = (
            x[C.kf211] * y[V.PIP3_Gab1P_RasGAP_RasGTP]
            - x[C.kr211] * y[V.PIP3_Gab1P_RasGAP] * y[V.RasGDP]
        )
        v[212] = (
            x[C.kf212] * y[V.E11P_Grb2_Gab1P] * y[V.PI3K] - x[C.kr212] * y[V.E11P_Grb2_Gab1P_PI3K]
        )
        v[213] = (
            x[C.kf213] * y[V.E12P_Grb2_Gab1P] * y[V.PI3K] - x[C.kr213] * y[V.E12P_Grb2_Gab1P_PI3K]
        )
        v[214] = (
            x[C.kf214] * y[V.E23P_Grb2_Gab1P] * y[V.PI3K] - x[C.kr214] * y[V.E23P_Grb2_Gab1P_PI3K]
        )
        v[215] = (
            x[C.kf215] * y[V.E24P_Grb2_Gab1P] * y[V.PI3K] - x[C.kr215] * y[V.E24P_Grb2_Gab1P_PI3K]
        )
        v[216] = (
            x[C.kf216] * y[V.E34P_Grb2_Gab1P] * y[V.PI3K] - x[C.kr216] * y[V.E34P_Grb2_Gab1P_PI3K]
        )
        v[217] = (
            x[C.kf217] * y[V.E44P_Grb2_Gab1P] * y[V.PI3K] - x[C.kr217] * y[V.E44P_Grb2_Gab1P_PI3K]
        )
        v[218] = (
            x[C.kf218] * y[V.E11P_ShcP_Grb2_Gab1P] * y[V.PI3K]
            - x[C.kr218] * y[V.E11P_ShcP_Grb2_Gab1P_PI3K]
        )
        v[219] = (
            x[C.kf219] * y[V.E12P_ShcP_Grb2_Gab1P] * y[V.PI3K]
            - x[C.kr219] * y[V.E12P_ShcP_Grb2_Gab1P_PI3K]
        )
        v[220] = (
            x[C.kf220] * y[V.E23P_ShcP_Grb2_Gab1P] * y[V.PI3K]
            - x[C.kr220] * y[V.E23P_ShcP_Grb2_Gab1P_PI3K]
        )
        v[221] = (
            x[C.kf221] * y[V.E24P_ShcP_Grb2_Gab1P] * y[V.PI3K]
            - x[C.kr221] * y[V.E24P_ShcP_Grb2_Gab1P_PI3K]
        )
        v[222] = (
            x[C.kf222] * y[V.E34P_ShcP_Grb2_Gab1P] * y[V.PI3K]
            - x[C.kr222] * y[V.E34P_ShcP_Grb2_Gab1P_PI3K]
        )
        v[223] = (
            x[C.kf223] * y[V.E44P_ShcP_Grb2_Gab1P] * y[V.PI3K]
            - x[C.kr223] * y[V.E44P_ShcP_Grb2_Gab1P_PI3K]
        )
        v[224] = x[C.kf224] * y[V.PIP3_Gab1P] * y[V.PI3K] - x[C.kr224] * y[V.PIP3_Gab1P_PI3K]
        v[225] = x[C.kf225] * y[V.E23P_PI3K] * y[V.PIP2] - x[C.kr225] * y[V.E23P_PI3K_PIP2]
        v[226] = x[C.kf226] * y[V.E24P_PI3K] * y[V.PIP2] - x[C.kr226] * y[V.E24P_PI3K_PIP2]
        v[227] = x[C.kf227] * y[V.E34P_PI3K] * y[V.PIP2] - x[C.kr227] * y[V.E34P_PI3K_PIP2]
        v[228] = x[C.kf228] * y[V.E44P_PI3K] * y[V.PIP2] - x[C.kr228] * y[V.E44P_PI3K_PIP2]
        v[229] = (
            x[C.kf229] * y[V.E11P_Grb2_Gab1P_PI3K] * y[V.PIP2]
            - x[C.kr229] * y[V.E11P_Grb2_Gab1P_PI3K_PIP2]
        )
        v[230] = (
            x[C.kf230] * y[V.E12P_Grb2_Gab1P_PI3K] * y[V.PIP2]
            - x[C.kr230] * y[V.E12P_Grb2_Gab1P_PI3K_PIP2]
        )
        v[231] = (
            x[C.kf231] * y[V.E23P_Grb2_Gab1P_PI3K] * y[V.PIP2]
            - x[C.kr231] * y[V.E23P_Grb2_Gab1P_PI3K_PIP2]
        )
        v[232] = (
            x[C.kf232] * y[V.E24P_Grb2_Gab1P_PI3K] * y[V.PIP2]
            - x[C.kr232] * y[V.E24P_Grb2_Gab1P_PI3K_PIP2]
        )
        v[233] = (
            x[C.kf233] * y[V.E34P_Grb2_Gab1P_PI3K] * y[V.PIP2]
            - x[C.kr233] * y[V.E34P_Grb2_Gab1P_PI3K_PIP2]
        )
        v[234] = (
            x[C.kf234] * y[V.E44P_Grb2_Gab1P_PI3K] * y[V.PIP2]
            - x[C.kr234] * y[V.E44P_Grb2_Gab1P_PI3K_PIP2]
        )
        v[235] = (
            x[C.kf235] * y[V.E11P_ShcP_Grb2_Gab1P_PI3K] * y[V.PIP2]
            - x[C.kr235] * y[V.E11P_ShcP_Grb2_Gab1P_PI3K_PIP2]
        )
        v[236] = (
            x[C.kf236] * y[V.E12P_ShcP_Grb2_Gab1P_PI3K] * y[V.PIP2]
            - x[C.kr236] * y[V.E12P_ShcP_Grb2_Gab1P_PI3K_PIP2]
        )
        v[237] = (
            x[C.kf237] * y[V.E23P_ShcP_Grb2_Gab1P_PI3K] * y[V.PIP2]
            - x[C.kr237] * y[V.E23P_ShcP_Grb2_Gab1P_PI3K_PIP2]
        )
        v[238] = (
            x[C.kf238] * y[V.E24P_ShcP_Grb2_Gab1P_PI3K] * y[V.PIP2]
            - x[C.kr238] * y[V.E24P_ShcP_Grb2_Gab1P_PI3K_PIP2]
        )
        v[239] = (
            x[C.kf239] * y[V.E34P_ShcP_Grb2_Gab1P_PI3K] * y[V.PIP2]
            - x[C.kr239] * y[V.E34P_ShcP_Grb2_Gab1P_PI3K_PIP2]
        )
        v[240] = (
            x[C.kf240] * y[V.E44P_ShcP_Grb2_Gab1P_PI3K] * y[V.PIP2]
            - x[C.kr240] * y[V.E44P_ShcP_Grb2_Gab1P_PI3K_PIP2]
        )
        v[241] = (
            x[C.kf241] * y[V.PIP3_Gab1P_PI3K] * y[V.PIP2] - x[C.kr241] * y[V.PIP3_Gab1P_PI3K_PIP2]
        )
        v[242] = x[C.kf242] * y[V.E23P_PI3K_PIP2] - x[C.kr242] * y[V.E23P_PI3K] * y[V.PIP3]
        v[243] = x[C.kf243] * y[V.E24P_PI3K_PIP2] - x[C.kr243] * y[V.E24P_PI3K] * y[V.PIP3]
        v[244] = x[C.kf244] * y[V.E34P_PI3K_PIP2] - x[C.kr244] * y[V.E34P_PI3K] * y[V.PIP3]
        v[245] = x[C.kf245] * y[V.E44P_PI3K_PIP2] - x[C.kr245] * y[V.E44P_PI3K] * y[V.PIP3]
        v[246] = (
            x[C.kf246] * y[V.E11P_Grb2_Gab1P_PI3K_PIP2]
            - x[C.kr246] * y[V.E11P_Grb2_Gab1P_PI3K] * y[V.PIP3]
        )
        v[247] = (
            x[C.kf247] * y[V.E12P_Grb2_Gab1P_PI3K_PIP2]
            - x[C.kr247] * y[V.E12P_Grb2_Gab1P_PI3K] * y[V.PIP3]
        )
        v[248] = (
            x[C.kf248] * y[V.E23P_Grb2_Gab1P_PI3K_PIP2]
            - x[C.kr248] * y[V.E23P_Grb2_Gab1P_PI3K] * y[V.PIP3]
        )
        v[249] = (
            x[C.kf249] * y[V.E24P_Grb2_Gab1P_PI3K_PIP2]
            - x[C.kr249] * y[V.E24P_Grb2_Gab1P_PI3K] * y[V.PIP3]
        )
        v[250] = (
            x[C.kf250] * y[V.E34P_Grb2_Gab1P_PI3K_PIP2]
            - x[C.kr250] * y[V.E34P_Grb2_Gab1P_PI3K] * y[V.PIP3]
        )
        v[251] = (
            x[C.kf251] * y[V.E44P_Grb2_Gab1P_PI3K_PIP2]
            - x[C.kr251] * y[V.E44P_Grb2_Gab1P_PI3K] * y[V.PIP3]
        )
        v[252] = (
            x[C.kf252] * y[V.E11P_ShcP_Grb2_Gab1P_PI3K_PIP2]
            - x[C.kr252] * y[V.E11P_ShcP_Grb2_Gab1P_PI3K] * y[V.PIP3]
        )
        v[253] = (
            x[C.kf253] * y[V.E12P_ShcP_Grb2_Gab1P_PI3K_PIP2]
            - x[C.kr253] * y[V.E12P_ShcP_Grb2_Gab1P_PI3K] * y[V.PIP3]
        )
        v[254] = (
            x[C.kf254] * y[V.E23P_ShcP_Grb2_Gab1P_PI3K_PIP2]
            - x[C.kr254] * y[V.E23P_ShcP_Grb2_Gab1P_PI3K] * y[V.PIP3]
        )
        v[255] = (
            x[C.kf255] * y[V.E24P_ShcP_Grb2_Gab1P_PI3K_PIP2]
            - x[C.kr255] * y[V.E24P_ShcP_Grb2_Gab1P_PI3K] * y[V.PIP3]
        )
        v[256] = (
            x[C.kf256] * y[V.E34P_ShcP_Grb2_Gab1P_PI3K_PIP2]
            - x[C.kr256] * y[V.E34P_ShcP_Grb2_Gab1P_PI3K] * y[V.PIP3]
        )
        v[257] = (
            x[C.kf257] * y[V.E44P_ShcP_Grb2_Gab1P_PI3K_PIP2]
            - x[C.kr257] * y[V.E44P_ShcP_Grb2_Gab1P_PI3K] * y[V.PIP3]
        )
        v[258] = (
            x[C.kf258] * y[V.PIP3_Gab1P_PI3K_PIP2] - x[C.kr258] * y[V.PIP3_Gab1P_PI3K] * y[V.PIP3]
        )
        v[259] = x[C.kf259] * y[V.PIP3] * y[V.Gab1] - x[C.kr259] * y[V.PIP3_Gab1]
        v[260] = x[C.V260] * y[V.PTEN] * y[V.PIP3] / (x[C.K260] + y[V.PIP3])
        v[261] = x[C.V261] * y[V.PIP3] * y[V.Akt] / (x[C.K261] + y[V.Akt])
        v[262] = x[C.V262] * y[V.AktP] / (x[C.K262] + y[V.AktP])
        v[263] = x[C.V263] * y[V.AktP] * y[V.Raf] / (x[C.K263] + y[V.Raf])
        v[264] = x[C.V264] * y[V.AktP] * y[V.RafP] / (x[C.K264] + y[V.RafP])
        v[265] = x[C.V265] * y[V.RafP_inactive] / (x[C.K265] + y[V.RafP_inactive])
        v[266] = x[C.V266] * y[V.RafPP_inactive] / (x[C.K266] + y[V.RafPP_inactive])
        v[267] = x[C.kf267] * y[V.AktP]
        v[268] = x[C.kf268] * y[V.E11P]
        v[269] = x[C.kf269] * y[V.E11P_Grb2]
        v[270] = x[C.kf270] * y[V.E11P_Shc]
        v[271] = x[C.kf271] * y[V.E11P_RasGAP]
        v[272] = x[C.kf272] * y[V.E11P_Shc_Grb2]
        v[273] = x[C.kf273] * y[V.E11P_ShcP_Grb2]
        v[274] = x[C.kf274] * y[V.E11P_ShcP_Grb2_PTP1B]
        v[275] = x[C.kf275] * y[V.E11P_Grb2_SOS]
        v[276] = x[C.kf276] * y[V.E11P_ShcP_Grb2_SOS]
        v[277] = x[C.kf277] * y[V.E11P_Grb2_Gab1]
        v[278] = x[C.kf278] * y[V.E11P_Grb2_SOS_RasGDP]
        v[279] = x[C.kf279] * y[V.E11P_ShcP_Grb2_SOS_RasGDP]
        v[280] = x[C.kf280] * y[V.E11P_Grb2_Gab1P]
        v[281] = x[C.kf281] * y[V.E11P_ShcP_Grb2_Gab1P]
        v[282] = x[C.kf282] * y[V.E11P_Grb2_Gab1P_PTP1B]
        v[283] = x[C.kf283] * y[V.E11P_ShcP_Grb2_Gab1P_PTP1B]
        v[284] = x[C.kf284] * y[V.E11P_Grb2_Gab1P_RasGAP]
        v[285] = x[C.kf285] * y[V.E11P_ShcP_Grb2_Gab1P_RasGAP]
        v[286] = x[C.kf286] * y[V.E11P_RasGAP_RasGTP]
        v[287] = x[C.kf287] * y[V.E11P_Grb2_Gab1P_RasGAP_RasGTP]
        v[288] = x[C.kf288] * y[V.E11P_ShcP_Grb2_Gab1P_RasGAP_RasGTP]
        v[289] = x[C.kf289] * y[V.E11P_Grb2_Gab1P_PI3K]
        v[290] = x[C.kf290] * y[V.E11P_ShcP_Grb2_Gab1P_PI3K]
        v[291] = (
            x[C.V291]
            * y[V.ERKPP] ** x[C.n291]
            / (x[C.K291] ** x[C.n291] + y[V.ERKPP] ** x[C.n291])
        )
        v[292] = x[C.kf292] * y[V.PreduspmRNA] - x[C.kr292] * (0.94 / 0.22) * y[V.duspmRNA]
        v[293] = x[C.kf293] * y[V.duspmRNA]
        v[294] = x[C.kf294] * y[V.duspmRNA]
        v[295] = x[C.V295] * y[V.ERKPP] * y[V.DUSP] / (x[C.K295] + y[V.DUSP])
        v[296] = x[C.kf296] * y[V.DUSP]
        v[297] = x[C.kf297] * y[V.DUSPP]
        v[298] = x[C.kf298] * y[V.ERKPP] * y[V.DUSP] - x[C.kr298] * y[V.DUSP_ERKPP]
        v[299] = x[C.kf299] * y[V.DUSP_ERKPP] - x[C.kr299] * y[V.ERKP] * y[V.DUSP]
        v[300] = x[C.kf300] * y[V.ERKP] * y[V.DUSP] - x[C.kr300] * y[V.DUSP_ERKP]
        v[301] = x[C.kf301] * y[V.DUSP_ERKP] - x[C.kr301] * y[V.ERK] * y[V.DUSP]
        v[302] = x[C.kf302] * y[V.ERK] * y[V.DUSP] - x[C.kr302] * y[V.DUSP_ERK]
        v[303] = x[C.kf303] * y[V.ERKPP] * y[V.DUSPP] - x[C.kr303] * y[V.DUSPP_ERKPP]
        v[304] = x[C.kf304] * y[V.DUSPP_ERKPP] - x[C.kr304] * y[V.ERKP] * y[V.DUSPP]
        v[305] = x[C.kf305] * y[V.ERKP] * y[V.DUSPP] - x[C.kr305] * y[V.DUSPP_ERKP]
        v[306] = x[C.kf306] * y[V.DUSPP_ERKP] - x[C.kr306] * y[V.ERK] * y[V.DUSPP]
        v[307] = x[C.kf307] * y[V.ERK] * y[V.DUSPP] - x[C.kr307] * y[V.DUSPP_ERK]
        v[308] = x[C.V308] * y[V.AktP] * y[V.GSK3b] / (x[C.K308] + y[V.GSK3b])
        v[309] = x[C.V309] * y[V.GSK3bP] / (x[C.K309] + y[V.GSK3bP])
        v[310] = (
            x[C.V310]
            * y[V.ERKPP] ** x[C.n310]
            / (
                x[C.K310] ** x[C.n310]
                + y[V.ERKPP] ** x[C.n310]
                + (y[V.cMycP] / x[C.KF310]) ** x[C.nF310]
            )
        )
        v[311] = x[C.kf311] * y[V.PrecmycmRNA] - x[C.kr311] * (0.94 / 0.22) * y[V.cmycmRNA]
        v[312] = x[C.kf312] * y[V.cmycmRNA]
        v[313] = x[C.kf313] * y[V.cmycmRNA]
        v[314] = x[C.V314] * y[V.ERKPP] * y[V.cMyc] / (x[C.K314] + y[V.cMyc])
        v[315] = x[C.V315] * y[V.cMycP] / (x[C.K315] + y[V.cMycP])
        v[316] = x[C.V316] * y[V.GSK3b] * y[V.cMycP] / (x[C.K316] + y[V.cMycP])
        v[317] = x[C.V317] * y[V.cMycPP] / (x[C.K317] + y[V.cMycPP])
        v[318] = x[C.kf318] * y[V.cMyc]
        v[319] = x[C.kf319] * y[V.cMycPP]

        if self.perturbation:
            for i, dv in self.perturbation.items():
                v[i] = v[i] * dv

        dydt = [0] * V.NUM
        dydt[V.EGF] = -v[1]
        dydt[V.ErbB1] = -v[1]
        dydt[V.EGF_ErbB1] = +v[1] - 2 * v[4] - v[5]
        dydt[V.HRG] = -v[2] - v[3]
        dydt[V.ErbB3] = -v[2]
        dydt[V.HRG_ErbB3] = +v[2] - v[6] - v[7]
        dydt[V.ErbB4] = -v[3]
        dydt[V.HRG_ErbB4] = +v[3] - v[7] - v[8] - 2 * v[9]
        dydt[V.E11] = +v[4] - v[10]
        dydt[V.ErbB2] = -v[5] - v[6] - v[8]
        dydt[V.E12] = +v[5] - v[11]
        dydt[V.E23] = +v[6] - v[12]
        dydt[V.E34] = +v[7] - v[14]
        dydt[V.E24] = +v[8] - v[13]
        dydt[V.E44] = +v[9] - v[15]
        dydt[V.E11P] = +v[10] - v[16] - v[17] - v[18] - v[268]
        dydt[V.E12P] = +v[11] - v[19] - v[20] - v[21]
        dydt[V.E23P] = +v[12] - v[22] - v[23] - v[24] - v[25]
        dydt[V.E24P] = +v[13] - v[26] - v[27] - v[28] - v[29]
        dydt[V.E34P] = +v[14] - v[30] - v[31] - v[32] - v[33]
        dydt[V.E44P] = +v[15] - v[34] - v[35] - v[36] - v[37]
        dydt[V.Grb2] = (
            -v[16]
            - v[19]
            - v[22]
            - v[26]
            - v[30]
            - v[34]
            - v[38]
            - v[39]
            - v[40]
            - v[41]
            - v[42]
            - v[43]
        )
        dydt[V.E11P_Grb2] = +v[16] - v[62] - v[74] - v[269]
        dydt[V.Shc] = -v[17] - v[20] - v[23] - v[27] - v[31] - v[35]
        dydt[V.E11P_Shc] = +v[17] - v[38] - v[270]
        dydt[V.RasGAP] = (
            -v[18]
            - v[21]
            - v[25]
            - v[29]
            - v[33]
            - v[37]
            - v[161]
            - v[162]
            - v[163]
            - v[164]
            - v[165]
            - v[166]
            - v[167]
            - v[168]
            - v[169]
            - v[170]
            - v[171]
            - v[172]
            - v[173]
        )
        dydt[V.E11P_RasGAP] = +v[18] - v[174] + v[193] - v[271]
        dydt[V.E12P_Grb2] = +v[19] - v[63] - v[75]
        dydt[V.E12P_Shc] = +v[20] - v[39]
        dydt[V.E12P_RasGAP] = +v[21] - v[175] + v[194]
        dydt[V.E23P_Grb2] = +v[22] - v[64] - v[76]
        dydt[V.E23P_Shc] = +v[23] - v[40]
        dydt[V.PI3K] = (
            -v[24]
            - v[28]
            - v[32]
            - v[36]
            - v[212]
            - v[213]
            - v[214]
            - v[215]
            - v[216]
            - v[217]
            - v[218]
            - v[219]
            - v[220]
            - v[221]
            - v[222]
            - v[223]
            - v[224]
        )
        dydt[V.E23P_PI3K] = +v[24] - v[225] + v[242]
        dydt[V.E23P_RasGAP] = +v[25] - v[176] + v[195]
        dydt[V.E24P_Grb2] = +v[26] - v[65] - v[77]
        dydt[V.E24P_Shc] = +v[27] - v[41]
        dydt[V.E24P_PI3K] = +v[28] - v[226] + v[243]
        dydt[V.E24P_RasGAP] = +v[29] - v[177] + v[196]
        dydt[V.E34P_Grb2] = +v[30] - v[66] - v[78]
        dydt[V.E34P_Shc] = +v[31] - v[42]
        dydt[V.E34P_PI3K] = +v[32] - v[227] + v[244]
        dydt[V.E34P_RasGAP] = +v[33] - v[178] + v[197]
        dydt[V.E44P_Grb2] = +v[34] - v[67] - v[79]
        dydt[V.E44P_Shc] = +v[35] - v[43]
        dydt[V.E44P_PI3K] = +v[36] - v[228] + v[245]
        dydt[V.E44P_RasGAP] = +v[37] - v[179] + v[198]
        dydt[V.E11P_Shc_Grb2] = +v[38] - v[44] + v[56] - v[272]
        dydt[V.E12P_Shc_Grb2] = +v[39] - v[45] + v[57]
        dydt[V.E23P_Shc_Grb2] = +v[40] - v[46] + v[58]
        dydt[V.E24P_Shc_Grb2] = +v[41] - v[47] + v[59]
        dydt[V.E34P_Shc_Grb2] = +v[42] - v[48] + v[60]
        dydt[V.E44P_Shc_Grb2] = +v[43] - v[49] + v[61]
        dydt[V.E11P_ShcP_Grb2] = +v[44] - v[50] - v[68] - v[80] - v[273]
        dydt[V.E12P_ShcP_Grb2] = +v[45] - v[51] - v[69] - v[81]
        dydt[V.E23P_ShcP_Grb2] = +v[46] - v[52] - v[70] - v[82]
        dydt[V.E24P_ShcP_Grb2] = +v[47] - v[53] - v[71] - v[83]
        dydt[V.E34P_ShcP_Grb2] = +v[48] - v[54] - v[72] - v[84]
        dydt[V.E44P_ShcP_Grb2] = +v[49] - v[55] - v[73] - v[85]
        dydt[V.PTP1B] = (
            -v[50]
            - v[51]
            - v[52]
            - v[53]
            - v[54]
            - v[55]
            + v[56]
            + v[57]
            + v[58]
            + v[59]
            + v[60]
            + v[61]
            - v[135]
            - v[136]
            - v[137]
            - v[138]
            - v[139]
            - v[140]
            - v[141]
            - v[142]
            - v[143]
            - v[144]
            - v[145]
            - v[146]
            - v[147]
            + v[148]
            + v[149]
            + v[150]
            + v[151]
            + v[152]
            + v[153]
            + v[154]
            + v[155]
            + v[156]
            + v[157]
            + v[158]
            + v[159]
            + v[160]
        )
        dydt[V.E11P_ShcP_Grb2_PTP1B] = +v[50] - v[56] - v[274]
        dydt[V.E12P_ShcP_Grb2_PTP1B] = +v[51] - v[57]
        dydt[V.E23P_ShcP_Grb2_PTP1B] = +v[52] - v[58]
        dydt[V.E24P_ShcP_Grb2_PTP1B] = +v[53] - v[59]
        dydt[V.E34P_ShcP_Grb2_PTP1B] = +v[54] - v[60]
        dydt[V.E44P_ShcP_Grb2_PTP1B] = +v[55] - v[61]
        dydt[V.SOS] = (
            -v[62]
            - v[63]
            - v[64]
            - v[65]
            - v[66]
            - v[67]
            - v[68]
            - v[69]
            - v[70]
            - v[71]
            - v[72]
            - v[73]
            - v[118]
            + v[119]
        )
        dydt[V.E11P_Grb2_SOS] = +v[62] - v[86] + v[98] - v[275]
        dydt[V.E12P_Grb2_SOS] = +v[63] - v[87] + v[99]
        dydt[V.E23P_Grb2_SOS] = +v[64] - v[88] + v[100]
        dydt[V.E24P_Grb2_SOS] = +v[65] - v[89] + v[101]
        dydt[V.E34P_Grb2_SOS] = +v[66] - v[90] + v[102]
        dydt[V.E44P_Grb2_SOS] = +v[67] - v[91] + v[103]
        dydt[V.E11P_ShcP_Grb2_SOS] = +v[68] - v[92] + v[104] - v[276]
        dydt[V.E12P_ShcP_Grb2_SOS] = +v[69] - v[93] + v[105]
        dydt[V.E23P_ShcP_Grb2_SOS] = +v[70] - v[94] + v[106]
        dydt[V.E24P_ShcP_Grb2_SOS] = +v[71] - v[95] + v[107]
        dydt[V.E34P_ShcP_Grb2_SOS] = +v[72] - v[96] + v[108]
        dydt[V.E44P_ShcP_Grb2_SOS] = +v[73] - v[97] + v[109]
        dydt[V.Gab1] = (
            -v[74]
            - v[75]
            - v[76]
            - v[77]
            - v[78]
            - v[79]
            - v[80]
            - v[81]
            - v[82]
            - v[83]
            - v[84]
            - v[85]
            - v[120]
            + v[121]
            - v[259]
        )
        dydt[V.E11P_Grb2_Gab1] = +v[74] - v[122] + v[148] - v[277]
        dydt[V.E12P_Grb2_Gab1] = +v[75] - v[123] + v[149]
        dydt[V.E23P_Grb2_Gab1] = +v[76] - v[124] + v[150]
        dydt[V.E24P_Grb2_Gab1] = +v[77] - v[125] + v[151]
        dydt[V.E34P_Grb2_Gab1] = +v[78] - v[126] + v[152]
        dydt[V.E44P_Grb2_Gab1] = +v[79] - v[127] + v[153]
        dydt[V.E11P_ShcP_Grb2_Gab1] = +v[80] - v[128] + v[154]
        dydt[V.E12P_ShcP_Grb2_Gab1] = +v[81] - v[129] + v[155]
        dydt[V.E23P_ShcP_Grb2_Gab1] = +v[82] - v[130] + v[156]
        dydt[V.E24P_ShcP_Grb2_Gab1] = +v[83] - v[131] + v[157]
        dydt[V.E34P_ShcP_Grb2_Gab1] = +v[84] - v[132] + v[158]
        dydt[V.E44P_ShcP_Grb2_Gab1] = +v[85] - v[133] + v[159]
        dydt[V.RasGDP] = (
            -v[86]
            - v[87]
            - v[88]
            - v[89]
            - v[90]
            - v[91]
            - v[92]
            - v[93]
            - v[94]
            - v[95]
            - v[96]
            - v[97]
            + v[193]
            + v[194]
            + v[195]
            + v[196]
            + v[197]
            + v[198]
            + v[199]
            + v[200]
            + v[201]
            + v[202]
            + v[203]
            + v[204]
            + v[205]
            + v[206]
            + v[207]
            + v[208]
            + v[209]
            + v[210]
            + v[211]
        )
        dydt[V.E11P_Grb2_SOS_RasGDP] = +v[86] - v[98] - v[278]
        dydt[V.E12P_Grb2_SOS_RasGDP] = +v[87] - v[99]
        dydt[V.E23P_Grb2_SOS_RasGDP] = +v[88] - v[100]
        dydt[V.E24P_Grb2_SOS_RasGDP] = +v[89] - v[101]
        dydt[V.E34P_Grb2_SOS_RasGDP] = +v[90] - v[102]
        dydt[V.E44P_Grb2_SOS_RasGDP] = +v[91] - v[103]
        dydt[V.E11P_ShcP_Grb2_SOS_RasGDP] = +v[92] - v[104] - v[279]
        dydt[V.E12P_ShcP_Grb2_SOS_RasGDP] = +v[93] - v[105]
        dydt[V.E23P_ShcP_Grb2_SOS_RasGDP] = +v[94] - v[106]
        dydt[V.E24P_ShcP_Grb2_SOS_RasGDP] = +v[95] - v[107]
        dydt[V.E34P_ShcP_Grb2_SOS_RasGDP] = +v[96] - v[108]
        dydt[V.E44P_ShcP_Grb2_SOS_RasGDP] = +v[97] - v[109]
        dydt[V.RasGTP] = (
            +v[98]
            + v[99]
            + v[100]
            + v[101]
            + v[102]
            + v[103]
            + v[104]
            + v[105]
            + v[106]
            + v[107]
            + v[108]
            + v[109]
            - v[174]
            - v[175]
            - v[176]
            - v[177]
            - v[178]
            - v[179]
            - v[180]
            - v[181]
            - v[182]
            - v[183]
            - v[184]
            - v[185]
            - v[186]
            - v[187]
            - v[188]
            - v[189]
            - v[190]
            - v[191]
            - v[192]
        )
        dydt[V.Raf] = -v[110] + v[111] - v[263] + v[265]
        dydt[V.RafP] = +v[110] - v[111] - v[264] + v[266]
        dydt[V.MEK] = -v[112] + v[113]
        dydt[V.MEKP] = +v[112] - v[113]
        dydt[V.ERK] = -v[114] + v[116] + v[301] - v[302] + v[306] - v[307]
        dydt[V.ERKP] = +v[114] - v[115] - v[116] + v[117] + v[299] - v[300] + v[304] - v[305]
        dydt[V.ERKPP] = +v[115] - v[117] - v[298] - v[303]
        dydt[V.SOSP] = +v[118] - v[119]
        dydt[V.Gab1P] = +v[120] - v[121]
        dydt[V.E11P_Grb2_Gab1P] = +v[122] - v[135] - v[161] - v[212] - v[280]
        dydt[V.E12P_Grb2_Gab1P] = +v[123] - v[136] - v[162] - v[213]
        dydt[V.E23P_Grb2_Gab1P] = +v[124] - v[137] - v[163] - v[214]
        dydt[V.E24P_Grb2_Gab1P] = +v[125] - v[138] - v[164] - v[215]
        dydt[V.E34P_Grb2_Gab1P] = +v[126] - v[139] - v[165] - v[216]
        dydt[V.E44P_Grb2_Gab1P] = +v[127] - v[140] - v[166] - v[217]
        dydt[V.E11P_ShcP_Grb2_Gab1P] = +v[128] - v[141] - v[167] - v[218] - v[281]
        dydt[V.E12P_ShcP_Grb2_Gab1P] = +v[129] - v[142] - v[168] - v[219]
        dydt[V.E23P_ShcP_Grb2_Gab1P] = +v[130] - v[143] - v[169] - v[220]
        dydt[V.E24P_ShcP_Grb2_Gab1P] = +v[131] - v[144] - v[170] - v[221]
        dydt[V.E34P_ShcP_Grb2_Gab1P] = +v[132] - v[145] - v[171] - v[222]
        dydt[V.E44P_ShcP_Grb2_Gab1P] = +v[133] - v[146] - v[172] - v[223]
        dydt[V.PIP3_Gab1] = -v[134] + v[160] + v[259]
        dydt[V.PIP3_Gab1P] = +v[134] - v[147] - v[173] - v[224]
        dydt[V.E11P_Grb2_Gab1P_PTP1B] = +v[135] - v[148] - v[282]
        dydt[V.E12P_Grb2_Gab1P_PTP1B] = +v[136] - v[149]
        dydt[V.E23P_Grb2_Gab1P_PTP1B] = +v[137] - v[150]
        dydt[V.E24P_Grb2_Gab1P_PTP1B] = +v[138] - v[151]
        dydt[V.E34P_Grb2_Gab1P_PTP1B] = +v[139] - v[152]
        dydt[V.E44P_Grb2_Gab1P_PTP1B] = +v[140] - v[153]
        dydt[V.E11P_ShcP_Grb2_Gab1P_PTP1B] = +v[141] - v[154] - v[283]
        dydt[V.E12P_ShcP_Grb2_Gab1P_PTP1B] = +v[142] - v[155]
        dydt[V.E23P_ShcP_Grb2_Gab1P_PTP1B] = +v[143] - v[156]
        dydt[V.E24P_ShcP_Grb2_Gab1P_PTP1B] = +v[144] - v[157]
        dydt[V.E34P_ShcP_Grb2_Gab1P_PTP1B] = +v[145] - v[158]
        dydt[V.E44P_ShcP_Grb2_Gab1P_PTP1B] = +v[146] - v[159]
        dydt[V.PIP3_Gab1P_PTP1B] = +v[147] - v[160]
        dydt[V.E11P_Grb2_Gab1P_RasGAP] = +v[161] - v[180] + v[199] - v[284]
        dydt[V.E12P_Grb2_Gab1P_RasGAP] = +v[162] - v[181] + v[200]
        dydt[V.E23P_Grb2_Gab1P_RasGAP] = +v[163] - v[182] + v[201]
        dydt[V.E24P_Grb2_Gab1P_RasGAP] = +v[164] - v[183] + v[202]
        dydt[V.E34P_Grb2_Gab1P_RasGAP] = +v[165] - v[184] + v[203]
        dydt[V.E44P_Grb2_Gab1P_RasGAP] = +v[166] - v[185] + v[204]
        dydt[V.E11P_ShcP_Grb2_Gab1P_RasGAP] = +v[167] - v[186] + v[205] - v[285]
        dydt[V.E12P_ShcP_Grb2_Gab1P_RasGAP] = +v[168] - v[187] + v[206]
        dydt[V.E23P_ShcP_Grb2_Gab1P_RasGAP] = +v[169] - v[188] + v[207]
        dydt[V.E24P_ShcP_Grb2_Gab1P_RasGAP] = +v[170] - v[189] + v[208]
        dydt[V.E34P_ShcP_Grb2_Gab1P_RasGAP] = +v[171] - v[190] + v[209]
        dydt[V.E44P_ShcP_Grb2_Gab1P_RasGAP] = +v[172] - v[191] + v[210]
        dydt[V.PIP3_Gab1P_RasGAP] = +v[173] - v[192] + v[211]
        dydt[V.E11P_RasGAP_RasGTP] = +v[174] - v[193] - v[286]
        dydt[V.E12P_RasGAP_RasGTP] = +v[175] - v[194]
        dydt[V.E23P_RasGAP_RasGTP] = +v[176] - v[195]
        dydt[V.E24P_RasGAP_RasGTP] = +v[177] - v[196]
        dydt[V.E34P_RasGAP_RasGTP] = +v[178] - v[197]
        dydt[V.E44P_RasGAP_RasGTP] = +v[179] - v[198]
        dydt[V.E11P_Grb2_Gab1P_RasGAP_RasGTP] = +v[180] - v[199] - v[287]
        dydt[V.E12P_Grb2_Gab1P_RasGAP_RasGTP] = +v[181] - v[200]
        dydt[V.E23P_Grb2_Gab1P_RasGAP_RasGTP] = +v[182] - v[201]
        dydt[V.E24P_Grb2_Gab1P_RasGAP_RasGTP] = +v[183] - v[202]
        dydt[V.E34P_Grb2_Gab1P_RasGAP_RasGTP] = +v[184] - v[203]
        dydt[V.E44P_Grb2_Gab1P_RasGAP_RasGTP] = +v[185] - v[204]
        dydt[V.E11P_ShcP_Grb2_Gab1P_RasGAP_RasGTP] = +v[186] - v[205] - v[288]
        dydt[V.E12P_ShcP_Grb2_Gab1P_RasGAP_RasGTP] = +v[187] - v[206]
        dydt[V.E23P_ShcP_Grb2_Gab1P_RasGAP_RasGTP] = +v[188] - v[207]
        dydt[V.E24P_ShcP_Grb2_Gab1P_RasGAP_RasGTP] = +v[189] - v[208]
        dydt[V.E34P_ShcP_Grb2_Gab1P_RasGAP_RasGTP] = +v[190] - v[209]
        dydt[V.E44P_ShcP_Grb2_Gab1P_RasGAP_RasGTP] = +v[191] - v[210]
        dydt[V.PIP3_Gab1P_RasGAP_RasGTP] = +v[192] - v[211]
        dydt[V.E11P_Grb2_Gab1P_PI3K] = +v[212] - v[229] + v[246] - v[289]
        dydt[V.E12P_Grb2_Gab1P_PI3K] = +v[213] - v[230] + v[247]
        dydt[V.E23P_Grb2_Gab1P_PI3K] = +v[214] - v[231] + v[248]
        dydt[V.E24P_Grb2_Gab1P_PI3K] = +v[215] - v[232] + v[249]
        dydt[V.E34P_Grb2_Gab1P_PI3K] = +v[216] - v[233] + v[250]
        dydt[V.E44P_Grb2_Gab1P_PI3K] = +v[217] - v[234] + v[251]
        dydt[V.E11P_ShcP_Grb2_Gab1P_PI3K] = +v[218] - v[235] + v[252] - v[290]
        dydt[V.E12P_ShcP_Grb2_Gab1P_PI3K] = +v[219] - v[236] + v[253]
        dydt[V.E23P_ShcP_Grb2_Gab1P_PI3K] = +v[220] - v[237] + v[254]
        dydt[V.E24P_ShcP_Grb2_Gab1P_PI3K] = +v[221] - v[238] + v[255]
        dydt[V.E34P_ShcP_Grb2_Gab1P_PI3K] = +v[222] - v[239] + v[256]
        dydt[V.E44P_ShcP_Grb2_Gab1P_PI3K] = +v[223] - v[240] + v[257]
        dydt[V.PIP3_Gab1P_PI3K] = +v[224] - v[241] + v[258]
        dydt[V.PIP2] = (
            -v[225]
            - v[226]
            - v[227]
            - v[228]
            - v[229]
            - v[230]
            - v[231]
            - v[232]
            - v[233]
            - v[234]
            - v[235]
            - v[236]
            - v[237]
            - v[238]
            - v[239]
            - v[240]
            - v[241]
            + v[260]
        )
        dydt[V.E23P_PI3K_PIP2] = +v[225] - v[242]
        dydt[V.E24P_PI3K_PIP2] = +v[226] - v[243]
        dydt[V.E34P_PI3K_PIP2] = +v[227] - v[244]
        dydt[V.E44P_PI3K_PIP2] = +v[228] - v[245]
        dydt[V.E11P_Grb2_Gab1P_PI3K_PIP2] = +v[229] - v[246]
        dydt[V.E12P_Grb2_Gab1P_PI3K_PIP2] = +v[230] - v[247]
        dydt[V.E23P_Grb2_Gab1P_PI3K_PIP2] = +v[231] - v[248]
        dydt[V.E24P_Grb2_Gab1P_PI3K_PIP2] = +v[232] - v[249]
        dydt[V.E34P_Grb2_Gab1P_PI3K_PIP2] = +v[233] - v[250]
        dydt[V.E44P_Grb2_Gab1P_PI3K_PIP2] = +v[234] - v[251]
        dydt[V.E11P_ShcP_Grb2_Gab1P_PI3K_PIP2] = +v[235] - v[252]
        dydt[V.E12P_ShcP_Grb2_Gab1P_PI3K_PIP2] = +v[236] - v[253]
        dydt[V.E23P_ShcP_Grb2_Gab1P_PI3K_PIP2] = +v[237] - v[254]
        dydt[V.E24P_ShcP_Grb2_Gab1P_PI3K_PIP2] = +v[238] - v[255]
        dydt[V.E34P_ShcP_Grb2_Gab1P_PI3K_PIP2] = +v[239] - v[256]
        dydt[V.E44P_ShcP_Grb2_Gab1P_PI3K_PIP2] = +v[240] - v[257]
        dydt[V.PIP3_Gab1P_PI3K_PIP2] = +v[241] - v[258]
        dydt[V.PIP3] = (
            +v[242]
            + v[243]
            + v[244]
            + v[245]
            + v[246]
            + v[247]
            + v[248]
            + v[249]
            + v[250]
            + v[251]
            + v[252]
            + v[253]
            + v[254]
            + v[255]
            + v[256]
            + v[257]
            + v[258]
            - v[259]
            - v[260]
        )
        dydt[V.Akt] = -v[261] + v[262]
        dydt[V.AktP] = +v[261] - v[262] - v[267]
        dydt[V.RafP_inactive] = +v[263] - v[265]
        dydt[V.RafPP_inactive] = +v[264] - v[266]
        dydt[V.PreduspmRNA] = +v[291] - v[292]
        dydt[V.duspmRNA] = +v[292] * (0.22 / 0.94) - v[294]
        dydt[V.DUSP] = +v[293] - v[295] - v[296] - v[298] + v[299] - v[300] + v[301] - v[302]
        dydt[V.DUSPP] = +v[295] - v[297] - v[303] + v[304] - v[305] + v[306] - v[307]
        dydt[V.DUSP_ERKPP] = +v[298] - v[299]
        dydt[V.DUSP_ERKP] = +v[300] - v[301]
        dydt[V.DUSP_ERK] = +v[302]
        dydt[V.DUSPP_ERKPP] = +v[303] - v[304]
        dydt[V.DUSPP_ERKP] = +v[305] - v[306]
        dydt[V.DUSPP_ERK] = +v[307]
        dydt[V.GSK3b] = -v[308] + v[309]
        dydt[V.GSK3bP] = +v[308] - v[309]
        dydt[V.PrecmycmRNA] = +v[310] - v[311]
        dydt[V.cmycmRNA] = +v[311] * (0.22 / 0.94) - v[312]
        dydt[V.cMyc] = +v[313] - v[314] + v[315] - v[318]
        dydt[V.cMycP] = +v[314] - v[315] - v[316] + v[317]
        dydt[V.cMycPP] = +v[316] - v[317] - v[319]

        return dydt


def param_values():
    """Parameter values"""
    x = [1] * C.NUM
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
    x[C.kr44] = 0
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
    x[C.kr56] = 0
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
    x[C.kr98] = 0
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
    x[C.kr122] = 0
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
    x[C.kr148] = 0
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
    x[C.kr193] = 0
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
    x[C.kr242] = 0
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
    x[C.n291] = 1
    x[C.kr292] = 0
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
    x[C.n310] = 1
    x[C.kr311] = 0

    return x


def initial_values():
    """Values of the initial condition"""
    y0 = [0] * V.NUM
    y0[V.ErbB1] = 1
    y0[V.ErbB3] = 1
    y0[V.ErbB4] = 1
    y0[V.ErbB2] = 1
    y0[V.Grb2] = 1
    y0[V.Shc] = 1
    y0[V.RasGAP] = 1
    y0[V.PI3K] = 1
    y0[V.PTP1B] = 1
    y0[V.SOS] = 1
    y0[V.Gab1] = 1
    y0[V.RasGDP] = 1
    y0[V.Raf] = 1
    y0[V.MEK] = 1
    y0[V.ERK] = 1
    y0[V.PIP2] = 1
    y0[V.PTEN] = 1
    y0[V.Akt] = 1
    y0[V.GSK3b] = 1

    return y0
