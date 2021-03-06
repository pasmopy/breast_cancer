module Exp
using StatsBase
using Statistics
include("./observable.jl")

const t = [0, 5, 15, 30, 45, 60, 90, 120]

function norm01(;egf1::Vector{Int}, hrg1::Vector{Int}, egf2::Vector{Int},
    hrg2::Vector{Int}, egf3::Vector{Int}, hrg3::Vector{Int})
    data1::Matrix{Float64} = hcat(egf1, hrg1)
    data2::Matrix{Float64} = hcat(egf2, hrg2)
    data3::Matrix{Float64} = hcat(egf3, hrg3)

    data1 .= data1 ./ maximum(data1)
    data2 .= data2 ./ maximum(data2)
    data3 .= data3 ./ maximum(data3)

    egf_ave::Vector{Float64} = zeros(length(t))
    hrg_ave::Vector{Float64} = zeros(length(t))
    egf_sem::Vector{Float64} = zeros(length(t))
    hrg_sem::Vector{Float64} = zeros(length(t))

    for i in eachindex(t)
        egf_ave[i] = mean([data1[i,1], data2[i,1], data3[i,1]])
        hrg_ave[i] = mean([data1[i,2], data2[i,2], data3[i,2]])
    end

    ave_vec::Vector{Float64} = vcat(egf_ave, hrg_ave)
    ave_min::Float64 = minimum(ave_vec)
    ave_max::Float64 = maximum(ave_vec)

    data1 .= (data1 .- ave_min) ./ (ave_max .- ave_min)
    data2 .= (data2 .- ave_min) ./ (ave_max .- ave_min)
    data3 .= (data3 .- ave_min) ./ (ave_max .- ave_min)

    for i in eachindex(t)
        egf_ave[i] = mean([data1[i,1], data2[i,1], data3[i,1]])
        hrg_ave[i] = mean([data1[i,2], data2[i,2], data3[i,2]])
        egf_sem[i] = std([data1[i,1], data2[i,1], data3[i,1]]) ./ sqrt(3)
        hrg_sem[i] = std([data1[i,2], data2[i,2], data3[i,2]]) ./ sqrt(3)
    end

    return egf_ave, hrg_ave, egf_sem, hrg_sem
end

experiments = Array{Dict{String,Array{Float64,1}},1}(undef, length(observables))
error_bars = Array{Dict{String,Array{Float64,1}},1}(undef, length(observables))

mcf7_data_pAkt = norm01(
    egf1=[32101, 156970, 90301, 76709, 63640, 52536, 46414, 57329],
    hrg1=[32101, 565508, 551901, 560064, 489678, 408802, 425323, 451502],
    egf2=[11612, 96189, 43622, 43238, 41007, 29902, 19255, 35079],
    hrg2=[11612, 397931, 432609, 417622, 434519, 509919, 361041, 292523],
    egf3=[66038, 208525, 102689, 117308, 125158, 92086, 68587, 78252],
    hrg3=[66038, 563079, 573540, 521062, 447462, 383774, 434807, 409615],
)
bt474_data_pAkt = norm01(
    egf1=[405198, 356865, 321475, 383445, 346872, 328052, 299123, 316633],
    hrg1=[405198, 357121, 419948, 488508, 495214, 443710, 402765, 451831],
    egf2=[432524, 619289, 581376, 481899, 429541, 399922, 376170, 334923],
    hrg2=[432524, 410919, 413878, 390581, 405359, 408471, 373108, 515120],
    egf3=[150446, 435897, 466378, 443105, 415827, 381441, 398841, 413906],
    hrg3=[150446, 556176, 560385, 539165, 589297, 552227, 540005, 539010],
)
mdamb231_data_pAkt = norm01(
    egf1=[86491, 826975, 575400, 354446, 143728, 107326, 88082, 108892],
    hrg1=[86491, 85990, 114224, 94636, 66139, 92359, 89540, 126099],
    egf2=[56004, 816102, 296684, 165282, 92711, 83090, 55002, 50605],
    hrg2=[56004, 50640, 64646, 55814, 70627, 66657, 68994, 46922],
    egf3=[36061, 729961, 218315, 166109, 68746, 53331, 45191, 21852],
    hrg3=[36061, 75029, 174050, 89546, 88262, 99161, 78673, 58412],
)
skbr3_data_pAkt = norm01(
    egf1=[185525, 566558, 675230, 679968, 669038, 596447, 547016, 439083],
    hrg1=[185525, 530483, 594649, 637933, 733390, 607193, 649356, 709869],
    egf2=[208726, 431186, 439068, 424423, 452676, 414664, 388795, 310282],
    hrg2=[208726, 467749, 509839, 497772, 470593, 429163, 340911, 314740],
    egf3=[290628, 699828, 586893, 512882, 460166, 432428, 442622, 383550],
    hrg3=[290628, 644516, 632090, 591780, 655548, 689898, 632575, 534658],
)
experiments[observables_index("Phosphorylated_Akt")] = Dict(
    "MCF7_EGF" => mcf7_data_pAkt[1],
    "MCF7_HRG" => mcf7_data_pAkt[2],
    "BT474_EGF" => bt474_data_pAkt[1],
    "BT474_HRG" => bt474_data_pAkt[2],
    "MDAMB231_EGF" => mdamb231_data_pAkt[1],
    "MDAMB231_HRG" => mdamb231_data_pAkt[2],
    "SKBR3_EGF" => skbr3_data_pAkt[1],
    "SKBR3_HRG" => skbr3_data_pAkt[2],
)
error_bars[observables_index("Phosphorylated_Akt")] = Dict(
    "MCF7_EGF" => mcf7_data_pAkt[3],
    "MCF7_HRG" => mcf7_data_pAkt[4],
    "BT474_EGF" => bt474_data_pAkt[3],
    "BT474_HRG" => bt474_data_pAkt[4],
    "MDAMB231_EGF" => mdamb231_data_pAkt[3],
    "MDAMB231_HRG" => mdamb231_data_pAkt[4],
    "SKBR3_EGF" => skbr3_data_pAkt[3],
    "SKBR3_HRG" => skbr3_data_pAkt[4],
)

mcf7_data_pERK = norm01(
    egf1=[65481, 446949, 221435, 283171, 265152, 266056, 204912, 188972],
    hrg1=[65481, 698717, 766252, 710005, 693622, 691856, 522173, 334410],
    egf2=[41927, 507623, 169918, 193671, 164088, 145916, 110844, 130362],
    hrg2=[41927, 605118, 699511, 654697, 579863, 490649, 299946, 229297],
    egf3=[118995, 807929, 338665, 267160, 253820, 230200, 157620, 153112],
    hrg3=[118995, 710436, 673318, 615206, 612686, 523198, 390301, 257664],
)
bt474_data_pERK = norm01(
    egf1=[358203, 550378, 633802, 632047, 500267, 394009, 339650, 221411],
    hrg1=[358203, 531893, 703437, 663640, 629213, 612525, 613871, 643056],
    egf2=[355065, 1421536, 1307969, 1101679, 939944, 689539, 507787, 468836],
    hrg2=[355065, 929915, 897601, 924274, 865529, 820386, 673456, 788623],
    egf3=[61593, 631017, 754722, 652440, 575812, 432406, 315961, 259708],
    hrg3=[61593, 480140, 487770, 463604, 438917, 452289, 470624, 531531],
)
mdamb231_data_pERK = norm01(
    egf1=[314472, 504819, 607786, 618492, 475195, 376035, 293988, 324600],
    hrg1=[314472, 156705, 183456, 277862, 141450, 199719, 253331, 407923],
    egf2=[458693, 1001334, 875594, 834259, 782639, 815888, 629576, 539187],
    hrg2=[458693, 322542, 403985, 331734, 263578, 262142, 276371, 313541],
    egf3=[365691, 747932, 937413, 945635, 870059, 706306, 510590, 451927],
    hrg3=[365691, 340876, 428510, 303543, 269653, 195660, 215972, 388050],
)
skbr3_data_pERK = norm01(
    egf1=[182424, 757623, 702797, 710750, 704235, 636329, 680400, 637799],
    hrg1=[182424, 662752, 754617, 725132, 637076, 590793, 610035, 541614],
    egf2=[335742, 721244, 699229, 542494, 486792, 449850, 446215, 452925],
    hrg2=[335742, 588327, 607123, 557027, 490253, 485330, 478285, 530225],
    egf3=[153334, 617406, 425659, 408737, 358565, 338492, 350159, 353366],
    hrg3=[153334, 636281, 559026, 436561, 389041, 360622, 333743, 355748],
)
experiments[observables_index("Phosphorylated_ERK")] = Dict(
    "MCF7_EGF" => mcf7_data_pERK[1],
    "MCF7_HRG" => mcf7_data_pERK[2],
    "BT474_EGF" => bt474_data_pERK[1],
    "BT474_HRG" => bt474_data_pERK[2],
    "MDAMB231_EGF" => mdamb231_data_pERK[1],
    "MDAMB231_HRG" => mdamb231_data_pERK[2],
    "SKBR3_EGF" => skbr3_data_pERK[1],
    "SKBR3_HRG" => skbr3_data_pERK[2],
)
error_bars[observables_index("Phosphorylated_ERK")] = Dict(
    "MCF7_EGF" => mcf7_data_pERK[3],
    "MCF7_HRG" => mcf7_data_pERK[4],
    "BT474_EGF" => bt474_data_pERK[3],
    "BT474_HRG" => bt474_data_pERK[4],
    "MDAMB231_EGF" => mdamb231_data_pERK[3],
    "MDAMB231_HRG" => mdamb231_data_pERK[4],
    "SKBR3_EGF" => skbr3_data_pERK[3],
    "SKBR3_HRG" => skbr3_data_pERK[4],
)

mcf7_data_pcMyc = norm01(
    egf1=[115975, 226001, 166894, 194150, 263331, 235172, 126949, 91142],
    hrg1=[115975, 62515, 81364, 155844, 390689, 664641, 848356, 856941],
    egf2=[185069, 276202, 204012, 234391, 290020, 360762, 325531, 242455],
    hrg2=[185069, 234416, 251732, 333993, 550670, 859790, 939956, 769616],
    egf3=[127244, 186118, 163387, 132053, 192949, 220987, 184381, 151547],
    hrg3=[127244, 110676, 152880, 277206, 461217, 637033, 908235, 712427],
)
bt474_data_pcMyc = norm01(
    egf1=[677642, 631653, 761828, 940801, 824185, 1063834, 1248333, 1126031],
    hrg1=[677642, 672273, 647219, 804473, 958061, 1105501, 1088048, 1323065],
    egf2=[382898, 389470, 359455, 414052, 501284, 540238, 621011, 445658],
    hrg2=[382898, 299678, 327189, 491156, 762165, 898141, 835015, 684669],
    egf3=[321190, 366660, 301538, 338201, 395081, 377184, 334796, 277554],
    hrg3=[321190, 208981, 211540, 332880, 434028, 539873, 598716, 603566],
)
mdamb231_data_pcMyc = norm01(
    egf1=[175088, 147000, 158258, 201315, 324002, 362181, 481797, 475243],
    hrg1=[175088, 239248, 242011, 307948, 345390, 408564, 454948, 443210],
    egf2=[398355, 426051, 420120, 513399, 483150, 516812, 529569, 506936],
    hrg2=[398355, 320668, 323452, 315337, 321620, 365778, 434971, 439945],
    egf3=[887432, 906168, 503726, 913939, 1236533, 1305121, 1000144, 955499],
    hrg3=[887432, 803746, 822021, 739253, 833273, 844401, 844497, 1064656],
)
skbr3_data_pcMyc = norm01(
    egf1=[216306, 260513, 398333, 450702, 457053, 431985, 360345, 234876],
    hrg1=[216306, 189928, 337528, 333291, 370834, 298729, 263847, 198003],
    egf2=[1100082, 874469, 1105692, 1121651, 1132333, 984411, 861916, 928244],
    hrg2=[1100082, 620915, 847173, 1070530, 849862, 753615, 516349, 486978],
    egf3=[624420, 636816, 875382, 972561, 1066021, 916545, 669543, 627495],
    hrg3=[624420, 511242, 678598, 826159, 933370, 839597, 600985, 610415],
)
experiments[observables_index("Phosphorylated_c-Myc")] = Dict(
    "MCF7_EGF" => mcf7_data_pcMyc[1],
    "MCF7_HRG" => mcf7_data_pcMyc[2],
    "BT474_EGF" => bt474_data_pcMyc[1],
    "BT474_HRG" => bt474_data_pcMyc[2],
    "MDAMB231_EGF" => mdamb231_data_pcMyc[1],
    "MDAMB231_HRG" => mdamb231_data_pcMyc[2],
    "SKBR3_EGF" => skbr3_data_pcMyc[1],
    "SKBR3_HRG" => skbr3_data_pcMyc[2],
)
error_bars[observables_index("Phosphorylated_c-Myc")] = Dict(
    "MCF7_EGF" => mcf7_data_pcMyc[3],
    "MCF7_HRG" => mcf7_data_pcMyc[4],
    "BT474_EGF" => bt474_data_pcMyc[3],
    "BT474_HRG" => bt474_data_pcMyc[4],
    "MDAMB231_EGF" => mdamb231_data_pcMyc[3],
    "MDAMB231_HRG" => mdamb231_data_pcMyc[4],
    "SKBR3_EGF" => skbr3_data_pcMyc[3],
    "SKBR3_HRG" => skbr3_data_pcMyc[4],
)

function get_timepoint(obs_name::String)::Vector{Float64}
    if obs_name in observables
        return t
    end
end
end # module