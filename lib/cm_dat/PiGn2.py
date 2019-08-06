
from matplotlib.colors import LinearSegmentedColormap
from numpy import nan, inf

# Used to reconstruct the colormap in pycam02ucs.cm.viscm
parameters = {'xp': [18.64284884140784, 18.64284884140784, 95.795068084779928, -77.797425212807227, -30.88053513237827, -10.810532153528101],
              'yp': [-19.924784033363096, -24.355823652070285, -37.648942508191809, 27.513404825737297, 38.981977956508814, 15.784182305630054],
              'min_Jp': 15,
              'max_Jp': 99.9574241618}

cm_data = [[ 0.20356191,  0.01159375,  0.27620087],
       [ 0.211824  ,  0.01431358,  0.2875344 ],
       [ 0.22071917,  0.01637861,  0.29826109],
       [ 0.22978039,  0.01838736,  0.30861388],
       [ 0.23891998,  0.02040174,  0.31873239],
       [ 0.24811834,  0.02240309,  0.32869698],
       [ 0.25736106,  0.02439289,  0.33854162],
       [ 0.26665551,  0.02632418,  0.34830521],
       [ 0.27599454,  0.02819877,  0.35799632],
       [ 0.28537659,  0.03000696,  0.3676242 ],
       [ 0.29480565,  0.03172411,  0.37720001],
       [ 0.30427924,  0.03334772,  0.3867247 ],
       [ 0.31379972,  0.03486142,  0.39620227],
       [ 0.32336758,  0.03625606,  0.40563284],
       [ 0.33298334,  0.03752391,  0.41501492],
       [ 0.34264679,  0.03866162,  0.42434472],
       [ 0.35236125,  0.03965456,  0.43362001],
       [ 0.3621215 ,  0.04052167,  0.44282877],
       [ 0.37193548,  0.04122178,  0.45196839],
       [ 0.38179375,  0.04181196,  0.46101821],
       [ 0.39170109,  0.04228995,  0.46996574],
       [ 0.40165495,  0.04268876,  0.47878681],
       [ 0.41164526,  0.04308253,  0.48744444],
       [ 0.42166753,  0.04353962,  0.49589434],
       [ 0.43171223,  0.04417439,  0.50406962],
       [ 0.44175593,  0.04519922,  0.51186297],
       [ 0.45174203,  0.04702377,  0.5190952 ],
       [ 0.46156946,  0.05032482,  0.52545733],
       [ 0.47096611,  0.05646592,  0.53039926],
       [ 0.47936462,  0.06731586,  0.53324927],
       [ 0.4862797 ,  0.08262543,  0.53429877],
       [ 0.49205833,  0.09926202,  0.53477146],
       [ 0.49722287,  0.11536398,  0.53530106],
       [ 0.50205176,  0.13052994,  0.53602884],
       [ 0.50669079,  0.14476353,  0.53699582],
       [ 0.51120642,  0.15820497,  0.53817387],
       [ 0.51562687,  0.17100511,  0.53951857],
       [ 0.51998716,  0.18324377,  0.5410346 ],
       [ 0.52429384,  0.19502618,  0.54268768],
       [ 0.52857657,  0.20638754,  0.54449892],
       [ 0.5328247 ,  0.21741899,  0.5464227 ],
       [ 0.53704804,  0.22815917,  0.54845706],
       [ 0.54125386,  0.23864204,  0.55059857],
       [ 0.54544897,  0.2488954 ,  0.55284521],
       [ 0.54963858,  0.25894398,  0.55519403],
       [ 0.55382015,  0.26881818,  0.55763163],
       [ 0.55799636,  0.27853701,  0.56015425],
       [ 0.56216946,  0.28811708,  0.56275839],
       [ 0.56634738,  0.29756595,  0.56544948],
       [ 0.57052712,  0.30690243,  0.56821752],
       [ 0.57470967,  0.31613847,  0.57105908],
       [ 0.57889595,  0.32528467,  0.5739711 ],
       [ 0.58308679,  0.33435043,  0.57695082],
       [ 0.58728301,  0.34334411,  0.57999584],
       [ 0.5914854 ,  0.35227319,  0.58310402],
       [ 0.5957078 ,  0.36113184,  0.5862903 ],
       [ 0.59994117,  0.36993621,  0.58953984],
       [ 0.60418332,  0.37869462,  0.59284731],
       [ 0.60843504,  0.387412  ,  0.59621142],
       [ 0.61270633,  0.39608475,  0.59964236],
       [ 0.61699822,  0.40471725,  0.60313906],
       [ 0.62130193,  0.41332183,  0.60668876],
       [ 0.62561839,  0.42190181,  0.6102909 ],
       [ 0.62996578,  0.43044618,  0.61396527],
       [ 0.63432983,  0.43897073,  0.61769315],
       [ 0.63870896,  0.44748013,  0.62147124],
       [ 0.64311686,  0.45596693,  0.62531387],
       [ 0.647549  ,  0.46443804,  0.62921427],
       [ 0.65199868,  0.47290157,  0.63316353],
       [ 0.65647941,  0.48135018,  0.63717556],
       [ 0.66098721,  0.48978974,  0.64124454],
       [ 0.66551518,  0.49822785,  0.64536177],
       [ 0.67007973,  0.5066552 ,  0.64954404],
       [ 0.67467079,  0.5150815 ,  0.65377923],
       [ 0.67928564,  0.52351088,  0.6580635 ],
       [ 0.68394325,  0.53193248,  0.66241601],
       [ 0.68862541,  0.54036101,  0.66681648],
       [ 0.69334146,  0.54879218,  0.67127386],
       [ 0.6980951 ,  0.55722559,  0.67579109],
       [ 0.70287651,  0.56566967,  0.68035717],
       [ 0.70770256,  0.57411533,  0.68498862],
       [ 0.71256052,  0.5825727 ,  0.68967159],
       [ 0.71745624,  0.59103981,  0.69441134],
       [ 0.72239417,  0.59951566,  0.69921155],
       [ 0.72736512,  0.60800769,  0.70406238],
       [ 0.7323858 ,  0.61650728,  0.70897973],
       [ 0.73744207,  0.62502479,  0.71394896],
       [ 0.74254368,  0.63355589,  0.71897905],
       [ 0.74768885,  0.64210328,  0.72406761],
       [ 0.75287612,  0.6506694 ,  0.72921268],
       [ 0.75811385,  0.65925087,  0.73442173],
       [ 0.7633925 ,  0.6678548 ,  0.73968507],
       [ 0.76872554,  0.6764749 ,  0.74501501],
       [ 0.77410273,  0.68511856,  0.75040131],
       [ 0.77953256,  0.69378241,  0.75585148],
       [ 0.78501232,  0.70246947,  0.76136248],
       [ 0.79054382,  0.71118015,  0.76693553],
       [ 0.79613039,  0.71991407,  0.77257322],
       [ 0.80176861,  0.72867454,  0.77827193],
       [ 0.80746651,  0.73745858,  0.78403855],
       [ 0.81321731,  0.74627134,  0.78986636],
       [ 0.81903089,  0.75510888,  0.79576386],
       [ 0.82490093,  0.76397604,  0.80172475],
       [ 0.83083466,  0.77287038,  0.80775496],
       [ 0.83682978,  0.78179458,  0.81385169],
       [ 0.84289039,  0.79074793,  0.82001783],
       [ 0.84901741,  0.79973142,  0.82625318],
       [ 0.85521285,  0.80874552,  0.83255826],
       [ 0.86148022,  0.81778998,  0.83893454],
       [ 0.86782046,  0.82686593,  0.8453806 ],
       [ 0.87423939,  0.83597214,  0.85189844],
       [ 0.88073805,  0.84510997,  0.85848443],
       [ 0.88732463,  0.85427738,  0.86513912],
       [ 0.89400204,  0.86347544,  0.8718548 ],
       [ 0.90078036,  0.87270216,  0.87862479],
       [ 0.90766534,  0.88195864,  0.88543018],
       [ 0.91466348,  0.89124663,  0.89224305],
       [ 0.92176378,  0.90057687,  0.8990171 ],
       [ 0.92891436,  0.90997627,  0.90570852],
       [ 0.93599949,  0.91949018,  0.91234095],
       [ 0.94292814,  0.92914063,  0.91906249],
       [ 0.94976852,  0.93888848,  0.92599282],
       [ 0.95664955,  0.94868791,  0.93308579],
       [ 0.96361876,  0.95853131,  0.94024436],
       [ 0.9706655 ,  0.9684299 ,  0.94741483],
       [ 0.97776646,  0.97839614,  0.95457898],
       [ 0.98490592,  0.98843757,  0.96173656],
       [ 0.99207124,  0.99856027,  0.96888909],
       [ 0.9900321 ,  0.99955205,  0.96692503],
       [ 0.97879778,  0.99140704,  0.95586087],
       [ 0.96760023,  0.98333455,  0.9448253 ],
       [ 0.95643452,  0.97533498,  0.9338172 ],
       [ 0.94530153,  0.96740602,  0.92284003],
       [ 0.93419817,  0.95954735,  0.91189267],
       [ 0.92312461,  0.95175718,  0.90097664],
       [ 0.91207969,  0.94403436,  0.89009187],
       [ 0.90106173,  0.93637807,  0.87923748],
       [ 0.89007213,  0.92878604,  0.86841547],
       [ 0.87910629,  0.92125881,  0.85762188],
       [ 0.86816907,  0.91379258,  0.84686184],
       [ 0.85725188,  0.90638977,  0.83612739],
       [ 0.84636293,  0.89904506,  0.82542679],
       [ 0.83549192,  0.89176164,  0.81475038],
       [ 0.82464657,  0.88453448,  0.8041059 ],
       [ 0.81381954,  0.8773654 ,  0.79348653],
       [ 0.80301385,  0.87025155,  0.78289544],
       [ 0.79222757,  0.86319231,  0.77233102],
       [ 0.78145735,  0.85618764,  0.76179037],
       [ 0.77070846,  0.84923378,  0.75127875],
       [ 0.75996929,  0.84233428,  0.7407854 ],
       [ 0.74925227,  0.83548231,  0.73032242],
       [ 0.73854154,  0.82868317,  0.71987505],
       [ 0.72784886,  0.82193037,  0.70945476],
       [ 0.71716557,  0.8152261 ,  0.69905369],
       [ 0.70649211,  0.80856867,  0.68867261],
       [ 0.69583251,  0.80195493,  0.67831565],
       [ 0.68517281,  0.79538911,  0.66797013],
       [ 0.67452871,  0.7888633 ,  0.65765111],
       [ 0.66388133,  0.78238361,  0.64734144],
       [ 0.65324244,  0.77594377,  0.63705253],
       [ 0.64260754,  0.76954408,  0.62678063],
       [ 0.63196744,  0.76318653,  0.61651778],
       [ 0.62133693,  0.75686386,  0.60627814],
       [ 0.61069262,  0.75058341,  0.59604077],
       [ 0.60005046,  0.74433751,  0.58582091],
       [ 0.58940626,  0.73812615,  0.57561533],
       [ 0.57874467,  0.73195331,  0.56541072],
       [ 0.56808449,  0.72581045,  0.55522488],
       [ 0.55740656,  0.71970278,  0.54504118],
       [ 0.54671388,  0.71362744,  0.53486317],
       [ 0.5360168 ,  0.70757908,  0.52470097],
       [ 0.52528743,  0.70156553,  0.51453046],
       [ 0.51454262,  0.69557908,  0.5043677 ],
       [ 0.50378609,  0.68961671,  0.49421699],
       [ 0.49298573,  0.68368703,  0.48405086],
       [ 0.4821635 ,  0.67778078,  0.47389   ],
       [ 0.47132031,  0.67189573,  0.46373637],
       [ 0.46042717,  0.66603885,  0.45356576],
       [ 0.44949689,  0.66020382,  0.44339091],
       [ 0.43853387,  0.65438712,  0.4332171 ],
       [ 0.42752351,  0.64859078,  0.42303316],
       [ 0.41644773,  0.6428175 ,  0.41282524],
       [ 0.4053242 ,  0.63705954,  0.40261046],
       [ 0.39414751,  0.63131584,  0.39238607],
       [ 0.38290164,  0.62558805,  0.38214041],
       [ 0.3715674 ,  0.61987815,  0.37185949],
       [ 0.36015888,  0.61417919,  0.36155777],
       [ 0.34866788,  0.60849005,  0.35123087],
       [ 0.33708517,  0.60280956,  0.34087376],
       [ 0.32540041,  0.59713655,  0.33048075],
       [ 0.31358502,  0.59147362,  0.3200312 ],
       [ 0.3016406 ,  0.58581589,  0.30952992],
       [ 0.28955836,  0.58016067,  0.2989732 ],
       [ 0.27732332,  0.57450619,  0.28835235],
       [ 0.26491898,  0.56885045,  0.27765755],
       [ 0.25232742,  0.56319115,  0.26687775],
       [ 0.23952942,  0.55752558,  0.2560005 ],
       [ 0.22650491,  0.55185047,  0.2450119 ],
       [ 0.21323389,  0.54616183,  0.23389649],
       [ 0.1996605 ,  0.54046058,  0.22260608],
       [ 0.18578167,  0.53473797,  0.21112986],
       [ 0.17160327,  0.52898403,  0.19945406],
       [ 0.1570897 ,  0.5231949 ,  0.18750841],
       [ 0.14227834,  0.51735769,  0.17524139],
       [ 0.12727245,  0.51145601,  0.16256709],
       [ 0.11240215,  0.50545936,  0.14940471],
       [ 0.09845334,  0.49931851,  0.13562236],
       [ 0.08741484,  0.49293175,  0.12111657],
       [ 0.08362166,  0.48608021,  0.10625853],
       [ 0.0904406 ,  0.47852424,  0.0931874 ],
       [ 0.10169467,  0.47051941,  0.08390117],
       [ 0.11198348,  0.46242782,  0.07739509],
       [ 0.12060725,  0.45435807,  0.07250507],
       [ 0.12776394,  0.44633774,  0.06862278],
       [ 0.13374884,  0.43836978,  0.06542873],
       [ 0.13877366,  0.43045419,  0.0627322 ],
       [ 0.14295863,  0.42259431,  0.0603918 ],
       [ 0.14647844,  0.41478183,  0.05834372],
       [ 0.14943588,  0.40701306,  0.05653535],
       [ 0.15183099,  0.39929614,  0.05487769],
       [ 0.1537765 ,  0.39162196,  0.05336495],
       [ 0.15531922,  0.38398881,  0.05196942],
       [ 0.15649346,  0.3763958 ,  0.05066515],
       [ 0.1573297 ,  0.36884179,  0.04943147],
       [ 0.15785491,  0.36132553,  0.04825184],
       [ 0.1580928 ,  0.35384568,  0.04711284],
       [ 0.15806422,  0.34640089,  0.04600348],
       [ 0.15778745,  0.33898981,  0.04491464],
       [ 0.15727855,  0.33161111,  0.04383863],
       [ 0.15655163,  0.32426348,  0.04276887],
       [ 0.15561907,  0.31694566,  0.04169967],
       [ 0.15449179,  0.3096564 ,  0.04062598],
       [ 0.15316687,  0.30239717,  0.03951782],
       [ 0.15166412,  0.29516448,  0.03840869],
       [ 0.14999378,  0.28795657,  0.03730232],
       [ 0.14816254,  0.28077223,  0.03619636],
       [ 0.14617631,  0.27361028,  0.03508879],
       [ 0.14403461,  0.26647088,  0.03396992],
       [ 0.14174101,  0.25935315,  0.03283663],
       [ 0.13930833,  0.25225382,  0.03170003],
       [ 0.13674003,  0.24517165,  0.03055924],
       [ 0.13403792,  0.23810569,  0.02941182],
       [ 0.13119685,  0.23105667,  0.0282457 ],
       [ 0.1282301 ,  0.22402045,  0.0270768 ],
       [ 0.12513946,  0.21699562,  0.02590494],
       [ 0.12192416,  0.20998134,  0.02472669],
       [ 0.11858276,  0.20297681,  0.02353821],
       [ 0.11512295,  0.19597857,  0.02235015],
       [ 0.11154527,  0.18898486,  0.02116301],
       [ 0.10784607,  0.181995  ,  0.01997135],
       [ 0.10402905,  0.17500594,  0.0187817 ],
       [ 0.10009587,  0.16801498,  0.01759808],
       [ 0.09604479,  0.16102014,  0.01641988],
       [ 0.0918742 ,  0.15401915,  0.01524727],
       [ 0.08758631,  0.14700822,  0.01408717],
       [ 0.08317967,  0.13998427,  0.01294147]]

PiGn2 = LinearSegmentedColormap.from_list('PiGn2', cm_data)
PiGn2_r = LinearSegmentedColormap.from_list('PiGn2_r', cm_data[::-1])

# Make a clean "from <colormap> import *"
__all__ = ['PiGn2', 'PiGn2_r']

# For use with viscm
test_cm = PiGn2

#
