
from matplotlib.colors import LinearSegmentedColormap
from numpy import nan, inf

# Used to reconstruct the colormap in pycam02ucs.cm.viscm
parameters = {'xp': [-3.7709220302195376, 1.8535040487520007, -12.07364624203467, -56.53339524723819, -27.339945599243094, 45.241933800818089, -18.501561760859261, -20.912030080418504, 15.244994712969941, 59.972573531457812, 25.422527617775557],
              'yp': [-33.673859810223433, -65.813437404346473, -9.8370064279154974, 6.2327823691460367, -17.87190082644625, -60.724670951943651, 10.785889194980115, 38.104530149984726, 46.675084175084208, 18.552953780226517, 13.464187327823709],
              'min_Jp': 15,
              'max_Jp': 99.9574241618}

cm_data = [[ 0.06314266,  0.03238354,  0.45534749],
       [ 0.06981308,  0.03453976,  0.47255294],
       [ 0.07632611,  0.03696225,  0.48950025],
       [ 0.08302377,  0.03853789,  0.50708237],
       [ 0.08969162,  0.0399033 ,  0.52474881],
       [ 0.09646965,  0.04048713,  0.54289491],
       [ 0.10322874,  0.04070096,  0.56118397],
       [ 0.11012715,  0.03980538,  0.58009346],
       [ 0.11688737,  0.03890581,  0.59885906],
       [ 0.12388217,  0.03613732,  0.61863326],
       [ 0.13078852,  0.0328362 ,  0.63855006],
       [ 0.13760705,  0.02880303,  0.65873656],
       [ 0.14432399,  0.02372747,  0.67938569],
       [ 0.1508247 ,  0.01749597,  0.70061559],
       [ 0.1566924 ,  0.01100357,  0.72215637],
       [ 0.16066581,  0.00706201,  0.74334331],
       [ 0.15721434,  0.02753911,  0.75587194],
       [ 0.14810885,  0.07260066,  0.75364428],
       [ 0.13883892,  0.10594286,  0.74779735],
       [ 0.13020017,  0.13190283,  0.74144736],
       [ 0.12199245,  0.15368405,  0.73502719],
       [ 0.1143349 ,  0.17258643,  0.72887225],
       [ 0.10727266,  0.18940687,  0.72310405],
       [ 0.10082849,  0.20466288,  0.71775668],
       [ 0.0950218 ,  0.21870643,  0.71282868],
       [ 0.08987465,  0.23178625,  0.70830492],
       [ 0.08541227,  0.24408399,  0.70416592],
       [ 0.08166101,  0.25573603,  0.70039174],
       [ 0.07864497,  0.26684729,  0.69696359],
       [ 0.0763822 ,  0.27750019,  0.69386422],
       [ 0.0748811 ,  0.2877607 ,  0.69107807],
       [ 0.07416111,  0.29767453,  0.68861314],
       [ 0.07417805,  0.30729516,  0.68643258],
       [ 0.07488851,  0.31666306,  0.68451341],
       [ 0.07625057,  0.32580857,  0.68284576],
       [ 0.07821455,  0.33475726,  0.68142158],
       [ 0.08072738,  0.34353039,  0.68023596],
       [ 0.08371479,  0.35215223,  0.67926425],
       [ 0.08712068,  0.36063908,  0.67849998],
       [ 0.0908912 ,  0.36900535,  0.67793673],
       [ 0.09496997,  0.37726591,  0.67756044],
       [ 0.09931058,  0.38543262,  0.67736273],
       [ 0.10387489,  0.39351507,  0.67733936],
       [ 0.10862863,  0.40152221,  0.6774849 ],
       [ 0.11353664,  0.40946404,  0.6777861 ],
       [ 0.11856849,  0.41734976,  0.67822928],
       [ 0.12371075,  0.42518341,  0.67881928],
       [ 0.12894604,  0.43297098,  0.67955069],
       [ 0.13426001,  0.440718  ,  0.68041775],
       [ 0.13963767,  0.44843077,  0.68140891],
       [ 0.14506305,  0.45611654,  0.68250591],
       [ 0.15054018,  0.46377577,  0.68371878],
       [ 0.15606519,  0.47141218,  0.68504127],
       [ 0.16163648,  0.47902915,  0.68646672],
       [ 0.16725481,  0.48662964,  0.68798806],
       [ 0.17292348,  0.49421622,  0.6895979 ],
       [ 0.17864872,  0.50179097,  0.69128847],
       [ 0.18444017,  0.50935542,  0.69305169],
       [ 0.19031159,  0.51691043,  0.69487921],
       [ 0.19628185,  0.52445603,  0.69676251],
       [ 0.20237627,  0.53199116,  0.69869296],
       [ 0.20862829,  0.53951345,  0.70066202],
       [ 0.21508185,  0.5470187 ,  0.70266144],
       [ 0.22178377,  0.55450775,  0.70463772],
       [ 0.22882325,  0.56196628,  0.70660431],
       [ 0.23630693,  0.56938216,  0.70853109],
       [ 0.24438925,  0.57673641,  0.71037649],
       [ 0.25330526,  0.58399828,  0.71207918],
       [ 0.26339957,  0.59111216,  0.71360961],
       [ 0.27520679,  0.5979908 ,  0.71489074],
       [ 0.28934663,  0.60449966,  0.71598712],
       [ 0.30606318,  0.61051603,  0.71718203],
       [ 0.32438849,  0.61608812,  0.71899147],
       [ 0.34294258,  0.62141124,  0.72161719],
       [ 0.36102149,  0.6266403 ,  0.72491544],
       [ 0.37843721,  0.63185534,  0.72868891],
       [ 0.39523443,  0.63708739,  0.73277428],
       [ 0.41148486,  0.64235215,  0.73706338],
       [ 0.42725609,  0.64765885,  0.74148323],
       [ 0.44261095,  0.65301254,  0.74598105],
       [ 0.45760207,  0.65841671,  0.75051804],
       [ 0.47225102,  0.66387877,  0.75507281],
       [ 0.48661467,  0.66939674,  0.75961572],
       [ 0.50072017,  0.67497329,  0.76412995],
       [ 0.5145674 ,  0.68061614,  0.76861226],
       [ 0.52819715,  0.68632286,  0.7730442 ],
       [ 0.54160859,  0.69209964,  0.77742665],
       [ 0.55482143,  0.69794729,  0.78175267],
       [ 0.56783811,  0.70387006,  0.78602533],
       [ 0.58067373,  0.70986869,  0.79024234],
       [ 0.59332441,  0.7159479 ,  0.79441384],
       [ 0.60581039,  0.72210626,  0.7985366 ],
       [ 0.61812084,  0.72834932,  0.80262879],
       [ 0.63027435,  0.73467522,  0.80669007],
       [ 0.6422787 ,  0.74108435,  0.81072773],
       [ 0.65412392,  0.74758102,  0.81476461],
       [ 0.66583243,  0.75416152,  0.81879809],
       [ 0.67741061,  0.76082575,  0.82283823],
       [ 0.68886233,  0.76757398,  0.82689774],
       [ 0.70018967,  0.77440672,  0.83099173],
       [ 0.71140906,  0.78132096,  0.83512069],
       [ 0.72252793,  0.78831574,  0.83929311],
       [ 0.73355377,  0.79538995,  0.84351695],
       [ 0.74449401,  0.80254246,  0.84779957],
       [ 0.75535603,  0.80977211,  0.85214772],
       [ 0.76614699,  0.81707775,  0.85656761],
       [ 0.77687287,  0.82445848,  0.86106621],
       [ 0.78754112,  0.83191308,  0.8656478 ],
       [ 0.79815913,  0.83944035,  0.87031576],
       [ 0.80873271,  0.84703947,  0.87507474],
       [ 0.81926719,  0.85470974,  0.87992919],
       [ 0.82976701,  0.86245073,  0.88488407],
       [ 0.84023519,  0.87026233,  0.8899455 ],
       [ 0.85067637,  0.8781442 ,  0.8951162 ],
       [ 0.86109303,  0.88609656,  0.90040071],
       [ 0.8714863 ,  0.89412012,  0.90580375],
       [ 0.88185541,  0.90221631,  0.91133012],
       [ 0.89219684,  0.91038764,  0.91698416],
       [ 0.90250302,  0.91863841,  0.92276848],
       [ 0.91276065,  0.92697568,  0.92868058],
       [ 0.92294946,  0.93541055,  0.93470524],
       [ 0.93304367,  0.94395899,  0.94080152],
       [ 0.9430306 ,  0.95263748,  0.94687378],
       [ 0.95294655,  0.96145015,  0.95275955],
       [ 0.96290671,  0.97037144,  0.95828236],
       [ 0.97306888,  0.97934976,  0.96335221],
       [ 0.983542  ,  0.98833915,  0.96800635],
       [ 0.99434837,  0.99732037,  0.97235092],
       [ 0.99621985,  0.99708721,  0.96735469],
       [ 0.98910203,  0.98764976,  0.95311834],
       [ 0.9821953 ,  0.97822468,  0.93881318],
       [ 0.97547116,  0.96881874,  0.92446184],
       [ 0.96890646,  0.95943701,  0.91008413],
       [ 0.96248898,  0.95008197,  0.89568332],
       [ 0.95620914,  0.94075497,  0.88126366],
       [ 0.95005925,  0.93145661,  0.86682988],
       [ 0.94403297,  0.92218697,  0.85238707],
       [ 0.93812496,  0.91294571,  0.8379405 ],
       [ 0.93233411,  0.9037315 ,  0.82348823],
       [ 0.92665763,  0.89454321,  0.80903321],
       [ 0.92109321,  0.88537952,  0.79457828],
       [ 0.91563842,  0.87623911,  0.7801272 ],
       [ 0.91029295,  0.86712008,  0.76567994],
       [ 0.90505549,  0.85802075,  0.75123894],
       [ 0.89992467,  0.84893937,  0.73680703],
       [ 0.89489978,  0.83987404,  0.72238598],
       [ 0.88998045,  0.83082269,  0.70797716],
       [ 0.88516517,  0.82178353,  0.69358415],
       [ 0.88045275,  0.81275467,  0.67921008],
       [ 0.87584238,  0.80373407,  0.66485739],
       [ 0.87133536,  0.79471906,  0.65052474],
       [ 0.86692842,  0.78570823,  0.63621923],
       [ 0.86262016,  0.77669961,  0.62194447],
       [ 0.85840907,  0.76769125,  0.60770418],
       [ 0.85429654,  0.75868029,  0.59349708],
       [ 0.85028007,  0.74966496,  0.57932868],
       [ 0.84635633,  0.74064381,  0.56520582],
       [ 0.84252327,  0.73161496,  0.55113296],
       [ 0.8387787 ,  0.72257655,  0.5371147 ],
       [ 0.83512518,  0.71352502,  0.52314782],
       [ 0.83155579,  0.70446008,  0.50924459],
       [ 0.82806728,  0.69538019,  0.49541101],
       [ 0.82465652,  0.6862838 ,  0.48165268],
       [ 0.82132018,  0.67716938,  0.46797534],
       [ 0.81805889,  0.66803387,  0.4543786 ],
       [ 0.81486489,  0.65887746,  0.44087471],
       [ 0.81173319,  0.64969936,  0.42747133],
       [ 0.80865922,  0.64049863,  0.41417518],
       [ 0.80563815,  0.63127449,  0.40099313],
       [ 0.80266611,  0.6220258 ,  0.38793048],
       [ 0.79973819,  0.61275187,  0.37499375],
       [ 0.7968461 ,  0.60345366,  0.36219373],
       [ 0.79398387,  0.59413121,  0.34953764],
       [ 0.79114533,  0.58478478,  0.3370326 ],
       [ 0.78832414,  0.5754148 ,  0.32468562],
       [ 0.78551416,  0.56602173,  0.31250312],
       [ 0.78270818,  0.55660672,  0.30049247],
       [ 0.7798987 ,  0.54717125,  0.28866097],
       [ 0.77707909,  0.53771654,  0.27701448],
       [ 0.77424278,  0.52824392,  0.2655584 ],
       [ 0.77138332,  0.51875486,  0.25429769],
       [ 0.76849381,  0.50925122,  0.2432374 ],
       [ 0.76556663,  0.49973557,  0.23238296],
       [ 0.76259701,  0.49020893,  0.22173642],
       [ 0.75957943,  0.48067294,  0.21130041],
       [ 0.75650866,  0.47112921,  0.20107703],
       [ 0.75337976,  0.4615793 ,  0.19106786],
       [ 0.7501863 ,  0.45202597,  0.18127558],
       [ 0.7469232 ,  0.44247129,  0.17170133],
       [ 0.74358878,  0.43291502,  0.16234321],
       [ 0.74017957,  0.42335827,  0.15320098],
       [ 0.73669234,  0.41380198,  0.14427416],
       [ 0.7331239 ,  0.40424718,  0.13556224],
       [ 0.72946771,  0.39469765,  0.12706684],
       [ 0.72572584,  0.38515042,  0.11878401],
       [ 0.72189621,  0.37560569,  0.11071302],
       [ 0.71797691,  0.3660635 ,  0.1028534 ],
       [ 0.71396596,  0.35652388,  0.09520513],
       [ 0.70985937,  0.34698868,  0.08776933],
       [ 0.70565928,  0.33745415,  0.08054562],
       [ 0.70136433,  0.32791943,  0.07353618],
       [ 0.69697318,  0.31838347,  0.06674438],
       [ 0.69248439,  0.30884517,  0.06017512],
       [ 0.68789757,  0.29930215,  0.05383523],
       [ 0.68321235,  0.28975156,  0.04773433],
       [ 0.67842749,  0.28019121,  0.04188515],
       [ 0.67354163,  0.27061862,  0.03636072],
       [ 0.66855512,  0.26102885,  0.03150742],
       [ 0.66346759,  0.2514171 ,  0.02731157],
       [ 0.65827645,  0.24178064,  0.02372197],
       [ 0.65297999,  0.23211508,  0.02069175],
       [ 0.64757633,  0.22241545,  0.01817765],
       [ 0.64206692,  0.21267093,  0.0161447 ],
       [ 0.63644967,  0.20287411,  0.0145582 ],
       [ 0.63071962,  0.1930207 ,  0.01338273],
       [ 0.62487374,  0.18310217,  0.01258896],
       [ 0.61890858,  0.17310874,  0.01215089],
       [ 0.61282014,  0.16302927,  0.01204588],
       [ 0.60660685,  0.15284485,  0.01226209],
       [ 0.60026273,  0.14254078,  0.01278509],
       [ 0.59377776,  0.13210708,  0.01359628],
       [ 0.58714329,  0.12152687,  0.01468699],
       [ 0.58034867,  0.11078224,  0.01605237],
       [ 0.57338073,  0.0998553 ,  0.01769137],
       [ 0.56622295,  0.08873018,  0.01960642],
       [ 0.5588547 ,  0.07739692,  0.02180245],
       [ 0.55125032,  0.06585838,  0.02428496],
       [ 0.54337852,  0.05414245,  0.02705628],
       [ 0.5352025 ,  0.04232293,  0.03010918],
       [ 0.52668183,  0.03110507,  0.0334188 ],
       [ 0.51777611,  0.02205008,  0.03693986],
       [ 0.50845819,  0.01523191,  0.04052421],
       [ 0.49872523,  0.01054831,  0.0438642 ],
       [ 0.48860559,  0.00776376,  0.04678873],
       [ 0.47815742,  0.00652746,  0.04918358],
       [ 0.4674569 ,  0.00643321,  0.05098008],
       [ 0.45658138,  0.0070938 ,  0.05216227],
       [ 0.44559546,  0.00819617,  0.05275893],
       [ 0.43454446,  0.00952137,  0.05282708],
       [ 0.42345487,  0.01093554,  0.05243525],
       [ 0.41233847,  0.01236758,  0.05165127],
       [ 0.40119796,  0.01378391,  0.0505346 ],
       [ 0.39014363,  0.0148316 ,  0.04899343],
       [ 0.37907097,  0.01582809,  0.04721344],
       [ 0.36801176,  0.01666972,  0.04518377],
       [ 0.35702727,  0.01718545,  0.04286751],
       [ 0.34600304,  0.01769357,  0.04040359],
       [ 0.3351019 ,  0.01776285,  0.03766479],
       [ 0.32413664,  0.01788554,  0.03494397],
       [ 0.3132996 ,  0.01758563,  0.03211375],
       [ 0.30240405,  0.01732319,  0.02934371],
       [ 0.29160425,  0.01674398,  0.02655319],
       [ 0.28077856,  0.0161301 ,  0.02384022],
       [ 0.26998562,  0.0153579 ,  0.02118968],
       [ 0.25922123,  0.01444636,  0.01862347],
       [ 0.24840331,  0.01355671,  0.01618812]]

RdYlCyBu = LinearSegmentedColormap.from_list('RdYlCyBu', cm_data)
RdYlCyBu_r = LinearSegmentedColormap.from_list('RdYlCyBu_r', cm_data[::-1])

# Make a clean "from <colormap> import *"
__all__ = ['RdYlCyBu', 'RdYlCyBu_r']

# For use with viscm
test_cm = RdYlCyBu

#
