def channel_generation():
    from numpy import sqrt
    import pandas as pd

    filename  = "short_channel"

    parameter = pd.read_csv("./channel_parameters.csv")
    box_length = float(parameter.box_length.values[0]) 
    box_width = float(parameter.box_width.values[0])
    box_height = float(parameter.box_height.values[0])
    channel_radius = float(parameter.channel_radius.values[0])
    channel_thickness = float(parameter.channel_thickness.values[0])
    channel_length = float(parameter.channel_length.values[0])
    membrane_thickness = float(parameter.membrane_thickness.values[0])
    Point = list()
    Point.append((0.0,0.0,0.0,0.0))
    Point.append((box_length/2, box_width/2, box_height/2, 1))
    Point.append((box_length/2, box_width/2, -box_height/2, 1))
    Point.append((box_length/2, -box_width/2, -box_height/2, 1))
    Point.append((box_length/2, -box_width/2, box_height/2, 1))
    Point.append((box_length/2, box_width/2, box_height/2, 1))
    Point.append((-box_length/2, box_width/2, box_height/2, 1))
    Point.append((-box_length/2, box_width/2, -box_height/2, 1))
    Point.append((-box_length/2, -box_width/2, -box_height/2, 1))
    Point.append((-box_length/2, -box_width/2, box_height/2, 1))
    Point.append((-channel_length/2, -(channel_thickness+channel_radius)*(1/sqrt(2)), (channel_thickness+channel_radius)*(1/sqrt(2)), 1))
    Point.append((-channel_length/2, -(channel_thickness+channel_radius)*(1/sqrt(2)), -(channel_thickness+channel_radius)*(1/sqrt(2)), 1))
    Point.append((-channel_length/2, (channel_thickness+channel_radius)*(1/sqrt(2)), -(channel_thickness+channel_radius)*(1/sqrt(2)), 1))
    Point.append((-channel_length/2, (channel_thickness+channel_radius)*(1/sqrt(2)), (channel_thickness+channel_radius)*(1/sqrt(2)), 1))
    Point.append((channel_length/2, (channel_thickness+channel_radius)*(1/sqrt(2)), (channel_thickness+channel_radius)*(1/sqrt(2)), 1))
    Point.append((channel_length/2, -(channel_thickness+channel_radius)*(1/sqrt(2)), (channel_thickness+channel_radius)*(1/sqrt(2)), 1))
    Point.append((channel_length/2, -(channel_thickness+channel_radius)*(1/sqrt(2)), -(channel_thickness+channel_radius)*(1/sqrt(2)), 1))
    Point.append((channel_length/2, (channel_thickness+channel_radius)*(1/sqrt(2)), -(channel_thickness+channel_radius)*(1/sqrt(2)), 1))
    Point.append((membrane_thickness/2, (channel_thickness+channel_radius)*(1/sqrt(2)), -(channel_thickness+channel_radius)*(1/sqrt(2)), 1))
    Point.append((membrane_thickness/2, -(channel_thickness+channel_radius)*(1/sqrt(2)), -(channel_thickness+channel_radius)*(1/sqrt(2)), 1))
    Point.append((membrane_thickness/2, -(channel_thickness+channel_radius)*(1/sqrt(2)), (channel_thickness+channel_radius)*(1/sqrt(2)), 1))
    Point.append((membrane_thickness/2, (channel_thickness+channel_radius)*(1/sqrt(2)), (channel_thickness+channel_radius)*(1/sqrt(2)), 1))
    Point.append((-membrane_thickness/2, (channel_thickness+channel_radius)*(1/sqrt(2)), (channel_thickness+channel_radius)*(1/sqrt(2)), 1))
    Point.append((-membrane_thickness/2, -(channel_thickness+channel_radius)*(1/sqrt(2)), (channel_thickness+channel_radius)*(1/sqrt(2)), 1))
    Point.append((-membrane_thickness/2, -(channel_thickness+channel_radius)*(1/sqrt(2)), -(channel_thickness+channel_radius)*(1/sqrt(2)), 1))
    Point.append((-membrane_thickness/2, (channel_thickness+channel_radius)*(1/sqrt(2)), -(channel_thickness+channel_radius)*(1/sqrt(2)), 1))
    Point.append((membrane_thickness/2, box_width/2, box_height/2, 1))
    Point.append((membrane_thickness/2, box_width/2, -box_height/2, 1))
    Point.append((membrane_thickness/2, -box_width/2, -box_height/2, 1))
    Point.append((-membrane_thickness/2, -box_width/2, -box_height/2, 1))
    Point.append((-membrane_thickness/2, box_width/2, -box_height/2, 1))
    Point.append((-membrane_thickness/2, box_width/2, box_height/2, 1))
    Point.append((-membrane_thickness/2, -box_width/2, box_height/2, 1))
    Point.append((membrane_thickness/2, -box_width/2, box_height/2, 1))
    Point.append((channel_length/2, -(channel_radius)*(1/sqrt(2)), (channel_radius)*(1/sqrt(2)), 1))
    Point.append((channel_length/2, -(channel_radius)*(1/sqrt(2)), -(channel_radius)*(1/sqrt(2)), 1))
    Point.append((channel_length/2, (channel_radius)*(1/sqrt(2)), -(channel_radius)*(1/sqrt(2)), 1))
    Point.append((channel_length/2, (channel_radius)*(1/sqrt(2)), (channel_radius)*(1/sqrt(2)), 1))
    Point.append((-channel_length/2, (channel_radius)*(1/sqrt(2)), (channel_radius)*(1/sqrt(2)), 1))
    Point.append((-channel_length/2, -(channel_radius)*(1/sqrt(2)), (channel_radius)*(1/sqrt(2)), 1))
    Point.append((-channel_length/2, -(channel_radius)*(1/sqrt(2)), -(channel_radius)*(1/sqrt(2)), 1))
    Point.append((-channel_length/2, (channel_radius)*(1/sqrt(2)), -(channel_radius)*(1/sqrt(2)), 1))
    Point.append((channel_length/2, 0, 0, 1))
    Point.append((membrane_thickness/2, 0, 0, 1))
    Point.append((-membrane_thickness/2, 0, 0, 1))
    Point.append((-channel_length/2, 0, 0, 1))
    
    # current cross-sections
    Point.append((0, -(channel_radius)*(1/sqrt(2)), (channel_radius)*(1/sqrt(2)), 1))
    Point.append((0, -(channel_radius)*(1/sqrt(2)), -(channel_radius)*(1/sqrt(2)), 1))
    Point.append((0, (channel_radius)*(1/sqrt(2)), -(channel_radius)*(1/sqrt(2)), 1))
    Point.append((0, (channel_radius)*(1/sqrt(2)), (channel_radius)*(1/sqrt(2)), 1))
    Point.append((0, 0, 0, 1))
    Point.append((channel_length/4, -(channel_radius)*(1/sqrt(2)), (channel_radius)*(1/sqrt(2)), 1))
    Point.append((channel_length/4, -(channel_radius)*(1/sqrt(2)), -(channel_radius)*(1/sqrt(2)), 1))
    Point.append((channel_length/4, (channel_radius)*(1/sqrt(2)), -(channel_radius)*(1/sqrt(2)), 1))
    Point.append((channel_length/4, (channel_radius)*(1/sqrt(2)), (channel_radius)*(1/sqrt(2)), 1))
    Point.append((channel_length/4, 0, 0, 1))
    Point.append((-channel_length/4, (channel_radius)*(1/sqrt(2)), (channel_radius)*(1/sqrt(2)), 1))
    Point.append((-channel_length/4, -(channel_radius)*(1/sqrt(2)), (channel_radius)*(1/sqrt(2)), 1))
    Point.append((-channel_length/4, -(channel_radius)*(1/sqrt(2)), -(channel_radius)*(1/sqrt(2)), 1))
    Point.append((-channel_length/4, (channel_radius)*(1/sqrt(2)), -(channel_radius)*(1/sqrt(2)), 1))
    Point.append((-channel_length/4, 0, 0, 1))
    
    #write to gmsh file
    fileString = "./{}.geo".format(filename)
    fo = open(fileString, "wb")
    fo.write(bytes("cl__1 = 1;\n", 'UTF-8'))
    for i in range(1,61):
        string = "Point({})".format(i)+" = {"+"{},{},{},1".format(Point[i][0],Point[i][1],Point[i][2])+"};"
        fo.write(bytes(string+'\n','UTF-8'))
    fo.write(bytes("""Line(2) = {1, 2};
Line(4) = {7, 6};
Line(5) = {6, 9};
Line(6) = {9, 8};
Line(7) = {8, 7};
Line(8) = {1, 4};
Line(9) = {4, 3};
Line(10) = {3, 2};
Line(11) = {26, 27};
Line(12) = {31, 26};
Line(13) = {6, 31};
Line(14) = {31, 30};
Line(15) = {30, 7};
Line(16) = {30, 27};
Line(17) = {27, 2};
Line(18) = {26, 1};
Line(19) = {31, 32};
Line(20) = {26, 33};
Line(21) = {32, 33};
Line(22) = {30, 29};
Line(23) = {29, 28};
Line(24) = {28, 27};
Line(25) = {28, 3};
Line(26) = {29, 8};
Circle(27) = {13, 45, 10};
Circle(28) = {10, 45, 11};
Circle(29) = {11, 45, 12};
Circle(30) = {12, 45, 13};
Circle(31) = {38, 45, 39};
Circle(32) = {39, 45, 40};
Circle(33) = {40, 45, 41};
Circle(34) = {41, 45, 38};
Circle(35) = {22, 44, 23};
Circle(36) = {23, 44, 24};
Circle(37) = {24, 44, 25};
Circle(38) = {25, 44, 22};
Circle(39) = {21, 43, 20};
Circle(40) = {20, 43, 19};
Circle(41) = {19, 43, 18};
Circle(42) = {18, 43, 21};
Circle(43) = {14, 42, 15};
Circle(44) = {15, 42, 16};
Circle(45) = {16, 42, 17};
Circle(46) = {17, 42, 14};
Circle(47) = {37, 42, 34};
Circle(48) = {34, 42, 35};
Circle(49) = {35, 42, 36};
Circle(50) = {36, 42, 37};
Line(51) = {13, 22};
Line(52) = {22, 21};
Line(53) = {21, 14};
Line(54) = {10, 23};
Line(55) = {23, 20};
Line(56) = {20, 15};
Line(57) = {37, 38};
Line(58) = {39, 34};
Line(59) = {41, 36};
Line(60) = {40, 35};
Line(61) = {12, 25};
Line(62) = {25, 18};
Line(63) = {18, 17};
Line(64) = {11, 24};
Line(65) = {24, 19};
Line(66) = {19, 16};
Line(107) = {9, 32};
Line(108) = {32, 29};
Line(109) = {28, 33};
Line(110) = {33, 4};
Line Loop(68) = {57, 31, 58, -47};
Ruled Surface(68) = {68};
Line Loop(70) = {34, -57, -50, -59};
Ruled Surface(70) = {70};
Line Loop(72) = {59, -49, -60, 33};
Ruled Surface(72) = {72};
Line Loop(74) = {58, 48, -60, -32};
Ruled Surface(74) = {74};
Line Loop(76) = {51, 35, -54, -27};
Ruled Surface(76) = {76};
Line Loop(78) = {52, 39, -55, -35};
Ruled Surface(78) = {78};
Line Loop(80) = {53, 43, -56, -39};
Ruled Surface(80) = {80};
Line Loop(82) = {51, -38, -61, 30};
Ruled Surface(82) = {82};
Line Loop(84) = {38, 52, -42, -62};
Ruled Surface(84) = {84};
Line Loop(86) = {42, 53, -46, -63};
Ruled Surface(86) = {86};
Line Loop(88) = {61, -37, -64, 29};
Ruled Surface(88) = {88};
Line Loop(90) = {37, 62, -41, -65};
Ruled Surface(90) = {90};
Line Loop(92) = {41, 63, -45, -66};
Ruled Surface(92) = {92};
Line Loop(94) = {54, 36, -64, -28};
Ruled Surface(94) = {94};
Line Loop(96) = {36, 65, -40, -55};
Ruled Surface(96) = {96};
Line Loop(98) = {56, 44, -66, -40};
Ruled Surface(98) = {98};
Line Loop(101) = {30, 27, 28, 29, -34, -33, -32, -31};
Plane Surface(101) = {101};
Line Loop(104) = {45, 46, 43, 44, -49, -48, -47, -50};
Plane Surface(104) = {104};
Line Loop(112) = {16, -24, -23, -22};
Plane Surface(112) = {112};
Line Loop(114) = {16, -11, -12, 14};
Plane Surface(114) = {114};
Line Loop(116) = {23, 109, -21, 108};
Plane Surface(116) = {116};
Line Loop(118) = {19, 108, -22, -14, 19, 108, -22, -14};
Plane Surface(118) = {118};
Line Loop(121) = {14, 22, -108, -19, -37, -36, -35, -38};
Plane Surface(121) = {121};
Line Loop(126) = {15, 4, 13, 14};
Plane Surface(126) = {126};
Line Loop(128) = {15, -7, -26, -22};
Plane Surface(128) = {128};
Line Loop(130) = {17, -10, -25, 24};
Plane Surface(130) = {130};
Line Loop(132) = {17, -2, -18, 11};
Plane Surface(132) = {132};
Line Loop(134) = {2, -10, -9, -8};
Plane Surface(134) = {134};
Line Loop(136) = {25, -9, -110, -109};
Plane Surface(136) = {136};
Line Loop(138) = {26, -6, 107, 108};
Plane Surface(138) = {138};
Line Loop(140) = {4, 5, 6, 7};
Plane Surface(140) = {140};
Line Loop(142) = {12, 20, -21, -19};
Plane Surface(142) = {142};
Line Loop(144) = {18, 8, -110, -20};
Plane Surface(144) = {144};
Line Loop(146) = {107, -19, -13, 5};
Plane Surface(146) = {146};
Line Loop(149) = {11, -24, 109, -20, -41, -40, -39, -42};
Plane Surface(149) = {149};
Surface Loop(151) = {76, 82, 84, 78, 80, 86, 104, 92, 90, 88, 94, 96, 98, 101, 70, 68, 74, 72};
Volume(151) = {151};
Surface Loop(153) = {114, 112, 116, 142, 149, 121, 84, 78, 96, 90};
Volume(153) = {153};
Surface Loop(155) = {146, 138, 128, 126, 140, 121, 82, 76, 94, 88, 101, 70, 68, 74, 104, 92, 86, 80, 98, 72, 149, 136, 130, 132, 134, 144};
Volume(155) = {155};

Circle(200) = {46,50,47};
Circle(201) = {47,50,48};
Circle(202) = {48,50,49};
Circle(203) = {49,50,46};
Line Loop(300) = {200,201,202,203};
Plane Surface(400) = {300};
Circle(204) = {51,55,52};
Circle(205) = {52,55,53};
Circle(206) = {53,55,54};
Circle(207) = {54,55,51};
Line Loop(301) = {204,205,206,207};
Plane Surface(401) = {301};
Circle(208) = {56,60,57};
Circle(209) = {57,60,58};
Circle(210) = {58,60,59};
Circle(211) = {59,60,56};
Line Loop(302) = {208,209,210,211};
Plane Surface(402) = {302};
Line Loop(303) = {31,32,33,34};
Plane Surface(403) = {303};
Line Loop(304) = {47,48,49,50};
Plane Surface(404) = {304};


Physical Surface(1) = {116,142,114,112};
Physical Surface(2) = {132,130,136,144};
Physical Surface(3) = {126,128,138,146};
Physical Surface(4) = {140};
Physical Surface(5) = {134};
Physical Surface(6) = {68, 70, 72, 74, 76, 80, 82, 86, 88, 92, 94, 98, 101, 104, 121, 149};
Physical Surface(7) = {300};
Physical Surface(8) = {301};
Physical Surface(9) = {302};
Physical Surface(10) = {303};
Physical Surface(11) = {304};

Physical Volume(1)  = {151};
Physical Volume(2)  = {153};
Physical Volume(3)  = {155};""",'UTF-8'))
if __name__=="__main__":
    channel_generation()
