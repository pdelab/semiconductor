def channel_generation():
    from numpy import sqrt
    import pandas as pd

    filename  = "gramA_channel"
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
    fo.write(bytes("""Line(1) = {7, 8};
    Line(2) = {8, 9};
    Line(3) = {9, 6};
    Line(4) = {6, 7};
    Line(5) = {7, 30};
    Line(6) = {30, 31};
    Line(7) = {31, 32};
    Line(8) = {32, 29};
    Line(9) = {29, 30};
    Line(10) = {8, 29};
    Line(11) = {9, 32};
    Line(12) = {6, 31};
    Line(13) = {27, 26};
    Line(14) = {26, 1};
    Line(15) = {1, 2};
    Line(16) = {2, 27};
    Line(17) = {27, 28};
    Line(18) = {28, 3};
    Line(19) = {3, 2};
    Line(20) = {1, 4};
    Line(21) = {4, 3};
    Line(22) = {28, 33};
    Line(23) = {33, 26};
    Line(24) = {33, 4};
    Line(25) = {32, 33};
    Line(26) = {29, 28};
    Line(27) = {31, 26};
    Line(28) = {27, 30};
    Circle(29) = {12, 45, 13};
    Circle(30) = {13, 45, 10};
    Circle(31) = {10, 45, 11};
    Circle(32) = {11, 45, 12};
    Circle(33) = {41, 45, 38};
    Circle(34) = {38, 45, 39};
    Circle(35) = {39, 45, 40};
    Circle(36) = {40, 45, 41};
    Circle(37) = {59, 60, 56};
    Circle(38) = {56, 60, 57};
    Circle(39) = {57, 60, 58};
    Circle(40) = {58, 60, 59};
    Circle(41) = {25, 44, 22};
    Circle(42) = {22, 44, 23};
    Circle(43) = {23, 44, 24};
    Circle(44) = {24, 44, 25};
    Circle(45) = {48, 50, 49};
    Circle(46) = {49, 50, 46};
    Circle(47) = {46, 50, 47};
    Circle(48) = {47, 50, 48};
    Circle(49) = {18, 43, 21};
    Circle(50) = {21, 43, 20};
    Circle(51) = {20, 43, 19};
    Circle(52) = {19, 43, 18};
    Circle(53) = {53, 55, 54};
    Circle(54) = {54, 55, 51};
    Circle(55) = {51, 55, 52};
    Circle(56) = {52, 55, 53};
    Circle(57) = {17, 42, 14};
    Circle(58) = {14, 42, 15};
    Circle(59) = {15, 42, 16};
    Circle(60) = {16, 42, 17};
    Circle(61) = {36, 42, 37};
    Circle(62) = {37, 42, 34};
    Circle(63) = {34, 42, 35};
    Circle(64) = {35, 42, 36};
    Line(65) = {12, 25};
    Line(66) = {41, 59};
    Line(67) = {13, 22};
    Line(68) = {38, 56};
    Line(69) = {11, 24};
    Line(70) = {40, 58};
    Line(71) = {10, 23};
    Line(72) = {39, 57};
    Line(73) = {59, 48};
    Line(74) = {56, 49};
    Line(75) = {58, 47};
    Line(76) = {57, 46};
    Line(77) = {48, 53};
    Line(78) = {49, 54};
    Line(79) = {47, 52};
    Line(80) = {46, 51};
    Line(81) = {25, 18};
    Line(82) = {22, 21};
    Line(83) = {24, 19};
    Line(84) = {23, 20};
    Line(85) = {18, 17};
    Line(86) = {53, 36};
    Line(87) = {21, 14};
    Line(88) = {54, 37};
    Line(89) = {19, 16};
    Line(90) = {52, 35};
    Line(91) = {20, 15};
    Line(92) = {51, 34};
    Line Loop(93) = {45, 46, 47, 48};
    Plane Surface(94) = {93};
    Line Loop(95) = {37, 38, 39, 40};
    Plane Surface(96) = {95};
    Line Loop(97) = {54, 55, 56, 53};
    Plane Surface(98) = {97};
    Line Loop(99) = {33, 34, 35, 36};
    Plane Surface(100) = {99};
    Line Loop(101) = {61, 62, 63, 64};
    Plane Surface(102) = {101};
    Line Loop(103) = {29, 30, 31, 32};
    Plane Surface(104) = {99, 103};
    Line Loop(105) = {57, 58, 59, 60};
    Plane Surface(106) = {101, 105};
    Line Loop(107) = {6, 7, 8, 9};
    Line Loop(108) = {41, 42, 43, 44};
    Plane Surface(109) = {107, 108};
    Line Loop(110) = {13, -23, -22, -17};
    Line Loop(111) = {49, 50, 51, 52};
    Plane Surface(112) = {110, 111};
    Line Loop(113) = {65, 41, -67, -29};
    Ruled Surface(114) = {113};
    Line Loop(115) = {65, -44, -69, 32};
    Ruled Surface(116) = {115};
    Line Loop(117) = {69, -43, -71, 31};
    Ruled Surface(118) = {117};
    Line Loop(119) = {30, 71, -42, -67};
    Ruled Surface(120) = {119};
    Line Loop(121) = {49, 87, -57, -85};
    Ruled Surface(122) = {121};
    Line Loop(123) = {60, -85, -52, 89};
    Ruled Surface(124) = {123};
    Line Loop(125) = {89, -59, -91, 51};
    Ruled Surface(126) = {125};
    Line Loop(127) = {91, -58, -87, 50};
    Ruled Surface(128) = {127};
    Line Loop(129) = {86, 61, -88, -53};
    Ruled Surface(130) = {129};
    Line Loop(131) = {86, -64, -90, 56};
    Ruled Surface(132) = {131};
    Line Loop(133) = {90, -63, -92, 55};
    Ruled Surface(134) = {133};
    Line Loop(135) = {92, -62, -88, 54};
    Ruled Surface(136) = {135};
    Line Loop(137) = {53, -78, -45, 77};
    Ruled Surface(138) = {137};
    Line Loop(139) = {78, 54, -80, -46};
    Ruled Surface(140) = {139};
    Line Loop(141) = {80, 55, -79, -47};
    Ruled Surface(142) = {141};
    Line Loop(143) = {79, 56, -77, -48};
    Ruled Surface(144) = {143};
    Line Loop(145) = {48, -73, -40, 75};
    Ruled Surface(146) = {145};
    Line Loop(147) = {75, -47, -76, 39};
    Ruled Surface(148) = {147};
    Line Loop(149) = {76, -46, -74, 38};
    Ruled Surface(150) = {149};
    Line Loop(151) = {74, -45, -73, 37};
    Ruled Surface(152) = {151};
    Line Loop(153) = {39, -70, -35, 72};
    Ruled Surface(154) = {153};
    Line Loop(155) = {70, 40, -66, -36};
    Ruled Surface(156) = {155};
    Line Loop(157) = {66, 37, -68, -33};
    Ruled Surface(158) = {157};
    Line Loop(159) = {38, -72, -34, 68};
    Ruled Surface(160) = {159};
    Line Loop(161) = {28, -9, 26, -17};
    Plane Surface(162) = {161};
    Line Loop(163) = {28, 6, 27, -13};
    Plane Surface(164) = {163};
    Line Loop(165) = {26, 22, -25, 8};
    Plane Surface(166) = {165};
    Line Loop(167) = {25, 23, -27, 7};
    Plane Surface(168) = {167};
    Line Loop(169) = {41, 82, -49, -81};
    Ruled Surface(170) = {169};
    Line Loop(171) = {81, -52, -83, 44};
    Ruled Surface(172) = {171};
    Line Loop(173) = {83, -51, -84, 43};
    Ruled Surface(174) = {173};
    Line Loop(175) = {50, -84, -42, 82};
    Ruled Surface(176) = {175};
    Line Loop(177) = {5, 6, -12, 4};
    Plane Surface(178) = {177};
    Line Loop(179) = {5, -9, -10, -1};
    Plane Surface(180) = {179};
    Line Loop(181) = {4, 1, 2, 3};
    Plane Surface(182) = {181};
    Line Loop(183) = {3, 12, 7, -11};
    Plane Surface(184) = {183};
    Line Loop(185) = {11, 8, -10, 2};
    Plane Surface(186) = {185};
    Line Loop(187) = {16, 17, 18, 19};
    Plane Surface(188) = {187};
    Line Loop(189) = {16, 13, 14, 15};
    Plane Surface(190) = {189};
    Line Loop(191) = {15, -19, -21, -20};
    Plane Surface(192) = {191};
    Line Loop(193) = {18, -21, -24, -22};
    Plane Surface(194) = {193};
    Line Loop(195) = {24, -20, -14, -23};
    Plane Surface(196) = {195};
    
    Physical Surface(1) = {164, 162, 166, 168};
    Physical Surface(2) = {178, 180, 186, 184};
    Physical Surface(3) = {190, 188, 194, 196};
    Physical Surface(4) = {182};
    Physical Surface(5) = {192};
    Physical Surface(6) = {109, 112};
    Physical Surface(7) = {114, 116, 118, 120, 104, 122, 124, 126, 128, 106, 158, 152, 138, 130, 132, 144, 146, 156, 154, 148, 142, 134, 136, 140, 150, 160};
    Physical Surface(8) = {94};
    Physical Surface(9) = {96};
    Physical Surface(10) = {100};
    Physical Surface(11) = {98};
    Physical Surface(12) = {102};
    
    Surface Loop(197) = {114, 116, 118, 120, 104, 158, 156, 154, 160, 148, 146, 152, 150, 142, 140, 138, 144, 134, 132, 130, 136, 106, 124, 122, 128, 126, 170, 176, 174, 172};
    Volume(198) = {197};
    Surface Loop(199) = {164, 162, 166, 168, 109, 170, 176, 174, 172, 112};
    Volume(200) = {199};
    Surface Loop(201) = {178, 180, 186, 184, 182, 109, 114, 116, 118, 120, 104, 100};
    Volume(202) = {201};
    Surface Loop(203) = {190, 188, 194, 192, 196, 112, 122, 128, 126, 124, 106, 102};
    Volume(204) = {203};
    Surface Loop(205) = {158, 156, 154, 160, 100, 96};
    Volume(206) = {205};
    Surface Loop(207) = {96, 152, 150, 148, 146, 94};
    Volume(208) = {207};
    Surface Loop(209) = {138, 140, 142, 144, 94, 98};
    Volume(210) = {209};
    Surface Loop(211) = {98, 130, 132, 134, 136, 102};
    Volume(212) = {211};
    
    Physical Volume(1) = {198};
    Physical Volume(2) = {200};
    Physical Volume(3) = {202,204,206,208,210,212};""",'UTF-8'))
if __name__=="__main__":
    channel_generation()
