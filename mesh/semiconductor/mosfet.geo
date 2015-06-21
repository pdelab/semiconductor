// DIMENSIONS
L      = 1.6;	// device length
W      = 0.2;   // device width
H      = 1.4;	// device height
emit_H = 0.3;	// emitter hieght/thickness
emit_L = 0.5;	// emitter length
coll_H = 0.3;	// collector hieght/thickness
coll_L = 0.5;	// collector length







// POINTS

// Box
Point(1) = { L/2, W/2, 0, 1.0};
Point(2) = {-L/2, W/2, 0, 1.0};
Point(3) = {-L/2,-W/2, 0, 1.0};
Point(4) = { L/2,-W/2, 0, 1.0};
Point(5) = { L/2, W/2, H, 1.0};
Point(6) = {-L/2, W/2, H, 1.0};
Point(7) = {-L/2,-W/2, H, 1.0};
Point(8) = { L/2,-W/2, H, 1.0};

// Emitter
Point(9)  = {-L/2, W/2, H-emit_H, 1.0};
Point(10) = {-L/2,-W/2, H-emit_H, 1.0};
Point(11) = {-L/2+emit_L-emit_H, W/2, H-emit_H, 1.0};
Point(12) = {-L/2+emit_L-emit_H,-W/2, H-emit_H, 1.0};
Point(13) = {-L/2+emit_L-emit_H, W/2, H, 1.0};
Point(14) = {-L/2+emit_L-emit_H,-W/2, H, 1.0};
Point(15) = {-L/2+emit_L, W/2, H, 1.0};
Point(16) = {-L/2+emit_L,-W/2, H, 1.0};

// Collector
Point(17) = { L/2, W/2, H-coll_H, 1.0};
Point(18) = { L/2,-W/2, H-coll_H, 1.0};
Point(19) = { L/2-coll_L+coll_H, W/2, H-coll_H, 1.0};
Point(20) = { L/2-coll_L+coll_H,-W/2, H-coll_H, 1.0};
Point(21) = { L/2-coll_L+coll_H, W/2, H, 1.0};
Point(22) = { L/2-coll_L+coll_H,-W/2, H, 1.0};
Point(23) = { L/2-coll_L, W/2, H, 1.0};
Point(24) = { L/2-coll_L,-W/2, H, 1.0};



// LINES
Line(1) = {6, 7};
Line(2) = {15, 16};
Line(3) = {23, 24};
Line(4) = {5, 8};
Line(5) = {13, 14};
Line(6) = {21, 22};
Line(7) = {9, 10};
Line(8) = {11, 12};
Line(9) = {19, 20};
Line(10) = {17, 18};
Line(11) = {6, 9};
Line(12) = {7, 10};
Line(13) = {13, 11};
Line(14) = {14, 12};
Line(15) = {7, 14};
Line(16) = {6, 13};
Line(17) = {13, 15};
Line(18) = {14, 16};
Line(19) = {10, 12};
Line(20) = {9, 11};
Line(21) = {23, 21};
Line(22) = {21, 5};
Line(23) = {21, 19};
Line(24) = {19, 17};
Line(25) = {17, 5};
Line(26) = {24, 22};
Line(27) = {22, 8};
Line(28) = {22, 20};
Line(29) = {20, 18};
Line(30) = {18, 8};
Line(31) = {15, 23};
Line(32) = {16, 24};
Circle(33) = {11, 13, 15};
Circle(34) = {12, 14, 16};
Circle(35) = {19, 21, 23};
Circle(36) = {20, 22, 24};
Line(37) = {9, 2};
Line(38) = {17, 1};
Line(39) = {18, 4};
Line(40) = {10, 3};
Line(41) = {3, 2};
Line(42) = {2, 1};
Line(43) = {1, 4};
Line(44) = {4, 3};





// SURFACES
Line Loop(45) = {16, 5, -15, -1};
Plane Surface(46) = {45};
Line Loop(47) = {5, 18, -2, -17};
Plane Surface(48) = {47};
Line Loop(49) = {16, 13, -20, -11};
Plane Surface(50) = {49};
Line Loop(51) = {13, 33, -17};
Plane Surface(52) = {51};
Line Loop(53) = {1, 12, -7, -11};
Plane Surface(54) = {53};
Line Loop(55) = {12, 19, -14, -15};
Plane Surface(56) = {55};
Line Loop(57) = {19, -8, -20, 7};
Plane Surface(58) = {57};
Line Loop(59) = {8, -14, -5, 13};
Plane Surface(60) = {59};
Line Loop(61) = {18, -34, -14};
Plane Surface(62) = {61};
Line Loop(63) = {21, 23, 35};
Plane Surface(64) = {63};
Line Loop(65) = {23, 24, 25, -22};
Plane Surface(66) = {65};
Line Loop(67) = {21, 6, -26, -3};
Plane Surface(68) = {67};
Line Loop(69) = {6, 28, -9, -23};
Plane Surface(70) = {69};
Line Loop(71) = {22, 4, -27, -6};
Plane Surface(72) = {71};
Line Loop(73) = {4, -30, -10, 25};
Plane Surface(74) = {73};
Line Loop(75) = {24, 10, -29, -9};
Plane Surface(76) = {75};
Line Loop(77) = {29, 30, -27, 28};
Plane Surface(78) = {77};
Line Loop(79) = {28, 36, 26};
Plane Surface(80) = {79};
Line Loop(81) = {32, -3, -31, 2};
Plane Surface(82) = {81};
Line Loop(83) = {34, -2, -33, 8};
Ruled Surface(84) = {83};
Line Loop(85) = {36, -3, -35, 9};
Ruled Surface(86) = {85};
Line Loop(87) = {40, 41, -37, 7};
Plane Surface(88) = {87};
Line Loop(89) = {37, 42, -38, -24, 35, -31, -33, -20};
Plane Surface(90) = {89};
Line Loop(91) = {39, -43, -38, 10};
Plane Surface(92) = {91};
Line Loop(93) = {43, 44, 41, 42};
Plane Surface(94) = {93};
Line Loop(95) = {44, -40, 19, 34, 32, -36, 29, 39};
Plane Surface(96) = {95};



// VOLUMES
Surface Loop(97) = {46, 50, 54, 56, 58, 84, 48, 62, 52};
Volume(98) = {97};
Surface Loop(99) = {68, 64, 80, 86, 76, 74, 72, 66, 78};
Volume(100) = {99};
Surface Loop(101) = {82, 96, 94, 92, 90, 88, 86, 76, 84, 58};
Volume(102) = {101};


// PHYSICAL ENTITIES
//Physical Surface(1) = {46, 48};
//Physical Surface(2) = {68, 72};
Physical Surface(1) = {46};
Physical Surface(2) = {72};
Physical Surface(3) = {82};
Physical Volume(1)  = {98};
Physical Volume(2)  = {100};
Physical Volume(3)  = {102};
