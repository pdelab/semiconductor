Point(1) = {2, 1, 1, 1.0};
Point(2) = {2, 1, -1, 1.0};
Point(3) = {2, -1, -1, 1.0};
Point(4) = {2, -1, 1, 1.0};
Point(5) = {2, -1, 1, 1.0};
Point(6) = {.25, 1, 1, 1.0};
Point(7) = {.25, -1, 1, 1.0};
Point(8) = {-.25, -1, 1, 1.0};
Point(9) = {-.25, 1, 1, 1.0};
Point(10) = {-2, 1, 1, 1.0};
Point(11) = {-2, -1, 1, 1.0};
Point(12) = {-2, -1, -1, 1.0};
Point(13) = {-2, 1, -1, 1.0};
Point(14) = {-1, 1, -1, 1.0};
Point(15) = {-1, -1, -1, 1.0};
Point(16) = {1, -1, -1, 1.0};
Point(17) = {1, 1, -1, 1.0};Point(18) = {0, 1, 1, 1.0};
Point(19) = {0, -1, 1, 1.0};
Point(20) = {0, -1, -1, 1.0};
Point(21) = {0, 1, -1, 1.0};
Line(1) = {10, 13};
Line(2) = {13, 12};
Line(3) = {12, 11};
Line(4) = {11, 10};
Line(5) = {1, 4};
Line(6) = {4, 3};
Line(7) = {3, 2};
Line(8) = {2, 1};
Line(9) = {1, 6};
Line(10) = {6, 18};
Line(11) = {18, 9};
Line(12) = {9, 10};
Line(13) = {13, 14};
Line(14) = {14, 21};
Line(15) = {21, 17};
Line(16) = {17, 2};
Line(17) = {12, 15};
Line(18) = {15, 20};
Line(19) = {20, 16};
Line(20) = {16, 3};
Line(21) = {14, 15};
Line(22) = {21, 20};
Line(23) = {17, 16};
Line(24) = {9, 8};
Line(25) = {8, 19};
Line(26) = {19, 7};
Line(27) = {7, 6};
Line(28) = {18, 19};
Line(29) = {8, 11};
Line(30) = {7, 4};
Line(31) = {18, 21};
Line(32) = {19, 20};
Line Loop(33) = {1, 2, 3, 4};
Plane Surface(34) = {33};
Line Loop(35) = {5, 6, 7, 8};
Plane Surface(36) = {35};
Line Loop(37) = {11, 24, 25, -28};
Plane Surface(38) = {37};
Line Loop(39) = {28, 26, 27, 10};
Plane Surface(40) = {39};
Line Loop(41) = {28, 32, -22, -31};
Plane Surface(42) = {41};
Line Loop(43) = {14, 22, -18, -21};
Plane Surface(44) = {43};
Line Loop(45) = {22, 19, -23, -15};
Plane Surface(46) = {45};
Line Loop(47) = {12, 1, 13, 14, -31, 11};
Plane Surface(48) = {47};
Line Loop(49) = {31, 15, 16, 8, 9, 10};
Plane Surface(50) = {49};
Line Loop(51) = {13, 21, -17, -2};
Plane Surface(52) = {51};
Line Loop(53) = {16, -7, -20, -23};
Plane Surface(54) = {53};
Line Loop(55) = {12, -4, -29, -24};
Plane Surface(56) = {55};
Line Loop(57) = {27, -9, 5, -30};
Plane Surface(58) = {57};
Line Loop(59) = {29, -3, 17, 18, -32, -25};
Plane Surface(60) = {59};
Line Loop(61) = {26, 30, 6, -20, -19, -32};
Plane Surface(62) = {61};
Surface Loop(67) = {56, 48, 34, 52, 44, 60, 38, 42};
Surface Loop(69) = {58, 40, 62, 36, 54, 50, 46, 42};
Volume(68) = {67};
Volume(70) = {69};

Physical Surface(3) = {34};
Physical Surface(4) = {36};
Physical Surface(1) = {38, 40};
Physical Surface(2) = {44, 46};

Physical Volume(1) = {68};
Physical Volume(2) = {70};