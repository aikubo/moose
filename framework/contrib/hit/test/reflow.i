[Block]
[./right]
type = TestAlignment
xy_data = "   0 0.0
1 1.0
10 3.0
100 2.0
1000 4.0"
[../]
[./left]
type = TestAlignment
xy_data = "100 1.0
200 12.0
300 313.0
400 4514.0
500 45615.0"
[../]
[./no_reflow]
type = TestReflow
long_vector = '10 1 20 2 30 3 40 4 50 5 60 6 70 7 80 8 90 9 110 11 120 12 130 13 140 14 150 15 160 16 170 17 180 18 190 19 210 21 220 22 230 23 240 24 250 25 260 26 270 27 280 28 290 29'
[../]
[./reflow]
type = TestReflow
long_vector = "10 1 20 2 30 3 40 4 50 5 60 6 70 7 80 8 90 9 110 11 120 12 130 13 140 14 150 15 160 16 170 17 180 18 190 19 210 21 220 22 230 23 240 24 250 25 260 26 270 27 280 28 290 29"
[../][]
