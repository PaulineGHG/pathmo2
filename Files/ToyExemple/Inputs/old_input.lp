%*
CHEMICALS
====================================================================================================
*%

% M1
atom("m1", c, 1..10).
atom("m1", o, 11..13).
bond("m1", 11, 1, simple).
bond("m1", 1, 12, double).
bond("m1", 1, 2, simple).
bond("m1", 2, 3, simple).
bond("m1", 3, 4, simple).
bond("m1", 4, 5, simple).
bond("m1", 5, 6, simple).
bond("m1", 5, 13, simple).
bond("m1", 6, 13, simple).
bond("m1", 6, 7, simple).
bond("m1", 7, 8, simple).
bond("m1", 8, 9, simple).
bond("m1", 9, 10, simple).

% M2
atom("m2", c, 1..10).
atom("m2", o, 11..13).
atom("m2", s, 14).
atom("m2", c, 15).
bond("m2", 11, 1, simple).
bond("m2", 1, 12, double).
bond("m2", 1, 2, simple).
bond("m2", 2, 3, simple).
bond("m2", 3, 4, simple).
bond("m2", 4, 5, simple).
bond("m2", 5, 6, simple).
bond("m2", 5, 13, simple).
bond("m2", 6, 7, simple).
bond("m2", 7, 8, simple).
bond("m2", 8, 9, simple).
bond("m2", 9, 10, simple).
bond("m2", 6, 14, simple).
bond("m2", 14, 15, simple).

% M5
atom("m5", c, 1..10).
atom("m5", o, 11..13).
bond("m5", 11, 1, simple).
bond("m5", 1, 12, double).
bond("m5", 1, 2, simple).
bond("m5", 2, 3, simple).
bond("m5", 3, 4, simple).
bond("m5", 4, 5, simple).
bond("m5", 5, 6, double).
bond("m5", 6, 7, simple).
bond("m5", 7, 8, simple).
bond("m5", 7, 13, simple).
bond("m5", 8, 13, simple).
bond("m5", 8, 9, simple).
bond("m5", 9, 10, simple).


%*
REACTIONS
====================================================================================================
*%

% R1
reaction(r1,"m1","m2").