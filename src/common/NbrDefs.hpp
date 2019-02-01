#ifndef __NBRDEFS_H__
#define __NBRDEFS_H__

BEGIN_EBI_NAMESPACE

#define NMRK -1
#define IN 1
#define OUT 0
#define CRV 2
#define ON_SURFACE 3

#define POSX 1
#define POSY 2
#define POSZ 3
#define NEGX -1
#define NEGY -2
#define NEGZ -3

#define NORM 0
#define FINE 1
#define CRSE 2

#define WNODE 3
#define XNODE 4

/* Normal Neighbors for Tables */

#define NRM0 0 /* x-, y-, z- */
#define NRM1 1 /* x-, y-, z */
#define NRM2 2 /* x-, y-, z+ */
#define NRM3 3 /* x-, y, z- */
#define NRM4 4 /* x-, y, z */
#define NRM5 5 /* x-, y, z+ */
#define NRM6 6 /* x-, y+, z- */
#define NRM7 7 /* x-, y+, z */
#define NRM8 8 /* x-, y+, z+ */
#define NRM9 9 /* x, y-, z- */
#define NRM10 10 /* x, y-, z */
#define NRM11 11 /* x, y-, z+ */
#define NRM12 12 /* x, y, z- */
#define NRM13 13 /* x, y, z */
#define NRM14 14 /* x, y, z+ */
#define NRM15 15 /* x, y+, z- */
#define NRM16 16 /* x, y+, z */
#define NRM17 17 /* x, y+, z+ */
#define NRM18 18 /* x+, y-, z- */
#define NRM19 19 /* x+, y-, z */
#define NRM20 20 /* x+, y-, z+ */
#define NRM21 21 /* x+, y, z- */
#define NRM22 22 /* x+, y, z */
#define NRM23 23 /* x+, y, z+ */
#define NRM24 24 /* x+, y+, z- */
#define NRM25 25 /* x+, y+, z */
#define NRM26 26 /* x+, y+, z+ */

/* For symmetries */
#define N_CORNER 0
#define N_EDGE 1
#define N_FACE 4
#define N_SELF 13

#define FC_CORNER 0
#define FC_EDGE 1
#define FC_FACE 5
 
#define N_SYM_CORNER 0
#define N_SYM_EDGE 1
#define N_SYM_FACE 2
#define N_SYM_SELF 3
#define F_SYM_CORNER 4
#define F_SYM_EDGE 5
#define F_SYM_FACE 6
#define C_SYM_CORNER 7
#define C_SYM_EDGE 8
#define C_SYM_FACE 9

/* For symmetries for W and X node evaluations */
#define WX_CORNER_0 0
#define WX_CORNER_1 1
#define WX_CORNER_2 7
#define WX_EDGE_0 2
#define WX_EDGE_1 8
#define WX_FACE_0 14

#define WX_SYM_C0 0
#define WX_SYM_C1 1
#define WX_SYM_C2 2
#define WX_SYM_E0 3
#define WX_SYM_E1 4
#define WX_SYM_F0 5


/* Fine Neighbors for Tables */

#define FN0 0
#define FN1 1
#define FN2 2
#define FN3 3
#define FN4 4
#define FN5 5
#define FN6 6
#define FN7 7
#define FN8 8
#define FN9 9
#define FN10 10
#define FN11 11
#define FN12 12
#define FN13 13
#define FN14 14
#define FN15 15
#define FN16 16
#define FN17 17
#define FN18 18
#define FN19 19
#define FN20 20
#define FN21 21
#define FN22 22
#define FN23 23
#define FN24 24
#define FN25 25
#define FN26 26
#define FN27 27
#define FN28 28
#define FN29 29
#define FN30 30
#define FN31 31
#define FN32 32
#define FN33 33
#define FN34 34
#define FN35 35
#define FN36 36
#define FN37 37
#define FN38 38
#define FN39 39
#define FN40 40
#define FN41 41
#define FN42 42
#define FN43 43
#define FN44 44
#define FN45 45
#define FN46 46
#define FN47 47
#define FN48 48
#define FN49 49
#define FN50 50
#define FN51 51
#define FN52 52
#define FN53 53
#define FN54 54
#define FN55 55

/* Coarse Neighbors for Tables 
 * These are not needed since Fine Neighbors
 * have a direct relationship with Coarse ones
 */

/*
#define CRS0 0
#define CRS1 1
#define CRS2 2
#define CRS3 3
#define CRS4 4
#define CRS5 5
#define CRS6 6
#define CRS7 7
#define CRS8 8
#define CRS9 9
#define CRS10 10
#define CRS11 11
#define CRS12 12
#define CRS13 13
#define CRS14 14
#define CRS15 15
#define CRS16 16
#define CRS17 17
#define CRS18 18
#define CRS19 19
#define CRS20 20
#define CRS21 21
#define CRS22 22
#define CRS23 23
#define CRS24 24
#define CRS25 25
#define CRS26 26
#define CRS27 27
#define CRS28 28
#define CRS29 29
#define CRS30 30
#define CRS31 31
#define CRS32 32
#define CRS33 33
#define CRS34 34
#define CRS35 35
#define CRS36 36
#define CRS37 37
#define CRS38 38
#define CRS39 39
#define CRS40 40
#define CRS41 41
#define CRS42 42
#define CRS43 43
#define CRS44 44
#define CRS45 45
#define CRS46 46
#define CRS47 47
#define CRS48 48
#define CRS49 49
#define CRS50 50
#define CRS51 51
#define CRS52 52
#define CRS53 53
#define CRS54 54
#define CRS55 55
*/

/* Listing of possible W neighbors
 *
 *						  ^ z
 *   					  |
 *						  |
 *          		  |
 *         y <-------
 *
 *
 *          x = x_B - 3*rad_B/2
 *
 *	|-----------------------------------|
 *	|     |     |     |     |     |     |
 *	|  35 |  29 |  23 |  17 |  11 |  5  |
 *	|-----------------------------------|
 *	|     |     |     |     |     |     |
 *	|  34 |  28 |  22 |  16 |  10 |  4  |
 *	|-----------------------------------|
 *	|     |     |     |     |     |     |
 *	|  33 |  27 |  21 |  15 |  9  |  3  |
 *	|-----------------------------------|
 *	|     |     |     |     |     |     |
 *	|  32 |  26 |  20 |  14 |  8  |  2  |
 *	|-----------------------------------|
 *	|     |     |     |     |     |     |
 *	|  31 |  25 |  19 |  13 |  7  |  1  |
 *	|-----------------------------------|
 *	|     |     |     |     |     |     |
 *	|  30 |  24 |  18 |  12 |  6  |  0  |
 *	|-----------------------------------|
 *
 * For symmetries, the above are defined as:
 * (See Syms.hpp for how the rest are mapped)
 *	|-----------------------------------|
 *	|     |     |     |     |     |     |
 *	|     |     |     |     |     |     |
 *	|-----------------------------------|
 *	|     |     |     |     |     |     |
 *	|     |     |     |     |     |     |
 *	|-----------------------------------|
 *	|     |     |     |     |     |     |
 *	|     |     |     |     |     |     |
 *	|-----------------------------------|
 *	|     |     |     |     |     |     |
 *	|     |     |     |WX_F0|WX_E1|WX_E0|
 *	|-----------------------------------|
 *	|     |     |     |     |     |     |
 *	|     |     |     |     |WX_C2|WX_C1|
 *	|-----------------------------------|
 *	|     |     |     |     |     |     |
 *	|     |     |     |     |     |WX_C0|
 *	|-----------------------------------|
 *
 *
 *           x = x_B - rad_B
 *
 *	|-----------------------------------|
 *	|     |     |     |     |     |     |
 *	|  55 |  49 |  47 |  45 |  43 |  41 |
 *	|-----------------------------------|
 *	|     |     |     |     |     |     |
 *	|  54 | FN15| FN11| FN7 | FN3 |  40 |
 *	|-----------------------------------|
 *	|     |     |     |     |     |     |
 *	|  53 | FN14| FN10| FN6 | FN2 |  39 |
 *	|-----------------------------------|
 *	|     |     |     |     |     |     |
 *	|  52 | FN13| FN9 | FN5 | FN1 |  38 |
 *	|-----------------------------------|
 *	|     |     |     |     |     |     |
 *	|  51 | FN12| FN8 | FN4 | FN0 |  37 |
 *	|-----------------------------------|
 *	|     |     |     |     |     |     |
 *	|  50 |  48 |  46 |  44 |  42 |  36 |
 *	|-----------------------------------|
 *
 *
 * Again for symmetries:
 * |-----------------------------------|
 *	|     |     |     |     |     |     |
 *	|     |     |     |     |     |     |
 *	|-----------------------------------|
 *	|     |     |     |     |     |     |
 *	|     | FN15| FN11| FN7 | FN3 |     |
 *	|-----------------------------------|
 *	|     |     |     |     |     |     |
 *	|     | FN14| FN10| FN6 | FN2 |     |
 *	|-----------------------------------|
 *	|     |     |     |     |     |     |
 *	|     | FN13| FN9 | FN5 | FN1 |     |
 *	|-----------------------------------|
 *	|     |     |     |     |     |     |
 *	|     | FN12| FN8 | FN4 | FN0 |     |
 *	|-----------------------------------|
 *	|     |     |     |     |     |     |
 *	|     |     |     |     |     |	   |
 *	|-----------------------------------|
 *
 *
 *          x = x_B - rad_B/2
 *
 *	|-----------------------------------|
 *	|     |     |     |     |     |     |
 *	|  75 |  69 |  67 |  65 |  63 |  61 |
 *	|-----------------------------------|
 *	|     |     |     |     |     |     |
 *	|  74 | FN27| FN23| FN21| FN19|  60 |
 *	|-----------------------------------|
 *	|     |     |           |     |     |
 *	|  73 | FN26|           | FN18|  59 |
 *	|------------     B     ------------|
 *	|     |     |           |     |     |
 *	|  72 | FN25|           | FN17|  58 |
 *	|-----------------------------------|
 *	|     |     |     |     |     |     |
 *	|  71 | FN24| FN22| FN20| FN16|  57 |
 *	|-----------------------------------|
 *	|     |     |     |     |     |     |
 *	|  70 |  68 |  66 |  64 |  62 |  56 |
 *	|-----------------------------------|
 *
 *          x = x_B + rad_B/2
 *
 *	|-----------------------------------|
 *	|     |     |     |     |     |     |
 *	|  95 |  89 |  87 |  85 |  83 |  81 |
 *	|-----------------------------------|
 *	|     |     |     |     |     |     |
 *	|  94 | FN39| FN35| FN33| FN31|  80 |
 *	|-----------------------------------|
 *	|     |     |           |     |     |
 *	|  93 | FN38|           | FN30|  79 |
 *	|------------     B     ------------|
 *	|     |     |           |     |     |
 *	|  92 | FN37|           | FN29|  78 |
 *	|-----------------------------------|
 *	|     |     |     |     |     |     |
 *	|  91 | FN36| FN34| FN32| FN28|  77 |
 *	|-----------------------------------|
 *	|     |     |     |     |     |     |
 *	|  90 |  88 |  86 |  84 |  82 |  76 |
 *	|-----------------------------------|
 *
 *          x = x_B + rad_B
 *
 *	|-----------------------------------|
 *	|     |     |     |     |     |     |
 *	| 115 | 109 | 107 | 105 | 103 | 101 |
 *	|-----------------------------------|
 *	|     |     |     |     |     |     |
 *	| 114 | FN55| FN51| FN47| FN43| 100 |
 *	|-----------------------------------|
 *	|     |     |     |     |     |     |
 *	| 113 | FN54| FN50| FN46| FN42|  99 |
 *	|-----------------------------------|
 *	|     |     |     |     |     |     |
 *	| 112 | FN53| FN49| FN45| FN41|  98 |
 *	|-----------------------------------|
 *	|     |     |     |     |     |     |
 *	| 111 | FN52| FN48| FN44| FN40|  97 |
 *	|-----------------------------------|
 *	|     |     |     |     |     |     |
 *	| 110 | 108 | 106 | 104 | 102 |  96 |
 *	|-----------------------------------|
 *
 *          x = x_B + 3*rad_B/2
 *
 *	|-----------------------------------|
 *	|     |     |     |     |     |     |
 *	| 151 | 145 | 139 | 133 | 127 | 121 |
 *	|-----------------------------------|
 *	|     |     |     |     |     |     |
 *	| 150 | 144 | 138 | 132 | 126 | 120 |
 *	|-----------------------------------|
 *	|     |     |     |     |     |     |
 *	| 149 | 143 | 137 | 131 | 125 | 119 |
 *	|-----------------------------------|
 *	|     |     |     |     |     |     |
 *	| 148 | 142 | 136 | 130 | 124 | 118 |
 *	|-----------------------------------|
 *	|     |     |     |     |     |     |
 *	| 147 | 141 | 135 | 129 | 123 | 117 |
 *	|-----------------------------------|
 *	|     |     |     |     |     |     |
 *	| 146 | 140 | 134 | 128 | 122 | 116 |
 *	|-----------------------------------|
 *
 */


/* W Nodes for Tables */

#define WN0 0
#define WN1 1
#define WN2 2
#define WN3 3
#define WN4 4
#define WN5 5
#define WN6 6
#define WN7 7
#define WN8 8
#define WN9 9
#define WN10 10
#define WN11 11
#define WN12 12
#define WN13 13
#define WN14 14
#define WN15 15
#define WN16 16
#define WN17 17
#define WN18 18
#define WN19 19
#define WN20 20
#define WN21 21
#define WN22 22
#define WN23 23
#define WN24 24
#define WN25 25
#define WN26 26
#define WN27 27
#define WN28 28
#define WN29 29
#define WN30 30
#define WN31 31
#define WN32 32
#define WN33 33
#define WN34 34
#define WN35 35
#define WN36 36
#define WN37 37
#define WN38 38
#define WN39 39
#define WN40 40
#define WN41 41
#define WN42 42
#define WN43 43
#define WN44 44
#define WN45 45
#define WN46 46
#define WN47 47
#define WN48 48
#define WN49 49
#define WN50 50
#define WN51 51
#define WN52 52
#define WN53 53
#define WN54 54
#define WN55 55
#define WN56 56
#define WN57 57
#define WN58 58
#define WN59 59
#define WN60 60
#define WN61 61
#define WN62 62
#define WN63 63
#define WN64 64
#define WN65 65
#define WN66 66
#define WN67 67
#define WN68 68
#define WN69 69
#define WN70 70
#define WN71 71
#define WN72 72
#define WN73 73
#define WN74 74
#define WN75 75
#define WN76 76
#define WN77 77
#define WN78 78
#define WN79 79
#define WN80 80
#define WN81 81
#define WN82 82
#define WN83 83
#define WN84 84
#define WN85 85
#define WN86 86
#define WN87 87
#define WN88 88
#define WN89 89
#define WN90 90
#define WN91 91
#define WN92 92
#define WN93 93
#define WN94 94
#define WN95 95
#define WN96 96
#define WN97 97
#define WN98 98
#define WN99 99
#define WN100 100
#define WN101 101
#define WN102 102
#define WN103 103
#define WN104 104
#define WN105 105
#define WN106 106
#define WN107 107
#define WN108 108
#define WN109 109
#define WN110 110
#define WN111 111
#define WN112 112
#define WN113 113
#define WN114 114
#define WN115 115
#define WN116 116
#define WN117 117
#define WN118 118
#define WN119 119
#define WN120 120
#define WN121 121
#define WN122 122
#define WN123 123
#define WN124 124
#define WN125 125
#define WN126 126
#define WN127 127
#define WN128 128
#define WN129 129
#define WN130 130
#define WN131 131
#define WN132 132
#define WN133 133
#define WN134 134
#define WN135 135
#define WN136 136
#define WN137 137
#define WN138 138
#define WN139 139
#define WN140 140
#define WN141 141
#define WN142 142
#define WN143 143
#define WN144 144
#define WN145 145
#define WN146 146
#define WN147 147
#define WN148 148
#define WN149 149
#define WN150 150
#define WN151 151

/* X Nodes for Tables */
/* X nodes have an inverse relationship with Wnodes.  So,
 * for example, XN0 == WN151
 * So, they need to be defined and are left commented out
 */

/*
#define XN0 0
#define XN1 1
#define XN2 2
#define XN3 3
#define XN4 4
#define XN5 5
#define XN6 6
#define XN7 7
#define XN8 8
#define XN9 9
#define XN10 10
#define XN11 11
#define XN12 12
#define XN13 13
#define XN14 14
#define XN15 15
#define XN16 16
#define XN17 17
#define XN18 18
#define XN19 19
#define XN20 20
#define XN21 21
#define XN22 22
#define XN23 23
#define XN24 24
#define XN25 25
#define XN26 26
#define XN27 27
#define XN28 28
#define XN29 29
#define XN30 30
#define XN31 31
#define XN32 32
#define XN33 33
#define XN34 34
#define XN35 35
#define XN36 36
#define XN37 37
#define XN38 38
#define XN39 39
#define XN40 40
#define XN41 41
#define XN42 42
#define XN43 43
#define XN44 44
#define XN45 45
#define XN46 46
#define XN47 47
#define XN48 48
#define XN49 49
#define XN50 50
#define XN51 51
#define XN52 52
#define XN53 53
#define XN54 54
#define XN55 55
#define XN56 56
#define XN57 57
#define XN58 58
#define XN59 59
#define XN60 60
#define XN61 61
#define XN62 62
#define XN63 63
#define XN64 64
#define XN65 65
#define XN66 66
#define XN67 67
#define XN68 68
#define XN69 69
#define XN70 70
#define XN71 71
#define XN72 72
#define XN73 73
#define XN74 74
#define XN75 75
#define XN76 76
#define XN77 77
#define XN78 78
#define XN79 79
#define XN80 80
#define XN81 81
#define XN82 82
#define XN83 83
#define XN84 84
#define XN85 85
#define XN86 86
#define XN87 87
#define XN88 88
#define XN89 89
#define XN90 90
#define XN91 91
#define XN92 92
#define XN93 93
#define XN94 94
#define XN95 95
#define XN96 96
#define XN97 97
#define XN98 98
#define XN99 99
#define XN100 100
#define XN101 101
#define XN102 102
#define XN103 103
#define XN104 104
#define XN105 105
#define XN106 106
#define XN107 107
#define XN108 108
#define XN109 109
#define XN110 110
#define XN111 111
#define XN112 112
#define XN113 113
#define XN114 114
#define XN115 115
#define XN116 116
#define XN117 117
#define XN118 118
#define XN119 119
#define XN120 120
#define XN121 121
#define XN122 122
#define XN123 123
#define XN124 124
#define XN125 125
#define XN126 126
#define XN127 127
#define XN128 128
#define XN129 129
#define XN130 130
#define XN131 131
#define XN132 132
#define XN133 133
#define XN134 134
#define XN135 135
#define XN136 136
#define XN137 137
#define XN138 138
#define XN139 139
#define XN140 140
#define XN141 141
#define XN142 142
#define XN143 143
#define XN144 144
#define XN145 145
#define XN146 146
#define XN147 147
#define XN148 148
#define XN149 149
#define XN150 150
#define XN151 151
*/

END_EBI_NAMESPACE

#endif
