/* Kernel Independent Fast Multipole Method
   Copyright (C) 2004 Lexing Ying, New York University

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2, or (at your option)
any later version.

This program is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
for more details.

You should have received a copy of the GNU General Public License
along with this program; see the file COPYING.  If not, write to the Free
Software Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
02111-1307, USA.  */

#ifndef _EXSOL3D_HPP_
#define _EXSOL3D_HPP_

#include "common/kernel3d.hpp"

BEGIN_EBI_NAMESPACE

enum {
  QNT_U = 0,
  QNT_P = 1,
  QNT_RHS = 2,
  QNT_MAX_U = 3,
  QNT_MAX_RHS = 4,
  QNT_ERR = -1
};

enum {
  CHS_EMPTY = 0,
  CHS_LAP_POLY_TEST = 1,
  CHS_LAP_ZRO = 10,
  CHS_LAP_CST = 11,
  CHS_LAP_XYZ = 12,
  CHS_LAP_X2  = 13,
  CHS_LAP_ESQ = 14,	
  CHS_LAP_EXTHAR = 15,
  CHS_LAP_X2Y2Z2 = 16,
  CHS_LAP_ANA = 17,
  CHS_LAP_HAR = 18,
  CHS_HARPER =  20,
  CHS_MODHEL_CST = 21,
  CHS_MODHEL_XYZ = 22,
  CHS_MODHEL_ESQ = 24,
  CHS_MODHEL_EXTHAR = 25,
  CHS_MODHEL_HAR = 28,
  CHS_STK_CST = 31,
  CHS_STK_STK = 32,
  CHS_STK_ROT = 33,
  CHS_STK_PB2 = 34,
  CHS_STK_RBR = 35, //rigid body rotation
  CHS_STK_HAR = 36,
  CHS_STK_ESQ = 37,
  CHS_STK_SIN = 38,
  CHS_STK_SPH = 39,
  CHS_NAV_CST = 51,
  CHS_NAV_NAV = 52,
  CHS_NAV_ROT = 53,
  CHS_NAV_PB3 = 54,
  CHS_LAP_ATAN1 = 81,
  CHS_LAP_ATAN2 = 82,
  CHS_LAP_ATAN5 = 85,
  CHS_LAP_ATAN10 = 810,
  CHS_LAP_TOR_ATAN10 = 8110,
  CHS_LAP_ATAN20 = 820,
  CHS_LAP_ATAN30 = 830,
  CHS_LAP_ATAN40 = 840,
  CHS_LAP_ATAN80 = 880,
  CHS_STK_ATAN1 = 91,
  CHS_STK_ATAN2 = 92,
  CHS_STK_ATAN5 = 95,
  CHS_STK_ATAN10 = 910,
  CHS_STK_ATAN20 = 920,
  CHS_STK_ATAN30 = 930,
  CHS_STK_ATAN40 = 940,
  CHS_STK_ATAN80 = 980,

  CHS_LAP_FREE_SPACE = 100,
  CHS_LAP_SPH = 101,
  CHS_LAP_COL = 102,
  CHS_LAP_PER = 103,
  CHS_LAP_DIR = 104,
  CHS_MODHEL_FREE_SPACE = 200,
  CHS_STK_FREE_SPACE = 300,
  
  //error
  CHS_ERR = -1
};

class Exsol3d
{
protected:
  int _et; //equation type
  vector<double> _coefs; //coefs
  int _ct; //choice type  //int _qt; //quantity type
  
public:
  Exsol3d() {;}
  Exsol3d(int et, const vector<double>& coefs, int ct): _et(et), _coefs(coefs), _ct(ct) {;}
  int& et() { return _et; }
  vector<double>& coefs() { return _coefs; }
  int& ct() { return _ct; }
  int tdof(int qt);
  int quantity(int qt, const DblNumMat& trgpos, DblNumVec& trgqnt);
};

END_EBI_NAMESPACE

#endif
