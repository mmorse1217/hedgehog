#ifndef __PERDIRMAPS_HPP__
#define __PERDIRMAPS_HPP__

#include "common/ebi.hpp"
#include <iostream>

BEGIN_EBI_NAMESPACE

using std::map;
using std::pair;

class PerDirMaps{
protected:
  map<int, DblNumMat> _PerUpwChk2UpwEqu;
  map<int, NumTns<DblNumMat> > _PerUpwEqu2UpwChk;
  map<int, DblNumMat> _PerDwnChk2DwnEqu;
  map<int, NumTns<DblNumMat> > _PerDwnEqu2DwnChk;
  map<int, OffTns<DblNumMat> > _PerUpwEqu2DwnChk;
  
  /* For Dirichlet Stuff */
  map<int, DblNumMat> _dwnEquGrdMap;

  DblNumVec _rtTrgDwnEquDen;
  //DblNumVec _rtTrgDwnChkVal;
  
public:
  PerDirMaps() {;}
  ~PerDirMaps(){;}
  map<int, DblNumMat>& dwnEquGrdMap(){ return _dwnEquGrdMap; }

  map<int, DblNumMat>& PerUpwChk2UpwEqu() { return _PerUpwChk2UpwEqu; }
  map<int, NumTns<DblNumMat> >& PerUpwEqu2UpwChk() { return _PerUpwEqu2UpwChk; }
  map<int, DblNumMat>& PerDwnChk2DwnEqu() { return _PerDwnChk2DwnEqu; }
  map<int, NumTns<DblNumMat> >& PerDwnEqu2DwnChk() { return _PerDwnEqu2DwnChk; }
  map<int, OffTns<DblNumMat> >& PerUpwEqu2DwnChk() { return _PerUpwEqu2DwnChk; }

  DblNumVec& rtTrgDwnEquDen() { return _rtTrgDwnEquDen; }
	 
};

END_EBI_NAMESPACE

#endif
