#ifndef __TREE_HPP__
#define __TREE_HPP__

#include "common/nummat.hpp"
#include "common/ebiobject.hpp"
#include "common/ebi.hpp"

BEGIN_EBI_NAMESPACE
/*
 * Required interface for any tree code to work properly with legacy markgrid()
 * TODO kill this entirely...
 */
class TreeNode {
   
    protected:
        int _color;
        int _numIntPts;
        double _dist;
        vector<char> _isIntDbl;
        vector<char> _isInt;
        int _numIntDblPts;

    public:
        vector<int> targets_in_box;

        //virtual ~TreeNode() {;}
        int& color(){ return _color; } 
        int& numIntPts(){ return _numIntPts; }
        double& dist() { return _dist; } 

        vector<char>& isInt() { return _isInt; } 
        vector<char>& isIntDbl() { return _isIntDbl; } 
        char& isIntDbl(int i) { return _isIntDbl[i]; } 
        char& isInt(int i) { return _isInt[i]; } 
        int& numIntDblPts(){ return _numIntDblPts; } 

        //MJM BUG use size_t instead of ints. Left as ints to enforce the 
        //template but could introduce bugs as number of points grow
        // TODO make virtual on template class
        vector<int>& trgOwnVecIdxs(){ return targets_in_box; }

};
class Tree {
  protected:
      uint _kSrcVal;
      DblNumMat _grdSrcSamPos;
      DblNumMat _grdDblSrcSamPos;

  public:
    // TODO once KIFMM is removed, RENAME EVERYTHING HERE
    //virtual Tree() = 0;
    virtual double radius() = 0;
    virtual double radius(int node_id) = 0;
    virtual Point3 center(int node_id) = 0;
    virtual bool terminal(int node_id) = 0;
    virtual TreeNode& node(int node_id) = 0;
    virtual int child(int node_id, Index3& index) = 0;
    
    virtual uint& kSrcVal() = 0;
    virtual uint srcGrdSze() = 0;
    virtual uint dblSrcGrdSze() = 0;

    virtual const DblNumMat& grdDblSrcSamPos() = 0;
    //virtual const DblNumMat& grdSrcSamPos() = 0;
    virtual DblNumMat grdDblSrcExaPos(int node_id, bool depth=false) = 0;
    virtual DblNumMat grdSrcExaPos(int node_id, bool depth=false) = 0;
    
    virtual int dwnOrderCollect(vector<int>& ret) = 0;
    virtual int upwOrderCollect(vector<int>& ret) = 0;
    virtual int revLvlOrderCollect(map<int, vector<int> >& ret) = 0;

};
END_EBI_NAMESPACE
#endif

