#ifndef TREEWRAPPER_H
#define TREEWRAPPER_H

#include "TTree.h"
#include "TLeaf.h"
#include "TBranch.h"

#include <map>
#include <string>
#include <vector>

using namespace std;

struct LeafAndBranch
{
  TLeaf* leaf;
  TBranch* branch;
};

/**
   Typical use:

   TTree *t = (TTree*)somefile.Get("sometree");
   TreeWrapper w(t);
   double energy=w.GetValue("evt.energy", 1);

*/
class TreeWrapper : public TObject {
public:
  /**
     Create a TreeWrapper for TTree t. TChain derives from TTree so
     can pass a TChain here too. On second thoughts that might not
     work. Who knows?
  */
  TreeWrapper(TTree* t);

  // Stupid default constructor forced by ROOT. Don't use
  TreeWrapper() {}

  virtual ~TreeWrapper() {}
  /**
     Register a branch that you want to use, a la
     TTree::SetBranchAddress. All you need here is the branch name
     though - no pointers to pointers etc

     It is now unneccessary to call this function. It'll happen automagically
     on first call of GetValue.
  */
  virtual bool AddBranch(const string& branchName);
  /**
     Get the value of branch \a branchName in tree entry \a ientry.
     Use non-zero \a leafVal for leaves that are arrays
  */
  virtual double GetValue(const string& branchName, int ientry, int leafVal=0);

  /**
     Number of entries in the tree
  */
  virtual int GetEntries() const { return tree->GetEntries(); }

  virtual TTree* GetTree() const { return tree; }

protected:
  TTree* tree;

  map<string, LeafAndBranch> leavesAndBranches;
  typedef map<string, LeafAndBranch>::iterator itLaB;

  // Are we wrapping a TChain?
  bool wrappingChain;
  // If we're on a TChain, the tree number of the last entry we read
  int lastReadTree;
  // The entry offset for the current tree
  int currentOffset;
  /**
     Reset the branch addresses because we're on a new TTree (if wrapping a chain)
  */
  virtual Bool_t Notify();

  // Find the tree number containing entry
  int GetTree(int entry, int guess = -1);
  // Set the branch and leaf addresses to the ones appropriate for tree number treeNum
  bool SetBranchAddresses();
  // Get the entry offset for tree treeNum
  int GetOffset(int treeNum);
};

#endif
