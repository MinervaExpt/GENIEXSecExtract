#define PR(x) cout << #x " = " << x << "\n";

#include "GENIEXSecExtract/TreeWrapper.h"

#include "TChain.h"
#include "TBranchElement.h"
#include "TStreamerInfo.h"
#include "TVirtualCollectionProxy.h"

#include <iostream>
#include <cstdlib>

TreeWrapper::TreeWrapper(TTree* t)
  : tree(t), wrappingChain(false), lastReadTree(-1), currentOffset(0)
{
  if(t->IsA()->InheritsFrom(TChain::Class())){
    wrappingChain=true;
    t->GetEntries();
    ((TChain*)t)->SetNotify(this);
  }
}

//======================================================================

bool TreeWrapper::AddBranch(const string& branchName)
{
  // Sometimes FindBranch works, sometimes GetBranch works. Don't know
  // why, so just try both
  TBranch* b=tree->FindBranch(branchName.c_str());
  if (!b) b=tree->GetBranch(branchName.c_str());
  TLeaf* l=tree->FindLeaf(branchName.c_str());
  if (!l) l=tree->GetLeaf(branchName.c_str());
  if (!b) {
    cerr << "Can't find branch " << branchName << endl;
    return false;
  }
  if (!l) {
    cerr << "Can't find leaf " << branchName << endl;
    return false;
  }
  tree->SetBranchStatus(branchName.c_str(), 1);

  LeafAndBranch lb;
  lb.leaf = l;
  lb.branch = b;
  leavesAndBranches[branchName] = lb;
  return true;
}

//======================================================================

double TreeWrapper::GetValue(const string& branchName, int ientry, int leafVal)
{
  //cout << "GetValue " << branchName << ", " << ientry << endl;
  itLaB it = leavesAndBranches.find(branchName);

  if (it==leavesAndBranches.end()) {
    if (!AddBranch(branchName)) {
      cerr << "Can't find branch " << branchName << ". Bailing" << endl;
      abort();
    }
    it = leavesAndBranches.find(branchName);
  }

  int localEntry;
  if(wrappingChain){
    localEntry=((TChain*)tree)->LoadTree(ientry);
  }
  else{
    localEntry=ientry;
  }
  //cout << "localEntry=" << localEntry << endl;
  // If we're trying to get a value from a variable sized array, we
  // need to make sure that the variable defining the size of the
  // array has been filled. The leafcount member of TLeaf is the leaf
  // that contains this number, so we fill it here. This is done by
  // TLeafI, but not TLeafF. I don't really understand why.
  //cout << it->second.leaf << endl;
  if(it->second.leaf->GetLeafCount())
    it->second.leaf->GetLeafCount()->GetBranch()->GetEntry(localEntry);
  
  it->second.branch->GetEntry(localEntry);
  return it->second.leaf->GetValue(leafVal);
}

//======================================================================

int TreeWrapper::GetTree(int entry, int guess)
{
  if(!wrappingChain) return 0;

  TChain* chain=(TChain*)tree;
  const int nTrees=chain->GetNtrees();
  Long64_t* offsets=chain->GetTreeOffset();

  if(guess >= 0 && guess < nTrees &&
     offsets[guess] <= entry && offsets[guess+1] > entry) return guess;

  // Dumb linear search. Too stupid to write a binary search
  // For chains of small numbers of files it wouldn't be worth it anyway
  for(int i=0; i<nTrees; ++i){
    if(offsets[i]<=entry && offsets[i+1]>entry) return i;
  }
  // Didn't find it
  return nTrees;
}

//======================================================================

Bool_t TreeWrapper::Notify()
{
  SetBranchAddresses();
  return kTRUE;
}

//======================================================================

bool TreeWrapper::SetBranchAddresses()
{
  // This should only be called when running on a TChain
  if(!wrappingChain) return false;

  const itLaB end = leavesAndBranches.end();
  for (itLaB it=leavesAndBranches.begin(); it!=end; ++it) {
    const string branchName=it->first;
    TBranch* b=tree->FindBranch(branchName.c_str());
    if (!b) b=tree->GetBranch(branchName.c_str());
    TLeaf* l=tree->FindLeaf(branchName.c_str());
    if (!l) l=tree->GetLeaf(branchName.c_str());
    if (!b) {
      cerr << "Can't find branch " << branchName << endl;
      return false;
    }
    if (!l) {
      cerr << "Can't find leaf " << branchName << endl;
      return false;
    }
    tree->SetBranchStatus(branchName.c_str(), 1);

    it->second.leaf = l;
    it->second.branch = b;
  }
  return true;
}

//======================================================================

int TreeWrapper::GetOffset(int treeNum)
{
  return ((TChain*)tree)->GetTreeOffset()[treeNum];
}

