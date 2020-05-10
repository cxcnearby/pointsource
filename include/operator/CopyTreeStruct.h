#ifndef POINTSOURCE_INCLUDE_OPERATOR_COPYTREESTRUCT_H_
#define POINTSOURCE_INCLUDE_OPERATOR_COPYTREESTRUCT_H_

#include "TLeaf.h"
#include <TChain.h>
#include <TFile.h>
#include <TTree.h>
#include <string>

class CopyTreeStruct {
public:
  CopyTreeStruct();
  CopyTreeStruct(const char *root_file_name);
  CopyTreeStruct(const char *root_file_name, const TString tree_name);
  ~CopyTreeStruct();
  int Init();
  int SetReadTree(TChain *&chain);
  int SetWriteTree(TTree *&tree);
  TString GetTreeName();
  bool IsExists(const TString leaf_name);
  double operator[](const TString leaf_name) const;

private:
  struct LeafInfo {
    TString leafname;
    TString datatype;
    int datanumber;
    TString leafword;
    Int_t data_int[1000];
    Long64_t data_long64[1000];
    Short_t data_short[1000];
    Float_t data_float[1000];
    Double_t data_double[1000];
    Bool_t data_bool[1000];
    // could add more!
  };
  std::map<TString, LeafInfo> mLeafInfo_;
  const char *template_name_;
  TString treename_;
  int init_status_ = -1;
  int setreadtree_status_ = -1;
  int setwritetree_status_ = -1;
};

#endif // COPYTREESTRUCT_H_
