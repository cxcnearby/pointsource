#ifndef POINTSOURCE_INCLUDE_REBASE_COPYTREESTRUCT_H_
#define POINTSOURCE_INCLUDE_REBASE_COPYTREESTRUCT_H_

#include "TLeaf.h"
#include <TChain.h>
#include <TFile.h>
#include <TTree.h>
#include <string>

class CopyTreeStruct {
public:
  struct LeafInfo_ {
    TString leafname_;
    TString datatype_;
    int datanumber_;
    TString leafword_;
    Int_t data_int_[1000];
    Long64_t data_long64_[1000];
    Short_t data_short_[1000];
    Float_t data_float_[1000];
    Double_t data_double_[1000];
    Bool_t data_bool_[1000];
    // could add more!
  };
  std::vector<LeafInfo_> vLeafInfo_;
  CopyTreeStruct();
  CopyTreeStruct(const char *root_file_name);
  CopyTreeStruct(const char *root_file_name, const TString tree_name);
  ~CopyTreeStruct();
  int Init();
  int SetReadTree(TChain *&chain);
  int SetWriteTree(TTree *&tree);
  template <typename T> int SetBond(const char *leaf, T &var);
  // int SetBond(const char *leaf, T *&var);
  TString GetTreeName();
  //  public :
private:
  const char *template_name_;
  TString treename_;
  int init_status_ = -1;
  int setreadtree_status_ = -1;
  int setwritetree_status_ = -1;
};

CopyTreeStruct::CopyTreeStruct() {}
CopyTreeStruct::CopyTreeStruct(const char *root_file_name)
    : template_name_(root_file_name) {}
CopyTreeStruct::CopyTreeStruct(const char *root_file_name,
                               const TString tree_name)
    : template_name_(root_file_name), treename_(tree_name) {}
CopyTreeStruct::~CopyTreeStruct() {}

#endif // COPYTREESTRUCT_H_
