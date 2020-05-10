#include "rebase/CopyTreeStruct.h"
#include <iostream>

CopyTreeStruct::CopyTreeStruct() {}
CopyTreeStruct::CopyTreeStruct(const char *root_file_name)
    : template_name_(root_file_name) {}
CopyTreeStruct::CopyTreeStruct(const char *root_file_name,
                               const TString tree_name)
    : template_name_(root_file_name), treename_(tree_name) {}
CopyTreeStruct::~CopyTreeStruct() {}

int CopyTreeStruct::Init() {
  TFile *file = TFile::Open(template_name_, "READ");
  if (!file) {
    init_status_ = 1;
    return init_status_;
  };
  // TObjArray *tobja_keys = (TObjArray *)file->GetListOfKeys();
  // if (!tobja_keys) {
  //   init_status_ = 2;
  //   return init_status_;
  // }
  // treename_ = tobja_keys->At(0)->GetName(); // An alternative way to get tree
  // name.

  // treename_ = file->GetListOfKeys()->At(0)->GetName(); // if you want to
  // check
  // more tree or hists,
  // use an TIter to get
  // all of them.
  //   treename_ = "trec"; // TODO: Here I use a deliberate name!!!
  if (treename_.IsNull()) {
    treename_ = file->GetListOfKeys()->At(0)->GetName();
    if (treename_.IsNull()) {
      init_status_ = 2;
      return init_status_;
    }
  }
  TTree *tree = (TTree *)file->Get(treename_);
  TIter next(tree->GetListOfLeaves());
  while (TLeaf *leaf = (TLeaf *)next()) {
    LeafInfo_ thisleaf;
    thisleaf.leafname_ = leaf->GetName();
    thisleaf.datatype_ = leaf->GetTypeName();
    thisleaf.datanumber_ = leaf->GetNdata();
    if (thisleaf.datatype_ == "Int_t") {
      if (thisleaf.datanumber_ == 1) {
        thisleaf.leafword_ = thisleaf.leafname_ + "/I";
      } else {
        thisleaf.leafword_ = thisleaf.leafname_ + "[" +
                             std::to_string(thisleaf.datanumber_) + "]/I";
      }
      init_status_ = 0;
    } else if (thisleaf.datatype_ == "Long64_t") {
      if (thisleaf.datanumber_ == 1) {
        thisleaf.leafword_ = thisleaf.leafname_ + "/L";
      } else {
        thisleaf.leafword_ = thisleaf.leafname_ + "[" +
                             std::to_string(thisleaf.datanumber_) + "]/L";
      }
      init_status_ = 0;
    } else if (thisleaf.datatype_ == "Short_t") {
      if (thisleaf.datanumber_ == 1) {
        thisleaf.leafword_ = thisleaf.leafname_ + "/S";
      } else {
        thisleaf.leafword_ = thisleaf.leafname_ + "[" +
                             std::to_string(thisleaf.datanumber_) + "]/S";
      }
      init_status_ = 0;
    } else if (thisleaf.datatype_ == "Float_t") {
      if (thisleaf.datanumber_ == 1) {
        thisleaf.leafword_ = thisleaf.leafname_ + "/F";
      } else {
        thisleaf.leafword_ = thisleaf.leafname_ + "[" +
                             std::to_string(thisleaf.datanumber_) + "]/F";
      }
      init_status_ = 0;
    } else if (thisleaf.datatype_ == "Double_t") {
      if (thisleaf.datanumber_ == 1) {
        thisleaf.leafword_ = thisleaf.leafname_ + "/D";
      } else {
        thisleaf.leafword_ = thisleaf.leafname_ + "[" +
                             std::to_string(thisleaf.datanumber_) + "]/D";
      }
      init_status_ = 0;
    } else if (thisleaf.datatype_ == "Bool_t") {
      if (thisleaf.datanumber_ == 1) {
        thisleaf.leafword_ = thisleaf.leafname_ + "/O";
      } else {
        thisleaf.leafword_ = thisleaf.leafname_ + "[" +
                             std::to_string(thisleaf.datanumber_) + "]/O";
      }
      init_status_ = 0;
    } else {
      init_status_ = 3;
    }
    vLeafInfo_.emplace_back(thisleaf);
    // std::cout << thisleaf.leafname_ << " " << thisleaf.datatype_ << " " <<
    // thisleaf.data_int_.size() << " " << thisleaf.data_float_.size() << " " <<
    // thisleaf.data_double_.size() << " " << thisleaf.data_short_.size() << " "
    // << thisleaf.data_long64_.size() << " " << thisleaf.data_bool_.size() <<
    // std::endl;
  }
  if (file)
    file->Close();
  return init_status_;
} // init_status: 0 == success; 1 == open file failed; 2 == getlistofkeys
  // failed; 3 == datatype not listed.

int CopyTreeStruct::SetReadTree(TChain *&chain) {
  while (init_status_ != 0) {
    if (init_status_ > 0) {
      setreadtree_status_ = 1;
      std::cerr << "CopyTreeStruct::SetReadTree: Init() failed! error code: "
                << init_status_ << std::endl;
      return setreadtree_status_;
    } else if (init_status_ < 0) {
      Init();
    }
  }
  chain = new TChain(treename_);
  for (auto &&v : vLeafInfo_) {
    if (v.datatype_ == "Int_t") {
      chain->SetBranchAddress(v.leafname_, &v.data_int_);
      setreadtree_status_ = 0;
    } else if (v.datatype_ == "Long64_t") {
      chain->SetBranchAddress(v.leafname_, &v.data_long64_);
      setreadtree_status_ = 0;
    } else if (v.datatype_ == "Short_t") {
      chain->SetBranchAddress(v.leafname_, &v.data_short_);
      setreadtree_status_ = 0;
    } else if (v.datatype_ == "Float_t") {
      chain->SetBranchAddress(v.leafname_, &v.data_float_);
      setreadtree_status_ = 0;
    } else if (v.datatype_ == "Double_t") {
      chain->SetBranchAddress(v.leafname_, &v.data_double_);
      setreadtree_status_ = 0;
    } else if (v.datatype_ == "Bool_t") {
      chain->SetBranchAddress(v.leafname_, &v.data_bool_);
      setreadtree_status_ = 0;
    } else {
      setreadtree_status_ = 2;
    }
  }
  return setreadtree_status_;
} // setreadtree_status_: 0 == success; 1 == init error; 2 == data type error.

int CopyTreeStruct::SetWriteTree(TTree *&tree) {
  while (init_status_ != 0) {
    if (init_status_ > 0) {
      setwritetree_status_ = 1;
      std::cerr << "CopyTreeStruct::SetWriteTree: Init() failed! error code: "
                << init_status_ << std::endl;
      return setwritetree_status_;
    } else if (init_status_ < 0) {
      Init();
    }
  }
  tree = new TTree(treename_, "an auto-generated tree based on template");
  for (auto &&v : vLeafInfo_) {
    if (v.datatype_ == "Int_t") {
      tree->Branch(v.leafname_, v.data_int_, v.leafword_);
      setwritetree_status_ = 0;
    } else if (v.datatype_ == "Long64_t") {
      tree->Branch(v.leafname_, v.data_long64_, v.leafword_);
      setwritetree_status_ = 0;
    } else if (v.datatype_ == "Short_t") {
      tree->Branch(v.leafname_, v.data_short_, v.leafword_);
      setwritetree_status_ = 0;
    } else if (v.datatype_ == "Float_t") {
      tree->Branch(v.leafname_, v.data_float_, v.leafword_);
      setwritetree_status_ = 0;
    } else if (v.datatype_ == "Double_t") {
      tree->Branch(v.leafname_, v.data_double_, v.leafword_);
      setwritetree_status_ = 0;
    } else if (v.datatype_ == "Bool_t") {
      tree->Branch(v.leafname_, v.data_bool_, v.leafword_);
      setwritetree_status_ = 0;
    } else {
      setwritetree_status_ = 2;
    }
  }
  return setwritetree_status_;
} // setwritetree_status_: 0 == success; 1 == init error; 2 == data type error.

TString CopyTreeStruct::GetTreeName() {
  while (init_status_ != 0) {
    if (init_status_ > 0) {
      std::cerr << "CopyTreeStruct::GetTreeName: Init() failed! error code: "
                << init_status_ << std::endl;
      exit(1);
    } else if (init_status_ < 0) {
      Init();
    }
  }
  return treename_;
}

template <typename T> int CopyTreeStruct::SetBond(const char *leaf, T &var) {
  // int CopyTreeStruct::SetBond(const char *leaf, T *&var) {
  int setbond_status_ = -1;
  while (init_status_ != 0) {
    if (init_status_ > 0) {
      setbond_status_ = 1;
      std::cerr << "CopyTreeStruct::SetBond: Init() failed! error code: "
                << init_status_ << std::endl;
      return setbond_status_;
    } else if (init_status_ < 0) {
      Init();
    }
  }
  for (auto &&v : vLeafInfo_) {
    if (leaf == v.leafname_) {
      if (v.datatype_ == "Int_t") {
        var.I = v.data_int_;
        setbond_status_ = 0;
      } else if (v.datatype_ == "Long64_t") {
        var.L = v.data_long64_;
        setbond_status_ = 0;
      } else if (v.datatype_ == "Short_t") {
        var.S = v.data_short_;
        setbond_status_ = 0;
      } else if (v.datatype_ == "Float_t") {
        var.F = v.data_float_;
        setbond_status_ = 0;
      } else if (v.datatype_ == "Double_t") {
        var.D = v.data_double_;
        setbond_status_ = 0;
      } else if (v.datatype_ == "Bool_t") {
        var.O = v.data_bool_;
        setbond_status_ = 0;
      } else
        setbond_status_ = 2;
      break;
    }
  }
  if (setbond_status_ == -1) {
    var.I = nullptr;
    var.L = nullptr;
    var.S = nullptr;
    var.F = nullptr;
    var.D = nullptr;
    var.O = nullptr;
  }
  return setbond_status_;
} // setbond_status_: 0 == success; 1 == init error; 2 == data type error; -1 ==
  // not set bond.
