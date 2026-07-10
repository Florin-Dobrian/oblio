#include "oblio/SparseMatrix.h"
#include "oblio/Permutation.h"
#include "test_util.h"
#include <cstdint>
#include <iostream>
using Oblio::SparseMatrix; using Oblio::Permutation;
int main(){ int pass=0,fail=0;
  auto ck=[&](bool ok,const char* n){ std::cout<<"  "<<n<<(ok?"  PASS":"  FAIL")<<"\n"; ok?++pass:++fail; };
  // full-symmetric tridiag n=4: col0{0,1} col1{0,1,2} col2{1,2,3} col3{2,3}
  SparseMatrix<double> A(4,{0,2,5,8,10},{0,1, 0,1,2, 1,2,3, 2,3},{2,-1, -1,2,-1, -1,2,-1, -1,2});
  ck(OblioTest::isStructurallySymmetric(A), "A structurally symmetric");
  ck(A.numCols()==4,"numCols == 4        ");
  ck(A.nnz()==10,   "nnz == 10 (full)    ");
  Permutation p(A.numCols());
  ck(p.size()==4 && p.validate(),"identity perm valid ");
  bool rt=true; for(std::size_t i=0;i<p.size();++i) rt=rt&&static_cast<std::size_t>(p.newToOld()[static_cast<std::size_t>(p.oldToNew()[i])])==i;
  ck(rt,"perm round-trip     ");
  std::cout<<"port smoke: "<<pass<<"/"<<(pass+fail)<<" passed\n"; return fail==0?0:1; }
