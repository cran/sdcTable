#include <Rcpp.h>
#include <string>
using namespace std;

RcppExport SEXP myPaste(SEXP stringvec, SEXP outvec, SEXP param) {
	SEXP rl = R_NilValue; // return list

	RcppStringVector stringVec(stringvec);
	RcppStringVector outVec(outvec);
	
	RcppParams rp(param);
	int keyVars = rp.getIntValue("nrKeyVars");
	int nrRows = stringVec.size();

	int by = nrRows / keyVars;
	for (int i=0; i < by; i++) {
		std::string str;
		for(int j=0; j < keyVars; j++) {
			str.append(stringVec(i+by*j));			
		}
		outVec(i) = str;
	}	
	RcppResultSet rs;
   	rs.add("str", outVec);
	rl = rs.getReturnList();
	return rl;
}
