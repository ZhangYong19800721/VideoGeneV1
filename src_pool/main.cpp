#include <stdio.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <iomanip>
#include "GenePool.h"
using namespace std;

int main() {
	GenePool pool = GenePool();
	pool.load(0,"Titanic.gene");
	pool.load(1,"xwlb.gene");

	GenePool::VDID id;
	GenePool::INDXPAIRSEQ s;
	long score = pool.find("xwlb-part001-part004.gene",id,s);

	cout<<"final match score = "<<score<<endl;
	for (int i = 0; i<s.size();++i){
		cout<<"x="<<s[i].first<<",y="<<s[i].second<<endl;
	}

	return 0;
}
