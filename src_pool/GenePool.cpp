/*
 * GenePool.cpp
 *
 *  Created on: Oct 23, 2017
 *      Author: root
 */

#include "GenePool.h"

vector<unsigned long> ReadGeneFile(const std::string filename);
unsigned int hamming(const unsigned long A);

GenePool::GenePool() {
	// TODO Auto-generated constructor stub

}

GenePool::~GenePool() {
	// TODO Auto-generated destructor stub
}

int GenePool::load(VDID id, const std::string filename) {
	std::vector<GENE> gene_sequence= ReadGeneFile(filename);
	m_VDID2GENE[id] = gene_sequence;

	hash_map<GENE,std::vector<INDX> > gene2indx;
	for (INDX i = 0; i < gene_sequence.size(); ++i) {
		GENE gene = gene_sequence[i];
		gene2indx[gene].push_back(i);
	}

	typedef hash_map<GENE,std::vector<INDX> > GENEMAP;
	for (GENEMAP::iterator it = gene2indx.begin(); it!=gene2indx.end(); ++it) {
		GENE gene = it->first;
		m_GENE2INDX[gene][id] = it->second;
	}

	return 0;
}

long GenePool::find(const std::string filename, VDID& id_, INDXPAIRSEQ &s_) {
	std::vector<GENE> gene1 = ReadGeneFile(filename);
	long score = 0;
	long highest_score = 0;

	typedef hash_map<VDID,std::vector<std::pair<INDX,INDX> > > VDID2INDXPAIR_MAP;
	VDID2INDXPAIR_MAP vdid_idxp_map;

	for (INDX x1 = 0; x1 < gene1.size(); ++x1) {
		GENE gene = gene1[x1];
		VDID2INDX_MAP genemap = m_GENE2INDX[gene];
		for (VDID2INDX_MAP::iterator it = genemap.begin(); it!=genemap.end(); ++it) {
			VDID id = it->first;

			std::vector<GENE> gene2 = m_VDID2GENE[id];
			for (unsigned long n = 0; n < it->second.size(); ++n) {
				INDX x2 = it->second[n];

				if (vdid_idxp_map[id].end() != std::find(vdid_idxp_map[id].begin(),vdid_idxp_map[id].end(),std::pair<INDX,INDX>(x1,x2))) {
					break;
				}

				INDXPAIRSEQ s; vector<int> dis_list;
				score = GreedyMatch(gene1,gene2,x1,x2,s,dis_list);
				vdid_idxp_map[id].insert(vdid_idxp_map[id].end(),s.begin(),s.end());

				if (score > highest_score) {
					highest_score = score;
					id_ = id;
					s_ = s;
				}
			}
		}
	}


	return highest_score;
}

long GenePool::GreedyMatch(const GENESEQ& gene1, const GENESEQ& gene2, INDX x1, INDX x2, INDXPAIRSEQ &s_, std::vector<int>& dis_vec_) {
	long score = 0;
	long highest_score = 0;

	// search forward
	for (int step=0; (x1+step)<gene1.size() && (x2+step)<gene2.size(); ++step) {
		unsigned int dis = hamming(gene1[x1+step]^gene2[x2+step]);
		s_.push_back(std::pair<INDX,INDX>(x1+step,x2+step));
		dis_vec_.push_back(dis);
		score += dis<=2?8:-5;
		if (score > highest_score) {
			highest_score = score;
		} else if (score <= highest_score - 15) {
			break;
		}
	}

	highest_score = score;
	// search backward
	for (int step=-1; (x1+step)>=0 && (x2+step)>=0; --step) {
		unsigned long dis = hamming(gene1[x1+step]^gene2[x2+step]);
		//s_.insert(s_.begin(),std::pair<INDX,INDX>(x1+step,x2+step)); //TODO to improve efficiency
		s_.push_back(std::pair<INDX,INDX>(x1+step,x2+step));
		dis_vec_.insert(dis_vec_.begin(),dis);
		score += dis<=2?8:-5;
		if (score > highest_score) {
			highest_score = score;
		} else if (score <= highest_score - 15) {
			break;
		}
	}

	return highest_score;
}
