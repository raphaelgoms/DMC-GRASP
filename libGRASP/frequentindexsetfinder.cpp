#include "frequentindexsetfinder.h"
#include "xmeans.h"
#include <cmath>

FrequentIndexSetFinder::FrequentIndexSetFinder(const vector<vector<double> >& _data, 
		const vector<double>& min,
		const vector<double>& max,
		double _minSupport)
{
	m_min = min;
	m_max = max;
	data = _data;
	minSupport = _minSupport;
}

vector< vector<int> > FrequentIndexSetFinder::getFrequentIndexSets() 
{
	vector< vector<int> > frequentIndexSets;

	for (size_t idxSetSize = 2; idxSetSize <= data.size(); idxSetSize++)
	{
		vector<vector<int>> indexSets = getCombinations(idxSetSize);
		for (vector<vector<int>>::iterator it = indexSets.begin(); it != indexSets.end(); ++it) {
			if (getSupport(*it) >= minSupport) {
				frequentIndexSets.push_back(*it);
			}
		}
	} 
	
	return frequentIndexSets;
}

vector< vector<int> > FrequentIndexSetFinder::getCombinations(size_t idxSetSize) {
	// Get all combinations of indexes of size idxSetSize
	vector< vector<int> > indexSets;
	int n = data[0].size();
	int t = idxSetSize; 

	std::vector<int> c(t+2);
	
	for (int j = 0; j < t; j++)
	{
		c[j] = j;
	}
	
	c[t] = n; c[t+1] = 0;

	int j;
	while (true) {
		indexSets.push_back(vector<int>(c.begin(), c.end()-2));
		 
		j = 1; 
		while (c[(j-1)] + 1 == c[(j-1)+1]) {
			c[(j-1)] = j - 1;
			j++;
		}	

		if (j > t) 
			break;

		c[(j-1)] = c[(j-1)] + 1;
	}

	return indexSets;
}



vector<vector<double>> FrequentIndexSetFinder::getSubFeatData(vector<int> indexSet){
	vector<vector<double>> filteredData;
	for (int i : indexSet) {
		filteredData.push_back(data[i]);
	}

	return filteredData;
	
}

double FrequentIndexSetFinder::distance(const vector<double> & v1, const vector<double> & v2) {
	double sum = 0;
	for (int i = 0; i < v1.size(); i++)	{
		sum += (v1[i] - v2[i]) * (v1[i] - v2[i]);
	}
	
	return sqrt(sum);
}

double FrequentIndexSetFinder::getSupport(const vector<int>& indexes) {
	// Get support of indexes
	// vector<vector<double>> filteredData = getSubFeatData(indexes);
	// vector<MCluster> clusters = Xmeans().xMeans(filteredData, data[0].size(), m_min, m_max);

	// int suport = 0;

	// int maxSize = -1;
	// int greatestCluster = -1;
	// for (int i = 0; i < clusters.size(); i++)
	// {
	// 	if (clusters[i].points.size() > maxSize) {
	// 		maxSize = clusters[i].points.size();
	// 		greatestCluster = i;
	// 	}
	// }

	// for (int i = 0; i < data.size(); i++) {
	// 	if (distance(clusters[greatestCluster].centroid, filteredData[i])) {
	// 		suport++;	
	// 	}
	// }
	


	return 0;
}