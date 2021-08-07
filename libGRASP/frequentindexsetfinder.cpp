#include <frequentindexsetfinder.h>

FrequentIndexSetFinder::FrequentIndexSetFinder(const std::vector<std::vector<double> >& _data, int _minSupport)
{
	data = _data;
	minSupport = _minSupport;
}

std::vector< std::vector<int> > FrequentIndexSetFinder::getFrequentIndexSets() 
{
	std::vector< std::vector<int> > frequentIndexSets;

	for (int idxSetSize = 2; idxSetSize < data.size(); idxSetSize++)
	{
		
	}
	
	return frequentIndexSets;
}