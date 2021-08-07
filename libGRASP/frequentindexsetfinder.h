#include <vector>
#include <set>
class FrequentIndexSetFinder {
private:
	double minSupport;
	std::vector<std::vector<double> > data;
public:
	
	FrequentIndexSetFinder(const std::vector<std::vector<double> >& _data, int _minSupport);
	std::vector< std::vector<int> > getFrequentIndexSets();
	std::vector< std::set<int> > getCombinations(int n);
}
