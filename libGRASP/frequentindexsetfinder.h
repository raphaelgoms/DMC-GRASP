#include <vector>
#include <set>
using namespace std;
class FrequentIndexSetFinder {
private:
	double minSupport;
	vector<double> m_min;
	vector<double> m_max;
	vector<vector<double> > data;
	vector<vector<double>> getSubFeatData(vector<int> indexSet);
	vector< vector<int> > getFrequentIndexSets();
	vector< vector<int> > getCombinations(size_t n);
	double getSupport(const vector<int>& indexSet);
	double distance(const vector<double> & v1, const vector<double> & v2);

public:
	FrequentIndexSetFinder(const vector<vector<double> >& _data, 
		const vector<double>& min,
		const vector<double>& max,
		double _minSupport);
};
