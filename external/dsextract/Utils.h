#ifndef UTILS 
#define UTILS 


#include <vector>
#include <map>

using namespace std;

// first contains the shingle matrix sorted by the first column Shingles 
// second contains a vector with the possitions of the first index in the
// matrix (implemented as a vector) 
typedef pair<vector<unsigned int>, vector<int> > PotentialClusters; //  potential clusters matrix

// pair that contains VNodeList with a flag = 0 indicating VN not used or 1 if it is used
typedef map<vector<unsigned int>, pair<unsigned int, bool > > VNodeList; //  potential virtual Nodes 

// vector of clusters of the hash matrix, this already have H1 and H2 column in common

typedef vector< PotentialClusters > HashClusters; //   Hash clusters with H1 and H2 ordered
typedef map<vector<unsigned int>, vector<unsigned int> > BicliqueVector; //  stores bicliques to write to a file later

#endif
