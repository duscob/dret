#ifndef ADJACENCY_MATRIX
#define ADJACENCY_MATRIX 
#include <vector>
#include <map>
#include <string>
#include <algorithm>

using namespace std;

class AdjacencyMatrix{
  
 public:
  
  static string int2String(int i);
  static vector<int> splitToInt(string str, string delim);

  unsigned getNumberROWS(){
	return matrix.size();
  }
  

// first item in pair is the outlink list in string format and the next element in pair is 
// the string that represent the vertex 
//  typedef vector<pair<string,string> > MatrixMap;  
  typedef pair<vector<unsigned int> , unsigned int>  MatrixNode;  
  typedef vector<MatrixNode *> MatrixMap;  

  /* 
   * Read the adjacency list file
   *
   * \param meta index for filename in francisco claude txt format
   *
   */
  void loadBinFileFranFormat(const char * filename);

  // stdout the matrix content by rows
  void dumpFile();

  // print matrix to a file content by rows
  void dumpFile(string filename, int numberBC, int minTH, int iter);
  void dumpRow(int row);

  /*
   * iterators
   */
  typedef MatrixMap::iterator matrix_iterator;
  typedef MatrixMap::const_iterator const_iterator;
  inline matrix_iterator begin() { return matrix.begin();};
  inline matrix_iterator end() { return matrix.end();};
  inline const_iterator cbegin() { return matrix.begin();};
  inline const_iterator cend() { return matrix.end();};
  inline long size() {return matrix.size();}


  AdjacencyMatrix(){totNodes = 0;};
  ~AdjacencyMatrix();

  int getNode(){return totNodes;}
  void nextNode(){totNodes++;}
  int totNodes;


// private:
  MatrixMap matrix;
  bool valid;
};


#endif
