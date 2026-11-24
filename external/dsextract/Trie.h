#ifndef TRIE_TREE_H
#define TRIE_TREE_H

#include <iostream>
#include <vector>
#include <algorithm>
#include <string>
#include <stack>


using namespace std;

class Trie;
class TrieNode;
        struct sortString{
                bool operator()(const unsigned int & rpStart, const unsigned int & rpEnd ){
			return(rpStart < rpEnd);
                }
        };

class TrieNode {
	unsigned int label;
	vector<unsigned int> vertexSet;
	int nChildren;
	int depth;
	TrieNode *firstChild;
	TrieNode *nextSibling; 
	TrieNode *parent;

	TrieNode(unsigned int elem, unsigned int posLabel, TrieNode *fChild, TrieNode *nChild){ vertexSet.push_back(elem); label = posLabel; firstChild = fChild; nextSibling = nChild; nChildren=0; parent = NULL; depth = 0; }
	friend class Trie;
};

class Trie {
       public: 
        Trie();
        ~Trie();
	unsigned int findKey(unsigned int key);
	void printTrie();
	void insert(unsigned int v, vector<unsigned int> &outlink);
	bool isEmpty();
	void makeEmpty();
	void generateVirtualNodeList(vector<pair<vector<unsigned int>, vector<unsigned int> > > & potentialVN);
	void generateListRec(TrieNode* ptr, vector<pair<vector<unsigned int>, vector<unsigned int> > > & potentialVN, vector<unsigned int> &vnOut);
	void generateListIt(TrieNode* ptr, vector<pair<vector<unsigned int>, vector<unsigned int> > > & potentialVN, vector<unsigned int> &vnOut);
	int ncount;

       private:
        TrieNode *root;
	void insert(unsigned int v, vector<unsigned int> &outlink, TrieNode * trieNode);
	void insert(unsigned int v, unsigned int posLabel, TrieNode * & t);
	TrieNode* find (unsigned int posLabel, TrieNode* t);
	void printTrieRec (TrieNode* t);
	void printNodeContent(TrieNode* t);
	void makeEmpty (TrieNode*& t);
	void addToSet(unsigned int v, unsigned int posLabel, TrieNode * & t);
	void printVector(vector<pair<vector<unsigned int>, vector<unsigned int> > > & v);
 	void printV( vector<unsigned int> v);

	void getChildren(TrieNode *ptr, vector<TrieNode *> &children);

        struct myComp{
                bool operator()(const unsigned int & rpStart, const unsigned int & rpEnd ){
                        if(rpStart ==  rpEnd)
				return true;
			else 
				return false;
                }
        };

        struct sortByLabel{
                bool operator()(const TrieNode * rpStart, const TrieNode * rpEnd ){
                        return (rpStart->label < rpEnd->label);
                }
        };


	struct CompareList{
    		bool operator()(const pair<vector<unsigned int>, vector<unsigned int> > & rpStart, const pair<vector<unsigned int>, vector<unsigned int> > & rpEnd ){
			vector<unsigned int>::iterator it; 
			vector<unsigned int> v1 = rpStart.first;
			vector<unsigned int> v2 = rpEnd.first;
			//cout<<" v1 size "<<v1.size()<<" v2 size "<<v2.size()<<"\n";
			if(v1.size() != v2.size()) return false;
			bool flag = equal(v1.begin(), v1.end(), v2.begin(), myComp());
			return flag;
    		}
	};

};

#endif





