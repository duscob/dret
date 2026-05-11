#include "AdjacencyMatrix.h"


#include <iostream>
#include <fstream>
#include <sstream>
#include <stdlib.h>
#include <cassert>

vector<string> split(string line, string delims){
   string::size_type bi, ei;
   vector<string> words;
   bi = line.find_first_not_of(delims);     
   while(bi != string::npos) {
  	ei = line.find_first_of(delims, bi);
    	if(ei == string::npos)
      		ei = line.length();
    		words.push_back(line.substr(bi, ei-bi));
    		bi = line.find_first_not_of(delims, ei);
  	}                 
  return words;
}

vector<int> AdjacencyMatrix::splitToInt(string line, string delims){
   string::size_type bi, ei;
   vector<int> words;
   bi = line.find_first_not_of(delims);     
   while(bi != string::npos) {
  	ei = line.find_first_of(delims, bi);
    	if(ei == string::npos)
      		ei = line.length();
		string aux = line.substr(bi,ei-bi);
    		words.push_back(atoi(aux.c_str()));
    		bi = line.find_first_not_of(delims, ei);
  	}                 
  return words;
}

string AdjacencyMatrix::int2String(int i){
   stringstream ss;
   string str;
   ss<<i;
   ss>>str;
   return str;
}

AdjacencyMatrix::~AdjacencyMatrix(){
	MatrixNode *mnode;	
	for(unsigned int i=0; i<matrix.size(); i++){
		mnode = matrix[i];
		delete mnode;
	}	
	matrix.clear();
}

//load bin file in Francisco format

void AdjacencyMatrix::loadBinFileFranFormat(const char * filename){
  FILE *fp = fopen64(filename,"r");
  if(fp==NULL)
        cout<<" error opening file "<<filename<<"\n";
  assert(fp!=NULL);

  uint rd = 0;

  uint nodes; uint edges;
  rd += fread(&nodes,sizeof(uint),1,fp);
  rd += fread(&edges,sizeof(uint),1,fp);

  unsigned int currentNode = 0;
  vector<unsigned int> aux;
  vector<unsigned int>::iterator vit;
  MatrixNode *mnode;

  cout<<" nodes "<<nodes<<" edges "<<edges<<"\n";
  for(uint i=0;i<nodes+edges;i++) {
    int k;
    rd += fread(&k,sizeof(unsigned int),1,fp);
    if(k<0) {
      if(aux.size() != 0){
	    // add to the outlink the same node to catch cliques
            /* Disabled, as we are only interested in blicliques.
               Jouni Sirén, 2012-09-27
	    vit = find(aux.begin(), aux.end(), currentNode);
	    if(vit == aux.end()){
		aux.push_back(currentNode);
	    	sort(aux.begin(), aux.end());
	    }
	    */
	    mnode = new pair< vector <unsigned int>, unsigned int>(aux,currentNode);
            matrix.push_back(mnode);
      }
      if(totNodes < -k)
	totNodes = -k;
      aux.clear();
      currentNode = -k;

    }
    else {
      aux.push_back(k);
      if(totNodes < k)
	totNodes = k;
    }
  }
  if(aux.size() != 0){
     /*
       Disabled, see above.
     vit = find(aux.begin(), aux.end(), currentNode);
     if(vit == aux.end()){
          aux.push_back(currentNode);
          sort(aux.begin(), aux.end());
     }
     */

     mnode = new pair< vector <unsigned int>, unsigned int>(aux,currentNode);
     matrix.push_back(mnode);
  }
  fclose(fp);
  totNodes++;
  cout<<"vns comienzan en "<<totNodes<<"\n";
  cout.flush();

}


// dump a single row 


void AdjacencyMatrix::dumpRow(int row){

  cout << "`" << matrix[row]->second << "'\t";
  vector<unsigned int> outlink = matrix[row]->first;
  vector<unsigned int>::iterator it;
  for(it = outlink.begin(); it != outlink.end(); it++) 
	cout<<*it<<" ";
  cout<<"\n";
}


// dump all the rows 

void AdjacencyMatrix::dumpFile(){

  MatrixMap::iterator end = matrix.end();

  for (MatrixMap::iterator it = matrix.begin();
       it != end;
       it ++){
	MatrixNode *p = *it;
	cout<<p->second<<":";
  	vector<unsigned int> outlink = p->first;
  	vector<unsigned int>::iterator vit;
  	for(vit = outlink.begin(); vit != outlink.end(); vit++) 
		cout<<*vit<<" ";
  	cout<<"\n";
  
  }
}

// dump all the rows to a file

void AdjacencyMatrix::dumpFile(string filename, int numberBC, int minTH, int iter) {
  
  vector<unsigned int>::iterator vit;
  ofstream *out = NULL; 
  //cout<<"iter "<<iter<<" numberBC "<<numberBC<<" minTH "<<minTH<<endl;
  // mod 100 to avoid to dump many files
  if(iter != 0 && (numberBC < minTH || (iter%100) == 0)){
  //if(numberBC < minTH ){
  	//ofstream out(filename.c_str());
	out = new ofstream(filename.c_str());
  }
  MatrixMap::iterator it=matrix.begin();
  while(it != matrix.end()){
           MatrixNode *p = *it;
	   //string str = int2String(p->second);
  	   vector<unsigned int> outlink = p->first;
	   
	   if(outlink.size() != 0 ){
  	     if(iter != 0 && (numberBC < minTH || (iter%100) == 0)){
  	     //if(numberBC < minTH ){
		*out<<p->second;
		*out<<": ";
  		vector<unsigned int>::iterator vit;
		//string s;
  		for(vit = outlink.begin(); vit != outlink.end(); vit++) {
			//s = s + " " + int2String(*vit);
			*out<<*vit;
			*out<<" ";
		}
		*out<<"\n";
	     }
	     it++;
 	   } else {
		//it=matrix.erase(it);
		it++;
		//cout<<" deleting p  second "<<p->second<<endl;
		//delete p;
  		//cout<<"erased row "<<str<<endl;
	   }
  }
  if(iter != 0 && (numberBC < minTH || (iter%100) == 0)){
  //if(numberBC < minTH ){
  	out->close();
	delete out;
  }
}


