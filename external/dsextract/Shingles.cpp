#include "Shingles.h"
#include "conf.h"

#include <iostream>
#include <algorithm>

#include <ctime>                           // for timing
#include <cmath>                           // for pow

#include <functional>                      // for hashing a word

#include <stdlib.h>

VNodeList Shingles::globalVNList;

void Shingles::cleanAll(){
        signNode.clear();
        prime_ =  static_cast<unsigned int>(pow(2.0,31.0)-1);
        srand(time(NULL));

        for (unsigned int i = 0; i < numberSignatures_; i++){

                A_.push_back( rand() % prime_ + 1);
                B_.push_back( rand() % prime_ + 1);
                MIN_.push_back(prime_ + 1);
        }


}


//
// compute shingles
//

void Shingles::computeShingles(AdjacencyMatrix & ds,
			       const unsigned short shingleSize){

  clock_t start, finish;
  double time;

  start = clock();

  // create the hasher
  std::hash<string> string_hash;

  // vector containing #shingleSize tokenized words 
  typedef vector<size_t> vShingles;
  
  // a single shingleID
  size_t shingleID;

  // a single shingleHASH
  unsigned int shingleHash;

  
  //  
  // for each row in the adjacency matrix 
  //
  AdjacencyMatrix::MatrixMap::const_iterator it_end=(ds.matrix).end();
  
  int di = 0;
  for(AdjacencyMatrix::MatrixMap::const_iterator it=(ds.matrix).begin(); it != it_end; 
      ++it, ++di){
    // REINIZIALIZE MIN
    for (unsigned int i = 0; i < numberSignatures_; i++)
      MIN_[i] = prime_;

    vShingles vHashedWords(shingleSize);

    // for each token
    unsigned int i=0;
      AdjacencyMatrix::MatrixNode *p = *it;
      vector<unsigned int> outlink = p->first;

	if (outlink.size() == 0)continue;

      vector<unsigned int>::iterator tok_iter;
     for(tok_iter = outlink.begin(); tok_iter != outlink.end(); tok_iter++){ 
      vHashedWords[ i % shingleSize ] = string_hash(AdjacencyMatrix::int2String(*tok_iter));

      if (i+1 >= shingleSize){  // skip to next if not yet a shingle-word
       
	// compute the shingle ID as sum of shingleSize words
	shingleID = 0;
	vShingles::const_iterator vs_it_end = vHashedWords.end();
	for(vShingles::const_iterator vs_it = vHashedWords.begin();
	    vs_it != vs_it_end;
	    vs_it++){
	  
	  shingleID = (shingleID + *vs_it);
	}   // each shingle
	
	// given a shingleID
	//   compute each hash function
	for (size_signatures j = 0; j < numberSignatures_; j++){
	  
	  // 32 + 32 bit = 64 mod prime_ -> 32 bit 
	  shingleHash = (((unsigned long long) A_[j] * 
			  (unsigned long long) shingleID)  + B_[j]) % prime_;
	
	  // take min_hash
	  if (shingleHash < MIN_[j]){
	    MIN_[j] = shingleHash;
	  }
	} // each signature
      } // if is a valid shingle-word

    } 
    

    for (size_signatures j = 0; j < numberSignatures_; j++){
      signNode.push_back(make_pair(MIN_[j], di));      
    }
    
  } // each vertex 

  finish = clock();
  time = double(finish - start) / CLOCKS_PER_SEC;

  if(debug){
   	cout<<"---- printSignNode\n";
  	printSignNode();
  }
}

void Shingles::writeBicliquesDisk(string outfile, unsigned currentPass){
        ofstream myfile;
        string filename = outfile + "-biclique-it-";
        filename = filename + AdjacencyMatrix::int2String(currentPass) + ".txt";
        myfile.open(filename.c_str());

        if (myfile.is_open()){
           BicliqueVector::iterator iter;
           for(iter = bicliqueOutVector.begin();  iter != bicliqueOutVector.end(); iter++){
                vector<unsigned int> vertexList = iter->first;
                vector<unsigned int> outlinks = iter->second;
                vector<unsigned int>::iterator viter;
		sort(vertexList.begin(), vertexList.end());
                for(viter = vertexList.begin(); viter != vertexList.end(); viter++){
                        myfile<<*viter;
                        myfile<<" ";
                }
                vector<unsigned int>::iterator oiter;
                myfile<<"- ";
                for(oiter = outlinks.begin(); oiter != outlinks.end(); oiter++){
                        //fprintf(fileout,"%s ", *oiter);
                        myfile<<*oiter;
                        myfile<<" ";
                }
                //fprintf(fileout,"\n", *oiter);
                myfile<<"\n";
           }
        }
	myfile.close();
}

void Shingles::clearBicliques(){
        bicliqueOutVector.clear();
}

int Shingles::numberBicliques(){
        return bicliqueOutVector.size();
}


// build m_rows which is a matrix of Shingles where 
// the first column is the adjacency matrix index row 

void Shingles::computePotentialClusters(PotentialClusters &potClusters, int numberROWS){

  unsigned int i=0;
  int m=0;
  vector<unsigned int> m_rows;
  vector<unsigned int> s_rows;
  //m_rows.reserve((numberSignatures_+1)*numberROWS);
  vector<unsigned int>::iterator iter;
  for(SignNodeVect::const_iterator it = signNode.begin();
      it != signNode.end();
     ++it, ++i){
  
    if((i % numberSignatures_) == 0){
	//cout<<it->second <<" "<<it->first;
	m_rows.push_back((unsigned int)it->second);
	m_rows.push_back(it->first);
	m = 0;
    } else {
	m++;
	//cout<<" "<<it->first;
	m_rows.push_back(it->first);
    }
    if(m == (numberSignatures_ - 1)){
	//cout<<"\n";
	;
    }
  }
  
  // it gets the hash column specified by second parameter sorted
  SignNodeVect col; 
  //col.reserve(numberROWS);
  col.clear();
  getSortColHash(m_rows, 1, col, 0, m_rows.size());

  // makes groups based on the hash value (key) associated to a vector of indexes asociated to rows of m_rows
  map<unsigned int, vector<int> > groups;
  //cout<<" sorted first column \n";
  makeGroupIes(groups, col);

  //make a new matrix now sorted by the hash column 1 in this case , after this m_rows can be eliminated
  // posibleClusters contains the row positions in s_rows that contains the same first hash columns
  vector<int> possibleClusters; 
  makeSortMatrix(s_rows, m_rows, groups, possibleClusters);

  possibleClusters.push_back(s_rows.size()/(numberSignatures_+1)*(numberSignatures_+1));
  
  m_rows.clear();
  potClusters.first = s_rows; 
  potClusters.second = possibleClusters;
}

void Shingles::genClusters(int column, PotentialClusters &pot, AdjacencyMatrix &ds){
  vector<unsigned int> cr = pot.first;
  vector<int> pos = pot.second;
  SignNodeVect col; 
  map<unsigned int, vector<int> > groups;

  for(unsigned k=0; k<pos.size()-1; k++){
     int num = pos[k+1] - pos[k];
     if(numberSignatures_ + 1 < num ){
  	vector<unsigned int> fc;
	//cout<<" fc begin \n";
	for(int i=pos[k]; i<pos[k+1]-1; i=i+numberSignatures_+1){
		pair<vector<unsigned int>, unsigned int> *adjit = (ds.matrix)[cr[i]];
		vector<unsigned int> outl = adjit->first;
		if(outl.size() == 0){
			//cout<<"fc outlink size 0 \n";
			continue;
		}
     		for(int j=0; j < (numberSignatures_ + 1); j++){
			fc.push_back(cr[i + j]);
     		}
  	}
	

	col.clear();
  	getSortColHash(fc, column, col, 0, fc.size());
	//cout<<" GROUPS \n";
	groups.clear();
  	makeGroupIes(groups, col);
  	vector<unsigned int> c_rows;
  	vector<int> posClusters; 
	posClusters.clear();
	
  	makeSortMatrix(c_rows, fc, groups, posClusters);
  	posClusters.push_back(c_rows.size());
	
	// stores rows that match previous column, so they are mine together
  	vector<unsigned int> fc2;

  	for(unsigned kk=0; kk<posClusters.size()-1; kk++){
     		int num = posClusters[kk+1] - posClusters[kk];
		vector<int> posc;
		posc.push_back(posClusters[kk]);
		posc.push_back(posClusters[kk + 1]);

		int numberEntries = num/(numberSignatures_ + 1);
		PotentialClusters pC;
		pC.first = c_rows;
		pC.second = posc;

  		//if(numberEntries > 10 && column < numberSignatures_ + 1){
  		if(numberEntries > 10 && column < numberSignatures_){
			genClusters(column + 1, pC, ds);
 		} else if(numberEntries > 1 ){
			//cout<<" can call miningPhase\n";
  			vector<unsigned int> fc1;
			//cout<<" fc1 begin \n";
			for(int i=posClusters[kk]; i<posClusters[kk+1]-1; i=i+numberSignatures_+1){
				pair<vector<unsigned int>, unsigned int> *adjit = (ds.matrix)[c_rows[i]];
				vector<unsigned int> outl = adjit->first;
				if(outl.size() == 0){
					cout<<"fc 1 outlink size 0 \n";
					continue;
				}
     				for(int j=0; j < (numberSignatures_ + 1); j++){
        				//cout<<cr[i + j]<<" ";
					fc1.push_back(c_rows[i + j]);
     				}
     				//cout<<"\n";
  			}
 			//cout<<" fc1 size "<<fc1.size()<<"\n";
	
		   	if(debug){			
				cout<<" For MiningPhase fc1 begin \n";
				for(unsigned i=0; i<fc1.size(); i=i+numberSignatures_+1){
     					for(int j=0; j < (numberSignatures_ + 1); j++){
						cout<<fc1[i + j]<<" ";
					} 
     					cout<<"\n";
				}
 				cout<<" For MiningPHase end fc1\n";
			}
			
			miningPhase(fc1, ds);
		} else if(numberEntries == 1){
			//cout<<" fc2 begin \n";
			for(int i=posClusters[kk]; i<posClusters[kk+1]-1; i=i+numberSignatures_+1){
				pair<vector<unsigned int>, unsigned int> *adjit = (ds.matrix)[c_rows[i]];
				vector<unsigned int> outl = adjit->first;
				if(outl.size() == 0){
					cout<<"fc 2 outlink size 0 \n";
					continue;
				}
     				for(int j=0; j < (numberSignatures_ + 1); j++){
        				//cout<<cr[i + j]<<" ";
					fc2.push_back(c_rows[i + j]);
     				}
     				//cout<<"\n";
  			}
		}
	}
	//cout<<" fc2 size "<<fc2.size()<<"\n";
	if(fc2.size() > 0){
	    if(debug){	
		cout<<" For MiningPhase fc2 begin column "<<column<<"\n";
		for(unsigned i=0; i<fc2.size(); i=i+numberSignatures_+1){
     			for(int j=0; j < (numberSignatures_ + 1); j++){
				cout<<fc2[i + j]<<" ";
			} 
     			cout<<"\n";
		}
 		cout<<" For MiningPHase end fc2 column "<<column<<"\n";
	    }
	    miningPhase(fc2,ds);
	}

	
     } // end if
   } // end for

}


// contructs the histogram for each node in the outlink lists
// counts is a vector that holds a pair with the edge number with its respective counter
// counts is sorted based on the counter of each edge using myfunction() compare method
// input: c_rows is the Hash matrix, and AdjacencyMatrix is the original matrix that
//        contains a vector of pairs in matrix field (where first element contains outlink in a string
//        and the second element contains the vertex) 
// output: mapVertexOut is a map with key as vertex and value as vector<int> as its outlink
//         according to counts (where edges with counter < 2 are discarded)

void Shingles::computeHistogram(vector<unsigned int> &c_rows, const AdjacencyMatrix &ds, 
				map<unsigned int, vector<unsigned int> > & mapVertexOut){
   map<unsigned int, vector<unsigned int> > mapVertexOutlink;
   map<unsigned int, vector<unsigned int> >::iterator mit;
   SignNodeVectString counts;
   SignNodeVectString::iterator it;
   for(unsigned k=0; k<c_rows.size(); k=k+numberSignatures_+1){
	AdjacencyMatrix::MatrixNode *adjit = (ds.matrix)[c_rows[k]];
	vector<unsigned int> tokens = adjit->first;
	vector<unsigned int>::iterator tok_iter;
        for (tok_iter = tokens.begin();
                        tok_iter != tokens.end();
                        tok_iter++){
                vector<unsigned int> nodes;
                unsigned int val = *tok_iter;
                unsigned int vert = adjit->second;
                mit = mapVertexOutlink.find(val);
                if(mit != mapVertexOutlink.end()){
                        nodes = mit->second;
                }
                nodes.push_back(vert);
                mapVertexOutlink[val] = nodes;
                SignNodeVectString::iterator pit;
                SignNodeVectString myvect;
                myvect.push_back(make_pair(val,1));
                pit = search(counts.begin(), counts.end(), myvect.begin(), myvect.end(), myfunction());

                if(pit != counts.end()){
                        //cout<<"val "<<val<<"  pit first "<<pit->first<<endl;
                        //cout<<" found edge "<<pit->first<<" with count "<<pit->second<<endl;
                        (pit->second)++;
                } else {
                        //cout<<" adding edge "<<val<<endl;
                        counts.push_back(make_pair(val, 1));
                }

        }

   }

   sort(counts.begin(), counts.end(), myCompare());
   //cout<<"count ends\n";

   SelectedEdges selectedEdges;
   for(mit=mapVertexOutlink.begin(); mit != mapVertexOutlink.end(); mit++){
        vector<unsigned int> nodes = mit->second;
	if(nodes.size() > 1){
		sort(nodes.begin(),nodes.end());
		selectedEdges.push_back(make_pair(mit->first, nodes));
	}
   }


   sort(selectedEdges.begin(), selectedEdges.end(), edgeCmp());
   SelectedEdges::iterator sit;
   for(sit = selectedEdges.begin(); sit != selectedEdges.end(); sit++){
        vector<unsigned int> nodes = sit->second;
	//cout<<"edge "<<sit->first<<": ";
        vector<unsigned int>::iterator it;
	for(it = nodes.begin(); it != nodes.end(); it++){
		vector<unsigned int> edges;
		//cout<<*it<<" ";
		mit = mapVertexOut.find(*it);
		if(mit != mapVertexOut.end()){
			edges = mit->second;
		}
		edges.push_back(sit->first);
		mapVertexOut[*it] = edges;
	}
	//cout<<"\n";
   }

   if(debug){
   	cout<<"counts \n";
   	for(it=counts.begin(); it != counts.end(); it++){
        	cout<<" edge "<<it->first<<" count "<<it->second<<"\n";
   	}
   }

}


// mapVertexOutlink contains the vertexes with their respective outlinks without edges with counter < 2

void Shingles::miningPhase(vector<unsigned int> &c_rows, AdjacencyMatrix &ds){
  map<unsigned int, vector<unsigned int> > mapVertexOutlink;
  map<unsigned int, vector<unsigned int> >::iterator mit;

  computeHistogram(c_rows,ds,mapVertexOutlink);
  
  Trie t;
  if(debug)
  	cout<<" map vertex outlink\n";
  for(mit= mapVertexOutlink.begin(); mit != mapVertexOutlink.end(); mit++){
	if(debug){
        	cout<<" vertex "<<mit->first<<" ";
        	for(vector<unsigned int>::iterator vit = (mit->second).begin(); vit != (mit->second).end(); vit++){
              		cout<<*vit<<" ";
        	}
        	cout<<"\n";
	}
	t.insert(mit->first,mit->second);
  }
  if(debug)
  	t.printTrie();
  if(t.isEmpty())return;

  VNodeList mapVNList;
  VNodeList::iterator vnit;

  vector<pair< vector<unsigned int> , vector<unsigned int> > > potentialVN;
  vector<pair< vector<unsigned int> , vector<unsigned int> > >::iterator it;
  t.generateVirtualNodeList(potentialVN);
  t.makeEmpty();
  sort(potentialVN.begin(),potentialVN.end(),CompareSaving());
  //cout<<"POTENTIAL VN\n";
  string vn = "v";
  for(it = potentialVN.begin(); it != potentialVN.end(); it++){
        pair< vector<unsigned int>, vector<unsigned int> > p = *it;
        vector<unsigned int> vnOutlinks = p.second;
        vector<unsigned int> vertexList = p.first;
	if((vnOutlinks.size() - 1) * (vertexList.size() - 1) == 1){
		if(debug){
			cout<<" saving = 1 vnOutlinks\n";
			printVector(vnOutlinks);
			cout<<"\n vertexs \n";
			printVector(vertexList);
		}
		continue;
	}
	//cout<<" vnid "<<uniqueid<<"\n";
	//printVector(vnOutlinks);
	
        unsigned int vnid = ds.getNode();

	pair<unsigned int,bool> myid;
	myid.first = vnid;
	myid.second = false;
	
	VNodeList::iterator itvnmap = globalVNList.find(vnOutlinks);
	if(itvnmap != globalVNList.end()){
		//cout<<" vnid found in global vnlist id \n";
		//cout<<(itvnmap->second).first<<"\n";
		mapVNList[vnOutlinks] = itvnmap->second;
	} else {
        	mapVNList[vnOutlinks] = myid;
		globalVNList[vnOutlinks] = myid;
	}
	ds.nextNode();
  }

  if(debug){
  	cout<<"vn dict \n";

  	for(vnit = mapVNList.begin(); vnit != mapVNList.end(); vnit++){
        	cout<<" vnid "<<(vnit->second).first<<" ";
        	vector<unsigned int> vnout = vnit->first;
        	for(vector<unsigned int>::iterator vit = vnout.begin(); vit != vnout.end(); vit++){
                	cout<<" "<<*vit;
        	}
        	cout<<"\n";
  	}
  }

  //cout<<"reemplazing vn dict in AdjacencyMatrix\n";
  reemplazing(c_rows,mapVNList, ds, potentialVN);
}

bool Shingles::transVector(AdjacencyMatrix::MatrixNode *adj, vector<unsigned int> outlink){
        vector<unsigned int> nodes;
        //cout<<"adj.first "<<adj.first<<"\n";
	vector<unsigned int> tokens = adj->first;
	vector<unsigned int>::iterator tok_iter;
        for (tok_iter = tokens.begin(); tok_iter != tokens.end(); tok_iter++){
             nodes.push_back(*tok_iter);
        }
        return includes(nodes.begin(), nodes.end(), outlink.begin(), outlink.end(), sortString());
}

// Reeplaces virtual nodes in AdjacencyMatrix, which contains the matrix in adjacency format

void Shingles::reemplazing(vector<unsigned int> &c_rows, 
  	VNodeList & mapVNList, AdjacencyMatrix &ds,
  	vector<pair< vector<unsigned int> , vector<unsigned int> > > &potentialVN){


  vector<pair< vector<unsigned int> , vector<unsigned int> > >::iterator it;
  VNodeList::iterator vnit;
  map<unsigned int, vector<unsigned int> >::iterator mit;

  map<unsigned int,int> mapVertexPosition;
  map<unsigned int,int>::iterator mapVPos;


  for(it = potentialVN.begin(); it != potentialVN.end(); it++){
        pair< vector<unsigned int> , vector<unsigned int> > p = *it;
        vector<unsigned int> vertexList = p.first;
        vector<unsigned int> vnOutlink = p.second;

        vector<unsigned int>::iterator vsit, vsit1;
        vector<int>::iterator indel;

	sort(vnOutlink.begin(), vnOutlink.end(), sortString());
	if(debug){
        	cout<<"\n vnOutlinks\n";
		printVector(vnOutlink);
        	cout<<"\n vertexList\n";
		printVector(vertexList);
	}

//eliminates vertex from vertexList if vertex already is used for another vn and
// its outlink has already changed with that previous vn

        vector<unsigned int> paborrar;
        // verifies whether a vertex in a vertexList was already modified by another virtual node
        if(mapVertexPosition.size() > 0){
                for(vsit = vertexList.begin(); vsit != vertexList.end(); vsit++){
                        mapVPos = mapVertexPosition.find(*vsit);
                        if(mapVPos != mapVertexPosition.end()){
                                //pair<vector<unsigned int>,unsigned int> nadj = (ds.matrix)[mapVPos->second];
                                AdjacencyMatrix::MatrixNode *nadj = (ds.matrix)[mapVPos->second];
			
                                if(!transVector(nadj, vnOutlink)){
                                        paborrar.push_back(*vsit);
                                }

                        }
                }

        }      
        for(vsit=paborrar.begin(); vsit != paborrar.end(); vsit++){
                vsit1 = find(vertexList.begin(), vertexList.end(), *vsit);
                if(vsit1 != vertexList.end()){
                        vertexList.erase(vsit1);
                }
        }
	paborrar.clear();

        if(vnOutlink.size() * vertexList.size() <= Shingles::bcsize || vertexList.size() < 2)
               continue;
	//cout<<" dumping ds "<<endl;
	//ds.dumpFile();

	vector<unsigned int> vertexBiclique;

	for(unsigned k=0; k<c_rows.size(); k=k+numberSignatures_+1){
        	AdjacencyMatrix::MatrixNode * adjit = (ds.matrix)[c_rows[k]];
                vector<unsigned int>::iterator vsit = find(vertexList.begin(), vertexList.end(), adjit->second);
                if(vsit == vertexList.end())continue;

	
		vector<unsigned int> tokens = adjit->first;
		vector<unsigned int>::iterator tok_iter;
                vector<unsigned int> nodes;
        	for (tok_iter = tokens.begin(); tok_iter != tokens.end(); tok_iter++){
                	unsigned int val = *tok_iter;
			if(val != 0)
				nodes.push_back(val);
		}

                vector<unsigned int>::iterator itdiff;
                vector<unsigned int>::iterator mitit;
		sort(nodes.begin(), nodes.end(), sortString());
		if(debug){
			cout<<" Nodes  for vertex "<<adjit->second<<"\n";
			printVector(nodes);
		}
		pair<unsigned int, bool> myid = mapVNList[vnOutlink];
		if(myid.first == adjit->second){
			//cout<<" found vnode igual a vertexList must not replace vnode outlink"<<adjit->second<<"\n";
			continue;
		}
		bool flag = false;
                if(includes(nodes.begin(), nodes.end(), vnOutlink.begin(), vnOutlink.end(), sortString())){
                	for(itdiff=vnOutlink.begin(); itdiff != vnOutlink.end(); itdiff++){
                		mitit = find(nodes.begin(), nodes.end(), *itdiff);
                		if(mitit != nodes.end()){
                			nodes.erase(mitit);
                        	}
                	}
                	// now must add vn id ... needs to change mapVNList. vector vertexes maps to vnid
			myid.second = true;
			mapVNList[vnOutlink] = myid;
			flag = true;
		}
		if(debug){
			cout<<" Nodes after reempl \n";
			printVector(nodes);
		}
		if(flag){
			
			vector<unsigned int> nnodes;
			for(unsigned int r=0; r<nodes.size(); r++){
				if(nodes[r] != 0)
					nnodes.push_back(nodes[r]);
			}
			
			(adjit->first).clear();
			adjit->first = nnodes;
			(ds.matrix)[c_rows[k]] = adjit;
			vertexBiclique.push_back(adjit->second);
                        mapVertexPosition[adjit->second] = c_rows[k];

		}
        }
       if(vertexBiclique.size() > 0){
                bicliqueOutVector[vertexBiclique] = vnOutlink;
        }

	//cout<<" after reempl in for dumping ds "<<endl;
	//ds.dumpFile();

  }

}

// print a vector
void Shingles::printVector(vector<unsigned int> &v){
	for(vector<unsigned int>::iterator it = v.begin(); it != v.end(); it++){
		cout<<*it<<" ";
	}
	cout<<"\n";
}

// Builds column sorted by hash value (shingle) (pair, first hash value, second row position)
// input: matrix matRows, number of column to sort, indexes ini and end for the matrix range to sort
// output: SignNodeVect col (sorted by shingle column)
// 

void Shingles::getSortColHash(vector<unsigned int> &matRows, int colNumber, SignNodeVect &col, int ini, int end){
  for(int i=ini; i < end; i=i+numberSignatures_+1){
	//cout<<matRows[i + 1]<<" ";
      	col.push_back(make_pair(matRows[i + colNumber], i));      
     	//cout<<"\n";
  }
  sort(col.begin(), col.end(), compareSignNode());
}

// Builds a map that associates a hash value from a column with the row positions in matrix,
// the idea is to group by column hash value with the set of row positions that contains the
// same col hash value
// input: vector col of pairs with first Hash value and second row position
// output: group (key : hash value) (value : vector of row positions)

void Shingles::makeGroupIes(map<unsigned int, vector<int> > &groups, SignNodeVect &col){
  int i=0;
  map<unsigned int, vector<int> >::iterator mit;
  for(SignNodeVect::const_iterator it = col.begin();
      it != col.end();
      ++it,++i){
	//cout<<"col "<<it->first<<" "<<it->second<<"\n";
	vector<int> ies;
	mit = groups.find(it->first);
	if(mit != groups.end()){
		ies = mit->second;			
	}
	ies.push_back(it->second);
	groups[it->first] = ies;
  }

}

// From a matrix build a new matrix sorted by groups defined in groups
// input: matRows (matrix to sort), groups (key: hash value, value: row positions)
// output: s_rows (sorted matRows based on groups), possibleClusters (vector of real positions in matrix 
//         of shingles implemented as a one dimension vector) , possibleClusters variable is used for next 
//	   clustering using next column 

void Shingles::makeSortMatrix(vector<unsigned int> &s_rows, vector<unsigned int> &matRows, map<unsigned int, vector<int> > &groups, vector<int> &possibleClusters){
  map<unsigned int, vector<int> >::iterator mit;
  int k = 0;
  for(mit=groups.begin(); mit!=groups.end(); mit++){
     	possibleClusters.push_back(k);
     	for(unsigned i=0; i<(mit->second).size(); i++){
     		for(int j = 0; j<(numberSignatures_ + 1); j++){
	 		s_rows.push_back(matRows[(mit->second)[i] + j]);
			k++;
		}
     	}
  }


}

void Shingles::printSignNodeMIN(){
  for(SignNodeVect::const_iterator it = signNode.begin(); it != signNode.end(); ++it){
	cout<<it->second <<" "<<it->first<<"\n";
  }
}

void Shingles::printSignNode(){
  unsigned int i=0;
  int m=0;
  for(SignNodeVect::const_iterator it = signNode.begin(); it != signNode.end(); ++it, ++i){
  
    if((i % numberSignatures_) == 0){
	cout<<it->second <<" "<<it->first;
	m = 0;
    } else {
	m++;
	cout<<" "<<it->first;
    }
    if(m == (numberSignatures_ - 1))
	cout<<"\n";
  }

}


