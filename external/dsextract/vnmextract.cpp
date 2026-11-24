#include <iostream>

//#include <boost/program_options.hpp>

#include "AdjacencyMatrix.h"
#include "Shingles.h"
#include "Utils.h"

#include <ctime>                           // for timing

using namespace std;

//static unsigned currentPass;

int main(int argc, char *argv[]){

  unsigned int shingle_size;
  unsigned int num_hashes = 8; 
  int format = 0;
  string outfile;

  if(argc < 8){
	cout<<"usage: ./vnmextract graph format[bin format] shingle_size iters bcsize(separated by ,) outputDir_and_file num_hashes\n";
	exit(1);
  }
  format = atoi(argv[2]);
  // shingle_size should be always 1
  shingle_size = atoi(argv[3]);
  unsigned iters = atoi(argv[4]);
  string bcsizeStr = argv[5];
  outfile = argv[6]; 
  num_hashes = atoi(argv[7]);
  
  vector<int> bcsizeVector = AdjacencyMatrix::splitToInt(bcsizeStr, ",");
  cout<<" filename input "<<argv[1]<<" format "<<format<<" shingle_size "<<shingle_size<<" outfile "<<outfile<<endl;

  AdjacencyMatrix ds;
  if(format == 0){
	cout<<"txt format not supported\n";
	exit(1);
  } else {
	// format is -1 2 3 4 -2 4 5 6 -5 1 2 3 4
	// format does not include nodes without outlinks 
	ds.loadBinFileFranFormat(argv[1]);
  }

  clock_t start, startini, finish;
  double time;
  start = clock();
  startini = clock();

  string filename;
  string itstr;
  int bcsize;

  unsigned i=0, k =0;
  while (true) {
	cout<<" num_hashes "<<num_hashes<<"\n";
	bcsize = bcsizeVector[k];
  	int numberRows = ds.getNumberROWS();
  	start = clock();
  	Shingles *sh = new Shingles(num_hashes, numberRows, bcsize);
  	finish = clock();
  	time = double(finish - start) / CLOCKS_PER_SEC;
  	cout << "\t elapsed time Shingles constructor " << time << endl;
  	start = clock();
  	sh->computeShingles(ds, shingle_size);
  	finish = clock();
  	time = double(finish - start) / CLOCKS_PER_SEC;
  	cout << "\t elapsed time Shingles computeShingles " << time << endl;
/*
  	start = clock();
	sh->streamClustersAndMine(numberRows, ds);
  	finish = clock();
  	time = double(finish - start) / CLOCKS_PER_SEC;
  	cout << "\t elapsed time Shingles streamClustersAndMine " << time << endl;
*/

  	PotentialClusters pClusters;
  	start = clock();
  	sh->computePotentialClusters(pClusters, numberRows);
  	finish = clock();
  	time = double(finish - start) / CLOCKS_PER_SEC;
  	cout << "\t elapsed time Shingles computePotentialClusters " << time << endl;
  	start = clock();
  	sh->genClusters(2,pClusters,ds);
  	finish = clock();

  	time = double(finish - start) / CLOCKS_PER_SEC;
  	cout << "\t elapsed time Shingles genClusters " << time << endl;

  	//cout<<"Adjacency Matrix i "<<i<<"\n";
	itstr = AdjacencyMatrix::int2String((int)i);
	filename = outfile + "-it-" + itstr;
	
  	start = clock();
  	sh->cleanAll();
  	finish = clock();
  	time = double(finish - start) / CLOCKS_PER_SEC;
  	cout << "\t elapsed time Shingles cleanAll " << time << endl;
  	numberRows = ds.getNumberROWS();
  	start = clock();
	sh->writeBicliquesDisk(outfile,i);	
  	finish = clock();
  	time = double(finish - start) / CLOCKS_PER_SEC;
  	cout << "\t elapsed time Shingles writeBicliquesDisk " << time << endl;
	unsigned threshold = sh->numberBicliques();
	cout<<" threshold "<<threshold<<endl;
  	start = clock();
  	ds.dumpFile(filename, threshold, iters, i);
  	finish = clock();
  	time = double(finish - start) / CLOCKS_PER_SEC;
  	cout << "\t elapsed time AdjacencyList dumpFile " << time << endl;
	sh->clearBicliques();
	delete sh;
	i++;
	if(threshold < iters){
		if(k < bcsizeVector.size() - 1){
			k++;
		}
		else {
			//cout<<" going out "<<endl;
			break;
		}
	} 


 }

  finish = clock();
  time = double(finish - startini) / CLOCKS_PER_SEC;
  cout << "\t elapsed time " << time << endl;


}
