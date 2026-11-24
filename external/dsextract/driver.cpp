#include <iostream>

#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <algorithm>

#include <sys/timeb.h>

#ifndef CLK_PER_SEC
#ifdef CLOCKS_PER_SEC
#define CLK_PER_SEC CLOCKS_PER_SEC
#endif
#endif

#include "HHH.h"

using namespace std;

double dblmainmax(double a, double b) {return (a >= b ? a : b);}

int main(int argc, char * argv[]) {
		int m; //number of heavy hitters in output
		int counters = 100;
		int threshold = 1000;
		int n = 100000;
		int k; // number of hashes for k-minhash
		int memory;
		double time, walltime;
		double epsil;
		vector<HeavyHitter *> ans;
		int i, j;
		clock_t begint, endt;
		struct timeb begintb, endtb;
		unsigned int **data;
		unsigned int *item;
		FILE * fp = NULL; //the file for the output

		if(argc != 6){
			printf("usage: ./HHH N counters (M) threshold k (minhash) output_file < input_file\n");
			return 0;
		}
		if (argc > 1) n = atoi(argv[1]);
		if (argc > 2) counters = atoi(argv[2]);
		if (argc > 3) threshold = atoi(argv[3]);
		if (argc > 4) k = atoi(argv[4]);
		if (argc > 5) fp = fopen(argv[5], "w");

		if(n/counters >= threshold) {
			printf("Unacceptable parameters: eps*n >= theshold\n");
			return 0;
		}

                data = (unsigned int **) malloc(sizeof(unsigned int) * n * k);

                //init((double)1/(double)counters);

		unsigned int w;
                for (i = 0; i < n; i++) {
                        item = (unsigned int *)malloc(sizeof(unsigned int) * k);
                        for(j=0; j<k; j++){
                                scanf("%d", &w);
                                item[j] = w;
                        }
                        data[i] = item;
                }

/*
                for(int i=0; i<n; i++){
                        for(int j=0; j<k; j++){
                                cout<<i<<" "<<j<<" "<<data[i][j]<<endl;
                        }
                }

*/

		double epsilon = (double)1/(double)counters;
		HHH h(epsilon, k);
		begint = clock();
		ftime(&begintb);
		for (i = 0; i < n; i++) {
			h.update(data[i], 1, k, 0, k, i);
		}
		endt = clock();
		ftime(&endtb);

		time = ((double)(endt-begint))/CLK_PER_SEC;
		walltime = ((double) (endtb.time-begintb.time))+((double)endtb.millitm-(double)begintb.millitm)/1000;
		//memory=maxmemusage();

		//data.clear();
		free(data);

		//h.printCounters(k);

		h.output(threshold, &m, k, ans);

		printf("%d HHHs\n", m);

		map<int, bool> hhrows;
		map<int, bool>::iterator mit;
		
		if (fp != NULL) {
			fprintf(fp, "argv[0] n counters threshold m time memory\n");
			fprintf(fp, "%s %d %d %d %d %lf\n", argv[0], n, counters, threshold, m, time);
			for (i = 0; i < m; i++) {
				HeavyHitter *ptr = ans[i];
				fprintf(fp, "%d ",ptr->mask);
				for(j=0; j<k; j++){
					fprintf(fp, "%d ",ptr->item[j]);
				}	
				fprintf(fp, "%d %d",ptr->lower, ptr->upper);
				//fprintf(fp, " rowlength %d ", ptr->rrows.size());
				fprintf(fp, " rowlength %d ", ptr->rowlength);
				fprintf(fp, " rows -- ");

                                LCLRows *rrows = ptr->rrows;
                                while(rrows != NULL){
					mit = hhrows.find(rrows->row);
					if(mit == hhrows.end()){
                                        	fprintf(fp,"%d ",rrows->row);
						hhrows[rrows->row] = true;
					} 
                                        rrows = rrows->next;
                                }
                                fprintf(fp,"\n");

				//vector<int> rows = ptr->rrows;
				//for(unsigned int pp=0; pp<rows.size(); pp++){
				//	fprintf(fp,"%d ",rows[pp]);
				//}
				//fprintf(fp, "\n");
				
			}
			fclose(fp);
		}
		
		epsil = -1;
		for (i = 0; i < m; i++) {
			HeavyHitter *ptr = ans[i];
			epsil = dblmainmax(epsil, ((double)(ptr->upper-ptr->lower))/n);
		}
		for(i=0; i<m; i++){
			delete ans[i];
		}

	
	return 0;
}


