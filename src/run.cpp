#include "parser.h"
#include "kmacs.h"
#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

int main(int argc, char **argv){
	int c;
	int rmq =0;
	while ((c = getopt (argc, argv, "r:")) != -1)
		switch (c){
		case 'r':
	        rmq = 1;
        	break;
		if (isprint (optopt))
                	fprintf (stderr, "Unknown option `-%c'.\n", optopt);
             	else
               		fprintf (stderr,"Unknown option character `\\x%x'.\n", optopt);
		return 1;
		default:
		abort ();
	}
    	FILE *stream;
	vector<string> taxa;
	vector<string> sequences;

    if (argc < 3) {
        fprintf(stderr, "usage: %s file.fasta k\n", argv[0]);
        return 1;
    }
	stream = fopen(argv[argc-2], "r");
    if (!stream) {
        perror("can't open input file");
	fprintf(stderr, "usage: %s file.fasta k\n", argv[0]);
		return 1;
    }
	int k = atoi(argv[argc-1]);
	double time=0.0;
	double start =clock();
	parse(stream, taxa, sequences);
	cout << "Parsing Done from FASTA file (Success)" << endl << endl;

	double** dmat = new double*[sequences.size()];
		for(int i=0;i<sequences.size();i++){
			dmat[i] = new double[sequences.size()];
	}
	
	calcDmat(dmat,  sequences, k, rmq);
	writeDmat(dmat, taxa);
	
	time += clock() - start;		
	time = time/CLOCKS_PER_SEC;
	cout << "Execution Time - " << time << "s" << endl;
	cout << "Output File - distanceMatrix.txt" << endl;

	for(int i=0;i<sequences.size();i++){
		delete [] dmat[i];
	}
	delete [] dmat;

	return 0;
}
