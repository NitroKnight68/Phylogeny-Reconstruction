#include "textParser.h"

// Parses the FASTA file into usable data
int parse(FILE *stream, vector<string>& taxa, vector<string>& sequences){
	vector<char> seqBuffer;
	vector<char> taxaBuffer;
	char c=getc(stream);
	while(c!=EOF){
		if(c=='>'){
			while((c=getc(stream))!='\n'){
				taxaBuffer.push_back(c);		
			}
		}
		else if(c==';')
			while((c=getc(stream))!='\n');
		else{
			while((c=getc(stream))!='>' && c!=EOF){
				if(!isalpha(c)){
					continue;
				}
				c=toupper(c);
				seqBuffer.push_back(c);	
			}	
		}
		if(seqBuffer.size()>0){
			string taxaString(taxaBuffer.begin(),taxaBuffer.end());
			string seqString(seqBuffer.begin(),seqBuffer.end());
			taxa.push_back(taxaString);
			sequences.push_back(seqString);
			seqBuffer.clear();
			taxaBuffer.clear();
		}
	}
	return 0;
}

// Writes the Distance Matrix to a TXT file
int writeDmat(double** DMat, vector<string>& taxa){
	ofstream outfile;
	outfile.open("distanceMatrix.txt");
		for(int i=0;i<taxa.size();i++){
				for(int j=0;j<taxa.size();j++){
					if(i!=j){
						outfile << DMat[j][i] << "  ";
						}
					else
						outfile << "0" << "  ";
				}
			outfile << endl;
		}
	outfile.close();
	return 0;
}
