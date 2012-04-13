/* This is  GPL */
/* Author: Sébastien Boisvert */

#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <map>
#include <stdlib.h>
#include <boost/random.hpp>
#include <boost/random/normal_distribution.hpp>
#include <ctime>

//#define VERBOSE

using namespace std;

char getRandomNucleotide(){
	int a=rand()%4;
	switch (a){
		case 0:
			return 'A';
		case 1:
			return 'T';
		case 2:
			return 'C';
		case 3:
			return 'G';
	}
	return 'N';
}

string addSubstitutionErrors(string*a,double rate){
	string b=*a;
	int resolution=100000000;
	int length=b.length();
	int threshold=(int)rate*resolution;

	for(int i=0;i<length;i++){
		int val=rand()%resolution;
		if(val<threshold){
			char oldNucleotide=b[i];
			char newNucleotide=getRandomNucleotide();
			while(newNucleotide==oldNucleotide){
				newNucleotide=getRandomNucleotide();
			}
			b[i]=newNucleotide;
		}
	}
	return b;
}

string reverseComplement(string*a){
	ostringstream b;
	int n=a->length()-1;
	while(n>=0){
		char symbol=a->at(n--);

		switch (symbol){
			case 'A':
				symbol='T';
				break;
			case 'T':
				symbol='A';
				break;
			case 'C':
				symbol='G';
				break;
			case 'G':
				symbol='C';
				break;
		}
		b<<symbol;
	}
	return b.str();
}

int main(int argc,char**argv){
	cout<<"Welcome to the VirtualNextGenSequencer, a free paired read simulator that generates free reads."<<endl;
	cout<<endl;

	if(argc!=9){

		cout<<argc<<" but needs 9"<<endl;

		cout<<"Usage: "<<endl;
		cout<<argv[0]<<" GENOME.FASTA SUBSTITUTION_RATE AVERAGE_OUTER_DISTANCE STANDARD_DEVIATION PAIRS READ_LENGTH OUT1.fasta OUT2.fasta"<<endl;
		cout<<endl;
		cout<<"Example: "<<endl;
		cout<<" VirtualNextGenSequencer Prevotella_melaninogenica_ATCC_25845_uid51377.fasta 0.0025 400 40 1000 101 joe_1.fasta joe_2.fasta"<<endl;

		return 0;
	}

	cout<<endl;
	cout<<"Parsing command."<<endl;

	string file=argv[1];
	double substitutionRate=atof(argv[2]);
	int average=atoi(argv[3]);
	int dev=atoi(argv[4]);
	int pairs=atoi(argv[5]);
	int length=atoi(argv[6]);
	string left=argv[7];
	string right=argv[8];

	cout<<endl;
	cout<<"<Sequencing simulation parameters>"<<endl;
	cout<<"AverageLength="<<average<<" StandardDeviation="<<dev<<endl;
	cout<<"Pairs="<<pairs<<" / "<<2*pairs<<" sequences"<<endl;
	cout<<"/ ReadLength="<<length<<endl;
	cout<<"SubstitutionRate="<<substitutionRate<<endl;

	int i=0;

	srand(time(NULL)+getpid()+4);
	ofstream leftOut(left.c_str());
	ofstream rightOut(right.c_str());
	cout<<"INPUT: "<<file<<endl;
	cout<<"OUTPUTS: "<<left<<" "<<right<<endl;



	/* the sampler for fragment length (not read length) */
	boost::mt19937 generator(time(NULL)+getpid());
	boost::normal_distribution<> distribution(average,dev);
	boost::variate_generator<boost::mt19937&,boost::normal_distribution<> > fragmentLengthSampler(generator,distribution);


	ostringstream genome;
	string header="";
	ifstream f(file.c_str());

	vector<string> genomeSequences;
	vector<string> headerNames;

	while(!f.eof()){
		char line[1000000];
		int maxSize=1000000;
		f.getline(line,maxSize);

		if(line[0]!='>'){
			genome<<line;
		}else{
			if(header!=""){
				string sequence=genome.str();
				genomeSequences.push_back(sequence);
				genome.str("");
				headerNames.push_back(header);
				header="";
			}

			header=line;
		}
	}
	f.close();

	string sequence=genome.str();
	genomeSequences.push_back(sequence);
	headerNames.push_back(header);

	cout<<endl;
	cout<<"Sequences: "<<genomeSequences.size()<<endl;

	vector<int> lengths;

	uint64_t totalLength=0;

	for(int i=0;i<(int)genomeSequences.size();i++){
		cout<<"  "<<i<<"  ----->  "<<headerNames[i]<<" / "<<genomeSequences[i].length()<<" nucleotides !"<<endl;

		int sequenceLength=genomeSequences[i].length();
		lengths.push_back(sequenceLength);

		totalLength+=sequenceLength;

	}
	cout<<endl;
	cout<<"TotalLength: "<<totalLength<<endl;
	cout<<endl;

	/* the sampler for position placement */
	boost::mt19937 generator2(time(NULL)+getpid()*3);
	boost::uniform_int<uint64_t> uniformDistribution(0,totalLength-1);
	boost::variate_generator<boost::mt19937&,boost::uniform_int<uint64_t> > fragmentPositionSampler(generator2,uniformDistribution);

	cout<<endl;
	cout<<"Commencing..."<<endl;
	cout<<"Sampling will be done on both DNA strands..."<<endl;
	cout<<endl;

	map<int,int> outputs;

	while(i<pairs){
		if(i%100000==0){
			cout<<i+1<<"/"<<pairs<<endl;
		}

		uint64_t globalStart=fragmentPositionSampler();
		int observedDistance=(int)fragmentLengthSampler();

		uint64_t start=-1;
		int sequenceId=-1;
		uint64_t sequenceLength=-1;

		// find the sequence that is hit */
	
		uint64_t sum=0;
		for(int k=0;k<(int)lengths.size();k++){
			sum+=lengths[k];

			if(globalStart<sum){
				start=globalStart-sum+lengths[k];
				sequenceId=k;
				sequenceLength=lengths[k];
				break;
			}
		}

		if(start==-1){
			cout<<"FATAL ERROR"<<endl;
		}

		/* check that it fits on the thing. */

		if(start+observedDistance>sequenceLength){
			continue;
		}

		int reverse=rand()%2;

		#ifdef VERBOSE
		cout<<"Sequence: "<<sequenceId<<" Length: "<<lengths[sequenceId]<<endl;
		cout<<"Global start= "<<globalStart<<" start= "<<start<<" length= "<<length<<endl;
		#endif

		string leftSequence=genomeSequences[sequenceId].substr(start,length);
		string rightSequence=genomeSequences[sequenceId].substr(start+observedDistance-length,length);

		if(reverse==1){
			rightSequence=reverseComplement(&rightSequence);
		}else{
			leftSequence=reverseComplement(&leftSequence);
		}

		leftSequence=addSubstitutionErrors(&leftSequence,substitutionRate);
		rightSequence=addSubstitutionErrors(&rightSequence,substitutionRate);

		leftOut<<">"<<i<<"/1"<<endl<<leftSequence<<endl;
		rightOut<<">"<<i<<"/2"<<endl<<rightSequence<<endl;

		outputs[sequenceId]++;
		outputs[sequenceId]++;

		i++;
	}

	cout<<i<<"/"<<pairs<<endl;
	cout<<endl;
	cout<<"Done."<<endl;
	cout<<endl;

	leftOut.close();
	rightOut.close();

	for(map<int,int>::iterator i=outputs.begin();i!=outputs.end();i++){
		cout<<" Sequence "<<headerNames[i->first]<<" -------------> "<<i->second<<" sequences written.."<<endl;
	}
	
	cout<<endl;
	cout<<"Have a nice day !"<<endl;
	
	return EXIT_SUCCESS;
}
