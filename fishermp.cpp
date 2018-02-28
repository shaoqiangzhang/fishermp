/*================================================================================
FisherMP: Motif Prediction using Fisher's exact test with openMP Parallel design
Author: Shaoqiang Zhang
Email: zhangshaoqiang@tjnu.edu.cn
Tianjin Normal University, Tianjin 300387, China
Date: 2018-02-27
==================================================================================*/
#include<map>
#include<string>
#include<iostream>
#include<vector>
#include<fstream>
#include<sstream>
#include<cmath>
#include<cctype>//for the "toupper" function
#include<ctime>
#include<cstdlib>
#include<algorithm>//for the "random_shuffle" functions
#include<numeric>//for the "accumulate" function
#include<omp.h> //openMP
#include<iomanip>
//==================================================================================
using namespace std;
typedef vector<vector<char> > MatrixChar;
typedef vector<vector<string> > MatrixStr;
typedef vector<vector<int> > MatrixInt;
typedef vector<vector<double> > MatrixDouble;
typedef pair<string, double> sdPAIR;
typedef map<string,int> hashmap;
typedef map<string,double> double_hash;
typedef vector<double_hash> sdhashvect;
typedef map<string,vector<int> > vecthash;
typedef vector<vecthash> vhashvect;
typedef pair<vector<string>,double> VectDoublePair;
typedef pair<vector<string>,int> strVectIntPair;
typedef map<vector<string>, int> VectIntHash;
typedef map<vector<string>,double> strVectDoubleHash;
//===================================================================================
string complementary(vector<char> charVect);//complementary DNA strand of a sequence vector
string strReverseComplement(string str);
double log_combin(int n, int m);
vector<int> vector_merge_unique(vector<int> v1, vector<int> v2);
struct CmpByValue {
	bool operator()(const sdPAIR& lhs, const sdPAIR& rhs) {
		return lhs.second < rhs.second;
	}
};
struct CmpVDPairByValue {
	bool operator()(const VectDoublePair& lhs, const VectDoublePair& rhs) {
		return lhs.second < rhs.second;
	}
};
//===================================================================================
int main(int argc, const char** argv){
	int MinSiteSize=5;//motif min length
	int MaxSiteSize=10;//motif max length
	int MotifNumber=10;//number of motifs to find
	int ThreadNum=MaxSiteSize-MinSiteSize+1;// number of threads
	int exist_backfile=0;//if exist background fasta file or not
	ifstream backfile;// define the background/negative file
	if(argc!=2 && argc!=4 && argc!=6 && argc!=8 && argc!=10 && argc!=12){
		cout<<"\n******\nUSAGE:\n******\n";
		cout<<argv[0]<<" <dataset> [OPTIONS]  > OutputFile\n";
		cout<<"\n<dataset>\tfile containing DNA sequences in FASTA format\n";
		cout<<"\nOPTIONS:\n";
		cout<<"-b\t\ta background data file in FASTA format(default=directly produced by the program itself)\n";
		cout<<"-m\t\tminimum size of binding sites to find(default="<<MinSiteSize<<")\n";
		cout<<"-M\t\tMaximum size of binding sites to find(default="<<MaxSiteSize<<")\n";
		cout<<"-n\t\tNumber of motifs to find(default="<<MotifNumber<<")\n";
		cout<<"-t\t\tnumber of threads to call(default="<<ThreadNum<<" i.e. =M-m+1)\n";
		cout<<"\n*******\nDesigned by Shaoqiang Zhang (zhangshaoqiang@tjnu.edu.cn), Feb/26/2018\n\n";
		cout<<"For example: \n"<<argv[0]<<" Klf1.fna -b Klf1.negative.seq -m 5 -M 9 -n 5 -t 10  > Klf1.motifs\n\n";
		cout<<endl;
		exit(1);
	}
	ifstream fastafile(argv[1]);
	for(int i=2;i<argc-1;i=i+2){
		string strindex(argv[i]);
		if(strindex=="-b"){
			backfile.open(argv[i+1],ios::out);
			exist_backfile=1;
		}else if(strindex=="-m"){
			string sminsitesize(argv[i+1]);
			istringstream isminsitesize(sminsitesize);
			isminsitesize>>MinSiteSize;
			if(MinSiteSize<3){
				cout<<"ERROR: the wrong parameter of '-m', please input a positive integer(>=3)."<<endl; exit(1);
			}
		}else if(strindex=="-M"){
			string smaxsitesize(argv[i+1]);
			istringstream ismaxsitesize(smaxsitesize);
			ismaxsitesize>>MaxSiteSize;
			if(MaxSiteSize<3){
				cout<<"ERROR: the wrong parameter of '-M', please input a positive integer(>='-m')."<<endl; exit(1);
			}
		}else if(strindex=="-n"){
			string smotifnum(argv[i+1]);
			istringstream ismotifnum(smotifnum);
			ismotifnum>>MotifNumber;
			if(MotifNumber<=0){
				cout<<"ERROR: the wrong parameter of '-n', please input a positive ineger."<<endl; exit(1);
			}
		}else if(strindex=="-t"){
			string sthreadnum(argv[i+1]);
			istringstream isthreadnum(sthreadnum);
			isthreadnum>>ThreadNum;
			if(ThreadNum<=0){
				cout<<"ERROR: the wrong parameter of '-t', please input a positive ineger."<<endl; exit(1);
			}
		}else{
			cout<<"ERROR: the wrong settings of parameters!"<<endl; exit(1);
		}
	}
	cout<<"The command line: ";
	for(int j=0;j<argc;j++){
		cout<<argv[j]<<" ";
	}
	cout<<"\n\n";
//=============================end of input command line========================================================

	time_t tstart,tend;
	time(&tstart);//start to keep running time
	long Acount=0; long Ccount=0; long Gcount=0; long Tcount=0;
	double Apercent, Cpercent, Gpercent, Tpercent;
	omp_set_num_threads(ThreadNum);

//===============begin===read fasta file into allSeqMat=========use the first thread===============================
	int line_count; char ch;
	double p_value;
	vector<string> seq_annotation;
	vector<string> back_annotation;
	vector<char> seqVector;
	MatrixChar allSeqMat;//store fasta file to a matrix
	MatrixChar backSeqMat;//store background file to a matrix
	/*if(ThreadNum==1){
		omp_set_num_threads(ThreadNum);
	}else{
		omp_set_num_threads(2);
	}*/
	//#pragma omp parallel sections default(shared) private(line_count, ch, seqVector)
	//{//begin parallel sections
		//#pragma omp section
		//{
			line_count=0;
			for(string s;getline(fastafile,s);){
				istringstream sin(s);
				sin>>ch;
				if(ch =='>'){
					line_count++;
					seq_annotation.push_back(s);
					if(!seqVector.empty()){
						allSeqMat.push_back(seqVector);
						seqVector.clear();
					}
				}else{
					if(ch>='a' && ch<='z'){
						ch=toupper(ch);//change letter into capital uppercase
					}
					if(ch =='A'){Acount++;}
					if(ch =='C'){Ccount++;}
					if(ch =='G'){Gcount++;}
					if(ch =='T'){Tcount++;}
					seqVector.push_back(ch);
					for(char base;sin>>base;){
						if(base>='a' && base<='z'){
							base=toupper(base);//change letter into capital
						}
						if(base =='A'){Acount++;}
						if(base =='C'){Ccount++;}
						if(base =='G'){Gcount++;}
						if(base =='T'){Tcount++;}
						seqVector.push_back(base);
					}
				}
			}
			if(line_count==0){
				cout<<"******\nERROR: the input file \"";
				cout<<argv[1]<<"\" isn't FASTA file!!\n******"<<endl;
				exit(1);
			}
			if(!seqVector.empty()){
				allSeqMat.push_back(seqVector);
			}
		//}
			//==================end=====read fasta file into allSeqMat========



//===begin===read the background file after the "-b " parameter into backSeqMat=====use the second thread===
		//#pragma omp section
		//{
			if(exist_backfile){//read the background file
				line_count=0;
				for(string bs;getline(backfile,bs);){
					istringstream sbin(bs);
					sbin>>ch;
					if(ch =='>'){
						line_count++;
						back_annotation.push_back(bs);
						if(!seqVector.empty()){
							backSeqMat.push_back(seqVector);
							seqVector.clear();
						}
					}else{
						if(ch>='a' && ch<='z')
							ch=toupper(ch);//change letter into capital uppercase
						seqVector.push_back(ch);
						for(char bbase;sbin>>bbase;){
							if(bbase>='a' && bbase<='z')
								bbase=toupper(bbase);//change letter into capital
							seqVector.push_back(bbase);
						}
					}
				}
				if(line_count==0){
					cout<<"******\nERROR: the input "<<backfile<<" file isn't FASTA file!!\n******"<<endl;
					exit(1);
				}
				if(!seqVector.empty()){
					backSeqMat.push_back(seqVector);
				}
			}
		//}
//===end==read the background file=====================
	//}//end parallel sections
	
	
	int total_ACGT=Acount+Ccount+Gcount+Tcount;
	//cout<<Acount<<" "<<Ccount<<" "<<Gcount<<" "<<Tcount<<endl;
	Apercent=static_cast<double>(Acount)/static_cast<double>(total_ACGT);
	Cpercent=static_cast<double>(Ccount)/static_cast<double>(total_ACGT);
	Gpercent=static_cast<double>(Gcount)/static_cast<double>(total_ACGT);
	Tpercent=static_cast<double>(Tcount)/static_cast<double>(total_ACGT);
	cout<<setprecision(3)<<"A="<<Apercent<<"\nC="<<Cpercent<<"\nG="<<Gpercent<<"\nT="<<Tpercent<<"\n"<<endl;
	vector<char> ACGT; ACGT.push_back('A');ACGT.push_back('C');ACGT.push_back('G');ACGT.push_back('T');
	vector<double> DNApercent; 
	DNApercent.push_back(Apercent);DNApercent.push_back(Cpercent);
	DNApercent.push_back(Gpercent);DNApercent.push_back(Tpercent);


	vector<char> curr_vector;
	double rp;
	line_count=allSeqMat.size();
	if(!exist_backfile){//produce a background matrix if no background file
		srand((unsigned)time(0));
		for(int i=0;i<line_count;++i){
			seqVector.clear();
			for(int x=0;x<allSeqMat[i].size();x++){
				rp=(rand()%100)/(double)101;//produce a random fraction 0<rp<1
				if(rp<=Apercent){
					seqVector.push_back('A');
				}else if(rp<=(Apercent+Cpercent)){
					seqVector.push_back('C');
				}else if(rp<=(Apercent+Cpercent+Gpercent)){
					seqVector.push_back('G');
				}else{
					seqVector.push_back('T');
				}
			}
			if(!seqVector.empty()){
				backSeqMat.push_back(seqVector);
			}
		}
	}

	//============================================================================================

	//read all k-mers into two maps for the input sequence set and the background sequence set
	vecthash kmerVectMap;//store k-mer to a hash with vector values
	vecthash backKmerVectMap;
	vhashvect kmervectmapvector;//******************************
	vhashvect kvector_private;
	vhashvect backkmervectmapvector;//******************************
	vhashvect backkvector_private;
	
	int kRange=MaxSiteSize-MinSiteSize+1;
	int bgk;
	vector<int> lineLabelVect;
	hashmap one_line_kmer_hash;
	omp_set_num_threads(ThreadNum);

	#pragma omp parallel for private(one_line_kmer_hash,bgk,lineLabelVect,backKmerVectMap,kmerVectMap, kvector_private,backkvector_private) \
	shared(kRange,kmervectmapvector, backkmervectmapvector,MinSiteSize,MaxSiteSize, allSeqMat,backSeqMat)
	for(int k=MinSiteSize;k<=MaxSiteSize+kRange;k++){
		if(k<=MaxSiteSize){
			for(int i=0;i<allSeqMat.size();i++){
				one_line_kmer_hash.clear();
				for(int j=0;j<=(allSeqMat[i].size()-k);j++){
					string kmerstr(allSeqMat[i].begin()+j, allSeqMat[i].begin()+j+k);
					string::size_type idx=kmerstr.find('N');
					if(idx== string::npos ){//if no 'N' in the string
						vector<char> kmer_vector(allSeqMat[i].begin()+j, allSeqMat[i].begin()+j+k);
						string rc_kmerstr=complementary(kmer_vector);
						if(one_line_kmer_hash.count(kmerstr)<1 && one_line_kmer_hash.count(rc_kmerstr)<1){
							one_line_kmer_hash.insert(hashmap::value_type(kmerstr,1));
							if(kmerstr!=rc_kmerstr){
								one_line_kmer_hash.insert(hashmap::value_type(rc_kmerstr,1));
							}
							if(kmerVectMap.count(kmerstr)<1 && kmerVectMap.count(rc_kmerstr)<1){
								lineLabelVect.clear();
								lineLabelVect.push_back(i);
								kmerVectMap.insert(vecthash::value_type(kmerstr,lineLabelVect));
							}else if(kmerVectMap.count(kmerstr)>0){
								kmerVectMap[kmerstr].push_back(i);
							}else if(kmerVectMap.count(rc_kmerstr)>0){
								kmerVectMap[rc_kmerstr].push_back(i);
							}
						}
					}
				}
			}
			kvector_private.push_back(kmerVectMap);
			#pragma omp critical (a)
			kmervectmapvector.insert(kmervectmapvector.end(),kvector_private.begin(),kvector_private.end());
		}
		if(k>MaxSiteSize){
			bgk=k-kRange;
			for(int i=0; i<backSeqMat.size();i++){
				one_line_kmer_hash.clear();
				for(int j=0;j<=(backSeqMat[i].size()-bgk);j++){
					string kmerstr(backSeqMat[i].begin()+j, backSeqMat[i].begin()+j+bgk);
					string::size_type idx=kmerstr.find('N');
					if(idx== string::npos ){
						vector<char> kmer_vector(backSeqMat[i].begin()+j, backSeqMat[i].begin()+j+bgk);
						string rc_kmerstr=complementary(kmer_vector);
						if(one_line_kmer_hash.count(kmerstr)<1 && one_line_kmer_hash.count(rc_kmerstr)<1){
							one_line_kmer_hash.insert(hashmap::value_type(kmerstr,1));
							if(kmerstr!=rc_kmerstr){
								one_line_kmer_hash.insert(hashmap::value_type(rc_kmerstr,1));
							}
							if(backKmerVectMap.count(kmerstr)<1 && backKmerVectMap.count(rc_kmerstr)<1){
								lineLabelVect.clear();
								lineLabelVect.push_back(i);
								backKmerVectMap.insert(vecthash::value_type(kmerstr,lineLabelVect));
							}else if(backKmerVectMap.count(kmerstr)>0){
								backKmerVectMap[kmerstr].push_back(i);
							}else if(backKmerVectMap.count(rc_kmerstr)>0){
								backKmerVectMap[rc_kmerstr].push_back(i);
							}
						}
					}
				}
			}
			backkvector_private.push_back(backKmerVectMap);
			#pragma omp critical (b)
			backkmervectmapvector.insert(backkmervectmapvector.end(),backkvector_private.begin(),backkvector_private.end());
		}		
	}//==========================end of reading k-mer into map=================================
	
	line_count=allSeqMat.size();
	int back_line_count=backSeqMat.size();
	allSeqMat.clear(); backSeqMat.clear();//clear two matrices allSeqMat,backSeqMat to free the memory
	kmerVectMap.clear();
	backKmerVectMap.clear();
	for(int i=0;i<kmervectmapvector.size();i++){
		vecthash kvm=kmervectmapvector[i];
		for(vecthash::const_iterator v=kvm.begin();v!=kvm.end();v++){
			kmerVectMap.insert(vecthash::value_type(v->first,v->second));
		}
	}
	kmervectmapvector.clear();
	for(int j=0;j<backkmervectmapvector.size();j++){
		vecthash bkvm=backkmervectmapvector[j];
		for(vecthash::const_iterator bv=bkvm.begin();bv!=bkvm.end();bv++){
			backKmerVectMap.insert(vecthash::value_type(bv->first,bv->second));
		}
	}
	backkmervectmapvector.clear();
	

	//================the following is to calcluate p-value for each k-mer======================
	
	vector<string> strVect;//store all k-mers into a vector
	MatrixInt strInLinesMatrix;//store sequence index into a matrix
	for (vecthash::const_iterator vh=kmerVectMap.begin(); vh!=kmerVectMap.end(); vh++){
		strVect.push_back(vh->first);
		strInLinesMatrix.push_back(vh->second);
	}
	int a,b,c,d;
	double_hash kmerPvalueHash;
	sdhashvect kmerphashvector;
	sdhashvect kpvector_private;
	int strvsize=strVect.size();
	omp_set_num_threads(ThreadNum);
	#pragma omp parallel default(shared) private(kmerPvalueHash,i,p_value, a, b, c, d, lineLabelVect,kpvector_private)
	{
		int thrid=omp_get_thread_num();
		int allthnum=omp_get_num_threads();
		for(int i=thrid;i<strvsize;i=i+allthnum){//compute p-value for each k-mer
			a=strInLinesMatrix[i].size();
			b=line_count-a;
			c=0;

			string kmer=strVect[i];
			string rcKmer=strReverseComplement(kmer);

			if(backKmerVectMap.count(kmer)>0){
				vecthash::iterator ch=backKmerVectMap.find(kmer);
				lineLabelVect=ch->second;
				c=lineLabelVect.size();
			}
			if(backKmerVectMap.count(rcKmer)>0){
				vecthash::iterator ch=backKmerVectMap.find(rcKmer);
				lineLabelVect=ch->second;
				c=lineLabelVect.size();
			}
			d=back_line_count-c;
			p_value=exp(log_combin(a+c, a)+log_combin(b+d,b)-log_combin(a+b+c+d,a+b));
			
			if((a>c) && (p_value<0.05)) //only keep the k-mers with #foreground>#background
				kmerPvalueHash.insert(double_hash::value_type(kmer,p_value));
		}
		kpvector_private.push_back(kmerPvalueHash);
		#pragma omp critical (p)
		kmerphashvector.insert(kmerphashvector.end(),kpvector_private.begin(),kpvector_private.end());
	}
	
//=========end of computing p-values==========================
	kmerPvalueHash.clear();
	for(int k=0;k<kmerphashvector.size();k++){
		double_hash kphash=kmerphashvector[k];
		for(double_hash::const_iterator ks=kphash.begin();ks!=kphash.end();ks++){
			kmerPvalueHash.insert(double_hash::value_type(ks->first, ks->second));
		}
	}//store k-mers and p-values into hash
	
	vector<sdPAIR> sdVector;
	for(double_hash::const_iterator ivm=kmerPvalueHash.begin();ivm!=kmerPvalueHash.end();ivm++){
		sdPAIR sdv=make_pair(ivm->first,ivm->second);
		sdVector.push_back(sdv);
	}//put these <kmer,pvalue> pairs into a vector
	
	sort(sdVector.begin(), sdVector.end(), CmpByValue());//sort k-mer with p-value from small to large
	
	/*for (int i=0;i<100;i++){
		cout<<sdVector[i].first<<"\t"<<sdVector[i].second<<"\n";
	}*/
	
	
//===The following is to form first-hand motif by computing difference between kmers======

	MatrixStr motif_matrix;
	int sdv=sdVector.size();
	hashmap used_kmer_hash;//define the kmers that have been used in the front motifs
	
	for (int i=0;i<100;i++){
		string currKmer=sdVector[i].first;
		if(used_kmer_hash.count(currKmer)<1){
			vector<string> oneMotifKmerVect;
			oneMotifKmerVect.clear();
			oneMotifKmerVect.push_back(currKmer);
			int strLength=currKmer.size();

			vector<char> currkmervect;
			currkmervect.clear();
			currkmervect.resize(strLength);
			currkmervect.assign(currKmer.begin(),currKmer.end());

			vector<int> identical_posit_vector;
			for(int p=0;p<strLength;p++){
				identical_posit_vector.push_back(0);//identical=0; different=1;
			}
			vector<int> bool_posit_vect=identical_posit_vector;
			vector<string> leftKmerVect;
			leftKmerVect.clear();

			for(int j=i+1;j<sdv;j++){//find exactly one position different with currKmer
				string nextKmer=sdVector[j].first;
				if(strLength==nextKmer.size()){
					vector<char> nextkmervect;
					nextkmervect.resize(strLength);
					nextkmervect.assign(nextKmer.begin(),nextKmer.end());
					string nextRCkmer=complementary(nextkmervect);
					vector<char> rckmervect;
					rckmervect.clear();
					rckmervect.resize(strLength);
					rckmervect.assign(nextRCkmer.begin(),nextRCkmer.end());

					int diff_position_count=0;
					int rc_diff_pos_count=0;
					int diff_pos=0; int rc_diff_pos=0;

					for(int k=0;k<strLength;k++){
						if(currkmervect[k]!=nextkmervect[k]){
							diff_position_count++;//hamming distance
							diff_pos=k;
						}
						if(currkmervect[k]!=rckmervect[k]){
							rc_diff_pos_count++;//hamming distance
							rc_diff_pos=k;
						}
					}

					if(diff_position_count==1){
						identical_posit_vector[diff_pos]=1;//revise value
						int total_diff_positions=accumulate(identical_posit_vector.begin(), identical_posit_vector.end(),0);//sum of vector's elements
						if(strLength>6 && total_diff_positions<=floor(static_cast<double>(strLength)*0.3)){
							oneMotifKmerVect.push_back(nextKmer);
							if(j<100){
								used_kmer_hash.insert(hashmap::value_type(nextKmer,1));
							}
							bool_posit_vect[diff_pos]=1;
						}else if(strLength<=6 && total_diff_positions==1){
							oneMotifKmerVect.push_back(nextKmer);
							if(j<100){
								used_kmer_hash.insert(hashmap::value_type(nextKmer,1));
							}
						}else{
							leftKmerVect.push_back(nextKmer);
							break;
						}
					}else if(rc_diff_pos_count==1){
						identical_posit_vector[rc_diff_pos]=1;
						int total_diff_positions=accumulate(identical_posit_vector.begin(),identical_posit_vector.end(),0);//sum of vector's elements
						if(strLength>6 && total_diff_positions<=floor(static_cast<double>(strLength)*0.3)){
							oneMotifKmerVect.push_back(nextRCkmer);
							if(j<100){
								used_kmer_hash.insert(hashmap::value_type(nextKmer,1));
							}
							bool_posit_vect[rc_diff_pos]=1;
						}else if(strLength<=6 && total_diff_positions==1){
							oneMotifKmerVect.push_back(nextRCkmer);
							if(j<100){
								used_kmer_hash.insert(hashmap::value_type(nextKmer,1));
							}
						}else{
							leftKmerVect.push_back(nextKmer);
							break;
						}
					}else{
						leftKmerVect.push_back(nextKmer);
					}
					oneMotifKmerVect.erase(unique(oneMotifKmerVect.begin(), oneMotifKmerVect.end()), oneMotifKmerVect.end());//unique the elements
				}
			}//=========find one-position non-conserved with the current first kmer

			leftKmerVect.erase(unique(leftKmerVect.begin(), leftKmerVect.end()),leftKmerVect.end());
			vector<string> tempMtfVect=oneMotifKmerVect;

			if(strLength>6){
				for(int j=0;j<leftKmerVect.size();j++){
					string nextKmer=leftKmerVect[j];
					vector<char> nextkmervect;
					nextkmervect.resize(strLength);
					nextkmervect.assign(nextKmer.begin(),nextKmer.end());
					string nextRCkmer=complementary(nextkmervect);
					vector<char> rckmervect;
					rckmervect.clear();
					rckmervect.resize(strLength);
					rckmervect.assign(nextRCkmer.begin(),nextRCkmer.end());
					for(int l=0;l<tempMtfVect.size();l++){
						vector<char> currkmer_vect;
						currkmer_vect.clear();
						currkmer_vect.resize(strLength);
						currkmer_vect.assign(tempMtfVect[l].begin(),tempMtfVect[l].end());

						int diff_position_count=0;
						int rc_diff_pos_count=0;
						int diff_pos=0; int rc_diff_pos=0;

						for(int k=0;k<strLength;k++){
							if(currkmer_vect[k]!=nextkmervect[k]){
								diff_position_count++;//hamming distance
								diff_pos=k;
							}
							if(currkmer_vect[k]!=rckmervect[k]){
								rc_diff_pos_count++;//hamming distance
								rc_diff_pos=k;
							}
						}

						if(diff_position_count==1){
							bool_posit_vect[diff_pos]=1;//revise value
							int total_diff_positions=accumulate(bool_posit_vect.begin(), bool_posit_vect.end(),0);//sum of vector's elements
							if(total_diff_positions<=floor(static_cast<double>(strLength)*0.3)){
								oneMotifKmerVect.push_back(nextKmer);
							}else{
								break;
							}
						}else if(rc_diff_pos_count==1){
							bool_posit_vect[rc_diff_pos]=1;
							int total_diff_positions=accumulate(bool_posit_vect.begin(),bool_posit_vect.end(),0);//sum of vector's elements
							if(total_diff_positions<=floor(static_cast<double>(strLength)*0.3)){
								oneMotifKmerVect.push_back(nextRCkmer);
							}else{
								break;
							}
						}	
					}
					oneMotifKmerVect.erase(unique(oneMotifKmerVect.begin(), oneMotifKmerVect.end()), oneMotifKmerVect.end());//unique the elements
				}
			}
			motif_matrix.push_back(oneMotifKmerVect);
		}
	}

//=================end of first-hand motif forming ==============================

	/*for(int i=0;i<motif_matrix.size();i++){
		for(int j=0;j<motif_matrix[i].size();j++){
			cout<<motif_matrix[i][j]<<"\t";
		}
		cout<<endl;
	}*/

	//---refine motifs by calculating new p-values--------------------------------
	vector<VectDoublePair> motifpvalue_pair_vect;
	// store all motifs and their p-values into a vector of (motif,p-value) pairs
	
	strVectDoubleHash motif_pvalue_hash;
	VectIntHash strvect_ma_hash;//store motif and ma
	VectIntHash strvect_mc_hash;//store motif and mc
	for(int i=0;i<motif_matrix.size();i++){
		vector<string> mv=motif_matrix[i];
		string first_kmer=mv[0];
		double_hash::iterator kk=kmerPvalueHash.find(first_kmer);
		double merged_Pvalue=kk->second;
		vector<string> oneMotifKmerVect;
		oneMotifKmerVect.clear();
		oneMotifKmerVect.push_back(first_kmer);
		vecthash::iterator ch=kmerVectMap.find(first_kmer);
		vector<int> mergedLineLabelVect=ch->second;
		string rc_string=strReverseComplement(first_kmer);
		vecthash::iterator bk;
		vector<int> mergedBackLineVect;
		mergedBackLineVect.clear();
		if(backKmerVectMap.count(first_kmer)>0){
			bk=backKmerVectMap.find(first_kmer);
			mergedBackLineVect=bk->second;
		}else if(backKmerVectMap.count(rc_string)>0){
			bk=backKmerVectMap.find(rc_string);
			mergedBackLineVect=bk->second;
		}
		int ma=0;
		int mc=0;
		for(int j=1;j<mv.size();j++){
			string nextKmer=mv[j];
			string nextRCkmer=strReverseComplement(nextKmer);
			if(kmerVectMap.count(nextKmer)>0){
				vecthash::iterator nk=kmerVectMap.find(nextKmer);
				mergedLineLabelVect=vector_merge_unique(mergedLineLabelVect,nk->second);
			}else if(kmerVectMap.count(nextRCkmer)>0){
				vecthash::iterator nk=kmerVectMap.find(nextRCkmer);
				mergedLineLabelVect=vector_merge_unique(mergedLineLabelVect,nk->second);
			}
			if(backKmerVectMap.count(nextKmer)>0){
				vecthash::iterator nk=backKmerVectMap.find(nextKmer);
				mergedBackLineVect=vector_merge_unique(mergedBackLineVect,nk->second);
			}else if(backKmerVectMap.count(nextRCkmer)>0){
				vecthash::iterator nk=backKmerVectMap.find(nextRCkmer);
				mergedBackLineVect=vector_merge_unique(mergedBackLineVect,nk->second);
			}
			//the following is to compute the p-value after merging k-mers
			a=mergedLineLabelVect.size();
			b=line_count-a;
			c=mergedBackLineVect.size();
			d=back_line_count-c;
			p_value=exp(log_combin(a+c, a)+log_combin(b+d,b)-log_combin(a+b+c+d,a+b));
			if(p_value<0.05){
				oneMotifKmerVect.push_back(nextKmer);
				merged_Pvalue=p_value;
				ma=a;
				mc=c;
			}else{
				break;
			}
		}
		motif_pvalue_hash.insert(strVectDoubleHash::value_type(oneMotifKmerVect, merged_Pvalue));
		strvect_ma_hash.insert(VectIntHash::value_type(oneMotifKmerVect, ma));
		strvect_mc_hash.insert(VectIntHash::value_type(oneMotifKmerVect, mc));
	}
	
	for(strVectDoubleHash::const_iterator svp=motif_pvalue_hash.begin();svp!=motif_pvalue_hash.end();svp++){
		VectDoublePair vdpair=make_pair(svp->first,svp->second);
		motifpvalue_pair_vect.push_back(vdpair);
	}//transfer hash to vector of pairs

	sort(motifpvalue_pair_vect.begin(), motifpvalue_pair_vect.end(), CmpVDPairByValue());
	//sort motifs with p-value from small to large
	
//=========================================================================================================
	for(int mm=0;mm<MotifNumber;mm++)
	{//----output the top number of motifs-------------------------------------
		vector<string> oneMotifKmerVect=motifpvalue_pair_vect[mm].first;
		VectIntHash::iterator via=strvect_ma_hash.find(oneMotifKmerVect);
		int ma=via->second;
		VectIntHash::iterator vic=strvect_mc_hash.find(oneMotifKmerVect);
		int mc=vic->second;
		cout<<"\nMOTIF TOP "<<mm+1<<" details:\n";
		cout<<"Postives    Negatives    P-value\n";
		cout<< ma<<"/"<<line_count<<"    "<<mc<<"/"<<back_line_count<<"    "<<motifpvalue_pair_vect[mm].second<<"\n";
		cout<<"Enriched Matching Words:\n";
		cout<<"Word      RC_word      Postives    Negatives    P-value\n";
		int totalaa=0;
		int mtflen=oneMotifKmerVect[0].size();
		MatrixInt PFM; PFM.clear();
		MatrixDouble PWM; PWM.clear();
		for(int qq=0;qq<4;qq++){
			vector<int> oneline;
			for(int pp=0;pp<mtflen;pp++){
				oneline.push_back(0);
			}
			PFM.push_back(oneline);
		}
		vector<int> seqNoVect;
		seqNoVect.clear();//
		for(int nn=0;nn<oneMotifKmerVect.size();nn++ ){
			string oneword=oneMotifKmerVect[nn];
			string rc_word=strReverseComplement(oneword);
			int aa=0; int cc=0;
			if(kmerVectMap.count(oneword)>0){
				vecthash::iterator ch=kmerVectMap.find(oneword);
				lineLabelVect=ch->second;
				seqNoVect.insert(seqNoVect.end(),lineLabelVect.begin(),lineLabelVect.end());//
				aa=lineLabelVect.size();
				totalaa=totalaa+aa;
				vector<char> charVect;
				charVect.resize(oneword.size());
				charVect.assign(oneword.begin(),oneword.end());
				for(int i=0;i<charVect.size();i++){
					if(charVect[i] == 'A'){
						PFM[0][i]=PFM[0][i]+aa;
					}
					if(charVect[i] == 'C'){
						PFM[1][i]=PFM[1][i]+aa;
					}
					if(charVect[i] == 'G'){
						PFM[2][i]=PFM[2][i]+aa;
					}
					if(charVect[i] == 'T'){
						PFM[3][i]=PFM[3][i]+aa;
					}
				}

			}
			if(kmerVectMap.count(rc_word)>0){
				vecthash::iterator ch=kmerVectMap.find(rc_word);
				lineLabelVect=ch->second;
				aa=lineLabelVect.size();
				totalaa=totalaa+aa;
				vector<char> charVect;
				charVect.resize(oneword.size());
				charVect.assign(oneword.begin(),oneword.end());
				for(int i=0;i<charVect.size();i++){
					if(charVect[i] == 'A'){
						PFM[0][i]=PFM[0][i]+aa;
					}
					if(charVect[i] == 'C'){
						PFM[1][i]=PFM[1][i]+aa;
					}
					if(charVect[i] == 'G'){
						PFM[2][i]=PFM[2][i]+aa;
					}
					if(charVect[i] == 'T'){
						PFM[3][i]=PFM[3][i]+aa;
					}
				}
			}
			if(backKmerVectMap.count(oneword)>0){
				vecthash::iterator ch=backKmerVectMap.find(oneword);
				lineLabelVect=ch->second;
				cc=lineLabelVect.size();
			}
			if(backKmerVectMap.count(rc_word)>0){
				vecthash::iterator ch=backKmerVectMap.find(rc_word);
				lineLabelVect=ch->second;
				cc=lineLabelVect.size();
			}
			if(kmerPvalueHash.count(oneword)>0){
				double_hash::iterator wd=kmerPvalueHash.find(oneword);
				p_value=wd->second;
			}
			if(kmerPvalueHash.count(rc_word)>0){
				double_hash::iterator wd=kmerPvalueHash.find(rc_word);
				p_value=wd->second;
			}
			cout<<oneword<<"    "<<rc_word<<"    "<<aa<<"/"<<line_count<<"    "<<cc<<"/"<<back_line_count<<"    "<<p_value<<"\n";
		}
		cout<<"These sites belong to the following sequence numbers: ";
		sort(seqNoVect.begin(),seqNoVect.end());
		seqNoVect.erase(unique(seqNoVect.begin(), seqNoVect.end()), seqNoVect.end());
		for (int i=0;i<seqNoVect.size();i++){
			if(i%10==0){
				cout<<"\n";
			}
			cout<<seqNoVect[i]+1<<" ";

		}
		cout<<endl;
		cout<<"The total number of sites: "<<totalaa<<"\nPosition Weight Matrix\n";
		vector<double> pwm_line;
		double pwm_score;
		for(int dna=0;dna<4;dna++){
			cout<<ACGT[dna];
			pwm_line.clear();
			for(int i=0;i<mtflen;i++){
				if(PFM[dna][i]>0){
					pwm_score=log(static_cast<double>(PFM[dna][i])/(static_cast<double>(totalaa)*DNApercent[dna]))/log(2);
					pwm_line.push_back(pwm_score);
					cout<<setprecision(4)<<" "<<pwm_score;
				}else{
					pwm_score=log(0.01/DNApercent[dna])/log(2);
					pwm_line.push_back(pwm_score);
					cout<<setprecision(4)<<" "<<pwm_score;
				}
			}
			PWM.push_back(pwm_line);
			cout<<endl;
		}
		cout<<"Information contents:\n";
		cout<<"I";
		for(int i=0;i<mtflen;i++){
			double ic=0.0;
			for(int dna=0;dna<4;dna++){
				ic=ic+(static_cast<double>(PFM[dna][i])/(static_cast<double>(totalaa)))*(PWM[dna][i]);
			}
			cout<<setprecision(4)<<" "<<ic;
		}
		cout<<endl;
		
		cout<<"Position Frequency Matrix:\n";
		cout<<"a";
		for(int i=0;i<mtflen;i++){
			cout<<" "<<PFM[0][i];
		}
		cout<<endl;
		cout<<"c";
		for(int i=0;i<mtflen;i++){
			cout<<" "<<PFM[1][i];
		}
		cout<<endl;
		cout<<"g";
		for(int i=0;i<mtflen;i++){
			cout<<" "<<PFM[2][i];
		}
		cout<<endl;
		cout<<"t";
		for(int i=0;i<mtflen;i++){
			cout<<" "<<PFM[3][i];
		}
		cout<<endl;

	}


	//=====below: record the running time==============
	time(&tend);
	double dif = difftime (tend,tstart);
	cout<<"\n\nTotal running time: "<<dif<<" seconds\n"<<endl;

}//end of main function

//========================functions in the following===============================================

string complementary(vector<char> charVect)
{//get the reverse complementary string of a sequence
	reverse(charVect.begin(),charVect.end());
	vector<char> complementVect;
	for(int j=0;j<charVect.size();j++){
		if(charVect[j]=='A'){
			complementVect.push_back('T');
		}else if(charVect[j]=='C'){
			complementVect.push_back('G');
		}else if(charVect[j]=='G'){
			complementVect.push_back('C');
		}else if(charVect[j]=='T'){
			complementVect.push_back('A');
		}else{
			complementVect.push_back('N');
		}
	}
	string rc_str(complementVect.begin(), complementVect.end());
	return rc_str;
}

string strReverseComplement(string str)
{//get the reverse complementary string of a sequence
	vector<char> charVect;
	charVect.resize(str.size());
	charVect.assign(str.begin(),str.end());
	reverse(charVect.begin(),charVect.end());
	vector<char> complementVect;
	for(int j=0;j<charVect.size();j++){
		if(charVect[j]=='A'){
			complementVect.push_back('T');
		}else if(charVect[j]=='C'){
			complementVect.push_back('G');
		}else if(charVect[j]=='G'){
			complementVect.push_back('C');
		}else if(charVect[j]=='T'){
			complementVect.push_back('A');
		}else{
			complementVect.push_back('N');
		}
	}
	string rc_str(complementVect.begin(), complementVect.end());
	return rc_str;
}

double log_combin(int n, int m)
{//calculate log combination C(n,m) using stirling formula
	double result;
	if(m >= n){
		result=0.0;
	}else if(m==0){
		result=0.0;
	}else{
		double N=static_cast<double>(n);
		double M=static_cast<double>(m);
		result=(N+0.5)*log(N)-0.5*log(2*3.14159265)-(M+0.5)*log(M)-(N-M+0.5)*log(N-M);
	}
	return result;
}

vector<int> vector_merge_unique(vector<int> v1, vector<int> v2)
{//merge two vectors v1 and v2, then delete the repeating elements
	vector<int> v3;
	v3.insert(v3.end(), v1.begin(), v1.end());
	v3.insert(v3.end(), v2.begin(), v2.end());
	sort(v3.begin(),v3.end());
	v3.erase(unique(v3.begin(), v3.end()), v3.end());
	return v3;
}

//===================end functions========================================
