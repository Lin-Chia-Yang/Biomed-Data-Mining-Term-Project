#include<iostream>
#include<fstream>
#include<string>
#include<vector>
#include<climits>
#include<math.h>
#include <omp.h>
using namespace std;
int find_motif_width(string str);
string convert(int num);
int hamdist(string motif, string DNA);
bool same_mutate_number(int mutate_number[],int size, int length, int mutate);

int find_motif_width(string str){
    int index=str.find_last_of(" ");
    str=str.substr(index+1);
    return stoi(str);
}

string convert(int num, int motif_width){
    string bits="";
    while(true){
        int remain = num & 3;
        bits = to_string(remain) + bits;
        num = num >> 2;
        if(num==0){
            break;
        }
    }
    while(bits.length()<motif_width){
        bits='0'+bits;
    }
    for(int i=0;i<bits.length();i++){
        if(bits[i]=='0'){
            bits[i]='A';
        }
        if(bits[i]=='1'){
            bits[i]='G';
        }
        if(bits[i]=='2'){
            bits[i]='C';
        }
        if(bits[i]=='3'){
            bits[i]='T';
        }
    }
    return bits;
}

int hamdist(string motif, string DNA){
    int score=0;
    for(int i=0;i<motif.length();i++){
        if(motif[i]!=DNA[i]){
            score += 1;
        }
    }
    return score;
}

bool same_mutate_number(int mutate_number[],int size, int length, int mutate){
    for(int i=0;i<size;i++){
        for(int j=0;j<length;j++){
            if(mutate_number[i * length + j] == mutate){
                break;
            }
            else if(j==length-1){
                return false;
            }
        }
    }
    return true;
}

int main(int argc, char *argv[]){
    if(argc != 2){
        cout << "error: need input file" << endl;
        return 0;
    }
    string input;
    vector<string> sequence;
    ifstream myfile;
    myfile.open(argv[1]);
    int motif_width=0;
    while(getline(myfile, input)){
        if(input[0]=='m'){
            motif_width=find_motif_width(input);
            //cout << "motif width is " << motif_width << endl;
        }
        if(input[0]=='>'){
            getline(myfile, input);
            sequence.push_back(input);
        }
    }
    int sequence_length=sequence[0].length();
    int bestdist=INT_MAX;
    string key;
    vector<vector<string>> bestword(sequence.size());
    vector<vector<int>> bestposition(sequence.size());
    for(int i=0;i<pow(4,motif_width);i++){
        string motif=convert(i, motif_width);
        //cout << motif <<endl;
        int mutate_number[sequence.size() * (sequence_length-motif_width)];
        #pragma omp parallel for num_threads(12)
        for(int j=0;j<sequence.size();j++){
            for(int k=0;k<sequence_length-motif_width;k++){
                mutate_number[j * (sequence_length-motif_width) + k] = hamdist(motif, sequence[j].substr(k, motif_width));
            }
        }
        for(int l=0;l<motif_width;l++){
            if(same_mutate_number(mutate_number, sequence.size(), sequence_length-motif_width, l) == true){
                if(l < bestdist){
                    key = motif;
                    bestdist = l;
                    for(int m=0;m<sequence.size();m++){
                        bestword[m].clear();
                        bestposition[m].clear();
                        for(int n=0;n<sequence_length-motif_width;n++){
                            if(mutate_number[m * (sequence_length-motif_width) + n] == l){
                                bestword[m].push_back(sequence[m].substr(n, motif_width));
                                bestposition[m].push_back(n+1);
                                n += motif_width - 1;
                            }
                        }
                    }
                }
                break;
            }
        }
    }
    ofstream outputfile;
    outputfile.open("output.txt");
    outputfile << "key = " << key << "(" << motif_width << "," << bestdist << ")"<<endl;
    cout << "key = " << key << "(" << motif_width << "," << bestdist << ")"<<endl;
    for(int i=0;i<bestword.size();i++){
        outputfile << "Sequence" << i+1 << ": ";
        cout << "Sequence" << i+1 << ": ";
        for(int j=0;j<bestword[i].size();j++){
            outputfile << "location: " << bestposition[i][j] << "~" << bestposition[i][j]+motif_width-1 << " ";
            cout << "location: " << bestposition[i][j] << "~" << bestposition[i][j]+motif_width-1 << " ";
            outputfile<< bestword[i][j] << " ";
            cout << bestword[i][j] << " ";
        }
        outputfile << endl;
        cout << endl;
    }
    outputfile.close();
    myfile.close();
    return 0;
}