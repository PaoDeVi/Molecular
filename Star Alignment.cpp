#include <algorithm>
#include <cctype>
#include <complex>
#include <fstream>
#include <iostream>
#include <list>
#include <map>
#include <math.h>
#include <queue>
#include <set>
#include <stack>
#include <stdio.h>
#include <string.h>
#include <string>
#include <time.h>
#include <vector>

#define N 10001

using namespace std;

int match = 1, mismatch_ = 1, gap_score = 2;
string sequencias[N];
int sum_score_index[N];
int dp[N][N]; 

int value(char a, char b){
  if(a==b){
    return match;
  }
  else{
    return -mismatch_;
  }
}

int needleman_wunsch(string Seq1, string Seq2,int n,int m) {
  for (int i = 0; i <= n; i++) dp[i][0] = dp[0][i] = -i * gap_score;
  for (int i = 1; i <= n; i++) {
    for (int j = 1; j <= m; j++) {
      int S = value(Seq1[i - 1], Seq2[j - 1]) ;
      dp[i][j] = max(dp[i - 1][j - 1] + S,
                     max(dp[i - 1][j] - gap_score, dp[i][j - 1] - gap_score));
    }
  }
  return dp[n][m];
}

pair<string, string> best_alignment(string Seq1, string Seq2,int n,int m) {
  string retSeq1, retSeq2;
  stack<char> SSeq1, SSeq2;
  int ii = n, jj = m;
  while (ii != 0 || jj != 0) {
    if (ii == 0) {
      SSeq1.push('-');
      SSeq2.push(Seq2[jj - 1]);
      jj--;
    } 
    else if (jj == 0) {
      SSeq1.push(Seq1[ii - 1]);
      SSeq2.push('-');
      ii--;
    } 
    else {
      int S = value(Seq1[ii - 1],Seq2[jj - 1]);
      if (dp[ii][jj] == dp[ii - 1][jj - 1] + S) {
        SSeq1.push(Seq1[ii - 1]);
        SSeq2.push(Seq2[jj - 1]);
        ii--;
        jj--;
      } else if (dp[ii - 1][jj] > dp[ii][jj - 1]) {
        SSeq1.push(Seq1[ii - 1]);
        SSeq2.push('-');
        ii--;
      } else {
        SSeq1.push('-');
        SSeq2.push(Seq2[jj - 1]);
        jj--;
      }
    }
  }

  while (!SSeq1.empty()) {
    retSeq1 += SSeq1.top();
    SSeq1.pop();
    retSeq2 += SSeq2.top();
    SSeq2.pop();
  }
  return make_pair(retSeq1, retSeq2);
}

void all_alignments(int i, int j, string align1, string align2,vector<pair<string, string>> &alignments,string Seq1, string Seq2) {
  if (i == 0 && j == 0) {
    reverse(align1.begin(), align1.end());
    reverse(align2.begin(), align2.end());
    alignments.push_back(make_pair(align1, align2));
    return;
  }

  if (i > 0 && j > 0) {
    int S = value(Seq1[i - 1], Seq2[j - 1]) ;
    if (dp[i][j] == dp[i - 1][j - 1] + S) {
      all_alignments(i - 1, j - 1, align1 + Seq1[i - 1], align2 + Seq2[j - 1],
                         alignments, Seq1,  Seq2);
    }
  }

  if (i > 0 && dp[i][j] == dp[i - 1][j] - gap_score) {
    all_alignments(i - 1, j, align1 + Seq1[i - 1], align2 + '-', alignments, Seq1,  Seq2);
  }

  if (j > 0 && dp[i][j] == dp[i][j - 1] - gap_score) {
    all_alignments(i, j - 1, align1 + '-', align2 + Seq2[j - 1], alignments, Seq1,  Seq2);
  }
}

void star_alignment(int num_seq){
  int max_score_index = 0, max_size = -1;
  string alignments[num_seq];

  //Selecciona aquella secuencia con mejor score
  for(int i = 0; i < num_seq; i++){
     sum_score_index[i] = 0;
    for(int j = 0; j < num_seq; j++){
      if(j!=i){
        sum_score_index[i] += needleman_wunsch(sequencias[i],sequencias[j],sequencias[i].length(),sequencias[j].length());
    }
    }
    if(i>0 && sum_score_index[i] > sum_score_index[max_score_index]){
      max_score_index = i;
    }
  }
  alignments[max_score_index] = sequencias[max_score_index];
  pair<string, string> optimal;
  cout<<"La secuencia "<<max_score_index+1<<" sera el centro de alineacion"<<endl;
  for(int i = 0; i < num_seq; i++){
    if(i!=max_score_index){
      needleman_wunsch(sequencias[max_score_index],sequencias[i],sequencias[max_score_index].length(),sequencias[i].length());
      vector<pair<string, string>> alignments_;
      all_alignments(sequencias[max_score_index].length(),sequencias[i].length(), "", "", alignments_,sequencias[max_score_index],sequencias[i]);
        for(int j = i; j >= 0; j--){
           int max_al = -999999;
           for (const auto &align : alignments_) {
          int score = needleman_wunsch(sequencias[j],align.second,sequencias[j].length(),align.second.length());
          if(score > max_al) alignments[i] = align.second;
        }
    }
  }
  }

  printf("El mejor score alineamiento es: %d \n",sum_score_index[max_score_index]);
  
  for(int i = 0; i < num_seq; i++){
    int n = alignments[i].length();
    if( n > max_size){
      max_size = n;
    }
  }
  
  for(int i = 0; i < num_seq; i++){
    int n = alignments[i].length();
    while(n < max_size){
      alignments[i] += '-';
      n = alignments[i].length();
    }
    cout<<alignments[i]<<endl;
  }
}

vector<string> obtenerSecuenciasADN(const string& nombreArchivo) {
    vector<string> secuencias;
    ifstream archivo(nombreArchivo);
    string linea;

    while (getline(archivo, linea)) {
        size_t posInicio = linea.find("F: 5'-") + 6;
        size_t posFin = linea.find("'", posInicio);
        if (posInicio != string::npos && posFin != string::npos) {
            string secuencia = linea.substr(posInicio+2, posFin - posInicio - 4);
            if(secuencia != "  R:") secuencias.push_back(secuencia);
        }
        posInicio = linea.find("R: 5'-") + 6;
        posFin = linea.find("'", posInicio);
        if (posInicio != string::npos && posFin != string::npos) {
            string secuencia = linea.substr(posInicio+2, posFin - posInicio - 4);
            if(secuencia != "  F:") secuencias.push_back(secuencia);
        }
    }

    archivo.close();
    return secuencias;
}

int front_sequences(vector<string> secuencias){
  int i = 0, index = 0;
  for (const auto& secuencia : secuencias) {
    if(i%2==0){
      sequencias[index] = secuencia;
      //cout << secuencia << endl;
      index++;
    }
     i++; 
  }
  return index;
}

int reverse_sequences(vector<string> secuencias){
  int i = 0, index = 0;
  for (const auto& secuencia : secuencias) {
    if(i%2==1){
      sequencias[index] = secuencia;
      //cout << secuencia << endl;
      index++;
    }
     i++; 
  }
  return index;
}

int all_sequences(vector<string> secuencias){
  int index = 0;
  for (const auto& secuencia : secuencias) {
      sequencias[index] = secuencia;
      index++;
    }
  return index;
}

int main() {
  //Ejemplo en clase
  int num_seq = 5;
  sequencias[0] = "ATTGCCATT"; 
  sequencias[1] = "ATGGCCATT"; 
  sequencias[2] = "ATCCAATTTT"; 
  sequencias[3] = "ATCTTCTT";
  sequencias[4] = "ACTGACC";
  
  string nombreArchivo = "BRCA1.txt";
  vector<string> secuencias = obtenerSecuenciasADN(nombreArchivo);
  //Secuencias directas
  num_seq = front_sequences(secuencias);
  //Secuencias inversas
  num_seq = reverse_sequences(secuencias);
  //Todos contra todos
  num_seq = all_sequences(secuencias);
  
  star_alignment(num_seq);

  return 0;
}
