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
int n, m; string Seq1, Seq2;

int dp[N][N]; string cadenas[N];

int value(char a, char b){
  if(a==b){
    return match;
  }
  else{
    return -mismatch_;
  }
}

int needleman_wunsch() {
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

pair<string, string> best_alignment() {
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

void all_alignments(int i, int j, string align1, string align2,vector<pair<string, string>> &alignments) {
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
                         alignments);
    }
  }

  if (i > 0 && dp[i][j] == dp[i - 1][j] - gap_score) {
    all_alignments(i - 1, j, align1 + Seq1[i - 1], align2 + '-', alignments);
  }

  if (j > 0 && dp[i][j] == dp[i][j - 1] - gap_score) {
    all_alignments(i, j - 1, align1 + '-', align2 + Seq2[j - 1], alignments);
  }
}

void print_dot_plot(const string &Seq1, const string &Seq2,
                    const vector<pair<int, int>> &matches) {

  int n = Seq1.length();
  int m = Seq2.length();
  vector<vector<char>> dot_plot(n, vector<char>(m, '.'));

  for (auto &match : matches) {
    int i = match.first;
    int j = match.second;
    dot_plot[i][j] = '*';
  }

  printf("   ");
  for (int j = 0; j < m; ++j) {
    printf("%c ", Seq2[j]);
  }
  printf("\n");

  for (int i = 0; i < n; ++i) {
    printf("%c ", Seq1[i]);
    for (int j = 0; j < m; ++j) {
      printf("%c ", dot_plot[i][j]);
    }
    printf("\n");
  }
}

string read(string file) {
  ifstream archivo(file);
  string contenido;
  char caracter;

  if (archivo.is_open()) {
    while (archivo.get(caracter)) {
      if (isalpha(caracter)) {
        contenido += caracter;
      }
    }
    archivo.close();
  } 
  return contenido;
}

void divide(string line) {
  int pos = line.find("Bacteria");
  int pos2 = line.find("SarsCov");
  int pos3 = line.find("Influenza");
  
  cadenas[0] = line.substr(pos + 8, pos2 - 8);
  cadenas[1] = line.substr(pos2 + 7, pos3 - 1075);
  cadenas[2] = line.substr(pos3 + 9);
  //printf("Posiciones %s\n\n %s\n\n %s\n\n", cadenas[0].c_str(),cadenas[1].c_str(), cadenas[2].c_str());
}
int main() {
  string archive;
  archive = read("Sequencias.txt");
  divide(archive);
  Seq1 = cadenas[0];
  //Seq1 = "AAAC";
  n = Seq1.length();
  Seq2 = cadenas [2];
  //Seq2 = "AGC";
  m = Seq2.length();
  printf("El mejor costo obtenido es: %d\n", needleman_wunsch());
//Mejor alineamiento
  pair<string, string> alignment = best_alignment();
  printf("El mejor alineamiento es:\n");
  printf("%s\n%s\n", alignment.first.c_str(), alignment.second.c_str());

  // Obtener todos los alineamientos posibles
  vector<pair<string, string>> alignments;
  all_alignments(n, m, "", "", alignments);

  printf("Todos los alineamientos encontrados:\n");
  for (const auto &align : alignments) {
    printf("%s\n%s\n\n", align.first.c_str(), align.second.c_str());
  }

  vector<pair<int, int>> matches;
  
  string &alignFirst = alignment.first;
  string &alignSecond = alignment.second;
  for (int i = 0; i < alignFirst.length(); ++i) {
    if (alignFirst[i] == alignSecond[i] && alignFirst[i] != '-') {
      matches.push_back(
        make_pair(i, i));
    }
  }
  
  printf("Matriz de Puntos:\n");
  print_dot_plot(Seq1, Seq2, matches);

  return 0;
}

