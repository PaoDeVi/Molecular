#include <iostream>
#include <vector>
#include <string>
#include <limits.h>
#include <algorithm>
#include <cmath>

using namespace std;

// Valores de emparejamiento
const int au = -4;
const int cg = -5;
const int gu = -1;
const int other = 0;

//Par
struct Pairing {
    int first;
    int second;
    bool emparejado;
};

//Alpha
int alpha_fun(char a, char b) {
    if ((a == 'A' && b == 'U') || (a == 'U' && b == 'A')) return au;
    else if ((a == 'C' && b == 'G') || (a == 'G' && b == 'C')) return cg;
    else if ((a == 'G' && b == 'U') || (a == 'U' && b == 'G')) return gu;
    else return other;
}

// Función para llenar la matriz de programación dinámica
void dp_fun(string seq, vector<vector<int>>& dp) {
    int n = seq.size();
    for (int l = 1; l < n; ++l) {
        for (int i = 0; i < n - l; ++i) {
            int j = i + l;
            dp[i][j] = dp[i][j - 1];
            for (int k = i; k < j; ++k) {
                int alpha = alpha_fun(seq[k], seq[j]);
                if (alpha < 0) {
                    int score = (k > 0 ? dp[i][k - 1] : 0) + (k + 1 <= j - 1 ? dp[k + 1][j - 1] : 0) + alpha;
                    dp[i][j] = min(dp[i][j], score);
                }
            }
        }
    }
}

void build_pairs(string seq, vector<vector<int>>& dp, vector<Pairing>& pairs, int i, int j,int threshold) {
    if (i >= j) return;
    
    if (dp[i][j] == dp[i][j - 1]) {
        build_pairs(seq, dp, pairs, i, j - 1,threshold);
    } 
    else {
        bool aux = false;
        for (int k = i; k < j; ++k) {
            int alpha = alpha_fun(seq[k], seq[j]);
            if (alpha < 0) {
                int score = (k > 0 ? dp[i][k - 1] : 0) + (k + 1 <= j - 1 ? dp[k + 1][j - 1] : 0) + alpha;
                if (dp[i][j] == score && abs(k-j)>threshold) {
                    aux = true;
                    pairs.push_back({k, j, true});
                    build_pairs(seq, dp, pairs, i, k - 1,threshold);
                    build_pairs(seq, dp, pairs, k + 1, j - 1,threshold);
                    break;
                }
            }
        }
    }
}

// Función para imprimir los pares emparejados en un formato similar a print_ARN
void print_ARN(vector<Pairing> pairs, string seq) {
    if (pairs.empty()) return;
    int last_index = pairs.size() - 1;
    //Linea 1
    cout << seq[pairs[0].first];
    
    for (int i = 1; i < pairs.size(); i++){
       cout << " - " << seq[pairs[i].first];
    }
       
    if (!pairs[last_index].emparejado) cout << " - ";
    cout << endl; 
    //Linea 2
    
    if (!pairs[0].emparejado) cout << " ";
    else cout << "|";
    
    for (int i = 1; i < pairs.size(); i++){
        cout << (pairs[i].emparejado ? "   |" : "    ");
    }
    
    if (!pairs[last_index].emparejado) cout << "  |";
    cout << endl;
    
    //Linea 3
    if (!pairs[0].emparejado) cout << "-";
    else cout << seq[pairs[0].second];
    
    for (int i = 1; i < pairs.size(); i++){
        if(pairs[i].second != pairs[i].first){
            cout << " - " << seq[pairs[i].second];
        }
        else{
            cout << " - " << "-";
        }
        
    }
    
    if (!pairs[last_index].emparejado) cout << " - ";
    cout << endl;
}

// Función para ordenar los pares por el campo `first`
void sortPairs(vector<Pairing>& pairs) {
    sort(pairs.begin(), pairs.end(), [](const Pairing& a, const Pairing& b) {
        return a.first < b.first;
    });
}

int main() {
    string seq;
    cout << "Ingrese la secuencia de ARN: "; cin >> seq;

    int n = seq.size();
    vector<vector<int>> dp(n, vector<int>(n, 0));
    vector<Pairing> pairs;

    // Llenar la matriz de DP
    dp_fun(seq, dp);

    cout << "Matriz de DP:\n";
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            cout << dp[i][j] << " ";
        }
        cout << "\n";
    }

    //Armar los pares
    build_pairs(seq, dp, pairs, 0, n - 1,2);
    
    sortPairs(pairs);
    cout << "Pares emparejados:\n";
    for (const auto& p : pairs) {
        cout << "(" << p.first << ", " << p.second << ") - (" << seq[p.first] << ", " << seq[p.second] << ")\n";
    }
    
    //Incluir los indices que no emparejaron
    for(int i = 0; i<n;i++){
        bool inpair = 0;
        for (const auto& p : pairs) {
            if(i == p.first || i == p.second) inpair = 1;
        }

        if(!inpair){
            pairs.push_back({i,i, false});
        }
    }
    sortPairs(pairs);
    print_ARN(pairs,seq);
    return 0;
}
