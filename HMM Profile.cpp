#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <stdlib.h> 

using namespace std;

vector<string> secuencias;
string aminolist = "STQYGAVLICMPFWDEKRH";
string background;
int DT = 0;

void print(){
    cout << "Secuencias" << endl;
    for (int i = 0; i < secuencias.size(); i++) {
        cout << secuencias[i] << endl;
    }
    cout << endl;
}

void delete_col(){
    float cant_secuencias = secuencias.size();
    for (int i = 0; i < secuencias[0].size();) {
        float aux = 0; 
        for (int j = 0; j < cant_secuencias; j++) {
            if (secuencias[j][i] != '-') aux++;
        }

        if (aux < cant_secuencias / 2) {
            for (int j = 0; j < cant_secuencias; j++) {
                secuencias[j].erase(i, 1);
            }
        }
        else i++; 
    }
}

class Estado{
public:
    int transciones[3] = { 0, 0, 0 };
    Estado* estados[3];
    
    void init(int prob1,int prob2,int prob3,Estado* est_1,Estado* est_2,Estado* est_3) {
        transciones[0] = prob1; estados[0] = est_1;
        transciones[1] = prob1+prob2; estados[1] = est_2;
        transciones[2] = prob1+prob2+prob3; estados[2] = est_3;
    }
};

class Match : public Estado{
public:
    vector<pair<int, char>> profileAcumulado;
    void calcularProfile(int col) {
        map<char, int> profile;
        for (int i = 0; i < secuencias.size(); i++) {
            char anminoacido = secuencias[i][col];
            if (anminoacido != '-') {
                if (profile.find(anminoacido) == profile.end()) 
                    profile[anminoacido] = 1;
                else
                    profile[anminoacido]++;
            }
        }
        
        for (auto& prof : profile) {
            char tipo = prof.first;
            int cantidad = prof.second;
            profileAcumulado.push_back(pair<int, char>(cantidad, tipo));
        }
    }

    void calcularProfAcum(int col){
        calcularProfile(col);
        for (int i = 1; i < profileAcumulado.size(); i++)
            profileAcumulado[i].first += profileAcumulado[i - 1].first;
    }

    void printProfile(){
        cout << profileAcumulado[0].second <<": "<< profileAcumulado[0].first +1<<"/"<< secuencias.size() + 20 << endl;
        for (int i = 1; i < profileAcumulado.size(); i++)
            cout << profileAcumulado[i].second<<": "<<profileAcumulado[i].first - profileAcumulado[i - 1].first+1<< "/" <<secuencias.size() + 20 << endl;
        cout << "Otros (AddOn) : 1/" << secuencias.size() + 20 << endl;
    }

    void emision(){
        int n = rand() % (20 + secuencias.size());
        if (n < 20) background.push_back(aminolist[n]);
        else {
            for (int i = 0; i < profileAcumulado.size(); i++) {
                if (n <= profileAcumulado[i].first)
                        background.push_back(profileAcumulado[i].second);
            }
        }
    }
};

vector<Match> match;
vector<Estado> insertion;
vector<Estado> deletes;
vector<Estado> aux; 

void calcularProbabilidadesTransiciones(int i, int& M, int& I, int& D) {
    M = I = D = 1; 
    for (int j = 0; j < secuencias.size(); j++) {
        if (secuencias[j][i + 1] == '-') D++;
        else M++;
    }
    cout <<"Match " << M << "/" << secuencias.size() + 3 << endl;
    cout <<"Insertion " << I << "/" << secuencias.size() + 3 << endl;
    cout <<"Delete " << D << "/" << secuencias.size() + 3 << endl;
    cout << endl;
}

void HMM(int col) {
    int p1, p2, p3;

    aux = vector<Estado>(2); 
    match = vector<Match>(col + 3);
    insertion = vector<Estado>(col + 1);
    deletes = vector<Estado>(col + 1);

    cout << "Transiciones " << " : "<<endl;
    calcularProbabilidadesTransiciones(-1, p1, p2, p3);
    aux[0].init(p1, p2, p3, &match[1], &insertion[0], &deletes[1]);
    insertion[0].init(p1, p2, p3, &match[1], &insertion[0], &deletes[1]);

    for (int i = 1; i < col; i++) {
        cout << "Transiciones en la columna " << i << " : "<<endl;
        calcularProbabilidadesTransiciones(i, p1, p2, p3);
        match[i].init(p1, p2, p3, &match[i + 1], &insertion[i], &deletes[i + 1]);
        insertion[i].init(p1, p2, p3, &match[i + 1], &insertion[i], &deletes[i + 1]);
        deletes[i].init(p1, p2, p3, &match[i + 1], &insertion[i], &deletes[i + 1]);

        // Calcula el profile
        match[i].calcularProfAcum(i - 1);
        cout << "Match Columna" << i << " : "<<endl;
        match[i].printProfile();
        cout << endl;
    }

    p1 = DT - 1;
    p2 = 1;

    match[col].init(p1, p2, 0, &aux[1], &insertion[col], nullptr);
    insertion[col].init(p1, p2, 0, &aux[1], &insertion[col], nullptr);
    deletes[col].init(p1, p2, 0, &aux[1], &insertion[col], nullptr);
}

int main() {
    secuencias = {
        "VGA--HAGEY",
        "V----NVDEV",
        "VEA--DVAGH",
        "VKG------D",
        "VYS--TYETS",
        "FNA--NIPKH",
        "IAGADNGAGY",
    };
    delete_col();
    int len = secuencias[0].size();
    cout << "Cantidad de estados Match relevantes: " << len << endl;
    print();
    DT = secuencias.size() + 3;
    HMM(len);
    return 0;
}
