#include <bits/stdc++.h>
#include <cstdint>
#include "kmersheader.hpp"
#include "murmurhash32.hpp"

using namespace std;
using u64 = uint64_t;

// ====== COUNT SKETCH ======
struct CountSketch {
    int d, w;
    vector<vector<int64_t>> table;
    vector<uint32_t> hashSeeds;
    vector<uint32_t> signSeeds;

    CountSketch(int d_, int w_) : d(d_), w(w_) {
        table.assign(d, vector<int64_t>(w, 0));
        hashSeeds.resize(d);
        signSeeds.resize(d);
        for (int i = 0; i < d; i++) {
            hashSeeds[i] = 1337 + i*17;
            signSeeds[i] = 2023 + i*31;
        }
    }

    int sign(u64 key, int j) const {
        uint32_t h = murmurhash(&key, signSeeds[j]);
        return (h & 1) ? 1 : -1;
    }

    void update(u64 key, int64_t delta) {
        for (int j = 0; j < d; j++) {
            uint32_t h = murmurhash(&key, hashSeeds[j]);
            int col = h % w;
            table[j][col] += delta * sign(key, j);
        }
    }

    int64_t estimate(u64 key) const {
        vector<int64_t> est(d);
        for (int j = 0; j < d; j++) {
            uint32_t h = murmurhash(&key, hashSeeds[j]);
            int col = h % w;
            est[j] = table[j][col] * sign(key, j);
        }
        nth_element(est.begin(), est.begin() + d/2, est.end());
        return est[d/2];
    }
};

void calibrar_countsketch_HH(const string& listaArchivos, int k, double phi, const string& outFile) {
    vector<int> ds = {9, 10, 11, 15, 20};
    vector<int> ws = {10000, 20000, 50000, 100000, 200000, 500000, 1000000, 2000000};

    ofstream out(outFile);
    if (!out) {
        cerr << "Error: no se pudo abrir " << outFile << "\n";
        return;
    }

    // Cabecera una sola vez
    out << "archivo,d,w,size_bytes,MAE,MRE,Precision,Recall,F1\n";

    ifstream in(listaArchivos);
    string filename;
    while (getline(in, filename)) {
        if (filename.empty()) continue;
        filename.erase(remove_if(filename.begin(), filename.end(),
                                 [](unsigned char c){ return isspace(c); }),
                       filename.end());
        if (filename.empty()) continue;

        cerr << "[Calibración] Procesando archivo: " << filename << "\n";

        unordered_map<u64, uint64_t> exact;
        uint64_t N = 0;

        ifstream f(filename);
        if (!f) {
            cerr << "  No pude abrir " << filename << "\n";
            continue;
        }

        string seq;
        while (getline(f, seq)) {
            transform(seq.begin(), seq.end(), seq.begin(), ::toupper);
            if (seq.empty() || seq[0] == '>') continue;

            string clean;
            for (char c : seq) if (c=='A'||c=='C'||c=='G'||c=='T') clean.push_back(c);
            if ((int)clean.size() < k) continue;

            for (size_t i = 0; i + k <= clean.size(); i++) {
                string sub = clean.substr(i, k);
                u64 val = encode_kmer(sub);
                u64 canon = canonical_kmer_bits(val, k);
                N++;
                exact[canon]++;
            }
        }

        if (N == 0) {
            cerr << "  No se encontraron k-mers válidos.\n";
            continue;
        }

        uint64_t threshold = (uint64_t)ceil(phi * (double)N);
        unordered_set<u64> gtHH;
        for (auto &p : exact) if (p.second >= threshold) gtHH.insert(p.first);

        for (int d : ds) {
            for (int w : ws) {
                CountSketch cs(d, w);
                for (auto &p : exact) cs.update(p.first, p.second);

                double mae = 0.0, mre = 0.0;
                size_t cnt = exact.size();
                for (auto &p : exact) {
                    int64_t est = cs.estimate(p.first);
                    double err = (double)est - (double)p.second;
                    mae += fabs(err);
                    if (p.second > 0) mre += fabs(err)/(double)p.second;
                }
                mae /= cnt;
                mre /= cnt;

                unordered_set<u64> csHH;
                for (auto &p : exact) {
                    int64_t est = cs.estimate(p.first);
                    if (max<int64_t>(0, est) >= threshold) csHH.insert(p.first);
                }

                size_t tp = 0;
                for (auto &hh : csHH) if (gtHH.count(hh)) tp++;
                size_t fp = csHH.size() - tp;
                size_t fn = gtHH.size() - tp;

                double precision = (csHH.empty() ? 0.0 : (double)tp / (double)(tp + fp));
                double recall = (gtHH.empty() ? 0.0 : (double)tp / (double)(tp + fn));
                double f1 = (precision + recall > 0 ? 2.0*precision*recall/(precision+recall) : 0.0);

                size_t sketch_size = cs.d * cs.w * sizeof(int64_t);

                out << filename << "," << d << "," << w << "," << sketch_size << ","
                    << mae << "," << mre << "," << precision << "," << recall << "," << f1 << "\n";
            }
        }
    }

    out.close();
    cerr << "Resultados acumulados en: " << outFile << "\n";
}



// ===================== MAIN =====================
int main(int argc, char* argv[]) {
    if (argc == 5) {
        string listaArchivos = argv[1];
        int k = stoi(argv[2]);
        double phi = stod(argv[3]);
        string outFile = argv[4];
        calibrar_countsketch_HH(listaArchivos, k, phi, outFile);
    } else {
        cerr << "Uso: " << argv[0] << " lista_archivos.txt k phi salida.csv\n";
        return 1;
    }
    return 0;
}