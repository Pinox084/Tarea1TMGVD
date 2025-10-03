#include <bits/stdc++.h>
#include "kmersheader.hpp"
#include "murmurhash32.hpp"

using namespace std;
using u64 = uint64_t;

// ===================== Count-Min with Conservative Update =====================
class CountMinCU {
public:
    size_t width;
    size_t rows;
    vector<uint32_t> seeds;
    vector<vector<uint32_t>> table;

    CountMinCU(size_t _width, size_t _rows, uint32_t seed0 = 0xBEEF) : width(_width), rows(_rows) {
        table.assign(rows, vector<uint32_t>(width, 0));
        seeds.resize(rows);
        for (size_t i = 0; i < rows; ++i)
            seeds[i] = seed0 + 0x9e3779b9u * (uint32_t)i;
    }

    void insert(u64 kmer_bits) {
        vector<size_t> idx(rows);
        for (size_t r = 0; r < rows; ++r) {
            uint32_t h = murmurhash(&kmer_bits, seeds[r]);
            idx[r] = (size_t)(h % (uint32_t)width);
        }
        uint32_t minv = UINT32_MAX;
        for (size_t r = 0; r < rows; ++r)
            minv = min(minv, table[r][idx[r]]);
        for (size_t r = 0; r < rows; ++r) {
            if (table[r][idx[r]] == minv) table[r][idx[r]]++;
        }
    }

    uint32_t estimate(u64 kmer_bits) const {
        uint32_t ans = UINT32_MAX;
        for (size_t r = 0; r < rows; ++r) {
            uint32_t h = murmurhash(&kmer_bits, seeds[r]);
            size_t idx = (size_t)(h % (uint32_t)width);
            ans = min(ans, table[r][idx]);
        }
        return ans == UINT32_MAX ? 0 : ans;
    }
};

// ===================== Tower Sketch =====================
class TowerSketch {
public:
    vector<CountMinCU> levels;

    TowerSketch(const vector<size_t>& widths,
                const vector<size_t>& rows_per_level,
                uint32_t seed0 = 0xBEEF) {
        if (widths.size() != rows_per_level.size())
            throw runtime_error("widths y rows_per_level deben coincidir");
        levels.reserve(widths.size());
        for (size_t i = 0; i < widths.size(); ++i)
            levels.emplace_back(widths[i], rows_per_level[i], seed0 + (uint32_t)i);
    }

    void insert(u64 canon_kmer_bits) {
        for (auto &lvl : levels) lvl.insert(canon_kmer_bits);
    }

    uint64_t estimate(u64 canon_kmer_bits) const {
        uint64_t best = UINT64_MAX;
        for (const auto &lvl : levels) {
            uint64_t e = lvl.estimate(canon_kmer_bits);
            best = min(best, e);
        }
        return best == UINT64_MAX ? 0 : best;
    }
};

// ===================== FUNCIÓN: Calibración =====================
void calibrar(const string& listaArchivos, int k, const string& outFile, double phi = 0.000002) {
    // ------------------ Configuraciones TowerSketch ------------------
    vector<vector<size_t>> configs_widths = {
        {1<<20, 1<<19, 1<<18, 1<<17, 1<<16}, 
        {1<<22, 1<<21, 1<<20, 1<<19, 1<<18}, 
        {1<<23, 1<<22, 1<<21, 1<<20, 1<<19},
        {1<<24, 1<<23, 1<<22, 1<<21, 1<<20, 1<<19, 1<<18}
    };
    vector<vector<size_t>> configs_rows = {
        {5, 5, 4, 4, 3},
        {6, 5, 5, 4, 4},
        {7, 6, 6, 5, 5}, 
        {10, 8, 7, 7, 6, 6, 5}
    };
    

    ofstream out(outFile);
    if (!out) {
        cerr << "Error: no se pudo abrir " << outFile << " para escribir.\n";
        return;
    }
    out << "archivo,config,N,MAE,RMSE,max_abs_error,max_rel_error,Precision,Recall,F1Score,TP,FP,FN\n";

    ifstream in(listaArchivos);
    string filename;
    while (in >> filename) {
        cerr << "[Calibración] Procesando archivo: " << filename << "\n";

        // ------------------ Conteo exacto ------------------
        unordered_map<u64,uint64_t> exact;
        uint64_t N = 0;

        vector<unique_ptr<TowerSketch>> sketches;
        for (size_t i = 0; i < configs_widths.size(); i++) {
            sketches.emplace_back(make_unique<TowerSketch>(configs_widths[i], configs_rows[i]));
        }

        ifstream f(filename);
        if (!f) continue;
        string seq;
        while (getline(f, seq)) {
            for (auto &c : seq) c = toupper(c);
            string clean;
            for (char c : seq)
                if (c=='A'||c=='C'||c=='G'||c=='T')
                    clean.push_back(c);
            if ((int)clean.size()<k) continue;

            for (size_t i=0;i+k<=clean.size();i++){
                string sub = clean.substr(i,k);
                u64 val = encode_kmer(sub);
                u64 canon = canonical_kmer_bits(val,k);
                N++;
                exact[canon]++;
                for (auto &sk: sketches) sk->insert(canon);
            }
        }

        if (N==0) continue;

        uint64_t threshold = (phi>0.0) ? (uint64_t)ceil(phi*(double)N) : 0;

        
        for (size_t cfg=0; cfg<sketches.size(); cfg++){
            double mae=0.0, rmse=0.0;
            uint64_t max_abs_error=0;
            double max_rel_error=0.0;
            size_t cnt=0;
            size_t TP=0, FP=0, FN=0;

            for (auto &p: exact){
                uint64_t est = sketches[cfg]->estimate(p.first);
                double err = (double)est - (double)p.second;
                mae += fabs(err);
                rmse += err*err;
                max_abs_error = max(max_abs_error,(uint64_t)fabs(err));
                if (p.second>0){
                    double rel = fabs(err)/(double)p.second;
                    max_rel_error = max(max_rel_error,rel);
                }

                // Evaluación para Heavy-Hitters
                if (threshold>0){
                    bool real_HH = p.second>=threshold;
                    bool est_HH  = est>=threshold;
                    if (real_HH && est_HH) TP++;
                    else if (!real_HH && est_HH) FP++;
                    else if (real_HH && !est_HH) FN++;
                }

                cnt++;
            }

            mae /= (double)cnt;
            rmse = sqrt(rmse/(double)cnt);
            double precision = (TP+FP>0)? (double)TP/(TP+FP) : 0.0;
            double recall    = (TP+FN>0)? (double)TP/(TP+FN) : 0.0;
            double F1        = (precision+recall>0)? 2.0*precision*recall/(precision+recall) : 0.0;

            cerr << "  Config " << cfg 
                 << " N=" << N 
                 << " MAE=" << mae 
                 << " RMSE=" << rmse
                 << " MaxAbs=" << max_abs_error
                 << " MaxRel=" << max_rel_error
                 << " Precision=" << precision
                 << " Recall=" << recall
                 << " F1=" << F1 << "\n";

            out << filename << ",cfg" << cfg << "," << N << ","
                << mae << "," << rmse << ","
                << max_abs_error << "," << max_rel_error << ","
                << precision << "," << recall << "," << F1 << ","
                << TP << "," << FP << "," << FN << "\n";
        }
    }

    out.close();
    cerr << "Resultados finales en: " << outFile << "\n";
}

// ===================== FUNCIÓN: Heavy Hitters =====================
void heavy_hitters(const string& listaArchivos, int k, double phi, const string& outFile) {
    vector<size_t> widths = {1<<22, 1<<21, 1<<20, 1<<19, 1<<18}; //conjunto de mejor rendimiento ya evaluado
    vector<size_t> rows   = {6, 6, 5, 5, 4};

    ofstream out(outFile);
    if (!out) {
        cerr << "Error: no se pudo abrir " << outFile << " para escribir.\n";
        return;
    }

    ifstream in(listaArchivos);
    string filename;
    while (in >> filename) {
        cerr << "Procesando archivo: " << filename << "\n";

        TowerSketch ts(widths, rows);
        uint64_t N = 0;
        unordered_set<u64> seen;

        ifstream f(filename);
        if (!f) continue;

        string seq;
        while (getline(f, seq)) {
            for (auto &c : seq) c = toupper(c);
            string clean;
            for (char c : seq)
                if (c == 'A' || c == 'C' || c == 'G' || c == 'T')
                    clean.push_back(c);
            if ((int)clean.size() < k) continue;

            for (size_t i = 0; i + k <= clean.size(); i++) {
                string sub = clean.substr(i, k);
                u64 val = encode_kmer(sub);
                u64 canon = canonical_kmer_bits(val, k);
                N++;
                ts.insert(canon);
                seen.insert(canon);
            }
        }

        uint64_t threshold = (uint64_t)ceil(phi * (double)N);
        cerr << "  Total k-mers válidos: " << N << "\n";
        cerr << "  Umbral absoluto: " << threshold << "\n";

        out << "### Archivo: " << filename << " ###\n";
        size_t bytes_written = 0;
        const size_t MAX_BYTES = 1 * 1024 * 1024;

        for (u64 canon : seen) {
            uint64_t est = ts.estimate(canon);
            if (est >= threshold) {
                string line = decode_kmer_digits(canon, k) + " " + to_string(est) + "\n";
                if (bytes_written + line.size() > MAX_BYTES) break;
                out << line;
                bytes_written += line.size();
            }
        }
        out << "\n";
        cerr << "  HHs guardados (" << bytes_written / 1024.0 << " KB)\n";
    }

    out.close();
    cerr << "Resultados HH en: " << outFile << "\n";
}

// ===================== MAIN =====================
int main(int argc, char* argv[]) {
    
    if (argc == 4) {
        // Modo calibración
        string listaArchivos = argv[1];
        int k = stoi(argv[2]);
        string outFile = argv[3];
        calibrar(listaArchivos, k, outFile);
    }
    else if (argc == 5) {
        // Modo heavy hitters
        string listaArchivos = argv[1];
        int k = stoi(argv[2]);
        double phi = stod(argv[3]);
        string outFile = argv[4];
        heavy_hitters(listaArchivos, k, phi, outFile);
    }
    else {
        cerr << "Uso calibración: " << argv[0] << " lista_archivos.txt k salida.csv\n";
        cerr << "Uso HH:          " << argv[0] << " lista_archivos.txt k phi salida.txt\n";
        return 1;
    }
    return 0;
}