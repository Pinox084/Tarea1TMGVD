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
            hashSeeds[i] = 1337 + i * 17;
            signSeeds[i] = 2023 + i * 31;
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
        nth_element(est.begin(), est.begin() + d / 2, est.end());
        return est[d / 2];
    }
};

// ===================== HEAVY HITTERS CON COUNT SKETCH =====================
void countsketch_HH(const string &listaArchivos, int k, double phi, const string &outFile) {
    const int d = 11;
    const int w = 2000000;

    ofstream out(outFile);
    if (!out) {
        cerr << "Error: no se pudo abrir " << outFile << " para escribir.\n";
        return;
    }

    ifstream in(listaArchivos);
    string filename;
    while (getline(in, filename)) {
        filename.erase(remove_if(filename.begin(), filename.end(),
                                 [](unsigned char c) { return isspace(c); }),
                       filename.end());
        if (filename.empty()) continue;

        cerr << "[HH] Procesando archivo: " << filename << "\n";

        ifstream f(filename);
        if (!f) {
            cerr << "  (se omite) no pude abrir " << filename << "\n";
            continue;
        }

        CountSketch cs(d, w);
        uint64_t N = 0;
        unordered_set<u64> seen;

        string line;
        while (getline(f, line)) {
            transform(line.begin(), line.end(), line.begin(), ::toupper);
            if (line.empty() || line[0] == '>') continue;

            string clean;
            for (char c : line)
                if (c == 'A' || c == 'C' || c == 'G' || c == 'T')
                    clean.push_back(c);

            if ((int)clean.size() < k) continue;

            for (size_t i = 0; i + k <= clean.size(); i++) {
                string sub = clean.substr(i, k);

                bool valid = true;
                for (char c : sub) {
                    if (c != 'A' && c != 'C' && c != 'G' && c != 'T') {
                        valid = false;
                        break;
                    }
                }
                if (!valid) continue;

                u64 val = encode_kmer(sub);
                if (val == UINT64_MAX) continue;

                u64 canon = canonical_kmer_bits(val, k);
                if (canon == UINT64_MAX) continue;

                cs.update(canon, 1);
                N++;
                seen.insert(canon);
            }
        }

        if (N == 0) {
            cerr << "  No se encontraron k-mers válidos.\n";
            continue;
        }

        uint64_t threshold = (uint64_t)ceil(phi * (double)N);
        cerr << "  Total k-mers: " << N << " | Umbral: " << threshold << "\n";

        out << "### Archivo: " << filename << " (N=" << N << ")\n";

        size_t bytes_written = 0;
        const size_t MAX_BYTES = 1 * 1024 * 1024; // 1 MB por archivo

        for (u64 canon : seen) {
            if (canon == UINT64_MAX) continue;

            int64_t est = cs.estimate(canon);
            if (est < 0) continue;  // ignorar estimaciones negativas
            if ((uint64_t)est == UINT64_MAX) continue; // descartar basura

            if ((uint64_t)est >= threshold) {
                string line_out = decode_kmer_digits(canon, k) + " " + to_string(est) + "\n";
                if (bytes_written + line_out.size() > MAX_BYTES) {
                    out << "# --- Límite 1MB alcanzado, truncando ---\n";
                    break;
                }
                out << line_out;
                bytes_written += line_out.size();
            }
        }
        out << "\n";
        cerr << "  HHs escritos: " << (bytes_written / 1024.0) << " KB\n";
    }

    out.close();
    cerr << "Resultados HH en: " << outFile << "\n";
}

// ===================== MAIN =====================
int main(int argc, char *argv[]) {
    if (argc == 5) {
        string listaArchivos = argv[1];
        int k = stoi(argv[2]);
        double phi = stod(argv[3]);
        string outFile = argv[4];
        countsketch_HH(listaArchivos, k, phi, outFile);
    } else {
        cerr << "Uso: " << argv[0] << " lista_archivos.txt k phi salida.txt\n";
        return 1;
    }
    return 0;
}
