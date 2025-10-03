#include <bits/stdc++.h>
#include "kmersheader.hpp"

using namespace std;
using u64 = uint64_t;

//Void de heavy hitters
void extraerHeavyHittersExactos(
    const string &listaArchivos,
    int k,
    double phi,
    const string &outHHFile
) {
    ifstream listin(listaArchivos);
    if (!listin) {
        cerr << "Error: no se pudo abrir lista de archivos: " << listaArchivos << "\n";
        return;
    }

    ofstream hhOut(outHHFile);
    if (!hhOut) {
        cerr << "Error: no se pudo abrir archivo de heavy hitters: " << outHHFile << "\n";
        return;
    }

    string filename;
    while (listin >> filename) {
        cerr << "Procesando archivo: " << filename << "\n";

        unordered_map<u64, uint64_t> exact;
        exact.reserve(1 << 20);
        uint64_t N = 0;

        ifstream fin(filename);
        if (!fin) {
            cerr << "  Warning: no pude abrir " << filename << " (se omite)\n";
            continue;
        }

        string line;
        while (getline(fin, line)) {
            for (char &c : line) c = toupper((unsigned char)c);
            string clean;
            clean.reserve(line.size());
            for (char c : line)
                if (c == 'A' || c == 'C' || c == 'G' || c == 'T')
                    clean.push_back(c);

            if ((int)clean.size() < k) continue;

            for (size_t i = 0; i + k <= clean.size(); ++i) {
                string sub = clean.substr(i, k);
                u64 val = encode_kmer(sub);
                if (val == UINT64_MAX) continue;
                u64 canon = canonical_kmer_bits(val, k);
                exact[canon]++;
                N++;
            }
        }
        fin.close();

        if (N == 0) {
            cerr << "  No se encontraron k-mers válidos en " << filename << "\n";
            hhOut << "### Archivo: " << filename << " ###\n\n";
            continue;
        }

        uint64_t threshold = (uint64_t)ceil(phi * (double)N);
        cerr << "  Total k-mers válidos: " << N << " | Umbral: " << threshold << "\n";

        vector<pair<u64, uint64_t>> HHs;
        HHs.reserve(128);
        for (const auto &p : exact) {
            if (p.second >= threshold) HHs.emplace_back(p.first, p.second);
        }

        hhOut << "### Archivo: " << filename << " ###\n";
        size_t bytes_written = 0;
        const size_t MAX_BYTES = 1 * 1024 * 1024; // 1 MB

        for (const auto &p : HHs) {
            string line_out = decode_kmer_digits(p.first, k) + " " + to_string(p.second) + "\n";
            if (bytes_written + line_out.size() > MAX_BYTES) {
                hhOut << "# --- Límite de 1MB alcanzado, truncando resultados ---\n";
                break;
            }
            hhOut << line_out;
            bytes_written += line_out.size();
        }
        hhOut << "\n";
        cerr << "  Heavy hitters encontrados: " << HHs.size()
             << " | escritos: " << (bytes_written / 1024.0) << " KB\n";
    }

    hhOut.close();
    cerr << "Archivo generado:\n  HH exactos: " << outHHFile << "\n";
}

//Void de frecuencias
void extraerFrecuenciasExactas(
    const string &listaArchivos,
    int k,
    const string &outFreqFile
) {
    ifstream listin(listaArchivos);
    if (!listin) {
        cerr << "Error: no se pudo abrir lista de archivos: " << listaArchivos << "\n";
        return;
    }

    ofstream freqOut(outFreqFile);
    if (!freqOut) {
        cerr << "Error: no se pudo abrir archivo de frecuencias: " << outFreqFile << "\n";
        return;
    }

    string filename;
    while (listin >> filename) {
        cerr << "Procesando archivo: " << filename << "\n";

        unordered_map<u64, uint64_t> exact;
        exact.reserve(1 << 20);
        uint64_t N = 0;

        ifstream fin(filename);
        if (!fin) {
            cerr << "  Warning: no pude abrir " << filename << " (se omite)\n";
            continue;
        }

        string line;
        while (getline(fin, line)) {
            for (char &c : line) c = toupper((unsigned char)c);
            string clean;
            clean.reserve(line.size());
            for (char c : line)
                if (c == 'A' || c == 'C' || c == 'G' || c == 'T')
                    clean.push_back(c);

            if ((int)clean.size() < k) continue;

            for (size_t i = 0; i + k <= clean.size(); ++i) {
                string sub = clean.substr(i, k);
                u64 val = encode_kmer(sub);
                if (val == UINT64_MAX) continue;
                u64 canon = canonical_kmer_bits(val, k);
                exact[canon]++;
                N++;
            }
        }
        fin.close();

        // Calcular memoria aproximada
        size_t mem_bytes = exact.size() * (sizeof(u64) + sizeof(uint64_t) + sizeof(size_t));
        double mem_mb = mem_bytes / (1024.0 * 1024.0);

        freqOut << "### Archivo: " << filename << " ###\n";
        freqOut << "# Total_kmers: " << N << "\n";
        freqOut << "# Memoria_MB: " << mem_mb << "\n";

        for (const auto &p : exact) {
            freqOut << decode_kmer_digits(p.first, k) << " " << p.second << "\n";
        }
        freqOut << "\n";

        cerr << "  Total k-mers válidos: " << N
             << " | Distintos: " << exact.size()
             << " | Memoria aprox: " << mem_mb << " MB\n";
    }

    freqOut.close();
    cerr << "Archivo generado:\n  Frecuencias exactas: " << outFreqFile << "\n";
}

// ===================== MAIN =====================
int main(int argc, char* argv[]) {
    if (argc < 5) {
        cerr << "Uso HH:   " << argv[0] << " hh lista_archivos.txt k phi out_hh.txt\n";
        cerr << "Uso Freq: " << argv[0] << " freq lista_archivos.txt k out_freq.txt\n";
        return 1;
    }

    string modo = argv[1];
    if (modo == "hh") {
        if (argc < 6) {
            cerr << "Faltan argumentos para modo hh\n";
            return 1;
        }
        string listaArchivos = argv[2];
        int k = stoi(argv[3]);
        double phi = stod(argv[4]);
        string outHHFile = argv[5];
        extraerHeavyHittersExactos(listaArchivos, k, phi, outHHFile);
    }
    else if (modo == "freq") {
        string listaArchivos = argv[2];
        int k = stoi(argv[3]);
        string outFreqFile = argv[4];
        extraerFrecuenciasExactas(listaArchivos, k, outFreqFile);
    }
    else {
        cerr << "Modo desconocido: " << modo << "\n";
        return 1;
    }

    return 0;
}
