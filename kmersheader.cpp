#include "kmersheader.hpp"

#include <unordered_map>
#include <string>
#include <vector>
#include <algorithm>
#include <iostream>
// PHI para 31 es 3.840459813e-07
// PHI para 21 es 3.840412721e-07
// PHI interesante 8.840459813e-07
int base2(char c) {
    switch(c) {
        case 'A': case 'a': return 0;
        case 'C': case 'c': return 1;
        case 'G': case 'g': return 2;
        case 'T': case 't': return 3;
        default: return -1; // N u otro carácter inválido
    }
}

// Complemento de un nucleótido en bits
int complementbit(int b) {
return 3 - b; // A<->T, C<->G
}

// Codificar k-mer a entero (2 bits por base)
u64 encode_kmer(const string &s) {
    u64 val = 0;
    for (char c : s) {
        int b = base2(c);
        if (b < 0) return UINT64_MAX; // inválido si hay N
        val = (val << 2) | b;
    }
    return val;
}

// Reverse complement en base-2
u64 revcomp_bits(u64 val, int k) {
    u64 rc = 0;
    for (int i = 0; i < k; i++) {
        int b = val & 0b11;
        b = complementbit(b);
        rc = (rc << 2) | b;
        val >>= 2;
    }
    return rc;
}

// Canonical k-mer en base-2
u64 canonical_kmer_bits(u64 val, int k) {
    u64 rc = revcomp_bits(val, k);
    return (val < rc) ? val : rc;
}

// Decodificar entero -> string en dígitos base-2
string decode_kmer_digits(u64 val, int k) {
    string s(k, '0');
    for (int i = k - 1; i >= 0; i--) {
        int b = val & 0b11;   // últimos 2 bits
        s[i] = char('0' + b); // convertir a '0','1','2','3'
        val >>= 2;
    }
    return s;
}
vector<pair<string, size_t>> encontrar_heavy_hitters(
    const string &seq, int k, double phi, size_t max_memory) 
{
    unordered_map<u64, size_t> freq;
    size_t total_kmers = 0;
    int n_kmers = (int)seq.size() - k + 1;
    if (n_kmers <= 0) return {};

    // Conteo exacto de todos los k-mers
    for (int i = 0; i < n_kmers; i++) {
        string sub = seq.substr(i, k);
        u64 enc = encode_kmer(sub);
        if (enc == UINT64_MAX) continue;
        u64 can = canonical_kmer_bits(enc, k);
        freq[can]++;
        total_kmers++;
    }

    if (total_kmers == 0) return {};

    double threshold = phi * total_kmers;

    // Filtrar heavy hitters
    vector<pair<string, size_t>> heavyhitters;
    heavyhitters.reserve(freq.size());

    for (const auto &kv : freq) {
        if (kv.second >= threshold) {
            heavyhitters.push_back({decode_kmer_digits(kv.first, k), kv.second});
        }
    }

    // Ordenar por frecuencia descendente
    sort(heavyhitters.begin(), heavyhitters.end(),
         [](auto &a, auto &b) { return a.second > b.second; });

    // Verificar memoria: si sobrepasa el límite, recorta
    size_t approx_size = heavyhitters.size() * (sizeof(pair<string, size_t>) + k);
    if (approx_size > max_memory) {
        // recorta los más frecuentes hasta que entre en memoria
        size_t max_items = max_memory / (sizeof(pair<string, size_t>) + k);
        if (max_items < heavyhitters.size()) {
            heavyhitters.resize(max_items);
        }
    }

    return heavyhitters;
}

std::unordered_map<u64, size_t> obtener_frecuencias_kmers(
    const std::string &seq,
    int k,
    size_t &total_kmers
) {
    std::unordered_map<u64, size_t> freq;
    total_kmers = 0;

    int n_kmers = (int)seq.size() - k + 1;
    if (n_kmers <= 0) return freq;

    for (int idx = 0; idx < n_kmers; ++idx) {
        std::string sub = seq.substr(idx, k);
        u64 enc = encode_kmer(sub);
        if (enc == UINT64_MAX) continue;  // ignora kmers con N
        u64 can = canonical_kmer_bits(enc, k);
        freq[can]++;
        total_kmers++;
    }

    return freq;
}