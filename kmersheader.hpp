#ifndef KMERSHEADER_HPP
#define KMERSHEADER_HPP

#include <bits/stdc++.h>
using namespace std;
using u64 = uint64_t;

// === Funciones de utilidades para k-mers ===

// Convertir letra -> base2 (0=A, 1=C, 2=G, 3=T)
int base2(char c);

// Complemento de un nucleótido en bits
int complementbit(int b);

// Codificar k-mer a entero (2 bits por base)
u64 encode_kmer(const string &s);

// Reverse complement en base-2
u64 revcomp_bits(u64 val, int k);

// Canonical k-mer en base-2
u64 canonical_kmer_bits(u64 val, int k);

// Decodificar entero -> string en dígitos base-2
string decode_kmer_digits(u64 val, int k);
// Convierte string ACGT a u64 codificado y obtiene el canonical k-mer codificado
inline u64 kmer_to_u64(const std::string &s) {
    return encode_kmer(s);
}
u64 canonical_kmer_bits(u64 val, int k);
// Devuelve la secuencia ACGT canónica
inline std::string canonical_kmer(const std::string &s) {
    u64 val = encode_kmer(s);
    if (val == UINT64_MAX) return ""; // contiene N u carácter inválido
    u64 can_bits = canonical_kmer_bits(val, s.size());

    // Decodificar de 2 bits/base a string ACGT
    std::string out;
    out.reserve(s.size());
    for (int i = (int)s.size() - 1; i >= 0; --i) {
        int b = can_bits & 0b11;
        switch(b) {
            case 0: out += 'A'; break;
            case 1: out += 'C'; break;
            case 2: out += 'G'; break;
            case 3: out += 'T'; break;
        }
        can_bits >>= 2;
    }
    std::reverse(out.begin(), out.end());
    return out;
}

#endif