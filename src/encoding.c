#include "encoding.h"
#include "rans_byte.h"
#include "params.h"

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>

#define SCALE_BITS 16
#define SCALE (1u << SCALE_BITS)

#if HAETAE_MODE == 2
#define M_H 252 // Alphabet size for h
static uint32_t f_h[M_H] = {25334, 15860, 3780, 330, 10, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 9, 325, 3772, 15873};
#define M_HB_Z1 37 // Alphabet size for HighBits(z1)
static uint32_t f_hb_z1[M_HB_Z1] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 4, 235, 3305, 15833, 26660, 15892, 3334, 240, 5, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};

#elif HAETAE_MODE == 3
#define M_H 252 // Alphabet size for h
static uint32_t f_h[M_H] = {17149, 13965, 7177, 2376, 496, 67, 5, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 4, 68, 500, 2375, 7185, 13930};
#define M_HB_Z1 71
static uint32_t f_hb_z1[M_HB_Z1] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 4, 55, 444, 2231, 7047, 14084, 17706, 14067, 7084, 2254, 443, 56, 3, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};

#elif HAETAE_MODE == 5
#define M_H 504 // Alphabet size for h
static uint32_t f_h[M_H] = {8061, 8117, 6881, 5275, 3639, 2241, 1249, 616, 281, 114, 43, 13, 3, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 4, 14, 41, 115, 281, 628, 1247, 2247, 3651, 5288, 6923, 8085};
#define M_HB_Z1 175
static uint32_t f_hb_z1[M_HB_Z1] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 8, 89, 589, 2561, 7297, 13713, 16776, 13740, 7347, 2569, 589, 88, 8, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
#endif // HAETAE_MODE


//TODO: insert hardcoded values:
static RansEncSymbol esyms_h[M_H];
static uint16_t symbol_h[SCALE] = {0};
static RansDecSymbol dsyms_h[M_H];

static RansEncSymbol esyms_hb_z1[M_HB_Z1];
static uint16_t symbol_hb_z1[SCALE] = {0};
static RansDecSymbol dsyms_hb_z1[M_HB_Z1];

// helper functions to generate lookup tables for encoding/decoding, 
// TODO: the output will be hardcoded once all parameters are fixed

void symbol_table(uint16_t *symbol, const uint32_t *freq, size_t alphabet_size) 
{
    int pos = 0;
    for (size_t sym=0; sym < alphabet_size; sym++) {
        for (uint32_t i=0; i < freq[sym]; i++)
            symbol[pos++] = sym;
    }
}

void cum_freq_table(uint32_t *cum_freq, const uint32_t *freq, size_t alphabet_size) {
    cum_freq[0] = 0;
    for (size_t i = 1; i < alphabet_size; i++) {
      cum_freq[i] = cum_freq[i-1] + freq[i-1];
    }
}

void encode_symbols(RansEncSymbol *esyms, const uint32_t *freq, size_t alphabet_size) {
    uint32_t cum_freq[alphabet_size];
    cum_freq_table(cum_freq, freq, alphabet_size);

    for (size_t i=0; i < alphabet_size; i++) 
    {
        RansEncSymbolInit(&esyms[i], cum_freq[i], freq[i], SCALE_BITS);
    }
}

void decode_symbols(RansDecSymbol *dsyms, uint16_t *symbol, const uint32_t *freq, size_t alphabet_size) {
    uint32_t cum_freq[alphabet_size];
    cum_freq_table(cum_freq, freq, alphabet_size);

    symbol_table(symbol, freq, alphabet_size);
    
    for (size_t i=0; i < alphabet_size; i++) 
    {
        RansDecSymbolInit(&dsyms[i], cum_freq[i], freq[i]);
    }
}

void precomputations_rans() {
    encode_symbols(esyms_h, f_h, M_H);
    decode_symbols(dsyms_h, symbol_h, f_h, M_H);

    encode_symbols(esyms_hb_z1, f_hb_z1, M_HB_Z1);
    decode_symbols(dsyms_hb_z1, symbol_hb_z1, f_hb_z1, M_HB_Z1);
}

/*************************************************
 * Name:        encode_h
 *
 * Description: rANS encode polynomial vector h
 *
 * Arguments:   - uint8_t *buf: pointer to output buffer
 *              - const int32_t *h: pointer to polynomial vector h
 **************************************************/
uint8_t* encode_h(uint8_t *buf, const int32_t *h) {
    RansState rans;
    uint8_t* ptr;
    size_t size_h = N*K;
    uint16_t s;

    RansEncInit(&rans);   
    ptr = buf + MAX_ENC_H_BYTES; // end of output buffer

    for (size_t i=size_h; i > 0; i--) 
    {
        s = h[i-1];
        RansEncPutSymbol(&rans, &ptr, &esyms_h[s]);
    }

    RansEncFlush(&rans, &ptr);
    return ptr; 
}

/*************************************************
 * Name:        decode_h
 *
 * Description: rANS decode polynomial vector h
 *
 * Arguments:   - int32_t *h: pointer to polynomial vector h
 *              - uint8_t *buf: pointer to output buffer
 **************************************************/
void decode_h(int32_t *h, uint8_t *buf) {
    RansState rans;
    size_t size_h = N*K;
    uint16_t s;

    RansDecInit(&rans, &buf);

    for (size_t i=0; i < size_h; i++) 
    {
        s = symbol_h[(uint16_t) RansDecGet(&rans, SCALE_BITS)];
        h[i] = s;
        RansDecAdvanceSymbol(&rans, &buf, &dsyms_h[s], SCALE_BITS);
    }
}

/*************************************************
 * Name:        encode_hb_z1
 *
 * Description: rANS encode polynomial vector HighBits(z1)
 *
 * Arguments:   - uint8_t *buf: pointer to output buffer
 *              - const int32_t *hb_z1: pointer to polynomial vector HighBits(z1)
 **************************************************/
uint8_t* encode_hb_z1(uint8_t *buf, const int32_t *hb_z1) {
    RansState rans;
    uint8_t* ptr;
    size_t size_hb_z1 = N*L;
    int16_t s;

    RansEncInit(&rans);   
    ptr = buf + MAX_ENC_HB_Z1_BYTES; // end of output buffer

    for (size_t i=size_hb_z1; i > 0; i--) 
    {
        s = hb_z1[i-1] + M_HB_Z1/2; // from centered to positive representation
        RansEncPutSymbol(&rans, &ptr, &esyms_hb_z1[s]);
    }

    RansEncFlush(&rans, &ptr);
    return ptr; 
}

/*************************************************
 * Name:        decode_hb_z1
 *
 * Description: rANS decode polynomial vector HighBits(z1)
 *
 * Arguments:   - int32_t *hb_z1: pointer to polynomial vector HighBits(z1)
 *              - uint8_t *buf: pointer to output buffer
 **************************************************/
void decode_hb_z1(int32_t *hb_z1, uint8_t *buf) {
    RansState rans;
    size_t size_hb_z1 = N*L;
    int16_t s;

    RansDecInit(&rans, &buf);

    for (size_t i=0; i < size_hb_z1; i++) 
    {
        s = symbol_hb_z1[(uint16_t) RansDecGet(&rans, SCALE_BITS)]; // cast to uint16_t ensures the index for symbol is in range if SCALE_BITS = 16 
        hb_z1[i] = s - M_HB_Z1/2; // from positive to centered representation
        RansDecAdvanceSymbol(&rans, &buf, &dsyms_hb_z1[s], SCALE_BITS);
    }
}