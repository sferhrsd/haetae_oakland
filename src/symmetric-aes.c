#include "aes256ctr.h"
#include "symmetric.h"
#include <stdint.h>

void haetae_aes256ctr_init(aes256ctr_ctx *state, const uint8_t key[32],
                           uint16_t nonce) {
    uint8_t expnonce[12] = {0};
    expnonce[0] = nonce;
    expnonce[1] = nonce >> 8;
    aes256ctr_init(state, key, expnonce);
}
