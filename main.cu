#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <cuda.h>
#include "sha256.cuh"
#include <dirent.h>
#include <ctype.h>
#include <math.h>
#include <gmp.h>

#define ROUNDS 50000

#define THREADS 1000
#define BLOCKS 1000

struct BLOCK {
    BYTE  version[04];
    BYTE prevHash[32];
    BYTE mrklRoot[32];
    BYTE timeStmp[04];
    BYTE bitsDiff[04];

    BYTE diffTrgt[32];

    long int x;
};
void printHashD(BYTE *hash) {

    printf("Hash : ");
    for (int i=0; i<32; i++) {
        printf("%02X", hash[i]);
    }
    printf("\n");
}

void d_int_to_little_endian_hex(long int num, BYTE *hex) {
    for (int i = 0; i < sizeof(int); i++) {
        hex[i] = (num >> (i * 8)) & 0xFF;
    }
}

void d_reverse(BYTE* arr, int size) {
    int start = 0;
    int end = size - 1;

    while (start < end) {
        BYTE temp = arr[start];
        arr[start] = arr[end];
        arr[end] = temp;
        start++; end--;
    }
}

__device__ void printHash(BYTE *hash) {

    printf("Hash : ");
    for (int i=0; i<32; i++) {
        printf("%02X", hash[i]);
    }
    printf("\n");
}

__device__ void reverse(BYTE* arr, int size) {
    int start = 0;
    int end = size - 1;

    while (start < end) {
        BYTE temp = arr[start];
        arr[start] = arr[end];
        arr[end] = temp;
        start++; end--;
    }
}

__device__ void int_to_little_endian_hex(long int num, BYTE *hex) {
    for (int i = 0; i < sizeof(int); i++) {
        hex[i] = (num >> (i * 8)) & 0xFF;
    }
}

__device__ bool is_valid(BYTE* blockHash, BYTE* targetHash) {

    for(int i = 0; i < 32; i++) {
        if(blockHash[i] > targetHash[i]) {
            return false;
        } else if(blockHash[i] < targetHash[i]) {
            return true;
        }
    }
    return false;
}

__device__ bool foundHash;

__global__ void startMining(BLOCK* block) {
    if (foundHash)
        return;

    BYTE digest[32];
    long int random;
    for (long int round=0; round<ROUNDS; round++) {

        random = blockIdx.x * (THREADS * ROUNDS) + threadIdx.x * ROUNDS + (round + 1);
        
        BYTE nonce[sizeof(int)];
        int_to_little_endian_hex(random, nonce);

        SHA256_CTX ctx;
        sha256_init(&ctx);
        sha256_update(&ctx, block->version, sizeof(block->version));
        sha256_update(&ctx, block->prevHash, sizeof(block->prevHash));
        sha256_update(&ctx, block->mrklRoot, sizeof(block->mrklRoot));
        sha256_update(&ctx, block->timeStmp, sizeof(block->timeStmp));
        sha256_update(&ctx, block->bitsDiff, sizeof(block->bitsDiff));
        sha256_update(&ctx, nonce, sizeof(nonce));
        sha256_final(&ctx, digest);

        sha256_init(&ctx);
        sha256_update(&ctx, digest, 32);
        sha256_final(&ctx, digest);
        reverse(digest, 32);

        if (is_valid(digest, block->diffTrgt)) {
            printHash(digest);
            printHash(block->diffTrgt);
            block->x = random;
            foundHash = true;
        }
    }
}

void hexstring_to_bytearray(const char *hexstring, BYTE *bytearray) {
    size_t i;
    size_t str_len = strlen(hexstring);

    for (i = 0; i < (str_len / 2); i++) {
        sscanf(hexstring + 2*i, "%2hhx", &bytearray[i]);
    }
}

void get_target_str(long int bits, BYTE* target) {

    long int exp = 8*((bits >> 24)-3);
    long int mant = bits & 0xffffff;

    mpz_t result, multiplied;
    mpz_inits(result, multiplied, NULL);
    mpz_ui_pow_ui(result, 2, exp);

    mpz_set_ui(multiplied, mant);
    mpz_mul(multiplied, multiplied, result);
    size_t count;
    BYTE *hex_bytes = (BYTE *)mpz_export(NULL, &count, 1, sizeof(BYTE), 1, 0, multiplied);

    for (size_t i=0; i<(32-count); ++i) {
        target[i] = 0x00;
    }
    for (size_t i=(32-count); i<32; ++i) {
        target[i] = hex_bytes[i-(32-count)];
    }
    mpz_clears(result, multiplied, NULL);
    free(hex_bytes);
}

int main() {

    long int version = 616259584;
    char prevBlock[] = "0000000000000000000286859eb09e9d6ab0bd39f9d6b1d2bce0a842851aa5db";
    char mrklRoot[] = "50b430dd0d84c3ffb2c5fd6d6b84f88dfc0aa1a92eafd469c549f71ad0cde202";
    long int time = 1713102559;
    long int bits = 486089497;
    
    BLOCK *d_currentBlock;
    BLOCK h_currentBlock;

    hexstring_to_bytearray(prevBlock, h_currentBlock.prevHash);
    hexstring_to_bytearray(mrklRoot, h_currentBlock.mrklRoot);
    d_reverse(h_currentBlock.prevHash, 32); d_reverse(h_currentBlock.mrklRoot, 32);

    d_int_to_little_endian_hex(version, h_currentBlock.version);
    d_int_to_little_endian_hex(time, h_currentBlock.timeStmp);
    d_int_to_little_endian_hex(bits, h_currentBlock.bitsDiff);
    get_target_str(bits, h_currentBlock.diffTrgt);

    cudaMalloc(&d_currentBlock, sizeof(BLOCK));
    cudaMemcpy(d_currentBlock, &h_currentBlock, sizeof(BLOCK), cudaMemcpyHostToDevice);

    startMining<<<BLOCKS, THREADS>>>(d_currentBlock);
    cudaDeviceSynchronize();
    
    cudaMemcpy(&h_currentBlock, d_currentBlock, sizeof(BLOCK), cudaMemcpyDeviceToHost);
    cudaFree(d_currentBlock);
    printf("x = %ld\n", h_currentBlock.x);

    return 0;
}
