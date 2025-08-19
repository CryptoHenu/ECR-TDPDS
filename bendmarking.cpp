/*
 * @Coding: UTF-8
 * @Author: Xiuling Li, Ziyi Dong
 * @Created: 05-14-2025
 * @Last Modified: 05-24-2025
 * @Copyright: © 2025 Ziyi Dong. All rights reserved.
 * @License: GPL v3.0
 * @Contact: ziyidong.cs@gmail.com
 */

#include <stdio.h>
#include <string.h>
#include <time.h>
#include <iostream>
#include "sha.h"
#include "pbc.h"
#include "wots.h"
#include "bendmarking.h"
#include "cpamaptozr.h"

#define RENUM 10000
#define SHA256_DIGEST_LENGTH 32

using namespace std;


void binary_string_to_G1(element_t g1, const char* binary_str, pairing_t pairing) {
    size_t len = strlen(binary_str);
    unsigned char bytes[len/8 + 1];
    memset(bytes, 0, sizeof(bytes));
    
    for(size_t i = 0; i < len; i++) {
        if(binary_str[i] == '1') {
            bytes[i/8] |= (1 << (7 - (i % 8)));
        }
    }
    
    element_init_G1(g1, pairing);
    element_from_bytes(g1, bytes);
}


void G2_to_Zr_via_hash(element_t z, element_t g2, pairing_t pairing) {
    unsigned char buffer[1024]; // 足够大的缓冲区
    int len = element_length_in_bytes(g2);
    
    element_init_Zr(z, pairing);
    element_to_bytes(buffer, g2);
    
    // 使用哈希函数（如SHA-256）处理buffer
    // 这里简化处理，实际应使用密码学哈希
    element_from_hash(z, buffer, len);
}


int bendmarking()
{

    FILE *file;
    file = fopen("bendmarking_output.txt", "w");
    if (!file)
    {
        perror("[Fail] Bendmarking_output.txt open fail.\n");
        exit(1);
    }
    fprintf(file, "=== Bendmarking Test Start === \n");
    fclose(file);

    int i;
    pairing_t pairing;
    element_t P;
    element_t Q, H, R, a, b, c;
    element_t BP;
    element_t a1, b1, c1;
    double relative_time;

    FILE *fp = fopen("../param/d201.param", "r");
    if (!fp)
    {
        printf("[Fail] Param file open fail.\n");
        return 1;
    }
    char param[10240];
    size_t count = fread(param, 1, sizeof(param), fp);
    fclose(fp);
    pairing_init_set_str(pairing, param);

    if (!pairing_is_symmetric(pairing))
    {
        printf("[Asymmetric] Pairing is an asymmetric pairing.\n");
    }
    else
    {
        printf("[Symmetric] Pairing is an symmetric pairing.\n");
    }

    element_init_G1(P, pairing);
    element_random(P);
    element_init_G1(Q, pairing);
    element_random(Q);
    element_init_G2(H, pairing);
    element_random(H);
    element_init_G1(R, pairing);
    element_init_Zr(a, pairing);
    element_init_Zr(b, pairing);
    element_init_Zr(c, pairing);
    element_random(a);
    element_random(b);
    element_init_GT(a1, pairing);
    element_init_GT(b1, pairing);
    element_init_GT(c1, pairing);
    element_random(a1);
    element_random(b1);
    element_init_GT(BP, pairing);

    clock_t start_time, end_time;

    // time_point_mul_G1
    start_time = clock();
    for (i = 1; i < RENUM; i++)
    {
        element_mul_zn(R, P, a);
    }
    end_time = clock();
    double time_point_mul_G1 = (double)(end_time - start_time) / CLOCKS_PER_SEC * 1000;
    file = fopen("bendmarking_output.txt", "a");
    if (!file)
    {
        perror("[Fail] Bendmarking_output.txt open fail.\n");
        exit(1);
    }
    fprintf(file, "time_point_mul_G1: %.6f ms, ", time_point_mul_G1);
    relative_time = time_point_mul_G1 / time_point_mul_G1;
    fprintf(file, "relative_time: %.6f \n", relative_time);
    fclose(file);

    // time_point_add_G1
    start_time = clock();
    for (i = 1; i < RENUM; i++)
    {
        element_add(R, P, Q);
    }
    end_time = clock();
    double time_point_add_G1 = (double)(end_time - start_time) / CLOCKS_PER_SEC * 1000;
    file = fopen("bendmarking_output.txt", "a");
    fprintf(file, "time_point_add_G1: %.6f ms, ", time_point_add_G1);
    relative_time = time_point_add_G1 / time_point_mul_G1;
    fprintf(file, "relative_time: %.6f \n", relative_time);
    fclose(file);

    // time_add_Zr
    start_time = clock();
    for (i = 1; i < RENUM; i++)
    {
        element_add(c, a, b);
    }
    end_time = clock();
    double time_add_Zr = (double)(end_time - start_time) / CLOCKS_PER_SEC * 1000;
    file = fopen("bendmarking_output.txt", "a");
    fprintf(file, "time_add_Zr: %.6f ms, ", time_add_Zr);
    relative_time = time_add_Zr / time_point_mul_G1;
    fprintf(file, "relative_time: %.6f \n", relative_time);
    fclose(file);

    // time_del_Zr
    start_time = clock();
    for (i = 1; i < RENUM; i++)
    {
        element_sub(c, a, b);
    }
    end_time = clock();
    double time_del_Zr = (double)(end_time - start_time) / CLOCKS_PER_SEC * 1000;
    file = fopen("bendmarking_output.txt", "a");

    fprintf(file, "time_del_Zr: %.6f ms, ", time_del_Zr);
    relative_time = time_del_Zr / time_point_mul_G1;
    fprintf(file, "relative_time: %.6f \n", relative_time);
    fclose(file);

    // time_mul_Zr
    start_time = clock();
    for (i = 1; i < RENUM; i++)
    {
        element_mul(c, a, b);
    }
    end_time = clock();
    double time_mul_Zr = (double)(end_time - start_time) / CLOCKS_PER_SEC * 1000;
    file = fopen("bendmarking_output.txt", "a");

    fprintf(file, "time_mul_Zr: %.6f ms, ", time_mul_Zr);
    relative_time = time_mul_Zr / time_point_mul_G1;
    fprintf(file, "relative_time:%.6f \n", relative_time);
    fclose(file);

    // time_div_Zr
    start_time = clock();
    for (i = 1; i < RENUM; i++)
    {
        element_div(c, a, b);
    }
    end_time = clock();
    double time_div_Zr = (double)(end_time - start_time) / CLOCKS_PER_SEC * 1000;
    file = fopen("bendmarking_output.txt", "a");

    fprintf(file, "time_div_Zr: %.6f ms, ", time_div_Zr);
    relative_time = time_div_Zr / time_point_mul_G1;
    fprintf(file, "relative_time: %.6f \n", relative_time);
    fclose(file);

    // time_inv_Zr
    start_time = clock();
    for (i = 1; i < RENUM; i++)
    {
        element_invert(c, a);
    }
    end_time = clock();
    double time_inv_Zr = (double)(end_time - start_time) / CLOCKS_PER_SEC * 1000;
    file = fopen("bendmarking_output.txt", "a");

    fprintf(file, "time_inv_Zr: %.6f ms, ", time_inv_Zr);
    relative_time = time_inv_Zr / time_point_mul_G1;
    fprintf(file, "relative_time: %.6f \n", relative_time);
    fclose(file);

    // time_mul_GT
    start_time = clock();
    for (i = 1; i < RENUM; i++)
    {
        element_mul(c1, a1, b1);
    }
    end_time = clock();
    double time_mul_GT = (double)(end_time - start_time) / CLOCKS_PER_SEC * 1000;
    file = fopen("bendmarking_output.txt", "a");

    fprintf(file, "time_mul_GT: %.6f ms, ", time_mul_GT);
    relative_time = time_mul_GT / time_point_mul_G1;
    fprintf(file, "relative_time: %.6f \n", relative_time);
    fclose(file);

    // time_div_GT
    start_time = clock();
    for (i = 1; i < RENUM; i++)
    {
        element_div(c1, a1, b1);
    }
    end_time = clock();
    double time_div_GT = (double)(end_time - start_time) / CLOCKS_PER_SEC * 1000;
    file = fopen("bendmarking_output.txt", "a");
    fprintf(file, "time_div_GT: %.6f ms, ", time_div_GT);
    relative_time = time_div_GT / time_point_mul_G1;
    fprintf(file, "relative_time: %.6f \n", relative_time);
    fclose(file);

    // time_pow_GT
    start_time = clock();
    for (i = 1; i < RENUM; i++)
    {
        element_pow_zn(c1, a1, b1);
    }
    end_time = clock();
    double time_pow_GT = (double)(end_time - start_time) / CLOCKS_PER_SEC * 1000;
    file = fopen("bendmarking_output.txt", "a");

    fprintf(file, "time_pow_GT: %.6f ms, ", time_pow_GT);
    relative_time = time_pow_GT / time_point_mul_G1;
    fprintf(file, "relative_time: %.6f \n", relative_time);
    fclose(file);

    // time_BP_G1_G1_GT
    start_time = clock();
    for (i = 1; i < RENUM; i++)
    {
        pairing_apply(BP, Q, H, pairing);
    }
    end_time = clock();
    double time_BP_G1_G1_GT = (double)(end_time - start_time) / CLOCKS_PER_SEC * 1000;
    file = fopen("bendmarking_output.txt", "a");

    fprintf(file, "time_BP_G1_G1_GT: %.6f ms, ", time_BP_G1_G1_GT);
    relative_time = time_BP_G1_G1_GT / time_point_mul_G1;
    fprintf(file, "relative_time: %.6f \n", relative_time);
    fclose(file);


    element_t g1, g1_inv;
    element_init_G1(g1, pairing);
    element_init_G1(g1_inv, pairing);

    start_time = clock();
    element_random(g1);

    for (i = 1; i < RENUM; i++)
    {
        // 计算逆元
    element_invert(g1_inv, g1);
    }
    end_time = clock();
    double time_neg_G1 = (double)(end_time - start_time) / CLOCKS_PER_SEC * 1000;
    file = fopen("bendmarking_output.txt", "a");

    fprintf(file, "time_neg_G1: %.6f ms, ", time_neg_G1);
    relative_time = time_neg_G1 / time_point_mul_G1;
    fprintf(file, "relative_time: %.6f \n", relative_time);
    fclose(file);

    
    // time_hash1
    unsigned char digest[SHA256_DIGEST_LENGTH];
    char Alice[] = "sender.alice@gmail.com";
    element_t user_Alice_Pub;
    element_init_Zr(user_Alice_Pub, pairing);
    start_time = clock();
    for (i = 1; i < RENUM; i++)
    {
        id_to_zr(pairing, Alice, user_Alice_Pub);
    }
    end_time = clock();
    double time_hash1 = (double)(end_time - start_time) / CLOCKS_PER_SEC * 1000;
    file = fopen("bendmarking_output.txt", "a");

    fprintf(file, "time_hash1: %.6f ms, ", time_hash1);
    relative_time = time_hash1 / time_point_mul_G1;
    fprintf(file, "relative_time: %.6f \n", relative_time);
    fclose(file);


    // time_hash2
    const char* binary_str = "0101010101010101"; // 示例字符串
    start_time = clock();
    for (i = 1; i < RENUM; i++)
    {
        size_t len = strlen(binary_str);
        unsigned char bytes[len/8 + 1];
        memset(bytes, 0, sizeof(bytes));
        
        for(size_t i = 0; i < len; i++) {
            if(binary_str[i] == '1') {
                bytes[i/8] |= (1 << (7 - (i % 8)));
            }
        }

        element_t g1;
        // 2. 将字节数组转换为G1元素
        element_init_G1(g1, pairing);
        element_from_bytes(g1, bytes);
        element_clear(g1);
    }
    end_time = clock();
    double time_hash2 = (double)(end_time - start_time) / CLOCKS_PER_SEC * 1000;
    file = fopen("bendmarking_output.txt", "a");

    fprintf(file, "time_hash2: %.6f ms, ", time_hash2);
    relative_time = time_hash2 / time_point_mul_G1;
    fprintf(file, "relative_time: %.6f \n", relative_time);
    fclose(file);



    // time_hash3
    start_time = clock();
    for (i = 1; i < RENUM; i++)
    {
        element_t g1, z;
        element_init_G2(g1, pairing);
        element_init_Zr(z, pairing);

        G2_to_Zr_via_hash(z, g1, pairing);

        element_clear(g1);
        element_clear(z);
    }
    end_time = clock();
    double time_hash3 = (double)(end_time - start_time) / CLOCKS_PER_SEC * 1000;
    file = fopen("bendmarking_output.txt", "a");

    fprintf(file, "time_hash3: %.6f ms, ", time_hash3);
    relative_time = time_hash3 / time_point_mul_G1;
    fprintf(file, "relative_time: %.6f \n", relative_time);
    fclose(file);





    // time_sign_key_gen
    uint8_t sk_seed[WOTS_N] = {1};
    uint8_t message[WOTS_N] = {0x12};

    uint8_t pk1[WOTS_LEN][WOTS_N];  
    uint8_t pk2[WOTS_LEN][WOTS_N];
    uint8_t sig[WOTS_LEN][WOTS_N];
    start_time = clock();
    for (i = 1; i < RENUM; i++)
    {
        wots_keygen(pk1, sk_seed);
    //print_hex("Public key (wots_keygen)", pk1, WOTS_LEN * WOTS_N);
    }
    end_time = clock();
    double time_sign_key_gen = (double)(end_time - start_time) / CLOCKS_PER_SEC * 1000;
    file = fopen("bendmarking_output.txt", "a");

    fprintf(file, "time_sign_key_gen: %.6f ms, ", time_sign_key_gen);
    relative_time = time_sign_key_gen / time_point_mul_G1;
    fprintf(file, "relative_time: %.6f \n", relative_time);
    fclose(file);


    // time_sign_gen
    start_time = clock();
    for (i = 1; i < RENUM; i++)
    {
        element_t elements[6];
        
        element_init_G1(elements[0], pairing);   // 第1个元素：G1
        element_init_GT(elements[1], pairing);   // 第2个元素：GT
        element_init_G1(elements[2], pairing);   // 第3个元素：G1
        element_init_GT(elements[3], pairing);   // 第4个元素：GT
        element_init_GT(elements[4], pairing);   // 第5个元素：GT
        element_init_G1(elements[5], pairing);   // 第6个元素：G1

        element_random(elements[0]);
        element_random(elements[1]);
        element_random(elements[2]);
        element_random(elements[3]);
        element_random(elements[4]);
        element_random(elements[5]);

        size_t total_len = 0;
        for (int i = 0; i < 6; i++) {
            total_len += element_length_in_bytes(elements[i]);
        }

        // 
        unsigned char *buffer = (unsigned char *)malloc(total_len);
        if (!buffer) {
            perror("[FAIL] Memory allocation failed.");
            exit(1);
        }

        // 
        size_t offset = 0;
        for (int i = 0; i < 6; i++) {
            int len = element_to_bytes(buffer + offset, elements[i]);
            if (len != element_length_in_bytes(elements[i])) {
                fprintf(stderr, "Serialization error: Element %d\n", i);
                free(buffer);
                exit(1);
            }
            offset += len;
        }
        SHA256(buffer, total_len, message); // hash to 256bit
        wots_sign(sig, message, sk_seed);
        free(buffer);
        element_clear(elements[0]);
        element_clear(elements[1]);
        element_clear(elements[2]);
        element_clear(elements[3]);
        element_clear(elements[4]);
        element_clear(elements[5]);
    }
    end_time = clock();
    double time_sign_gen = (double)(end_time - start_time) / CLOCKS_PER_SEC * 1000;
    file = fopen("bendmarking_output.txt", "a");

    fprintf(file, "time_sign_gen: %.6f ms, ", time_sign_gen);
    relative_time = time_sign_gen / time_point_mul_G1;
    fprintf(file, "relative_time: %.6f \n", relative_time);
    fclose(file);


    // time_sign_verify
    start_time = clock();
    for (i = 1; i < RENUM; i++)
    {
        wots_pk_from_sig(pk2, sig, message);
        //print_hex("Recovered public key (wots_pk_from_sig)", pk2, WOTS_LEN * WOTS_N);

        int receiversuccess = 1;
        for (int i = 0; i < WOTS_LEN; i++) {
            if (memcmp(pk1[i], pk2[i], WOTS_N) != 0) {
                receiversuccess = 0;
                break;
            }
        }
        //printf("WOTS+ verification %s\n", receiversuccess ? "passed" : "failed");
    }
    end_time = clock();
    double time_sign_verify = (double)(end_time - start_time) / CLOCKS_PER_SEC * 1000;
    file = fopen("bendmarking_output.txt", "a");

    fprintf(file, "ttime_sign_verify: %.6f ms, ", time_sign_verify);
    relative_time = time_sign_verify / time_point_mul_G1;
    fprintf(file, "relative_time: %.6f \n", relative_time);
    fclose(file);



    file = fopen("bendmarking_output.txt", "a");
    if (!file)
    {
        perror("[Fail] Bendmarking_output.txt open fail.\n");
        exit(1);
    }
    fprintf(file, "=== Bendmarking Test End === \n");
    fclose(file);


















    // clear memory
    element_clear(P);
    element_clear(Q);
    element_clear(R);
    element_clear(a);
    element_clear(b);
    element_clear(c);
    element_clear(a1);
    element_clear(b1);
    element_clear(c1);
    element_clear(BP);

    pairing_clear(pairing);

    return 1;
}