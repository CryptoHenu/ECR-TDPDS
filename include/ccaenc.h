#ifndef CCAENC_H
#define CCAENC_H

#include "ccastruct.h"
#include "cpastruct.h"
#include "pbc.h"

void ccaEnc(pairing_t pairing, pkg_params pkg_params, ts_params ts_params, element_t user_Alice_Pub, UserPrivateKey User_Alice_Priv, element_t Time_Pub, element_t vk, element_t PT, ccaCiphertext &PCT);

void ccaReEnc(pairing_t pairing, ccaCiphertext PCT, element_t rk, pkg_params pkg_params, element_t vk, ccaReCiphertext &RCT);


#endif