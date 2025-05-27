/*
 * @Coding: UTF-8
 * @Author: Ziyi Dong, Xiuling Li
 * @Created: 05-14-2025
 * @Last Modified: 05-24-2025
 * @Copyright: © 2025 Ziyi Dong. All rights reserved.
 * @License: GPL v3.0
 * @Contact: ziyidong.cs@gmail.com
 */

#include "cpakeygen.h"


// User private key generation
void PrivatekeyGen(pairing_t pairing, element_t pkg_priv, pkg_params pkg_params, element_t user_Alice_Pub, UserPrivateKey &privatekey)
{
    element_t diff, inv;
    element_random(privatekey.r);
    element_init_Zr(diff, pairing);
    element_init_Zr(inv, pairing);

    element_sub(diff, pkg_priv, user_Alice_Pub);
    element_invert(inv, diff);
    element_neg(privatekey.K, pkg_params.g);
    element_pow_zn(privatekey.K, privatekey.K, privatekey.r);
    element_add(privatekey.K, privatekey.K, pkg_params.h);
    element_pow_zn(privatekey.K, privatekey.K, inv);

    element_clear(diff);
    element_clear(inv);
}


// TimeTrapDoor generation
void TimeTrapDoorGen(pairing_t pairing, element_t ts_priv, ts_params ts_params, element_t Time_Pub, TimeTrapDoor &Time_St)
{
    element_t diff, inv;            
    element_random(Time_St.r); 
    element_init_Zr(diff, pairing);
    element_init_Zr(inv, pairing);
    
    element_sub(diff, ts_priv, Time_Pub); // diff �� a - b
    element_invert(inv, diff);
    element_neg(Time_St.K, ts_params.g);
    element_pow_zn(Time_St.K, Time_St.K, Time_St.r);
    element_add(Time_St.K, Time_St.K, ts_params.h);
    element_pow_zn(Time_St.K, Time_St.K, inv);


    element_clear(diff);
    element_clear(inv);
}

// RK, X generation function
void RkGen(pairing_t pairing, pkg_params pkg_params, element_t user_Alice_Pub, UserPrivateKey User_Alice_Priv, Ciphertext PCT, element_t &rk, element_t &X)
{
    element_t Q, temp;
    element_init_G1(Q, pairing);
    element_init_G1(temp, pairing);
    
    element_random(Q);

    element_pow_zn(temp, Q, User_Alice_Priv.r);
    element_add(rk, temp, User_Alice_Priv.K);
    pairing_apply(X, PCT.C3, temp, pairing);

    element_clear(Q);
    element_clear(temp);

}

// Rj generation function;
void RjGen(pairing_t pairing, pkg_params pkg_params, UserPrivateKey User_Alice_Priv, element_t user_Pub, element_t rk, element_t X, element_t k3, Rj &rj)
{
    element_t temp1, temp2;
    element_init_Zr(temp1, pairing);
    element_init_G1(temp2, pairing);

    // u
    element_mul(temp1, k3, user_Pub);
    element_neg(rj.u, pkg_params.g);
    element_pow_zn(rj.u, rj.u, temp1);
    element_pow_zn(temp2, pkg_params.g1, k3);
    element_add(rj.u, rj.u, temp2);

    // v
    element_pow_zn(rj.v, pkg_params.e_g_g, k3);

    // w
    element_invert(rj.w, pkg_params.e_g_h);
    element_pow_zn(rj.w, rj.w, k3);
    element_mul(rj.w, rj.w, X);

    element_clear(temp1);
    element_clear(temp2);

}