/********************************************************************************************
* SIDH: an efficient supersingular isogeny cryptography library
*
* Abstract: ephemeral supersingular isogeny Diffie-Hellman key exchange (SIDH)
*********************************************************************************************/ 

#include "random/random.h"


static void init_basis(digit_t *gen, f2elm_t XP, f2elm_t XQ, f2elm_t XR)
{ // Initialization of basis points
    
    fpcopy(gen,                  XP[0]);
    fpcopy(gen +   NWORDS_FIELD, XP[1]);
    fpcopy(gen + 2*NWORDS_FIELD, XQ[0]);
    fpcopy(gen + 3*NWORDS_FIELD, XQ[1]);
    fpcopy(gen + 4*NWORDS_FIELD, XR[0]);
    fpcopy(gen + 5*NWORDS_FIELD, XR[1]);
}

static void init_basis_base_for_Bob(digit_t* gen, f2elm_t XPB, f2elm_t XQB, f2elm_t YQB)
{ // Initialization of basis points

    fpcopy(gen, XPB[0]);
    fpzero(XPB[1]);
    fpcopy(gen, XQB[0]);
    fpneg(XQB[0]);
    fpzero(XQB[1]);
    fpzero(YQB[0]);
    fpcopy(gen + NWORDS_FIELD, YQB[1]);
}


void random_mod_order_A(unsigned char* random_digits)
{  // Generation of Alice's secret key  
   // Outputs random value in [0, 2^eA - 1]

    randombytes(random_digits, SECRETKEY_A_BYTES);
    random_digits[SECRETKEY_A_BYTES-1] &= MASK_ALICE;    // Masking last byte 
}


void random_mod_order_B(unsigned char* random_digits)
{  // Generation of Bob's secret key  
   // Outputs random value in [0, 2^Floor(Log(2, oB)) - 1]

    randombytes(random_digits, SECRETKEY_B_BYTES);
    random_digits[SECRETKEY_B_BYTES-1] &= MASK_BOB;     // Masking last byte 
}


int EphemeralKeyGeneration_A(const unsigned char* PrivateKeyA, unsigned char* PublicKeyA)
{ // Alice's ephemeral public key generation
  // Input:  a private key PrivateKeyA in the range [0, 2^eA - 1]. 
  // Output: the public key PublicKeyA consisting of 3 elements in GF(p^2) which are encoded by removing leading 0 bytes.
    point_proj_t R, R0 = { 0 }, R1 = { 0 }, R2 = { 0 }, phiP = { 0 }, phiQ = { 0 }, phiR = { 0 },pts[MAX_INT_POINTS_ALICE];
    point_full_proj_t tR = { 0 }, tR0 = { 0 }, tR1 = { 0 };
    f2elm_t XPA, XQA, XRA, coeff[3], A24plus = { 0 }, A24minus = { 0 },C24 = { 0 };
    unsigned int i, row, m, index = 0, pts_index[MAX_INT_POINTS_ALICE], npts = 0, ii = 0;
    int bit1 = 0;
    digit_t SecretKeyA[NWORDS_ORDER] = {0};

    // Initialize basis points
    init_basis((digit_t*)A_gen, XPA, XQA, XRA);
    init_basis((digit_t*)B_gen, phiP->X, phiQ->X, phiR->X);
    fpcopy((digit_t*)&Montgomery_one, (phiP->Z)[0]);
    fpcopy((digit_t*)&Montgomery_one, (phiQ->Z)[0]);
    fpcopy((digit_t*)&Montgomery_one, (phiR->Z)[0]);

    // Initialize constants: A24plus = A+2C, C24 = 4C, where A=6, C=1
    fpcopy((digit_t*)&Montgomery_one, A24plus[0]);
    mp2_add(A24plus, A24plus, A24plus);
    mp2_add(A24plus, A24plus, C24);
    fp2copy(C24, A24minus);
    mp2_add(C24, C24, A24plus);

    // Retrieve kernel point
    decode_to_digits(PrivateKeyA, SecretKeyA, SECRETKEY_A_BYTES, NWORDS_ORDER);
    LADDER3PT_for_Alice((digit_t*)&DBL_QA[0], SecretKeyA, R0, R1);

    fp2copy(R0->X, R2->X);
    fp2copy(R0->Z, R2->Z);
    fpcopy((digit_t*)&Montgomery_one, (tR->X)[0]);
    fpzero((tR->X)[1]);
    fpcopy((digit_t*)&YQA4, (tR->Y)[0]);
    fpzero((tR->Y)[1]);
    fpcopy((digit_t*)&Montgomery_one, (tR->Z)[0]);
    fpzero((tR->Z)[1]);
    xADD2((tR->X)[0], R2, R1->X, R1->Z);
    RecoverY3PT(R0, R1, R2, tR);
    plus_for_Alice2(SecretKeyA, tR);
    bit1 = (SecretKeyA[0 >> LOG2RADIX] >> (0 & (RADIX - 1))) & 1;
    bit1 = -bit1;
    fp2copy(XPA, tR0->X);
    fp2copy(YPA, tR0->Y);
    fpcopy((digit_t*)&Montgomery_one, tR0->Z[0]);
    fpzero(tR0->Z[1]);
    fpcopy((digit_t*)&PplusQA[0], tR1->X[0]);
    fpcopy((digit_t*)&PplusQA[NWORDS_FIELD], tR1->X[1]);
    fpcopy((digit_t*)&PplusQA[2 * NWORDS_FIELD], tR1->Y[0]);
    fpcopy((digit_t*)&PplusQA[3 * NWORDS_FIELD], tR1->Y[1]);
    fpcopy(&Montgomery_one, tR1->Z[0]);
    fpzero(tR1->Z[1]);
    swap_points2(tR0, tR1, bit1);
    plus_for_Alice3(tR, tR0, R);
    xTPL(R, R, A24minus, A24plus);

#if (OALICE_BITS % 2 == 1)
    point_proj_t S;

    xDBLe(R, S, A24plus, C24, (int)(OALICE_BITS-1));
    get_2_isog(S, A24plus, C24); 
    eval_2_isog(phiP, S); 
    eval_2_isog(phiQ, S); 
    eval_2_isog(phiR, S);
    eval_2_isog(R, S);
#endif

    // Traverse tree
    index = 0;        
    for (row = 1; row < MAX_Alice; row++) {
        while (index < MAX_Alice-row) {
            fp2copy(R->X, pts[npts]->X);
            fp2copy(R->Z, pts[npts]->Z);
            pts_index[npts++] = index;
            m = strat_Alice[ii++];
            xDBLe(R, R, A24plus, C24, (int)(2*m));
            index += m;
        }
        get_4_isog(R, A24plus, C24, coeff);        

        for (i = 0; i < npts; i++) {
            eval_4_isog(pts[i], coeff);
        }
        eval_4_isog(phiP, coeff);
        eval_4_isog(phiQ, coeff);
        eval_4_isog(phiR, coeff);

        fp2copy(pts[npts-1]->X, R->X); 
        fp2copy(pts[npts-1]->Z, R->Z);
        index = pts_index[npts-1];
        npts -= 1;
    }

    get_4_isog(R, A24plus, C24, coeff); 
    eval_4_isog(phiP, coeff);
    eval_4_isog(phiQ, coeff);
    eval_4_isog(phiR, coeff);

    inv_3_way(phiP->Z, phiQ->Z, phiR->Z);
    fp2mul_mont(phiP->X, phiP->Z, phiP->X);
    fp2mul_mont(phiQ->X, phiQ->Z, phiQ->X);
    fp2mul_mont(phiR->X, phiR->Z, phiR->X);
                
    // Format public key                   
    fp2_encode(phiP->X, PublicKeyA);
    fp2_encode(phiQ->X, PublicKeyA + FP2_ENCODED_BYTES);
    fp2_encode(phiR->X, PublicKeyA + 2*FP2_ENCODED_BYTES);

    return 0;
}

int EphemeralKeyGeneration_B(const unsigned char* PrivateKeyB, unsigned char* PublicKeyB)
{ // Bob's ephemeral public key generation
  // Input:  a private key PrivateKeyB in the range [0, 2^Floor(Log(2,oB)) - 1]. 
  // Output: the public key PublicKeyB consisting of 3 elements in GF(p^2) which are encoded by removing leading 0 bytes.
    point_proj_t R = { 0 }, R0 = { 0 }, R1 = { 0 }, R2 = { 0 }, phiP = { 0 }, phiQ = { 0 }, phiR = { 0 }, pts[MAX_INT_POINTS_BOB];
    point_full_proj_t tR;
    felm_t* pre = 0;
    pre = (felm_t*)pre_for_Bob + OBOB_BITS + 2;
    f2elm_t XPB, XQB, YQB, coeff[3], A24plus = { 0 }, A24minus = { 0 }, A24 = { 0 }, C24 = { 0 }, A = { 0 };
    unsigned int i, row, m, index = 0, pts_index[MAX_INT_POINTS_BOB], npts = 0, ii = 0;
    digit_t SecretKeyB[NWORDS_ORDER] = {0};

    // Initialize basis points
    init_basis_base_for_Bob((digit_t*)pre_for_Bob, XPB, XQB, YQB);
    init_basis((digit_t*)A_gen, phiP->X, phiQ->X, phiR->X);
    fpcopy((digit_t*)&Montgomery_one, (phiP->Z)[0]);
    fpcopy((digit_t*)&Montgomery_one, (phiQ->Z)[0]);
    fpcopy((digit_t*)&Montgomery_one, (phiR->Z)[0]);

    
    // Initialize constants: A24 = A+2C, C24 = 4C, where A=0, C=1
    fpcopy((digit_t*)&Montgomery_one, A24plus[0]);
    mp2_add(A24plus, A24plus, A24plus);
    fp2copy(A24plus, A24);
    mp2_add(A24plus, A24plus, C24);
    // Initialize constants: A24minus = A-2C, A24plus = A+2C, where A=6, C=1
    mp2_add(A24plus, A24plus, A24minus);
    mp2_add(A24plus, A24minus, A);
    mp2_add(A24minus, A24minus, A24plus);

    // Retrieve kernel point
    decode_to_digits(PrivateKeyB, SecretKeyB, SECRETKEY_B_BYTES, NWORDS_ORDER);
    LADDER3PT_for_Bob(XPB[0], SecretKeyB, R0, R1);
    xDBL(R0, R0, A24, C24);
    xDBL(R0, R0, A24, C24);
    xDBL(R1, R1, A24, C24);
    xDBL(R1, R1, A24, C24);
    fp2copy(R0->X, R2->X);
    fp2copy(R0->Z, R2->Z);
    pre = pre + 1;
    fpcopy(*pre, (tR->X)[0]);
    fpzero((tR->X)[1]);
    pre = pre + 1;
    fpcopy(*pre, (tR->Y)[0]);
    fpzero((tR->Y)[1]);
    fpcopy((digit_t*)&Montgomery_one, (tR->Z)[0]);
    fpzero((tR->Z)[1]);
    xADD2((tR->X)[0], R2, R1->X, R1->Z);
    RecoverY3PT(R0, R1, R2, tR);
    plus_for_Bob(tR, XQB, YQB, R);
    iso_for_Bob(R, R);
    
    // Traverse tree
    index = 0;  
    for (row = 1; row < MAX_Bob; row++) {
        while (index < MAX_Bob-row) {
            fp2copy(R->X, pts[npts]->X);
            fp2copy(R->Z, pts[npts]->Z);
            pts_index[npts++] = index;
            m = strat_Bob[ii++];
            xTPLe(R, R, A24minus, A24plus, (int)m);
            index += m;
        } 
        get_3_isog(R, A24minus, A24plus, coeff);

        for (i = 0; i < npts; i++) {
            eval_3_isog(pts[i], coeff);
        }     
        eval_3_isog(phiP, coeff);
        eval_3_isog(phiQ, coeff);
        eval_3_isog(phiR, coeff);

        fp2copy(pts[npts-1]->X, R->X); 
        fp2copy(pts[npts-1]->Z, R->Z);
        index = pts_index[npts-1];
        npts -= 1;
    }
    
    get_3_isog(R, A24minus, A24plus, coeff);
    eval_3_isog(phiP, coeff);
    eval_3_isog(phiQ, coeff);
    eval_3_isog(phiR, coeff);

    inv_3_way(phiP->Z, phiQ->Z, phiR->Z);
    fp2mul_mont(phiP->X, phiP->Z, phiP->X);
    fp2mul_mont(phiQ->X, phiQ->Z, phiQ->X);
    fp2mul_mont(phiR->X, phiR->Z, phiR->X);
    // Format public key
    fp2_encode(phiP->X, PublicKeyB);
    fp2_encode(phiQ->X, PublicKeyB + FP2_ENCODED_BYTES);
    fp2_encode(phiR->X, PublicKeyB + 2*FP2_ENCODED_BYTES);

    return 0;
}

int EphemeralSecretAgreement_A(const unsigned char* PrivateKeyA, const unsigned char* PublicKeyB, unsigned char* SharedSecretA)
{ // Alice's ephemeral shared secret computation
  // It produces a shared secret key SharedSecretA using her secret key PrivateKeyA and Bob's public key PublicKeyB
  // Inputs: Alice's PrivateKeyA is an integer in the range [0, oA-1]. 
  //         Bob's PublicKeyB consists of 3 elements in GF(p^2) encoded by removing leading 0 bytes.
  // Output: a shared secret SharedSecretA that consists of one element in GF(p^2) encoded by removing leading 0 bytes.  
    point_proj_t R, pts[MAX_INT_POINTS_ALICE];
    f2elm_t coeff[3], PKB[3], jinv;
    f2elm_t A24plus = {0}, C24 = {0}, A = {0};
    unsigned int i, row, m, index = 0, pts_index[MAX_INT_POINTS_ALICE], npts = 0, ii = 0;
    digit_t SecretKeyA[NWORDS_ORDER] = {0};
      
    // Initialize images of Bob's basis
    fp2_decode(PublicKeyB, PKB[0]);
    fp2_decode(PublicKeyB + FP2_ENCODED_BYTES, PKB[1]);
    fp2_decode(PublicKeyB + 2*FP2_ENCODED_BYTES, PKB[2]);

    // Initialize constants: A24plus = A+2C, C24 = 4C, where C=1
    get_A(PKB[0], PKB[1], PKB[2], A);
    mp_add((digit_t*)&Montgomery_one, (digit_t*)&Montgomery_one, C24[0], NWORDS_FIELD);
    mp2_add(A, C24, A24plus);
    mp_add(C24[0], C24[0], C24[0], NWORDS_FIELD);

    // Retrieve kernel point
    decode_to_digits(PrivateKeyA, SecretKeyA, SECRETKEY_A_BYTES, NWORDS_ORDER);
    LADDER3PT(PKB[0], PKB[1], PKB[2], SecretKeyA, ALICE, R, A);    

#if (OALICE_BITS % 2 == 1)
    point_proj_t S;

    xDBLe(R, S, A24plus, C24, (int)(OALICE_BITS-1));
    get_2_isog(S, A24plus, C24);
    eval_2_isog(R, S);
#endif

    // Traverse tree
    index = 0;        
    for (row = 1; row < MAX_Alice; row++) {
        while (index < MAX_Alice-row) {
            fp2copy(R->X, pts[npts]->X);
            fp2copy(R->Z, pts[npts]->Z);
            pts_index[npts++] = index;
            m = strat_Alice[ii++];
            xDBLe(R, R, A24plus, C24, (int)(2*m));
            index += m;
        }
        get_4_isog(R, A24plus, C24, coeff);        

        for (i = 0; i < npts; i++) {
            eval_4_isog(pts[i], coeff);
        }

        fp2copy(pts[npts-1]->X, R->X); 
        fp2copy(pts[npts-1]->Z, R->Z);
        index = pts_index[npts-1];
        npts -= 1;
    }

    get_4_isog(R, A24plus, C24, coeff); 
    mp2_add(A24plus, A24plus, A24plus);                                                
    fp2sub(A24plus, C24, A24plus); 
    fp2add(A24plus, A24plus, A24plus);                    
    j_inv(A24plus, C24, jinv);
    fp2_encode(jinv, SharedSecretA);    // Format shared secret

    return 0;
}

int EphemeralSecretAgreement_B(const unsigned char* PrivateKeyB, const unsigned char* PublicKeyA, unsigned char* SharedSecretB)
{ // Bob's ephemeral shared secret computation
  // It produces a shared secret key SharedSecretB using his secret key PrivateKeyB and Alice's public key PublicKeyA
  // Inputs: Bob's PrivateKeyB is an integer in the range [0, 2^Floor(Log(2,oB)) - 1]. 
  //         Alice's PublicKeyA consists of 3 elements in GF(p^2) encoded by removing leading 0 bytes.
  // Output: a shared secret SharedSecretB that consists of one element in GF(p^2) encoded by removing leading 0 bytes.  
    point_proj_t R, pts[MAX_INT_POINTS_BOB];
    f2elm_t coeff[3], PKB[3], jinv;
    f2elm_t A24plus = {0}, A24minus = {0}, A = {0};
    unsigned int i, row, m, index = 0, pts_index[MAX_INT_POINTS_BOB], npts = 0, ii = 0;
    digit_t SecretKeyB[NWORDS_ORDER] = {0};
      
    // Initialize images of Alice's basis
    fp2_decode(PublicKeyA, PKB[0]);
    fp2_decode(PublicKeyA + FP2_ENCODED_BYTES, PKB[1]);
    fp2_decode(PublicKeyA + 2*FP2_ENCODED_BYTES, PKB[2]);

    // Initialize constants: A24plus = A+2C, A24minus = A-2C, where C=1
    get_A(PKB[0], PKB[1], PKB[2], A);
    mp_add((digit_t*)&Montgomery_one, (digit_t*)&Montgomery_one, A24minus[0], NWORDS_FIELD);
    mp2_add(A, A24minus, A24plus);
    mp2_sub_p2(A, A24minus, A24minus);


    // Retrieve kernel point
    decode_to_digits(PrivateKeyB, SecretKeyB, SECRETKEY_B_BYTES, NWORDS_ORDER);
    mp_add(SecretKeyB, SecretKeyB, SecretKeyB, NWORDS_ORDER);
    mp_add(SecretKeyB, SecretKeyB, SecretKeyB, NWORDS_ORDER);
    LADDER3PT2(PKB[1], PKB[0], PKB[2], SecretKeyB, BOB, R, A);
    
    // Traverse tree
    index = 0;  
    for (row = 1; row < MAX_Bob; row++) {
        while (index < MAX_Bob-row) {
            fp2copy(R->X, pts[npts]->X);
            fp2copy(R->Z, pts[npts]->Z);
            pts_index[npts++] = index;
            m = strat_Bob[ii++];
            xTPLe(R, R, A24minus, A24plus, (int)m);
            index += m;
        }
        get_3_isog(R, A24minus, A24plus, coeff);

        for (i = 0; i < npts; i++) {
            eval_3_isog(pts[i], coeff);
        } 

        fp2copy(pts[npts-1]->X, R->X); 
        fp2copy(pts[npts-1]->Z, R->Z);
        index = pts_index[npts-1];
        npts -= 1;
    }
     
    get_3_isog(R, A24minus, A24plus, coeff);    
    fp2add(A24plus, A24minus, A);                 
    fp2add(A, A, A);
    fp2sub(A24plus, A24minus, A24plus);                   
    j_inv(A, A24plus, jinv);
    fp2_encode(jinv, SharedSecretB);    // Format shared secret

    return 0;
}
