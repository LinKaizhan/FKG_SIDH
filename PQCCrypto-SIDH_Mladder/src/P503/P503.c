/********************************************************************************************
* SIDH: an efficient supersingular isogeny cryptography library
*
* Abstract: supersingular isogeny parameters and generation of functions for P503
*********************************************************************************************/  

#include "P503_api.h" 
#include "P503_internal.h"
#include "../internal.h"


// Encoding of field elements, elements over Z_order, elements over GF(p^2) and elliptic curve points:
// --------------------------------------------------------------------------------------------------
// Elements over GF(p) and Z_order are encoded with the least significant octet (and digit) located at the leftmost position (i.e., little endian format). 
// Elements (a+b*i) over GF(p^2), where a and b are defined over GF(p), are encoded as {a, b}, with a in the least significant position.
// Elliptic curve points P = (x,y) are encoded as {x, y}, with x in the least significant position. 
// Internally, the number of digits used to represent all these elements is obtained by approximating the number of bits to the immediately greater multiple of 32.
// For example, a 503-bit field element is represented with Ceil(503 / 64) = 8 64-bit digits or Ceil(503 / 32) = 16 32-bit digits.

//
// Curve isogeny system "SIDHp503". Base curve: Montgomery curve By^2 = Cx^3 + Ax^2 + Cx defined over GF(p503^2), where A=6, B=1, C=1 and p503 = 2^250*3^159-1
//
         
const uint64_t p503[NWORDS64_FIELD]              = { 0xFFFFFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFF, 0xABFFFFFFFFFFFFFF, 
                                                     0x13085BDA2211E7A0, 0x1B9BF6C87B7E7DAF, 0x6045C6BDDA77A4D0, 0x004066F541811E1E };
const uint64_t p503x2[NWORDS64_FIELD]            = { 0xFFFFFFFFFFFFFFFE, 0xFFFFFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFF, 0x57FFFFFFFFFFFFFF,
                                                     0x2610B7B44423CF41, 0x3737ED90F6FCFB5E, 0xC08B8D7BB4EF49A0, 0x0080CDEA83023C3C }; 
const uint64_t p503x4[NWORDS64_FIELD]            = { 0xFFFFFFFFFFFFFFFC, 0xFFFFFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFF, 0xAFFFFFFFFFFFFFFF, 
                                                     0x4C216F6888479E82, 0x6E6FDB21EDF9F6BC, 0x81171AF769DE9340, 0x01019BD506047879 };
const uint64_t p503p1[NWORDS64_FIELD]            = { 0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0xAC00000000000000,
                                                     0x13085BDA2211E7A0, 0x1B9BF6C87B7E7DAF, 0x6045C6BDDA77A4D0, 0x004066F541811E1E };
const uint64_t p503p1x64[NWORDS64_FIELD/2]       = { 0xC216F6888479E82B, 0xE6FDB21EDF9F6BC4, 0x1171AF769DE93406, 0x1019BD5060478798 };  
const uint64_t p503x16p[2*NWORDS64_FIELD]        = { 0x0000000000000010, 0x0000000000000000, 0x0000000000000000, 0x8000000000000000, 
                                                     0x9EF484BBBDC30BEA, 0x8C8126F090304A1D, 0xF7472844B10B65FC, 0x30F32157CFDC3C33, 
                                                     0x1463AB4329A333F7, 0xDFC933977C47D3A4, 0x338A3767F6F2520B, 0x4F8CB7565CCC13FA, 
                                                     0xDE43B73AACD2189B, 0xBCF845CAC5405FBD, 0x516D02A09E684B7A, 0x0001033A4091BB86 }; 
// Order of Alice's subgroup
const uint64_t Alice_order[NWORDS64_ORDER]       = { 0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0400000000000000 }; 
// Order of Bob's subgroup
const uint64_t Bob_order[NWORDS64_ORDER]         = { 0xC216F6888479E82B, 0xE6FDB21EDF9F6BC4, 0x1171AF769DE93406, 0x1019BD5060478798 };
// Alice's generator values {XPA0 + XPA1*i, XQA0 + XQA1*i, XRA0 + XRA1*i} in GF(p503^2), expressed in Montgomery representation
const uint64_t A_gen[6 * NWORDS64_FIELD]         = { 0x3353B596D45A95A6, 0xDF7E0A94A39B96C0, 0x715DC90A72A3223F, 0xCB73F56E5AD9430F,
                                                     0xE4B5DA591AEE475D, 0x322F1CE730413BD7, 0x4EEA4028D168DAD2, 0xB254087875FFA,    // XPA0
                                                     0xAC3985C5BB18D89D, 0x45F2445C680A1E40, 0xF59454B450FBAB11, 0x95DC27D8152A0DAE,
                                                     0x42A4FD439715E500, 0xB958FBA1CD4CC505, 0xC4E5AB2ABB732FC5, 0x268ED322F62ACA,   // XPA1
                                                     0xDD27E5ADF7F57AB4, 0x7C1379D2B09F0434, 0x6E267408F1C8C89F, 0xC3BB383C07B60035,
                                                     0x9268C9183A95ECD5, 0x9327EC043E0F021F, 0xE63D2D907A9DE5A5, 0x3110B6B4E0CD93,   // XQA0
                                                     0x40B6BC5F5C2675E6, 0x62AD4B61EEDC2C5C, 0xA1CCA6B5091EF540, 0xC6273D4E1D8FC7FE,
                                                     0x266D8B99EE63A78F, 0x39604E6927906566, 0xAB8BA8F2C6A977F8, 0xCD759EE7AB739,    // XQA1
                                                     0x1482EA2C7A8F5FA0, 0xB42C8B9C007E5FE5, 0xCFCFF2625C69E7FD, 0x8334C3F384C268F5,
                                                     0xD71E78E25FA4DB2F, 0x64BECFBE41708879, 0x103FF021EF7BF9, 0x2695BB8221E83B,     // XRA0
                                                     0xA08787E922A1030, 0x8D34581F64BCE547, 0x2FA5BED41306271A, 0xEC24812ABD206DCF,
                                                     0x978FA888C3CC6366, 0x2BFF991CDB7CE058, 0xA0BCCC1A447CF056, 0x2425429A072D82}; // XRA1
const uint64_t DBL_QA[2 * NWORDS64_FIELD]        = { 0x273f6c464cb9ab1a, 0x83722dbae9836b36, 0x7236dd158d1a1bbf, 0xbe84ed2fd6fc9b11,
                                                     0xf4fac85bba91e9b4, 0x783d71c36f23ae76, 0x6fc94cf24bda330a, 0x3929a6320c9596,
                                                     0x3e3209393cb32a2, 0x689964ccab348a84, 0x74471effced8819b, 0x661d7240b28e2790,
                                                     0x71aae7baae2179ca, 0x5da64f579d150d5b, 0x17919259b69ea954, 0x77328aa89bbc7};
const uint64_t YPA[2 * NWORDS64_FIELD]           = { 0xbc88bb85404378e5, 0x61071195bc44bf8f, 0xc92d13994ce9b8b3, 0x9ed615392dcf6ca2, 
                                                     0xc4a95165fb25bfda, 0xeeea8545ebeaec62, 0xac09c1c3e91b41fd, 0xb43ac79a90a0c,
                                                     0xe8b38a79e90eaadb, 0x840b284661ccfc39, 0x6d5091432c311ad, 0xdafbd9cd646033b,
                                                     0x3faf77bc98339af0, 0x75f0c7a7aa5d03a3, 0xd188da98de124c6a, 0x2780b2b7b1c9cc};
const uint64_t PplusQA[4 * NWORDS64_FIELD]       = { 0x1482ea2c7a8f5fa0, 0xb42c8b9c007e5fe5, 0xcfcff2625c69e7fd, 0x8334c3f384c268f5,
                                                     0xd71e78e25fa4db2f, 0x64becfbe41708879, 0x0103ff021ef7bf9, 0x2695bb8221e83b,
                                                     0xf5f787816dd5efcf, 0x72cba7e09b431ab8, 0xd05a412becf9d8e5, 0xbfdb7ed542df9230,
                                                     0x7b78b3515e458439, 0xef9c5daba0019d56, 0xbf88faa395fab479, 0x1c41b2a779f09b,
                                                     0x4f151d6b2697df41, 0xf2286438aadddb71, 0x62378cc5be23004f, 0x822807933e84ad42, 
                                                     0x36db6c363e3d2500, 0x95941f4db77237b2, 0xfd917b6f231a9e7c, 0x1e1e9b5aa4f411,
                                                     0xb3d457620c43d607, 0x96c95df412038dbe, 0xf4e14d69d3ef397c, 0x3ff63724a560957, 
                                                     0xbde046fe40105145, 0x5762a00b1b2c71f, 0x97e2dfed61620d74, 0x1c4742befd3c};
// Bob's generator values {XPB0, XQB0, XRB0 + XRB1*i} in GF(p503^2), expressed in Montgomery representation
const uint64_t B_gen[6 * NWORDS64_FIELD]         = { 0xB810321963CF561F, 0xACA612873FBC647F, 0xE5C29CB78215B634, 0xB277ACABE764F907,
                                                     0x76DBA8FCCDFF4721, 0x1B4E6541441EB543, 0xDAAB92E8B2DD0517, 0x1ECAA65407E4C9,   // XPB0
                                                     0xF7EEE8D8D30365E6, 0x48F0AF97691E0303, 0xA8AC75108BFDA627, 0x7C0F65DCF8450F1,
                                                     0xCD74E9CA0E92BECA, 0x342E232149CA1DFA, 0x8E841EC6D7725DE3, 0x2429A4E9A12CB0,   // XPB1
                                                     0xB810321963CF561F, 0xACA612873FBC647F, 0xE5C29CB78215B634, 0xB277ACABE764F907,
                                                     0x76DBA8FCCDFF4721, 0x1B4E6541441EB543, 0xDAAB92E8B2DD0517, 0x1ECAA65407E4C9,   // XQB0
                                                     0x81117272CFC9A19, 0xB70F506896E1FCFC, 0x57538AEF740259D8, 0xA43F09A2307BAF0E,
                                                     0x45937210137F28D6, 0xE76DD3A731B45FB4, 0xD1C1A7F7030546EC, 0x1C3D5057DFF16D,   // XQB1
                                                     0x6E3DEF7C8A5A47D2, 0x12D9AF90F92FC868, 0xCE33D50FC931894B, 0x2927354E05ED037C,
                                                     0x4864AD1D8B6E4E56, 0x2C6BB7E4CD4284DD, 0x50A30A93843DDC28, 0x38195667C39958,   // XRB0
                                                     0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0 }; // XRB1
/* Basis for Bob on E_0, expressed in Montgomery representation */
const uint64_t P3[2 * NWORDS64_FIELD]            = { 0x4256C520FB388820, 0x744FD7C3BAAF0A13, 0x4B6A2DDDB12CBCB8, 0xE46826E27F427DF8,
                                                     0xFE4A663CD505A61B, 0xD6B3A1BAF025C695, 0x7C3BB62B8FCC00BD, 0x3AFDDE4A35746C,
                                                     0x440192590061240E, 0x60C942451EC3E20D, 0x2195638E3B7632CA, 0xBA84AC322AA59D16,
                                                     0x3751CBF97048E02D, 0x6A583E4C816EAC44, 0x7A984D4F477762C1, 0x27B5AB2E503D63 };


                                                                                                                                       // Montgomery constant Montgomery_R2 = (2^512)^2 mod p503
const uint64_t Montgomery_R2[NWORDS64_FIELD]     = { 0x5289A0CF641D011F, 0x9B88257189FED2B9, 0xA3B365D58DC8F17A, 0x5BC57AB6EFF168EC,
                                                     0x9E51998BD84D4423, 0xBF8999CBAC3B5695, 0x46E9127BCE14CDB6, 0x003F6CFCE8B81771 };                                                   
// Value one in Montgomery representation 
const uint64_t Montgomery_one[NWORDS64_FIELD]    = { 0x00000000000003F9, 0x0000000000000000, 0x0000000000000000, 0xB400000000000000, 
                                                     0x63CB1A6EA6DED2B4, 0x51689D8D667EB37D, 0x8ACD77C71AB24142, 0x0026FBAEC60F5953 };


// Fixed parameters for isogeny tree computation
const unsigned int strat_Alice[MAX_Alice-1] = { 
61, 32, 16, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 8, 4, 2, 1, 1, 2, 1, 1, 
4, 2, 1, 1, 2, 1, 1, 16, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 8, 4, 2, 1, 
1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 29, 16, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 
1, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 13, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 
1, 1, 2, 1, 1, 5, 4, 2, 1, 1, 2, 1, 1, 2, 1, 1, 1 };

const unsigned int strat_Bob[MAX_Bob-1] = { 
71, 38, 21, 13, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 5, 4, 2, 1, 1, 2, 1, 
1, 2, 1, 1, 1, 9, 5, 3, 2, 1, 1, 1, 1, 2, 1, 1, 1, 4, 2, 1, 1, 1, 2, 1, 1, 17, 9, 
5, 3, 2, 1, 1, 1, 1, 2, 1, 1, 1, 4, 2, 1, 1, 1, 2, 1, 1, 8, 4, 2, 1, 1, 1, 2, 1, 
1, 4, 2, 1, 1, 2, 1, 1, 33, 17, 9, 5, 3, 2, 1, 1, 1, 1, 2, 1, 1, 1, 4, 2, 1, 1, 1, 
2, 1, 1, 8, 4, 2, 1, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 16, 8, 4, 2, 1, 1, 1, 2, 
1, 1, 4, 2, 1, 1, 2, 1, 1, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1 };
           
// Setting up macro defines and including GF(p), GF(p^2), curve, isogeny and kex functions
#define fpcopy                        fpcopy503
#define fpzero                        fpzero503
#define fpadd                         fpadd503
#define fpsub                         fpsub503
#define fpneg                         fpneg503
#define fpdiv2                        fpdiv2_503
#define fpcorrection                  fpcorrection503
#define fpmul_mont                    fpmul503_mont
#define fpsqr_mont                    fpsqr503_mont
#define fpinv_mont                    fpinv503_mont
#define fpinv_chain_mont              fpinv503_chain_mont
#define fpinv_mont_bingcd             fpinv503_mont_bingcd
#define fp2copy                       fp2copy503
#define fp2zero                       fp2zero503
#define fp2add                        fp2add503
#define fp2sub                        fp2sub503
#define mp_sub_p2                     mp_sub503_p2
#define mp_sub_p4                     mp_sub503_p4
#define sub_p4                        mp_sub_p4
#define fp2neg                        fp2neg503
#define fp2div2                       fp2div2_503
#define fp2correction                 fp2correction503
#define fp2mul_mont                   fp2mul503_mont
#define fp2sqr_mont                   fp2sqr503_mont
#define fp2inv_mont                   fp2inv503_mont
#define fp2inv_mont_bingcd            fp2inv503_mont_bingcd
#define fpequal_non_constant_time     fpequal503_non_constant_time
#define mp_add_asm                    mp_add503_asm
#define mp_subaddx2_asm               mp_subadd503x2_asm
#define mp_dblsubx2_asm               mp_dblsub503x2_asm
#define crypto_kem_keypair            crypto_kem_keypair_SIKEp503
#define crypto_kem_enc                crypto_kem_enc_SIKEp503
#define crypto_kem_dec                crypto_kem_dec_SIKEp503
#define random_mod_order_A            random_mod_order_A_SIDHp503
#define random_mod_order_B            random_mod_order_B_SIDHp503
#define EphemeralKeyGeneration_A      EphemeralKeyGeneration_A_SIDHp503
#define EphemeralKeyGeneration_B      EphemeralKeyGeneration_B_SIDHp503
#define EphemeralSecretAgreement_A    EphemeralSecretAgreement_A_SIDHp503
#define EphemeralSecretAgreement_B    EphemeralSecretAgreement_B_SIDHp503

#include "../fpx.c"
#include "../ec_isogeny.c"
#include "../sidh.c"    
#include "../sike.c"