/********************************************************************************************
* SIDH: an efficient supersingular isogeny cryptography library
*
* Abstract: supersingular isogeny parameters and generation of functions for P434
*********************************************************************************************/  

#include "P434_api.h" 
#include "P434_internal.h"
#include "../internal.h"


// Encoding of field elements, elements over Z_order, elements over GF(p^2) and elliptic curve points:
// --------------------------------------------------------------------------------------------------
// Elements over GF(p) and Z_order are encoded with the least significant octet (and digit) located at the leftmost position (i.e., little endian format). 
// Elements (a+b*i) over GF(p^2), where a and b are defined over GF(p), are encoded as {a, b}, with a in the least significant position.
// Elliptic curve points P = (x,y) are encoded as {x, y}, with x in the least significant position. 
// Internally, the number of digits used to represent all these elements is obtained by approximating the number of bits to the immediately greater multiple of 32.
// For example, a 434-bit field element is represented with Ceil(434 / 64) = 7 64-bit digits or Ceil(434 / 32) = 14 32-bit digits.

//
// Curve isogeny system "SIDHp434". Base curve: Montgomery curve By^2 = Cx^3 + Ax^2 + Cx defined over GF(p434^2), where A=6, B=1, C=1 and p434 = 2^216*3^137-1
//
         
const uint64_t p434[NWORDS64_FIELD]              = { 0xFFFFFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFF, 0xFDC1767AE2FFFFFF, 
                                                     0x7BC65C783158AEA3, 0x6CFC5FD681C52056, 0x0002341F27177344 };
const uint64_t p434x2[NWORDS64_FIELD]            = { 0xFFFFFFFFFFFFFFFE, 0xFFFFFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFF, 0xFB82ECF5C5FFFFFF,
                                                     0xF78CB8F062B15D47, 0xD9F8BFAD038A40AC, 0x0004683E4E2EE688 }; 
const uint64_t p434x4[NWORDS64_FIELD]            = { 0xFFFFFFFFFFFFFFFC, 0xFFFFFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFF, 0xF705D9EB8BFFFFFF, 
                                                     0xEF1971E0C562BA8F, 0xB3F17F5A07148159, 0x0008D07C9C5DCD11 }; 
const uint64_t p434p1[NWORDS64_FIELD]            = { 0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0xFDC1767AE3000000,
                                                     0x7BC65C783158AEA3, 0x6CFC5FD681C52056, 0x0002341F27177344 };  
const uint64_t p434x16p[2*NWORDS64_FIELD]        = { 0x0000000000000010, 0x0000000000000000, 0x0000000000000000, 0x47D130A3A0000000, 
                                                     0x873470F9D4EA2B80, 0x6074052FC75BF530, 0x54497C1B1D119772, 0xC55F373D2CDCA412, 
                                                     0x732CA2221C664B96, 0x6445AB96AF6359A5, 0x221708AB42ABE1B4, 0xAE3D3D0063244F01, 
                                                     0x18B920F2ECF68816, 0x0000004DB194809D }; 
// Order of Alice's subgroup
const uint64_t Alice_order[NWORDS64_ORDER]       = { 0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000001000000 }; 
// Order of Bob's subgroup
const uint64_t Bob_order[NWORDS64_ORDER]         = { 0x58AEA3FDC1767AE3, 0xC520567BC65C7831, 0x1773446CFC5FD681, 0x0000000002341F27 };
// Alice's generator values {XPA0 + XPA1*i, XQA0 + xQA1*i, XRA0 + XRA1*i} in GF(p434^2), expressed in Montgomery representation
const uint64_t A_gen[6 * NWORDS64_FIELD]         = { 0x6E18D3A63313A738, 0x1DCC496DD6DDE298, 0xA35F3F7DAFBE2B43, 0xC6B9A5CC670071EB,
                                                     0x2EA3DB085283675A, 0xFDFE173A0297F36,  0x0002200804EB824D, // XPA0
                                                     0xB999E9E259F7BFA8, 0x2584D67D0C2EEAA9, 0x80AB07D4E9625724, 0x781DA616A7A76E54,
                                                     0x9BE449736374F491, 0x8C6F86E8B0C4D74A, 0x0001C1D4812CBD98, // XPA1
                                                     0x257DBD53095FD263, 0xBBB3C7A7B4EDB1D4, 0xA817B7FDDD5BB8DA, 0xF5DE963B242B7AB3,
                                                     0x7F51B5362FC94CB6, 0xE7D2496B526DFF16, 0x0001E962CF69118C, // XQA0
                                                     0xED9DC89467FB039D, 0x17C71E114B5803D0, 0x816C3379BE9647BF, 0xB07F441A15434B64,
                                                     0xCC65C1804AF4CBD1, 0xF06BF5F074032C77, 0x0001A251F94CF02C, // XQA1
                                                     0xA26194AB4BD1A16F, 0xCFCD9F7F04D5AB10, 0x1BB4A7C04C37482C, 0x71DEE733632DA36D,
                                                     0x7335784B5ECF957F, 0x66AE2381533A7F09, 0x000232BFFE6FA42F, // XRA0
                                                     0x60ACBE5D899CFA6A, 0x82AC55A556E5A22F, 0x437D8C2AC83FDC6B, 0x620A8DA602543EDE,
                                                     0xD19ABA8092A1E8C2, 0xAFF1AA61981C95D3, 0x0001A7232B0C035E }; // XRA1
const uint64_t DBL_QA[2 * NWORDS64_FIELD]        = { 0x6448cb5bd976250e, 0x3092cf8c8803d8b4, 0x2542331c81c2c2f8, 0xc086e30db24aa32,
                                                     0xd1b7f3e5532ee315, 0x226067da42cd56db, 0x1dea86eb48d8a,
                                                     0x414795f6a70d543e, 0xdc7d2e0b9229a814, 0x679cd711c5b2ac8e, 0xcc92a943030f0d18,
                                                     0xc2bb3cda074e0076, 0x19332e71dc423ba4, 0x1c95a2fc045dd };
const uint64_t YPA[2 * NWORDS64_FIELD]           = { 0x9b989be60cff0d15, 0x8b80a32171813f53, 0xf4f067606a56228e, 0x48f8237e159577b0,
                                                     0x42529574b9e74156, 0xd8d26313f4aa9f9c, 0x1279ac6bc876c,
                                                     0x9597544cbe9d88df, 0x13801f440df32748, 0xe4ecaff9c15d0ceb, 0x7867d92eb045a646,
                                                     0x2399062ba8c64ef, 0xe9258c0bdf8bbff7, 0x1ce4bbf872205 };
const uint64_t PplusQA[4 * NWORDS64_FIELD]       = { 0xa26194ab4bd1a16f, 0xcfcd9f7f04d5ab10, 0x1bb4a7c04c37482c, 0x71dee733632da36d,
                                                     0x7335784b5ecf957f, 0x66ae2381533a7f09, 0x232bffe6fa42f,
                                                     0x9f5341a276630595, 0x7d53aa5aa91a5dd0, 0xbc8273d537c02394, 0x9bb6e8d4e0abc121,
                                                     0xaa2ba1f79eb6c5e1, 0xbd0ab574e9a88a82, 0x8cfbfc0b6fe5,
                                                     0xda1e960ae3e5c4e8, 0x6effc0350686260f, 0xbd6eaccf62467b5, 0x65939cfb1161e478,
                                                     0x5d0ed5901e82ddcb, 0xda8be7ee6d455d94, 0x19017b8ce77b2 ,
                                                     0xd30ece1ea3e19f40, 0x3abb724e9467b8fd, 0xc34cec4a1f9f85d4, 0xe0b40f984e683dc0,
                                                     0x684c9b19b4180b6e, 0x7314c90c41f2842e, 0xe0745aab36b4 };
// Bob's generator values {XPB0, XQB0, XRB0 + XRB1*i} in GF(p434^2), expressed in Montgomery representation
const uint64_t B_gen[6 * NWORDS64_FIELD]         = { 0xE172658571249BA8, 0x9D8F52CB15829DA0, 0xE3A7C7F9F0E3F832, 0x8B825DD0B9410D30,
                                                     0xF42F815734752EDA, 0xCB35DD9160997586, 0x00018B3AAAAD0F79, // XPB0
                                                     0xCF0B435C40C1375D, 0x58AC8A63992B36EF, 0x416D0B3DFB0C1DF5, 0xB257E9CFE8985F15,
                                                     0xA493D98A7A1D6DF2, 0x6D6781A5B3FDE61F, 0x000179AC0D886A3F, // XPB1
                                                     0xE172658571249BA8, 0x9D8F52CB15829DA0, 0xE3A7C7F9F0E3F832, 0x8B825DD0B9410D30,
                                                     0xF42F815734752EDA, 0xCB35DD9160997586, 0x00018B3AAAAD0F79, // XQB0
                                                     0x30F4BCA3BF3EC8A2, 0xA753759C66D4C910, 0xBE92F4C204F3E20A, 0x4B698CAAFA67A0EA,
                                                     0xD73282EDB73B40B1, 0xFF94DE30CDC73A36, 0x0000BA73198F0904, // XQB1
                                                     0x9F7367022EFDF650, 0xA8C21C687A91D6BC, 0xDDB909C497C4BFED, 0x66FD362A30232EBF,
                                                     0x84AC5026408590E1, 0x5378004CB74DA4ED, 0x00008AA46B9E55B2, // XRB0
                                                     0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0 }; // XRB1
/* Basis for Bob on E_0, expressed in Montgomery representation */
const uint64_t P3[2 * NWORDS64_FIELD]            = { 0x214C34BB192F67A0, 0xDD49D3D02115D30, 0x700652C1A7B66ED, 0x1F856B48F4FF0024,
                                                     0xFBDE6F4E6A705221, 0xB951A3D6C93D87B8, 0xAE8ADB818ED6,
                                                     0x51D889FE197209C1, 0x191BCD9DBE4FE0EF, 0x447818CF5E54DD8A, 0x3F42710E8562A583,
                                                     0x647BDBB01C66DCB5, 0xF402D36C15EA12E1, 0xA1E1D287C14C};


// Montgomery constant Montgomery_R2 = (2^448)^2 mod p434
const uint64_t Montgomery_R2[NWORDS64_FIELD]     = { 0x28E55B65DCD69B30, 0xACEC7367768798C2, 0xAB27973F8311688D, 0x175CC6AF8D6C7C0B,
                                                     0xABCD92BF2DDE347E, 0x69E16A61C7686D9A, 0x000025A89BCDD12A };                                                   
// Value one in Montgomery representation 
const uint64_t Montgomery_one[NWORDS64_FIELD]    = { 0x000000000000742C, 0x0000000000000000, 0x0000000000000000, 0xB90FF404FC000000, 
                                                     0xD801A4FB559FACD4, 0xE93254545F77410C, 0x0000ECEEA7BD2EDA };


// Fixed parameters for isogeny tree computation
const unsigned int strat_Alice[MAX_Alice-1] = { 
48, 28, 16, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 13, 7, 4, 2, 1, 1, 2, 1, 1, 3, 2, 1, 1, 
1, 1, 5, 4, 2, 1, 1, 2, 1, 1, 2, 1, 1, 1, 21, 12, 7, 4, 2, 1, 1, 2, 1, 1, 3, 2, 1, 1, 1, 1, 5, 3, 2, 1, 1, 1, 1, 2, 1, 1, 1, 9, 5, 3, 2, 1, 1, 
1, 1, 2, 1, 1, 1, 4, 2, 1, 1, 1, 2, 1, 1 };

const unsigned int strat_Bob[MAX_Bob-1] = { 
66, 33, 17, 9, 5, 3, 2, 1, 1, 1, 1, 2, 1, 1, 1, 4, 2, 1, 1, 1, 2, 1, 1, 8, 4, 2, 1, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 16, 8, 4, 2, 1, 1, 1, 
2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 32, 16, 8, 4, 3, 1, 1, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 8, 4, 2, 
1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 16, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1 };
           
// Setting up macro defines and including GF(p), GF(p^2), curve, isogeny and kex functions
#define fpcopy                        fpcopy434
#define fpzero                        fpzero434
#define fpadd                         fpadd434
#define fpsub                         fpsub434
#define fpneg                         fpneg434
#define fpdiv2                        fpdiv2_434
#define fpcorrection                  fpcorrection434
#define fpmul_mont                    fpmul434_mont
#define fpsqr_mont                    fpsqr434_mont
#define fpinv_mont                    fpinv434_mont
#define fpinv_chain_mont              fpinv434_chain_mont
#define fpinv_mont_bingcd             fpinv434_mont_bingcd
#define fp2copy                       fp2copy434
#define fp2zero                       fp2zero434
#define fp2add                        fp2add434
#define fp2sub                        fp2sub434
#define mp_sub_p2                     mp_sub434_p2
#define mp_sub_p4                     mp_sub434_p4
#define sub_p4                        mp_sub_p4
#define fp2neg                        fp2neg434
#define fp2div2                       fp2div2_434
#define fp2correction                 fp2correction434
#define fp2mul_mont                   fp2mul434_mont
#define fp2sqr_mont                   fp2sqr434_mont
#define fp2inv_mont                   fp2inv434_mont
#define fp2inv_mont_bingcd            fp2inv434_mont_bingcd
#define fpequal_non_constant_time     fpequal434_non_constant_time
#define mp_add_asm                    mp_add434_asm
#define mp_subaddx2_asm               mp_subadd434x2_asm
#define mp_dblsubx2_asm               mp_dblsub434x2_asm
#define crypto_kem_keypair            crypto_kem_keypair_SIKEp434
#define crypto_kem_enc                crypto_kem_enc_SIKEp434
#define crypto_kem_dec                crypto_kem_dec_SIKEp434
#define random_mod_order_A            random_mod_order_A_SIDHp434
#define random_mod_order_B            random_mod_order_B_SIDHp434
#define EphemeralKeyGeneration_A      EphemeralKeyGeneration_A_SIDHp434
#define EphemeralKeyGeneration_B      EphemeralKeyGeneration_B_SIDHp434
#define EphemeralSecretAgreement_A    EphemeralSecretAgreement_A_SIDHp434
#define EphemeralSecretAgreement_B    EphemeralSecretAgreement_B_SIDHp434

#include "../fpx.c"
#include "../ec_isogeny.c"
#include "../sidh.c"    
#include "../sike.c"