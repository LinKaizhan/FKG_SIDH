/********************************************************************************************
* SIDH: an efficient supersingular isogeny cryptography library
*
* Abstract: supersingular isogeny parameters and generation of functions for P610
*********************************************************************************************/  

#include "P610_api.h" 
#include "P610_internal.h"
#include "../internal.h"


// Encoding of field elements, elements over Z_order, elements over GF(p^2) and elliptic curve points:
// --------------------------------------------------------------------------------------------------
// Elements over GF(p) and Z_order are encoded with the least significant octet (and digit) located at the leftmost position (i.e., little endian format). 
// Elements (a+b*i) over GF(p^2), where a and b are defined over GF(p), are encoded as {a, b}, with a in the least significant position.
// Elliptic curve points P = (x,y) are encoded as {x, y}, with x in the least significant position. 
// Internally, the number of digits used to represent all these elements is obtained by approximating the number of bits to the immediately greater multiple of 32.
// For example, a 610-bit field element is represented with Ceil(610 / 64) = 10 64-bit digits or Ceil(610 / 32) = 20 32-bit digits.

//
// Curve isogeny system "SIDHp610". Base curve: Montgomery curve By^2 = Cx^3 + Ax^2 + Cx defined over GF(p610^2), where A=6, B=1, C=1 and p610 = 2^305*3^192-1
//
         
const uint64_t p610[NWORDS64_FIELD]              = { 0xFFFFFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFF, 0x6E01FFFFFFFFFFFF, 
                                                     0xB1784DE8AA5AB02E, 0x9AE7BF45048FF9AB, 0xB255B2FA10C4252A, 0x819010C251E7D88C, 0x000000027BF6A768 };
const uint64_t p610x2[NWORDS64_FIELD]            = { 0xFFFFFFFFFFFFFFFE, 0xFFFFFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFF, 0xDC03FFFFFFFFFFFF,
                                                     0x62F09BD154B5605C, 0x35CF7E8A091FF357, 0x64AB65F421884A55, 0x03202184A3CFB119, 0x00000004F7ED4ED1 };
const uint64_t p610x4[NWORDS64_FIELD]            = { 0xFFFFFFFFFFFFFFFC, 0xFFFFFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFF, 0xB807FFFFFFFFFFFF, 
                                                     0xC5E137A2A96AC0B9, 0x6B9EFD14123FE6AE, 0xC956CBE8431094AA, 0x06404309479F6232, 0x00000009EFDA9DA2 };
const uint64_t p610p1[NWORDS64_FIELD]            = { 0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x6E02000000000000,
                                                     0xB1784DE8AA5AB02E, 0x9AE7BF45048FF9AB, 0xB255B2FA10C4252A, 0x819010C251E7D88C, 0x000000027BF6A768 };   
const uint64_t p610x16p[2*NWORDS64_FIELD]        = { 0x0000000000000010, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x3FC0000000000000, 
                                                     0xD0F642EAB4A9FA32, 0xA308175F6E00CA89, 0xB549A0BDE77B5AAC, 0xCDFDE7B5C304EE69, 0x7FDB7FF0812B12EF, 
                                                     0xE09BA529B9FE1167, 0xD249C196DAB8CD7F, 0xD4E22754A3F20928, 0x97825638B19A7CCE, 0x05E04550FC4CCE0D, 
                                                     0x8FB5DA1152CDE50C, 0xF9649BA3EA408644, 0x4473C93E6441063D, 0xBE190269D1337B7B, 0x0000000000000062 }; 
// Order of Alice's subgroup
const uint64_t Alice_order[NWORDS64_ORDER]       = { 0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0002000000000000 }; 
// Order of Bob's subgroup
const uint64_t Bob_order[NWORDS64_ORDER]         = { 0x26F4552D58173701, 0xDFA28247FCD5D8BC, 0xD97D086212954D73, 0x086128F3EC46592A, 0x00013DFB53B440C8 };
// Alice's generator values {XPA0 + XPA1*i, XQA0 + xQA1*i, XRA0 + XRA1*i} in GF(p610^2), expressed in Montgomery representation
const uint64_t A_gen[6 * NWORDS64_FIELD]         = { 0x31C8AF7FFC0DE9FA, 0x8A8AD55D2AC8A709, 0x95A4DC49B64E5B2C, 0xF08C77AAE90ABE83, 0x675E4FF97C95845D,
                                                     0xF8A22591248401F0, 0x73F573A4FF34A84A, 0x37D18A6C3D989158, 0xEE73973862A3E95, 0x24084FCCB, // XPA0
                                                     0x4B8C9CED6DEF0B8B, 0x652C800D926AB992, 0x3DFA6D6B8FD37D80, 0xA30C578CD98EFD79, 0x9FC067E58CCBD32E,
                                                     0x2B0599AEAF150FDB, 0xBA321B31886F3292, 0xE0011F56247547A1, 0x28CA0747910BFAE2, 0xFC020A14, // XPA1
                                                     0x2728178178DEAFBD, 0xD377C4656DBC71F0, 0x968642007B807932, 0xB8B04B1039062A21, 0xF824771B468A977C,
                                                     0x260F1C50354F46AB, 0x78A3D37CDBBD4DC5, 0x1FB1BAC6851BA175, 0xA73444F1CAC4A10, 0xF3A5C2BB, // XQA0
                                                     0x4F828B752E825BB4, 0x82CEA210AC766C69, 0x8B1BBC87DAD8BEDD, 0x9BFC5B9CE215B423, 0xF7E1BCC0C541177C,
                                                     0x7727E3A0F1A1AF24, 0xFBCFE4177D2B0221, 0xBB15BDCC160D902A, 0x3FE1467B4A911446, 0x1A495CB35, // XQA1
                                                     0x38687702D78D1A93, 0x58C09FD23B1E1B56, 0xC54917327D5C0FAB, 0xB6D55B7BE801A3C, 0xEB3AE21C8B93E9E9,
                                                     0xECB45AD6D24FF76A, 0x850645B4F39EC5F2, 0xE6F78202586C9B3A, 0x2923209A250F7F66, 0x26FB150F, // XRA0
                                                     0x5AC7B27F9096F718, 0x487DDD2820132C83, 0x6B21AC48569E12D8, 0x57B54E5A827D1CD9, 0xDB7C4BEB143E4130,
                                                     0xB6781CA1DA245EAD, 0xCC09878A2A6D7C45, 0x980726C5232C75E5, 0x50D3A7350792C35F, 0x172B595DB }; // XRA1
const uint64_t DBL_QA[2 * NWORDS64_FIELD]        = { 0x2c9e52fa31b9b76, 0xef4088ad3e54c6dd, 0xc18e7055d2cff348, 0x24b3268c87d5f690, 0xdd80d94ade7b0a93,
                                                     0x5ee075b1e9a6c6bd, 0x27f68f76241404bd, 0x2b267148416a9627, 0x27270dafd0dd30ff, 0xd8b7e841,
                                                     0xb84b1e242a63879e, 0xd3a74c3d2770fa06, 0x49df32c277de73a3, 0xca452cb04eba1741, 0xe36ec74b21763cf,
                                                     0xee808c414124f7b3, 0xcdbc7c4c7fa2f565, 0x6ec6a04436a3b6dd, 0x655153fcac56e490, 0x1c8ae36cc };
const uint64_t YPA[2 * NWORDS64_FIELD]           = { 0x3ca84837d69d8728, 0xb2bdfe3304cb7401, 0x8c840937950ad3e9, 0xce8094a539aa6c49, 0xf0802aae490f29a0, 
                                                     0x5458a8e61bb9d01f, 0x3592a73de4758511, 0x7dea75b85a60f316, 0xf835eeac9b12cc1d, 0x11c4e0162,
                                                     0x87a90900552b058, 0xf34899fe9411dc6a, 0x3807cf5b95b0168, 0xc986baf1e3ffded4, 0x1d10eac33aa0781a,
                                                     0xd9569230f9a2d512, 0xf8295f6189dbaaf3, 0x26b44d4cecb1a5e8, 0x9ca4ce754143daa4, 0xaf517dfc };
const uint64_t PplusQA[4 * NWORDS64_FIELD]       = { 0x3ae8bb4f4de5e21a, 0x3c646dfde429c031, 0x9dc5916c37a21fc4, 0x1754faf5d9dc1ba3, 0x53f7022de9e07850, 
                                                     0x97ac6836c73d072e, 0x26e37b2502a716d7, 0x3e643c9018eca8e5, 0xc796641a6ee9017f, 0x12e6a48ec,
                                                     0x6bff31eccb1e2092, 0xd916e73e07769500, 0xedf799cc675ee22b, 0xb3c36ed05b36434, 0x629758b74e92643e, 
                                                     0x3e35456235455243, 0x87624a13997758c7, 0xfde1837097e7d59b, 0x6eeffed35309078b, 0xb12d2a52,
                                                     0xeabb251c79a581c9, 0x1ddadde7d0c4adbf, 0x3979ea0e826c6034, 0xefeb3adf3ea1a68c, 0x174c6cd565164f3,
                                                     0x971a26fbfb9544bc, 0x83ce13424dc2d699, 0xaeb453e747a11622, 0x23dcc826e38ff746, 0xd2346570,
                                                     0x2ec71192464d8b22, 0x3fd75abed41d8c72, 0x2e206d4f17f372db, 0xd91f67a83c6616ec, 0x268b00035db0c31,
                                                     0xccb96bd0238db8cf, 0x71cc72e3696eb8e7, 0x83599d21e5430d78, 0x55416a92cdf519d0, 0x1c23f7a19 };
// Bob's generator values {XPB0, XQB0, XRB0 + XRB1*i} in GF(p610^2), expressed in Montgomery representation
const uint64_t B_gen[6 * NWORDS64_FIELD]         = { 0xD4A2CF040BC56F2C, 0x58F1D1D2B190EDE7, 0x2229F10D3BC7BA47, 0x769AB0F0EDD86AA4, 0x97F1214B80D8463,
                                                     0x9B23774D13ED3EEE, 0x9A182E846DAA95C6, 0x343741369B273442, 0x61FB37462569D4BB, 0x1815EF8B9, // XPB0
                                                     0xF380CA27C26BF32E, 0xD594C3EA0698D298, 0x21D388E632D1CA2E, 0xDD1E0B34330E0AB0, 0xEA7B89CAD59CA8C2,
                                                     0x28C129BFC584BEC1, 0x48D1E802FC7418CF, 0x11F3A548C5DFFDF7, 0xDB0E9AF98D314F67, 0x219918D2B, // XPB1
                                                     0xD4A2CF040BC56F2C, 0x58F1D1D2B190EDE7, 0x2229F10D3BC7BA47, 0x769AB0F0EDD86AA4, 0x97F1214B80D8463,
                                                     0x9B23774D13ED3EEE, 0x9A182E846DAA95C6, 0x343741369B273442, 0x61FB37462569D4BB, 0x1815EF8B9, // XQB0
                                                     0xC7F35D83D940CD1, 0x2A6B3C15F9672D67, 0xDE2C7719CD2E35D1, 0x22E1F4CBCCF1F54F, 0x838676352A63573D,
                                                     0x88B72428E4D5F16C, 0x5215D742081BE0DC, 0xA0620DB14AE42733, 0xA68175C8C4B68925, 0x62651A3C, // XQB1
                                                     0x4F62205A5DAFB369, 0xA2B75D5BC06C691F, 0x6B82C9B893D51C38, 0x2C2467D7AB7DAA2C, 0x8A8D5AC13C2C5ADD,
                                                     0xBC3AEC544F8953F5, 0xBC43C1BE1B1DC069, 0xB8CDA0908AEBCD84, 0xA213356DB0FBFCFF, 0x15F063030, // XRB0
                                                     0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0 }; // XRB1
/* Basis for Bob on A = 0, expressed in Montgomery representation */
const uint64_t P3[2 * NWORDS64_FIELD]            = { 0x203596CF0245B227, 0xFE7D4CB978F11517, 0xEC79574E9D7DD13A, 0xD24627B69D4DFF63, 0x85B4D3B2B5426BBF,
                                                     0xFF0237C357683FCA, 0x2C3E0FE7792534B1, 0x8B68DB1AFC3F9CDE, 0x5AFD2B5021786921, 0x16CFF1918,
                                                     0xDFE1CAFF47350FFB, 0x7F6641B5806DBD07, 0xD558CE2B43292C47, 0x28EB4A4147C77BD6, 0x143218EB29F5FB6C,
                                                     0x5F457BD167A2260F, 0x26D9639E9DD4A15D, 0xEC9DFA3764433777, 0x9D8C59E2D257CACF, 0x1D2D65779};
// Montgomery constant Montgomery_R2 = (2^640)^2 mod p610
const uint64_t Montgomery_R2[NWORDS64_FIELD]     = { 0xE75F5D201A197727, 0xE0B85963B627392E, 0x6BC1707818DE493D, 0xDC7F419940D1A0C5, 0x7358030979EDE54A,
                                                     0x84F4BEBDEED75A5C, 0x7ECCA66E13427B47, 0xC5BB4E65280080B3, 0x7019950F516DA19A, 0x000000008E290FF3 };                                                    
// Value one in Montgomery representation 
const uint64_t Montgomery_one[NWORDS64_FIELD]    = { 0x00000000670CC8E6, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x9A34000000000000,
                                                     0x4D99C2BD28717A3F, 0x0A4A1839A323D41C, 0xD2B62215D06AD1E2, 0x1369026E862CAF3D, 0x000000010894E964 };


// Fixed parameters for isogeny tree computation
const unsigned int strat_Alice[MAX_Alice-1] = { 
67, 37, 21, 12, 7, 4, 2, 1, 1, 2, 1, 1, 3, 2, 1, 1, 1, 1, 5, 3, 2, 1, 1, 1, 1, 
2, 1, 1, 1, 9, 5, 3, 2, 1, 1, 1, 1, 2, 1, 1, 1, 4, 2, 1, 1, 1, 2, 1, 1, 16, 9, 
5, 3, 2, 1, 1, 1, 1, 2, 1, 1, 1, 4, 2, 1, 1, 1, 2, 1, 1, 8, 4, 2, 1, 1, 2, 1, 
1, 4, 2, 1, 1, 2, 1, 1, 33, 16, 8, 5, 2, 1, 1, 1, 2, 1, 1, 1, 4, 2, 1, 1, 2, 1, 
1, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 16, 8, 4, 2, 1, 1, 1, 2, 1, 1, 
4, 2, 1, 1, 2, 1, 1, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1 };

const unsigned int strat_Bob[MAX_Bob-1] = { 
86, 48, 27, 15, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 7, 4, 2, 1, 1, 2, 1, 
1, 3, 2, 1, 1, 1, 1, 12, 7, 4, 2, 1, 1, 2, 1, 1, 3, 2, 1, 1, 1, 1, 5, 3, 2, 1, 1, 
1, 1, 2, 1, 1, 1, 21, 12, 7, 4, 2, 1, 1, 2, 1, 1, 3, 2, 1, 1, 1, 1, 5, 3, 2, 1, 1, 
1, 1, 2, 1, 1, 1, 9, 5, 3, 2, 1, 1, 1, 1, 2, 1, 1, 1, 4, 2, 1, 1, 1, 2, 1, 1, 38, 
21, 12, 7, 4, 2, 1, 1, 2, 1, 1, 3, 2, 1, 1, 1, 1, 5, 3, 2, 1, 1, 1, 1, 2, 1, 1, 1, 
9, 5, 3, 2, 1, 1, 1, 1, 2, 1, 1, 1, 4, 2, 1, 1, 1, 2, 1, 1, 17, 9, 5, 3, 2, 1, 1, 
1, 1, 2, 1, 1, 1, 4, 2, 1, 1, 1, 2, 1, 1, 8, 4, 2, 1, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 
1, 1 };

// Setting up macro defines and including GF(p), GF(p^2), curve, isogeny and kex functions
#define fpcopy                        fpcopy610
#define fpzero                        fpzero610
#define fpadd                         fpadd610
#define fpsub                         fpsub610
#define fpneg                         fpneg610
#define fpdiv2                        fpdiv2_610
#define fpcorrection                  fpcorrection610
#define fpmul_mont                    fpmul610_mont
#define fpsqr_mont                    fpsqr610_mont
#define fpinv_mont                    fpinv610_mont
#define fpinv_chain_mont              fpinv610_chain_mont
#define fpinv_mont_bingcd             fpinv610_mont_bingcd
#define fp2copy                       fp2copy610
#define fp2zero                       fp2zero610
#define fp2add                        fp2add610
#define fp2sub                        fp2sub610
#define mp_sub_p2                     mp_sub610_p2
#define mp_sub_p4                     mp_sub610_p4
#define sub_p4                        mp_sub_p4
#define fp2neg                        fp2neg610
#define fp2div2                       fp2div2_610
#define fp2correction                 fp2correction610
#define fp2mul_mont                   fp2mul610_mont
#define fp2sqr_mont                   fp2sqr610_mont
#define fp2inv_mont                   fp2inv610_mont
#define fp2inv_mont_bingcd            fp2inv610_mont_bingcd
#define fpequal_non_constant_time     fpequal610_non_constant_time
#define mp_add_asm                    mp_add610_asm
#define mp_subaddx2_asm               mp_subadd610x2_asm
#define mp_dblsubx2_asm               mp_dblsub610x2_asm
#define crypto_kem_keypair            crypto_kem_keypair_SIKEp610
#define crypto_kem_enc                crypto_kem_enc_SIKEp610
#define crypto_kem_dec                crypto_kem_dec_SIKEp610
#define random_mod_order_A            random_mod_order_A_SIDHp610
#define random_mod_order_B            random_mod_order_B_SIDHp610
#define EphemeralKeyGeneration_A      EphemeralKeyGeneration_A_SIDHp610
#define EphemeralKeyGeneration_B      EphemeralKeyGeneration_B_SIDHp610
#define EphemeralSecretAgreement_A    EphemeralSecretAgreement_A_SIDHp610
#define EphemeralSecretAgreement_B    EphemeralSecretAgreement_B_SIDHp610

#include "../fpx.c"
#include "../ec_isogeny.c"
#include "../sidh.c"
#include "../sike.c"