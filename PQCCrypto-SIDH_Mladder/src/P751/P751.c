/********************************************************************************************
* SIDH: an efficient supersingular isogeny cryptography library
*
* Abstract: supersingular isogeny parameters and generation of functions for P751
*********************************************************************************************/  

#include "P751_api.h" 
#include "P751_internal.h"
#include "../internal.h"


// Encoding of field elements, elements over Z_order, elements over GF(p^2) and elliptic curve points:
// --------------------------------------------------------------------------------------------------
// Elements over GF(p) and Z_order are encoded with the least significant octet (and digit) located at the leftmost position (i.e., little endian format). 
// Elements (a+b*i) over GF(p^2), where a and b are defined over GF(p), are encoded as {a, b}, with a in the least significant position.
// Elliptic curve points P = (x,y) are encoded as {x, y}, with x in the least significant position. 
// Internally, the number of digits used to represent all these elements is obtained by approximating the number of bits to the immediately greater multiple of 32.
// For example, a 751-bit field element is represented with Ceil(751 / 64) = 12 64-bit digits or Ceil(751 / 32) = 24 32-bit digits.

//
// Curve isogeny system "SIDHp751". Base curve: Montgomery curve By^2 = Cx^3 + Ax^2 + Cx defined over GF(p751^2), where A=6, B=1, C=1 and p751 = 2^372*3^239-1
//
         
const uint64_t p751[NWORDS64_FIELD]              = { 0xFFFFFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFF, 0xEEAFFFFFFFFFFFFF,
                                                     0xE3EC968549F878A8, 0xDA959B1A13F7CC76, 0x084E9867D6EBE876, 0x8562B5045CB25748, 0x0E12909F97BADC66, 0x00006FE5D541F71C };
const uint64_t p751x2[NWORDS64_FIELD]            = { 0xFFFFFFFFFFFFFFFE, 0xFFFFFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFF, 0xDD5FFFFFFFFFFFFF, 
                                                     0xC7D92D0A93F0F151, 0xB52B363427EF98ED, 0x109D30CFADD7D0ED, 0x0AC56A08B964AE90, 0x1C25213F2F75B8CD, 0x0000DFCBAA83EE38 }; 
const uint64_t p751x4[NWORDS64_FIELD]            = { 0xFFFFFFFFFFFFFFFC, 0xFFFFFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFF, 0xFFFFFFFFFFFFFFFF, 0xBABFFFFFFFFFFFFF, 
                                                     0x8FB25A1527E1E2A3, 0x6A566C684FDF31DB, 0x213A619F5BAFA1DB, 0x158AD41172C95D20, 0x384A427E5EEB719A, 0x0001BF975507DC70 }; 
const uint64_t p751p1[NWORDS64_FIELD]            = { 0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0xEEB0000000000000,
                                                     0xE3EC968549F878A8, 0xDA959B1A13F7CC76, 0x084E9867D6EBE876, 0x8562B5045CB25748, 0x0E12909F97BADC66, 0x00006FE5D541F71C };   
const uint64_t p751x16p[2*NWORDS64_FIELD]        = { 0x0000000000000010, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x2A00000000000000, 
                                                     0x826D2F56C0F0EAE2, 0xAD4C9CBD81067123, 0xF62CF3052282F124, 0x53A95F7469B516FE, 0x3DADEC0D08A4732F, 0x58AD934557C11C7E, 
                                                     0x7F731B89B2DA43F2, 0x51AE9F5F5F6AFF3B, 0xD74319A6C9BCA375, 0x5BAB790796CF84D4, 0xA421554FE2E49CA8, 0x20AD617C8DF437CF, 
                                                     0x3AB06E7A12F5FF7B, 0x70A25E037E40347E, 0x51F1D323FB4C1151, 0xAE0D99AA4835FED9, 0xDF5429960D2536B6, 0x000000030E91D466 };
// Order of Alice's subgroup
const uint64_t Alice_order[NWORDS64_ORDER]       = { 0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0010000000000000 }; 
// Order of Bob's subgroup
const uint64_t Bob_order[NWORDS64_ORDER]         = { 0xC968549F878A8EEB, 0x59B1A13F7CC76E3E, 0xE9867D6EBE876DA9, 0x2B5045CB25748084, 0x2909F97BADC66856, 0x06FE5D541F71C0E1 };
// Alice's generator values {XPA0 + XPA1*i, XQA0 + xQA1*i, XRA0 + XRA1*i} in GF(p751^2), expressed in Montgomery representation
const uint64_t A_gen[6 * NWORDS64_FIELD]         = { 0x2584350E0C33C304, 0x51E9C29E234DC61E, 0xC6E65A7BF90ACC05, 0xB1333E2E19B3A930, 0xA4F7CA2F7F66909F, 0xE01E9E6F6704BF9E,
                                                     0xE2345D48C0219D6D, 0x70F37AD9933FC182, 0x7B9D4D5870CFACA3, 0x3B8DAF20190D460D, 0xB02D6FF9AAFA0C7, 0x15A435D19526, // XPA0
                                                     0xE85E3F2B4EDDAF22, 0x4824EDCA0A253CB2, 0x65C70852876C50A0, 0x917389F0D88B919,  0x93FBE011EFA068E5, 0x72703759A4651388,
                                                     0xA266A6AEE1213EE0, 0xC496ABC50E388B6E, 0x564CB9FE0EBD72B5, 0x88B483157D3BADC8, 0x326D337A76B5317, 0x440F6F4F2D5A, // XPA1
                                                     0xCDD55D2646A1DE32, 0xAA056CDD8B80E53E, 0xAA87189B3A885C53, 0x9F6D9809057564A1, 0xC59794A13E1D38B8, 0x97F8ED39F3FA7DE5,
                                                     0xFC0CAF68C8B95129, 0x393F28B240A42FFE, 0xCD99B2F9792DEF96, 0xF1036825CBF416B9, 0x877B835F0533F2AE, 0xCDFFE238E18, // XQA0
                                                     0xAC4EF1B17010B136, 0xEC411E1B5AD8A667, 0x7737372EDB66A1AF, 0x43593ECED672CF87, 0x1E418547C7B8A975, 0x8CC78DCB18BD469A,
                                                     0x6C9FB93FD2EF8496, 0x8A4AC42666AB8545, 0x8A973B8387C15F1D, 0xC1657503D4BB4ADA, 0x22F49E4311D7BBF0, 0x1299B8FDA94C, // XQA1
                                                     0xC04B8957D3A4748F, 0xF3FB80F19063629F, 0x595434555D4EBE94, 0x8E1FEF11BFD1E0DA, 0xE31E3377248C0BB4, 0x9A05DEFF75EA51BA,
                                                     0x398686FBB343398A, 0x20331307B470DA54, 0x964FA62AD10005C5, 0x9EA5CC4D64E5D9EE, 0xC84675CF9B96060F, 0x1DECCB78CFAC, // XRA0
                                                     0x6B20FF684759DDC2, 0xD50EB91730DEAFBF, 0xAA5CA048E2DAF488, 0xE29708E28654FC18, 0x542928AD1F445359, 0xA311B83D79E73FF6,
                                                     0x850B7F5926826B22, 0x2D46731863BDB99D, 0x467A80CD8320B69D, 0xC046B12F05BFD513, 0x35D9B2FF794BDB40, 0x633276495B85 }; // XRA1
const uint64_t DBL_QA[2 * NWORDS64_FIELD]        = { 0xac2b2d74f883dfe4, 0xda9b5d82caa27d78, 0xf8656ebc40d57f4c, 0x5e1cd5bdbf041897, 0x1a30c6a718d110c8, 0x3c8def0dc70d6806,
                                                     0x91ab2c2c9282d88c, 0x8b6aedd25d129720, 0xaa92dd198282d20d, 0x976b9255cb297eb, 0xf6d8ab5c106ebe7a, 0xc5fb17b0515,
                                                     0xd5592babbcc2584e, 0xe0547c84cd5e0c1, 0xfbe528cb2d17b51b, 0x2168cc83a03036bc, 0x46149ea13591e9e3, 0xff230f71abe6a6e0,
                                                     0xd4a9a33beebb78fa, 0x63627d7cdda2d559, 0x601fcfc408949785, 0xcde4532f5618bcf6, 0xbc83162a741e1d9b, 0x6f443172fd95 };
const uint64_t YPA[2 * NWORDS64_FIELD]           = { 0xcf298a24ab4eecc8, 0x426be362d17b58e1, 0xdec9e4ab0c0813e4, 0xbb213d92b1a23dec, 0x73f9337ebf1afb1e, 0x22a428421e3f369d,
                                                     0x4c504fba8d4c1f41, 0x97c03c026b64e556, 0x94524150e5242247, 0x8f397d005b7d0f3, 0x2eeefb40c2d1d40f, 0x49ffa7880cc0,
                                                     0xc735ae6a6d9ab879, 0x4431fcf02097bf97, 0xcd6c8982e0b17062, 0xd38791b330dbf671, 0xbbe57bf59a8d9150, 0x58f968f33f45a7d8, 
                                                     0x903068b77ec83b26, 0x7dadadd772211f21, 0x2a8dea498d8a12bb, 0xae73b6ae3e7657f, 0x11695a4a18565fcd, 0x303fe0d52cdb};
const uint64_t PplusQA[4 * NWORDS64_FIELD]       = { 0xa7e27390cd89ed0, 0xf359eb3682d601b3, 0x7d893292e008d357, 0xca8297ede777799d, 0xaf494679324a5427, 0xc30a8772971b92b1,
                                                     0x291a6a8f56c222be, 0xd0d7f09ad8d323fa, 0xfa385cdcf0d0c4d5, 0x22c76490c77b6efb, 0x2466ad8bf91afd5d, 0x43e734531e6,
                                                     0xd733e958cb9c582b, 0x19b03b0647850342, 0x31be64c55b229ccd, 0xe25b77d735d49cce, 0x6555570cab893df8, 0xa452b814fe47d118,
                                                     0x8791843a4b61b101, 0xe761b1d5e99f89fe, 0x2b227b1d56a0931d, 0x6a279550bc1a989b, 0xdd9f4643c9ddd6df, 0x436657c7481e,
                                                     0x18718a5b58de448d, 0x44678e528b714548, 0x3bac89684e17847f, 0xd03d5e8a7a093d5a, 0xcf07e039c76e6f3e, 0x5ae4a7f32526fa36,
                                                     0x24d18348e9a45d10, 0x45d3164a37a0d0e6, 0xf221a442947e4bd7, 0x6bfbf5ae6db2c791, 0x1abb91ec57aeaac6, 0x5432433db9ad,
                                                     0x33293f2db350111f, 0xcedff7e53611ec93, 0x2d739b88b42a7c75, 0x4edb6b4121ae0dd0, 0x6b32ae397dd99f95, 0xd0ac8b36d0e24c89, 
                                                     0x23d6ae11a2b1c61a, 0x8a05380734cd9e89, 0xdcb14cd3c9f292f1, 0x3e24282abab56ebd, 0x69cc3fa3be707915, 0x5dafe89bb9f2 };

// Bob's generator values {XPB0, XQB0, XRB0 + XRB1*i} in GF(p751^2), expressed in Montgomery representation
const uint64_t B_gen[6 * NWORDS64_FIELD]         = { 0x110F4508C6634CCB, 0x31910BC05E296F4C, 0xED17AB0D6C029EA6, 0x9C863AB6172B9974, 0x5C15236CDB216F99, 0xDC025064818EC7D7,
                                                     0xC2180F387487EBF0, 0x946B1D0F025CBC3B, 0x5AE34395A520CB46, 0xB52034F98A879F2C, 0x3D2FAE10A22AB7C7, 0x174CD090DA3D, // XPB0
                                                     0xC3C6A839776171F2, 0x5883AFB529C8E50A, 0xDE1622BBD192925, 0x64CCE86B1826A21, 0x441AF1ABE9F6568E, 0x3F29EEC0BC6F962D,
                                                     0xA7845A0127159975, 0x109DCD6D92B0C3F2, 0x462438CD0100EE2E, 0xFB7869F2B1DF80EB, 0x563B0C55F0EEDC53, 0x1958C37D4721, // XPB1
                                                     0x110F4508C6634CCB, 0x31910BC05E296F4C, 0xED17AB0D6C029EA6, 0x9C863AB6172B9974, 0x5C15236CDB216F99, 0xDC025064818EC7D7,
                                                     0xC2180F387487EBF0, 0x946B1D0F025CBC3B, 0x5AE34395A520CB46, 0xB52034F98A879F2C, 0x3D2FAE10A22AB7C7, 0x174CD090DA3D, // XQB0
                                                     0x3C3957C6889E8E0D, 0xA77C504AD6371AF5, 0xF21E9DD442E6D6DA, 0xF9B331794E7D95DE, 0xBBE50E541609A971, 0xAF86113F439069D2,
                                                     0x3C683C8422E2DF33, 0xC9F7CDAC81470884, 0xC22A5F9AD5EAFA48, 0x89EA4B11AAD2D65C, 0xB7D78449A6CC0012, 0x568D11C4AFFA, // XQB1
                                                     0x31BB0964DFBDC34F, 0xFDC65CF4959AB106, 0xA3071E4B8B04D8FF, 0x9B68CFCE270DE486, 0x2339E590896E0095, 0xFC753508AD83E33E,
                                                     0x73A274E4A6908387, 0x88D1B207BBE8E2DC, 0xA6D0583233DC71F, 0xCF7F2ECC609DE5BE, 0xB8AF0669FBD1CF01, 0x1F3EF25DD512, // XRB0
                                                     0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0 }; // XRB1
/* Basis for Bob on A = 0, expressed in Montgomery representation */
const uint64_t P3[2 * NWORDS64_FIELD]            = { 0xF1A8C9ED7B96C4AB, 0x299429DA5178486E, 0xEF4926F20CD5C2F4, 0x683B2E2858B4716A, 0xDDA2FBCC3CAC3EEB, 0xEC055F9F3A600460,
                                                     0xD5A5A17A58C3848B, 0x4652D836F42EAED5, 0x2F2E71ED78B3A3B3, 0xA771C057180ADD1D, 0xC780A5D2D835F512, 0x114EA3B55AC1,
                                                     0x2E1EB8ED8C1C8C94, 0x6CFE456B25DBE01, 0x1EB54C3E8010F57A, 0x4B222D95FC81619D, 0xF99EBD204D501496, 0xC18348F9B629361,
                                                     0xC29E9A16BEDE6F96, 0x3B39F30163DAD41D, 0x807D3D1ECF2AC04E, 0xE088443F222A4988, 0x61B49A7524F1EA12, 0x41BF31133104};

// Montgomery constant Montgomery_R2 = (2^768)^2 mod p751
const uint64_t Montgomery_R2[NWORDS64_FIELD]     = { 0x233046449DAD4058, 0xDB010161A696452A, 0x5E36941472E3FD8E, 0xF40BFE2082A2E706, 0x4932CCA8904F8751 ,0x1F735F1F1EE7FC81, 
                                                     0xA24F4D80C1048E18, 0xB56C383CCDB607C5, 0x441DD47B735F9C90, 0x5673ED2C6A6AC82A, 0x06C905261132294B, 0x000041AD830F1F35 };                                                    
// Value one in Montgomery representation 
const uint64_t Montgomery_one[NWORDS64_FIELD]    = { 0x00000000000249ad, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x8310000000000000,
                                                     0x5527b1e4375c6c66, 0x697797bf3f4f24d0, 0xc89db7b2ac5c4e2e, 0x4ca4b439d2076956, 0x10f7926c7512c7e9, 0x00002d5b24bce5e2 };


// Fixed parameters for isogeny tree computation
const unsigned int strat_Alice[MAX_Alice-1] = { 
80, 48, 27, 15, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 7, 4, 2, 1, 1, 2, 1, 
1, 3, 2, 1, 1, 1, 1, 12, 7, 4, 2, 1, 1, 2, 1, 1, 3, 2, 1, 1, 1, 1, 5, 3, 2, 1, 1, 
1, 1, 2, 1, 1, 1, 21, 12, 7, 4, 2, 1, 1, 2, 1, 1, 3, 2, 1, 1, 1, 1, 5, 3, 2, 1, 
1, 1, 1, 2, 1, 1, 1, 9, 5, 3, 2, 1, 1, 1, 1, 2, 1, 1, 1, 4, 2, 1, 1, 1, 2, 1, 1, 
33, 20, 12, 7, 4, 2, 1, 1, 2, 1, 1, 3, 2, 1, 1, 1, 1, 5, 3, 2, 1, 1, 1, 1, 2, 1, 
1, 1, 8, 5, 3, 2, 1, 1, 1, 1, 2, 1, 1, 1, 4, 2, 1, 1, 2, 1, 1, 16, 8, 4, 2, 1, 1, 
1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1 };

const unsigned int strat_Bob[MAX_Bob-1] = { 
112, 63, 32, 16, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 8, 4, 2, 1, 1, 2, 1, 
1, 4, 2, 1, 1, 2, 1, 1, 16, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 8, 4, 2, 
1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 31, 16, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 
1, 1, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 15, 8, 4, 2, 1, 1, 2, 1, 1, 4, 
2, 1, 1, 2, 1, 1, 7, 4, 2, 1, 1, 2, 1, 1, 3, 2, 1, 1, 1, 1, 49, 31, 16, 8, 4, 2, 
1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 
15, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 7, 4, 2, 1, 1, 2, 1, 1, 3, 2, 1, 
1, 1, 1, 21, 12, 8, 4, 2, 1, 1, 2, 1, 1, 4, 2, 1, 1, 2, 1, 1, 5, 3, 2, 1, 1, 1, 1, 
2, 1, 1, 1, 9, 5, 3, 2, 1, 1, 1, 1, 2, 1, 1, 1, 4, 2, 1, 1, 1, 2, 1, 1 };

// Setting up macro defines and including GF(p), GF(p^2), curve, isogeny and kex functions
#define fpcopy                        fpcopy751
#define fpzero                        fpzero751
#define fpadd                         fpadd751
#define fpsub                         fpsub751
#define fpneg                         fpneg751
#define fpdiv2                        fpdiv2_751
#define fpcorrection                  fpcorrection751
#define fpmul_mont                    fpmul751_mont
#define fpsqr_mont                    fpsqr751_mont
#define fpinv_mont                    fpinv751_mont
#define fpinv_chain_mont              fpinv751_chain_mont
#define fpinv_mont_bingcd             fpinv751_mont_bingcd
#define fp2copy                       fp2copy751
#define fp2zero                       fp2zero751
#define fp2add                        fp2add751
#define fp2sub                        fp2sub751
#define mp_sub_p2                     mp_sub751_p2
#define mp_sub_p4                     mp_sub751_p4
#define sub_p4                        mp_sub_p4
#define fp2neg                        fp2neg751
#define fp2div2                       fp2div2_751
#define fp2correction                 fp2correction751
#define fp2mul_mont                   fp2mul751_mont
#define fp2sqr_mont                   fp2sqr751_mont
#define fp2inv_mont                   fp2inv751_mont
#define fp2inv_mont_bingcd            fp2inv751_mont_bingcd
#define fpequal_non_constant_time     fpequal751_non_constant_time
#define mp_add_asm                    mp_add751_asm
#define mp_subaddx2_asm               mp_subadd751x2_asm
#define mp_dblsubx2_asm               mp_dblsub751x2_asm
#define crypto_kem_keypair            crypto_kem_keypair_SIKEp751
#define crypto_kem_enc                crypto_kem_enc_SIKEp751
#define crypto_kem_dec                crypto_kem_dec_SIKEp751
#define random_mod_order_A            random_mod_order_A_SIDHp751
#define random_mod_order_B            random_mod_order_B_SIDHp751
#define EphemeralKeyGeneration_A      EphemeralKeyGeneration_A_SIDHp751
#define EphemeralKeyGeneration_B      EphemeralKeyGeneration_B_SIDHp751
#define EphemeralSecretAgreement_A    EphemeralSecretAgreement_A_SIDHp751
#define EphemeralSecretAgreement_B    EphemeralSecretAgreement_B_SIDHp751

#include "../fpx.c"
#include "../ec_isogeny.c"
#include "../sidh.c"
#include "../sike.c"