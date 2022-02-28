# FKG_SIDH
Supersingular isogeny Diffie-Hellman (SIDH) is attractive for its small public key size, but it is still unsatisfactory due to its efficiency, compared to other post-quantum proposals. In this paper, we focus on the performance of SIDH when the starting curve is 𝐸_6 : 𝑦^2 = 𝑥^3 + 6𝑥^2 + 𝑥, which is fixed in Round-3 SIKE implementation. Inspired by the previous work, we present several tricks to accelerate key generation of SIDH and each process of SIKE. Our experimental results show that the
performance of this work is at least 6.09% faster than that of the current SIKE implementation, and we can further improve the performance when large storage is available.

Our code is based on SIDH v3.4 (https://github.com/Microsoft/PQCrypto-SIDH), and the usage of our code is as same as the latter.

See README.md in the file PQCrypto-SIDH_Mladder (or PQCrypto-SIDH_3PTladder_precomp) for details.
