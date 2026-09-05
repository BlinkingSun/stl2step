// Public-domain SHA-256 (D-140-1(6) featureId). No new dependency.
#ifndef GRADE_SHA256_HPP
#define GRADE_SHA256_HPP

#include <cstdint>
#include <cstdio>
#include <cstring>
#include <string>
#include <vector>

namespace grade {

class Sha256 {
public:
    Sha256() { reset(); }

    void reset() {
        len_ = 0;
        h_[0] = 0x6a09e667u;
        h_[1] = 0xbb67ae85u;
        h_[2] = 0x3c6ef372u;
        h_[3] = 0xa54ff53au;
        h_[4] = 0x510e527fu;
        h_[5] = 0x9b05688cu;
        h_[6] = 0x1f83d9abu;
        h_[7] = 0x5be0cd19u;
        nblock_ = 0;
    }

    void update(const void* data, size_t n) {
        const auto* p = static_cast<const uint8_t*>(data);
        while (n > 0) {
            const size_t room = 64 - nblock_;
            const size_t take = n < room ? n : room;
            std::memcpy(block_ + nblock_, p, take);
            nblock_ += take;
            p += take;
            n -= take;
            if (nblock_ == 64) {
                compress();
                nblock_ = 0;
            }
            len_ += take;
        }
    }

    void update(const std::string& s) { update(s.data(), s.size()); }

    std::string hex() {
        uint8_t out[32];
        finish(out);
        static const char* kHex = "0123456789abcdef";
        std::string r(64, '0');
        for (int i = 0; i < 32; ++i) {
            r[2 * i] = kHex[out[i] >> 4];
            r[2 * i + 1] = kHex[out[i] & 15];
        }
        return r;
    }

    std::string hex12() { return hex().substr(0, 12); }

private:
    static uint32_t rotr(uint32_t x, int n) { return (x >> n) | (x << (32 - n)); }

    void compress() {
        static const uint32_t K[64] = {
            0x428a2f98u, 0x71374491u, 0xb5c0fbcfu, 0xe9b5dba5u, 0x3956c25bu, 0x59f111f1u,
            0x923f82a4u, 0xab1c5ed5u, 0xd807aa98u, 0x12835b01u, 0x243185beu, 0x550c7dc3u,
            0x72be5d74u, 0x80deb1feu, 0x9bdc06a7u, 0xc19bf174u, 0xe49b69c1u, 0xefbe4786u,
            0x0fc19dc6u, 0x240ca1ccu, 0x2de92c6fu, 0x4a7484aau, 0x5cb0a9dcu, 0x76f988dau,
            0x983e5152u, 0xa831c66du, 0xb00327c8u, 0xbf597fc7u, 0xc6e00bf3u, 0xd5a79147u,
            0x06ca6351u, 0x14292967u, 0x27b70a85u, 0x2e1b2138u, 0x4d2c6dfcu, 0x53380d13u,
            0x650a7354u, 0x766a0abbu, 0x81c2c92eu, 0x92722c85u, 0xa2bfe8a1u, 0xa81a664bu,
            0xc24b8b70u, 0xc76c51a3u, 0xd192e819u, 0xd6990624u, 0xf40e3585u, 0x106aa070u,
            0x19a4c116u, 0x1e376c08u, 0x2748774cu, 0x34b0bcb5u, 0x391c0cb3u, 0x4ed8aa4au,
            0x5b9cca4fu, 0x682e6ff3u, 0x748f82eeu, 0x78a5636fu, 0x84c87814u, 0x8cc70208u,
            0x90befffau, 0xa4506cebu, 0xbef9a3f7u, 0xc67178f2u};
        uint32_t w[64];
        for (int i = 0; i < 16; ++i) {
            w[i] = (uint32_t(block_[4 * i]) << 24) | (uint32_t(block_[4 * i + 1]) << 16) |
                   (uint32_t(block_[4 * i + 2]) << 8) | uint32_t(block_[4 * i + 3]);
        }
        for (int i = 16; i < 64; ++i) {
            const uint32_t s0 = rotr(w[i - 15], 7) ^ rotr(w[i - 15], 18) ^ (w[i - 15] >> 3);
            const uint32_t s1 = rotr(w[i - 2], 17) ^ rotr(w[i - 2], 19) ^ (w[i - 2] >> 10);
            w[i] = w[i - 16] + s0 + w[i - 7] + s1;
        }
        uint32_t a = h_[0], b = h_[1], c = h_[2], d = h_[3];
        uint32_t e = h_[4], f = h_[5], g = h_[6], h = h_[7];
        for (int i = 0; i < 64; ++i) {
            const uint32_t S1 = rotr(e, 6) ^ rotr(e, 11) ^ rotr(e, 25);
            const uint32_t ch = (e & f) ^ ((~e) & g);
            const uint32_t t1 = h + S1 + ch + K[i] + w[i];
            const uint32_t S0 = rotr(a, 2) ^ rotr(a, 13) ^ rotr(a, 22);
            const uint32_t maj = (a & b) ^ (a & c) ^ (b & c);
            const uint32_t t2 = S0 + maj;
            h = g;
            g = f;
            f = e;
            e = d + t1;
            d = c;
            c = b;
            b = a;
            a = t1 + t2;
        }
        h_[0] += a;
        h_[1] += b;
        h_[2] += c;
        h_[3] += d;
        h_[4] += e;
        h_[5] += f;
        h_[6] += g;
        h_[7] += h;
    }

    void finish(uint8_t out[32]) {
        const uint64_t bitlen = len_ * 8;
        uint8_t pad = 0x80;
        update(&pad, 1);
        pad = 0;
        while (nblock_ != 56) update(&pad, 1);
        uint8_t lenb[8];
        for (int i = 0; i < 8; ++i) lenb[7 - i] = static_cast<uint8_t>((bitlen >> (8 * i)) & 0xffu);
        update(lenb, 8);
        for (int i = 0; i < 8; ++i) {
            out[4 * i] = static_cast<uint8_t>(h_[i] >> 24);
            out[4 * i + 1] = static_cast<uint8_t>(h_[i] >> 16);
            out[4 * i + 2] = static_cast<uint8_t>(h_[i] >> 8);
            out[4 * i + 3] = static_cast<uint8_t>(h_[i]);
        }
    }

    uint32_t h_[8];
    uint8_t block_[64];
    size_t nblock_;
    uint64_t len_;
};

inline std::string sha256File(const std::string& path) {
    std::FILE* f = std::fopen(path.c_str(), "rb");
    if (!f) return std::string(64, '0');
    Sha256 s;
    uint8_t buf[4096];
    size_t n;
    while ((n = std::fread(buf, 1, sizeof buf, f)) > 0) s.update(buf, n);
    std::fclose(f);
    return s.hex();
}

inline uint64_t fileSize(const std::string& path) {
    std::FILE* f = std::fopen(path.c_str(), "rb");
    if (!f) return 0;
    if (std::fseek(f, 0, SEEK_END) != 0) {
        std::fclose(f);
        return 0;
    }
    const long sz = std::ftell(f);
    std::fclose(f);
    return sz < 0 ? 0 : static_cast<uint64_t>(sz);
}

}  // namespace grade

#endif
