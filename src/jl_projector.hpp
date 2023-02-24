#ifndef EMP_ZK_JL_PROJECTOR_H
#define EMP_ZK_JL_PROJECTOR_H

#include "utils.hpp"

#include <utility>

#include "flint/fmpz.h"
#include "flint/flint.h"
#include "flint/fq.h"

#include <openssl/sha.h>

uint64_t random_oracle(int i, int j) {
	static unsigned char key[] = {0x01, 0x02, 0x03, 0x04, 0x05, 0x06, 0x07, 0x08, 0x09, 0x10, 0x11, 0x12, 0x13, 0x14, 0x15};

	SHA256_CTX ctx;
	SHA256_Init(&ctx);
	SHA256_Update(&ctx, key, 16);

	SHA256_Update(&ctx, (unsigned char*) &i, sizeof(int));
	SHA256_Update(&ctx, (unsigned char*) &j, sizeof(int));

	unsigned char result[32];
	SHA256_Final(result, &ctx);

	uint64_t res;
	memcpy(&res, result, sizeof(uint64_t));

	return mod(res);
}

typedef struct _jl_info {
	int k{};
	vector<uint64_t> keys_0;
	vector<uint64_t> keys_1;
	vector<uint64_t> keys_2;
	vector<uint64_t> keys_3;

	vector<uint64_t> limit;
	fmpz_t p;
} jl_info;

jl_info jl_projector_init(int k, vector<uint64_t> limit) {
	jl_info res;
	res.k = k;

	
	std::random_device rd;
    std::mt19937 gen(rd());

	std::uniform_int_distribution<unsigned long long> dis(0, 4611686018427322367ULL);

	
	/*
	for(int i = 0; i < k; i++) {
		res.keys_0.push_back(random_oracle(i, 0));
		res.keys_1.push_back(random_oracle(i, 1));
		res.keys_2.push_back(random_oracle(i, 2));
		res.keys_3.push_back(random_oracle(i, 3));
	}*/

	
	for(int i = 0; i < k; i++) {
		res.keys_0.push_back(dis(gen));
		res.keys_1.push_back(dis(gen));
		res.keys_2.push_back(dis(gen));
		res.keys_3.push_back(dis(gen));
	}

	res.limit = limit;
	fmpz_init_set_ui(res.p, p);

	return res;
}

void jl_projector_clean(jl_info &info) {
	fmpz_clear(info.p);
}

void jl_query_prf_with_prepared_data(int i, IntFp &zk_j, jl_info &info, uint64_t idx, bool** prepared_bits, uint64_t** prepared_witnesses, uint64_t &res, IntFp &zk_res, vector<IntFp> &zk_zero_checking) {
	auto a_start = clock_start();
	IntFp zk_j_cubic = zk_j * zk_j;
	zk_j_cubic = zk_j_cubic * zk_j;

	IntFp zk_val = IntFp(info.keys_0[i], PUBLIC);
	zk_val = zk_val + zk_j_cubic * info.keys_3[i];
	zk_val = zk_val + zk_j * info.keys_1[i];
	zk_val = zk_val + zk_j * zk_j * info.keys_2[i];

	const uint64_t QNR = 7;

	uint64_t inp_a = prepared_witnesses[i][idx];
	bool is_qr = prepared_bits[i][idx];
	IntFp zk_is_qr = IntFp(is_qr, ALICE);
	auto a_total = time_from(a_start);

	//printf("Aux 2A total: %0.8f\n", a_total); 

	auto b_start = clock_start();
	// perform bit testing of zk_is_qr
	{
		IntFp tmp;
		tmp = zk_is_qr * zk_is_qr;
		tmp = tmp + zk_is_qr.negate();

		zk_zero_checking.emplace_back(tmp);
	}

	IntFp zk_a = IntFp(inp_a, ALICE);

	// Compute f(T) * b + n * f(T) * (1 - b) mod p
	IntFp zk_left_sum = zk_val * zk_is_qr;
	IntFp zk_right_sum = IntFp(1, PUBLIC) + zk_is_qr.negate();
	zk_right_sum = zk_right_sum * QNR;
	zk_right_sum = zk_right_sum * zk_val;
	IntFp zk_check_sq_a = zk_left_sum + zk_right_sum;

	// Check that f(T) * b + n * f(T) * (1 - b) - a^2 mod p = 0
	IntFp zk_check_sq_a_tmp = zk_a * zk_a;
	zk_check_sq_a = zk_check_sq_a + zk_check_sq_a_tmp.negate();
	zk_zero_checking.emplace_back(zk_check_sq_a);

	auto b_total = time_from(b_start);

	//printf("Aux 2B total: %0.8f\n", b_total); 

	// return the QR test result
	res = is_qr;
	zk_res = zk_is_qr;
}

void jl_projector_prepare(uint64_t inp_1, uint64_t inp_2, uint64_t inp_3, uint64_t inp_4, jl_info &info, uint64_t idx, bool** prepared_bits, uint64_t** prepared_witnesses) {
	assert(inp_1 <= info.limit[0]);
	assert(inp_2 <= info.limit[1]);
	assert(inp_3 <= info.limit[2]);
	assert(inp_4 <= info.limit[3]);

	// compute j
	uint64_t j = 0;
	j = add_mod(j, mult_mod(inp_1, info.limit[1] * info.limit[2] * info.limit[3]));
	j = add_mod(j, mult_mod(inp_2, info.limit[2] * info.limit[3]));
	j = add_mod(j, mult_mod(inp_3, info.limit[3]));
	j = add_mod(j, mult_mod(inp_4, 1));

	const uint64_t QNR = 7;

	// compute if val is a quadratic residue
	fmpz_t fT;
	fmpz_t a;
	fmpz_init(fT);
	fmpz_init(a);

	for(int i = 0; i < info.k; i++) {
		uint64_t j_cubic = mult_mod(j, j);
		j_cubic = mult_mod(j_cubic, j);

		uint64_t val = info.keys_0[i];
		val += mult_mod(info.keys_1[i], j);
		val = mod(val);
		val += mult_mod(info.keys_2[i], mult_mod(j, j));
		val = mod(val);
		val += mult_mod(info.keys_3[i], j_cubic);
		val = mod(val);

		fmpz_set_ui(fT, val);
		uint64_t inp_a;
		// find the square root of fT (or the square root of QNR * fT)
		bool is_qr = fmpz_jacobi(fT, info.p) == 1;
		if(is_qr) {
			if (fmpz_sqrtmod(a, fT, info.p) == 1) {
				inp_a = fmpz_get_ui(a); // Input witness for the quadratic residue of f(T)*b + n * f(T) * (1-b) mod p
			} else {
				cerr << "Error with computing square root of f(T)" << endl;
				exit(1);
			}
		} else {
			fmpz_mul_ui(fT, fT, QNR);
			if (fmpz_sqrtmod(a, fT, info.p) == 1) {
				inp_a = fmpz_get_ui(a); // Input witness for the quadratic residue of f(T)*b + n * f(T) * (1-b) mod p
			} else {
				cerr << "Error with computing square root of n * f(T)" << endl;
				exit(1);
			}
		}

		prepared_bits[i][idx] = is_qr;
		prepared_witnesses[i][idx] = inp_a;
	}

	fmpz_clear(fT);
	fmpz_clear(a);
}

void jl_projector(IntFp zk_inp_1, IntFp zk_inp_2, IntFp zk_inp_3, IntFp zk_inp_4, jl_info &info, uint64_t idx, bool** prepared_bits, uint64_t** prepared_witnesses, vector<uint64_t> &res_projected, vector<IntFp> &zk_res_projected, vector<IntFp> &zk_zero_checking) {
	// compute j
	IntFp zk_j = IntFp(0, PUBLIC);
	zk_j = zk_j + zk_inp_1 * info.limit[1] * info.limit[2] * info.limit[3];
	zk_j = zk_j + zk_inp_2 * info.limit[2] * info.limit[3];
	zk_j = zk_j + zk_inp_3 * info.limit[3];
	zk_j = zk_j + zk_inp_4 * 1;

	vector<uint64_t> res;
	vector<IntFp> zk_res;
	for(int i = 0; i < info.k; i++) {
		res.push_back(0);
		zk_res.push_back(IntFp(0, PUBLIC));
	}

	for(int i = 0; i < info.k; i++) {
		jl_query_prf_with_prepared_data(i, zk_j, info, idx, prepared_bits, prepared_witnesses, res[i], zk_res[i], zk_zero_checking);
	}

	for(int i = 0; i < info.k; i++) {
		res_projected[i] = add_mod(res_projected[i], res[i]);
		zk_res_projected[i] = zk_res_projected[i] + zk_res[i];
	}
}

#endif //EMP_ZK_JL_PROJECTOR_H
