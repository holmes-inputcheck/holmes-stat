#include "../src/jl_projector.hpp"
#include <map>
#include <valarray>
#include "../src/bench_dataset_1.hpp"
#include <algorithm>
#include <random>


const int threads = 32;

// the k parameter for JL, which is related to accuracy
 // k is the same as r in the paper
int k = 200;

// number of entries to go through the range checks
// set this number to be high enough so offline phase batching pattern does not affect too much
int num_records = 41118;

typedef struct _data_entry {
	vector<uint64_t> inp;
	vector<IntFp> zk_inp;
	vector<int> sample;
} data_entry;


void count_dataset(vector<financial_data_entry> &dataset, vector<uint64_t> &counter_dataset, vector<uint64_t> &dataset_encoding, vector<uint64_t> limit) {
	for (int i = 0; i < dataset.size(); i++) {
		uint64_t j = 0;
		j = add_mod(j, dataset[i].inp_education);
		j = mult_mod(j, limit[3]);
		j = add_mod(j, dataset[i].inp_marital);
		j = mult_mod(j, limit[2]);
		j = add_mod(j, dataset[i].inp_job);
		j = mult_mod(j, limit[1]);
		j = add_mod(j, dataset[i].inp_age);

		dataset_encoding.emplace_back(j);

		counter_dataset[j] += 1;
	}
}

void compute_legendre(vector<financial_data_entry> &dataset, vector<uint64_t> &dataset_encoding, jl_info info, vector<double> &projection_vector) {
	for(int i = 0; i < dataset.size(); i++) {
		fmpz_t fT;
		fmpz_init(fT);
		
		uint64_t j = dataset_encoding[i];

		for(int l = 0; l < k; l++) {
			uint64_t j_cubic = mult_mod(j, j);
			j_cubic = mult_mod(j_cubic, j);

			uint64_t val = info.keys_0[l];
			val += mult_mod(info.keys_1[l], j);
			val = mod(val);
			val += mult_mod(info.keys_2[l], mult_mod(j, j));
			val = mod(val);
			val += mult_mod(info.keys_3[l], j_cubic);
			val = mod(val);

			fmpz_set_ui(fT, val);

			// find the square root of fT (or the square root of QNR * fT)
			bool is_qr = fmpz_jacobi(fT, info.p) == 1;

			if (is_qr) {
				projection_vector[l] += 1 / sqrt(info.k);
			} else {
				projection_vector[l] -= 1 / sqrt(info.k);
			}
		}

		fmpz_clear(fT);
	}
}

int main(int argc, char** argv) {
	std::cout << std::endl << "------------ Chi-squared test real dataset graph ------------" << std::endl;

	{
		FILE * fp = fopen("./benchmark_input.txt", "r");
		fscanf(fp, "%d%d", &k, &num_records);
		fclose(fp);
	}

	vector<uint64_t> precomputed_dataset;
	vector<double> normalized_chi_sq;
	vector<double> unnormalized_chi_sq;
	vector<double> jl_chi_sq;


	/************************************************************************************/

	/* Initialize training and public dataset */
	int num_records = 41188;
	auto public_dataset = load_dataset(precomputed_dataset, num_records);
	const uint64_t QNR = 7;

	auto rng = std::linear_congruential_engine<uint_fast32_t, 48271, 0, 2147483647> {};
	std::shuffle(std::begin(public_dataset), std::end(public_dataset), rng);

	vector<financial_data_entry> train_dataset;
	for (int i = 0; i < num_records / 2; i++) {
		train_dataset.emplace_back(public_dataset[i]);
	}

	/* Initialize JL */
	vector<uint64_t> limit = {9, 12, 4, 8};
	auto info = jl_projector_init(k, limit);
	int D = 3456;
	
	/* Count the public dataset */
	vector<uint64_t> counter_public_dataset;
	for (int i = 0; i < D; i++) {
		counter_public_dataset.emplace_back(0);
	}
	vector<uint64_t> public_dataset_encoding;
	count_dataset(public_dataset, counter_public_dataset, public_dataset_encoding, limit);
	
	/* Approximate A * expected with JL + Legendre PRF */
	vector<double> expec_proj_vector;
	for(int i = 0; i < k; i++) {
		expec_proj_vector.push_back(0);
	}
	compute_legendre(public_dataset, public_dataset_encoding, info, expec_proj_vector);
	for (int i = 0; i < info.k; i++) {
		expec_proj_vector[i] = (double(train_dataset.size()) / public_dataset.size()) * double(expec_proj_vector[i]);
	}

	/* Initialize RNG */
	std::random_device rd;
    std::mt19937_64 gen(rd());
	std::uniform_int_distribution<uint64_t> index_dis(0, train_dataset.size() - 1);
	std::uniform_int_distribution<uint64_t> dim_dis(0, 3);
	
	/* Initialize poisoning the dataset */
	int num_poisons = 12000; // Total number of poisons
	int e_iters = 1; // Granularity of each poison
	
	for (int e = 0; e < num_poisons; e++) {
		std::cout << e << endl;
		/* Poison the dataset */
		int epp = 0;
		while(epp < e_iters) {
			uint64_t rand_ind = index_dis(gen);
			uint64_t rand_dim = dim_dis(gen);
			uint64_t inp;

			switch (rand_dim) {
				case 0:
					inp = train_dataset[rand_ind].inp_age;
					if (inp == info.limit[rand_dim] - 1) {
						train_dataset[rand_ind].inp_age -= 1;
					} else {
						train_dataset[rand_ind].inp_age += 1;
					}
				case 1:
					inp = train_dataset[rand_ind].inp_job;
					if (inp == info.limit[rand_dim] - 1) {
						train_dataset[rand_ind].inp_job -= 1;
					} else {
						train_dataset[rand_ind].inp_job += 1;
					}
				case 2:
					inp = train_dataset[rand_ind].inp_marital;
					if (inp == info.limit[rand_dim] - 1) {
						train_dataset[rand_ind].inp_marital -= 1;
					} else {
						train_dataset[rand_ind].inp_marital += 1;
					}
				case 3:
					inp = train_dataset[rand_ind].inp_education;
					if (inp == info.limit[rand_dim] - 1) {
						train_dataset[rand_ind].inp_education -= 1;
					} else {
						train_dataset[rand_ind].inp_education += 1;
					}
			}

			epp++;
		}
		
		/* Count the poisoned training dataset */
		vector<uint64_t> counter_train_dataset;
		for (int i = 0; i < D; i++) {
			counter_train_dataset.emplace_back(0);
		}
		vector<uint64_t> train_dataset_encoding;
		count_dataset(train_dataset, counter_train_dataset, train_dataset_encoding, limit);
		
		/* Compute the normalized and unnormalized chi-square statistic */
		double naive_unnormalized_chi = 0.0;
		double naive_normalized_chi = 0.0;
		int nonzero_D = 0;
		for (int i = 0; i < D; i++) {
			double O_minus_E = (double(counter_train_dataset[i]) - double(counter_public_dataset[i]) * double(train_dataset.size()) / public_dataset.size());
			naive_unnormalized_chi += O_minus_E * O_minus_E / train_dataset.size();

			if (counter_public_dataset[i] != 0) {
				naive_normalized_chi += O_minus_E * O_minus_E / (double(counter_public_dataset[i]) * double(train_dataset.size()) / public_dataset.size());
				nonzero_D++;
			}
		}
		normalized_chi_sq.emplace_back(naive_normalized_chi);
		unnormalized_chi_sq.emplace_back(naive_unnormalized_chi);

		/* Approximate A * count with JL + Legendre PRF */		
		vector<double> hist_proj_vector;
		for(int i = 0; i < k; i++) {
			hist_proj_vector.push_back(0);
		}
		compute_legendre(train_dataset, train_dataset_encoding, info, hist_proj_vector);

		/* Approximate A * (count - expected) */
		double jl_unnormalized_chi = 0.0;
		for (int i = 0; i < k; i++) {
			jl_unnormalized_chi += (hist_proj_vector[i] - expec_proj_vector[i]) * (hist_proj_vector[i] - expec_proj_vector[i]);
		}
		jl_chi_sq.emplace_back(jl_unnormalized_chi / train_dataset.size());
	}

	cout << "Finished all the JL projections" << endl;

	{
		FILE * fp = fopen("./benchmark_result.txt", "w");
		FILE * fp_expec = fopen("./chi_expectation.txt", "w");
		fprintf(fp, "Changes,actualchi,unchi,jlchi\n");
		for (int i = 0; i < num_poisons; i++) {
			 fprintf(fp, "%llu,%0.8f,%0.8f,%0.8f\n", i * e_iters, normalized_chi_sq[i], unnormalized_chi_sq[i], jl_chi_sq[i]);
		}

		fprintf(fp_expec, "Expectation\n");
		for (int i = 0; i < D; i++) {
			 fprintf(fp_expec, "%llu\n", counter_public_dataset[i]);
		}
		
		fclose(fp);
		fclose(fp_expec);
	}
	
	return 0;
}
