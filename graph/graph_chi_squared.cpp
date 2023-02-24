#include "../src/jl_projector.hpp"
#include "../src/bench_dataset_1.hpp"
#include <map>
#include <valarray>

const int threads = 32;

// number of dimensions
int num_dimensions = 1;

// the size of each dimension
int size_of_each_dimension = 35; // for simplicity, we assume each dimension has the same size

// the k parameter for JL, which is related to accuracy
// k is the same as r in the paper
int k = 200;

// number of entries to go through the range checks
// set this number to be high enough so offline phase batching pattern does not affect too much
int num_records = 40000;

typedef struct _data_entry {
	vector<uint64_t> inp;
	vector<int> sample;
} data_entry;

std::map<vector<int>, int> compute_sample_map(vector<data_entry> &dataset) {
	std::map<vector<int>, int> sample_map;
	for (int i = 0; i < num_records; i++) {
		vector<int> buckets;
		for (int j = 0; j < num_dimensions; j++) {
			buckets.emplace_back(dataset[i].sample[j]);
			sample_map[buckets]++;
		}
	}

	return sample_map;
}

vector<data_entry> create_fake_dataset(vector<vector<int>> indices, vector<double> expectation) {
	vector<data_entry> dataset;

	for (int i=0; i<indices.size(); i++) {
        vector<int> index = indices[i];
        int num_repetitions = round(expectation[i]);

        for (int n = 0; n < num_repetitions; n++) {
            data_entry new_entry;
            for (int j = 0; j < index.size(); j++) {
                int number = index[j];
                new_entry.sample.emplace_back(number);
                new_entry.inp.emplace_back(number);
            }
            dataset.push_back(new_entry);
        }        
    }

	return dataset;	
}

vector<double> permute_prob(vector<vector<double>> &prob_set, int end) {
	// permutes in this fashion:
	// does it in a lexicographic ordering fashion (of cartesian products)
	// 0 0 0 0 0
	// 0 0 0 0 1
	// 0 0 0 0 2
	// 0 0 0 1 0
	// 0 0 0 1 1 
	// 0 0 0 1 2 
	// 0 0 0 2 0
	// ...
	
	if (end == num_dimensions - 1) {
		vector<double> vec;
		for (int i = 0; i < size_of_each_dimension; i++)
			vec.push_back(prob_set[end][i] * num_records);
		return vec;
	}

	vector<double> vec = permute_prob(prob_set, end + 1);
	int removals = vec.size();

	for (int j = 0; j < size_of_each_dimension; j++) {
		for (int i = 0; i < removals; i++) {
			double mult = vec[i];
			vec.push_back(prob_set[end][j] * mult);
		}
	}

	for (int i = 0; i < removals; i++) {
		vec.erase(vec.begin());
	}

	return vec;
}

vector<vector<int>> permute_indices(int end) {
	// permutes in this fashion:
	// does it in a lexicographic ordering fashion (of cartesian products)
	// 0 0 0 0 0
	// 0 0 0 0 1
	// 0 0 0 0 2
	// 0 0 0 1 0
	// 0 0 0 1 1 
	// 0 0 0 1 2 
	// 0 0 0 2 0
	// ...
	
	if (end == num_dimensions - 1) {
		vector<vector<int>> vec;
		for (int i = 0; i < size_of_each_dimension; i++){
			vector<int> starting;
			starting.emplace_back(i);
			vec.emplace_back(starting);
		}
		return vec;
	}

	vector<vector<int>> vec = permute_indices(end + 1);
	int removals = vec.size();

	for (int j = 0; j < size_of_each_dimension; j++) {
		for (int i = 0; i < removals; i++) {
			vector<int> og = vec[i];
			vector<int> new_coordinate;
			new_coordinate.emplace_back(j);

			for (int m = 0; m < og.size(); m++) {
				new_coordinate.emplace_back(og[m]);
			}

			vec.push_back(new_coordinate);
		}
	}

	for (int i = 0; i < removals; i++) {
		vec.erase(vec.begin());
	}
		
	return vec;
}

vector<std::discrete_distribution<int>> create_distribution_set() {
	// create the distributions of each dimension
	vector<std::discrete_distribution<int>> dist_set;

	std::random_device rd{};
    std::mt19937 gen{rd()};

	// place the distributions in a vector set
	double mean = size_of_each_dimension / 2;
	std::normal_distribution<> cat_dist{mean, 5};
	std::uniform_int_distribution<> random_weight_dist(0, size_of_each_dimension - 1);
	
	std::vector<int> weights;

	for (int j = 0; j < size_of_each_dimension; j++) {
		weights.push_back(0);
	}

	for(int j=0; j < 50000; j++) {
		int weight = std::round(cat_dist(gen));
		if (weight < 0 || weight >= size_of_each_dimension)
			weight = random_weight_dist(gen);

		weights[weight] += 1;

	}
		
	// push in distributions
	for (int i = 0; i < num_dimensions; i++) {
		std::discrete_distribution<int> real_dist(weights.begin(), weights.end());
		dist_set.emplace_back(real_dist);
	}

	return dist_set;
}

uint64_t compute_idx_encoding(vector<int> data, jl_info info) {
	uint64_t j = 0;
	for (int beta = 0; beta < num_dimensions; beta ++) {
		uint64_t multiplier = 1;
		for (int nu = beta + 1; nu < num_dimensions; nu ++) {
			multiplier *= info.limit[nu];
		}
		j = add_mod(j, mult_mod(data[beta], multiplier));
	}

	return j;
}

void compute_legendre(size_t dataset_size, vector<uint64_t> &dataset_encoding, jl_info info, vector<double> &projection_vector) {
	#pragma omp parallel for default(shared)
	for(int i = 0; i < dataset_size; i++) {
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
	std::cout << std::endl << "------------ Chi-squared test simulated dataset graph ------------" << std::endl;

	{
		FILE * fp = fopen("./benchmark_input.txt", "r");
		fscanf(fp, "%d%d%d%d", &num_dimensions, &size_of_each_dimension, &k, &num_records);
		fclose(fp);
	}

	// our dimension size
	k = 200;
	vector<uint64_t> precomputed_dataset;
	vector<double> normalized_chi_sq;
    vector<double> unnormalized_chi_sq;
	vector<double> jl_chi_sq;

	/************************************************************************************/

	/* Create the expectation vector N*prob(dim_i), with lexicographically ordered indices */
	vector<std::discrete_distribution<int>> dist_set = create_distribution_set();
	vector<vector<double>> prob_set;
	for (auto dist : dist_set) {
		prob_set.emplace_back(dist.probabilities());
	}
	vector<double> dq = permute_prob(prob_set, 0); 

	/* Create the labeling of the lexicographically ordered indices */
	vector<vector<int>> indices = permute_indices(0); 

	/* Initialize the simulated dataset */
    auto dataset = create_fake_dataset(indices, dq);
    num_records = dataset.size();

	/* Initialize JL */
	vector<uint64_t> limit;
	{
		for(int i = 0; i < num_dimensions; i++) {
			limit.push_back(size_of_each_dimension);
		}
	}
	auto info = jl_projector_init(k, limit);
    const uint64_t QNR = 7;
	
	/* Initialize RNG */
	std::random_device rd;
    std::mt19937_64 gen(rd());
	std::uniform_int_distribution<uint64_t> index_dis(0, dataset.size() - 1);
	std::uniform_int_distribution<uint64_t> dim_dis(0, num_dimensions);

	/* Compute A * expected with JL + Legendre PRF */
	vector<uint64_t> expected_encoding;
	for (int i = 0; i < dq.size(); i++) {
		uint64_t rounded = int(round(dq[i]));
		
		// compute the one-hot encoding and store it for simplicity
		uint64_t j = compute_idx_encoding(indices[i], info);
		for (int l = 0; l < rounded; l++) {
			expected_encoding.emplace_back(j);
		}
    }
	vector<double> expec_proj_vector;
	for (int i = 0; i < info.k; i++) {
		expec_proj_vector.push_back(0);
	}
	compute_legendre(expected_encoding.size(), expected_encoding, info, expec_proj_vector);

	/* Initialize poisoning the dataset */
	double previous_chi_square = -1;
    int num_poisons = 12000; // Total number of poisons
	int e_iters = 1; // Granularity of each poison

	for (int e = 0; e < num_poisons; e++) {
		std::cout << e << endl;
		/* Poison the dataset */
		int epp = 0;
		while (epp < e_iters) {
			uint64_t rand_ind = index_dis(gen);
	
			for (int dim = 0; dim < num_dimensions; dim++) {
				uint64_t samp = (uint64_t) dataset[rand_ind].sample[dim];
				if (samp == info.limit[dim] - 1) {
					dataset[rand_ind].sample[dim] -= 1;
					dataset[rand_ind].inp[dim] -= 1;
				} else {
					dataset[rand_ind].sample[dim] += 1;
					dataset[rand_ind].inp[dim] += 1;
				}
				
				double naive_unnormalized_chi = 0.0;
				double naive_normalized_chi = 0.0;
				std::map<vector<int>, int> sample_map = compute_sample_map(dataset);
				for (int i = 0; i < indices.size(); i++) {
					double E_minus_O = double(sample_map[indices[i]]) - dq[i];
					naive_unnormalized_chi += E_minus_O * E_minus_O / num_records;
					naive_normalized_chi += E_minus_O * E_minus_O / dq[i];
				}
				
				if (naive_normalized_chi < previous_chi_square) {
					continue;
				} else {
					previous_chi_square = naive_normalized_chi;
					epp++;
				}
			}
		}

		/* Compute the normalized and unnormalized chi-square statistic */
		double naive_unnormalized_chi = 0.0;
		double naive_normalized_chi = 0.0;
		std::map<vector<int>, int> sample_map = compute_sample_map(dataset);
		for (int i = 0; i < indices.size(); i++) {
			double E_minus_O = double(sample_map[indices[i]]) - dq[i];
			naive_unnormalized_chi += E_minus_O * E_minus_O / num_records;
			naive_normalized_chi += E_minus_O * E_minus_O / dq[i];
		}
		unnormalized_chi_sq.emplace_back(naive_unnormalized_chi);
		normalized_chi_sq.emplace_back(naive_normalized_chi);

        /* Approximate A * count using JL + Legendre */
		vector<uint64_t> dataset_encoding;
		for (int i = 0; i < dataset.size(); i++) {
			uint64_t j = compute_idx_encoding(dataset[i].sample, info);
			dataset_encoding.emplace_back(j);
		}
        vector<double> hist_proj_vector;
        for(int i = 0; i < info.k; i++) {
            hist_proj_vector.push_back(0);
        }
		compute_legendre(dataset.size(), dataset_encoding, info, hist_proj_vector);

        /* Approximate A * (count - expected) */
        double jl_unnormalized_chi = 0;
        for (int i = 0; i < info.k; i++) {
            jl_unnormalized_chi += (hist_proj_vector[i] - expec_proj_vector[i]) * (hist_proj_vector[i] - expec_proj_vector[i]);
        }
        jl_chi_sq.emplace_back(jl_unnormalized_chi / dataset.size());
    }

	cout << "Finished all the JL projections" << endl;
	
	{
		FILE * fp = fopen("./benchmark_result.txt", "w");
		fprintf(fp, "Changes,actualchi,unchi,jlchi\n");
		for (int i = 0; i < num_poisons; i++) {
			 fprintf(fp, "%llu,%0.8f,%0.8f,%0.8f\n", i * e_iters, normalized_chi_sq[i], unnormalized_chi_sq[i], jl_chi_sq[i]);
		}

        FILE * fp_expec = fopen("./chi_expectation.txt", "w");
        fprintf(fp_expec, "Expectation\n");
		for (int i = 0; i < dq.size(); i++) {
			 fprintf(fp_expec, "%0.8f\n", dq[i]);
		}
		
		fclose(fp);
        fclose(fp_expec);
	}
	
	return 0;
}