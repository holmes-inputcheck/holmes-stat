#include "../src/utils.hpp"

const int threads = 32;

// number of entries to add to our dataset
// set this number to be high enough so offline phase batching pattern does not affect too much
int num_records;

// fixed point shift accuracy, default is 1000
int fixed_point_shift;
uint64_t population_mean = 100;
uint64_t population_sd = 10;

typedef struct _data_entry {
	uint64_t inp;
	double inp_dec;
} data_entry;

vector<data_entry> create_random_dataset(vector<uint64_t> &precomputed_dataset) {
	// randomly generate a dataset where all elements are all in the right ranges

	std::default_random_engine data_generator; // set up gaussian distribution
    std::normal_distribution<double> data_distribution(population_mean, population_sd); // mean centered at 10, sd = 1
	
	vector<data_entry> dataset;
	for(int i = 0; i < num_records; i++) {
		data_entry new_entry;
		double random_element = data_distribution(data_generator);
		new_entry.inp_dec = random_element;

		int64_t raw_entry = (int64_t) (round(random_element * fixed_point_shift));

		if (raw_entry < 0)
			new_entry.inp = mult_mod(p - 1, -raw_entry);
		else
			new_entry.inp = raw_entry;
		
		/*
		if (new_entry.inp_dec < 0) 
			new_entry.zk_inp = IntFp(-raw_entry, ALICE).negate();
		else
			new_entry.zk_inp = IntFp(raw_entry);*/

		precomputed_dataset.push_back(new_entry.inp);
		dataset.push_back(new_entry);
	}

	return dataset;
}

void compute_sum(vector<data_entry> &dataset, uint64_t &inp_sum, double &inp_dec_sum) {
    for(int i = 0; i < dataset.size(); i++) {
		inp_sum = add_mod(inp_sum, dataset[i].inp);
		inp_dec_sum += dataset[i].inp_dec;
	}
}

double compute_sample_mean(vector<data_entry> &dataset, vector<uint64_t> &indices_changed) {
    uint64_t sum = 0;
    for(int i = 0; i < indices_changed.size(); i++) {
		sum += dataset[indices_changed[i]].inp;
	}
    return double(sum) / indices_changed.size();
}

double compute_variance(vector<data_entry> &dataset, double pop_mean) {
    double var = 0;
    for(int i = 0; i < dataset.size(); i++ ) {
        var += (dataset[i].inp_dec - pop_mean) * (dataset[i].inp_dec - pop_mean) / (num_records - 1);
    }
    return var;
}

double compute_squared_sum(vector<data_entry> &dataset) {
	double squared_sum = 0;
	for (int i = 0; i < dataset.size(); i++) {
		squared_sum += dataset[i].inp_dec * dataset[i].inp_dec;
	}
	return squared_sum;
}

int main(int argc, char** argv) {
	std::cout << std::endl << "------------ T-test simulated dataset graph ------------" << std::endl;

	{
		FILE * fp = fopen("./benchmark_input.txt", "r");
		fscanf(fp, "%d%d", &num_records, &fixed_point_shift);
		fclose(fp);
	}

    std::random_device rd;  //Will be used to obtain a seed for the random number engine
    std::mt19937 index_gen(rd());
	std::uniform_int_distribution<> index_distribution(0, num_records);

	vector<uint64_t> precomputed_dataset;

	/************************************************************************************/

	auto dataset = create_random_dataset(precomputed_dataset);

	// Initialize the sum to be PUBLIC because both parties must agree on the same input sum
    uint64_t total_changes = 12000; // to change the mean by 1: x * fixed_point_shift / num_records = 1;
    uint64_t number_changed[total_changes];
    uint64_t output_mean[total_changes];
	
    uint64_t output_t_test[total_changes]; // outputted t statistic according to student's t-test
    double output_t_dec_test[total_changes]; //output t statistical with double precision

    for (int i = 0; i < total_changes; i++) {
        uint64_t rand_index = index_distribution(index_gen);
        dataset[rand_index].inp = add_mod(dataset[rand_index].inp, fixed_point_shift / 10); // deviate an entry by 0.1 (fixed_point_shift / 10)
		dataset[rand_index].inp_dec += 1;

        uint64_t inp_sum = 0;
		double inp_dec_sum = 0.0;
        compute_sum(dataset, inp_sum, inp_dec_sum);

        output_mean[i] = inp_sum / num_records; //add the mean changed;
		double output_dec_mean = inp_dec_sum / num_records;

		double sq_sum = compute_squared_sum(dataset);

		// compute the sample variance
		double variance = (sq_sum * num_records - inp_dec_sum * inp_dec_sum) / num_records / (num_records - 1);
		double test_variance = compute_variance(dataset, output_dec_mean);
		//cout << "variance : " << variance << " , test variance: " << test_variance << endl;

		
        // compute t test statistic
        number_changed[i] = i + 1; // add the number changed;
		uint64_t scaled_pop_mean = population_mean * fixed_point_shift;

		uint64_t scaled_pop_sd = population_sd * fixed_point_shift;
		uint64_t scaled_sqroot_n = (int64_t) (sqrt(num_records) * fixed_point_shift);
		
		//cout << output_mean[i] - scaled_pop_mean << " " << output_dec_mean - population_mean << endl;

		
		// TO-DO: do the proper negation squaring for computing the t-test statistic
        output_t_test[i] = mult_mod(add_mod(output_mean[i], mult_mod(p - 1, scaled_pop_mean)), (scaled_sqroot_n / scaled_pop_sd)); // multiply since we don't want divide by 0
		//output_t_dec_test[i] =  (output_dec_mean - population_mean) / (sqrt(variance) / sqrt(num_records));
		output_t_dec_test[i] =  (output_dec_mean - population_mean) / (population_sd / sqrt(num_records));  
		//cout << "output t-test: " << output_t_test[i] << ", output t-test precision: " << output_t_dec_test[i] << endl;
		
	}

	{
		FILE * fp = fopen("./benchmark_result.txt", "w");
		fprintf(fp, "Changes,t,precisiont\n");
        for (int i = 0; i < total_changes; i++) {
			printf("%llu,%llu,%0.8f\n", number_changed[i], output_t_test[i], output_t_dec_test[i]);
            fprintf(fp, "%llu,%llu,%0.8f\n", number_changed[i], output_t_test[i], output_t_dec_test[i]);
        }
		fclose(fp);
	}

	return 0;
}
