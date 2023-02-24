#include "../src/bench_dataset_1.hpp"

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

void compute_sum(vector<financial_data_entry> &dataset, uint64_t &inp_sum) {
    for(int i = 0; i < dataset.size(); i++) {
		inp_sum = add_mod(inp_sum, dataset[i].inp_duration);
	}
}

double compute_sample_mean(vector<financial_data_entry> &dataset, vector<uint64_t> &indices_changed) {
    uint64_t sum = 0;
    for(int i = 0; i < indices_changed.size(); i++) {
		sum += dataset[indices_changed[i]].inp_duration;
	}
    return double(sum) / indices_changed.size();
}

double compute_variance(vector<financial_data_entry> &dataset, double pop_mean, int num_records) {
    double var = 0;
    for(int i = 0; i < dataset.size(); i++ ) {
        var += (dataset[i].inp_duration - pop_mean) * (dataset[i].inp_duration - pop_mean) / (num_records - 1);
    }
    return var;
}

double compute_squared_sum(vector<financial_data_entry> &dataset) {
	double squared_sum = 0;
	for (int i = 0; i < dataset.size(); i++) {
		squared_sum += double(dataset[i].inp_duration) * double(dataset[i].inp_duration);
	}
	return squared_sum;
}

int main(int argc, char** argv) {
	std::cout << std::endl << "------------ T-test real dataset graph ------------" << std::endl;

	{
		FILE * fp = fopen("./benchmark_input.txt", "r");
		fscanf(fp, "%d%d", &num_records, &fixed_point_shift);
		fclose(fp);
	}


    std::random_device rd;  //Will be used to obtain a seed for the random number engine
    std::mt19937 index_gen(rd());
	std::uniform_int_distribution<> index_distribution(0, num_records);

	auto total_time_start = clock_start();

	vector<uint64_t> precomputed_dataset;

	/************************************************************************************/

	//auto dataset = create_random_dataset(precomputed_dataset);

    auto public_dataset = load_dataset(precomputed_dataset, num_records);
	auto rng = std::default_random_engine {};
	std::shuffle(std::begin(public_dataset), std::end(public_dataset), rng);

    int public_num_records = public_dataset.size();

	vector<financial_data_entry> first_dataset;
	for (int i = 0; i < public_num_records / 2; i++) {
		first_dataset.emplace_back(public_dataset[i]);
	}

    num_records = first_dataset.size();

    uint64_t public_inp_sum = 0;
    compute_sum(public_dataset, public_inp_sum);
    double public_sum = (double) public_inp_sum;
    double public_mean = public_sum / public_num_records;
    //cout << std::setprecision (12) << "public mean: " << public_mean << endl;
    double public_sq_sum = compute_squared_sum(public_dataset);
    double public_variance = (public_sq_sum * public_num_records - public_sum * public_sum) / public_num_records / (public_num_records - 1);
    double public_test_variance = compute_variance(public_dataset, public_mean, public_num_records);
    //cout << std::setprecision (12) << "public variance : " << public_variance << " , public test variance: " << public_test_variance << endl;
	
    //cout << std::setprecision (12) << "public size: " << public_num_records << ", first size: " << num_records << endl;
    // Initialize the sum to be PUBLIC because both parties must agree on the same input sum
    uint64_t total_changes = 12000; // to change the mean by 1: x * fixed_point_shift / num_records = 1;
    uint64_t number_changed[total_changes];
    uint64_t output_mean[total_changes];
	
    uint64_t output_t_test[total_changes]; // outputted t statistic according to student's t-test
    double original_degree_of_freedom;
    uint64_t degree_of_freedom[total_changes];
    double output_t_dec_test[total_changes]; //output t statistical with double precision

    for (int i = 0; i < total_changes; i++) {
        uint64_t inp_sum = 0;
        compute_sum(first_dataset, inp_sum);
        double first_sum = (double) inp_sum;

		double first_mean = first_sum / num_records;
        //cout << std::setprecision (12) << "first dataset mean: " << first_mean << endl;
		double first_sq_sum = compute_squared_sum(first_dataset);

		// compute the sample variance
		double first_variance = (first_sq_sum * num_records - first_sum * first_sum) / num_records / (num_records - 1);
		double first_test_variance = compute_variance(first_dataset, first_mean, num_records);
		//cout << std::setprecision (12) << "first dataset variance : " << first_variance << " , first test variance: " << first_test_variance << endl;

        double t1 = first_variance / num_records;
        double t2 = public_variance / public_num_records;
       
		
        // compute t test statistic
        number_changed[i] = i; // add the number changed;
		
        double t_stat = (first_mean - public_mean) / sqrt(t1 + t2);
        output_t_dec_test[i] =  t_stat;

        
        double num = (t1 + t2) * (t1 + t2);
        double denom = (t1 * t1) / (num_records - 1) + (t2 * t2) / (public_num_records - 1);

        if (i == 0) {
            original_degree_of_freedom = num / denom;
        }

        degree_of_freedom[i] = (uint64_t) round(num/denom);

		//cout << std::setprecision (12) << "Output t-test precision: " << output_t_dec_test[i] << endl;
        //cout << std::setprecision (12) << "Degree of freedom: " << degree_of_freedom[i] << endl;

        // pertube our first dataset
        uint64_t rand_index = index_distribution(index_gen);
        first_dataset[rand_index].inp_duration = add_mod(first_dataset[rand_index].inp_duration, 25);
		
	}
    
	{
		FILE * fp = fopen("./benchmark_result.txt", "w");
		fprintf(fp, "Changes,df,precisiont\n");
        for (int i = 0; i < total_changes; i++) {
			//printf("%llu,%0.8f,%0.8f\n", number_changed[i], degree_of_freedom[i], output_t_dec_test[i]);
            fprintf(fp, "%llu,%llu,%0.8f\n", number_changed[i], degree_of_freedom[i], output_t_dec_test[i]);
        }
		fclose(fp);
	}

	return 0;
}
