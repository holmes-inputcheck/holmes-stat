#ifndef EMP_ZK_BENCH_DATASET_1_H
#define EMP_ZK_BENCH_DATASET_1_H

#include "csv_reader.hpp"

typedef struct _financial_data_entry {
	uint64_t inp_age;
	IntFp zk_age;
	uint64_t inp_job;
	IntFp zk_job;
	uint64_t inp_marital;
	IntFp zk_marital;
	uint64_t inp_education;
	IntFp zk_education;
	uint64_t inp_default;
	IntFp zk_default;
	uint64_t inp_housing;
	IntFp zk_housing;
	uint64_t inp_loan;
	IntFp zk_loan;
	uint64_t inp_contact;
	IntFp zk_contact;
	uint64_t inp_month;
	IntFp zk_month;
	uint64_t inp_day_of_week;
	IntFp zk_day_of_week;
	uint64_t inp_duration;
	IntFp zk_duration;
	uint64_t inp_campaign;
	IntFp zk_campaign;
	uint64_t inp_pdays;
	IntFp zk_pdays;
	uint64_t inp_previous;
	IntFp zk_previous;
	uint64_t inp_poutcome;
	IntFp zk_poutcome;
	uint64_t inp_emp_var_rate;
	IntFp zk_emp_var_rate;
	uint64_t inp_cons_price_idx;
	IntFp zk_cons_price_idx;
	uint64_t inp_cons_conf_idx;
	IntFp zk_cons_conf_idx;
	uint64_t inp_euribor3m;
	IntFp zk_euribor3m;
	uint64_t inp_nr_employed;
	IntFp zk_nr_employed;
	uint64_t inp_y;
	IntFp zk_y;

	uint64_t derived_age_group;
	IntFp zk_derived_age_group;
} financial_data_entry;

vector<financial_data_entry> load_dataset(vector<uint64_t> &precomputed_dataset, int num_records) {
	assert(num_records <= 41188);

	vector<vector<uint64_t>> dataset_csv = read_csv("./src/dataset_1.csv", 41188);

	vector<financial_data_entry> dataset;
	dataset.reserve(num_records);

	
	for(int i = 0; i < dataset_csv.size(); i++) {
		vector<uint64_t> cur = dataset_csv[i];

		for(int j = 0; j < 21; j++) {
			precomputed_dataset.push_back(cur[j]);
		}

		financial_data_entry new_entry;
		new_entry.inp_age = cur[0] / 10 - 1;
		new_entry.inp_job = cur[1] - 1;
		new_entry.inp_marital = cur[2] - 1;
		new_entry.inp_education = cur[3] - 1;
		new_entry.inp_default = cur[4];
		new_entry.inp_housing = cur[5];
		new_entry.inp_loan = cur[6];
		new_entry.inp_contact = cur[7];
		new_entry.inp_month = cur[8];
		new_entry.inp_day_of_week = cur[9];
		new_entry.inp_duration = cur[10];
		new_entry.inp_campaign = cur[11];
		new_entry.inp_pdays = cur[12];
		new_entry.inp_previous = cur[13];
		new_entry.inp_poutcome = cur[14];
		new_entry.inp_emp_var_rate = cur[15];
		new_entry.inp_cons_price_idx = cur[16];
		new_entry.inp_cons_conf_idx = cur[17];
		new_entry.inp_euribor3m = cur[18];
		new_entry.inp_nr_employed = cur[19];
		new_entry.inp_y = cur[20];

		dataset.push_back(new_entry);
	}
	

	return dataset;
}

#endif //EMP_ZK_BENCH_DATASET_1_H