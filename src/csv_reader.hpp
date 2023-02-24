#ifndef EMP_ZK_CSV_READER_H
#define EMP_ZK_CSV_READER_H

#include "utils.hpp"

// Reference: https://www.gormanalysis.com/blog/reading-and-writing-csv-files-with-cpp/

vector<vector<uint64_t>> read_csv(string file_name, int row_limit = INT32_MAX) {
	// Open CSV file
	ifstream csv_file;
	csv_file.open(file_name, ios::in);
	if (!csv_file || !csv_file.is_open() || !csv_file.good()) {
		cerr << "File " << file_name << " does not exist or is ill-formatted" << endl;
		exit(1);
	}

	// Simplify to value-only
	// due to challenges in accessing randomly with a vector of pairs
	vector<vector<uint64_t>> entries;

	string first_line;
	getline(csv_file, first_line);

	string buffer_line;
	int num_rows = 0;
	while (getline(csv_file, buffer_line) && num_rows < row_limit) {
		vector<uint64_t> line_numbers;

		stringstream ss(buffer_line);
		string cell;
		while (getline(ss, cell, ',')) {
			line_numbers.push_back((uint64_t) stoi(cell));
		}

		entries.push_back(line_numbers);
		num_rows++;
	}

	csv_file.close();
	return entries;
}

#endif //EMP_ZK_CSV_READER_H
