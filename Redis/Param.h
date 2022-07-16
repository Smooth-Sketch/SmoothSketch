#ifndef _PARAM_SMOOTH_H_
#define _PARAM_SMOOTH_H_
#include <bits/stdc++.h>
#include "Common/Matrix.h"

const uint32_t window_size = 10000;

// memory ratio of stage 1 to stage 2
const double stage_ratio = 0.2;

// the minimum value for a_K
const double var_thres = 1;

// number of consecutive windows
const int P = 7;

// degree of polynomial
const int K = 2;

// number of recorded windows in stage 1
const int S = 4;

// threshold for mean square error
const double error_thres = 4;

// number of cells in each bucket
const int bucket_size = 4;

const int KEY_LEN = 8;

const int _NumOfArray = 3;

// threshold for potential
const double potential_thres = 0.5;

double calcu_matrix[K + 1][P] = {};
double calcu_matrix_try[K + 1][S] = {};

void init_matrix() {
	// by normal equation, the LSE satisfies X'Xb = X'Y, hence b = (X'X)^{-1}X'Y
	// to get b for an arbitrary Y, we only need to find (X'X)^{-1}X'
	double X[P][K + 1] = {}, XX[K + 1][K + 1] = {}, temp[(K + 1) * (K + 1)];
	for (int i = 0; i <= P - 1; ++i) {
		for (int j = 0; j <= K; ++j) {
			X[i][j] = pow(i, j);
		}
	}
	for (int i = 0; i <= K; ++i) {
		for (int j = 0; j <= K; ++j) {
			for (int k = 0; k <= P - 1; ++k) {
				XX[i][j] += X[k][i] * X[k][j];
			}
		}
	}

	for (int i = 0; i <= K; ++i) {
		for (int j = 0; j <= K; ++j) {
			temp[i * (K + 1) + j] = XX[i][j];
		}
	}
	// find the inverse matrix of X'X
	inverse(K + 1, temp);

	for (int i = 0; i <= K; ++i) {
		for (int j = 0; j <= K; ++j) {
			XX[i][j] = temp[i * (K + 1) + j];
		}
	}

	// get (X'X)^{-1}X'
	for (int i = 0; i <= K; ++i) {
		for (int j = 0; j <= P - 1; ++j) {
			for (int k = 0; k <= K; ++k) {
				calcu_matrix[i][j] += XX[i][k] * X[j][k];
			}
		}
	}


	double Z[S][K + 1] = {};
	for (int i = 0; i <= S - 1; ++i) {
		for (int j = 0; j <= K; ++j) {
			X[i][j] = pow(i, j);
		}
	}
	for (int i = 0; i <= K; ++i) {
		for (int j = 0; j <= K; ++j) {
			for (int k = 0; k <= S - 1; ++k) {
				XX[i][j] += X[k][i] * X[k][j];
			}
		}
	}
	for (int i = 0; i <= K; ++i) {
		for (int j = 0; j <= K; ++j) {
			temp[i * (K + 1) + j] = XX[i][j];
		}
	}
	inverse(K + 1, temp);
	for (int i = 0; i <= K; ++i) {
		for (int j = 0; j <= K; ++j) {
			XX[i][j] = temp[i * (K + 1) + j];
		}
	}
	for (int i = 0; i <= K; ++i) {
		for (int j = 0; j <= K; ++j) {
			XX[i][j] = temp[i * (K + 1) + j];
		}
	}
	for (int i = 0; i <= K; ++i) {
		for (int j = 0; j <= S - 1; ++j) {
			for (int k = 0; k <= K; ++k) {
				calcu_matrix_try[i][j] += XX[i][k] * Z[j][k];
			}
		}
	}
}

void linear_regressing(double* y, double* b) {
	memset(b, 0, (K + 1) * sizeof(double));
	for (int i = 0; i <= K; ++i) {
		for (int j = 0; j <= P - 1; ++j) {
			b[i] += calcu_matrix[i][j] * y[j];
		}
	}
}

void calcu_variation(uint32_t c[], int d[], uint32_t size, uint32_t k) {
	int e[size];
	memcpy(e, c, size * sizeof(uint32_t));
	while (k--) {
		for (int i = 0; i < size - 1; ++i)
			e[i] = e[i + 1] - e[i];
		--size;
	}
	memcpy(d, e, size * sizeof(int));
}

void linear_regressing_try(double* y, double* b) {
	// only use S number of windows for linear regressing
	memset(b, 0, (K + 1) * sizeof(double));
	for (int i = 0; i <= K; ++i) {
		for (int j = 0; j <= S - 1; ++j) {
			b[i] += calcu_matrix[i][j] * y[j];
		}
	}
}

#endif
