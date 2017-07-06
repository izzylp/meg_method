/* Compile and run: gcc -W -Wall -o meg meg.c
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>

typedef enum {false, true} bool;

struct linear_system
{
	int32_t n;
	double *matrix;
	double *b;
};

uint32_t index_of(const uint32_t i, const uint32_t j, const uint32_t n)
{
	return i * n + j;
}

struct linear_system *parse_file(const char *file)
{
	FILE *f = fopen(file, "r");
	struct linear_system *ret = malloc(sizeof(struct linear_system));
	uint32_t n;

	fscanf(f, "%u", &n);
	if (n == 0) {
		return NULL;
	}
	ret->n = n;
	ret->matrix = malloc(n * n * sizeof(double));
	ret->b = malloc(n * sizeof(double));

	size_t i, j;
	for (i = 0; i < n; ++i) {
		for (j = 0; j < n; ++j) {
			fscanf(f, "%lf", &ret->matrix[index_of(i, j, ret->n)]);
		}
		fscanf(f, "%lf", &ret->b[i]);
	}
	return ret;
}

void ls_free(struct linear_system *ls)
{
	free(ls->matrix);
	free(ls->b);
	free(ls);
}

void ls_to_string(const struct linear_system *ls)
{
	int32_t i, j;

	for (i = 0; i < ls->n; ++i) {
		for (j = 0; j < ls->n; ++j) {
			printf("%lf ", ls->matrix[index_of(i, j, ls->n)]);
			if (j + 1 == ls->n) {
				printf("%lf ", ls->b[i]);
				printf("\n");
			}
		}
	}
}

double _det(const struct linear_system *ls, int32_t i, int32_t j, int32_t n)
{
	double ret;
	if (n == 2) {
		ret = ls->matrix[index_of(i, j, ls->n)] * ls->matrix[index_of(i + 1, (j + 1) % ls->n, ls->n)] -
				ls->matrix[index_of(i + 1, j, ls->n)] * ls->matrix[index_of(i, (j + 1) % ls->n, ls->n)];
		if (j + 1 == ls->n) {
			ret = -ret;
		}
		return ret;
	}

	int32_t jj;
	int32_t k = 1;

	for (jj = j; jj < n; ++jj) {
		ret += k * ls->matrix[index_of(i, jj, ls->n)] * _det(ls, i + 1, (jj + 1) % ls->n, n - 1);
		k *= -1;
	}
	return ret;
}

double det(const struct linear_system *ls)
{
	return _det(ls, 0, 0, ls->n);
}

bool is_invertible(const struct linear_system *ls)
{
	printf("%lf\n", det(ls));
	return det(ls) != 0.0;
}

void invert_lines(const struct linear_system *ls, int32_t x, int32_t y)
{
	int32_t j;
	double tmp;

	for (j = 0; j < ls->n; ++j) {
		tmp = ls->matrix[index_of(x, j, ls->n)];
		ls->matrix[index_of(x, j, ls->n)] = ls->matrix[index_of(y, j, ls->n)];
		ls->matrix[index_of(y, j, ls->n)] = tmp;
	}
	tmp = ls->b[x];
	ls->b[x] = ls->b[y];
	ls->b[y] = tmp;
}

void solve(const struct linear_system *ls)
{
	int32_t i, j, k, step = 0;
	double pivot, m;

	/*if (!is_invertible(ls)) {
		printf("Matrix is not invertible\n");
		return;
	}*/

	printf("\n\n");
	/* Iterate through lines */
	for (k = 0; k < ls->n - 1; ++k) {
		for (i = k; i < ls->n - 1; ++i) {
			if (ls->matrix[index_of(i, k, ls->n)] != 0) {
				pivot = ls->matrix[index_of(i, k, ls->n)];
				break;
			}
		}

		printf("STEP #%d\n", step);
		printf("PIVOT: %lf\n\n", pivot);

		if (i != k) {
			invert_lines(ls, i, k);
		}

		for (i = k + 1; i < ls->n; ++i) {
			m = ls->matrix[index_of(i, k, ls->n)] / ls->matrix[index_of(k, k, ls->n)];
			ls->b[i] -= m * ls->b[k];
			for (j = k; j < ls->n; ++j) {
				ls->matrix[index_of(i, j, ls->n)] -= m * ls->matrix[index_of(k, j, ls->n)];
			}
			printf("m = %lf\n", m);
			ls_to_string(ls);
			printf("\n\n");
		}
		++step;
	}

	double *result = malloc(ls->n * sizeof(double));
	int32_t i_result = ls->n - 1;
	for (k = ls->n - 1; k >= 0; --k) {
		result[k] = ls->b[k];
		for (j = k + 1; j < ls->n; ++j) {
			result[k] -= result[j] * ls->matrix[index_of(k, j, ls->n)];
		}
		result[k] /= ls->matrix[index_of(k, k, ls->n)];
		i_result--;
	}

	char v = 'a';

	for (i = 0; i < ls->n; ++i) {
		printf("%c = %lf\n", v, result[i]);
		++v;
	}
}

int main(int argc, char *argv[])
{
	if (argc != 2) {
		printf("Wrong number of parameters!\n");
		return 1;
	}

	struct linear_system *ls = parse_file(argv[1]);

	if (!ls) {
		printf("Input error\n");
		return 2;
	}
	ls_to_string(ls);
	solve(ls);
	ls_free(ls);

	return 0;
}
