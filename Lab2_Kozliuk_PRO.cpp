#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>

void input_matrix(gsl_matrix* matrix)
{
	for (int i = 0; i < matrix->size1; i++) 
	{
		printf("Row %d: ", i + 1);
		for (int j = 0; j < matrix->size2; j++) 
		{
			scanf_s("%lf", gsl_matrix_ptr(matrix, i, j)); // введення елементу матриці
		}
	}
}

void output_matrix(gsl_matrix* matrix) 
{
	for (int i = 0; i < matrix->size1; i++) {
		for (int j = 0; j < matrix->size1; j++) {
			printf("%g ", gsl_matrix_get(matrix, i, j));
		}
		printf("\n");
	}
}

gsl_vector* calculate_b(int *matrix_dimension) 
{

	gsl_vector* v = gsl_vector_alloc(*matrix_dimension);
	for (int i = 1; i <= v->size; i++) 
	{
		if (i % 2 == 0) 
		{
			gsl_vector_set(v, i - 1, pow(i, 1)/12);
		}
		else 
		{
			gsl_vector_set(v, i - 1, i);
		}
	}
	
	return v;
}

void output_vector(gsl_vector* vector) 
{
	for (int i = 0; i < vector->size; i++) {
		printf("%0.2lf\n", gsl_vector_get(vector, i));
	}
}

gsl_vector* calculate_y1(gsl_matrix* matrix, gsl_vector* vector, int* matrix_dimension) {

	gsl_vector* result = gsl_vector_alloc(*matrix_dimension);
	double alpha = 1.0;
	double beta = 0.0;
	gsl_blas_dgemv(CblasNoTrans, alpha, matrix, vector, beta, result);
	return result;
}



int main() {
	int matrix_dimension = 0;
	
	while (matrix_dimension < 1)
	{
		printf("Enter the dimension of the square matrix:");
		scanf_s("%d", &matrix_dimension);
		if (matrix_dimension < 1)
			printf("Incorect matrix dimension!\n");

	}

	gsl_matrix* A = gsl_matrix_alloc(matrix_dimension, matrix_dimension); 
	gsl_vector* b = calculate_b(&matrix_dimension);					 
	gsl_vector* y1 = NULL;

	gsl_matrix* A1 = gsl_matrix_alloc(matrix_dimension, matrix_dimension);
	gsl_vector* b1 = gsl_vector_alloc(matrix_dimension);
	gsl_vector* c1 = gsl_vector_alloc(matrix_dimension);
	gsl_vector* y2 = NULL;


	gsl_matrix* A2 = gsl_matrix_alloc(matrix_dimension, matrix_dimension);
	gsl_matrix* B2 = gsl_matrix_alloc(matrix_dimension, matrix_dimension);
	gsl_matrix* C2 = NULL;
	gsl_matrix* Y3 = NULL;


	printf("**********************************");
	printf("Enter the values of the matrix A:\n");
	input_matrix(A);
	printf("**********************************");
	printf("Formula b: bi = i^2 / 12 for even and bi = i for odd\n");
	printf("Vector b:\n");
	output_vector(b);
	printf("**********************************");
	 y1 = calculate_y1(A, b, &matrix_dimension);
	printf("Formula y1 = A * b\n y1:\n");
	output_vector(y1);
	printf("**********************************");





	return 0;
}
