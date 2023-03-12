#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_vector.h>

char input_validator(const char* message = NULL) 
{
	if (message != NULL)
		printf(message);
	static char choice = 0;
	while (choice != 1 && choice != 2 && choice != 3)
	{
		printf("Enter '1' to fill randomly"
			   "\nEnter '2' to fill from keyboard"
			   "\nEnter '3' to exit"
			   "\nYour choice is ");
		scanf_s("%d", &choice);
	
		if (choice < 1 || choice > 3) 
			printf("Incorect input!\n");
	}

	return choice;
}
void generate_random_num_vector(gsl_vector* vector) 
{
	for (int i = 0; i < vector->size; i++) {
		gsl_vector_set(vector, i, rand() % 100 + 1);
	}
}
void generate_random_num_matrix(gsl_matrix* matrix) {
	for (int i = 0; i < matrix->size1; i++) {
		for (int j = 0; j < matrix->size1; j++) {
			gsl_matrix_set(matrix, i, j, rand() % 100 + 1);
		}
	}
}
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
void input_vector(gsl_vector* vector) {
	
	for (int j = 0; j < vector->size; j++) 
	{
		printf("Row %d: ", j + 1);
		scanf_s("%lf", gsl_vector_ptr(vector, j)); // введення елементу матриці
	}
}
void output_matrix(gsl_matrix* matrix) 
{
	for (int i = 0; i < matrix->size1; i++) {
		for (int j = 0; j < matrix->size1; j++) {
			printf("%g\t", gsl_matrix_get(matrix, i, j));
		}
		printf("\n");
	}
}
void output_vector(gsl_vector* vector) 
{
	for (int i = 0; i < vector->size; i++) {
		printf("%0.2lf\n", gsl_vector_get(vector, i));
	}
}
void control_matrix_input(gsl_matrix* matrix, char choice) {
	switch (choice) {
	case 1:
		generate_random_num_matrix(matrix);
		output_matrix(matrix);
		break;
	case 2:
		input_matrix(matrix);
		break;
	case 3:
		exit(1);
	};
};
void control_vector_input(gsl_vector* vector, char choice) {
	switch (choice) {
	case 1:
		generate_random_num_vector(vector);
		output_vector(vector);
		break;
	case 2:
		input_vector(vector);
		break;
	case 3:
		exit(1);
	};
};
gsl_vector* calculate_y1(gsl_matrix* matrix, gsl_vector* vector, int* matrix_dimension) {

	gsl_vector* result = gsl_vector_alloc(*matrix_dimension);
	double alpha = 1.0;
	double beta = 0.0;
	gsl_blas_dgemv(CblasNoTrans, alpha, matrix, vector, beta, result);
	return result;
}
gsl_vector* calculate_b(int* matrix_dimension) {

	gsl_vector* v = gsl_vector_alloc(*matrix_dimension);
	for (int i = 1; i <= v->size; i++) {
		if (i % 2 == 0) {
			gsl_vector_set(v, i - 1, pow(i, 1) / 12);
		}
		else {
			gsl_vector_set(v, i - 1, i);
		}
	}

	return v;
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
	control_matrix_input(A, input_validator("\nEnter the values of the matrix A:\n"));
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
