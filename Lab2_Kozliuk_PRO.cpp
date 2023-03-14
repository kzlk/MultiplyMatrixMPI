#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_vector.h>

#include "Timer.h"

typedef void (*log_callback)(const char*);
log_callback logger = NULL;

enum VECTOR_TYPE
{
	V_ROW,
	V_COL
};
enum RESULTS
{
	R_FILE = 1,
	R_CONSOLE,
	R_BOTH,
	R_NONE,
	R_INIT
};
enum FILLING
{
	F_RANDOM = 1,
	F_KEYBOARD,
	F_EXIT,
	F_DEFAULT
};
void output_matrix(const gsl_matrix* matrix)
{
	for (int i = 0; i < matrix->size1; i++)
	{
		for (int j = 0; j < matrix->size1; j++)
		{
			printf("%0.2lf\t", gsl_matrix_get(matrix, i, j));
		}
		printf("\n");
	}
}
void output_vector(const gsl_vector* vector, const VECTOR_TYPE type = V_COL)
{
		for (int i = 0; i < vector->size; i++)
		{
			if(type == V_COL)
				printf("%0.2lf\n", gsl_vector_get(vector, i));
			else
				printf("%0.2lf\t", gsl_vector_get(vector, i));
		}
}

FILLING input_validator(const char* message = NULL)
{
	if (message != NULL)
		printf(message);


	static FILLING choice;
	choice = F_DEFAULT;

	while (choice != F_RANDOM && choice != F_KEYBOARD && choice != F_EXIT)
	{
		printf("Enter '1' to fill randomly"
			   "\nEnter '2' to fill from keyboard"
			   "\nEnter '3' to exit"
			   "\nYour choice is ");
		scanf_s("%d", &choice);
	
		if (choice != F_RANDOM && choice != F_KEYBOARD && choice != F_EXIT)
			printf("Incorrect input!\n");
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
			scanf_s("%lf", gsl_matrix_ptr(matrix, i, j)); 
		}
	}
} 
void input_vector(gsl_vector* vector)
{
	
	for (int j = 0; j < vector->size; j++) 
	{
		printf("Row %d: ", j + 1);
		scanf_s("%lf", gsl_vector_ptr(vector, j)); 
	}
}


void control_matrix_input(gsl_matrix* matrix, const char choice)
{
	switch (choice)
	{
	case F_RANDOM:
		generate_random_num_matrix(matrix);
		output_matrix(matrix);
		break;
	case F_KEYBOARD:
		input_matrix(matrix);
		break;
	case F_EXIT:
		exit(1);
	default: ;
	}
};
void control_vector_input(gsl_vector* vector, const char choice)
{
	switch (choice) {
	case F_RANDOM:
		generate_random_num_vector(vector);
		output_vector(vector);
		break;
	case F_KEYBOARD:
		input_vector(vector);
		break;
	case F_EXIT:
		exit(1);
	default: ;
	}
};

gsl_vector* mult_matrix_by_vector( gsl_matrix* matrix, const gsl_vector* vector)
{
	gsl_vector* result = gsl_vector_alloc(matrix->size1);
	double alpha = 1.0;
	double beta = 0.0;
	gsl_blas_dgemv(CblasNoTrans, alpha, matrix, vector, beta, result);
	return result;
}
gsl_vector* calculate_b(int* matrix_dimension)
{
	gsl_vector* v = gsl_vector_alloc(*matrix_dimension);
	for (int i = 1; i <= v->size; i++) {
		if (i % 2 == 0) {
			gsl_vector_set(v, i - 1, pow(i, 2) / 12);
		}
		else {
			gsl_vector_set(v, i - 1, i);
		}
	}

	return v;
}
gsl_vector* calculate_12b1_minus_c1(const gsl_vector* b1, const gsl_vector* c1)
{
	gsl_vector* result = gsl_vector_alloc(b1->size); // where n is the size of the vectors

	double alpha = 12.0; // coefficient for b1

	gsl_vector_memcpy(result, b1); // copy b1 to result vector

	gsl_vector_scale(result, alpha); // multiply result by 12

	gsl_vector_sub(result, c1); // subtract c1 from result

	return result;
}
gsl_vector* calculate_y2( gsl_matrix* matrix, const gsl_vector* vector)
{
	gsl_vector* result = gsl_vector_alloc(matrix->size1);
	double alpha = 1.0;
	double beta = 0.0;
	gsl_blas_dgemv(CblasNoTrans, alpha, matrix, vector, beta, result);
	return result;
}
gsl_matrix* calculate_C2(const int* matrix_dimension)
{
	gsl_matrix* matrix = gsl_matrix_alloc(*matrix_dimension, *matrix_dimension);
	double temp = NULL;
	for (int i = 0; i < *matrix_dimension; i++)
	{
		for(int j = 0; j < *matrix_dimension; j++)
		{
			temp = 1 / (i + 1 + pow((j + 1),2));
			gsl_matrix_set(matrix, i, j, temp);
		}
	}
	return matrix;
}
gsl_matrix* calculate_Y3(gsl_matrix* A2, gsl_matrix* B2, gsl_matrix* C2)
{
	gsl_matrix* Y3 = gsl_matrix_alloc(A2->size1, A2->size1);

	gsl_matrix_memcpy(Y3, B2); // Y3 <- B2

	gsl_matrix_sub(Y3, C2); // Y3 = Y3 - C2 

	gsl_matrix_mul_elements(Y3, A2); // Y3 = Y3 * B2

	return Y3;

}



int get_matrix_dimension()
{
	int matrix_dimension = 0;
	while (matrix_dimension < 1)
	{
		printf("Enter the dimension of the square matrix:");
		scanf_s("%d", &matrix_dimension);
		if (matrix_dimension < 1)
			printf("Incorrect matrix dimension!\n");
	}

	return matrix_dimension;

}
RESULTS get_intermediate_result()
{
	RESULTS print_intermediate_result = R_INIT;

	printf("**********************************");

	while (print_intermediate_result != R_FILE
		&& print_intermediate_result != R_CONSOLE
		&& print_intermediate_result != R_BOTH
		&& print_intermediate_result != R_NONE)
	{
		printf("\nDo you wanna to output intermediate result?\n");
		printf("Enter '1' to write to file "
			"\nEnter '2' to write to console "
			"\nEnter '3' to write to both console and file"
			"\nEnter '4' to cancel write intermediate results"
			"\nYour choice is ");
		scanf_s("%d", &print_intermediate_result);

		if (print_intermediate_result < 1 || print_intermediate_result > 4)
			printf("Incorrect input!\n");
	}
	return print_intermediate_result;
}

gsl_vector* main_calculation()
{


	return NULL;
}


int main() {
	
	int matrix_dimension = get_matrix_dimension();
	RESULTS print_intermediate_result = get_intermediate_result();

	//function = NULL
	//switch result

	gsl_matrix* A = gsl_matrix_alloc(matrix_dimension, matrix_dimension); 
	gsl_vector* b = NULL;					 
	gsl_vector* y1 = NULL;

	gsl_matrix* A1 = gsl_matrix_alloc(matrix_dimension, matrix_dimension);
	gsl_vector* b1 = gsl_vector_alloc(matrix_dimension);
	gsl_vector* c1 = gsl_vector_alloc(matrix_dimension);
	gsl_vector* y2 = NULL;

	gsl_matrix* A2 = gsl_matrix_alloc(matrix_dimension, matrix_dimension);
	gsl_matrix* B2 = gsl_matrix_alloc(matrix_dimension, matrix_dimension);
	gsl_matrix* C2 = NULL;
	gsl_matrix* Y3 = NULL;

	/********INPUT VALUES*********/
	printf("\n**********************************");
	control_matrix_input(A, input_validator("\nEnter the values of the matrix A:\n"));
	printf("\n**********************************");

	control_matrix_input(A1, input_validator("\nEnter the values of the matrix A1:\n"));
	printf("\n**********************************");

	control_vector_input(b1, input_validator("\nEnter column-vector b1:\n"));
	printf("\n**********************************");

	control_vector_input(c1, input_validator("\nEnter column-vector c1:\n"));
	printf("\n**********************************");

	control_matrix_input(A2, input_validator("\nEnter the values of the matrix A2:\n"));
	printf("\n**********************************");

	control_matrix_input(B2, input_validator("\nEnter the values of the matrix B2:\n"));
	printf("\n**********************************");

	/********CALCULATE VALUES FOR MAIN FORMULA*********/
	printf("\nCalculated values for main formula\n");

	//calculate b1
	printf("\nFormula b: bi = i^2 / 12 for even and bi = i for odd\n");
	printf("Vector b:\n");
	b = calculate_b(&matrix_dimension);
	output_vector(b);
	printf("\n**********************************");

	//calculate y1
	y1 = mult_matrix_by_vector(A, b);
	printf("\nFormula y1 = A * b\n y1:\n");
	output_vector(y1);
	printf("\n**********************************");

	//for intermediate results
	gsl_vector* b12_min_c1 = calculate_12b1_minus_c1(b1, c1);

	//calculate y2
	y2 = mult_matrix_by_vector(A1, b12_min_c1);
	printf("\nFormula y2 = A1 * (12b1 - c1)\n y2=:\n");
	output_vector(y2);
	printf("\n**********************************");

	//calculate C2; Cij = 1 / (i+j^2)
	printf("\nFormula Cij = 1 / (1+j^2)\n C2=:\n");
	C2 = calculate_C2(&matrix_dimension);
	output_matrix(C2);
	printf("\n**********************************");

	//calculate Y3 = A2 * (B2 - C2)
	printf("\nFormula Y3 = A2 * (B2 - C2)\n Y2=:\n");
	Y3 = calculate_Y3(A2, B2, C2);
	output_matrix(Y3);
	printf("\n**********************************");


	/********CALCULATE MAIN FORMULA STEP BY STEP*********/

	//gsl_blas_dger

		

		

	/*Free memory*/


	return 0;
}
