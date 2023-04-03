#include <iostream>
#include <gsl/gsl_blas.h>
#include  "io.h"
#include "mpi.h"

#include "../Tools/Timer/Timer.h"
#include "../Tools/Logger/Logger.h"

#define MAX_FILENAME_LENGTH 256
#define DEFAULT_OUTPUT_FILENAME "output.txt"

#define DEFAULT_OUTPUT_FILENAME_PROC1 "process1_output.txt"
#define DEFAULT_OUTPUT_FILENAME_PROC2 "process2_output.txt"
#define DEFAULT_OUTPUT_FILENAME_PROC3 "process3_output.txt"

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
			printf("%lf\t", gsl_matrix_get(matrix, i, j));
		}
		printf("\n");
	}
	fflush(stdout);
}
void output_vector(const gsl_vector* vector, const VECTOR_TYPE type = V_COL)
{
	for (int i = 0; i < vector->size; i++)
	{
		if (type == V_COL)
			printf("%lf\n", gsl_vector_get(vector, i));
		else
			printf("%lf\t", gsl_vector_get(vector, i));
	}
	fflush(stdout);
}

FILLING input_validator(const char* message = NULL)
{
	if (message != NULL)
		printf(message);
	fflush(stdout);

	static FILLING choice;
	choice = F_DEFAULT;

	while (choice != F_RANDOM && choice != F_KEYBOARD && choice != F_EXIT)
	{
		fflush(stdout);
		printf("Enter '1' to fill randomly"
			"\nEnter '2' to fill from keyboard"
			"\nEnter '3' to exit"
			"\nYour choice is ");
		fflush(stdout);
		scanf_s("%d", &choice);

		if (choice != F_RANDOM && choice != F_KEYBOARD && choice != F_EXIT)
			printf("Incorrect input!\n");
	}

	return choice;
}

void generate_random_num_vector(gsl_vector* vector)
{
	for (int i = 0; i < vector->size; i++) {
		//gsl_vector_set(vector, i, rand() % 100 + 1);
		gsl_vector_set(vector, i, (i + 1)* 13 + i);
	}
}
void generate_random_num_matrix(gsl_matrix* matrix) {
	for (int i = 0; i < matrix->size1; i++) {
		for (int j = 0; j < matrix->size1; j++) {
			//gsl_matrix_set(matrix, i, j, rand() % 100 + 1);
			gsl_matrix_set(matrix, i, j, (i + 1 * j + 1 * j + 1)*13+i);
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
			fflush(stdout);
			scanf_s("%lf", gsl_matrix_ptr(matrix, i, j));
		}
		//fflush(stdout);
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
		//output_matrix(matrix);
		break;
	case F_KEYBOARD:
		input_matrix(matrix);
		break;
	case F_EXIT:
		exit(1);
	default:;
	}
};
void control_vector_input(gsl_vector* vector, const char choice)
{
	switch (choice) {
	case F_RANDOM:
		generate_random_num_vector(vector);
		//output_vector(vector);
		break; 
	case F_KEYBOARD:
		//input_vector(vector);
		break;
	case F_EXIT:
		exit(1);
	default:;
	}
};

gsl_vector* mult_matrix_by_vector(gsl_matrix* matrix, const gsl_vector* vector)
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
gsl_vector* calculate_y2(gsl_matrix* matrix, const gsl_vector* vector)
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
		for (int j = 0; j < *matrix_dimension; j++)
		{
			temp = 1 / (i + 1 + pow((j + 1), 2));
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
gsl_matrix* mult_row_by_col_vector(const gsl_vector* row_vector, const gsl_vector* col_vector)
{
	gsl_matrix* result = gsl_matrix_alloc(row_vector->size, row_vector->size);
	gsl_matrix_set_zero(result);
	//1 - col vector, 2 - row vector
	gsl_blas_dger(1.0, col_vector, row_vector, result);
	return  result;
}
gsl_matrix* mult_matrix_by_matrix(const gsl_matrix* matrix1, const gsl_matrix* matrix2)
{
	gsl_matrix* result = gsl_matrix_alloc(matrix1->size1, matrix1->size1);
	gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, matrix1, matrix2, 0.0, result);
	return  result;
}
gsl_matrix* add_matrix_to_matrix(const gsl_matrix* matrix1, const gsl_matrix* matrix2)
{
	gsl_matrix* result = gsl_matrix_alloc(matrix1->size1, matrix1->size1);
	gsl_matrix_memcpy(result, matrix2);
	gsl_matrix_add(result, matrix1);
	return  result;
}
gsl_vector* add_vector_to_vector(const gsl_vector* vector1, const gsl_vector* vector2)
{
	gsl_vector* result = gsl_vector_alloc(vector1->size);
	gsl_vector_memcpy(result, vector1); // копіюємо значення вектора v1 в result
	gsl_vector_add(result, vector2);	// додаємо до result вектор v2
	return result;
}
gsl_matrix* mult_scalar_by_matrix(const double scalar, const gsl_matrix* matrix)
{
	gsl_matrix* result = gsl_matrix_alloc(matrix->size1, matrix->size1);
	gsl_matrix_memcpy(result, matrix);
	gsl_matrix_scale(result, scalar);
	return result;
}
int get_matrix_dimension()
{
	int matrix_dimension = 0;
	while (matrix_dimension < 1)
	{
		printf("Enter the dimension of the square matrix:");
		fflush(stdout);
		scanf_s("%d", &matrix_dimension);
		if (matrix_dimension < 1)
			printf("Incorrect matrix dimension!\n");
	}
	fflush(stdout);
	return matrix_dimension;

}
RESULTS get_intermediate_result()
{
	RESULTS print_intermediate_result = R_INIT;

	while (print_intermediate_result != R_FILE
		&& print_intermediate_result != R_CONSOLE
		&& print_intermediate_result != R_BOTH
		&& print_intermediate_result != R_NONE)
	{
		fflush(stdout);
		printf("\nDo you wanna to output intermediate result?\n");
		printf("Enter '1' to write to file "
			"\nEnter '2' to write to console "
			"\nEnter '3' to write to both console and file"
			"\nEnter '4' to cancel write intermediate results"
			"\nYour choice is ");
		fflush(stdout);
		scanf_s("%d", &print_intermediate_result);
		fflush(stdout);
		if (print_intermediate_result < 1 || print_intermediate_result > 4)
			printf("Incorrect input!\n");
	}
	return print_intermediate_result;
}

char* get_input_file()
{
	char* filename = (char*)malloc(MAX_FILENAME_LENGTH * sizeof(char));

	printf("Enter the input filename: ");
	fflush(stdout);
	scanf_s("%s", filename, MAX_FILENAME_LENGTH);

	// Ensure the filename doesn't exceed the maximum length
	size_t filename_length = strlen(filename);
	if (filename_length >= MAX_FILENAME_LENGTH) {
		printf("Filename is too long.\n");
		free(filename);
		return NULL;
	}

	return filename;

}

gsl_vector* main_calculation(const gsl_vector* y1, const gsl_vector* y2, const gsl_matrix* Y3, logger* my_logger = NULL)
{
	/********CALCULATE MAIN FORMULA STEP BY STEP*********/

	//y1 * y2'
	gsl_matrix* y1_mul_y2_trans = mult_row_by_col_vector(y2, y1);

	//Y3^2 -> Y3 * Y3
	gsl_matrix* Y3_square = mult_matrix_by_matrix(Y3, Y3);

	//Y3^2 * y1
	gsl_vector* Y3_square_mul_y1 = mult_matrix_by_vector(Y3_square, y1);

	//[Y3^2 * y1] + y2
	gsl_vector* Y3_square_mul_y1_add_y2 = add_vector_to_vector(Y3_square_mul_y1, y2);

	//y1' * ([Y3^2 * y1] + y2)
	double y1_trans_mul_Y3_square_mul_y1_add_y2 = NULL;
	gsl_blas_ddot(y1, Y3_square_mul_y1_add_y2, &y1_trans_mul_Y3_square_mul_y1_add_y2);

	//y1' * ([Y3^2 * y1] + y2) * Y3
	gsl_matrix* y1_trans_mul_Y3_square_mul_y1_add_y2_mul_Y3 = mult_scalar_by_matrix(y1_trans_mul_Y3_square_mul_y1_add_y2, Y3);

	//(y1' * ([Y3^2 * y1] + y2) * Y3 + y1*y2')
	gsl_matrix* pre_final_res1 = add_matrix_to_matrix(y1_trans_mul_Y3_square_mul_y1_add_y2_mul_Y3, y1_mul_y2_trans);

	//(y1' * ([Y3^2 * y1] + y2) * Y3 + y1*y2') * y2
	gsl_vector* pre_final_res2 = mult_matrix_by_vector(pre_final_res1, y2);

	//y2' + [(y1' * ([Y3^2 * y1] + y2) * Y3 + y1*y2') * y2]'
	gsl_vector* final_result = add_vector_to_vector(y2, pre_final_res2);


	if (my_logger->func != NULL)
	{

		log_vector(my_logger, y2, "y2'", V_ROW);
		log_matrix(my_logger, y1_mul_y2_trans, "y1 * y2'");

		log_matrix(my_logger, Y3_square, "Y3^2");

		log_vector(my_logger, Y3_square_mul_y1, "Y3^2 * y1");

		log_vector(my_logger, Y3_square_mul_y1_add_y2, "[Y3^2 * y1] + y2");

		char buf[100];
		sprintf_s(buf, 100, "\nNumber y1' * ([Y3^2 * y1] + y2): %lf \n", y1_trans_mul_Y3_square_mul_y1_add_y2);
		log_result(my_logger, buf);

		log_matrix(my_logger, y1_trans_mul_Y3_square_mul_y1_add_y2_mul_Y3, "y1' * ([Y3^2 * y1] + y2) * Y3");
		log_matrix(my_logger, pre_final_res1, "(y1' * ([Y3^2 * y1] + y2) * Y3 + y1*y2')");

		log_vector(my_logger, pre_final_res2, "(y1' * ([Y3^2 * y1] + y2) * Y3 + y1*y2') * y2");

		log_vector(my_logger, pre_final_res2, "((y1' * ([Y3^2 * y1] + y2) * Y3 + y1*y2') * y2)'", V_ROW);
	}

	gsl_matrix_free(y1_mul_y2_trans);
	gsl_matrix_free(Y3_square);
	gsl_vector_free(Y3_square_mul_y1);
	gsl_vector_free(Y3_square_mul_y1_add_y2);
	gsl_matrix_free(y1_trans_mul_Y3_square_mul_y1_add_y2_mul_Y3);
	gsl_matrix_free(pre_final_res1);
	gsl_vector_free(pre_final_res2);

	return final_result;
}

int main(int argc, char* argv[])
{
	MPI_Status Status;
	int matrix_dimension;
	RESULTS print_intermediate_result;
	logger* my_logger = NULL;

	MPI_Init(&argc, &argv);
	int rank, size;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);

	if(rank == 0)
	{
		matrix_dimension = get_matrix_dimension();
		print_intermediate_result = get_intermediate_result();
	}

	MPI_Bcast(&matrix_dimension, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&print_intermediate_result, 1, MPI_INT, 0, MPI_COMM_WORLD);


	if (rank == 0)
	{
		if (print_intermediate_result == R_FILE || print_intermediate_result == R_BOTH)
			my_logger = create_logger(print_intermediate_result, DEFAULT_OUTPUT_FILENAME_PROC1);
		else my_logger = create_logger(print_intermediate_result);

		//Allocate matrix and vector
		gsl_matrix* A = gsl_matrix_alloc(matrix_dimension, matrix_dimension);
		gsl_matrix* A1 = gsl_matrix_alloc(matrix_dimension, matrix_dimension);
		gsl_vector* b1 = gsl_vector_alloc(matrix_dimension);
		gsl_vector* c1 = gsl_vector_alloc(matrix_dimension);
		gsl_matrix* A2 = gsl_matrix_alloc(matrix_dimension, matrix_dimension);
		gsl_matrix* B2 = gsl_matrix_alloc(matrix_dimension, matrix_dimension);

		/********INPUT VALUES*********/
		printf("\n**********************************");
		control_matrix_input(A, input_validator("\nEnter the values of the matrix A:\n"));
		printf("\n**********************************");

		control_matrix_input(A1, input_validator("\nEnter the values of the matrix A1:\n"));
		printf("\n**********************************");
		//gsl_matrix_free(A1);

		control_vector_input(b1, input_validator("\nEnter column-vector b1:\n"));
		printf("\n**********************************");
		//gsl_vector_free(b1);

		control_vector_input(c1, input_validator("\nEnter column-vector c1:\n"));
		printf("\n**********************************");
		//gsl_vector_free(c1);

		control_matrix_input(A2, input_validator("\nEnter the values of the matrix A2:\n"));
		printf("\n**********************************");
		MPI_Send(A2->data, matrix_dimension * matrix_dimension, MPI_DOUBLE, 1, 0, MPI_COMM_WORLD);
		//gsl_matrix_free(A2);

		control_matrix_input(B2, input_validator("\nEnter the values of the matrix B2:\n"));
		printf("\n**********************************");
		fflush(stdout);
		MPI_Send(B2->data, matrix_dimension * matrix_dimension, MPI_DOUBLE, 1, 0, MPI_COMM_WORLD);
		//gsl_matrix_free(B2);
		/*******************************************************************************************/
		gsl_vector* Y3_square_mul_y1 = gsl_vector_alloc(matrix_dimension);
		gsl_matrix* Y3 = gsl_matrix_alloc(matrix_dimension, matrix_dimension);
		/*calculation in first process*/

		timerify* timer = createTimer();
		startTimer(timer);

		gsl_vector* b = calculate_b(&matrix_dimension);//1
		gsl_vector* y1 = mult_matrix_by_vector(A, b);//2
		//send y1 to 1 proc
		MPI_Send(y1->data, matrix_dimension, MPI_DOUBLE, 1, 0, MPI_COMM_WORLD);
		gsl_vector* b12_min_c1 = calculate_12b1_minus_c1(b1, c1);//3
		gsl_vector* y2 = mult_matrix_by_vector(A1, b12_min_c1);//4
		//y1 * y2'
		gsl_matrix* y1_mul_y2_trans = mult_row_by_col_vector(y2, y1);//5

		MPI_Recv(Y3_square_mul_y1->data, matrix_dimension, MPI_DOUBLE, 1, 0, MPI_COMM_WORLD, &Status);

		//
			//[Y3^2 * y1] + y2
		gsl_vector* Y3_square_mul_y1_add_y2 = add_vector_to_vector(Y3_square_mul_y1, y2);

		//y1' * ([Y3^2 * y1] + y2)
		double y1_trans_mul_Y3_square_mul_y1_add_y2 = NULL;
		gsl_blas_ddot(y1, Y3_square_mul_y1_add_y2, &y1_trans_mul_Y3_square_mul_y1_add_y2);

		MPI_Recv(Y3->data, matrix_dimension * matrix_dimension, MPI_DOUBLE, 1, 0, MPI_COMM_WORLD, &Status);

		//y1' * ([Y3^2 * y1] + y2) * Y3
		gsl_matrix* y1_trans_mul_Y3_square_mul_y1_add_y2_mul_Y3 = mult_scalar_by_matrix(y1_trans_mul_Y3_square_mul_y1_add_y2, Y3);

		//(y1' * ([Y3^2 * y1] + y2) * Y3 + y1*y2')
		gsl_matrix* pre_final_res1 = add_matrix_to_matrix(y1_trans_mul_Y3_square_mul_y1_add_y2_mul_Y3, y1_mul_y2_trans);

		//(y1' * ([Y3^2 * y1] + y2) * Y3 + y1*y2') * y2
		gsl_vector* pre_final_res2 = mult_matrix_by_vector(pre_final_res1, y2);

		//y2' + [(y1' * ([Y3^2 * y1] + y2) * Y3 + y1*y2') * y2]'
		gsl_vector* final_result = add_vector_to_vector(y2, pre_final_res2);

		stopTimer(timer);

		char buf2[100];
		sprintf_s(buf2, 100, "Elapsed time: %f seconds\n\n", getElapsedSeconds(timer));
		
		if (my_logger->func != NULL) log_result(my_logger, buf2);
		printf("%s", buf2);

		/********************Logging and release memory**********************/
		if (my_logger->func != NULL)
		{
			log_matrix(my_logger, A, "input A:");
			log_vector(my_logger, b1, "input b1:");
			log_matrix(my_logger, A1, "input A1:");
			log_vector(my_logger, c1, "input c1:");
			//log_matrix(my_logger, A2, "input A2:");
			//log_matrix(my_logger, B2, "input B2:");
			log_vector(my_logger, b, "b: bi = i^2 / 12 for even and bi = i for odd\n b:");
			log_vector(my_logger, y1, "y1 = A * b\n y1:");
			log_vector(my_logger, b12_min_c1, "b12 minus c1 =:");
			log_vector(my_logger, y2, "y2 = A1 * (12b1 - c1)\n y2=:");
			log_vector(my_logger, y2, "y2'", V_ROW);
			log_matrix(my_logger, y1_mul_y2_trans, "y1 * y2'");
			log_vector(my_logger, Y3_square_mul_y1_add_y2, "[Y3^2 * y1] + y2");

			char buf[100];
			sprintf_s(buf, 100, "\nNumber y1' * ([Y3^2 * y1] + y2): %lf \n", y1_trans_mul_Y3_square_mul_y1_add_y2);
			log_result(my_logger, buf);

			log_matrix(my_logger, y1_trans_mul_Y3_square_mul_y1_add_y2_mul_Y3, "y1' * ([Y3^2 * y1] + y2) * Y3");
			log_matrix(my_logger, pre_final_res1, "(y1' * ([Y3^2 * y1] + y2) * Y3 + y1*y2')");

			log_vector(my_logger, pre_final_res2, "(y1' * ([Y3^2 * y1] + y2) * Y3 + y1*y2') * y2");

			log_vector(my_logger, pre_final_res2, "((y1' * ([Y3^2 * y1] + y2) * Y3 + y1*y2') * y2)'", V_ROW);

			log_vector(my_logger, final_result, "result y2' + [(y1' * ([Y3^2 * y1] + y2) * Y3 + y1*y2') * y2]':", V_ROW);
		}
		gsl_vector_free(b1);
		gsl_vector_free(c1);
		gsl_matrix_free(A1);
		gsl_vector_free(b12_min_c1);
		gsl_vector_free(y2);
		gsl_vector_free(y1);
		gsl_matrix_free(y1_mul_y2_trans);
		gsl_vector_free(Y3_square_mul_y1);
		gsl_vector_free(Y3_square_mul_y1_add_y2);
		gsl_matrix_free(y1_trans_mul_Y3_square_mul_y1_add_y2_mul_Y3);
		gsl_matrix_free(pre_final_res1);
		gsl_vector_free(pre_final_res2);
		gsl_vector_free(final_result);
		gsl_matrix_free(Y3);
		gsl_vector_free(b);
		gsl_matrix_free(A);
		destroy_logger(my_logger);
	}
	if(rank == 1)
	{
		if (print_intermediate_result == R_FILE || print_intermediate_result == R_BOTH)
			my_logger = create_logger(print_intermediate_result, DEFAULT_OUTPUT_FILENAME_PROC2);
		else my_logger = create_logger(print_intermediate_result);

		gsl_vector* y1 = gsl_vector_alloc(matrix_dimension);
		gsl_matrix* A2 = gsl_matrix_alloc(matrix_dimension, matrix_dimension);
		gsl_matrix* B2 = gsl_matrix_alloc(matrix_dimension, matrix_dimension);
		MPI_Recv(A2->data, matrix_dimension * matrix_dimension, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &Status);
		MPI_Recv(B2->data, matrix_dimension * matrix_dimension, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &Status);

		gsl_matrix* C2 = calculate_C2(&matrix_dimension);//1
		gsl_matrix* B2_min_C2 = gsl_matrix_alloc(A2->size1, A2->size1);//2
		MPI_Recv(y1->data, matrix_dimension, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &Status);
		gsl_matrix_memcpy(B2_min_C2, B2); // B2_min_C2 <- B2
		gsl_matrix_sub(B2_min_C2, C2); // B2_min_C2 = B2_min_C2 - C2 
		gsl_matrix* Y3 = mult_matrix_by_matrix(A2, B2_min_C2);//3

		//Y3^2 -> Y3 * Y3
		gsl_matrix* Y3_square = mult_matrix_by_matrix(Y3, Y3);//4

		//Y3^2 * y1
		gsl_vector* Y3_square_mul_y1 = mult_matrix_by_vector(Y3_square, y1);//5
		MPI_Send(Y3_square_mul_y1->data, matrix_dimension, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
		MPI_Send(Y3->data, matrix_dimension * matrix_dimension, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);

		/********************Logging and release memory**********************/
	if (my_logger->func != NULL)
		{
			log_matrix(my_logger, A2, "input A2:");
			log_matrix(my_logger, B2, "input B2:");
			log_matrix(my_logger, C2, "Cij = 1 / (i+j^2)\n C2=:");
			log_matrix(my_logger, B2_min_C2, "B2 - C2 =:");
			log_matrix(my_logger, Y3, "Y3 = A2 * (B2 - C2)\n Y3=:");
			log_matrix(my_logger, Y3_square, "Y3^2");
			log_vector(my_logger, Y3_square_mul_y1, "Y3^2 * y1");
		}
		gsl_matrix_free(A2);
		gsl_matrix_free(B2);
		gsl_matrix_free(C2);
		gsl_matrix_free(B2_min_C2);
		gsl_matrix_free(Y3);
		gsl_matrix_free(Y3_square);
		gsl_vector_free(Y3_square_mul_y1);
		destroy_logger(my_logger);


		/*
		log_matrix(my_logger, A1, "input A1:");
		log_vector(my_logger, b1, "input b1:");
		log_vector(my_logger, c1, "input c1:");
		log_vector(my_logger, b12_min_c1, "b12 minus c1 =:");
		log_vector(my_logger, y2, "y2 = A1 * (12b1 - c1)\n y2=:");
		log_vector(my_logger, y2, "y2'", V_ROW);
		log_matrix(my_logger, y1_mul_y2_trans, "y1 * y2'");
		log_vector(my_logger, Y3_square_mul_y1_add_y2, "[Y3^2 * y1] + y2");

		char buf[100];
		sprintf_s(buf, 100, "\nNumber y1' * ([Y3^2 * y1] + y2): %lf \n", y1_trans_mul_Y3_square_mul_y1_add_y2);
		log_result(my_logger, buf);

		log_matrix(my_logger, y1_trans_mul_Y3_square_mul_y1_add_y2_mul_Y3, "y1' * ([Y3^2 * y1] + y2) * Y3");
		log_matrix(my_logger, pre_final_res1, "(y1' * ([Y3^2 * y1] + y2) * Y3 + y1*y2')");

		log_vector(my_logger, pre_final_res2, "(y1' * ([Y3^2 * y1] + y2) * Y3 + y1*y2') * y2");

		log_vector(my_logger, pre_final_res2, "((y1' * ([Y3^2 * y1] + y2) * Y3 + y1*y2') * y2)'", V_ROW);

		log_vector(my_logger, final_result, "result y2' + [(y1' * ([Y3^2 * y1] + y2) * Y3 + y1*y2') * y2]':", V_ROW);

		gsl_vector_free(b1);
		gsl_vector_free(c1);
		gsl_matrix_free(A1);
		gsl_vector_free(b12_min_c1);
		gsl_vector_free(y2);
		gsl_vector_free(y1);
		gsl_matrix_free(y1_mul_y2_trans);
		gsl_vector_free(Y3_square_mul_y1);
		gsl_vector_free(Y3_square_mul_y1_add_y2);
		gsl_matrix_free(y1_trans_mul_Y3_square_mul_y1_add_y2_mul_Y3);
		gsl_matrix_free(pre_final_res1);
		gsl_vector_free(pre_final_res2);
		gsl_vector_free(final_result);
		gsl_matrix_free(Y3);
		destroy_logger(my_logger);
		*/
	}

	MPI_Finalize();
	return 0;
}