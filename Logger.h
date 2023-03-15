#pragma once
#include <cstdio>
#include <gsl/gsl_matrix_double.h>
#include <gsl/gsl_vector.h>

typedef void (*log_func)(const char*, FILE* file);

enum RESULTS
{
    R_FILE = 1,
    R_CONSOLE,
    R_BOTH,
    R_NONE,
    R_INIT
};

typedef struct {
    log_func func;
    FILE* file;
} logger;

void log_to_file(const char* message, FILE* file);
void log_to_console(const char* message, FILE* file = NULL);
void log_to_both(const char* message, FILE* file);

logger* create_logger(RESULTS log_type, const char* filename = NULL);
void log_result(logger* log, const char* message);
void destroy_logger(logger* log);

void log_matrix(logger* log, const gsl_matrix* matrix, const char* name);
void log_vector(logger* log, const gsl_vector* vector, const char* name);
