
#include <cstring>

#include "logger.h"

void log_to_file(const char* message, FILE* file)
{
    fprintf(file, "%s", message);
}

void log_to_console(const char* message, FILE* file)
{
    printf("%s", message);
}

void log_to_both(const char* message, FILE* file)
{
    log_to_file(message, file);
    log_to_console(message);
}

logger* create_logger(RESULTS log_type, const char* filename)
{
    logger* log = (logger*)malloc(sizeof(logger));

    switch (log_type)
    {
    case R_FILE:
        log->func = log_to_file;
        fopen_s(&log->file, filename, "a");
        break;
    case R_CONSOLE:
        log->func = log_to_console;
        log->file = NULL;
        break;
    case R_BOTH:
        log->func = log_to_both;
        fopen_s(&log->file, filename, "a");
        break;
    default:
        log->func = NULL;
        log->file = NULL;
        break;
    case R_NONE:
        break;
    case R_INIT:
        break;
    }

    return log;
}

void log_result(logger* log, const char* message)
{
    log->func(message, log->file);
    if (log->file)
    {
        fflush(log->file);
    }
}

void destroy_logger(logger* log)
{
    if (log->file)
    {
        fclose(log->file);
    }
    free(log);
}

void log_matrix(logger* log, const gsl_matrix* matrix, const char* name)
{
    const size_t buf_size = 1000;
    char* message = (char*)malloc(buf_size);
    sprintf_s(message, buf_size, "Matrix %s:\n", name);
    log_result(log, message);
    for (size_t i = 0; i < matrix->size1; i++)
    {
        sprintf_s(message, buf_size, "  [ ");
        for (size_t j = 0; j < matrix->size2; j++)
        {
            sprintf_s(message + strlen(message), buf_size - strlen(message), "%.2f ", gsl_matrix_get(matrix, i, j));
        }
        sprintf_s(message + strlen(message), buf_size - strlen(message), "]\n");
        log_result(log, message);
    }
    free(message);
}

void log_vector(logger* log, const gsl_vector* vector, const char* name) {
    const size_t buf_size = 1000;
    char* message = (char*)malloc(buf_size);
    sprintf_s(message, buf_size, "Vector %s: [ ", name);
    for (size_t i = 0; i < vector->size; i++) {
        sprintf_s(message + strlen(message), buf_size - strlen(message), "%.2f ", gsl_vector_get(vector, i));
    }
    sprintf_s(message + strlen(message), buf_size - strlen(message), "]\n");
    log_result(log, message);
    free(message);
}