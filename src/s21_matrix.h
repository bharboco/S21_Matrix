#ifndef S21_MATRIX_H
#define S21_MATRIX_H

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define SUCCESS 1
#define FAILURE 0

typedef struct matrix_struct {
  double **matrix;
  int rows;     // вниз
  int columns;  // вправо
} matrix_t;

typedef enum {
  OK = 0,          // OK
  INCORRECT = 1,   // Ошибка, некорректная матрица
  ERROR_CALC = 2,  // Ошибка вычисления
} res;

int s21_create_matrix(int rows, int columns, matrix_t *result);     // +
void s21_remove_matrix(matrix_t *A);                                // +
int s21_eq_matrix(matrix_t *A, matrix_t *B);                        //+
int s21_sum_matrix(matrix_t *A, matrix_t *B, matrix_t *result);     // +
int s21_sub_matrix(matrix_t *A, matrix_t *B, matrix_t *result);     // +
int s21_mult_number(matrix_t *A, double number, matrix_t *result);  //+
int s21_mult_matrix(matrix_t *A, matrix_t *B, matrix_t *result);    //+
int s21_transpose(matrix_t *A, matrix_t *result);                   //+
int s21_calc_complements(matrix_t *A, matrix_t *result);            // +
int s21_determinant(matrix_t *A, double *result);                   // +
int s21_inverse_matrix(matrix_t *A, matrix_t *result);              //+

void s21_initialization(matrix_t *A);
int s21_check_matrix(matrix_t *check);  // +
int s21_minor(matrix_t *A, matrix_t *result, int index_rows,
              int index_columns);  // +
#endif                             // S21_MATRIX_H