#include "s21_matrix.h"

int s21_create_matrix(int rows, int columns, matrix_t *result) {
  if (result == NULL || rows < 1 || columns < 1) return INCORRECT;
  int res = OK;
  result->rows = rows;
  result->columns = columns;
  result->matrix =
      malloc(rows * columns * sizeof(double) + rows * sizeof(double));
  if (result->matrix != NULL) {
    double *pointer = (double *)(result->matrix + rows);
    for (int i = 0; i < rows; i++) {
      result->matrix[i] = pointer + columns * i;
    }
  } else
    res = INCORRECT;
  return res;
}

void s21_remove_matrix(matrix_t *A) {
  free(A->matrix);
  A->rows = 0, A->columns = 0, A->matrix = NULL;
}

int s21_eq_matrix(matrix_t *A, matrix_t *B) {
  if (s21_check_matrix(A) || s21_check_matrix(B) || A->rows != B->rows ||
      A->columns != B->columns) {
    return FAILURE;
  }
  int res = SUCCESS;
  for (int i = 0; i < A->rows; i++) {
    for (int j = 0; j < A->columns; j++) {
      if (fabs(A->matrix[i][j] - B->matrix[i][j]) > 1e-7) {  // 0.0000001
        res = FAILURE;
      }
    }
  }
  return res;
}

int s21_sum_matrix(matrix_t *A, matrix_t *B, matrix_t *result) {
  int res = OK;
  if (s21_check_matrix(A) || s21_check_matrix(B)) return INCORRECT;
  if (A->rows == B->rows && A->columns == B->columns) {
    if (s21_create_matrix(A->rows, A->columns, result) == OK) {
      for (int i = 0; i < A->rows; i++) {
        for (int j = 0; j < A->columns; j++) {
          result->matrix[i][j] = A->matrix[i][j] + B->matrix[i][j];
        }
      }
    } else
      res = INCORRECT;
  } else
    res = ERROR_CALC;
  return res;
}

int s21_sub_matrix(matrix_t *A, matrix_t *B, matrix_t *result) {
  int res = OK;
  if (s21_check_matrix(A) || s21_check_matrix(B)) return INCORRECT;
  if (A->rows == B->rows && A->columns == B->columns) {
    if (s21_create_matrix(A->rows, A->columns, result) == OK) {
      for (int i = 0; i < A->rows; i++) {
        for (int j = 0; j < A->columns; j++) {
          result->matrix[i][j] = A->matrix[i][j] - B->matrix[i][j];
        }
      }
    } else
      res = INCORRECT;
  } else
    res = ERROR_CALC;
  return res;
}

int s21_mult_number(matrix_t *A, double number, matrix_t *result) {
  if (s21_check_matrix(A) != OK) return INCORRECT;
  int res = OK;
  if (s21_create_matrix(A->rows, A->columns, result) == OK) {
    for (int i = 0; i < A->rows; i++) {
      for (int j = 0; j < A->columns; j++) {
        result->matrix[i][j] = A->matrix[i][j] * number;
      }
    }
  } else
    res = INCORRECT;
  return res;
}

int s21_mult_matrix(matrix_t *A, matrix_t *B, matrix_t *result) {
  int res = OK;
  if (s21_check_matrix(A) || s21_check_matrix(B))
    res = INCORRECT;
  else if (A->columns != B->rows)
    res = ERROR_CALC;
  else if (s21_create_matrix(A->rows, B->columns, result) != OK)
    res = INCORRECT;
  else
    for (int i = 0; i < result->rows; i++)
      for (int j = 0; j < result->columns; j++)
        for (int k = 0; k < A->columns; k++)
          result->matrix[i][j] += A->matrix[i][k] * B->matrix[k][j];
  return res;
}

int s21_transpose(matrix_t *A, matrix_t *result) {
  int res = OK;
  if (s21_check_matrix(A))
    res = INCORRECT;
  else if (s21_create_matrix(A->columns, A->rows, result) != OK)
    res = INCORRECT;
  else
    for (int i = 0; i < A->rows; i++)
      for (int j = 0; j < A->columns; j++)
        result->matrix[j][i] = A->matrix[i][j];
  return res;
}

int s21_calc_complements(matrix_t *A, matrix_t *result) {
  if (s21_check_matrix(A) || result == NULL) return INCORRECT;
  int res = OK;
  if (A->columns != A->rows) {
    res = ERROR_CALC;
  } else if (A->rows == 1) {
    if (s21_create_matrix(A->rows, A->columns, result) == OK) {
      result->matrix[0][0] = A->matrix[0][0];
    }
  } else if (s21_create_matrix(A->rows, A->columns, result) == OK) {
    for (int i = 0; i < A->rows; i++) {
      for (int j = 0; j < A->columns; j++) {
        matrix_t makeshift = {0};
        double te = 0;
        s21_minor(A, &makeshift, i, j);
        if (s21_determinant(&makeshift, &te) == OK) {
          result->matrix[i][j] = pow(-1, (i + j)) * te;
          s21_remove_matrix(&makeshift);
        } else
          res = ERROR_CALC;
      }
    }
  }
  return res;
}

int s21_determinant(matrix_t *A, double *result) {
  if (s21_check_matrix(A) || result == NULL) return INCORRECT;
  int res = OK;
  *result = 0;
  if (A->columns != A->rows)
    res = ERROR_CALC;
  else if (A->rows == 1)
    *result = A->matrix[0][0];
  else if (A->rows == 2)
    *result =
        A->matrix[0][0] * A->matrix[1][1] - A->matrix[0][1] * A->matrix[1][0];
  else {
    int sign = 1;
    for (int i = 0; i < A->columns; i++) {
      matrix_t makeshift = {0};
      int check_minor = s21_minor(A, &makeshift, 0, i);
      double makeshift_d = 0;
      if (s21_determinant(&makeshift, &makeshift_d) == OK &&
          check_minor == OK) {
        *result += A->matrix[0][i] * makeshift_d * sign;
        sign = -sign;
        s21_remove_matrix(&makeshift);
      } else
        res = ERROR_CALC;
      s21_remove_matrix(&makeshift);
    }
  }
  return res;
}

int s21_inverse_matrix(matrix_t *A, matrix_t *result) {
  int res = OK;
  if (s21_check_matrix(A) != OK)
    res = INCORRECT;
  else if (A->rows != A->columns)
    res = ERROR_CALC;
  matrix_t calc = {0}, transpons = {0};
  double determinant = 0.0;
  s21_determinant(A, &determinant);
  if (fabs(determinant) < 1e-7)
    res = ERROR_CALC;
  else
    s21_calc_complements(A, &calc);
  s21_transpose(&calc, &transpons);
  s21_mult_number(&transpons, 1 / determinant, result);
  s21_remove_matrix(&calc);
  s21_remove_matrix(&transpons);
  return res;
}

// tools
int s21_check_matrix(matrix_t *check) {
  int res = OK;
  if (check == NULL || check->matrix == NULL || check->rows < 1 ||
      check->columns < 1) {
    res = INCORRECT;
  }
  return res;
}

int s21_minor(matrix_t *A, matrix_t *matrix_minor, int index_rows,
              int index_columns) {
  int res = OK;
  if (s21_create_matrix(A->rows - 1, A->columns - 1, matrix_minor) == OK) {
    int minor_row = 0;
    for (int i = 0; i < A->rows; i++) {
      if (i == index_rows) continue;
      int minor_columns = 0;
      for (int j = 0; j < A->columns; j++) {
        if (j == index_columns) continue;
        matrix_minor->matrix[minor_row][minor_columns] = A->matrix[i][j];
        minor_columns++;
      }
      minor_row++;
    }
  } else
    res = INCORRECT;
  return res;
}

void s21_initialization(matrix_t *A) {
  int count = 1;
  for (int i = 0; i < A->rows; i++) {
    for (int j = 0; j < A->columns; j++) {
      A->matrix[i][j] = count++;
    }
  }
}