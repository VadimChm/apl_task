//check if input matrix is orthogonal
//gcc ort_mat.c | ./a.out <n - dim matrix> <file, if needed>
#include<stdio.h>
#include <stdlib.h>
#include <math.h>
#define MAX_PRINT 10
#define EPS 1e-14

void print_matrix (double *mat, int n, FILE *f);
double ij_func (int i, int j);
void fill_matrix (double *mat, int n);
int read_matrix (double *mat, int n, char *file_name);
double diff_with_id_matrix (double *mat, int n);
double scale_mult_ij_row (double *mat, int n, int I, int J);
double mTm_diff_id (double *mat, int n);


void print_matrix (double *mat, int n, FILE *f)
{
  int i = 0;
  int j = 0;
  if (n <= 0)
    return;
  for (i = 0; i < n; i++)
    {
      for (j = 0; j < n; j++)
        {
          fprintf (f, "%f ", mat[i * n + j]);
          if (j >= (MAX_PRINT - 1))
            break;
        }
      fprintf (f, "\n");
      if (i >= (MAX_PRINT - 1))
        break;
    }
}

double ij_func (int i, int j)
{
  if ((i == 0) && (j == 0))
    return 1. / sqrt (2);
  if ((i == 0) && (j == 2))
    return -1. / sqrt (2);
  if ((i == 2) && (j == 0))
    return 1. / sqrt (2);
  if ((i == 2) && (j == 2))
    return 1. / sqrt (2);
  if (i == j)
    return 1;

  return 0;
}

void fill_matrix (double *mat, int n)
{
  int i = 0;
  int j = 0;
  for (i = 0; i < n; i++)
    {
      for (j = 0; j < n; j++)
        mat[i * n + j] = ij_func (i, j);
    }
}

int read_matrix (double *mat, int n, char *file_name)
{
  int i = 0;
  int j = 0;
  if (!mat || !file_name || n < 0)
    return -1;

  FILE *f = fopen (file_name, "r");
  if (!f)
    {
      fprintf (stderr, "Can not open file \"%s\"!\n", file_name);
      return -2;
    }

  for (i = 0; i < n; i++)
    {   
      for (j = 0; j < n; j++)
        {
          if (fscanf (f, "%lf", mat + (i * n + j)) != 1)
            {
              fprintf (stderr, "Wrong data in file \"%s\"!\n", file_name);
              fclose (f);
              return -3;
            }
        }
    }

  fclose (f);
  return 0;
}

double diff_with_id_matrix (double *mat, int n)
{
  double res = 0;
  int i = 0;
  int j = 0;
  if (!mat)
    return 0;

  for (i = 0; i < n; i++)
    {
      for (int j = 0; j < n; j++)
        {
          if (i == j)
            res += fabs (mat[i * n + j] - 1.);
          else 
            res += fabs (mat[i * n + j]);
        }
    }

  return res;
}

double scale_mult_ij_row (double *mat, int n, int I, int J)
{
  int i = 0;
  double res = 0;
  if (!mat || I < 0 || I > n || J < 0 || J > n)
    return 0;

  for (i = 0; i < n; i++)
    res += mat[I * n + i] * mat[J * n + i];

  return res;
}

double mTm_diff_id (double *mat, int n)
{
  int i = 0;
  int j = 0;
  double res = 0;
  double mTm_ij = 0;
  if (!mat)
    return 0;

  for (i = 0; i < n; i++)
    {
      for (j = 0; j < n; j++)
        {
          mTm_ij = scale_mult_ij_row (mat, n, i, j);

          if (i == j)
            res += fabs (mTm_ij - 1);
          else
            res += fabs (mTm_ij);
        }
    }

  return res;
}

int main(int argc, char *argv[])
{
  int n = 0;
  char *file_name = 0;
  double *matr = 0;
  double res = 0;
  if (((argc != 2) && (argc != 3)) || 
      ((n = atoi (argv[1])) <= 0)) 
    {
      printf("Ошибка ввода...\n");
      printf("[размерность матрицы] [имя файла (если нужно)]\n");
      return 1;
    }
  if (argc == 3)
    file_name = argv[2];

  matr=(double*) malloc (n * n * sizeof(double));
  if (!matr)
    {
      fprintf (stderr, "Not enough memory!\n");
      return 2;
    }

  if (argc == 3)
    {
      if (read_matrix (matr, n, file_name) < 0)
        {
          free (matr);
          return 3;
        }
    }
  else 
    fill_matrix (matr, n);

  printf ("Input matrix:\n");
  print_matrix (matr, n, stdout);

  res = mTm_diff_id (matr, n);
  if (res > EPS)
    printf ("\n-> False\n");
  else 
    printf ("\n-> True\n");

  free (matr);
  return 0;
}