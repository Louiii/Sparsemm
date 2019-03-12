#ifndef _UTILS_H
#define _UTILS_H

#include <stdio.h>

struct coord {
    int i, j;
};

struct _p_COO {
    int m, n, NZ;
    struct coord *coords;
    double *data;
};

typedef struct _p_COO *COO;


struct _p_CSR {
    int m, n, NZ;
    double *data;
    int *IA;
    int *JA;
};

typedef struct _p_CSR *CSR;

// typedef struct {
//   double *data;
//   struct coord *coords;
//   size_t used;
//   size_t size;
// } Array;
typedef struct {
  double *data;
  struct coord *coords;
  size_t used;
  size_t size;
  // struct Array* next;
} Array;

// struct _p_DYN {
//     int size, act_size;
//     double *data;
//     struct coord *coords;
//     struct _p_DYN *next;
// };

// typedef struct _p_DYN *DYN;

void alloc_sparse(int, int, int, COO*);
void free_sparse(COO*);
// void alloc_csr(int, int, int, int, CSR*);
// void free_csr(CSR*);
void alloc_dense(int, int, double **);
void free_dense(double **);
void zero_dense(int, int, double *);

void convert_sparse_to_dense(const COO, double **);
void convert_dense_to_sparse(const double *, int, int, COO *);

void read_sparse(const char *, COO *);
void write_sparse(FILE *, COO);
void read_sparse_binary(const char *, COO *);
void write_sparse_binary(FILE *, COO);
void print_sparse(COO);
void random_matrix(int, int, double, COO *);

#endif
