#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "utils.h"


// #ifdef LIKWID_PERFMON
// #include <likwid.h>
// #else
// #define LIKWID_MARKER_INIT
// #define LIKWID_MARKER_THREADINIT
// #define LIKWID_MARKER_SWITCH
// #define LIKWID_MARKER_REGISTER(regionTag)
// #define LIKWID_MARKER_START(regionTag)
// #define LIKWID_MARKER_STOP(regionTag)
// #define LIKWID_MARKER_CLOSE
// #define LIKWID_MARKER_GET(regionTag, nevents, events, time, count)
// #endif

/* Compute C = C + A*B in dense, column major, format. */
static void dgemm(int m, int n, int k, const double *a, const double *b, double *c)
{
    int i, j, p;
    int lda = m;
    int ldb = k;
    int ldc = m;
    for (j = 0; j < n; j++) {
        for (p = 0; p < k; p++) {
            for (i = 0; i < m; i++) {
                c[j*ldc + i] = c[j*ldc + i] + a[p*lda + i] * b[j*ldb + p];
            }
        }
    }
}

/* Computes C = A*B by converting A and B to dense column major
 * format, performing the matrix-matrix multiplication, and then
 * converting the result back to sparse.
 * C will be allocated by this routine.
 */
void basic_sparsemm(const COO A, const COO B, COO *C){
    // LIKWID_MARKER_START("basic_profile");

    double *a = NULL;
    double *b = NULL;
    double *c = NULL;
    int m, n, k;
    convert_sparse_to_dense(A, &a);
    convert_sparse_to_dense(B, &b);

    *C = NULL;
    m = A->m;
    k = A->n;
    n = B->n;
    if (k != B->m) {
        fprintf(stderr, "Invalid matrix sizes, got %d x %d and %d x %d\n", A->m, A->n, B->m, B->n);
        free(a);
        free(b);
        exit(1);
    }
    alloc_dense(m, n, &c);
    zero_dense(m, n, c);

    dgemm(m, n, k, a, b, c);
    free_dense(&a);
    free_dense(&b);
    convert_dense_to_sparse(c, m, n, C);
    free_dense(&c);

    // LIKWID_MARKER_STOP("basic_profile");
}

////////////////////////////////////////////////////
// void cnvrt_COO_to_CSR(const COO sparse, CSR *A_csr) {
//     CSR sr = calloc(1, sizeof(struct _p_CSR));
//     sr->m = sparse->m;
//     sr->n = sparse->n;
//     sr->NZ = sparse->NZ;
//     // copy the y coords , JA
//     //make IA
//     int *JA = malloc(sparse->NZ * sizeof(int));
//     int *IA = malloc((1 + sparse->m) * sizeof(int));
//     double *data = malloc(sparse->NZ * sizeof(double));
    
//     int current_row = -1;
//     int acc = 0;
//     int tot = 0;
//     int x, nz, rc;
//     for (nz = 0; nz < sparse->NZ; nz++) {//(  (x, y) in coords) {
//         x = sparse->coords[nz].i;
//         JA[nz] = sparse->coords[nz].j;
//         data[nz] = sparse->data[nz];
//         if (current_row != x) {//we are now on a new row. it may be greater than 1 row down.
//             tot += acc;
//             IA[current_row + 1] = tot;//MAY NOT BE CURRENT ROW
//             acc = 0;
//             if (x - current_row > 1) {// there are rows that have 0 elems on: current_row+1,...,x-1
//                 for (rc = current_row + 1; rc < x; rc++) {//i in range(current_row+1, x) {
//                     IA[rc + 1] = tot;
//                 }
//             }
//         }
//         current_row = x;
//         acc ++;
//     }
//     IA[current_row + 1] = tot + acc;

//     sr->data = data;

//     sr->IA = IA;
//     sr->JA = JA;

//     *A_csr = sr;
// }
//////////////////////////////////////////

/* Computes O = (A + B + C) (D + E + F) by converting to dense column
 * major, performing the matrix matrix multiplication, and converting
 * back to sparse.  This routine allocates O.*/
void basic_sparsemm_sum(const COO A, const COO B, const COO C,
                        const COO D, const COO E, const COO F,
                        COO *O)
{
    double *a = NULL;
    double *b = NULL;
    double *c = NULL;
    double *d = NULL;
    double *e = NULL;
    double *f = NULL;
    double *o = NULL;
    int i, j, m, n, k;

    m = A->m;
    k = A->n;
    n = D->n;
    if (A->m != B->m || A->n != B->n) {
        fprintf(stderr, "A (%d x %d) and B (%d x %d) are not the same shape\n",
                A->m, A->n, B->m, B->n);
        exit(1);
    }
    if (A->m != C->m || A->n != C->n) {
        fprintf(stderr, "A (%d x %d) and C (%d x %d) are not the same shape\n",
                A->m, A->n, C->m, C->n);
        exit(1);
    }
    if (D->m != E->m || D->n != E->n) {
        fprintf(stderr, "D (%d x %d) and E (%d x %d) are not the same shape\n",
                D->m, D->n, E->m, E->n);
        exit(1);
    }
    if (D->m != F->m || D->n != F->n) {
        fprintf(stderr, "D (%d x %d) and F (%d x %d) are not the same shape\n",
                D->m, D->n, F->m, F->n);
        exit(1);
    }

    if (A->n != D->m) {
        fprintf(stderr, "Invalid matrix sizes, got %d x %d and %d x %d\n",
                A->m, A->n, D->m, D->n);
        exit(1);
    }
        
    convert_sparse_to_dense(A, &a);
    convert_sparse_to_dense(B, &b);
    convert_sparse_to_dense(C, &c);
    convert_sparse_to_dense(D, &d);
    convert_sparse_to_dense(E, &e);
    convert_sparse_to_dense(F, &f);

    /* Compute sums */
    for (j = 0; j < k; j++) {
        for (i = 0; i < m; i++) {
            a[j*m + i] += b[j*m + i] + c[j*m + i];
        }
    }
    /////////////////////////////////////////
    //test:
    // int v;
    // COO *tes;
    // convert_dense_to_sparse(a, m, k, tes);
    // CSR *test;
    // cnvrt_COO_to_CSR(*tes, test);
    // printf("\n\nchecking A:\n\n");
    // printf("A_csr = %p\n", (*test));
    // printf("FINAL : A_csr->m = %d, A_csr->n = %d, A_csr->NZ = %d\n", (*test)->m, (*test)->n, (*test)->NZ);

    // for (v=0; v<(*test)->m; v++){
    //     printf("IA[%d] = %d\n", v, (*test)->IA[v]);
    // }
    // for (v=0; v<(*test)->NZ; v++){
    //     printf("JA[%d] = %d data = %f\n", v, (*test)->JA[v], (*test)->data[v]);
    // }
    ////////////////////////////////////////
    for (j = 0; j < n; j++) {
        for (i = 0; i < k; i++) {
            d[j*k + i] += e[j*k + i] + f[j*k + i];
        }
    }
    free_dense(&b);
    free_dense(&c);
    free_dense(&e);
    free_dense(&f);
    alloc_dense(m, n, &c);
    zero_dense(m, n, c);
    dgemm(m, n, k, a, d, c);
    free_dense(&a);
    free_dense(&d);
    convert_dense_to_sparse(c, m, n, O);
    free_dense(&c);
}
