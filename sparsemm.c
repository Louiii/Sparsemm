#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

#include "utils.h"

// #if __STDC_VERSION__ < 199901L
// #define restrict  
// #endif
// #include <immintrin.h>
// #include <time.h>

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


void basic_sparsemm(const COO, const COO, COO*);
void basic_sparsemm_sum(const COO, const COO, const COO,
                            const COO, const COO, const COO,
                            COO *);
void optimised_sparsemm(const COO, const COO, COO*);
void optimised_sparsemm_sum(const COO, const COO, const COO,
                            const COO, const COO, const COO,
                            COO *);

void check_BASIC_OPT(COO *Cbasic, COO *Copt) {
    double *basic, *opt;
    int pass = 0;
    convert_sparse_to_dense(*Cbasic, &basic);
    convert_sparse_to_dense(*Copt, &opt);
    int i,j;
    int m = (*Cbasic)->m;
    int n = (*Cbasic)->n;
    for (j = 0; j < n; j++) {
        for (i = 0; i < m; i++) {
            double diff = fabs(opt[j*m + i] - basic[j*m + i]);
            if (diff != diff || diff > 1e-3) {
                fprintf(stderr, "MM Failed check at entry (%d, %d), basic value %g, opt value %g\n", i, j, basic[j*m + i], opt[j*m + i]);
                pass = 1;
            }
        }
    }
    if (pass == 0){printf("AB SUCCESS!!\n");}else{printf("AB FAIL!!\n");}

    free(basic);
    free(opt);
    free_sparse(&Cbasic);
    free_sparse(&Copt);
}

static int check_sparsemm()
{
    COO A, B, Cbasic, Copt;
    double *basic, *opt;
    int i, j, m, n, k;
    int pass = 0;

    m = 100;//m - number of rows
    k = 220;
    n = 432;
    random_matrix(m, k, 0.01, &A);
    random_matrix(k, n, 0.02, &B);

    basic_sparsemm(A, B, &Cbasic);
    optimised_sparsemm(A, B, &Copt);

    convert_sparse_to_dense(Cbasic, &basic);
    convert_sparse_to_dense(Copt, &opt);

    for (j = 0; j < n; j++) {
        for (i = 0; i < m; i++) {
            double diff = fabs(opt[j*m + i] - basic[j*m + i]);
            if (diff != diff || diff > 1e-3) {
                fprintf(stderr, "sparseMM Failed check at entry (%d, %d), basic value %g, opt value %g\n", i, j, basic[j*m + i], opt[j*m + i]);
                pass = 1;
            }
        }
    }
    if (pass == 0) {
        printf("SUCCESS check_sparsemm\n");
    }

    free(basic);
    free(opt);
    free_sparse(&A);
    free_sparse(&B);
    free_sparse(&Cbasic);
    free_sparse(&Copt);

    return pass;
}

static int check_sparsemm_sum()
{
    COO A, B, C, D, E, F, Obasic, Oopt;
    double *basic, *opt;
    int i, j, m, n, k;
    int pass = 0;

    m = 128;
    k = 301;
    n = 322;
    random_matrix(m, k, 0.01, &A);
    random_matrix(m, k, 0.04, &B);
    random_matrix(m, k, 0.01, &C);
    random_matrix(k, n, 0.02, &D);
    random_matrix(k, n, 0.03, &E);
    random_matrix(k, n, 0.15, &F);
    basic_sparsemm_sum(A, B, C, D, E, F, &Obasic);
    printf("completed basic ssmm\n");
    optimised_sparsemm_sum(A, B, C, D, E, F, &Oopt);

    convert_sparse_to_dense(Obasic, &basic);
    convert_sparse_to_dense(Oopt, &opt);

    for (j = 0; j < n; j++) {
        for (i = 0; i < m; i++) {
            // printf("basic:%f,  opt:%f\n", basic[j*m + i], opt[j*m + i]);
            double diff = fabs(opt[j*m + i] - basic[j*m + i]);
            if (diff != diff || diff > 1e-3) {
                fprintf(stderr, "sum sparseMM-SUM Failed check at entry (%d, %d), basic value %g, opt value %g\n", i, j, basic[j*m + i], opt[j*m + i]);
                pass = 1;
            }
        }
    }
    if (pass == 0) {
        printf("SUCCESS check_sparsemm_sum\n");
    }
    

    free(basic);
    free(opt);
    free_sparse(&A);
    free_sparse(&B);
    free_sparse(&C);
    free_sparse(&D);
    free_sparse(&E);
    free_sparse(&F);
    free_sparse(&Obasic);
    free_sparse(&Oopt);

    return pass;
}

int main(int argc, char **argv)
{
    COO O;
    FILE *f;

    // LIKWID_MARKER_INIT;
    // LIKWID_MARKER_THREADINIT;

    void (*reader)(const char *, COO *) = &read_sparse;
    void (*writer)(FILE *, COO) = &write_sparse;
    if (!(argc == 2 || argc == 4 || argc == 8 || argc == 5 || argc == 9)) {
        fprintf(stderr, "Invalid arguments.\n");
        fprintf(stderr, "Usage: %s CHECK\n", argv[0]);
        fprintf(stderr, "  Check the implemented routines using randomly generated matrices.\n");
        fprintf(stderr, "Alternate usage: %s [--binary] O A B\n", argv[0]);
        fprintf(stderr, "  Computes O = A B\n");
        fprintf(stderr, "  Where A and B are filenames of matrices to read.\n");
        fprintf(stderr, "  O is the filename of the matrix to be written.\n\n");
        fprintf(stderr, "Alternate usage: %s [--binary] O A B C D E F\n", argv[0]);
        fprintf(stderr, "  Computes O = (A + B + C) (D + E + F)\n");
        fprintf(stderr, "  Where A-F are the files names of matrices to read.\n");
        fprintf(stderr, "  O is the filename of the matrix to be written.\n\n");
        fprintf(stderr, "If the --binary flag is given use binary reading and writing of matrices.\n\n");
        return 1;
    }

    if (argc == 2) {
        int pass = 0;
        if (strcmp(argv[1], "CHECK")) {
            fprintf(stderr, "Invalid mode, expecting CHECK, got %s\n", argv[1]);
            return 1;
        }
        // a |= b, assigns a to (a OR b)
        pass |= check_sparsemm();// a = 0 || (0/1), if check_sparsemm() == 1 it failed.
        pass |= check_sparsemm_sum();
        if (pass == 0){
            printf("SUCCESS!\n");
        } else {
            printf("FAILED!\n");
        }
        return pass;
    } 
    if (argc == 5 || argc == 9) {
      if (strcmp(argv[1], "--binary") != 0) {
        fprintf(stderr, "Expecting flag --binary, not '%s'\n", argv[1]);
        return 1;
      }
      reader = &read_sparse_binary;
      writer = &write_sparse_binary;
      argc--;
      argv++;
    }
    if (argc == 4) {
        COO A, B;
        read_sparse(argv[2], &A);
        read_sparse(argv[3], &B);

        
        // basic_sparsemm(A, B, &O);
        // optimised_sparsemm(A, B, &O);
        //////////MY CHECK!/////////
        COO O1;
        COO O2;
        printf("opt starting...\n");
        optimised_sparsemm(A, B, &O2);
        printf("opt complete.\nbasic starting...\n");
        basic_sparsemm(A, B, &O1);
        printf("basic complete.\n");
        check_BASIC_OPT(&O1, &O2);
        ////////////////////////////
        printf("NON-ZEROS : %d, n = %d, m = %d\n", O->NZ, O->n, O->m);
        
        free_sparse(&A);
        free_sparse(&B);
    } else {
        COO A, B, C, D, E, F;
        read_sparse(argv[2], &A);
        read_sparse(argv[3], &B);
        read_sparse(argv[4], &C);
        read_sparse(argv[5], &D);
        read_sparse(argv[6], &E);
        read_sparse(argv[7], &F);

        optimised_sparsemm_sum(A, B, C, D, E, F, &O);

        free_sparse(&A);
        free_sparse(&B);
        free_sparse(&C);
        free_sparse(&D);
        free_sparse(&E);
        free_sparse(&F);
    }
    // LIKWID_MARKER_CLOSE;

    f = fopen(argv[1], "w");
    write_sparse(f, O);
    free_sparse(&O);
    fclose(f);
    return 0;
}
