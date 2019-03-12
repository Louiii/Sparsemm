#include "utils.h"
#include <stdlib.h>
#include <stdio.h>

#include <math.h>
#include <string.h>

#ifdef LIKWID_PERFMON
#include <likwid.h>
#else
#define LIKWID_MARKER_INIT
#define LIKWID_MARKER_THREADINIT
#define LIKWID_MARKER_SWITCH
#define LIKWID_MARKER_REGISTER(regionTag)
#define LIKWID_MARKER_START(regionTag)
#define LIKWID_MARKER_STOP(regionTag)
#define LIKWID_MARKER_CLOSE
#define LIKWID_MARKER_GET(regionTag, nevents, events, time, count)
#endif

void basic_sparsemm(const COO, const COO, COO *);
void basic_sparsemm_sum(const COO, const COO, const COO,
                        const COO, const COO, const COO,
                        COO *);

void initArray(Array *a, size_t initialSize) {
  a->coords = (struct coord *)malloc(initialSize * sizeof(struct coord));
  a->data = (double *)malloc(initialSize * sizeof(double));
  a->used = 0;
  a->size = initialSize;
}

void expandArray(Array *a, int new_size) {
    a->size = new_size;
    a->data = (double *)realloc(a->data, new_size * sizeof(double));
    a->coords = (struct coord *)realloc(a->coords, new_size * sizeof(struct coord));
    printf("HAD TO REALLOC\n");
}

void freeArray(Array *a) {
  free(a->coords);
  free(a->data);
  a->used = a->size = 0;
}

int cmp_y(const void *x, const void *y){
    int (*a)[3] = (int(*)[3])x;
    int (*b)[3] = (int(*)[3])y;
    if((*a)[1]==(*b)[1])
        return (*a)[0] - (*b)[0];
    else
        return (*a)[1] - (*b)[1];
}

int cmp_x(const void *x, const void *y){
    int (*a)[3] = (int(*)[3])x;
    int (*b)[3] = (int(*)[3])y;
    if((*a)[0]==(*b)[0])
        return (*a)[1] - (*b)[1];
    else
        return (*a)[0] - (*b)[0];
}

void convert_COO_to_CSR(const COO sparse, CSR *A_csr) {
    int *JA = malloc(sparse->NZ * sizeof(int));
    int *IA = malloc((1 + sparse->m) * sizeof(int));
    double *data = malloc(sparse->NZ * sizeof(double));

    int coords_data[sparse->NZ][3];
    int idx;
    for (idx = 0; idx < sparse->NZ; idx++) {
        coords_data[idx][0] = sparse->coords[idx].i;
        coords_data[idx][1] = sparse->coords[idx].j;
        coords_data[idx][2] = idx;
    }
    qsort(coords_data, sparse->NZ, sizeof(coords_data[0]), cmp_x);


    int current_col = -1;
    int acc = 0;
    int tot = 0;
    int x, nz, rc;
    for (nz = 0; nz < sparse->NZ; nz++) {//(  (x, y) in coords) {
        JA[nz] = coords_data[nz][1];
        x = coords_data[nz][0];
        data[nz] = sparse->data[coords_data[nz][2]];//index in data
        if (current_col != x) {//we are now on a new row. it may be greater than 1 row down.
            tot += acc;
            IA[current_col + 1] = tot;//MAY NOT BE CURRENT ROW
            acc = 0;
            if (x - current_col > 1) {// there are rows that have 0 elems on: current_col+1,...,x-1
                for (rc = current_col + 1; rc < x; rc++) {//i in range(current_col+1, x) {
                    IA[rc + 1] = tot;
                }
            }
        }
        current_col = x;
        acc ++;
    }
    IA[current_col + 1] = tot + acc;

    CSR sr = calloc(1, sizeof(struct _p_CSR));
    sr->m = sparse->m;
    sr->n = sparse->n;
    sr->NZ = sparse->NZ;
    sr->data = data;
    sr->IA = IA;
    sr->JA = JA;
    *A_csr = sr;
}

void convert_COO_to_CSC(const COO sparse, CSR *B_csc){//, int **coords_data) {
    CSR sr = malloc(sizeof(struct _p_CSR));
    sr->m = sparse->m;
    sr->n = sparse->n;
    sr->NZ = sparse->NZ;

    int coords_data[sparse->NZ][3];
    int idx;
    for (idx = 0; idx < sparse->NZ; idx++) {
        coords_data[idx][0] = sparse->coords[idx].i;
        coords_data[idx][1] = sparse->coords[idx].j;
        coords_data[idx][2] = idx;
    }
    qsort(coords_data, sparse->NZ, sizeof(coords_data[0]), cmp_y);

    int *JA = malloc(sparse->NZ * sizeof(int));
    int *IA = malloc((1 + sparse->n) * sizeof(int));
    double *data = malloc(sparse->NZ * sizeof(double));

    int current_col = -1;
    int acc = 0;
    int tot = 0;
    int x, y, nz, rc, data_idx;
    for (nz = 0; nz < sparse->NZ; nz++) {//(  (x, y) in coords) {
        JA[nz] = coords_data[nz][0];
        y = coords_data[nz][1];
        data[nz] = sparse->data[coords_data[nz][2]];//index in data
        if (current_col != y) {//we are now on a new row. it may be greater than 1 row down.
            tot += acc;
            IA[current_col + 1] = tot;//MAY NOT BE CURRENT ROW
            acc = 0;
            if (y - current_col > 1) {// there are rows that have 0 elems on: current_col+1,...,x-1
                for (rc = current_col + 1; rc < y; rc++) {//i in range(current_col+1, x) {
                    IA[rc + 1] = tot;
                }
            }
        }
        current_col = y;
        acc ++;
    }
    IA[current_col + 1] = tot + acc;

    sr->data = data;
    sr->IA = IA;
    sr->JA = JA;

    *B_csc = sr;
}

void optimised_sparsemm(const COO A, const COO B, COO *C) {
    // LIKWID_MARKER_START("optimised_profile");
	// A -> CSR
    LIKWID_MARKER_START("A_to_CSR");
	CSR A_csr;// = calloc(1, sizeof(struct _p_CSR));
    convert_COO_to_CSR(A, &A_csr);
    LIKWID_MARKER_STOP("A_to_CSR");
    // // B -> CSC
    LIKWID_MARKER_START("B_to_CSB");
    CSR B_csr;
    convert_COO_to_CSC(B, &B_csr);//, coords_data);//
    // LIKWID_MARKER_STOP("B_to_CSB");
	//m - number of rows
	// (mA x nA) * (mB x nB )
	// nA = mB
    LIKWID_MARKER_START("A_times_B");
    int A_m = A->m;
    int B_n = B->n;

    int improved_predict = ceil( pow((A_m * B_n * 0.0067),0.5) );//experimentally found
    
    Array a;
    int new_size;
    initArray(&a, improved_predict);

	int arow, bcol, i, j;
	float total;
	int nz = 0;


    // parallelize this loop;
    // if there are p GPUs, do:
    // split B_csr->IA into p parts (v = 0:p-1)
    // record each col in part of B (JA[IA[v]:IA[v+p]])
    // ----- send to each GPU: -----
    //      loops over section, add to its own array
    //      send back the pointer to the array and the number of elements
    // merge all arrays into one and quicksort.

    #pragma GCC ivdep
    
    #pragma acc data copyin(A, B) copyout(C)
	
    for ( bcol = 0; bcol < B_n; bcol++ ) {//for col in B;
        //make zero vector length of column
        
        #pragma acc kernels loop gang
		
        for (arow = 0; arow < A_m; arow++) {// for row in A:
            
            #pragma acc worker
            
            total = 0;
            i = A_csr->IA[arow], j = B_csr->IA[bcol];// i is the col index in a row of A, j is the row index of a col in B
            
            #pragma acc vector reduction (+:sum)

            while (i < A_csr->IA[arow+1] && j < B_csr->IA[bcol+1]) { // m
				if (A_csr->JA[i] < B_csr->JA[j])
					i++;
				else if (B_csr->JA[j] < A_csr->JA[i])
					j++;
				else {/* if arr1[i] == arr2[j] */
					total += (A_csr->data[i])*(B_csr->data[j]);
					j++;
					i++;
				}
			}
            //find out how long it spends in this if statement
			if (total != 0) {
                if (a.used == a.size) {
                    new_size = a.size * ceil( B_n/(bcol+1) + 0.05 );
                    printf("new size = %d ", new_size);
                    expandArray(&a, new_size);
                }
                a.coords[a.used].i = arow;
                a.coords[a.used].j = bcol;
                a.data[a.used] = total;
                a.used++;
                nz++; 
			}
		}
	}
    alloc_sparse(A_m, B_n, nz, C);
    (*C)->coords = a.coords;
    (*C)->data = a.data;

	free(A_csr);
	free(B_csr);
    LIKWID_MARKER_STOP("A_times_B");
    // LIKWID_MARKER_STOP("optimised_profile");
}

void convert_3COO_to_1CSR(const COO sparseA, const COO sparseB, const COO sparseC, CSR *A2_csr) {
    CSR sr = (struct _p_CSR *)malloc(sizeof(struct _p_CSR));
    sr->m = sparseA->m;
    sr->n = sparseA->n;
    int predNZ = sparseA->NZ + sparseB->NZ + sparseC->NZ;
    // printf("NZ - A, B, C = %d, %d, %d\n", sparseA->NZ, sparseB->NZ, sparseC->NZ);
    // printf("predNZ %d\n", predNZ);
    // copy the y coords , JA
    //make IA
    int *JA = (int *)malloc(predNZ * sizeof(int));
    int *IA = (int *)malloc((1 + sparseA->m) * sizeof(int));
    double *data = (double *)calloc(predNZ , sizeof(double));

    int A_NZ = sparseA->NZ;
    int B_NZ = sparseB->NZ;
    int C_NZ = sparseC->NZ;

    int idx;
    int coords_dataA[A_NZ][3];
    for (idx = 0; idx < A_NZ; idx++) {
        coords_dataA[idx][0] = sparseA->coords[idx].i;
        coords_dataA[idx][1] = sparseA->coords[idx].j;
        coords_dataA[idx][2] = idx;
    }
    qsort(coords_dataA, A_NZ, sizeof(coords_dataA[0]), cmp_x);

    int coords_dataB[B_NZ][3];
    for (idx = 0; idx < B_NZ; idx++) {
        coords_dataB[idx][0] = sparseB->coords[idx].i;
        coords_dataB[idx][1] = sparseB->coords[idx].j;
        coords_dataB[idx][2] = idx;
    }
    qsort(coords_dataB, B_NZ, sizeof(coords_dataB[0]), cmp_x);

    int coords_dataC[C_NZ][3];
    for (idx = 0; idx < C_NZ; idx++) {
        coords_dataC[idx][0] = sparseC->coords[idx].i;
        coords_dataC[idx][1] = sparseC->coords[idx].j;
        coords_dataC[idx][2] = idx;
    }
    qsort(coords_dataC, C_NZ, sizeof(coords_dataC[0]), cmp_x);

    // printf("A_NZ %d %d %d\n", A_NZ, B_NZ, C_NZ );
    int current_row = 0;
    // int no_options;
    int a_idx = 0, b_idx = 0, c_idx = 0;
    int acc = 0;
    int tot = 0;
    int x, rc, i, min_y, min_idx;
    // double d;
    int nz = 0;
    // for (nz = 0; nz < sparse->NZ; nz++) {//(  (x, y) in coords) {
    int options[3];
    // printf("sparseA->n %d\n", sparseA->n);
    // int v;
    // for (v=0;v<A_NZ;v++){printf("row = %d, JA[%d] = %d\n", sparseA->coords[v].i, v, sparseA->coords[v].j);}
    // printf("[a_idx].i, b, c = %d, %d, %d\n", sparseA->coords[a_idx].i, sparseB->coords[b_idx].i, sparseC->coords[c_idx].i);
    while (a_idx < A_NZ || b_idx < B_NZ || c_idx < C_NZ) {
    	// printf("idx: a = %d/%d, b = %d/%d, c = %d/%d : row: a = %d/%d,  b = %d/%d,  c = %d/%d\n", a_idx, A_NZ, b_idx, B_NZ, c_idx, C_NZ, sparseA->coords[a_idx].i, sparseA->m, sparseB->coords[b_idx].i, sparseB->m, sparseC->coords[c_idx].i, sparseC->m);
    	options[0] = 0, options[1] = 0, options[2] = 0;
    	if (a_idx < A_NZ && coords_dataA[a_idx][0] == current_row) { options[0] = 1; }
    	if (b_idx < B_NZ && coords_dataB[b_idx][0] == current_row) { options[1] = 1; }
    	if (c_idx < C_NZ && coords_dataC[c_idx][0] == current_row) { options[2] = 1; }

        // if (a_idx >= A_NZ || A_NZ == 0) { options[0] = 0; }
        // if (b_idx >= B_NZ || B_NZ == 0) { options[1] = 0; }
        // if (c_idx >= C_NZ || C_NZ == 0) { options[2] = 0; }

    	if (coords_dataA[a_idx][0] != current_row && coords_dataB[b_idx][0] != current_row && coords_dataC[c_idx][0] != current_row) {
    		//INCREMENT LEVEL!
    		// printf("---------------increment row from %d, to ", current_row);
    		// printf("a_idx = %d, b_idx = %d, c_idx = %d\n", a_idx, b_idx, c_idx);
    		tot += acc;
            IA[current_row + 1] = tot;
            // printf("[a_idx].i, b, c = %d, %d, %d\n", sparseA->coords[a_idx].i, sparseB->coords[b_idx].i, sparseC->coords[c_idx].i);
    		// if (current_row == sparseA->m) {break;}
            // prev_row = current_row;
            // int v;
            // for (v=0;v<current_row;v++) {printf("%d\n", IA[v]);}

            if (coords_dataA[a_idx][0] == current_row + 1 || coords_dataB[b_idx][0] == current_row + 1 || coords_dataC[c_idx][0] == current_row + 1) {
                current_row++;
            } else { 
                current_row++;
                IA[current_row + 1] = tot;
                printf("\n\nmissing cols need to be added to IA\n\n\n");
        
                // x = sparseA->coords[a_idx].i;
                // if (sparseB->coords[b_idx].i < x) { x = sparseB->coords[b_idx].i; }
                // if (sparseC->coords[c_idx].i < x) { x = sparseC->coords[c_idx].i; }
            
                // if (x - current_row > 1) {// there are rows that have 0 elems on: current_row+1,...,x-1
                //     int rc;
                //     for (rc = current_row + 1; rc < x; rc++) {//i in range(current_row+1, x) {
                //         IA[rc + 1] = tot;
                //     }
                // }
                // current_row = x;
                // printf("\n\nmissing rows need to be added to IA\n\n\n");
            }

      //       current_row = sparseA->coords[a_idx].i;
    		// if (sparseB->coords[b_idx].i < current_row) { current_row = sparseB->coords[b_idx].i; }
    		// if (sparseC->coords[c_idx].i < current_row) { current_row = sparseC->coords[c_idx].i; }
    		// printf("%d, of %d\n", current_row, sparseA->m);
    		acc = 0;
    	} else {
            // if (options[0] == 0 && options[1] == 0 && options[2] == 0) {
            //     // IA[current_row + 1] = 0;
            //     // current_row++;
            //     printf("\nfuaacccj\n\n");
            //     break;
            // } else {
                // printf("options:");
    	    	for (i = 0; i<3; i++){
                    // printf("%d,", options[i]);
    	    		if (options[i] == 1) {
    	    			if (i == 0) { options[0] = coords_dataA[a_idx][1]; }
    	    			if (i == 1) { options[1] = coords_dataB[b_idx][1]; }
    	    			if (i == 2) { options[2] = coords_dataC[c_idx][1]; }
    	    		} else {
    	    			options[i] = sparseA->n;
    	    		}
    	    	}
    	    	min_y = options[0];
    	    	min_idx = 0;
    	    	if (options[1] < min_y){ min_y = options[1]; min_idx = 1;}
    	    	if (options[2] < min_y){ min_y = options[2]; min_idx = 2;}
    	    	// increment appropriate index and place in appropriate data.
    	    	//IF MIN IS REPEATED ADD TOGETHER AND INCREMENT TOGTHER
                // printf(" > min y = %d > ", min_y);
                // printf("options:");
    	    	for (i=0;i<3;i++){
                    // printf("%d,", options[i]);
    	    		if (options[i] == min_y) {
    			    	if (i==0) { data[nz] += sparseA->data[coords_dataA[a_idx][2]]; a_idx++; }
    			    	if (i==1) { data[nz] += sparseB->data[coords_dataB[b_idx][2]]; b_idx++; }
    			    	if (i==2) { data[nz] += sparseC->data[coords_dataC[c_idx][2]]; c_idx++; }
    			    }
    			}
                // printf(" a_idx = %d \n",a_idx);
    		    JA[nz] = min_y;
    		    nz++;
    		    acc++;
            // }
		}
    }
    // printf("Actual NZ = %d\n", nz);
    IA[current_row + 1] = tot + acc;
    IA[0]=0;
    // printf("current_row %d\n", current_row);
    // data = (double *)realloc(data, nz * sizeof(double));
    // JA = (int *)realloc(JA, nz * sizeof(int));

    sr->NZ = nz;
    sr->data = data;
    sr->IA = IA;
    sr->JA = JA;

    *A2_csr = sr;
}

void convert_3COO_to_1CSC(const COO sparseD, const COO sparseE, const COO sparseF, CSR *B2_csr) {
    CSR sr = (struct _p_CSR *)malloc(sizeof(struct _p_CSR));
    sr->m = sparseD->m;
    sr->n = sparseD->n;

    int predNZ = sparseD->NZ + sparseE->NZ + sparseF->NZ;

    int *JA = (int *)malloc(predNZ * sizeof(int));
    int *IA = (int *)malloc((1 + sparseD->n) * sizeof(int));
    double *data = (double *)calloc(predNZ, sizeof(double));

    int D_NZ = sparseD->NZ;
    int E_NZ = sparseE->NZ;
    int F_NZ = sparseF->NZ;


    // MAKE TRANSFORM DICTS FOR D,E,F
    int idx;
    int coords_dataD[D_NZ][3];
    for (idx = 0; idx < D_NZ; idx++) {
    	coords_dataD[idx][0] = sparseD->coords[idx].i;
    	coords_dataD[idx][1] = sparseD->coords[idx].j;
    	coords_dataD[idx][2] = idx;
    }
    qsort(coords_dataD, D_NZ, sizeof(coords_dataD[0]), cmp_y);

    int coords_dataE[E_NZ][3];
    for (idx = 0; idx < E_NZ; idx++) {
    	coords_dataE[idx][0] = sparseE->coords[idx].i;
    	coords_dataE[idx][1] = sparseE->coords[idx].j;
    	coords_dataE[idx][2] = idx;
    }
    qsort(coords_dataE, E_NZ, sizeof(coords_dataE[0]), cmp_y);

    int coords_dataF[F_NZ][3];
    for (idx = 0; idx < F_NZ; idx++) {
    	coords_dataF[idx][0] = sparseF->coords[idx].i;
    	coords_dataF[idx][1] = sparseF->coords[idx].j;
    	coords_dataF[idx][2] = idx;
    }
    qsort(coords_dataF, F_NZ, sizeof(coords_dataF[0]), cmp_y);


    int current_col = 0;
    int d_idx = 0, e_idx = 0, f_idx = 0;
    int acc = 0;
    int tot = 0;
    int x, rc, i, min_x, min_idx;
    // double d;
    int nz = 0;
    // for (nz = 0; nz < sparse->NZ; nz++) {//(  (x, y) in coords) {
    int options[3];
    // printf("NZ : %d %d %d\n", D_NZ, E_NZ, F_NZ);
    int v;
    // printf("[d_idx].i, b, c = %d, %d, %d\n", sparseD->coords[d_idx].i, sparseD->coords[e_idx].i, sparseF->coords[f_idx].i);
    while (d_idx < D_NZ || e_idx < E_NZ || f_idx < F_NZ) {
    	// printf("d_idx = %d, e_idx = %d, f_idx = %d\n", d_idx, e_idx, f_idx);
    	options[0] = 0, options[1] = 0, options[2] = 0;
    	if (d_idx < D_NZ && coords_dataD[d_idx][1] == current_col) { options[0] = 1; }
    	if (e_idx < E_NZ && coords_dataE[e_idx][1] == current_col) { options[1] = 1; }
    	if (f_idx < F_NZ && coords_dataF[f_idx][1] == current_col) { options[2] = 1; }
        // printf("current col: %d, acc: %d, total: %d, didx=%d, eidx=%d, fidx=%d\n",current_col, acc, tot, d_idx, e_idx, f_idx);
    	if (coords_dataD[d_idx][1] != current_col && coords_dataE[e_idx][1] != current_col && coords_dataF[f_idx][1] != current_col) {
    		// printf("cols: %d, %d, %d\n", coords_dataD[d_idx][1], coords_dataE[e_idx][1], coords_dataF[f_idx][1]);
            //INCREMENT LEVEL!
    		tot += acc;
            IA[current_col + 1] = tot;
            // for (v=0;v<current_col;v++){printf("IA[%d] = %d\n", v, IA[v]);}
            // printf("\n");
            if (coords_dataD[d_idx][1] == current_col + 1 || coords_dataE[e_idx][1] == current_col + 1 || coords_dataF[f_idx][1] == current_col + 1) {
                current_col++;
            } else { 
                current_col++;
                IA[current_col + 1] = tot;
                // printf("\n\nmissing cols need to be added to IA\n\n\n");
            }
    		acc = 0;
    	} else {
            // printf("options:");
            // for (v=0;v<3;v++){printf("%d, ", options[v]);}
            // printf("\n");
	    	for (i = 0; i<3; i++){
	    		if (options[i] == 1) {
	    			if (i == 0) { options[0] = coords_dataD[d_idx][0]; }
	    			if (i == 1) { options[1] = coords_dataE[e_idx][0]; }
	    			if (i == 2) { options[2] = coords_dataF[f_idx][0]; }
	    		} else {
	    			options[i] = sparseD->m;
	    		}
	    	}
            // printf("options:");
            // for (v=0;v<3;v++){printf("%d, ", options[v]);}
            // printf("\n");
	    	min_x = options[0];
	    	min_idx = 0;
	    	if (options[1] < min_x){ min_x = options[1]; min_idx = 1;}
	    	if (options[2] < min_x){ min_x = options[2]; min_idx = 2; }
	    	// increment appropriate index and place in appropriate data.
	    	//IF MIN IS REPEATED ADD TOGETHER AND INCREMENT TOGTHER
	    	for (i=0;i<3;i++){
	    		if (options[i] == min_x) {
			    	if (i==0) { data[nz] += sparseD->data[coords_dataD[d_idx][2]]; d_idx++; }
			    	if (i==1) { data[nz] += sparseE->data[coords_dataE[e_idx][2]]; e_idx++; }
			    	if (i==2) { data[nz] += sparseF->data[coords_dataF[f_idx][2]]; f_idx++; }
                    // printf("data:%f\n", data[nz]);
			    }
			}
            // for (v=0;v<nz;v++){printf("%f, ", data[v]);}
		    JA[nz] = min_x;
		    nz++;
		    acc++;
		}
    }
    // printf("Actual NZ = %d\n", nz);
    IA[current_col + 1] = tot + acc;
    IA[0]=0;
    // printf("final_col index %d\n", current_col+1);
    // for (v=0;v<nz;v++){printf("%f, ", data[v]);}
    // data = (double *)realloc(data, nz * sizeof(double));
    // JA = (int *)realloc(JA, nz * sizeof(int));

    sr->NZ = nz;
    sr->data = data;
    sr->IA = IA;
    sr->JA = JA;

    *B2_csr = sr;
}

void optimised_sparsemm2(const CSR *A_csr, const CSR *B_csr, COO *C) {
	//m - number of rows
	// (mA x nA) * (mB x nB )
	// nA = mB
    int A_m = (*A_csr)->m;
    int B_n = (*B_csr)->n;

    int A_NZ = (*A_csr)->NZ;
    int B_NZ = (*B_csr)->NZ;

    int improved_predict = ceil( pow((A_m * B_n * 0.0067),0.5) );//experimentally found

    Array a;
    int new_size;
    initArray(&a, improved_predict);

    int arow, bcol, i, j;
    float total;
    int nz = 0;
	for ( bcol = 0; bcol < B_n; bcol++ ) {//for col in B;
		for (arow = 0; arow < A_m; arow++) {// for row in A:
			total = 0;
			i = (*A_csr)->IA[arow], j = (*B_csr)->IA[bcol];// i is the col index in a row of A, j is the row index of a col in B
			while (i < (*A_csr)->IA[arow+1] && j < (*B_csr)->IA[bcol+1]) { // m
				if ((*A_csr)->JA[i] < (*B_csr)->JA[j])
					i++;
				else if ((*B_csr)->JA[j] < (*A_csr)->JA[i])
					j++;
				else {/* if arr1[i] == arr2[j] */
					total += ((*A_csr)->data[i])*((*B_csr)->data[j]);
					j++;
					i++;
				}
			}
            if (total != 0) {
                if (a.used == a.size) {
                    new_size = a.size * ceil( B_n/(bcol+1) + 0.05 );
                    printf("new size = %d ", new_size);
                    expandArray(&a, new_size);
                }
                a.coords[a.used].i = arow;
                a.coords[a.used].j = bcol;
                a.data[a.used] = total;
                a.used++;
                nz++;
            }
		}
	}
    alloc_sparse(A_m, B_n, nz, C);
    (*C)->coords = a.coords;
    (*C)->data = a.data;

	free(*A_csr);
	free(*B_csr);
}

void optimised_sparsemm_sum(const COO A, const COO B, const COO C,
                            const COO D, const COO E, const COO F,
                            COO *O)
{
    LIKWID_MARKER_START("ABC_to_CSR");
    CSR A2_csr;
    convert_3COO_to_1CSR(A, B, C, &A2_csr);

    LIKWID_MARKER_STOP("ABC_to_CSR");
    LIKWID_MARKER_START("DEF_to_CSC");
    CSR B2_csr;
    convert_3COO_to_1CSC(D, E, F, &B2_csr);
    LIKWID_MARKER_STOP("DEF_to_CSC");

    LIKWID_MARKER_START("A_times_B");
    optimised_sparsemm2(&A2_csr, &B2_csr, O);
    LIKWID_MARKER_STOP("A_times_B");
    printf("NZ = %d, n = %d\n", (*O)->NZ, A->n);
}