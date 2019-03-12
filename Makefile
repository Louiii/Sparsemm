CFLAGS = -O3 -march=native  -D_GNU_SOURCE
LDFLAGS = -lm
CC = gcc
PROF = -llikwid -pthread -DLIKWID_PERFMON -I/ddn/apps/Cluster-Apps/likwid/4.1/include -L/ddn/apps/Cluster-Apps/likwid/4.1/lib
VEC = -fopt-info-vec

MFL = small_matrices
SZ = small
TYP = 2D
OUTDIR = outop

OBJ = optimised-sparsemm.o basic-sparsemm.o utils.o
HEADER = utils.h

.PHONY: clean help check

all: sparsemm

help:
	@echo "Available targets are"
	@echo "  clean: Remove all build artifacts"
	@echo "  check: Perform a simple test of your optimised routines"
	@echo "  sparsemm: Build the sparse matrix-matrix multiplication binary"
	@echo "  profile all matrices."
	@echo "  -- edit vars in vim to change between small and large (TYP 2D -> 3D) and out dir."
	@echo "  -- to change basic to opt, edit sparsemm.c"

clean:
	-rm -f sparsemm $(OBJ)

check: sparsemm
	./sparsemm CHECK

profile:
	sbatch -o $(OUTDIR)/dg1-mass-mass sparse.slurm $(MFL)/DG1-mass-$(TYP) $(MFL)/DG1-mass-$(TYP)
	sbatch -o $(OUTDIR)/dg1-ip-mass sparse.slurm $(MFL)/DG1-ip-laplace-$(TYP) $(MFL)/DG1-mass-$(TYP)
	sbatch -o $(OUTDIR)/dg1-mass-ip sparse.slurm $(MFL)/DG1-mass-$(TYP) $(MFL)/DG1-ip-laplace-$(TYP)
	sbatch -o $(OUTDIR)/dg1-ip-ip sparse.slurm $(MFL)/DG1-ip-laplace-$(TYP) $(MFL)/DG1-ip-laplace-$(TYP)
	sbatch -o $(OUTDIR)/dg2-mass-mass sparse.slurm $(MFL)/DG2-mass-$(TYP) $(MFL)/DG2-mass-$(TYP)
	sbatch -o $(OUTDIR)/dg2-ip-mass sparse.slurm $(MFL)/DG2-ip-laplace-$(TYP) $(MFL)/DG2-mass-$(TYP)
	sbatch -o $(OUTDIR)/dg2-mass-ip sparse.slurm $(MFL)/DG2-mass-$(TYP) $(MFL)/DG2-ip-laplace-$(TYP)
	sbatch -o $(OUTDIR)/dg2-ip-ip sparse.slurm $(MFL)/DG2-ip-laplace-$(TYP) $(MFL)/DG2-ip-laplace-$(TYP)
	sbatch -o $(OUTDIR)/dg3-mass-mass sparse.slurm $(MFL)/DG3-mass-$(TYP) $(MFL)/DG3-mass-$(TYP)
	sbatch -o $(OUTDIR)/dg3-ip-mass sparse.slurm $(MFL)/DG3-ip-laplace-$(TYP) $(MFL)/DG3-mass-$(TYP)
	sbatch -o $(OUTDIR)/dg3-mass-ip sparse.slurm $(MFL)/DG3-mass-$(TYP) $(MFL)/DG3-ip-laplace-$(TYP)
	sbatch -o $(OUTDIR)/dg3-ip-ip sparse.slurm $(MFL)/DG3-ip-laplace-$(TYP) $(MFL)/DG3-ip-laplace-$(TYP)
	sbatch -o $(OUTDIR)/dg4-mass-mass sparse.slurm $(MFL)/DG4-mass-$(TYP) $(MFL)/DG4-mass-$(TYP)
	sbatch -o $(OUTDIR)/dg4-ip-mass sparse.slurm $(MFL)/DG4-ip-laplace-$(TYP) $(MFL)/DG4-mass-$(TYP)
	sbatch -o $(OUTDIR)/dg4-mass-ip sparse.slurm $(MFL)/DG4-mass-$(TYP) $(MFL)/DG4-ip-laplace-$(TYP)
	sbatch -o $(OUTDIR)/dg4-ip-ip sparse.slurm $(MFL)/DG4-ip-laplace-$(TYP) $(MFL)/DG4-ip-laplace-$(TYP)
	sbatch -o $(OUTDIR)/cg1-mass-mass sparse.slurm $(MFL)/$(SZ)-CG1-mass-3D $(MFL)/$(SZ)-CG1-mass-3D
	sbatch -o $(OUTDIR)/cg1-lp-mass sparse.slurm $(MFL)/$(SZ)-CG1-laplace-3D $(MFL)/$(SZ)-CG1-mass-3D
	sbatch -o $(OUTDIR)/cg1-mass-lp sparse.slurm $(MFL)/$(SZ)-CG1-mass-3D $(MFL)/$(SZ)-CG1-laplace-3D
	sbatch -o $(OUTDIR)/cg1-lp-lp sparse.slurm $(MFL)/$(SZ)-CG1-laplace-3D $(MFL)/$(SZ)-CG1-laplace-3D
	sbatch -o $(OUTDIR)/cg2-mass-mass sparse.slurm $(MFL)/$(SZ)-CG2-mass-3D $(MFL)/$(SZ)-CG2-mass-3D
	sbatch -o $(OUTDIR)/cg2-lp-mass sparse.slurm $(MFL)/$(SZ)-CG2-laplace-3D $(MFL)/$(SZ)-CG2-mass-3D
	sbatch -o $(OUTDIR)/cg2-mass-lp sparse.slurm $(MFL)/$(SZ)-CG2-mass-3D $(MFL)/$(SZ)-CG2-laplace-3D
	sbatch -o $(OUTDIR)/cg2-lp-lp sparse.slurm $(MFL)/$(SZ)-CG2-laplace-3D $(MFL)/$(SZ)-CG2-laplace-3D

vec_report: optimised-sparsemm.c $(OBJ)
	$(CC) $(CFLAGS) $(PROF) -c -o $@ $<
	$(CC) $(CFLAGS) $(VEC) $(PROF) -o $@ $< $(OBJ) $(LDFLAGS)

sparsemm: sparsemm.c $(OBJ)
	$(CC) $(CFLAGS) -std=c99 -o $@ $< $(OBJ) $(LDFLAGS)

%.o: %.c $(HEADER)
	$(CC) $(CFLAGS) -std=c99 -c -o $@ $<

# sparsemm: sparsemm.c $(OBJ)
# 	$(CC) $(CFLAGS) $(PROF) -std=c99 -o $@ $< $(OBJ) $(LDFLAGS)

# %.o: %.c $(HEADER)
# 	$(CC) $(CFLAGS) $(PROF) -std=c99 -c -o $@ $<

profile2:
	sbatch -o $(OUTDIR)/dg1-mmm-mmm sparse2.slurm $(MFL)/DG1-mass-$(TYP) $(MFL)/DG1-mass-$(TYP) $(MFL)/DG1-mass-$(TYP) $(MFL)/DG1-mass-$(TYP) $(MFL)/DG1-mass-$(TYP) $(MFL)/DG1-mass-$(TYP)
	sbatch -o $(OUTDIR)/dg1-iii-iii sparse2.slurm $(MFL)/DG1-ip-laplace-$(TYP) $(MFL)/DG1-ip-laplace-$(TYP) $(MFL)/DG1-ip-laplace-$(TYP) $(MFL)/DG1-ip-laplace-$(TYP) $(MFL)/DG1-ip-laplace-$(TYP) $(MFL)/DG1-ip-laplace-$(TYP)
	sbatch -o $(OUTDIR)/dg1-imi-mim sparse2.slurm $(MFL)/DG1-ip-laplace-$(TYP) $(MFL)/DG1-mass-$(TYP) $(MFL)/DG1-ip-laplace-$(TYP) $(MFL)/DG1-mass-$(TYP) $(MFL)/DG1-ip-laplace-$(TYP) $(MFL)/DG1-mass-$(TYP)
	sbatch -o $(OUTDIR)/dg1-mim-imi sparse2.slurm $(MFL)/DG1-mass-$(TYP) $(MFL)/DG1-ip-laplace-$(TYP) $(MFL)/DG1-mass-$(TYP) $(MFL)/DG1-ip-laplace-$(TYP) $(MFL)/DG1-mass-$(TYP) $(MFL)/DG1-ip-laplace-$(TYP) 

	sbatch -o $(OUTDIR)/dg2-mmm-mmm sparse2.slurm $(MFL)/DG2-mass-$(TYP) $(MFL)/DG2-mass-$(TYP) $(MFL)/DG2-mass-$(TYP) $(MFL)/DG2-mass-$(TYP) $(MFL)/DG2-mass-$(TYP) $(MFL)/DG2-mass-$(TYP)
	sbatch -o $(OUTDIR)/dg2-iii-iii sparse2.slurm $(MFL)/DG2-ip-laplace-$(TYP) $(MFL)/DG2-ip-laplace-$(TYP) $(MFL)/DG2-ip-laplace-$(TYP) $(MFL)/DG2-ip-laplace-$(TYP) $(MFL)/DG2-ip-laplace-$(TYP) $(MFL)/DG2-ip-laplace-$(TYP)
	sbatch -o $(OUTDIR)/dg2-imi-mim sparse2.slurm $(MFL)/DG2-ip-laplace-$(TYP) $(MFL)/DG2-mass-$(TYP) $(MFL)/DG2-ip-laplace-$(TYP) $(MFL)/DG2-mass-$(TYP) $(MFL)/DG2-ip-laplace-$(TYP) $(MFL)/DG2-mass-$(TYP)
	sbatch -o $(OUTDIR)/dg2-mim-imi sparse2.slurm $(MFL)/DG2-mass-$(TYP) $(MFL)/DG2-ip-laplace-$(TYP) $(MFL)/DG2-mass-$(TYP) $(MFL)/DG2-ip-laplace-$(TYP) $(MFL)/DG2-mass-$(TYP) $(MFL)/DG2-ip-laplace-$(TYP) 

	sbatch -o $(OUTDIR)/dg3-mmm-mmm sparse2.slurm $(MFL)/DG3-mass-$(TYP) $(MFL)/DG3-mass-$(TYP) $(MFL)/DG3-mass-$(TYP) $(MFL)/DG3-mass-$(TYP) $(MFL)/DG3-mass-$(TYP) $(MFL)/DG3-mass-$(TYP)
	sbatch -o $(OUTDIR)/dg3-iii-iii sparse2.slurm $(MFL)/DG3-ip-laplace-$(TYP) $(MFL)/DG3-ip-laplace-$(TYP) $(MFL)/DG3-ip-laplace-$(TYP) $(MFL)/DG3-ip-laplace-$(TYP) $(MFL)/DG3-ip-laplace-$(TYP) $(MFL)/DG3-ip-laplace-$(TYP)
	sbatch -o $(OUTDIR)/dg3-imi-mim sparse2.slurm $(MFL)/DG3-ip-laplace-$(TYP) $(MFL)/DG3-mass-$(TYP) $(MFL)/DG3-ip-laplace-$(TYP) $(MFL)/DG3-mass-$(TYP) $(MFL)/DG3-ip-laplace-$(TYP) $(MFL)/DG3-mass-$(TYP)
	sbatch -o $(OUTDIR)/dg3-mim-imi sparse2.slurm $(MFL)/DG3-mass-$(TYP) $(MFL)/DG3-ip-laplace-$(TYP) $(MFL)/DG3-mass-$(TYP) $(MFL)/DG3-ip-laplace-$(TYP) $(MFL)/DG3-mass-$(TYP) $(MFL)/DG3-ip-laplace-$(TYP)

	sbatch -o $(OUTDIR)/dg4-mmm-mmm sparse2.slurm $(MFL)/DG4-mass-$(TYP) $(MFL)/DG4-mass-$(TYP) $(MFL)/DG4-mass-$(TYP) $(MFL)/DG4-mass-$(TYP) $(MFL)/DG4-mass-$(TYP) $(MFL)/DG4-mass-$(TYP)
	sbatch -o $(OUTDIR)/dg4-iii-iii sparse2.slurm $(MFL)/DG4-ip-laplace-$(TYP) $(MFL)/DG4-ip-laplace-$(TYP) $(MFL)/DG4-ip-laplace-$(TYP) $(MFL)/DG4-ip-laplace-$(TYP) $(MFL)/DG4-ip-laplace-$(TYP) $(MFL)/DG4-ip-laplace-$(TYP)
	sbatch -o $(OUTDIR)/dg4-imi-mim sparse2.slurm $(MFL)/DG4-ip-laplace-$(TYP) $(MFL)/DG4-mass-$(TYP) $(MFL)/DG4-ip-laplace-$(TYP) $(MFL)/DG4-mass-$(TYP) $(MFL)/DG4-ip-laplace-$(TYP) $(MFL)/DG4-mass-$(TYP)
	sbatch -o $(OUTDIR)/dg4-mim-imi sparse2.slurm $(MFL)/DG4-mass-$(TYP) $(MFL)/DG4-ip-laplace-$(TYP) $(MFL)/DG4-mass-$(TYP) $(MFL)/DG4-ip-laplace-$(TYP) $(MFL)/DG4-mass-$(TYP) $(MFL)/DG4-ip-laplace-$(TYP)

	sbatch -o $(OUTDIR)/cg1-mmm-mmm sparse2.slurm $(MFL)/$(SZ)-CG1-mass-3D $(MFL)/$(SZ)-CG1-mass-3D $(MFL)/$(SZ)-CG1-mass-3D $(MFL)/$(SZ)-CG1-mass-3D $(MFL)/$(SZ)-CG1-mass-3D $(MFL)/$(SZ)-CG1-mass-3D
	sbatch -o $(OUTDIR)/cg1-iii-iii sparse2.slurm $(MFL)/$(SZ)-CG1-laplace-3D $(MFL)/$(SZ)-CG1-laplace-3D $(MFL)/$(SZ)-CG1-laplace-3D $(MFL)/$(SZ)-CG1-laplace-3D $(MFL)/$(SZ)-CG1-laplace-3D $(MFL)/$(SZ)-CG1-laplace-3D
	sbatch -o $(OUTDIR)/cg1-imi-mim sparse2.slurm $(MFL)/$(SZ)-CG1-laplace-3D $(MFL)/$(SZ)-CG1-mass-3D $(MFL)/$(SZ)-CG1-laplace-3D $(MFL)/$(SZ)-CG1-mass-3D $(MFL)/$(SZ)-CG1-laplace-3D $(MFL)/$(SZ)-CG1-mass-3D
	sbatch -o $(OUTDIR)/cg1-mim-imi sparse2.slurm $(MFL)/$(SZ)-CG1-mass-3D $(MFL)/$(SZ)-CG1-laplace-3D $(MFL)/$(SZ)-CG1-mass-3D $(MFL)/$(SZ)-CG1-laplace-3D $(MFL)/$(SZ)-CG1-mass-3D $(MFL)/$(SZ)-CG1-laplace-3D

	sbatch -o $(OUTDIR)/cg2-mmm-mmm sparse2.slurm $(MFL)/$(SZ)-CG2-mass-3D $(MFL)/$(SZ)-CG2-mass-3D $(MFL)/$(SZ)-CG2-mass-3D $(MFL)/$(SZ)-CG2-mass-3D $(MFL)/$(SZ)-CG2-mass-3D $(MFL)/$(SZ)-CG2-mass-3D
	sbatch -o $(OUTDIR)/cg2-iii-iii sparse2.slurm $(MFL)/$(SZ)-CG2-laplace-3D $(MFL)/$(SZ)-CG2-laplace-3D $(MFL)/$(SZ)-CG2-laplace-3D $(MFL)/$(SZ)-CG2-laplace-3D $(MFL)/$(SZ)-CG2-laplace-3D $(MFL)/$(SZ)-CG2-laplace-3D
	sbatch -o $(OUTDIR)/cg2-imi-mim sparse2.slurm $(MFL)/$(SZ)-CG2-laplace-3D $(MFL)/$(SZ)-CG2-mass-3D $(MFL)/$(SZ)-CG2-laplace-3D $(MFL)/$(SZ)-CG2-mass-3D $(MFL)/$(SZ)-CG2-laplace-3D $(MFL)/$(SZ)-CG2-mass-3D
	sbatch -o $(OUTDIR)/cg2-mim-imi sparse2.slurm $(MFL)/$(SZ)-CG2-mass-3D $(MFL)/$(SZ)-CG2-laplace-3D $(MFL)/$(SZ)-CG2-mass-3D $(MFL)/$(SZ)-CG2-laplace-3D $(MFL)/$(SZ)-CG2-mass-3D $(MFL)/$(SZ)-CG2-laplace-3D
