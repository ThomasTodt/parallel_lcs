#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>

#ifndef max
#define max(a, b) (((a) > (b)) ? (a) : (b))
#endif

#ifndef min
#define min(a, b) (((a) < (b)) ? (a) : (b))
#endif

#define CACHE_LINE 64
#define PAD (CACHE_LINE / sizeof(mtype))

typedef unsigned short mtype;

void printMatrix(char * seqA, char * seqB, mtype ** scoreMatrix, int sizeA, int sizeB)
{
    int i, j;

	//print header
	printf("Score Matrix:\n");
	printf("========================================\n");

	//print LCS score matrix allong with sequences

	printf("    ");
	printf("%5c   ", ' ');

	for (j = 0; j < sizeA; j++)
		printf("%5c   ", seqA[j]);
	printf("\n");
	for (i = 0; i < sizeB + 1; i++) {
		if (i == 0)
			printf("    ");
		else
			printf("%c   ", seqB[i - 1]);
		for (j = 0; j < sizeA + 1; j++) {
			printf("%5d   ", scoreMatrix[i][j]);
		}
		printf("\n");
	}
	printf("========================================\n");
}

char* read_seq(char *fname) {
    FILE *fseq = fopen(fname, "rt");
    if (!fseq) {
        printf("Error reading file %s\n", fname);
        exit(1);
    }

    fseek(fseq, 0L, SEEK_END);
    long size = ftell(fseq);
    rewind(fseq);

    char *seq = (char *)calloc(size + 1, sizeof(char));
    if (!seq) {
        printf("Error allocating memory for sequence %s.\n", fname);
        exit(1);
    }

    int i = 0;
    while (!feof(fseq)) {
        seq[i] = fgetc(fseq);
        if ((seq[i] != '\n') && (seq[i] != EOF))
            i++;
    }
    seq[i] = '\0';
    fclose(fseq);
    return seq;
}

// Allocate diagonal representation of DP matrix
mtype** allocateDiagMatrix(int m, int n) {
    int total = m + n + 1; // numero de linhas
    mtype **diag = (mtype **)malloc(total * sizeof(mtype *));
    for (int k = 0; k < total; k++) {
        int len = (k < m + 1 && k < n + 1) ? k + 1 :
                  (k < max(m, n) + 1) ? min(m, n) + 1 :
                  m + n + 1 - k;
        // printf("len: %d\n", len);
        diag[k] = (mtype *)calloc((len + PAD), sizeof(mtype));

        // alignment
        // size_t size = (len + PAD) * sizeof(mtype);
        // diag[k] = (mtype *)aligned_alloc(CACHE_LINE, size);
        // memset(diag[k], 0, size);  // Zero-initialize, like calloc
    }
    return diag;
}

void freeDiagMatrix(mtype **diag, int m, int n) {
    int total_diags = m + n + 1;
    for (int k = 0; k < total_diags; k++)
        free(diag[k]);
    free(diag);
}

int LCS_antidiag_o(char *A, char *B, int m, int n) {
    int max_len = min(m, n) + 1;

    mtype *D_prev2 = calloc(max_len + PAD, sizeof(mtype));
    mtype *D_prev = calloc(max_len + PAD, sizeof(mtype));
    mtype *D_curr = calloc(max_len + PAD, sizeof(mtype));

    int total_diags = m + n + 1;

    #ifdef DEBUGMATRIX
        mtype **scoreMatrix = malloc((m + 1) * sizeof(mtype *));
        for (int i = 0; i <= m; i++) {
            scoreMatrix[i] = calloc(n + 1, sizeof(mtype));
        }
    #endif

    for (int k = 2; k < total_diags; k++) {
        int i_start = max(1, k - n);
        int i_end = min(k - 1, m);

        #pragma omp parallel for
        for (int i = i_start; i <= i_end; i++) {
            int j = k - i;
            int idx = i - i_start;

            if (A[i - 1] == B[j - 1]) {
                int prev_idx = i - 1 - max(1, k - 2 - n);
                D_curr[idx] = D_prev2[prev_idx] + 1;
            } else {
                int left_idx = i - 1 - max(1, k - 1 - n);
                int up_idx = i - max(1, k - 1 - n);
                mtype left = (left_idx >= 0) ? D_prev[left_idx] : 0;
                mtype up = (up_idx >= 0) ? D_prev[up_idx] : 0;
                D_curr[idx] = max(left, up);
            }

            #ifdef DEBUGMATRIX
                scoreMatrix[j][i] = D_curr[idx];
            #endif
        }

        // rotate buffers
        mtype *temp = D_prev2;
        D_prev2 = D_prev;
        D_prev = D_curr;
        D_curr = temp;
    }

    int final_k = m + n;
    int idx = m - max(1, final_k - n);
    int score = D_prev[idx];


    #ifdef DEBUGMATRIX
        printMatrix(A, B, scoreMatrix, m, n);

        for (int i = 0; i <= m; i++) {
            free(scoreMatrix[i]);
        }
        free(scoreMatrix);
    #endif


    free(D_prev2);
    free(D_prev);
    free(D_curr);
    return score;
}

int LCS_antidiag(char *A, char *B, int m, int n) {
    mtype **diag = allocateDiagMatrix(m, n);
    int total_diags = m + n + 1;

    for (int k = 2; k < total_diags; k++) {
        int i_start = max(1, k - n);
        int i_end = min(k - 1, m);
        
        int prev_k = k - 2;
        int prev_k2 = k - 1;

        #pragma omp parallel for default(none) shared(diag, A, B, k, m, n, i_start, i_end, prev_k, prev_k2)
        for (int i = i_start; i <= i_end; i++)
        {
            int j = k - i;
            if (j > n) continue; // ?

            // Check match
            if (A[i - 1] == B[j - 1]) {
                int prev_i = i - 1;
                int prev_i_start = max(1, prev_k - n);
                int prev_idx = prev_i - prev_i_start;

                diag[k][i - i_start] = diag[prev_k][prev_idx] + 1;
            } else {
                int left_i = i - 1;
                int left_i_start = max(1, prev_k2 - n);
                int left_idx = left_i - left_i_start;
                mtype left = (left_i >= left_i_start && left_idx >= 0) ? diag[prev_k2][left_idx] : 0;

                int up_i = i;
                // int up_idx = up_i - up_i_start;
                int up_idx = up_i - left_i_start;
                // mtype up = (up_i >= up_i_start && up_idx >= 0) ? diag[prev_k2][up_idx] : 0;
                mtype up = (up_i >= left_i_start && up_idx >= 0) ? diag[prev_k2][up_idx] : 0;


                diag[k][i - i_start] = max(left, up);
            }
        }

    }

    int final_k = m + n;
    int i_start = max(1, final_k - n);
    int idx = m - i_start;
    int score = diag[final_k][idx];


    freeDiagMatrix(diag, m, n);
    return score;
}


int main(int argc, char **argv) {
    omp_set_num_threads(4);

    char *seqA = read_seq("fileA.in");
    char *seqB = read_seq("fileB.in");

    int sizeA = strlen(seqA);
    int sizeB = strlen(seqB);

    double start = omp_get_wtime();
    int score = LCS_antidiag_o(seqA, seqB, sizeA, sizeB);
    double end = omp_get_wtime();

    printf("\nScore: %d\n", score);
    printf("Time: %.6f seconds\n", end - start);

    free(seqA);
    free(seqB);
    return 0;
}
