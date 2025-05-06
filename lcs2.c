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

typedef unsigned short mtype;

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
        diag[k] = (mtype *)calloc(len, sizeof(mtype));
    }
    return diag;
}

void freeDiagMatrix(mtype **diag, int m, int n) {
    int total_diags = m + n + 1;
    for (int k = 0; k < total_diags; k++)
        free(diag[k]);
    free(diag);
}

int LCS_antidiag(char *A, char *B, int m, int n) {
    mtype **diag = allocateDiagMatrix(m, n);
    int total_diags = m + n + 1;

    for (int k = 2; k < total_diags; k++) {
        int i_start = max(1, k - n);
        int i_end = min(k - 1, m);
        int i_start_prev = max(1, (k - 1) - n);

        omp_set_num_threads(4);
        #pragma omp parallel for
        for (int i = i_start; i <= i_end; i++)
        {
            int j = k - i;
            if (j > n) continue; // ?

            // printf("A[%d]: %c == B[%d]: %c?\n", i-1, A[i-1], j-1, B[j-1]);
            // printf("i_start: %d\n", i_start);

            // if (A[i - 1] == B[j - 1])
            // {
            //     printf("sim!\n\n");
            //     printf("guardando em diag[%d][%d]: diag[%d][%d] + 1 -> %d\n\n", k, i-i_start, k-2, i-i_start, diag[k - 2][i - i_start_prev] + 1);
            //     diag[k][i - i_start] = diag[k - 2][i - i_start_prev] + 1;
            // }
            // else
            // {
            //     printf("não!\n\n");
            //     // mtype left = diag[k - 1][i - i_start - 1];
            //     int left_idx = i - 1 - i_start_prev;
            //     mtype left = (i - 1 >= i_start_prev && left_idx >= 0) ? diag[k - 1][left_idx] : 0;
            //     printf("diag[%d][%d]: left: %d\n", k-1, left_idx, left);
            //     // mtype up   = diag[k - 1][i - i_start];
            //     int up_idx = i - i_start_prev;
            //     mtype up = (i >= i_start_prev && up_idx >= 0) ? diag[k - 1][up_idx] : 0;
            //     printf("diag[%d][%d]: up: %d\n\n\n", k-1, up_idx, up);
            //     diag[k][i - i_start] = max(left, up);
            //     printf("guardando em diag[%d][%d]: max(left, up) -> %d\n\n", k, i-i_start, max(left, up));
            // }

            // Check match
            if (A[i - 1] == B[j - 1]) {
                int prev_k = k - 2;
                int prev_i = i - 1;
                int prev_i_start = max(1, prev_k - n);
                int prev_idx = prev_i - prev_i_start;

                diag[k][i - i_start] = diag[prev_k][prev_idx] + 1;
            } else {
                int prev_k = k - 1;

                // Left: (i-1, j)
                int left_i = i - 1;
                int left_i_start = max(1, prev_k - n);
                int left_idx = left_i - left_i_start;
                mtype left = (left_i >= left_i_start && left_idx >= 0) ? diag[prev_k][left_idx] : 0;

                // Up: (i, j-1)
                int up_i = i;
                int up_i_start = max(1, prev_k - n);
                int up_idx = up_i - up_i_start;
                mtype up = (up_i >= up_i_start && up_idx >= 0) ? diag[prev_k][up_idx] : 0;

                diag[k][i - i_start] = max(left, up);
            }
        }

        // Print current antidiagonal
        // printf("Diagonal k=%d: \n", k);
        // for (int i = i_start; i <= i_end; i++) {
            // int idx = i - i_start;
            // printf("%d ", diag[k][idx]);
        // }
        // printf("\n\n\n");
    }

    // int score = diag[m + n][0]; // last diagonal, only one value
    int final_k = m + n;
    int i_start = max(1, final_k - n);
    int idx = m - i_start;
    int score = diag[final_k][idx];

    // // imprime a matriz
    // mtype **scoreMatrix = (mtype **)malloc((m + 1) * sizeof(mtype *));
    // for (int i = 0; i <= m; i++)
    //     scoreMatrix[i] = (mtype *)calloc(n + 1, sizeof(mtype));

    // for (int k = 2; k < total_diags; k++) {
    //     int i_start = max(1, k - n);
    //     int i_end = min(k - 1, m);
    //     for (int i = i_start; i <= i_end; i++) {
    //         int j = k - i;
    //         if (j > n) continue;

    //         int idx = i - i_start;
    //         scoreMatrix[i][j] = diag[k][idx];
    //     }
    // }

    // // Print matrix
    // printf("\nScore Matrix:\n    ");
    // for (int j = 0; j <= n; j++)
    //     printf("%5c", B[j]);
    // printf("\n");

    // for (int i = 1; i <= m; i++) {
    //     printf("%c ", A[i-1]);
    //     for (int j = 1; j <= n; j++) {
    //         printf("%5d", scoreMatrix[i][j]);
    //     }
    //     printf("\n");
    // }

    // // Free the reconstructed matrix
    // for (int i = 0; i <= m; i++)
    //     free(scoreMatrix[i]);
    // free(scoreMatrix);


    freeDiagMatrix(diag, m, n);
    return score;
}

int main(int argc, char **argv) {
    char *seqA = read_seq("fileA.in");
    char *seqB = read_seq("fileB.in");

    int sizeA = strlen(seqA);
    int sizeB = strlen(seqB);

    double start = omp_get_wtime();
    int score = LCS_antidiag(seqA, seqB, sizeA, sizeB);
    double end = omp_get_wtime();

    printf("\nScore: %d\n", score);
    printf("Time: %.6f seconds\n", end - start);

    free(seqA);
    free(seqB);
    return 0;
}
