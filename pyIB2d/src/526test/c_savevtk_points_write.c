/*Rewritten of Savevtk_points, add int N as additional input */
#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
void c_savevtk_points_write(int N, double *x, char* filename, char* vectorName){
    int i,j;
    FILE *fptr = fopen(filename, "w");
    if (fptr == NULL){
        puts("no file exists");
    }
    else{
        fprintf(fptr, "# vtk DataFile Version 2.0\n");
        fprintf(fptr, "%s\n",vectorName);
        fprintf(fptr, "ASCII\n");
        fprintf(fptr, "DATASET UNSTRUCTURED_GRID\n\n");
        fprintf(fptr,"POINTS %d float\n",N);
        for (int i = 0; i < N; i++){
            fprintf(fptr, "%.15f %.15f %.15f\n", x[i*3+0],x[i*3+1],x[i*3+2]);
        }
        fprintf(fptr, "\n");
        fprintf(fptr,"CELLS %d %d\n",N,2*N);
        for (int i = 0; i < N; i++){
            fprintf(fptr, "%d %d\n", 1, i);
        }
        fprintf(fptr,"\n");
        for (int i = 0; i < N; i++){
            fprintf(fptr, "%d",1);
        }
        fprintf(fptr,"1");
        fprintf(fptr,"\n");
        fclose(fptr);
    }
}
/*Savevtk_points ends */