/* Test wrtFile.c in C here! */
//
//  main.cpp
//  Savetk
//
//  Created by Ao Zeng on 5/17/16.
//  Copyright © 2016 Ao Zeng. All rights reserved.
//
#include <assert.h>
#include <stdlib.h>
#include <stdio.h>

/* Rewritten of savevtk_scalar(array, filename, colorMap,dx,dy) with int arrayRow,int arrayCol as additional inputs */

void c_savevtk_scalar(int arrayRow,int arrayCol,double* array, char* filename, char* colorMap,double dx, double dy){
    FILE *fptr = fopen(filename, "w");
    if (fptr == NULL){
        puts("no file exists");
    }
    else{
    fprintf(fptr, "# vtk DataFile Version 2.0\n");
    fprintf(fptr, "Comment goes here\n");
    fprintf(fptr, "ASCII\n");
    fprintf(fptr,"\n");
    fprintf(fptr, "DATASET STRUCTURED_POINTS\n");
    fprintf(fptr, "DIMENSIONS    %d   %d   %d\n",arrayRow ,arrayCol,1);
    fprintf(fptr,"\n");
    fprintf(fptr, "ORIGIN    0.000   0.000   0.000\n");
    fprintf(fptr,"SPACING   %f%f   1.000\n",dx, dy);
    fprintf(fptr,"\n");
    fprintf(fptr,"POINT_DATA   %d\n",arrayCol*arrayCol*1);
    fprintf(fptr,"SCALARS %s double\n",colorMap);
    fprintf(fptr,"LOOKUP_TABLE default\n");
    fprintf(fptr,"\n");
    for(int i = 0; i < arrayRow; i++){
        for (int j = 0; j < arrayCol; j++){
            fprintf(fptr, "%f ", array[i*arrayCol+j]);
        }
        fprintf(fptr,"\n");
    }
    fclose(fptr);
}
}
/* savevtk_scalar(array, filename, colorMap,dx,dy) ends*/



/* Rewritten of savevtk_vector(X, Y, filename, vectorName,dx,dy), int Xrow, int Xcol, int Yrow, int Ycol are added as additional inputs*/
void c_savevtk_vector(int Xrow, int Xcol, int Yrow, int Ycol, double* X, double* Y,char* filename,char* vectorName,double dx,double dy){
    FILE *fptr = fopen(filename, "w");
    if (fptr == NULL){
        puts("no file exists");
    }
    else{
    assert(Xrow == Yrow && Xcol == Ycol && "Error: velocity arrays of unequal size");
    fprintf(fptr,"# vtk DataFile Version 2.0\n");
    fprintf(fptr,"Comment goes here\n");
    fprintf(fptr,"ASCII\n");
    fprintf(fptr,"\n");
    fprintf(fptr,"DATASET STRUCTURED_POINTS\n");
    fprintf(fptr,"DIMENSIONS    %d   %d   %d\n", Xrow, Xcol,1); /*Is this correct? */
    fprintf(fptr,"\n");
    fprintf(fptr,"ORIGIN    0.000   0.000   0.000\n");
    /*fid.write('SPACING   1.000   1.000   1.000\n') #if want [1,32]x[1,32] rather than [0,Lx]x[0,Ly] */
    fprintf(fptr, "SPACING   %f%f   1.000\n",dx,dy);
    fprintf(fptr,"\n");
    fprintf(fptr,"POINT_DATA   %d\n",Xrow*Xcol);
    fprintf(fptr,"VECTORS %s double\n",vectorName);
    fprintf(fptr, "\n");
    for (int i = 0; i< Xrow;i++){
        for (int j = 0; j<Xcol;j++){
            fprintf(fptr,"%f ",X[i*Xcol+j]);
            fprintf(fptr,"%f ",Y[i*Xcol+j]); /*The length of the number should e-15 !!!!!!!!! */
            fprintf(fptr,"0 ");
        }
        fprintf(fptr,"\n");
    }
    fclose(fptr);
}
}
/*savevtk_vector(X, Y, filename, vectorName,dx,dy) ends */


/*Rewritten of Savevtk_points, add int N as additional input */
void c_savevtk_points_write(int N, double *x, char* filename, char* vectorName){
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
            fprintf(fptr, "%.15e %.15e %.15e\n", x[i*3+0],x[i*3+1],x[i*3+2]);
        }
        fprintf(fptr, "\n");
        fprintf(fptr,"CELLS %d %d\n",N,2*N);
        for (int i = 0; i < N; i++){
            fprintf(fptr, "%d %d\n", 1, i);
        }
        fprintf(fptr,"\n");
        fprintf(fptr,"CELL_TYPES %d\n",N);
        for (int i = 0; i < N; i++){
            fprintf(fptr, "%d ",1);
        }
        fprintf(fptr,"\n");
        fclose(fptr);
    }
}
/*Savevtk_points ends */
 




/*Rewritten of Savevtk_points_connects, add int N, int Nc as additional input */
void c_savevtk_points_connects_write(int N, int Nc, double* x,char* filename, char* vectorname, double* connectsMat){
    FILE *fptr = fopen(filename, "w");
    if (fptr == NULL){
        puts("no file exists");
    }
    else{
        fprintf(fptr,"# vtk DataFile Version 2.0 \n");
        fprintf(fptr,"%s",vectorname);
        fprintf(fptr,"\n");
        fprintf(fptr,"ASCII\n");
        fprintf(fptr,"DATASET UNSTRUCTURED_GRID\n\n");
        fprintf(fptr, "POINTS %d float\n",N);
        for (int i = 0; i < N; i++){
            fprintf(fptr, "%.15e %.15e %.15e \n", x[i*3+0],x[i*3+1],x[i*3+2]);
        }
        fprintf(fptr, "\n");
        
        fprintf(fptr,"CELLS %d %d\n", Nc, Nc*3);
        
        for(int i = 0; i < Nc; i++){
            fprintf(fptr, "%d %d %d\n",2,(int) connectsMat[i*2+0],(int) connectsMat[i*2+1]);
        }
        fprintf(fptr, "\n");
        
        fprintf(fptr,"CELL_TYPES %d\n", Nc);
        for (int i = 0; i<Nc; i++){
            fprintf(fptr, "3 ");
        }
        fprintf(fptr,"\n");
        fclose(fptr);
    }
};
