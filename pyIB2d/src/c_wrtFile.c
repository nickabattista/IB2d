//
//  main.cpp
//  Savetk
//
//  Created by Ao Zeng on 5/17/16.
//  Copyright © 2016 Ao Zeng. All rights reserved.
//
#include <assert.h>
#include <iostream>
#include <stdio.h>
using namespace std;
/* Rewritten of savevtk_scalar(array, filename, colorMap,dx,dy) with int arrayRow,int arrayCol as additional inputs */

void savevtk_scalar(int arrayRow,int arrayCol,double** array, char* filename, char* colorMap,double dx, double dy){
    int i,j;
    FILE *fptr = fopen(filename, "w");
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
            fprintf(fptr, "%f ", array[i][j]);
        }
        fprintf(fptr,"\n");
    }
    
}
/* savevtk_scalar(array, filename, colorMap,dx,dy) ends*/



/* Rewritten of savevtk_vector(X, Y, filename, vectorName,dx,dy), int Xrow, int Xcol, int Yrow, int Ycol are added as additional inputs*/
void savevtk_vector(int Xrow, int Xcol, int Yrow, int Ycol, double** X, double** Y,char* filename,char* vectorName,double dx,double dy){
    int i,j;
    FILE *fptr = fopen(filename, "w");
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
    fprintf(fptr,"Vectors %s double\n",vectorName);
    fprintf(fptr, "\n");
    for (int i = 0; i< Xcol;i++){
        for (int j = 0; j<Ycol;j++){
            fprintf(fptr,"%f ",X[i][j]);
            fprintf(fptr,"%f ",Y[i][j]); /*The length of the number should e-15 !!!!!!!!! */
            fprintf(fptr,"1");
        }
        fprintf(fptr,"\n");
    }
    fclose(fptr);
}
/*savevtk_vector(X, Y, filename, vectorName,dx,dy) ends */


/*Rewritten of Savevtk_points, add int N as additional input */
void savevtk_points_write(int N, double **x, char* filename, char* vectorName){
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
            fprintf(fptr, "%.15f %.15f %.15f\n", x[i][0],x[i][1],x[i][2]);
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





/*Rewritten of Savevtk_points_connects, add int N, int Nc as additional input */
void savevtk_points_connects_write(int N, int Nc, double** x,char* filename, char* vectorname, double** connectsMat){
    FILE *fptr = fopen(filename, "w");
    if (fptr == NULL){
            puts("no file exists");
        }
        else{
            fprintf(fptr,"This is the header for the Test!\n"); /*Declare this is a test */
            fprintf(fptr,"# vtk DataFile Version 2.0 \n");
            fprintf(fptr,"%s",vectorname);
            fprintf(fptr,"\n");
            fprintf(fptr,"ASCII\n");
            fprintf(fptr,"DATASET UNSTRUCTED_GRID\n \n");
            fprintf(fptr, "POINTS %d float\n",N);
            for (int i = 0; i < N; i++){
                fprintf(fptr, "%.15f %.15f %.15f \n", x[i][0],x[i][1],x[i][2]);
            }
            fprintf(fptr, "\n");
            
            fprintf(fptr,"CELL %d %d\n", Nc, Nc*3);
            
            for(int i = 0; i < Nc; i++){
                fprintf(fptr, "%d %f %f\n",2,connectsMat[i][0],connectsMat[i][1]);
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
/*Savevtk_points_connects ends */

int main(){
    int i,j;
    char filename[] = "/Users/aozeng/Documents/Comp411/MiniProject/test.txt";
    char vectorname[] =  "vector";
    double** x1;
    double** connectsMat1;
    x1 = new double *[100000];
    for(int i = 0; i <100000; i++){
        x1[i] = new double[3];
    }
    for (int i=0;i<100000;i++){
        for (int j = 0; j<3;j++)
            x1[i][j] = i*3+j;
    }
    
    connectsMat1 = new double *[100000];
    for(int i = 0; i <100000; i++){
        connectsMat1[i] = new double[2];
    }
    for (int i=0;i<100000;i++){
        for (int j = 0; j<2;j++)
            connectsMat1[i][j] = i*2+j;
    }
    savevtk_points_connects_write(100000, 100000, x1, filename, vectorname, connectsMat1);
    return 0;
}







