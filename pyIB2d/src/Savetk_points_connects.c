//
//  main.cpp
//  Savevtk_points_connects
//
//  Created by Ao Zeng on 5/13/16.
//  Copyright © 2016 Ao Zeng. All rights reserved.
//

#include <iostream>
#include <stdio.h>
using namespace std;

class Savevtk_points_connects
{
    double x[100000][3];
    char filename[100];
    char vectorname[100];
    double connectsMat[100000][2];
    int N;
    int Nc;
    FILE *fptr;
    int i,j;

public:
    
    Savevtk_points_connects (int newN, int newNc, double newX[100000][3],char newFilename[100], char newVectorname[100], double newConnectsMat[100000][2]);

    void write(){
        if (fptr == NULL){
            puts("no file exists");
        }
        else{
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
    }

};

    
Savevtk_points_connects::Savevtk_points_connects (int newN, int newNc, double newX[100000][3],char newFilename[100], char newVectorname[100], double newConnectsMat[100000][2]){
    for (int i = 0; i< newN; i++){
        for (int j =0; j<2; j++){
            x[i][j] = newX[i][j];
        }
    }
    for (int i = 0; i< newN; i++){
        for (int j =0; j<2; j++){
            connectsMat[i][j] = newConnectsMat[i][j];
        }
    }
    N = newN;
    Nc = newNc;
    strcpy(filename, newFilename);
    strcpy(vectorname, newVectorname);
    fptr = fopen(filename, "w");
        };


int main(){
    int i,j;
    char filename[] = "/Users/aozeng/Documents/Comp411/MiniProject/test.txt";
    char vectorname [] =  "vector\0";
    double x1[100000][3];
    double connectsMat1[100000][2];
    for (int i=0;i<100000;i++){
        for (int j = 0; j<3;j++)
            x1[i][j] = i*3+j;
    }
    for (int i=0;i<100000;i++){
        for (int j = 0; j<2;j++)
            connectsMat1[i][j] = i*2+j;
    }
    Savevtk_points_connects sample = Savevtk_points_connects(100000, 100000, x1, filename, vectorname, connectsMat1);
    sample.write();
    return 0;
}







