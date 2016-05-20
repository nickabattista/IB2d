#include <stdio.h>


void add(double* array1, double* array2, int row1, int col1){
    int index =0;
    for (int i=0;i<row1;i++){
        for (int j = 0; j<col1;j++){
            array1[index] = array1[index]+array2[index];
            index++;
        }
    }
}

int main(){
    double x1[30];
    double x2[30];
    int index1 = 0;
    int index2 = 0;
    for (int i=0; i< 10; i++){
        for (int j=0; j<3;j++){
            x1[index1] = index1;
            index1++;
        }
    }
    
    for (int i=0; i< 10; i++){
        for (int j=0; j<3;j++){
            x2[index2] = index2;
            index2++;
        }
    }
    
    add(x1,x2,10,3);
    for (int i=0; i< 10;i++){
        for (int j=0;j<3;j++){
            printf("%f", x1[i*3+j]);
        }
        printf("\n");
    }
    
    return 0;
}

