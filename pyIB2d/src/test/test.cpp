#include <stdio.h>


double** add(double** array1, double** array2, int row1, int col1, int row2, int col2){
    for (int i=0;i<row1;i++){
        for (int j = 0; j<col1;j++){
            array1[i][j] = array2[i][j]+array1[i][j];
        }
        printf("\n");
    }
    return array1;
}

int main(){
    double** result;
    double** x1;
    double** x2;
    x1= new double *[10];
    for(int i = 0; i <10; i++){
        x1[i] = new double[3];
    }
    for (int i=0;i<10;i++){
        for (int j = 0; j<3;j++)
            x1[i][j] = i*3+j;
    }
    x2 = new double *[10];
    for(int i = 0; i <10; i++){
        x2[i] = new double[3];
    }
    for (int i=0;i<10;i++){
        for (int j = 0; j<3;j++)
            x2[i][j] = i*3+j;
    }
    
    result = add(x1,x2,10,3,10,3);
    
    for (int i=0;i<10;i++){
        for (int j = 0; j<3;j++){
            printf("%f ",result[i][j]);
        }
        printf("\n");
    }
    return 0;
}

