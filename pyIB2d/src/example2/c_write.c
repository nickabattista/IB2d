/*
c_write.c

simple C function that alters data passed in via a pointer

    used to see how we can do this with Cython/numpy

*/

void c_write (double* array) {

    int i ,j;

    for (i = 0; i < 2; i++) {
        for(j = 0; j<2;j++){
            printf("%f",array[i*2+j]);
        }
        }

    return ;
}