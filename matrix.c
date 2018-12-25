/* 
 * auth = peic 
 * 对矩阵乘法进行优化， 对矩阵进行转置 
 */  
  
#include <stdio.h>  
#include <stdlib.h>  
#include <time.h>
#include <omp.h>  
  
#define N 1024  
void generate(int *A,int *B,int *C,int size)
{
        for(int i = 0; i < size; i++)  
        for(int j = 0; j < size; j++)  
        {  
            *(A+i*size+j) = 1;  
            *(B+i*size+j) = 2;  
            *(C+i*size+j) = 0;  
        }
}
void save(char name[],int *C,int size)
{
    FILE *fp;
    if ((fp = fopen(name,"w"))==NULL)
    {
    printf("the file can not open..");
    exit(0);
    }//出错处理
    for(int i = 0; i < size; i++) 
    { 
        for(int j = 0; j < size; j++) 
        { 
            fprintf(fp,"%5d", *(C+i*size+j)); 
     
        } 
        fprintf(fp,"\n"); 
     
    } 
    fclose(fp);//OK就这么简单
}  
void check()
{

} 
void MartixMultip(int* A, int* B, int* C, int size)  
{ 
    generate(A,B,C,size);    
    int i, j, k;  
    //先对矩阵进行转置  
    int *tmp = (int *)malloc(sizeof(int) * size * size);  
    for(i=0; i<size; i++)  
        for(j=0; j<size; j++)  
        {  
            *(tmp+i*size+j) = *(B+j*size+i);  
        }  
    //对转置后的矩阵做乘法，具有空间局域性 
    double begin_time = (double)clock() ;   
    for(i=0; i < size; i++)  
    {     
        for(j=0; j<size; j++)  
        {     
            for(k=0; k<size; k++)  
            {  
                //i表示行,j表示列,k表示游标  
                *(C+i*size+j) += (*(A+i*size+k)) * (*(tmp+j*size+k));  
            }  
        }  
    }
 
    double end_time = (double)clock() ;  
    double cost_time = end_time - begin_time;  
    printf("Multip cost time: %fms\n", cost_time); 
    save("multp.txt",C,size);  
    free(tmp);  
}  

void compute_normal(int *A,int *B,int *C,int size)
{
    generate(A,B,C,size);   
    
    int i, j, k;
    double begin_time = (double)clock() ;  
    for(i=0; i < size; i++)  
    {     
        for(j=0; j<size; j++)  
        {     
            for(k=0; k<size; k++)  
            {  
                //i表示行,j表示列,k表示游标  
                *(C+i*size+j) += (*(A+i*size+k)) * (*(B+k*size+j));  
            }  
        }  
    }  
    /*输出运算结果 
    for(int i = 0; i < size; i++) 
    { 
        for(int j = 0; j < size; j++) 
        { 
            printf("%4d", *(C+i*size+j)); 
     
        } 
        printf("\n"); 
     
    } 
    */  
    double end_time = (double)clock();  // / CLOCKS_PER_SEC
    double cost_time = end_time - begin_time;  
    printf("normal cost time: %fms\n", cost_time);
    save("normal.txt",C,size);   
}


void MartixBlock(int  *A, int *B, int *C, int n, int bsize) 
{
    int r, c, k, kk, cc;
    int sum;
    int en = bsize * (n/bsize); /* Amount that frts evenly into blocks */
  
    generate(A,B,C,n);  
    double begin_time = (double)clock() ;  
    for (kk = 0; kk < en; kk += bsize) { 
        for (cc = 0; cc < en; cc += bsize) {
            for (r = 0; r < n; r++) {
                for (c = cc; c < cc + bsize; c++) {
                    sum = *(C+r*n+c);
                    for (k = kk; k < kk + bsize; k++) {
                         sum += (*(A+r*n+k))*(*(B+k*n+c));
                    }
                    *(C+r*n+c) = sum;
                }
            }
        }
    }
    double end_time = (double)clock();  // / CLOCKS_PER_SEC
    double cost_time = end_time - begin_time;  
    printf("Block cost time: %fms\n", cost_time);
    save("Block.txt",C,n);  
}


void MartixOpenmp(int* A, int* B, int* C, int size)  
{  
    int i, j, k;  
    //先对矩阵进行转置 
    generate(A,B,C,size);  
    int *tmp = (int *)malloc(sizeof(int) * size * size);  
    for(i=0; i<size; i++)  
        for(j=0; j<size; j++)  
        {  
            *(tmp+i*size+j) = *(B+j*size+i);  
        }  
      
    //对转置后的矩阵做乘法，具有空间局域性 
    double begin_time = omp_get_wtime() ; 
    #pragma omp parallel shared(A,B,C) private(i,j,k)  
    {  
        #pragma omp for schedule(dynamic)  
        for(i=0; i < size; i++)  
        {     
            for(j=0; j<size; j++)  
            {     
                for(k=0; k<size; k++)  
                {  
                    //i表示行,j表示列,k表示游标  
                    *(C+i*size+j) += (*(A+i*size+k)) * (*(tmp+j*size+k));  
                }  
                
            }  
        } 
    }
    double end_time = omp_get_wtime();  // / CLOCKS_PER_SEC
    double cost_time = end_time - begin_time;  
    printf("openmp cost time: %fms\n", cost_time*1000000);
    save("openmp.txt",C,size);  
}  
int main()  
{ 
    int size = N; 
    int* A = (int *)malloc(sizeof(int) * size * size);  
    int* B = (int *)malloc(sizeof(int) * size * size);  
    int* C = (int *)malloc(sizeof(int) * size * size); 
    compute_normal(A,B,C,size); 
    MartixMultip(A,B,C,size);
    MartixBlock(A,B,C,size,4); 
    MartixOpenmp(A,B,C,size);
    free(A);  
    free(B);  
    free(C); 
    return 0; 
}
