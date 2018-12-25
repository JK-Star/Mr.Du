#include <mpi.h>
#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#define M 1024
#define N 1024 
void generate(int *A,int *B,int size)
{
        for(int i = 0; i < size; i++)  
        for(int j = 0; j < size; j++)  
        {  
            *(A+i*size+j) = 1;  
            *(B+i*size+j) = 2;    
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
    fclose(fp);
} 
int main(int argv, char *argc[])
{

	int my_rank;/*My process rank*/ 
	int comm_sz;/*Number of processes*/ 
	int local_M; int i,j,k; 
	double start,finish;/*timer*/ 
	int tem; 
	//初始化
	MPI_Init(NULL,NULL); 
	MPI_Comm_rank(MPI_COMM_WORLD,&my_rank); 
	MPI_Comm_size(MPI_COMM_WORLD,&comm_sz); 
	//每个矩阵分配到的行数 
	local_M=M/comm_sz; //分配到每个进程的矩阵 
	int *local_Matrix_one=(int*)malloc(local_M*N*sizeof(int)); //定义两个矩阵 
	int *Matrix_one=NULL; 
	int *Matrix_two=(int*)malloc(M*N*sizeof(int)); //每个进程里的结果矩阵 
	int *local_result=(int*)malloc(local_M*N*sizeof(int)); //结果矩阵 
	int *result_Matrix=NULL; 
	start=MPI_Wtime();
	if(my_rank==0) { 

		Matrix_one=(int*)malloc(M*N*sizeof(int));
		generate(Matrix_one,Matrix_two,M);
		//数据分发 
		MPI_Scatter(Matrix_one,local_M*N,MPI_INT,local_Matrix_one,local_M*N,MPI_INT,0,MPI_COMM_WORLD); 
		//数据广播 
		MPI_Bcast(Matrix_two,M*N,MPI_INT,0,MPI_COMM_WORLD); //计算local结果 
		for(i=0;i<local_M;i++) 
		for(j=0;j<M;j++)
		{ 
			tem=0; 
			for(k=0;k<N;k++) 
			tem +=local_Matrix_one[i*M+k]*Matrix_two[j*M+k]; 
			local_result[i*M+j]=tem; 
		} 
		free(local_Matrix_one); 
		result_Matrix=(int*)malloc(M*N*sizeof(int)); 
		//结果聚集 
		MPI_Gather(local_result,local_M*N,MPI_INT,result_Matrix,local_M*N,MPI_INT,0,MPI_COMM_WORLD); //剩余行处理（处理不能整除的情况） 
		int rest=M%comm_sz; 
		if(rest!=0) 
		for(i=M-rest-1;i<M;i++) 
		for(j=0;j<M;j++)
		{ 
			tem=0; 
			for(k=0;k<N;k++) 
			tem +=Matrix_one[i*M+k]*Matrix_two[j*M+k]; 
			result_Matrix[i*M+j]=tem; 
		} 
		finish=MPI_Wtime(); 
		free(Matrix_one); 
		free(Matrix_two); 
		free(local_result);
		printf("rank:%d \t time:%f ms\n",my_rank,(finish-start)*1000000);  
		//将结果写入文件 
		save("mpi_ScatGat.txt",result_Matrix,M);
	}
			 
	else{ //printf("process %d of %d\n",my_rank,comm_sz); 
	//数据分发 
	MPI_Scatter(Matrix_one,local_M*N,MPI_INT,local_Matrix_one,local_M*N,MPI_INT,0,MPI_COMM_WORLD); 
	//数据广播 
	MPI_Bcast(Matrix_two,M*N,MPI_INT,0,MPI_COMM_WORLD);
	 //计算local结果 
	 for(i=0;i<local_M;i++) 
	 for(j=0;j<M;j++)
	 { 
	 	tem=0; 
	 	for(k=0;k<N;k++) 
	 	tem +=local_Matrix_one[i*M+k]*Matrix_two[j*M+k]; 
	 	local_result[i*M+j]=tem; 
	 } 
	 free(local_Matrix_one); 
	 free(Matrix_two); 
	 //结果聚集 
	 MPI_Gather(local_result,local_M*N,MPI_INT,result_Matrix,local_M*N,MPI_INT,0,MPI_COMM_WORLD); 
	 free(local_result); //printf("%d %d\n",local_M,my_rank); 
	} 
	MPI_Finalize(); 
	return 0; 
}
