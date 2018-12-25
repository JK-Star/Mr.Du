#include <mpi.h>
#include<stdio.h>
#include<math.h>
#include<stdlib.h>
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
int main(int argv, char *argc[])
{

	double start, stop; 
	int i, j, k, l; 
	int *a, *b, *c, *buffer, *ans; 
	int size = N; 
	int rank, numprocs, line; 
	MPI_Init(NULL,NULL);//MPI Initialize 
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);//获得当前进程号 
	MPI_Comm_size(MPI_COMM_WORLD,&numprocs);//获得进程个数 
	line = size/numprocs;//将数据分为(进程数)个块,主进程也要处理数据 
	a = (int*)malloc(sizeof(int)*size*size); 
	b = (int*)malloc(sizeof(int)*size*size); 
	c = (int*)malloc(sizeof(int)*size*size); 
	//缓存大小大于等于要处理的数据大小，大于时只需关注实际数据那部分 
	buffer = (int*)malloc(sizeof(int)*size*line);
	//数据分组大小 
	ans = (int*)malloc(sizeof(int)*size*line);//保存数据块计算的结果 
	//主进程对矩阵赋初值，并将矩阵N广播到各进程,将矩阵M分组广播到各进程
	start = MPI_Wtime();  
	if (rank==0) { 
		generate(a,b,c,size);
		//将矩阵N发送给其他从进程 
		for (i=1;i<numprocs;i++) 
		{ 
			MPI_Send(b,size*size,MPI_INT,i,0,MPI_COMM_WORLD);
		} 
		 //依次将a的各行发送给各从进程 
		for(l=1; l<numprocs; l++) 
		{ 
		 MPI_Send(a+(l-1)*line*size,size*line,MPI_INT,l,1,MPI_COMM_WORLD);
		} 
		 //接收从进程计算的结果 
		for (k=1;k<numprocs;k++) 
		{ 
			MPI_Recv(ans,line*size,MPI_INT,k,3,MPI_COMM_WORLD,MPI_STATUS_IGNORE); 
		//将结果传递给数组c 
			for (i=0;i<line;i++) 
			{ 
				for (j=0;j<size;j++) 
				{ 
					c[((k-1)*line+i)*size+j] = ans[i*size+j]; 
				} 
			} 
		} 
		//计算a剩下的数据 
		for (i=(numprocs-1)*line;i<size;i++) 
		{ 
			for (j=0;j<size;j++) 
			{ 
				int temp=0; 
				for (k=0;k<size;k++) 
				temp += a[i*size+k]*b[k*size+j]; c[i*size+j] = temp; 
			} 
		}
		////统计时间
		stop = MPI_Wtime(); 

		save("mpi_SendRcev.txt",c,size);
		//结果测试 
		
		printf("rank:%d \t time:%f ms\n",rank,(stop-start)*1000000); 
		free(a); 
		free(b); 
		free(c); 
		free(buffer); 
		free(ans); 
		} 
		//其他进程接收数据，计算结果后，发送给主进程 
		else 
		{ 
			//接收广播的数据(矩阵b) 
			MPI_Recv(b,size*size,MPI_INT,0,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE); 
			MPI_Recv(buffer,size*line,MPI_INT,0,1,MPI_COMM_WORLD,MPI_STATUS_IGNORE); 
			//计算乘积结果，并将结果发送给主进程 
			for (i=0;i<line;i++) 
			{
			 for (j=0;j<size;j++) 
			 {
			  int temp=0; 
			  for(k=0;k<size;k++) 
			  temp += buffer[i*size+k]*b[k*size+j]; 
			  ans[i*size+j]=temp; 
			 } 
			} 
			//将计算结果传送给主进程 
			MPI_Send(ans,line*size,MPI_INT,0,3,MPI_COMM_WORLD); 
		} 
		MPI_Finalize();
		//结束 
		return 0; 
}
