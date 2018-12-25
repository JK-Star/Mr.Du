/* 
 * auth = peic 
 * 对矩阵乘法进行优化， 对矩阵进行转置 
 */  
  
#include <stdio.h>  
#include <stdlib.h>  
#include <time.h>
#include <stddef.h>  
#include <x86intrin.h>  
#include <avx2intrin.h>
#include<omp.h>

#define N 1024  
  
void computer_pi_native(size_t dt)
{
    
    double begin_time1 = (double)clock() ; 
    double pi=0.0;
    double delta=1.0/dt;
    for(size_t i=0;i<dt;i++)
    {
        double x=(double)i/dt;
        pi+=delta/(1.0+x*x);
    }
    double end_time1 = (double)clock() ; 
    double cost_time1 = end_time1 - begin_time1;  
    printf("pi=%.10lf\t",pi*4.0); 
    printf("serial cost time: %fus\n", cost_time1);   
}  
void computer_pi_avx_f(size_t dt)
{
    double begin_time1 = (double)clock() ; 
    float pi=0.0;
    float delta=1.0f/dt;
    __m256 ymm0,ymm1,ymm2,ymm3,ymm4;
    ymm0=_mm256_set1_ps(1.0);
    ymm1=_mm256_set1_ps(delta);
    ymm2=_mm256_set_ps(delta*7,delta*6,delta*5,delta*4,delta*3,delta*2,delta,0.0);
    ymm4=_mm256_setzero_ps();
    for(int i=0;i<=dt-8;i+=8)
    {
        ymm3=_mm256_set1_ps(i*delta);
        ymm3=_mm256_add_ps(ymm3,ymm2);
        ymm3=_mm256_mul_ps(ymm3,ymm3);
        ymm3=_mm256_add_ps(ymm0,ymm3);
        ymm3=_mm256_div_ps(ymm1,ymm3);
        ymm4=_mm256_add_ps(ymm4,ymm3);
    }

    float tmp[8] __attribute__((aligned(16)));
    _mm256_store_ps(tmp,ymm4);
    pi+=tmp[0]+tmp[1]+tmp[2]+tmp[3]+tmp[4]+tmp[5]+tmp[6]+tmp[7];
    double end_time1 = (double)clock() ; 
    double cost_time1 = end_time1 - begin_time1;  
    printf("pi=%.10lf\t",pi*4.0); 
    printf("avx-float cost time: %fus\n", cost_time1);  
    
} 
void computer_pi_avx_d(size_t dt)
{
    double begin_time1 = (double)clock() ;  
    double pi=0.0;
    double delta=1.0/dt;
    __m256d ymm0,ymm1,ymm2,ymm3,ymm4;
    ymm0=_mm256_set1_pd(1.0);
    ymm1=_mm256_set1_pd(delta);
    ymm2=_mm256_set_pd(delta*3,delta*2,delta,0.0);
    ymm4=_mm256_setzero_pd();
    for(int i=0;i<=dt-4;i+=4)
    {
        ymm3=_mm256_set1_pd(i*delta);
        ymm3=_mm256_add_pd(ymm3,ymm2);
        ymm3=_mm256_mul_pd(ymm3,ymm3);
        ymm3=_mm256_add_pd(ymm0,ymm3);
        ymm3=_mm256_div_pd(ymm1,ymm3);
        ymm4=_mm256_add_pd(ymm4,ymm3);
    }

    double tmp[4] __attribute__((aligned(32)));
    _mm256_store_pd(tmp,ymm4);
    pi+=tmp[0]+tmp[1]+tmp[2]+tmp[3];
    double end_time1 = (double)clock() ; 
    double cost_time1 = end_time1 - begin_time1;  
    printf("pi=%.10lf\t",pi*4.0); 
    printf("avx-double cost time: %fus\n", cost_time1);   
    
}  
void computer_pi_sse_f(size_t dt)
{
   double begin_time1 = (double)clock() ;  
    float pi=0.0;
    float delta=1.0f/dt;
    __m128 xmm0,xmm1,xmm2,xmm3,xmm4;
    xmm0=_mm_set1_ps(1.0);
    xmm1=_mm_set1_ps(delta);
    xmm2=_mm_set_ps(delta*3,delta*2,delta,0.0);
    xmm4=_mm_setzero_ps();
    for(int i=0;i<=dt-4;i+=4)
    {
        xmm3=_mm_set1_ps(i*delta);
        xmm3=_mm_add_ps(xmm3,xmm2);
        xmm3=_mm_mul_ps(xmm3,xmm3);
        xmm3=_mm_add_ps(xmm0,xmm3);
        xmm3=_mm_div_ps(xmm1,xmm3);
        xmm4=_mm_add_ps(xmm4,xmm3);
    }

    float tmp[4] __attribute__((aligned(16)));
    _mm_store_ps(tmp,xmm4);
    pi+=tmp[0]+tmp[1]+tmp[2]+tmp[3];
    double end_time1 = (double)clock() ; 
    double cost_time1 = end_time1 - begin_time1;  
    printf("pi=%.10lf\t",pi*4.0); 
    printf("sse-float cost time: %fus\n", cost_time1);   
     
}
void computer_pi_sse_d(size_t dt)
{
   double begin_time1 = (double)clock() ; 
    double pi=0.0;
    double delta=1.0/dt;
    __m128d xmm0,xmm1,xmm2,xmm3,xmm4;
    xmm0=_mm_set1_pd(1.0);
    xmm1=_mm_set1_pd(delta);
    xmm2=_mm_set_pd(delta,0.0);
    xmm4=_mm_setzero_pd();
    for(int i=0;i<=dt-2;i+=2)
    {
        xmm3=_mm_set1_pd(i*delta);
        xmm3=_mm_add_pd(xmm3,xmm2);
        xmm3=_mm_mul_pd(xmm3,xmm3);
        xmm3=_mm_add_pd(xmm0,xmm3);
        xmm3=_mm_div_pd(xmm1,xmm3);
        xmm4=_mm_add_pd(xmm4,xmm3);
    }

    double tmp[2] __attribute__((aligned(32)));
    _mm_store_pd(tmp,xmm4);
    pi+=tmp[0]+tmp[1];
    double end_time1 = (double)clock() ; 
    double cost_time1 = end_time1 - begin_time1;  
    printf("pi=%.10lf\t",pi*4.0); 
    printf("sse-double cost time: %fus\n", cost_time1);  
     
}
void computer_pi_native_o(size_t dt)
{
    double begin_time1 = (double)clock() ;
    double pi=0.0f;
    double temp;
    #pragma omp parallel for simd default(none) shared(dt) reduction(+:pi)
    for(int i=0;i<dt;i++)
    {
        float temp=(i+0.5f)/dt;
        
        temp=4/(1+temp*temp);
        pi+=temp;
        
    }
    double end_time1 = (double)clock() ; 
    double cost_time1 = end_time1 - begin_time1;  
    printf("pi=%.10lf\t",pi/dt); 
    printf("openmp native cost time: %fus\n", cost_time1); 
  
    
} 
void computer_pi_native_openmp(size_t dt)
{
    
  
    double begin_time1 = omp_get_wtime() ; 
    double pi=0.0;
    double delta=1.0/dt;
 
    #pragma omp parallel for default(none) shared(dt) reduction(+:pi)
    for(int i=0;i<dt;i++)
    {
        //printf("%d",omp_get_num_threads());
        float temp=(i+0.5f)/dt;
        
        temp=4/(1+temp*temp);
        pi+=temp;
        
    }
    double end_time1 = omp_get_wtime() ; 
    double cost_time1 = end_time1 - begin_time1;  
    printf("pi=%.10lf\t",pi/dt); 
    
    printf("openmp cost time: %lfus\n", cost_time1*1000000);  
  
    
}


void computer_pi_native_openmp_simd(size_t dt)
{
    
  
    double begin_time1 = omp_get_wtime() ; 
    double pi=0.0;
    double delta=1.0/dt;
    double tmp_tt;
    __m256d ymm0,ymm1,ymm2,ymm3,ymm4;
    
 
    #pragma omp parallel for reduction(+:pi)  
    for(int i=0;i<dt;i+=dt/4)
        {   ymm0=_mm256_set1_pd(1.0);
            ymm1=_mm256_set1_pd(delta);
            ymm2=_mm256_set_pd(delta*3,delta*2,delta,0.0);
            ymm4=_mm256_setzero_pd();//firstprivate(ymm0,ymm1,ymm2,ymm3,ymm4) 
            #pragma omp firstprivate(ymm0,ymm1,ymm2,ymm3,ymm4)
            for(int j=i;j<i+dt/4-4;j+=4)
            {
            ymm3=_mm256_set1_pd(j*delta);
            ymm3=_mm256_add_pd(ymm3,ymm2);
            ymm3=_mm256_mul_pd(ymm3,ymm3);
            ymm3=_mm256_add_pd(ymm0,ymm3);
            ymm3=_mm256_div_pd(ymm1,ymm3);
            ymm4=_mm256_add_pd(ymm4,ymm3);
            //printf("%d %d %lf %d\n",i,j,pi,omp_get_thread_num());
            } 
            double tmp[4] __attribute__((aligned(32)));
            _mm256_store_pd(tmp,ymm4);
            tmp_tt=tmp[0]+tmp[1]+tmp[2]+tmp[3];
            pi+=tmp_tt;               
        }
    double end_time1 = omp_get_wtime() ; 
    double cost_time1 = end_time1 - begin_time1;  
    printf("pi=%.10lf\t",pi*4); 
    
    printf("openmp_simd cost time: %lfus\n", cost_time1*1000000);        
}
void compute_pi_omp_avx_double_cycle(int dt){

    double begin_time1 = omp_get_wtime() ; 
	double pi=0.0;
	double delta=1.0/dt;
	__m256d yum0,yum1,yum2,yum3,yum4,yum5,yum6;
	yum0 = _mm256_set1_pd(1.0);
	yum1 = _mm256_set1_pd(delta);
	yum2 = _mm256_set_pd(delta*3,delta*2,delta,0.0);
	yum4 = _mm256_setzero_pd();
	yum5 = _mm256_setzero_pd();
	for(int i=0;i<=dt-8;i+=8){
		yum3 = _mm256_set1_pd(i*delta);
		yum3 = _mm256_add_pd(yum3,yum2);
		yum3 = _mm256_mul_pd(yum3,yum3);
		yum3 = _mm256_add_pd(yum0,yum3);
		yum3 = _mm256_div_pd(yum1,yum3);
		yum4 = _mm256_add_pd(yum3,yum4);
		
		yum6 = _mm256_set1_pd((i+4)*delta);
		yum6 = _mm256_add_pd(yum6,yum2);
		yum6 = _mm256_mul_pd(yum6,yum6);
		yum6 = _mm256_add_pd(yum0,yum6);
		yum6 = _mm256_div_pd(yum1,yum6);
		yum5 = _mm256_add_pd(yum5,yum6);
	}
	yum4 = _mm256_add_pd(yum4,yum5);
	double tmp[4] __attribute__((aligned(32)));
	_mm256_store_pd(tmp,yum4);
	pi = tmp[0]+tmp[1]+tmp[2]+tmp[3];
    double end_time1 = omp_get_wtime() ; 
    double cost_time1 = end_time1 - begin_time1;  
    printf("pi=%.10lf\t",pi*4); 
    
    printf("cycle cost time: %lfus\n", cost_time1*1000000);   
}
int main()  
{  size_t dt=1280000;   
   computer_pi_native(dt);
   computer_pi_avx_f(dt);
   computer_pi_avx_d(dt);
   computer_pi_sse_f(dt);
   computer_pi_sse_d(dt);
   compute_pi_omp_avx_double_cycle(dt);
  // computer_pi_native_o(dt);
   computer_pi_native_openmp(dt);
   computer_pi_native_openmp_simd(dt);
   
  

    return 0; 
}