//This code simulates the mixed species Viscek Model - no repulsion, no turning angle, scalar noise and angle based implementation.

// Angles are chosen over [-PI,PI] and are maintained so throughout.
// The size of spheres in the animations are not to scale.
#include<stdio.h>
#include<iostream>
#include<stdlib.h>
#include<math.h>
#include<time.h>
#include<unistd.h>
using namespace std;

#define PI 3.14159


const int N=100; // Number of individuals in each flock


float L,r=1,delta_r,v=.03,delta_t=1,eta;


float simulate(float x1[2][N],float y1[2][N],float theta1[2][N],float x2[2][N],float y2[2][N],float theta2[2][N],long int T, FILE *fp); //Runs simulation
void initialize(float x1[2][N],float y1[2][N],float theta1[2][N],float x2[2][N],float y2[2][N],float theta2[2][N]); // Gives random initial values to position and direction of each particle
void update_pos(float x1[2][N],float y1[2][N],float theta1[2][N],float x2[2][N],float y2[2][N],float theta2[2][N]); // Updates position of each particle
void update_vel(float x1[2][N],float y1[2][N],float theta1[2][N],float x2[2][N],float y2[2][N],float theta2[2][N]); // Updates velocity of each particle
float distance(float x1,float x2,float y1,float y2); // calculates distance along torus between two points
float arctan(float sin,float cos); //generates the angle between [-PI PI] for given sin and cos values
float limit(float x); // limits angle to [-PI PI] for some angle x
void swap(float x1[2][N],float y1[2][N],float theta1[2][N],float x2[2][N],float y2[2][N],float theta2[2][N]); // x[0][i] = x[1][i] and so for other variables
float Orderparameter(float theta1[2][N],float theta2[2][N]); 
float *correlation(float x[2][N],float y[2][N],float theta[2][N],float l,float delta_l);


void timeseriesplot();
void plot(FILE *fp);
void write_to_file(float x[2][N],float y[2][N]);



int main()
{
     float x1[2][N],y1[2][N],theta1[2][N],x2[2][N],y2[2][N],theta2[2][N];
     int n=1;
     float v_a[n],mean,stdev;
     float density;
     char name[100];
     long int T = 5000;
     FILE *data;
     
     
     density = 5.00;
     
     L = sqrt(N/density);
     
     for(eta = 0;eta <= 7;eta+= .25)
     {
	  sprintf(name,"Correlations/Density_%.2fNoise_%.2f.csv",density,eta);
	  printf("%s\n",name);
	  data = fopen(name,"w");
	  
	  v_a[0] = simulate(x1,y1,theta1,x2,y2,theta2,T,data); 
	  
	  fclose(data);
	  
     }
     
     
     return 0;
}


float simulate(float x1[2][N],float y1[2][N],float theta1[2][N],float x2[2][N],float y2[2][N],float theta2[2][N],long int T,FILE *fp)
{
     
     // FILE *gnupipe;
     //  gnupipe = popen("gnuplot -persistent","w");
     
     
     initialize(x1,y1,theta1,x2,y2,theta2);   
     long int t;
     int i;
     for(t=0;t<T;t++)
     {
	  update_pos(x1,y1,theta1,x2,y2,theta2);
	  update_vel(x1,y1,theta1,x2,y2,theta2);
	  
	  // fprintf(fp,"%ld\t%f\n",t,Orderparameter(theta));
	  // write_to_file(x,y);
	  
	  //plot(gnupipe);
	  
	  swap(x1,y1,theta1,x2,y2,theta2);
     }
     
     // pclose(gnupipe);
     /*
     for(i = 0;i < (N-1);i++)
     {
	  fprintf(fp,"x%d,y%d,theta%d,",i,i,i); 
	  
     }
     fprintf(fp,"x%d,y%d,theta%d\n",i,i,i);
     
     for(t=0;t<500;t++)
     {
	  update_pos(x,y,theta);
	  update_vel(x,y,theta);
	  
	  for(i = 0; i < (N-1); i++)
	  {
	       fprintf(fp,"%.2f,%.2f,%.2f,",x[1][i],y[1][i],theta[1][i]);
	       
	  }
	  
	  fprintf(fp,"%.2f,%.2f,%.2f\n",x[1][i],y[1][i],theta[1][i]);
	  
	  
	  swap(x,y,theta);
     }*/
     
     
     return (Orderparameter(theta1,theta2));
}


void initialize(float x1[2][N],float y1[2][N],float theta1[2][N],float x2[2][N],float y2[2][N],float theta2[2][N])
{
     srand(time(NULL));
     
     for(int i=0;i<N;i++)
     {
	  x1[0][i]= ( (float)rand()/(float)RAND_MAX )*L;
	  y1[0][i]= ( ((float)rand())/((float)RAND_MAX) )*L;
	  theta1[0][i] = -PI + ( ((float)rand())/ ((float)RAND_MAX) )*2*PI;   
	  
	  x2[0][i]= ( (float)rand()/(float)RAND_MAX )*L;
	  y2[0][i]= ( ((float)rand())/((float)RAND_MAX) )*L;
	  theta2[0][i] = -PI + ( ((float)rand())/ ((float)RAND_MAX) )*2*PI;  
     }
}

void update_pos(float x1[2][N],float y1[2][N],float theta1[2][N],float x2[2][N],float y2[2][N],float theta2[2][N])
{
     for(int i=0;i<N;i++)
     {
	  x1[1][i] = x1[0][i] + v*cos(theta1[0][i])*delta_t;
	  y1[1][i] = y1[0][i] + v*sin(theta1[0][i])*delta_t;
	  
	  x2[1][i] = x2[0][i] + v*cos(theta2[0][i])*delta_t;
	  y2[1][i] = y2[0][i] + v*sin(theta2[0][i])*delta_t;
	  
	  if(x1[1][i] > L)
	  {
	       x1[1][i] -= L;
	  }
	  
	  if(x1[1][i] < 0)
	  {
	       x1[1][i] += L;
	  }
	  
	  if(y1[1][i] > L)
	  {
	       y1[1][i]-= L;
	  }
	  
	  if(y1[1][i] < 0)
	  {
	       y1[1][i] += L;
	  }
	  
	  if(x2[1][i] > L)
	  {
	       x2[1][i] -= L;
	  }
	  
	  if(x2[1][i] < 0)
	  {
	       x2[1][i] += L;
	  }
	  
	  if(y2[1][i] > L)
	  {
	       y2[1][i]-= L;
	  }
	  
	  if(y2[1][i] < 0)
	  {
	       y2[1][i] += L;
	  }
	  
     }  
}
void update_vel(float x1[2][N],float y1[2][N],float theta1[2][N],float x2[2][N],float y2[2][N],float theta2[2][N])
{
     for(int i=0; i<N;i++)
     {
	  float sumsin1,sumsin2;
	  float sumcos1,sumcos2;
	  float count1,count2;
	  float d_same,d_diff;
	  sumsin1 = 0;
	  sumcos1 = 0;
	  count1 = 0;
	  
	  sumsin1 = 0;
	  sumcos1 = 0;
	  count1 = 0;
	  
	  for(int j =0 ; j < N; j++)
	  {
	       d_same = distance(x1[0][i],x1[0][j],y1[0][i],y1[0][j]);
	       
	       if(d_same < r)
	       {
		    sumsin1 += sin( theta1[0][j]);
		    sumcos1 += cos( theta1[0][j]);
		    count1 += 1;
	       }
	       
	       d_diff = distance(x1[0][i],x2[0][j],y1[0][i],y2[0][j]);
	       
	       if(d_diff < (r-delta_r))
	       {
		    sumsin2 += sin( theta2[0][j]);
		    sumcos2 += cos( theta2[0][j]);
		    count2 += 1;
	       }
	  }
	  
	  float sinavg;
	  float cosavg;
	  float thetaavg;
	  float delta_theta; 
	  
	 
	  sinavg = (sumsin1+sumsin2)/(count1 + count2);
	  cosavg = (sumcos1 + sumcos2)/(count1 + count2);
	  thetaavg = arctan(sinavg,cosavg); 
	  delta_theta = eta*( (float)rand()/(float)RAND_MAX ) - eta/2;
	  
	  theta1[1][i] = thetaavg + delta_theta;
	  
	  theta1[1][i] = limit(theta1[1][i]);
     }
     
     
     for( int i=0; i<N;i++)
     {
	  float sumsin1,sumsin2;
	  float sumcos1,sumcos2;
	  float count1,count2;
	  float d_same,d_diff;
	  sumsin1 = 0;
	  sumcos1 = 0;
	  count1 = 0;
	  
	  sumsin1 = 0;
	  sumcos1 = 0;
	  count1 = 0;
	  
	  for(int j =0 ; j < N; j++)
	  {
	       d_same = distance(x2[0][i],x2[0][j],y2[0][i],y2[0][j]);
	       
	       if(d_same < r)
	       {
		    sumsin1 += sin( theta2[0][j]);
		    sumcos1 += cos( theta2[0][j]);
		    count1 += 1;
	       }
	       
	       d_diff = distance(x2[0][i],x1[0][j],y2[0][i],y1[0][j]);
	       
	       if(d_diff < (r-delta_r))
	       {
		    sumsin2 += sin( theta1[0][j]);
		    sumcos2 += cos( theta1[0][j]);
		    count2 += 1;
	       }
	  }
	  
	  float sinavg;
	  float cosavg;
	  float thetaavg;
	  float delta_theta; 
	  
	  
	  sinavg = (sumsin1+sumsin2)/(count1 + count2);
	  cosavg = (sumcos1 + sumcos2)/(count1 + count2);
	  thetaavg = arctan(sinavg,cosavg); 
	  delta_theta = eta*( (float)rand()/(float)RAND_MAX ) - eta/2;
	  
	  theta2[1][i] = thetaavg + delta_theta;
	  
	  theta2[1][i] = limit(theta2[1][i]);
     }

}

void swap(float x1[2][N],float y1[2][N],float theta1[2][N],float x2[2][N],float y2[2][N],float theta2[2][N])
{
     for(int i=0;i<N;i++)
     {
	  x1[0][i] = x1[1][i];
	  y1[0][i] = y1[1][i];
	  theta1[0][i] = theta1[1][i]; 
	  
	  x2[0][i] = x2[1][i];
	  y2[0][i] = y2[1][i];
	  theta2[0][i] = theta2[1][i]; 
     }
}

float Orderparameter(float theta1[2][N],float theta2[2][N])
{
     float sumcos,sumsin;
     sumcos =0;
     sumsin =0;
     
     for(int i=0;i<N;i++)
     {
	  sumsin += sin(theta1[1][i]) + sin(theta2[1][i]);
	  sumcos += cos(theta1[1][i]) + cos(theta2[1][i]) ;
     }
     
     return((sqrt(sumsin*sumsin + sumcos*sumcos))/(2*N));
}

float distance(float x1,float x2,float y1,float y2)
{
     float dx,dy;
     dx = fabs(x1-x2);
     dy = fabs(y1-y2);
     
     if(dx > L/2)
	  dx = L-dx;
     if( dy > L/2)
	  dy = L-dy;
     
     return( sqrt(dx*dx + dy*dy) );
}

float arctan(float sin,float cos)
{
     if(sin > 0  && cos > 0) return(atan(sin/cos));
     if(sin < 0  && cos > 0) return(atan(sin/cos));
     if(sin > 0  && cos < 0) return(PI + atan(sin/cos));
     if(sin < 0  && cos < 0) return(-PI + atan(sin/cos));
     
     if(sin == 0 && cos > 0) return(0);
     if(sin == 0 && cos < 0) return(PI);
     if(sin > 0 && cos == 0) return(PI/2);
     if(sin < 0 && cos == 0) return(-PI/2);
     
}

float limit(float x)
{
     if(x <= PI && x >= -PI) return(x);
     
     while(!(x <=PI && x >= -PI))
     {
	  
	  if(x > PI && x <= 2*PI) x = x - 2*PI;
	  else if(x < -PI && x >= -2*PI) x = x + 2*PI;
	  else
	  {
	       if(x > 2*PI)
	       {
		    while(x >= 2*PI) x = x - 2*PI;
	       }
	       
	       if(x < 0)
	       {
		    while(x <= 0) x = x + 2*PI;
	       }
	  }  
     }
     
     return x;
     
}

void write_to_file(float x[2][N],float y[2][N])
{
     FILE *fp;
     
     fp = fopen("data1.temp","w");
     
     for(int i=0;i<N;i++)
     {
	  fprintf(fp,"%f\t%f\t",x[1][i],y[1][i]);
     }
     
     fprintf(fp,"\n");
     
     fclose(fp);
}

void plot(FILE *fp)
{
     int i=1;
     double l=1;
     float x;
     x= 4.9*10/L;
     
     fprintf(fp,"set terminal wxt size 600,600\n");
     fprintf(fp,"set size square 1,1\n");
     fprintf(fp,"set key off\n");
     fprintf(fp,"unset xtic \n");
     fprintf(fp,"unset ytic \n");
     fprintf(fp,"set grid xtic\n");
     fprintf(fp,"set grid ytic\n");
     fprintf(fp,"set xrange [0:%lf]\n",L);
     fprintf(fp,"set yrange [0:%lf]\n",L);
     
     fprintf(fp,"plot \'data1.temp\' u %d:%d w points pointtype 6 ps %.1f,",1,2,x);
     for(i=3;i<(2*N-1);i+=2)
     {
	  fprintf(fp,"\'data1.temp\' u %d:%d w points pointtype 6 ps %.1f,",i,i+1,x);  
     }
     
     
     fprintf(fp,"\'data1.temp\' u %d:%d w points pointtype 6 ps %.1f\n",i,i+1,x);
     
     
     usleep(10000);
     
     //pclose(gnuplotPipe);
     
}

float *correlation(float x[2][N],float y[2][N],float theta[2][N],float l = 0,float delta_l = L/100)
{
     float* p;
     p = (float *)malloc( sizeof(float)*(int)( (L-l)/delta_l ) );
     int k = 0;
     while(l<=L && k<= sizeof(float)*(int)( (L-l)/delta_l ))
     {
	  int count = 0;
	  float sum = 0;
	  for(int i=0;i<N;i++)
	  {
	       for(int j=0;j<N;j++)
	       {
		    float d = distance(x[1][i],x[1][j],y[1][i],y[1][j]);
		    if(d > l && d < (l+delta_l) )
		    {
			 count++;
			 sum += cos(theta[1][i])*cos(theta[1][j]) + sin(theta[1][i])*sin(theta[1][j]);
		    }
	       }
	  }
	  
	  p[k] = sum/count;
	  k++;
	  l += delta_l;
     }
     
     return p;
     
}

void timeseriesplot()
{
     char name[50],output[50];
     FILE *fp;
     fp = popen("gnuplot","w");
     
     
     
     for(float density = 1;density<=5;density+=0.5)
     {  
	  for(float eta =0;eta<=5;eta+=.5)
	  {
	       for(long int T = 5000; T<=15000; T+= 5000)
	       {
		    
		    FILE *fp;
		    fp = popen("gnuplot","w");
		    
		    
		    
		    
		    fprintf(fp,"set terminal png\n");
		    fprintf(fp,"set xrange [0:%ld]\n",T);
		    fprintf(fp,"set yrange [0:1.2]\n");
		    fprintf(fp,"set xlabel 'T'\n");
		    fprintf(fp,"set ylabel 'Order Parameter'\n");
		    fprintf(fp,"set title \'D_%.1f,eta_%.1f,T_%ld\'",density,eta,T);
		    sprintf(output,"time_%.1f_%.1f_%ld.png",density,eta,T);
		    fprintf(fp,"set output \"%s\"",output);
		    sprintf(name,"time_%.1f_%.1f_%ld_0.dat",density,eta,T);
		    
		    fprintf(fp,"plot \'%s\' u 1:2,",name);
		    for(int i=1;i<9;i++)
		    {
			 sprintf(name,"time_%.1f_%.1f_%ld_%d.dat",density,eta,T,i);
			 fprintf(fp,"\'%s\' u 1:2,",name);
			 
		    }
		    sprintf(name,"time_%.1f_%.1f_%ld_9.dat",density,eta,T);
		    fprintf(fp,"\'%s\' u 1:2 \n",name);
		    
		    fclose(fp);
		    
	       } 
	  }
	  
	  
     }
     
     
     
}


