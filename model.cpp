//This code simulates the original Viscek Model as is in the 1995 PRL - no repulsion, no turning angle, scalar noise and angle based implementation
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


const int N=100;


float L,r=1,v=.03,delta_t=1,eta;


float simulate(float x[2][N],float y[2][N],float theta[2][N],long int T, FILE *fp); //Runs simulation
void initialize(float x[2][N],float y[2][N],float theta[2][N]); // Gives random initial values to position and direction of each particle
void update_pos(float x[2][N],float y[2][N],float theta[2][N]); // Updates position of each particle
void update_vel(float x[2][N],float y[2][N],float theta[2][N]); // Updates velocity of each particle
float distance(float x1,float x2,float y1,float y2); // calculates distance along torus between two points
float arctan(float sin,float cos); //generates the angle between [-PI PI] for given sin and cos values
float limit(float x); // limits angle to [-PI PI] for some angle x
void swap(float x[2][N],float y[2][N],float theta[2][N]); // x[0][i] = x[1][i] and so for other variables
float Orderparameter(float theta[2][N]); 
float *correlation(float x[2][N],float y[2][N],float theta[2][N],float l,float delta_l);


void timeseriesplot();
void plot(FILE *fp);
void write_to_file(float x[2][N],float y[2][N]);



int main()
{
     float x[2][N],y[2][N],theta[2][N];
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
	  
	  v_a[0] = simulate(x,y,theta,T,data); 
	  
	  fclose(data);
	  
     }
     
     
     return 0;
}


float simulate(float x[2][N],float y[2][N],float theta[2][N],long int T,FILE *fp)
{
     
     // FILE *gnupipe;
     //  gnupipe = popen("gnuplot -persistent","w");
     
     
     initialize(x,y,theta);   
     long int t;
     int i;
     for(t=0;t<T;t++)
     {
	  update_pos(x,y,theta);
	  update_vel(x,y,theta);
	  
	  // fprintf(fp,"%ld\t%f\n",t,Orderparameter(theta));
	  // write_to_file(x,y);
	  
	  //plot(gnupipe);
	  
	  swap(x,y,theta);
     }
     
     // pclose(gnupipe);
     
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
     }
     
     
     return (Orderparameter(theta));
}


void initialize(float x[2][N],float y[2][N],float theta[2][N])
{
     srand(time(NULL));
     
     for(int i=0;i<N;i++)
     {
	  x[0][i]= ( (float)rand()/(float)RAND_MAX )*L;
	  y[0][i]= ( ((float)rand())/((float)RAND_MAX) )*L;
	  theta[0][i] = -PI + ( ((float)rand())/ ((float)RAND_MAX) )*2*PI;   
     }
}

void update_pos(float x[2][N],float y[2][N],float theta[2][N])
{
     for(int i=0;i<N;i++)
     {
	  x[1][i] = x[0][i] + v*cos(theta[0][i])*delta_t;
	  y[1][i] = y[0][i] + v*sin(theta[0][i])*delta_t;
	  
	  if(x[1][i] > L)
	  {
	       x[1][i] -= L;
	  }
	  
	  if(x[1][i] < 0)
	  {
	       x[1][i] += L;
	  }
	  
	  if(y[1][i] > L)
	  {
	       y[1][i]-= L;
	  }
	  
	  if(y[1][i] < 0)
	  {
	       y[1][i] += L;
	  }
	  
     }  
}
void update_vel(float x[2][N],float y[2][N],float theta[2][N])
{
     for(int i=0; i<N;i++)
     {
	  float sumsin;
	  float sumcos;
	  float count;
	  float d;
	  sumsin = 0;
	  sumcos = 0;
	  count = 0;
	  
	  for(int j =0 ; j < N; j++)
	  {
	       d = distance(x[0][i],x[0][j],y[0][i],y[0][j]);
	       
	       if(d < r)
	       {
		    sumsin += sin( theta[0][j]);
		    sumcos += cos( theta[0][j]);
		    count += 1;
	       }
	  }
	  
	  float sinavg;
	  float cosavg;
	  float thetaavg;
	  float delta_theta; 
	  
	  if(count == 0) {
	       count = 1;
	       sumsin = sin(theta[0][i]);
	       sumcos = cos(theta[0][i]);
	  }
	  
	  sinavg = sumsin/count;
	  cosavg = sumcos/count;
	  thetaavg = arctan(sinavg,cosavg); 
	  delta_theta = eta*( (float)rand()/(float)RAND_MAX ) - eta/2;
	  
	  theta[1][i] = thetaavg + delta_theta;
	  
	  theta[1][i] = limit(theta[1][i]);
     }
}

void swap(float x[2][N],float y[2][N],float theta[2][N])
{
     for(int i=0;i<N;i++)
     {
	  x[0][i] = x[1][i];
	  y[0][i] = y[1][i];
	  theta[0][i] = theta[1][i]; 
     }
}

float Orderparameter(float theta[2][N])
{
     float sumcos,sumsin;
     sumcos =0;
     sumsin =0;
     
     for(int i=0;i<N;i++)
     {
	  sumsin += sin(theta[1][i]);
	  sumcos += cos(theta[1][i]);
     }
     
     return((sqrt(sumsin*sumsin + sumcos*sumcos))/N);
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


