// Autor: Antonio David Subires Santana
// Estudiante del Grado en Física por la Universdad de Granada
// Programa que calcula la propagación de caminantes aleatorios en hypercubes



#include "complex.h"
#include <stdio.h>
#include <math.h>
#include <gsl/gsl_rng.h>

int main(){
    int t,i,j,dimension,walk,x[50000];
    double norma, QRW[51],AUX[51],p,media,desviacion;
    //Inicialización de la gsl
    const gsl_rng_type * T;
    gsl_rng * r;
    
    gsl_rng_env_setup();
    gsl_rng_default_seed=3; 

    T = gsl_rng_mt19937; 
    r = gsl_rng_alloc (T);
    gsl_rng_set(r,5);

    //Finish gsl_rng
    
    FILE *fdatos,*fdatos2;
    fdatos=fopen("Hitting2.txt","w");
    fdatos2=fopen("Hittingpromedio3.txt","w");
    for(dimension=2;dimension<=20;dimension++){
        for(i=0;i<50000;i++){
            x[i]=0;
        }
        desviacion=0.;
        media=0.;
        printf("dimension=%i\n",dimension);
        for(walk=0;walk<50000;walk++){
            
            for(i=0;i<dimension+1;i++){
                    QRW[i]=0;
                    AUX[i]=0;
            }
            
            QRW[0]=1;
            AUX[0]=1;
            t=0;
            
            while(t>=0){
                for(i=0;i<dimension+1;i++){
                    AUX[i]=QRW[i];
                }
                // Propago al caminante en el hipercubo
                for(i=0;i<dimension+1;i++){
                    if(QRW[i]!=0){
                        if(i==0){
                            AUX[i+1]=AUX[i+1]+QRW[i];
                        }
                        else if(i==dimension){
                            AUX[i-1]=AUX[i-1]+QRW[i];
                        }
                        else{
                            AUX[i+1]=AUX[i+1]+QRW[i];
                            AUX[i-1]=AUX[i-1]+QRW[i];
                        }
                        AUX[i]=AUX[i]-QRW[i];
                        //AUX[i]=0.0;
                    
                    }
                }
                t=t+1;
                x[walk]=t;
                media=media+1;
                // Imprimo el auxiliar sobre el caminante
                for(i=0;i<dimension+1;i++){
                    QRW[i]=AUX[i];
                }
                norma=0.0;
                // Normalizo
                for(i=0;i<dimension+1;i++){
                    norma=norma+QRW[i]*QRW[i];
                }
                for(i=0;i<dimension+1;i++){
                    QRW[i]=QRW[i]/sqrt(norma);
                }
                
                /*
                printf("\n\nmedia = %i\n\n",t);
                for(i=0;i<dimension+1;i++){
                    printf("%lf\t%i\n",QRW[i],i);
                }
                //*/
                // Aquí mido
                p=gsl_rng_uniform (r);
                //printf("\nprobabilidad=%lf\tprobabilidad objetivo=%lf\n",p,QRW[dimension]);
                if(p<QRW[dimension]){
                    fprintf(fdatos,"%i\t%i\n",dimension,t);
                    break;
                }
                else{
                    //printf("renormalizo\n");
                    QRW[dimension]=0;
                    norma=0.0;
                        for(i=0;i<dimension+1;i++){
                            norma=norma+QRW[i];
                        }
                        for(i=0;i<dimension+1;i++){
                            QRW[i]=QRW[i]/sqrt(norma);
                        }
                }
            }
        }
        media=media/50000.;
        for(i=0;i<50000;i++){
            desviacion+=(x[i]-media)*(x[i]-media);
        }
        desviacion=desviacion/50000.;
        desviacion=sqrt(desviacion);
        fprintf(fdatos2,"%i\t%lf\t0\t%lf\n",dimension,media,desviacion);
    }
    fclose(fdatos);
    fclose(fdatos2);
    gsl_rng_free (r);
    return 0;
}
