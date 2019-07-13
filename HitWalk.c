// Autor: Antonio David Subires Santana
// Estudiante del Grado en Física por la Universdad de Granada
// Caminante aleatorio cuántico en 1 dimensión
// Todas las funciones que se utilizan son complejas
// La moneda se define al principio C (coin) y también la función de onda inicial (QRW)
// QRW representa el estado de mi red (monodimensional) la primera columna se refiere a los estados de spin +
// mientras que la segunda a los de spin -


#include "complex.h"
#include <stdio.h>
#include <math.h>
#include <gsl/gsl_rng.h>


int main()
{
    int L,medida;
    int k, j;
    double p, probabibilidad;
    double NORMA;
    fcomplex C[2][2], QRW[2002][2], AUX[2002][2];
    
    // I define the COIN
    // we use the Hadamard's operator
    
    C[0][0]=Complex(1.0/sqrt(2),0);
    C[0][1]=Complex(1.0/sqrt(2),0);
    C[1][0]=Complex(1.0/sqrt(2),0);
    C[1][1]=Complex(-1.0/sqrt(2),0);
    
    FILE *fprobabilidad;
    fprobabilidad=fopen("QuantumHitting.txt","w");
    
    //Inicialización de la gsl
    const gsl_rng_type * T;
    gsl_rng * r;
    
    gsl_rng_env_setup();
    gsl_rng_default_seed=3; 

    T = gsl_rng_mt19937; 
    r = gsl_rng_alloc (T);
    gsl_rng_set(r,5);

    //Finish gsl_rng
    L=1000;
    
    for(medida=10;medida<500;medida++){
  /*
    if(L%2==0){
        medida=L/2;
    }
    else{
        medida=(L+1)/2;
    }
   // */
    // The place on I want to measure the wave
    
    //I initialize the matrix QRW and AUX matrix as null
    for(k=0;k<2002;k++){
        for(j=0;j<2;j++){
            QRW[k][j]=Complex(0,0);
            AUX[k][j]=Complex(0,0);
        }
    }
    // Estado inicial ANTISIMÉTRICO
    
    /*
    QRW[L][1]=Complex(1.0,0);
    AUX[L][1]=Complex(1.0,0);
    //*/
    
    
    // Estado inicial SIMÉTRICO
    ///*
    QRW[L][0]=Complex(1.0/sqrt(2),0);
    AUX[L][0]=Complex(1.0/sqrt(2),0);
    
    QRW[L][1]=Complex(0,1.0/sqrt(2));
    AUX[L][1]=Complex(0,1.0/sqrt(2));
    //*/
    // Let's go write the algorithm
    int t;
    t=1;
    while(t<(L+1)){
        for(k=0;k<(2*L+1);k++){
            // Now, we must test if the spin is + or -
            // First +
                if(Cabs(QRW[k][0])!=0){
                if(k>0){
                    AUX[k-1][1]=Cadd(AUX[k-1][1],Cmul(C[0][1],QRW[k][0]));
                }
                if(k<2*L){
                    
                    AUX[k+1][0]=Cadd(AUX[k+1][0],Cmul(C[0][0],QRW[k][0]));
                }
                AUX[k][0]=Complex(0,0);
                }
            
            // Now -
                if(Cabs(QRW[k][1])!=0){
                if(k>0){
                    AUX[k-1][1]=Cadd(AUX[k-1][1],Cmul(C[1][1],QRW[k][1]));
                }
                if(k<2*L){
                    AUX[k+1][0]=Cadd(AUX[k+1][0],Cmul(C[1][0],QRW[k][1]));                    
                }
                AUX[k][1]=Complex(0,0);
                }
            }

        for(k=0;k<(2*L+1);k++){
            for(j=0;j<2;j++){
                QRW[k][j]=AUX[k][j];
            }
        }
    // Hacemos una medida y vemos qué pasa.
    
        t=t+1;
        p=gsl_rng_uniform (r);
        probabibilidad=Cabs(QRW[L-1+medida][0])*Cabs(QRW[L-1+medida][0])+Cabs(QRW[L-1+medida][1])*Cabs(QRW[L-1+medida][1]);
        if(p<probabibilidad){
            fprintf(fprobabilidad,"%i\t%i\n",t, medida+1);
            break;
        }
        p=gsl_rng_uniform (r);
        probabibilidad=Cabs(QRW[L-1-medida][0])*Cabs(QRW[L-1-medida][0])+Cabs(QRW[L-1-medida][1])*Cabs(QRW[L-1-medida][1]);
        if(p<probabibilidad){
            fprintf(fprobabilidad,"%i\t%i\n",t, medida);
            break;
        }
    
        else{
            NORMA=0.0;
            QRW[L-1+medida][1]=Complex(0,0);
            QRW[L-1+medida][0]=Complex(0,0);
            QRW[L-1-medida][1]=Complex(0,0);
            QRW[L-1-medida][0]=Complex(0,0);
            for(k=0;k<(2*L+1);k++){
                for(j=0;j<2;j++){
                    NORMA = NORMA + (Cabs(QRW[k][j])*Cabs(QRW[k][j]));
                }
            }
            NORMA=sqrt(NORMA);

    
            for(k=0;k<(2*L+1);k++){
                for(j=0;j<2;j++){
                    QRW[k][j]=RCmul(1.0/(1.0*NORMA),QRW[k][j]);
                }
            }
        }
    }
    }

    fclose(fprobabilidad);
    gsl_rng_free (r);
    
    
    return 0;
}
