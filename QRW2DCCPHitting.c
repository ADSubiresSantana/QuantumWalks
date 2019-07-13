// Autor: Antonio David Subires Santana
// Estudiante del Grado en Física por la Universdad de Granada
//Caminante aleatorio cuántico en 2 dimensiones
// Incluye condiciones periódicas de contorno

// Las direcciones de propagación del caminante son diagonales




#include "complex.h"
#include <stdio.h>
#include <math.h>
#include <gsl/gsl_rng.h>


int main()
{
    // Una fila tiene 2L+1 elementos
    // Como C inicializa las filas y columnas en 0, un índice va a ir desde 0 hasta 2L (incluido)
    // L es el centro de una fila o columna
    // Hay un total de 4L^2 elementos
    int L,k,j,i,x,y,M1,M2,M3,M4,M5,M6,M7,M8,caminante;
    // Para trabajar con la red, en lugar de usar una matriz uso un vector de 4L^2 elementos
    // La primera fila (0) llegará hasta 2L el siguiente elemento sería el 2L+1 que se corresponderia en la red con el elemento (0,1)
    // Con esta notación un elemento (i,j) se expresaría como (fila*(2*L+1)+columna)
    // Llamaremos i=fila*(2*L+1)
    fcomplex C[4][4], Hadamard[4][4], Grover[4][4],QRW[40401][4], AUX[40401][4];
    double NORMA, RaizNorma,probabilidad,p,promediotemporal;

    
    // Ficheros
    FILE *fmapa2d,*fmapa3d,*fhitting;
    fmapa2d=fopen("Datos.txt","w");
    fmapa3d=fopen("Evolucion.txt","w");
    fhitting=fopen("QuantumHitting2DCond29onlyHadamard.txt","w");
    //Inicialización de la gsl
    const gsl_rng_type * T;
    gsl_rng * r;
    
    gsl_rng_env_setup();
    gsl_rng_default_seed=3; 

    T = gsl_rng_mt19937; 
    r = gsl_rng_alloc (T);
    gsl_rng_set(r,5);

    //Finish gsl_rng
                                                                                        // Distintas monedas
    // Operador de Hadamard "separable"
    /*
    for(k=0;k<4;k++){
        for(j=0;j<4;j++){
            C[k][j]=Complex(0,0);
        }
    }
    
    C[0][0]=Complex(1.0/sqrt(2),0);
    C[0][1]=Complex(1.0/sqrt(2),0);
    C[1][0]=Complex(1.0/sqrt(2),0);
    C[1][1]=Complex(-1.0/sqrt(2),0);
    C[2][2]=Complex(1.0/sqrt(2),0);
    C[2][3]=Complex(1.0/sqrt(2),0);
    C[3][2]=Complex(1.0/sqrt(2),0);
    C[3][3]=Complex(-1.0/sqrt(2),0);
    //*/
    // Operador de Hadamard construido mediante producto tensorial
    ///*
    //FILA 0 (relacionado con LEFT - DOWN)
    
    Hadamard[0][0]=Complex(1.0/2,0);
    Hadamard[0][1]=Complex(1.0/2,0);
    Hadamard[0][2]=Complex(1.0/2,0);
    Hadamard[0][3]=Complex(1.0/2,0);
    
    // FILA 1 (relacionado con LEFT - UP)
    
    Hadamard[1][0]=Complex(1.0/2,0);
    Hadamard[1][1]=Complex(-1.0/2,0);
    Hadamard[1][2]=Complex(1.0/2,0);
    Hadamard[1][3]=Complex(-1.0/2,0);
    
    // FILA 2 (relacionado con RIGHT - UP)
    
    Hadamard[2][0]=Complex(1.0/2,0);
    Hadamard[2][1]=Complex(1.0/2,0);
    Hadamard[2][2]=Complex(-1.0/2,0);
    Hadamard[2][3]=Complex(-1.0/2,0);
    
    // FILA 3 (relacionado con RIGHT - DOWN)
    
    Hadamard[3][0]=Complex(1.0/2,0);
    Hadamard[3][1]=Complex(-1.0/2,0);
    Hadamard[3][2]=Complex(-1.0/2,0);
    Hadamard[3][3]=Complex(1.0/2,0);
    //*/
        // Operador de Grover
    ///*
    //FILA 0 (relacionado con LEFT - DOWN)
    
    Grover[0][0]=Complex(-1.0/2,0);
    Grover[0][1]=Complex(1.0/2,0);
    Grover[0][2]=Complex(1.0/2,0);
    Grover[0][3]=Complex(1.0/2,0);
    
    // FILA 1 (relacionado con LEFT - UP)
    
    Grover[1][0]=Complex(1.0/2,0);
    Grover[1][1]=Complex(-1.0/2,0);
    Grover[1][2]=Complex(1.0/2,0);
    Grover[1][3]=Complex(1.0/2,0);
    
    // FILA 2 (relacionado con RIGHT - UP)
    
    Grover[2][0]=Complex(1.0/2,0);
    Grover[2][1]=Complex(1.0/2,0);
    Grover[2][2]=Complex(-1.0/2,0);
    Grover[2][3]=Complex(1.0/2,0);
    
    // FILA 3 (relacionado con RIGHT - DOWN)
    
    Grover[3][0]=Complex(1.0/2,0);
    Grover[3][1]=Complex(1.0/2,0);
    Grover[3][2]=Complex(1.0/2,0);
    Grover[3][3]=Complex(-1.0/2,0);
    //*/
            // Operador de DFT
    /*
    //FILA 0 (relacionado con LEFT - DOWN)
    
    C[0][0]=Complex(1.0/2,0);
    C[0][1]=Complex(1.0/2,0);
    C[0][2]=Complex(1.0/2,0);
    C[0][3]=Complex(1.0/2,0);
    
    // FILA 1 (relacionado con LEFT - UP)
    
    C[1][0]=Complex(1.0/2,0);
    C[1][1]=Complex(0,1.0/2);
    C[1][2]=Complex(-1.0/2,0);
    C[1][3]=Complex(0,-1.0/2);
    
    // FILA 2 (relacionado con RIGHT - UP)
    
    C[2][0]=Complex(1.0/2,0);
    C[2][1]=Complex(-1.0/2,0);
    C[2][2]=Complex(1.0/2,0);
    C[2][3]=Complex(-1.0/2,0);
    
    // FILA 3 (relacionado con RIGHT - DOWN)
    
    C[3][0]=Complex(1.0/2,0);
    C[3][1]=Complex(0,-1.0/2);
    C[3][2]=Complex(-1.0/2,0);
    C[3][3]=Complex(0,1.0/2);
    //*/
                                                                                    // Tamaño del sistema y medidores
    
    for(L=3;L<=30;L++){
        printf("Tamaño del sistema L=%i\n",L);
        promediotemporal=0;
        
        M1=2*L;
        M2=L;
        M3=0;
        M4=2*L*L+L;
        M5=4*L*L+2*L;
        M6=4*L*L+3*L;
        M7=4*L*L+4*L;
        M8=2*L*L+3*L;

                  
        for(caminante=1;caminante<=1000;caminante++){
                                                                  // Condiciones iniciales del caminante
        //I initialize the matrix QRW and AUX matrix as null
        
        for(k=0;k<((2*L+1)*(2*L+1));k++){
            for(j=0;j<4;j++){
                QRW[k][j]=Complex(0,0);
                AUX[k][j]=Complex(0,0);
            }
        }
        ///*
                    // Estado inicial simétrico 1
                
        QRW[L*(2*L+1)+L][0]=Complex(1.0/2,0);
        AUX[L*(2*L+1)+L][0]=Complex(1.0/2,0);
        
        QRW[L*(2*L+1)+L][1]=Complex(0,1.0/2);
        AUX[L*(2*L+1)+L][1]=Complex(0,1.0/2);
            
        QRW[L*(2*L+1)+L][2]=Complex(0,1.0/2);
        AUX[L*(2*L+1)+L][2]=Complex(0,1.0/2);

        QRW[L*(2*L+1)+L][3]=Complex(-1.0/2,0);
        AUX[L*(2*L+1)+L][3]=Complex(-1.0/2,0);
        //*/
                    // Estado inicial simétrico 2
        /*        
        QRW[L*(2*L+1)+L][0]=Complex(1.0/2,0);
        AUX[L*(2*L+1)+L][0]=Complex(1.0/2,0);
        
        QRW[L*(2*L+1)+L][1]=Complex(-1.0/2,0);
        AUX[L*(2*L+1)+L][1]=Complex(-1.0/2,0);
            
        QRW[L*(2*L+1)+L][2]=Complex(-1.0/2,0);
        AUX[L*(2*L+1)+L][2]=Complex(-1.0/2,0);

        QRW[L*(2*L+1)+L][3]=Complex(1.0/2,0);
        AUX[L*(2*L+1)+L][3]=Complex(1.0/2,0);
        //*/
                                                                                        // Aplico el operador U
        int t;
        t=1;
        while(t>0){
                                                                                        // Combinación de monedas
            if(t<0){
                for(k=0;k<4;k++){
                    for(j=0;j<4;j++){
                        C[k][j]=Grover[k][j];
                    }
                }
            }
            else{
                for(k=0;k<4;k++){
                    for(j=0;j<4;j++){
                        C[k][j]=Hadamard[k][j];
                    }
                }
            }
                                                                                            // Para sacar un gif
            // Escribo el mapa de calor
            /*
            i=0;
            for(k=0;k<=40400;k++){
                fprintf(fprobabilidad,"%lf\t",(Cabs(QRW[k][0])*Cabs(QRW[k][0])+Cabs(QRW[k][1])*Cabs(QRW[k][1])+Cabs(QRW[k][2])*Cabs(QRW[k][2])+Cabs(QRW[k][3])*Cabs(QRW[k][3])));
                i+=1;
                if(i%(2*L+1)==0){
                    fprintf(fprobabilidad,"\n");
                }
            }
            fprintf(fprobabilidad,"\n\n");
            //*/  
                                                                                        // Evolución del caminante
            k=0;
            while(k<((2*L+1)*(2*L+1))){
                                        /* 
                                            Primera fila: referida a spin LEFT - DOWN
                                            Segunda fila: referida a spin LEFT - UP
                                            Tercera fila: referida a spin RIGHT - DOWN
                                            Cuarta fila: referida a spin RIGHT - UP*/ 
                for(i=0;i<=3;i++){
                    if(Cabs(QRW[k][i])!=0){
                            // Left-down
                            
                            if(k==0){
                                AUX[4*L*L+4*L][0]=Cadd(AUX[4*L*L+4*L][0],Cmul(C[0][i],QRW[k][i]));    
                            }
                            else if(k%(2*L+1)==0){
                                AUX[k-1][0]=Cadd(AUX[k-1][0],Cmul(C[0][i],QRW[k][i])); 
                            }
                            else if(k<=2*L){
                                AUX[k+4*L*L+2*L-1][0]=Cadd(AUX[k+4*L*L+2*L-1][0],Cmul(C[0][i],QRW[k][i])); 
                            }
                            else{
                            AUX[k-2*L-2][0]=Cadd(AUX[k-2*L-2][0],Cmul(C[0][i],QRW[k][i]));
                            }
                            
                            
                            // Left-up
                            if(k==(4*L*L+2*L)){
                                AUX[2*L][1]=Cadd(AUX[2*L][1],Cmul(C[1][i],QRW[k][i]));
                            }
                            else if(k%(2*L+1)==0){
                            AUX[k+4*L+1][1]=Cadd(AUX[k+4*L+1][1],Cmul(C[1][i],QRW[k][i])); 
                            }
                            else if(k>(4*L*L+2*L)){
                                AUX[k-4*L*L-2*L-1][1]=Cadd(AUX[k-4*L*L-2*L-1][1],Cmul(C[1][i],QRW[k][i]));
                            }
                            else{
                                AUX[k+2*L][1]=Cadd(AUX[k+2*L][1],Cmul(C[1][i],QRW[k][i]));
                            }
                            
                            
                            //Right-down
                            if(k==2*L){
                                AUX[4*L*L+2*L][2]=Cadd(AUX[4*L*L+2*L][2],Cmul(C[2][i],QRW[k][i]));
                            }
                            else if(((k+1)%(2*L+1))==0){
                                AUX[k-4*L-1][2]=Cadd(AUX[k-4*L-1][2],Cmul(C[2][i],QRW[k][i]));
                            }
                            else if(k<2*L){
                                AUX[k+4*L*L+2*L+1][2]=Cadd(AUX[k+4*L*L+2*L+1][2],Cmul(C[2][i],QRW[k][i]));
                            }
                            else{
                                AUX[k-2*L][2]=Cadd(AUX[k-2*L][2],Cmul(C[2][i],QRW[k][i]));
                            }
                            
                            
                            // Right-up
                            if(k==(4*L*L+4*L)){
                                AUX[0][3]=Cadd(AUX[0][3],Cmul(C[3][i],QRW[k][i]));
                            }
                            else if(((k+1)%(2*L+1))==0){
                                AUX[k+1][3]=Cadd(AUX[k+1][3],Cmul(C[3][i],QRW[k][i]));
                            }
                            else if(k>=(4*L*L+2*L)){
                                AUX[k-4*L*L-2*L+1][3]=Cadd(AUX[k-4*L*L-2*L+1][3],Cmul(C[3][i],QRW[k][i]));
                            }
                            else{
                            AUX[k+2*L+2][3]=Cadd(AUX[k+2*L+2][3],Cmul(C[3][i],QRW[k][i]));
                            }
                            
                            
                        AUX[k][i]=Csub(AUX[k][i],QRW[k][i]);
                    }
                }
                    k=k+1;
                }

        for(k=0;k<((2*L+1)*(2*L+1));k++){
                for(j=0;j<4;j++){
                    QRW[k][j]=AUX[k][j];
                }
        }
        
        t=t+1;
                                                                        /* Esta parte del código se encarga de
                                                                            detectar la partícula*/
                                                                            
            // Mido arriba-derecha
            p=gsl_rng_uniform (r);
            probabilidad=Cabs(QRW[M1][0])*Cabs(QRW[M1][0])+Cabs(QRW[M1][1])*Cabs(QRW[M1][1])+Cabs(QRW[M1][2])*Cabs(QRW[M1][2])+Cabs(QRW[M1][3])*Cabs(QRW[M1][3]);
            if(p<probabilidad){
                promediotemporal=promediotemporal+t;
                break;
            }
            
            for(j=0;j<4;j++){
                    QRW[M1][j]=Complex(0,0);
                }
                
            NORMA=0.0;
            for(k=0;k<((2*L+1)*(2*L+1));k++){
                for(j=0;j<4;j++){
                    NORMA = NORMA + (Cabs(QRW[k][j])*Cabs(QRW[k][j]));
                }
            }
            
            RaizNorma=sqrt(NORMA);
            
            for(k=0;k<((2*L+1)*(2*L+1));k++){
                for(j=0;j<4;j++){
                    QRW[k][j]=RCmul(1.0/RaizNorma,QRW[k][j]);
                }
            }
                // Mido arriba-centro
            p=gsl_rng_uniform (r);
            probabilidad=Cabs(QRW[M2][0])*Cabs(QRW[M2][0])+Cabs(QRW[M2][1])*Cabs(QRW[M2][1])+Cabs(QRW[M2][2])*Cabs(QRW[M2][2])+Cabs(QRW[M2][3])*Cabs(QRW[M2][3]);
            if(p<probabilidad){
                promediotemporal=promediotemporal+t;
                break;
            }
            
            for(j=0;j<4;j++){
                    QRW[M2][j]=Complex(0,0);
                }
                
            NORMA=0.0;
            for(k=0;k<((2*L+1)*(2*L+1));k++){
                for(j=0;j<4;j++){
                    NORMA = NORMA + (Cabs(QRW[k][j])*Cabs(QRW[k][j]));
                }
            }
            
            RaizNorma=sqrt(NORMA);
            
            for(k=0;k<((2*L+1)*(2*L+1));k++){
                for(j=0;j<4;j++){
                    QRW[k][j]=RCmul(1.0/RaizNorma,QRW[k][j]);
                }
            }
                // Mido arriba-izquierda
            p=gsl_rng_uniform (r);
            probabilidad=Cabs(QRW[M3][0])*Cabs(QRW[M3][0])+Cabs(QRW[M3][1])*Cabs(QRW[M3][1])+Cabs(QRW[M3][2])*Cabs(QRW[M3][2])+Cabs(QRW[M3][3])*Cabs(QRW[M3][3]);
            if(p<probabilidad){
                promediotemporal=promediotemporal+t;
                break;
            }
            
            for(j=0;j<4;j++){
                    QRW[M3][j]=Complex(0,0);
                }
                
            NORMA=0.0;
            for(k=0;k<((2*L+1)*(2*L+1));k++){
                for(j=0;j<4;j++){
                    NORMA = NORMA + (Cabs(QRW[k][j])*Cabs(QRW[k][j]));
                }
            }
            
            RaizNorma=sqrt(NORMA);
            
            for(k=0;k<((2*L+1)*(2*L+1));k++){
                for(j=0;j<4;j++){
                    QRW[k][j]=RCmul(1.0/RaizNorma,QRW[k][j]);
                }
            }
                // Mido mitad-izquierda
            p=gsl_rng_uniform (r);
            probabilidad=Cabs(QRW[M4][0])*Cabs(QRW[M4][0])+Cabs(QRW[M4][1])*Cabs(QRW[M4][1])+Cabs(QRW[M4][2])*Cabs(QRW[M4][2])+Cabs(QRW[M4][3])*Cabs(QRW[M4][3]);
            if(p<probabilidad){
                promediotemporal=promediotemporal+t;
                break;
            }
            
            for(j=0;j<4;j++){
                    QRW[M4][j]=Complex(0,0);
                }
                
            NORMA=0.0;
            for(k=0;k<((2*L+1)*(2*L+1));k++){
                for(j=0;j<4;j++){
                    NORMA = NORMA + (Cabs(QRW[k][j])*Cabs(QRW[k][j]));
                }
            }
            
            RaizNorma=sqrt(NORMA);
            
            for(k=0;k<((2*L+1)*(2*L+1));k++){
                for(j=0;j<4;j++){
                    QRW[k][j]=RCmul(1.0/RaizNorma,QRW[k][j]);
                }
            }
            // Mido abajo-izquierda
            p=gsl_rng_uniform (r);
            probabilidad=Cabs(QRW[M5][0])*Cabs(QRW[M5][0])+Cabs(QRW[M5][1])*Cabs(QRW[M5][1])+Cabs(QRW[M5][2])*Cabs(QRW[M5][2])+Cabs(QRW[M5][3])*Cabs(QRW[M5][3]);
            if(p<probabilidad){
                promediotemporal=promediotemporal+t;
                break;
            }
            
            for(j=0;j<4;j++){
                    QRW[M5][j]=Complex(0,0);
                }
                
            NORMA=0.0;
            for(k=0;k<((2*L+1)*(2*L+1));k++){
                for(j=0;j<4;j++){
                    NORMA = NORMA + (Cabs(QRW[k][j])*Cabs(QRW[k][j]));
                }
            }
            
            RaizNorma=sqrt(NORMA);
            
            for(k=0;k<((2*L+1)*(2*L+1));k++){
                for(j=0;j<4;j++){
                    QRW[k][j]=RCmul(1.0/RaizNorma,QRW[k][j]);
                }
            }
                // Mido abajo-centro
            p=gsl_rng_uniform (r);
            probabilidad=Cabs(QRW[M6][0])*Cabs(QRW[M6][0])+Cabs(QRW[M6][1])*Cabs(QRW[M6][1])+Cabs(QRW[M6][2])*Cabs(QRW[M6][2])+Cabs(QRW[M6][3])*Cabs(QRW[M6][3]);
            if(p<probabilidad){
                promediotemporal=promediotemporal+t;
                break;
            }
            
            for(j=0;j<4;j++){
                    QRW[M6][j]=Complex(0,0);
                }
                
            NORMA=0.0;
            for(k=0;k<((2*L+1)*(2*L+1));k++){
                for(j=0;j<4;j++){
                    NORMA = NORMA + (Cabs(QRW[k][j])*Cabs(QRW[k][j]));
                }
            }
            
            RaizNorma=sqrt(NORMA);
            
            for(k=0;k<((2*L+1)*(2*L+1));k++){
                for(j=0;j<4;j++){
                    QRW[k][j]=RCmul(1.0/RaizNorma,QRW[k][j]);
                }
            }
                // Mido abajo-derecha
            p=gsl_rng_uniform (r);
            probabilidad=Cabs(QRW[M7][0])*Cabs(QRW[M7][0])+Cabs(QRW[M7][1])*Cabs(QRW[M7][1])+Cabs(QRW[M7][2])*Cabs(QRW[M7][2])+Cabs(QRW[M7][3])*Cabs(QRW[M7][3]);
            if(p<probabilidad){
                promediotemporal=promediotemporal+t;
                break;
            }
            
            for(j=0;j<4;j++){
                    QRW[M7][j]=Complex(0,0);
                }
                
            NORMA=0.0;
            for(k=0;k<((2*L+1)*(2*L+1));k++){
                for(j=0;j<4;j++){
                    NORMA = NORMA + (Cabs(QRW[k][j])*Cabs(QRW[k][j]));
                }
            }
            
            RaizNorma=sqrt(NORMA);
            
            for(k=0;k<((2*L+1)*(2*L+1));k++){
                for(j=0;j<4;j++){
                    QRW[k][j]=RCmul(1.0/RaizNorma,QRW[k][j]);
                }
            }
                // Mido mitad-centro
            p=gsl_rng_uniform (r);
            probabilidad=Cabs(QRW[M8][0])*Cabs(QRW[M8][0])+Cabs(QRW[M8][1])*Cabs(QRW[M8][1])+Cabs(QRW[M8][2])*Cabs(QRW[M8][2])+Cabs(QRW[M8][3])*Cabs(QRW[M8][3]);
            if(p<probabilidad){
                promediotemporal=promediotemporal+t;
                break;
            }
            
            for(j=0;j<4;j++){
                    QRW[M8][j]=Complex(0,0);
                }
                
            NORMA=0.0;
            for(k=0;k<((2*L+1)*(2*L+1));k++){
                for(j=0;j<4;j++){
                    NORMA = NORMA + (Cabs(QRW[k][j])*Cabs(QRW[k][j]));
                }
            }
            
            RaizNorma=sqrt(NORMA);
            
            for(k=0;k<((2*L+1)*(2*L+1));k++){
                for(j=0;j<4;j++){
                    QRW[k][j]=RCmul(1.0/RaizNorma,QRW[k][j]);
                }
            }
        }
            }
    
        fprintf(fhitting,"%i\t%lf\n",L,promediotemporal/1000.);
    }
                                                                                                //NORMA

    NORMA=0.0;
    for(k=0;k<((2*L+1)*(2*L+1));k++){
        for(j=0;j<4;j++){
            NORMA = NORMA + (Cabs(QRW[k][j])*Cabs(QRW[k][j]));
        }
    }
    
    RaizNorma=sqrt(NORMA);
    printf("Norma=%.25lf\n",RaizNorma);
        
    
    for(k=0;k<((2*L+1)*(2*L+1));k++){
        for(j=0;j<4;j++){
            QRW[k][j]=RCmul(1.0/RaizNorma,QRW[k][j]);
        }
    }
                                                                    // Imprimo datos en los ficheros para representarlo en 2D y en 3D
            // 2D
            i=0;
            for(k=0;k<((2*L+1)*(2*L+1));k++){
                fprintf(fmapa2d,"%lf\t",(Cabs(QRW[k][0])*Cabs(QRW[k][0])+Cabs(QRW[k][1])*Cabs(QRW[k][1])+Cabs(QRW[k][2])*Cabs(QRW[k][2])+Cabs(QRW[k][3])*Cabs(QRW[k][3])));
                i+=1;
                if(i%(2*L+1)==0){
                    fprintf(fmapa2d,"\n");
                }
            }
            // 3D
            x=0;
            y=0;
            for(k=0;k<((2*L+1)*(2*L+1));k++){
                    if((Cabs(QRW[k][0])*Cabs(QRW[k][0])+Cabs(QRW[k][1])*Cabs(QRW[k][1])+Cabs(QRW[k][2])*Cabs(QRW[k][2])+Cabs(QRW[k][3])*Cabs(QRW[k][3])!=0)){
                    fprintf(fmapa3d,"%i\t%i\t%lf\n",x,y,Cabs(QRW[k][0])*Cabs(QRW[k][0])+Cabs(QRW[k][1])*Cabs(QRW[k][1])+Cabs(QRW[k][2])*Cabs(QRW[k][2])+Cabs(QRW[k][3])*Cabs(QRW[k][3]));
                    }
                    x=x+1;
                    if(x%(2*L+1)==0){
                        y=y+1;
                        x=0;
                    }
                
            }
                                                                                // Cierro ficheros y libero memoria
    fclose(fmapa2d);
    fclose(fmapa3d);
    fclose(fhitting);
    gsl_rng_free (r);
    return 0;
}

