/* Daianna Gonzalez Padilla.
   Computo cientifico, 2do Semestre 2021, LCG: Proyecto Final
   Programa 4: Aproximacion funcional por el metodo de minimos cuadrados 
   
   Diccionario de variables:
    empezar permite reiniciar el programa segun su valor
    puntos es una matriz para almacenar las coordenadas de los puntos
    matriz es una matriz para resolver el sistema de ecuaciones 
    px es un vector para almacenar los P(x) con P el polinomio obtenido y x los valores en x de los puntos dados
    n y t es el numero de puntos
    m es el grado del polinomio al cual se quieren ajustar los puntos
    tam_filas es el numero de filas que tiene matriz
    tam_columnas es el numero de columnas que tiene matriz
    i,j,k son contadores y representan indices 
    s,s1,s2,s3,c, ppx y py son acumuladores
    coef es el coeficiente de correlacion
    aii es una entrada de la matriz en la diagonal en la coordenada i,i
    aki es un entrada de la matriz en la coordenada k,i         */



#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

//Se define el tipo de dato byte
typedef unsigned char byte;

int main(int argc, char *argv[])
{
    
    // Se declaran las variables 
    	float **matriz = NULL; //Matriz dada por un puntero a punteros a float
    	float **puntos = NULL; 
    	float *px = NULL; 
      	byte i,j,k,n,m,tam_filas, tam_columnas, t, empezar; 
      	float s,s1,s2,s3,py,c,aii,aki, ppx, coef;
      	
    //Se pregunta al usuario si quiere comenzar
    puts ("Para comenzar oprima 1\n");
    scanf("%d", &empezar);
    while (empezar==1) {
 
         
        
    	puts("Para cuantos puntos desea ajustar la curva?\n "); //Se pide el numero de puntos n
    	scanf("%d", &n);
    	fflush(stdin);
    	
    
       // Se solicita la memoria para la matriz puntos de tamaño 2n para almacenar las coordenadas:
        
      	// Se reserva primero el espacio para un  vector de tamaño n de punteros a float, se evalua que la peticion fue exitosa	
      	puntos = (float**) malloc (n * sizeof(float *));
      	if(!puntos) {
    		puts("Error: Memoria insuficiente");
            exit(1);  }
      	//Se reserva el espacio al que apunta cada puntero del vector: dos floats para cada punto (x,y) dado
      	for (i=0 ; i < n ; i++) {
        	*(puntos + i) = (float *) malloc (2 * sizeof(float));
        	//Si la peticion no fue exitosa se libera el bloque del vector puntos y los valores previos guardados
        	if(!*(puntos + i)) {
        		puts("Error: Memoria insuficiente");
        		for(j= 0; j< i-1; j++)
    		          free(puntos + i);
                free(puntos);
        		exit(1);  
          }
      	}
       //Se piden los n puntos y se guardan en la matriz puntos
    	for (i=0 ; i < n ; i++) {
    	    printf("Ingrese la coordenada en x del punto %d",i+1);
    	    printf("\n");
            scanf("%f", &puntos[i][0]);
            fflush(stdin);
    	    printf("Ingrese la coordenada en y del punto %d",i+1);
    	    printf("\n");
    	    scanf("%f", &puntos[i][1]);
    	    fflush(stdin);
    
          }
        
      	
        t=n;
        //Se pide el grado de P(x) para ajustar una curva, se guarda en m
    	puts("Ingrese el grado del polinomio para ajustar una curva:\n ");
    	scanf("%d", &m);
    	fflush(stdin);
    	//Las dimensiones de la matriz de los coeficientes de P(x) son m+1 filas y  x m+2 columnas
    	tam_filas=m+1;
    	tam_columnas=m+2;
    	
    	
    	// Se dimensiona la matriz de los coeficientes aij para resolver el polinomio:
      	// Se reserva espacio para un vector ( de tamaño = numero de filas) de punteros a float, cada uno de estos apuntando a otro vector de floats
        //  de tamaño = numero de columnas 	
      	matriz = (float **) malloc (tam_filas * sizeof(float *));
      	//Si no hay memoria suficiente se liberan los bloques de las coordenadas y su matriz
      	if(!matriz) {
    		puts("Error: Memoria insuficiente");
    		for(j= 0; j<t; j++)
                   free(puntos + j);
            free(puntos);
            exit(1);  }
      	
      	//Ya que se tiene el vector con el numero de celdas como el de filas, cada uno contiene un puntero a un vector de tantas celdas como numero de columnas
        for (i=0 ; i < tam_filas ; i++) {
        	*(matriz+i) = (float *) malloc (tam_columnas * sizeof(float));
        	// Si la peticion de memoria no fue exitosa se liberan los bloques ya reservados
        	if(!matriz[i]) {
        		puts("Error: Memoria insuficiente");
        		for(j= 0; j< i-1; j++)
    		          free(matriz+j);
                for(k= 0; k<t; k++)
                   free((puntos + k));
                free(puntos);
                free(matriz);
        		exit(1);  
          }
      	}
      	
       //Se llena la matriz de coeficientes donde cada entrada ij es la suma de los X**(i+j) y las entradas de la ultima columna son la suma de Y*(X**i),
       // con i las filas y j las columnas y (X,Y) los puntos dados
    	for (i=0 ; i < tam_filas ; i++) {
          for (j=0 ; j < tam_columnas ; j++) {
              s=0;
              c=0;
              for (k=0; k<t; k++){
                s=s+(pow(puntos[k][0],(i+j))); 
                c=c+((pow(puntos[k][0],i)*puntos[k][1]));   }
    
              if (j!=tam_columnas-1) {
               matriz[i][j]=s; 
               }    
              else {
               matriz[i][j]=c;
               }
            }
         }
    
       
       //Obtención de la matriz escalonada por Gauss-Jordan a partir de la matriz de coeficientes 
       // i es la fila o renglon Ri, j la columna Cj
       i=0;
       while (i<tam_filas) {
          aii=matriz[i][i];
          for (j=0 ; j < tam_columnas ; j++) {
            //En cada iteracion se dividen las entradas del renglon Ri entre la entrada de la diagonal aii para obtener 1 en esa posicion 
            matriz[i][j]= matriz[i][j]/aii; 
            }
          //Los demas renglones Rk con k!=i resultan de restarlos con los elementos del renglon Ri*aki para obtener ceros en las entradas 
          // restantes de la columna Cj
          for (k=0 ; k < tam_filas ; k++) {
              aki=matriz[k][i];
              for (j=0 ; j < tam_columnas ; j++) {
                if (k!=i) {
                   matriz[k][j]= matriz[k][j]-(matriz[i][j]*aki);  }
              }
          }
          i++;
       }
       
       //Los valores de los coeficientes son los valores de la ultima columna de la matriz, se imprimen
       puts("El polinomio que se ajusta a sus puntos es:\n");
       for (k=0 ; k < tam_filas ; k++) {
           if (k==tam_filas-1){
      	    printf("%.4f", matriz[k][tam_columnas-1]);  printf("%s", " x^"); printf("%d\n", k); }
           else {
            printf("%.4f", matriz[k][tam_columnas-1]);  printf("%s", " x^"); printf("%d", k); printf("%s", "    +   ");
            }
           } 
       
       //Calcular el coef de Pearson:
       
       // Se reserva el espacio para un  vector 	
      	px = (float*) malloc (t * sizeof(float));
       //Si la peticion no fue exitosa se libera el bloque del vector puntos y los valores previos guardados
      	if(!px) {
    		puts("Error: Memoria insuficiente");
    		for(j= 0; j< tam_filas; j++)
              free(matriz+j);
            for(k= 0; k<t; k++)
              free((puntos + k));
            free(puntos);
            free(matriz);
            exit(1);  
       }
       
       //Promedio de los elementos del vector con los Y's de los puntos dados
       py=0;
       for (i=0; i<t; i++) {
           py=py+puntos[i][1]; } 
       py=py/t;
      
       
       //Elementos del vector con los P(x): cada x se evalua en el polinomio obtenido, se calcula el promedio de tales elementos
       ppx=0;
       for (i=0; i<t; i++) {
           s=0;
           for (j=0; j<tam_filas; j++) {
               s=s+ (pow(puntos[i][0],j)* matriz[j][tam_columnas-1]); } 
           px[i]=s;
           ppx=ppx+px[i]; } 
       ppx=ppx/t;
      
          
      
       //Se calcula el numerador s1: suma de (y-py)(px-ppx) para y las coordenadas en y de los puntos dados, y px los y's obtenidos al evaluar las x en P(x)
       // y con py y ppx el  promedio de los y's en cada vector
       // Con s2 y s3 se calculan las sumatorias del denominador: los y-py   y los  p(x)-ppx
       s1=0;
       s2=0;
       s3=0;
       for (i=0; i<t; i++) {
           s1=s1+ ((puntos[i][1]-py)*(px[i]-ppx));
           s2=s2+(pow(puntos[i][1]-py,2));
           s3=s3+(pow(px[i]-ppx,2)); }
          
       //Expresion final: se imprime el coeficiente de determinacion 
       coef=s1/(sqrt(s2)* sqrt(s3));
       printf("R^2: %.4f\n ", pow(coef,2)); puts("\n");
       
       
       //Al salir se liberan todos los bloques reservados
       for(j= 0; j< tam_filas; j++)
         free(matriz[j]);
       for(k= 0; k<t; k++)
        free((puntos[k]));
       free(puntos);
       free(matriz);
       free(px);  
       
       //Se vuelve a preguntar al usuario si desea ajustar otro conjunto de puntos
       puts("Para ajustar otro conjunto de puntos oprima 1, para salir 2\n");
       scanf ("%d", &empezar);
    
       
        
    }
   puts ("Gracias por usar este programa \n");
   return 0;
         
    

} 

