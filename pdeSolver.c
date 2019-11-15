/**
 * \file pdeSolver.c
 * \author Jackson Borguezani & Roberta Tomigian
 * \mainpage Solução Discreta para uma Equação Diferencial Parcial de Duas Variáveis
 * \section introSec ReadMe
 * Dada uma Equação Diferencial Parcial com duas variáveis independentes, discretizamos a malha utilizando Diferenças Finitas Centrais e o método de Gauss-Seidel para encontrar a possível solução do Sistema Linear. 
 *
 * \brief **Introdução à Computação Científica CI1164** \n Trabalho 1 -- Prof. Armando Nicolui \n **Alunos :** \n Jackson Rossi Borguezani GRR20176573 \n Roberta Tomigian GRR20171631
 * 
 * \date 22 out 2019
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <sys/time.h>
#include <math.h>
#include "utils.h"
#include "SistemasLineares.h"
#include "pdeSolver.h"
#include <getopt.h>


#define M_PI 3.14159265358979323846
#define fronteira_d 0
#define fronteira_e 0

/*!
	\fn void gera_matriz(SistLinear_t *SL)
	*
	\brief Monta a matriz pentadiagonal do SL.
	*
	\param SL: Ponteiro do Sistema Linear
	*
	\details Essa função gera a matriz de coeficientes a partir da aproximação de valores da Equação Diferencial Parcial
	*
*/
void gera_matriz(SistLinear_t *SL)
{
    real_t hx = M_PI / (SL->nx+1);
    real_t hy = M_PI / (SL->ny+1);

	// Diagonal Principal
	SL->dp = ((4*M_PI*M_PI)+(2/(hx*hx))+(2/(hy*hy)));

	// Diagonal Superior 
	SL->ds = -((1/(hy*hy)) - (1/(2*hy)));

	// Diagonal Inferior 
	SL->di = -((1/(hx*hx)) + (1/(2*hx)));

	// Diagonal Superior Afastada 
	SL->dsa = -((1/(hx*hx)) - (1/(2*hx)));

	// Diagonal Inferior Afastada 
	SL->dia = -((1/(hy*hy)) + (1/(2*hy)));
	return;
}

/*!
	\fn void gera_vetor_b(SistLinear_t *SL)
	*
	\brief Calcula o vetor dos termos independentes B do SL.
	*
	\param SL: Ponteiro do Sistema Linear
	*
	\details Essa função gera o vetor de termos independentes b a partir da aproximação de valores da Equação Diferencial Parcial
	*
*/
void gera_vetor_b(SistLinear_t *SL)
{
    real_t hx = M_PI / (SL->nx+1);
    real_t hy = M_PI / (SL->ny+1);  
	real_t x = 0.0;
	real_t y = 0.0;

    for (int i = 0; i < (SL->nx+2); ++i, x += hx)
    {
		y = 0.0;
        for (int j = 0; j < (SL->ny+2); ++j, y += hy)
        {
            SL->b[(i*(SL->ny+2)) + j] = 4*M_PI*M_PI*(sin(2*M_PI*x)*sinh(M_PI*y)+sin(2*M_PI*(M_PI-x))*sinh(M_PI*(M_PI-y)));
        }
    }
	return;
}

/*!
    \fn void contorno_x(SistLinear_t *SL)
    *
    \brief Monta o vetor X com as constantes correspondentes da discretização
    *
    \param SL: Sistema Linear
    *
*/
void contorno_x(SistLinear_t *SL){
    real_t hx = M_PI / (SL->nx+1);
    real_t hy = M_PI / (SL->ny+1);	
	real_t x = 0;
	real_t y = 0;
	for(int i=0;i < ((SL->nx)+2);++i, x += hx){
		y = 0;
		for(int j=0;j < ((SL->ny)+2);++j, y += hy){
			if(i == 0 || i == ((SL->nx)+1)){
				SL->x[(i*((SL->ny)+2))+j] = 0;
			}else if(j == 0){
				SL->x[(i*((SL->ny)+2))+j] = sin(2*M_PI*(M_PI-x))*sinh(M_PI*M_PI);
			}else if(j == ((SL->ny)+1)){
				SL->x[(i*((SL->ny)+2))+j] = sin(2*M_PI*x)*sinh(M_PI*M_PI);
			}else{
				SL->x[(i*((SL->ny)+2))+j] = 0;
			}
		}
	}
	return;
}


/*!
	\fn void saida_gnuplot(Metrica *P, int flagArq, char *arqOut, SistLinear_t *SL)
	*
	\brief Gera a saída da solução no formato aceito pelo gnuplot
	*
	\param *P: Estrutura que armazena o tempo, número de iterações e o resíduo de cada iteração de Gauss
	*
	\param flagArq: Flag sinalizando se há um arquivo de saída
    *
    \param *arqOut: Arquivo de saída (pode ser nulo)
	*
    \param *SL: Sistema Linear
    *
*/
void saida_gnuplot(Metrica *P, int flagArq, char *arqOut, SistLinear_t *SL)
{
    real_t hx = M_PI / (SL->nx+1);
    real_t hy = M_PI / (SL->ny+1);
	real_t x = 0.0;
	real_t y = 0.0;

    FILE* fp;
    if(flagArq)
    {
        fp = fopen(arqOut ,"w+");
        if (fp == NULL)
            fprintf(stderr, "Um erro ocorreu ao tentar criar o arquivo.\n");
    }else
        fp = stdout;
        
    fprintf(fp,"###########\n");
    fprintf(fp,"# Tempo Método GS: %f\n", P->mediaTempo);
    fprintf(fp,"#\n");
    fprintf(fp,"# Norma L2 do Residuo\n");

    for (int i=1; i <= P->iter; ++i)
        fprintf(fp, "# i = %d: %f\n",i,P->norma[i]);

    fprintf(fp,"#\n");
    fprintf(fp,"###########\n");

    for (int i=0;i<((SL->nx)+2);i++, x += hx){
		y = 0.0;
		for(int j=0; j<((SL->ny)+2);j++, y += hy){
			fprintf(fp, "%lf %lf %lf\n",x,y,SL->x[i*((SL->ny)+2)+j]);
		}
		fprintf(fp,"\n");
    }
    fclose(fp);   
	return;
}

int main(int argc, char *argv[])
{
    /*=========================== Le os parametros ===========================*/
    int opt;
    char *arqOut = NULL;
    int flagArq = 0, x, y, maxIter;

    const struct option stopcoes[] = {
        {"nx", required_argument, 0, 'x'},
        {"ny", required_argument, 0, 'y'},
        {"maxIt", required_argument, 0, 'i'},
        {"arqOut", optional_argument, 0, 'o'},
        {0, 0, 0, 0},
    };

    while ((opt = getopt_long(argc, argv, "nx:ny:i:o:", stopcoes, NULL)) > 0)
    {
        switch (opt)
        {
        case 'n':
            break;
        case 'x': /* opção -nx */
            x = atoi(optarg);
            break;
        case 'y': /* opção -ny */
            y = atoi(optarg);
            break;
        case 'i': /* opção -i */
            maxIter = atoi(optarg);
            break;
        case 'o': /* opção -o */
            arqOut = optarg;
            flagArq = 1;
            break;
        default:
            fprintf(stderr, "Opção inválida ou faltando argumento\n");
            return -1;
        }
    }
    /*=====================================================================*/


    /*=========================== Aloca e Monta SL ========================*/
    SistLinear_t *SL;

    SL = alocaSistLinear(x, y);
    inicializaSistLinear(SL, x, y);
    gera_matriz(SL);
    gera_vetor_b(SL);
    /*=====================================================================*/

    int i, tamLin = SL->nx * SL->ny;

    /*=========================== Resolve o SL ========================*/
    Metrica *P;
    P = alocaMetrica(x,y,maxIter);
    int ret = gaussSeidel(SL,maxIter,P);

    /*=========================== Monta saida do programa ========================*/
    saida_gnuplot(P,flagArq,arqOut,SL);

    /*===================================================================================*/
    liberaSistLinear(SL);
    return 0;
}

