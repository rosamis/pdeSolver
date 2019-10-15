/**
 * \file pdeSolver.c
 * \author Jackson Borguezani & Roberta Tomigian
 * \mainpage Solução Discreta para uma Equação Diferencial Parcial
 * \section introSec ReadMe
 * Dada uma Equação Diferencial Parcial com duas variáveis independentes, discretizamos a malha utilizando Diferenças Finitas Centrais e o método de Gauss-Seidel para encontrar a possível solução do Sistema Linear. 
 *
 * \brief **Introdução à Computação Científica CI1164** \n Trabalho 1 -- Prof. Armando Nicolui \n **Alunos :** \n Jackson Rossi Borguezani GRR20176573 \n Roberta Tomigian GRR20171631
 * 
 * \date 03 out 2019
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
	\fn real_t quadrado(real_t p)
	*
	\brief Função que calcula o quadrado de um número
	*
	\param p: o valor a ser calculado
	*
*/

real_t quadrado(real_t p)
{
    return p*p;
}

/*!
	\fn real_t gera_d_principal(real_t hx_quadrado, real_t hy_quadrado, real_t pi_quadrado)
	*
	\brief Função que retorna a fórmula usada no calculo dos coeficientes da diagonal principal
	*
	\param hx_quadrado: Deslocamento no eixo x ao quadrado 
	*
	\param hy_quadrado: Deslocamento no eixo y ao quadrado 
    *
	\param pi_quadrado: o valor de pi ao quadrado
*/

real_t gera_d_principal(real_t hx_quadrado, real_t hy_quadrado, real_t pi_quadrado)
{

    return 4 * hy_quadrado + 4 * hx_quadrado + 8 * pi_quadrado * hx_quadrado* hy_quadrado;
}

/*!
	\fn real_t gera_d_superior(real_t hx, real_t hy_quadrado)
	*
	\brief Função que retorna a fórmula usada no calculo dos coeficientes da diagonal superior
	*
	\param hx: Deslocamento no eixo x 
	*
	\param hy_quadrado: Deslocamento no eixo y ao quadrado 
*/

real_t gera_d_superior(real_t hx, real_t hy_quadrado)
{
    return hx * hy_quadrado - 2 * hy_quadrado;
}

/*!
	\fn real_t gera_d_inferior(real_t hx, real_t hy_quadrado)
	*
	\brief Função que retorna a fórmula usada no calculo dos coeficientes da diagonal inferior
	*
	\param hx: Deslocamento no eixo x 
	*
	\param hy_quadrado: Deslocamento no eixo y ao quadrado 
*/

real_t gera_d_inferior(real_t hx, real_t hy_quadrado)
{
    return -(2 * hy_quadrado) - hx * hy_quadrado;
}

/*!
	\fn real_t gera_ds_afastada(real_t hx_quadrado, real_t hy)
	*
	\brief Função que retorna a fórmula usada no calculo dos coeficientes da diagonal superior afastada
	*
	\param hy: Deslocamento no eixo y  
	*
	\param hx_quadrado: Deslocamento no eixo x ao quadrado 
*/

real_t gera_ds_afastada(real_t hx_quadrado, real_t hy)
{
    return hx_quadrado * hy - 2 * hx_quadrado;
}

/*!
	\fn real_t gera_ds_afastada(real_t hx_quadrado, real_t hy)
	*
	\brief Função que retorna a fórmula usada no calculo dos coeficientes da diagonal inferior afastada
	*
	\param hy: Deslocamento no eixo y  
	*
	\param hx_quadrado: Deslocamento no eixo x ao quadrado 
*/

real_t gera_di_afastada(real_t hx_quadrado, real_t hy)
{
    return -2 * hx_quadrado - hx_quadrado * hy;
}

/*!
	\fn real_t funcao_f(int i, int j, real_t pi_quadrado)
	*
	\brief Função que retorna a fórmula usada no calculo da função principal
	*
	\param i: o valor de i 
	*
	\param j: o valor de j 
	*
	\param pi_quadrado: o valor de pi_quadrado ao quadrado 
*/

real_t funcao_f(int i, int j, real_t pi_quadrado)
{
    return (4 * pi_quadrado * (sin(2 * M_PI * i) * sinh(M_PI * j) + sin(2 * M_PI * (M_PI - i)) * sinh(M_PI * (M_PI - j))));
}

/*!
	\fn real_t b_principal(int i, int j, real_t hx_quadrado, real_t hy_quadrado, real_t pi_quadrado)
	*
	\brief Função que retorna o calculo da função das constantes nos termos independentes de b
	*
	\param i: o valor de i 
	*
	\param j: o valor de j 
	*
	\param pi_quadrado: o valor de pi_quadrado ao quadrado 
	*
	\param hx_quadrado: Deslocamento no eixo x ao quadrado  
	*
	\param hy_quadrado: Deslocamento no eixo y ao quadrado 
*/

real_t b_principal(int i, int j, real_t hx_quadrado, real_t hy_quadrado, real_t pi_quadrado)
{
    real_t f = funcao_f(i,j,pi_quadrado);

    return 2 * hx_quadrado * hy_quadrado * f; //-------------------------------------------

}

/*!
	\fn real_t b_limite_superior(real_t hx_quadrado, real_t hy, int i, real_t pi_quadrado)
	*
	\brief Função que retorna o calculo da função das constantes dos termos independentes de b no limite superior
	*
	\param i: o valor de i 
	*
	\param pi_quadrado: o valor de pi_quadrado ao quadrado 
	*
	\param hx_quadrado: Deslocamento no eixo x ao quadrado 
	*
	\param hy: Deslocamento no eixo y   
*/

real_t b_limite_superior(real_t hx_quadrado, real_t hy, int i, real_t pi_quadrado)
{
    return (hx_quadrado * hy - 2 * hx_quadrado) * sin(2 * M_PI * i) * sinh(pi_quadrado);
}

/*!
	\fn real_t b_limite_superior(real_t hx_quadrado, real_t hy, int i, real_t pi_quadrado)
	*
	\brief Função que retorna o calculo da função das constantes dos termos independentes de b no limite inferior
	*
	\param i: o valor de i 
	*
	\param pi_quadrado: o valor de pi_quadrado ao quadrado 
	*
	\param hx_quadrado: Deslocamento no eixo x ao quadrado  
	*
	\param hy: Deslocamento no eixo y  
*/

real_t b_limite_inferior(real_t hx_quadrado, real_t hy, int i, real_t pi_quadrado)
{
    return (-2 * hx_quadrado - hx_quadrado * hy) * sin(2 * M_PI * (M_PI - i)) * sinh(pi_quadrado);
}

/*!
	\fn real_t b_limite_esquerda(real_t hx, real_t hy_quadrado)
	*
	\brief Função que retorna o calculo da função das constantes dos termos independentes de b no limite à esquerda
	*
	\param hy_quadrado: Deslocamento no eixo y ao quadrado 
	*
	\param hx: Deslocamento no eixo x
*/

real_t b_limite_esquerda(real_t hx, real_t hy_quadrado)
{
    return (-2 * hy_quadrado - hx * hy_quadrado) * fronteira_e;
}

/*!
	\fn real_t b_limite_direita(real_t hx, real_t hy_quadrado)
	*
	\brief Função que retorna o calculo da função das constantes dos termos independentes de b no limite à direita
	*
	\param hy_quadrado: Deslocamento no eixo y ao quadrado 
	*
	\param hx: Deslocamento no eixo x
*/

real_t b_limite_direita(real_t hx, real_t hy_quadrado)
{
    return (hx * hy_quadrado - 2 * hy_quadrado) * fronteira_d;
}

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
    int k;
    real_t hx = M_PI / SL->nx;
    real_t hy = M_PI / SL->ny;
    real_t hx_quadrado = quadrado(hx);
    real_t hy_quadrado = quadrado(hy);
    real_t pi_quadrado = quadrado(M_PI);

    int tamLin = SL->nx * SL->ny;

    for (k = 1; k <= tamLin; k++) // Diagonal Principal
        SL->dp[k] = gera_d_principal(hx_quadrado,hy_quadrado,pi_quadrado);

    for (k = 1; k <= tamLin - 1; k++) // Diagonal Superior
        SL->ds[k] = gera_d_superior(hx,hy_quadrado);

    for (k = 2; k <= tamLin; k++) // Diagonal Inferior
        SL->di[k] = gera_d_inferior(hx,hy_quadrado);

    for (k = 1; k <= tamLin - SL->nx; k++) // Diagonal Superior Afastada
        SL->dsa[k] = gera_ds_afastada(hx_quadrado,hy);

    for (k = SL->ny + 1; k <= tamLin; k++) // Diagonal Inferior Afastada
        SL->dia[k] = gera_di_afastada(hx_quadrado,hy);
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
    int k = 1;
    int tamLin = SL->nx * SL->ny;
    real_t hx = M_PI / SL->nx;
    real_t hy = M_PI / SL->ny;
    real_t hx_quadrado = quadrado(hx);
    real_t hy_quadrado = quadrado(hy);   
    real_t pi_quadrado = quadrado(M_PI);

    for (int j = 1; j <= SL->ny; j++)
    {
        for (int i = 1; i <= SL->nx; i++)
        {
            SL->b[k] = b_principal(i,j,hx_quadrado,hy_quadrado,pi_quadrado);
            if (i == 1)
            {
                SL->b[k] = SL->b[k] - b_limite_esquerda(hx,hy_quadrado);
                if (j > 1)
                    SL->di[k - 1] = 0.0;
            }
            if (i == SL->nx - 1)
            {
                SL->b[k] = SL->b[k] - b_limite_direita(hx,hy_quadrado);
                if (j < SL->ny)
                    SL->ds[k + 1] = 0.0;
            }
            if (j == 1)
                SL->b[k] = SL->b[k] - b_limite_inferior(hx_quadrado,hy,i,pi_quadrado);

            if (j == SL->ny)
                SL->b[k] = SL->b[k] - b_limite_superior(hx_quadrado,hy,i,pi_quadrado);
            k++;
        }
    }
}

/*!
	\fn void solucao(real_t hx, real_t hy, real_t *x, real_t *y, int ny, int nx)
	*
	\brief Essa função calcula os valores dos eixos x e y da malha
	*
	\param nx: Número de pontos a serem calculados no eixo x
	*
	\param ny: Número de pontos a serem calculados no eixo y
	*
	\param hx: Deslocamento no eixo x
	*
	\param hy: Deslocamento no eixo y
	*
	\param y: Vetor de valores em y
	*
	\param x: Vetor de valores em x
	*
*/
void solucao(real_t hx, real_t hy, real_t *x, real_t *y, int ny, int nx)
{
    for(int j = 0; j<ny;++j)
    { 
        for(int i = 0; i<nx; ++i)
        {
            y[j] = hy*j;
            x[i] = hx*i;
        }
    }
}

/*!
	\fn real_t* aloca_vetor(int nx, int ny)
	*
	\brief Aloca e inicializa vetor 
	*
	\param nx: Número de pontos a serem calculados no eixo x
	*
	\param ny: Número de pontos a serem calculados no eixo y
	*
	\details Essa função gera vetor de dimensão nx*ny
	*
*/
real_t* aloca_vetor(int nx, int ny, int flag)
{
	int t;
	if (flag)
		t = nx;
	else 
		t = ny;
		
	real_t *v = (real_t*)malloc(t * sizeof(real_t));
	
	if(v)
		for(int i = 0; i < t; ++i)
			v[i] = 0.0;
	
	return v;
}

/*!
	\fn void saida_gnuplot(real_t *R, real_t n, int flagArq, char *arqOut, double tempo, int iter)
	*
	\brief Gera a saída da solução no formato aceito pelo gnuplot
	*
	\param nx: Número de pontos a serem calculados no eixo x
	*
	\param ny: Número de pontos a serem calculados no eixo y
	*
	\details Essa função gera vetor de dimensão nx*ny
	*
*/
void saida_gnuplot(Metrica *P, int flagArq, char *arqOut, SistLinear_t *SL)
{
    real_t hx = M_PI / SL->nx;
    real_t hy = M_PI / SL->ny;
    real_t *vx = aloca_vetor(SL->nx,SL->ny,1);
    real_t *vy = aloca_vetor(SL->nx,SL->ny,0);
    solucao(hx,hy,vx,vy,SL->ny,SL->nx);

    FILE* fp;
    if(flagArq)
    {
        fp = fopen(arqOut ,"w+");
        if (fp == NULL)
            printf("Um erro ocorreu ao tentar criar o arquivo.\n");
    }else
        fp = stdout;
        
    fprintf(fp,"###########\n");
    fprintf(fp,"# Tempo Método GS: %f\n", P->mediaTempo);
    fprintf(fp,"#\n");
    fprintf(fp,"# Norma L2 do Residuo\n");

    for (int i=1; i<=P->iter; ++i)
        fprintf(fp, "# i = %d: %f\n",i,P->norma[i]);

    fprintf(fp,"#\n");
    fprintf(fp,"###########\n");

    for (int i=0;i<SL->nx;i++){
		for(int j=0; j<SL->ny;j++){
			fprintf(fp, "%f %f %f\n",vx[j],vy[i],SL->x[j+i+1]);
	}
        
    }

    fclose(fp);   
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

