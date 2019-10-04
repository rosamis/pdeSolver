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
//#define EPS 1.0e-4

real_t quadrado(real_t p)
{
    return p*p;
}
real_t gera_d_principal(real_t hx, real_t hy)
{
    real_t hx_quadrado = quadrado(hx);
    real_t hy_quadrado = quadrado(hy);
    real_t pi_quadrado = quadrado(M_PI);
    return 4 * hy_quadrado + 4 * hx_quadrado + 8 * pi_quadrado * hx_quadrado* hy_quadrado;
}

real_t gera_d_superior(real_t hx, real_t hy)
{
    real_t hy_quadrado = quadrado(hy);
    return hx * hy_quadrado - 2 * hy_quadrado;
}

real_t gera_d_inferior(real_t hx, real_t hy)
{
    real_t hy_quadrado = quadrado(hy);
    return -(2 * hy_quadrado) - hx * hy_quadrado;
}

real_t gera_ds_afastada(real_t hx, real_t hy)
{
    real_t hx_quadrado = quadrado(hx);
    return hx_quadrado * hy - 2 * hx_quadrado;
}

real_t gera_di_afastada(real_t hx, real_t hy)
{
    real_t hx_quadrado = quadrado(hx);
    return -2 * hx_quadrado - hx_quadrado * hy;
}

real_t funcao_f(int i, int j, real_t pi_quadrado)
{
    return (4 * pi_quadrado * (sin(2 * M_PI * i) * sinh(M_PI * j) + sin(2 * M_PI * (M_PI - i)) * sinh(M_PI * (M_PI - j))));
}

real_t b_principal(int i, int j, real_t hx_quadrado, real_t hy_quadrado, real_t pi_quadrado)
{
    real_t f = funcao_f(i,j,pi_quadrado);

    return 2 * hx_quadrado * hy_quadrado * f;

}

real_t b_limite_superior(real_t hx_quadrado, real_t hy, int i, real_t pi_quadrado)
{
    return (hx_quadrado * hy - 2 * hx_quadrado) * sin(2 * M_PI * i) * sinh(pi_quadrado);
}

real_t b_limite_inferior(real_t hx_quadrado, real_t hy, int i, real_t pi_quadrado)
{
    return (-2 * hx_quadrado - hx_quadrado * hy) * sin(2 * M_PI * (M_PI - i)) * sinh(pi_quadrado);
}

real_t b_limite_esquerda(real_t hx, real_t hy_quadrado)
{
    return (-2 * hy_quadrado - hx * hy_quadrado) * fronteira_e;
}

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
	\details A função libera todos os ponteiros da estrutura do Sistema Linear.
	*
*/
void gera_matriz(SistLinear_t *SL)
{
    int k;
    real_t hx = M_PI / SL->nx;
    real_t hy = M_PI / SL->ny;

    int tamLin = SL->nx * SL->ny;

    for (k = 1; k <= tamLin; k++) // Diagonal Principal
        SL->dp[k] = gera_d_principal(hx,hy);

    for (k = 1; k <= tamLin - 1; k++) // Diagonal Superior
        SL->ds[k] = gera_d_superior(hx,hy);

    for (k = 2; k <= tamLin; k++) // Diagonal Inferior
        SL->di[k] = gera_d_inferior(hx,hy);

    for (k = 1; k <= tamLin - SL->nx; k++) // Diagonal Superior Afastada
        SL->dsa[k] = gera_ds_afastada(hx,hy);

    for (k = SL->ny + 1; k <= tamLin; k++) // Diagonal Inferior Afastada
        SL->dia[k] = gera_di_afastada(hx,hy);
}

/*!
	\fn void gera_vetor_b(SistLinear_t *SL)
	*
	\brief Calcula o vetor dos termos independentes B do SL.
	*
	\param SL: Ponteiro do Sistema Linear
	*
	\details A função libera todos os ponteiros da estrutura do Sistema Linear.
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

void saida_gnuplot(real_t *R, real_t n, int flagArq, char *arqOut, double tempo, int iter)
{
    if(flagArq)
    {
        FILE* fp;
        fp = fopen(arqOut ,"w+");

        if (fp == NULL)
            printf("Um erro ocorreu ao tentar criar o arquivo.\n");
        
        fprintf(fp,"###########\n");
        fprintf(fp,"# Tempo Método GS: %f\n", tempo/iter);
        fprintf(fp,"#\n");
        fprintf(fp,"# Norma L2 do Residuo\n");
        for (int i=1; i<=iter; ++i)
        {
            fprintf(fp, "# i = %d: %f\n",i,R[i]);
        }
        fprintf(fp,"###########\n");

        fclose(fp);

    }else{
        printf("###########\n");
        printf("# Tempo Método GS: %f\n", tempo/iter);
        printf("#\n");
        printf("# Norma L2 do Residuo\n");
        for (int i=1; i<=iter; ++i)
        {
            printf("# i = %d: %f\n",i,R[i]);
        }      
        printf("###########\n");
    }
}

int main(int argc, char *argv[])
{

    /*=========================== Le os parametros ===========================*/
    int opt;
    char *arqOut = NULL;
    int flagArq = 0, x, y, maxIter;
    FILE *fp;

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
    printf("\nDP:\n");
    for (i = 1; i <= tamLin; i++)
    {
        printf("%.2f ", SL->dp[i]);
    }
    printf("\n\n");

    printf("\nDS:\n");
    for (i = 1; i <= tamLin; i++)
    {
        printf("%.2f ", SL->ds[i]);
    }
    printf("\n\n");

    printf("\nDI:\n");
    for (i = 1; i <= tamLin; i++)
    {
        printf("%.2f ", SL->di[i]);
    }
    printf("\n\n");

    printf("\nDSA:\n");
    for (i = 1; i <= tamLin; i++)
    {
        printf("%.2f ", SL->dsa[i]);
    }
    printf("\n\n");

    printf("\nDIA:\n");
    for (i = 1; i <= tamLin; i++)
    {
        printf("%.2f ", SL->dia[i]);
    }
    printf("\n\n");

    printf("\nB:\n");
    for (i = 1; i <= tamLin; i++)
    {
        printf("%.2f ", SL->b[i]);
    }
    printf("\n\n");

    int iter;
    real_t *R = (real_t *)malloc(maxIter * sizeof(real_t));
    int ret = gaussSeidel(R,SL,maxIter,&iter);

    printf("\nX:\n");
    for (i = 1; i <= tamLin; i++)
    {
        printf("%.2f ", SL->x[i]);
    }
    printf("\n\n");
    real_t n = normaL2Residuo(SL);
    printf("Norma L2: %f\n", n);

    saida_gnuplot(R,n,flagArq,arqOut,tempo,iter);

    /*===================================================================================*/

 /*
  printf("\n--------------------------------------------\n\n");

  printf("#-----------------------------------------\n# Ordem: %d\n", tamLin );
  printf("# Metodo    , ret, norma    ,  iteracoes, tempo\n#-----------------------------------------\n");
*/
    /*double tempo_inicial = timestamp();
  int ret = gaussSeidel(SL,maxIter);
  double tempo_final = timestamp();

  double norma = normaL2Residuo(SL);
/*  
  printf("Gauss Seidel, Norma: %.6g,Tempo: %.5g, Erro: %d\n", norma, tempo_final-tempo_inicial, ret);
*/
    liberaSistLinear(SL);
    return 0;
}
