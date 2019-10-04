/**
 * \file SistemasLineares.c
 * \author Jackson Borguezani & Roberta Tomigian
 * \mainpage Funções de alocação, inicialiazação e resolução de um sistema linear com o Cálculo da NormaL2
 * \section introSec ReadMe
 * Trabalho 1 | ICC - Prof. Armando Nicolui. |
 * Aluno: Jackson Rossi Borguezani GRR20176573.
 * \brief Arquivo que define os métodos de calcular uma solução para um sistema linear e calcular a norma do resíduo L2.
 *\ date 03 out 2019
 *
 */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "SistemasLineares.h"
#include "utils.h"

/*!
	\fn SistLinear_t *alocaSistLinear(unsigned int nx, unsigned int ny)
	*
	\brief Alocação do Sistema Linear Pentadiagonal
	*
	\param nx: Número de pontos a serem calculados no eixo x
	*
	\param ny: Número de pontos a serem calculados no eixo y
	*
	\details Função responsável por alocar dinamicamente os vetores do SL.
	*
*/
SistLinear_t *alocaSistLinear(unsigned int nx, unsigned int ny)
{
	SistLinear_t *SL = (SistLinear_t *)malloc(sizeof(SistLinear_t));
	int t = nx * ny;
	if (SL)
	{
		SL->dp = (real_t *)malloc(t * sizeof(real_t));
		SL->di = (real_t *)malloc(t * sizeof(real_t));
		SL->ds = (real_t *)malloc(t * sizeof(real_t));
		SL->dia = (real_t *)malloc(t * sizeof(real_t));
		SL->dsa = (real_t *)malloc(t * sizeof(real_t));
		SL->b = (real_t *)malloc((t + 1) * sizeof(real_t));
		SL->x = (real_t *)malloc((t + 1) * sizeof(real_t));
		if (!(SL->dp) || !(SL->b) || !(SL->x) || !(SL->di) || !(SL->dia) || !(SL->dsa) || !(SL->ds))
		{
			liberaSistLinear(SL);
		}
	}
	return SL;
}

Metrica *alocaMetrica(unsigned int nx, unsigned int ny, int maxIter){
	Metrica *P = (Metrica *)malloc(sizeof(Metrica));
	int t = maxIter;
	if(P){
		P->norma = (double*)malloc(t * sizeof(double));
		if(!(P->norma)){
			free(P);
		} 
		for(int i = 0; i <= t; ++i){
			P->norma[i] = 0.0;
		}
		P->iter = 0;
		P->mediaTempo = 0.0;
	}
	return P;
}

/*!
	\fn void inicializaSistLinear(SistLinear_t *SL, int x, int y)
	*
	\brief Inicialização do Sistema Linear com zeros
	*
	\param SL: Ponteiro do Sistema Linear
	*
	\param x: Número de pontos a serem calculados no eixo x
	*
	\param y: Número de pontos a serem calculados no eixo y
	*
*/
void inicializaSistLinear(SistLinear_t *SL, int x, int y)
{
	SL->nx = x;
	SL->ny = y;
	int tam = x * y;
	for (unsigned int i = 0; i < tam; ++i)
	{
		SL->dp[i] = 0.0;
		SL->ds[i] = 0.0;
		SL->di[i] = 0.0;
		SL->dia[i] = 0.0;
		SL->dsa[i] = 0.0;
		SL->b[i] = 0.0;
		SL->x[i] = 0.0;
	}
}

/*!
	\fn void liberaSistLinear(SistLinear_t *SL)
	*
	\brief Liberação dos Ponteiros do Sistema Linear.
	*
	\param SL: Ponteiro do Sistema Linear
	*
	\details A função libera todos os ponteiros da estrutura do Sistema Linear.
	*
*/
void liberaSistLinear(SistLinear_t *SL)
{
	free(SL->dp);
	free(SL->ds);
	free(SL->di);
	free(SL->dia);
	free(SL->dsa);
	free(SL->b);
	free(SL->x);
	free(SL);
}

/*!
	\fn int gaussSeidel(SistLinear_t *SL, int maxIter)
	*
	\brief Método da Gauss-Seidel
	*
	\param SL: Ponteiro do Sistema Linear
	*
	\param maxIter: Número máximo de iterações
	*
	\details Esta função calcula a solução de um sistema linear pentadiagonal, com estrutura de dados em vetores.
	*
	\details Erros possíveis no metodo de Gauss-Seidel
	*
	· Erro (-1) DIVISAO_POR_0: Caso a diagonal principal seja zero, não é possível calcular a solução pois teremos uma divisao por 0
	*
	· Erro (-2) ATINGIU_MAX_ITERACOES: Ocorre quando o número máximo de iterações, passado por parâmetro, é atingido
	*
*/
int gaussSeidel(SistLinear_t *SL, int maxIter, Metrica *P)
{
	int i, j, k, l, tam, nx, ny;
	real_t xk, *Xi, *Bi, *Dp, *Di, *Ds, *Dia, *Dsa, somatorio, norma, diff, normaL2;

	nx = SL->nx;
	ny = SL->ny;

	Dp = SL->dp;
	Di = SL->di;
	Ds = SL->ds;
	Dia = SL->dia;
	Dsa = SL->dsa;
	Xi = SL->x;
	Bi = SL->b;
	tam = SL->nx * SL->ny;

	double tempo=0.0;
	double somaTempo = 0.0;

	k = 1;
	do
	{
		tempo = timestamp();
		//primeira equação
		i = 1;
		if (Dp[i] == 0)
		{
			fprintf(stderr, "Erro (-1)\n");
			return -1;
		}
		xk = (Bi[i] - Ds[i] * Xi[i + 1] - Dsa[i] * Xi[i + nx]) / Dp[i];
		norma = fabs(xk - Xi[1]);
		//printf("N prim%d: %.3f\n", i, norma);
		Xi[i] = xk;
		//printf("\nN1: %.3f\n",norma);

		//equações centrais
		for (i = 2; i < tam; ++i)
		{
			/*if (Dp[i] == 0)
			{
				fprintf(stderr, "Erro (-1)\n");
				return -1;
			}*/
			if (i > nx)
			{
				xk = (Bi[i] - Dia[i] * Xi[i - nx] - Di[i] * Xi[i - 1] - Ds[i] * Xi[i + 1] - Dsa[i] * Xi[i + nx]) / Dp[i];
			}
			else
			{
				xk = (Bi[i] - Di[i] * Xi[i - 1] - Ds[i] * Xi[i + 1] - Dsa[i] * Xi[i + nx]) / Dp[i];
			}
			diff = fabs(xk - Xi[i]);
			norma = (diff > norma) ? (diff) : (norma);
			Xi[i] = xk;
			//printf("N%d: %.3f\n",i,norma);
		}

		//ultima equação
		xk = (Bi[i] - Dia[i] * Xi[i - nx] - Di[i] * Xi[i - 1]) / Dp[i];
		diff = fabs(xk - Xi[i]);
		norma = (diff > norma) ? (diff) : (norma);
		Xi[i] = xk;
		//printf("Nfora: %.3f\n", norma);

		++k;
		double tempo1 = timestamp();
		double tempoFim = tempo1 - tempo;
		somaTempo += tempoFim;
		//P->norma[k] = normaL2Residuo(SL);

	} while (k < maxIter && norma > EPS);
	//printf("K: %d\n", k);
	//printf("N fim: %.3f\n", norma);
	double mTempo = somaTempo / k;
	P->mediaTempo = mTempo;
	P->iter = k;
	printf("Tempo média: %f\n",mTempo );


	if (k >= maxIter)
	{
		fprintf(stderr, "Erro (-2)\n");
		return -2;
	}

	return 0;
}

/*!
	\fn double normaL2Residuo(SistLinear_t * SL)
	*
	\brief Cálculo da NormaL2
	*
	\param SL: Ponteiro do Sistema Linear
	*
*/
double normaL2Residuo(SistLinear_t *SL)
{
	real_t *R;

	int i, j, k, l, tam, nx, ny;
	real_t xk, *Xk1, *Xi, *Bi, *Dp, *Di, *Ds, *Dia, *Dsa, somatorio, norma, diff;

	nx = SL->nx;
	ny = SL->ny;

	Dp = SL->dp;
	Di = SL->di;
	Ds = SL->ds;
	Dia = SL->dia;
	Dsa = SL->dsa;
	Xi = SL->x;
	Bi = SL->b;
	tam = SL->nx * SL->ny;

	R = (real_t *)malloc(tam * sizeof(real_t));

	for (unsigned int i = 1; i <= tam; ++i)
	{
		R[i] = 0.0;
	}
	//primeira equação
	i = 1;
	R[i] = Bi[i] - (Dp[i] * Xi[i] + Ds[i] * Xi[i + 1] + Dsa[i] * Xi[i + nx]);
	//printf("\nR1: %.3f\n",R[i]);

	//equações centrais
	for (i = 2; i < tam; ++i)
	{
		if (i > nx)
		{
			R[i] = Bi[i] - (Dp[i] * Xi[i] + Dia[i] * Xi[i - nx] + Di[i] * Xi[i - 1] + Ds[i] * Xi[i + 1] + Dsa[i] * Xi[i + nx]);
		}
		else
		{
			R[i] = Bi[i] - (Dp[i] * Xi[i] + Di[i] * Xi[i - 1] + Ds[i] * Xi[i + 1] + Dsa[i] * Xi[i + nx]);
		}
		//printf("R%d: %.3f\n",i,R[i]);
	}

	//ultima equação
	R[i] = Bi[i] - (Dp[i] * Xi[i] + Dia[i] * Xi[i - nx] + Di[i] * Xi[i - 1]);
	//printf("\nR%d: %.3f\n",i,R[i]);

	real_t soma = 0.0;
	for (i = 1; i <= tam; ++i)
	{
		soma += R[i] * R[i];
	}
	//free(R);

	return sqrt(soma);
}