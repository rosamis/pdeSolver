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
#include "pdeSolver.h"

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
	if (SL)
	{
		SL->b = (real_t *)malloc(((nx+2)*(ny+2)) * sizeof(real_t));
		SL->x = (real_t *)malloc(((nx+2)*(ny+2)) * sizeof(real_t));
		if (!(SL->b) || !(SL->x))
		{
			liberaSistLinear(SL);
		}
	}
	return SL;
}

/*!
	\fn Metrica *alocaMetrica(unsigned int nx, unsigned int ny, int maxIter)
	*
	\brief Aloca e inicia a estrutura de dados métrica
	*
	\param nx: Número de pontos a serem calculados no eixo x
	*
	\param ny: Número de pontos a serem calculados no eixo y
	*
	\param maxIter: Número máximo de iterações
	*
*/
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
	int tam = (x+2) * (y+2);
	for (unsigned int i = 0; i < tam; ++i)
	{
		SL->b[i] = 0.0;
	}
	return;
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
	free(SL->b);
	free(SL->x);
	free(SL);
	return;
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
	\param P: É um ponteiro para uma estrutura de dados que guarda o vetor norma, a média do tempo e o número de iterações de cada iteração do Gauss-Seidel.
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
	real_t xk, *Xi, norma, diff, normaL2;

	nx = SL->nx;
	ny = SL->ny;

	Xi = SL->x;

	double tempo = 0.0;
	double tempof = 0.0;
	double tempoFim = 0.0;
	double somaTempo = 0.0;

	contorno_x(SL);
	k = 1;
	do
	{
		tempo = timestamp();
		//primeira equação
		for(unsigned int i = 1; i < (nx+1); ++i){
			for(unsigned int j = 1; j < (ny+1); ++j){
				xk = SL->b[(i*(ny+2)) + j];
				xk -= ((SL->di) * Xi[(i*(ny+2)) + j - 1]);
				xk -= ((SL->dia) * Xi[((i-1)*(ny+2)) + j]);
				xk -= ((SL->ds) * Xi[(i*(ny+2)) + j + 1]);
				xk -= ((SL->dsa) * Xi[((i+1)*(ny+2)) + j]);
				xk /= (SL->dp);

				Xi[(i*(ny+2)) + j] = xk;
			}
		}

		++k;
		tempof = timestamp();
		tempoFim = tempof - tempo;
		somaTempo += tempoFim;
		//P->norma[k] = normaL2Residuo(SL);

	} while (k < maxIter && norma>EPS);
	P->mediaTempo = somaTempo / k;
	P->iter = k;

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
	real_t xk, *Xk1, *Xi, *Bi, Dp, Di, Ds, Dia, Dsa, somatorio, norma, diff;

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
		R[i] = 0.0;
	//primeira equação
	i = 1;
	R[i] = Bi[i] - (Dp * Xi[i] + Ds * Xi[i + 1] + Dsa * Xi[i + nx]);

	//equações centrais
	for (i = 2; i < tam; ++i)
	{
		if (i > nx)
			R[i] = Bi[i] - (Dp * Xi[i] + Dia * Xi[i - nx] + Di * Xi[i - 1] + Ds * Xi[i + 1] + Dsa * Xi[i + nx]);
		else
			R[i] = Bi[i] - (Dp * Xi[i] + Di * Xi[i - 1] + Ds * Xi[i + 1] + Dsa * Xi[i + nx]);
	}

	//ultima equação
	R[i] = Bi[i] - (Dp * Xi[i] + Dia * Xi[i - nx] + Di * Xi[i - 1]);

	real_t soma = 0.0;
	for (i = 1; i <= tam; ++i)
		soma += R[i] * R[i];
	
	free(R);
	return sqrt(soma);
}
