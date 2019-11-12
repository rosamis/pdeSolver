/**
 * \file SistemasLineares.c
 * \author Jackson Borguezani & Roberta Tomigian
 * \mainpage Funções de alocação, inicialiazação e resolução de um sistema linear com o Cálculo da NormaL2
 * \section introSec ReadMe
 * Trabalho 1 | ICC - Prof. Armando Nicolui. |
 * Aluno: Jackson Rossi Borguezani GRR20176573. \n Roberta Tomigian GRR20171631.
 * \brief Arquivo que define os métodos de calcular uma solução para um sistema linear e calcular a norma do resíduo L2.
 *\ date 22 out 2019
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
		tempof = timestamp();
		tempoFim = tempof - tempo;
		somaTempo += tempoFim;
		P->norma[k] = normaL2Residuo(SL);
		++k;
	} while (k <= maxIter);
	P->mediaTempo = somaTempo / k;
	P->iter = k-1;

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
	real_t *Xi, *Bi, Dp, Di, Ds, Dia, Dsa, somatorio, norma, diff, xk;

	nx = SL->nx;
	ny = SL->ny;

	Dp = SL->dp;
	Di = SL->di;
	Ds = SL->ds;
	Dia = SL->dia;
	Dsa = SL->dsa;
	Xi = SL->x;
	Bi = SL->b;
	tam = (SL->nx + 2) * (SL->ny + 2);

	R = (real_t *)malloc(tam * sizeof(real_t));
	
	for (unsigned int i = 0; i < tam; ++i)
		R[i] = 0.0;

	real_t soma = 0.0;

	for(unsigned int i = 1; i < (nx+1); ++i){
		for(unsigned int j = 1; j < (ny+1); ++j){				
			xk = Bi[(i*(ny+2)) + j];
			xk -= (Di * Xi[(i*(ny+2)) + j - 1]);				
			xk -= (Dia * Xi[((i-1)*(ny+2)) + j]);				
			xk -= (Ds * Xi[(i*(ny+2)) + j + 1]);				
			xk -= (Dsa * Xi[((i+1)*(ny+2)) + j]);				
			xk -= (Dp * Xi[(i*(ny+2)) + j]);
			R[(i*(ny+2)) + j] = xk;
			soma += xk*xk;
		}
	}
	free(R);
	return sqrt(soma);
}
