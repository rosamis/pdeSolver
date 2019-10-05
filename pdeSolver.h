#ifndef __PDESOLVER_H__
#define __PDESOLVER_H__
/**
 * \file pdeSolver.h
 * \author Jackson Borguezani & Roberta Tomigian
 * \mainpage Biblioteca de Funções e Estruturas da Solução Discreta para uma Equação Diferencial Parcial
 * \section introSec ReadMe
 *
 * \brief **Introdução à Computação Científica CI1164** \n Trabalho 1 -- Prof. Armando Nicolui \n **Alunos :** \n Jackson Rossi Borguezani GRR20176573 \n Roberta Tomigian GRR20171631
 * 
 * \date 03 out 2019
 *
 */

/*======================================== Funções que geram a matriz ======================================*/
real_t quadrado(real_t p);
real_t gera_d_principal(real_t hx_quadrado, real_t hy_quadrado, real_t pi_quadrado);
real_t gera_d_superior(real_t hx, real_t hy_quadrado);
real_t gera_d_inferior(real_t hx, real_t hy_quadrado);
real_t gera_ds_afastada(real_t hx_quadrado, real_t hy);
real_t gera_di_afastada(real_t hx_quadrado, real_t hy);
void gera_matriz(SistLinear_t * SL);

/*=========================== Funções que geram o vetor de termos independentes ===========================*/
real_t funcao_f(int i, int j, real_t pi_quadrado);
real_t b_principal(int i, int j, real_t hx_quadrado, real_t hy_quadrado, real_t pi_quadrado);
real_t b_limite_superior(real_t hx_quadrado, real_t hy, int i, real_t pi_quadrado);
real_t b_limite_inferior(real_t hx_quadrado, real_t hy, int i, real_t pi_quadrado);
real_t b_limite_esquerda(real_t hx, real_t hy_quadrado);
real_t b_limite_direita(real_t hx, real_t hy_quadrado);
void gera_vetor_b(SistLinear_t * SL);

/*================================= Funções que geram a solução da equação =================================*/
void solucao(real_t hx, real_t hy, real_t *x, real_t *y, int ny, int nx);
real_t* aloca_vetor(int nx, int ny);

/*================================== Funções que geram a saída da solução ==================================*/
void saida_gnuplot(Metrica *P, int flagArq, char *arqOut, SistLinear_t * SL);

#endif //__PDESOLVER_H__