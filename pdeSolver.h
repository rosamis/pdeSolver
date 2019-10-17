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
void gera_matriz(SistLinear_t * SL);

/*=========================== Funções que geram o vetor de termos independentes ===========================*/
void gera_vetor_b(SistLinear_t * SL);
void contorno_x(SistLinear_t *SL);

/*================================== Funções que geram a saída da solução ==================================*/
void saida_gnuplot(Metrica *P, int flagArq, char *arqOut, SistLinear_t * SL);

#endif //__PDESOLVER_H__