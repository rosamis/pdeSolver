/**
 * \file SistemasLineares.h
 * \author Jackson Borguezani & Roberta Tomigian
 * \mainpage Biblioteca de Funções e Estruturas do Sistema Linear
 * \section introSec ReadMe
 * Trabalho 1 | ICC - Prof. Armando Nicolui. |
 * Aluno: Jackson Borguezani GRR20176573 | Roberta Tomigian GRR20171631
 * \brief Arquivo que define as funções e estruturas do sistema linear.
 *\ date 03 out 2019
 *
 */
#include <time.h>
#include <sys/time.h>
#include <math.h>
#include "utils.h"

typedef double real_t;

#define EPS 1.0e-4

/*! \struct SistLinear_t
    \brief Estrutura de dados do Sistema Linear Pentadiagonal
*/
typedef struct {
  real_t *dp; /**< Vetor da diagonal principal. */
  real_t *ds; /**< Vetor da diagonal superior. */
  real_t *di; /**< Vetor da diagonal inferior. */
  real_t *dia; /**< Vetor da diagonal inferior afastada. */
  real_t *dsa; /**< Vetor da diagonal superior afastada. */
  real_t *b; /**< Vetor de termos independentes. */
  real_t *x; /**< Vetor solução. */
  unsigned int nx; /**< Quantidade de pontos em x. */
  unsigned int ny; /**< Quantidade de pontos em y. */
} SistLinear_t;

SistLinear_t* alocaSistLinear (unsigned int nx, unsigned int ny);
void inicializaSistLinear (SistLinear_t *SL, int x, int y);
void liberaSistLinear (SistLinear_t *SL);
int gaussSeidel (SistLinear_t *SL, int maxIter);
double normaL2Residuo(SistLinear_t *SL);