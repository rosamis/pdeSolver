#ifndef __PDESOLVER_H__
#define __PDESOLVER_H__

/*!
    \fn void gera_matriz(SistLinear_t * SL)
	\brief Função que gera a matriz de coeficientes A
	\param SL é o Sistema Linear
	\details Essa função gera a matriz de coeficientes a partir da aproximação de valores da Equação Diferencial Parcial
	*
*/
void gera_matriz(SistLinear_t * SL);

/*!
    \fn void gera_vetor_b(SistLinear_t * SL)
	\brief Função que gera o vetor de termos independentes b
	\param SL é o Sistema Linear
	\details Essa função gera o vetor de termos independentes b a partir da aproximação de valores da Equação Diferencial Parcial
	*
*/
void gera_vetor_b(SistLinear_t * SL);

#endif //__PDESOLVER_H__