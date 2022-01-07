#ifndef POLYMUL_H
#define POLYMUL_H

#include "SABER_params.h"
#include <stdint.h>

void toom_cook_4way_interpol(const uint16_t *a_i, uint16_t *res_i);
void toom_cook_4way_evaluate(const uint16_t *a1, const uint16_t sw[7][64], uint16_t *res);
void evaluation_single(const uint16_t *b, uint16_t sw[7][64]);

#endif
