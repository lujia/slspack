
#ifndef SLSPACK_UTILS_H
#define SLSPACK_UTILS_H

#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include <assert.h>
#include <math.h>

#include "config.h"

/* timer */
double slspack_get_time(void);

/* output to stdout */
int slspack_printf(const char *fmt, ...);

/* output error */
void slspack_error(int code, const char *fmt, ...);

/* print warning info */
void slspack_warning(const char *fmt, ...);

void * slspack_malloc(size_t n);
void * slspack_calloc(size_t n);
void slspack_free(void *p);

void slspack_memcpy(void *dst, const void *src, size_t n);
void * slspack_mem_copy(const void *src, size_t n);

#endif
