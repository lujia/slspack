
#include "utils.h"

#if TIME_WITH_SYS_TIME
#include <sys/time.h>
#include <time.h>
#else
#if HAVE_SYS_TIME_H
#include <sys/time.h>
#else
#include <time.h>
#endif
#endif

#if USE_UNIX /* linux, unix */
double slspack_get_time(void)
{
    struct timeval tv;

    gettimeofday(&tv, (struct timezone *)0);
    return tv.tv_sec + (double)tv.tv_usec * 1e-6;
}

#else /* windows */
#include <windows.h>

double slspack_get_time()
{
    LARGE_INTEGER timer;
    static LARGE_INTEGER fre;
    static int init = 0;

    if (init != 1) {
        QueryPerformanceFrequency(&fre);
        init = 1;
    }

    QueryPerformanceCounter(&timer);

    return (double)1. * timer.QuadPart / fre.QuadPart;
}
#endif

int slspack_printf(const char *fmt, ...)
{
    va_list ap;
    int ret;

    va_start(ap, fmt);
    ret = vfprintf(stdout, fmt, ap);
    va_end(ap);

    fflush(stdout);

    return ret;
}

void slspack_error(int code, const char *fmt, ...)
{
    va_list ap;
    char s[2048];

    sprintf(s, "error: %s", fmt);
    va_start(ap, fmt);
    vfprintf(stdout, s, ap);
    va_end(ap);

    fflush(stdout);

    if (code == 0) return;
    exit(code);
}

void slspack_warning(const char *fmt, ...)
{
    va_list ap;
    char s[2028];

    sprintf(s, "warning: %s", fmt);
    va_start(ap, fmt);
    vfprintf(stdout, s, ap);
    va_end(ap);

    fflush(stdout);

    return;
}

void * slspack_malloc(size_t n) 
{ 
    void *ptr;

    if (n == 0) return NULL;

    ptr = malloc(n);

    if (ptr == NULL) {
        slspack_error(1, "slspack: failed to malloc memory\n");
    }

    return ptr;
}

void * slspack_calloc(size_t n) 
{ 
    void *ptr;

    if (n == 0) return NULL;

    ptr = calloc(1, n);

    if (ptr == NULL) {
        slspack_error(1, "slspack: failed to malloc memory\n");
    }

    return ptr;
}

void slspack_free(void *p) 
{ 
    if (p == NULL) return;

    free(p);
}

void slspack_memcpy(void *dst, const void *src, size_t n)
{
    if (n == 0) return;

    assert(dst != NULL && src != NULL);
    memcpy(dst, src, n);
}

void * slspack_mem_copy(const void *src, size_t n)
{
    void *dst;

    if (n == 0) return NULL;

    assert(src != NULL);

    dst = slspack_malloc(n);
    slspack_memcpy(dst, src, n);

    return dst;
}
