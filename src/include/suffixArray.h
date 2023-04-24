#ifndef _SAIS_H
#define _SAIS_H 1

#ifdef __cplusplus
extern "C"
{
#endif /* __cplusplus */

    // Creating a Suffix Array
    int sais(const unsigned char *T, int *SA, int *LCP, int n);

#ifdef __cplusplus
} /* extern "C" */
#endif /* __cplusplus */

#endif /* _SAIS_H */