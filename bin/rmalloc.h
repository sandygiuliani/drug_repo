#ifndef __RMALLOC_H
#define __RMALLOC_H

/*
 *  File      : rmalloc.h
 *  Author    : Ronald Kriemann
 *  Purpose   : malloc library 
 */

#include <malloc.h>

/*
 * supported arguments for mallopt
 */

/* extra amount of memory to allocate from system in
   each sbrk/mmap/malloc-call */
#ifndef M_TOP_PAD
#  define M_TOP_PAD  -2
#endif

#ifdef __cplusplus
extern "C" {
#endif

/**************************************
 *
 * 64bit clean definition of mallinfo
 *
 **************************************/

struct ul_mallinfo  {
        unsigned long arena;    /* total space in arena */
        unsigned long ordblks;  /* number of ordinary blocks */
        unsigned long smblks;   /* number of small blocks */
        unsigned long hblks;    /* number of holding blocks */
        unsigned long hblkhd;   /* space in holding block headers */
        unsigned long usmblks;  /* space in small blocks in use */
        unsigned long fsmblks;  /* space in free small blocks */
        unsigned long uordblks; /* space in ordinary blocks in use */
        unsigned long fordblks; /* space in free ordinary blocks */
        unsigned long keepcost; /* cost of enabling keep option */
};
    
/**************************************
 *
 * allocation and deallocation
 *
 **************************************/

extern void * malloc  ( size_t size );
extern void * calloc  ( size_t nmemb, size_t size);
extern void   free    ( void * ptr  );
extern void * realloc ( void * ptr, size_t size );
    
/* set malloc options */
extern int mallopt ( int param, int val );
    
/* report memory usage */
extern struct mallinfo    mallinfo ( void );
extern struct ul_mallinfo rmallinfo ( void );

#ifdef __cplusplus
}
#endif

#endif  /* __RMALLOC_H */
