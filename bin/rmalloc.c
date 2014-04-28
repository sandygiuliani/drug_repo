/*
 * RMALLOC by Ronald Kriemann
 * --------------------------
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation
 * (version 2.1 of the License).
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 *
 *
 * Changelog
 * ---------
 *
 * v1.1
 *         - implemented deallocation of empty data segments if enough
 *           memory is available
 *         - non-empty wilderness are now handled by middle-sized
 *           management ensuring good-fit allocation on wilderness
 *           (greatly reduces internal fragmentation in certain situations)
 *         - support for 16 byte alignment added (for SSE applications)
 *         - fixed bug when setting size of split blocks
 *         - modified message and error output
 *         - changed Rmalloc tracing to write more data
 *         - removed STAT_FRAGMENTATION (was not used anymore)
 *         - guarded munmap by mutex in mt-mode (seems to be not mt-safe)
 *         - non-handled (large) blocks are now part of statistics
 * v1.0.3
 *         - changed license to LGPL
 *         - ported to Win32
 *         - guarded special functions by defines (mmap, getpagesize, etc.)
 *         - moved trace.fd to trace.file (file IO instead of read/write)
 *         - C++ operators now really throw std::bad_alloc
 * v1.0.2
 *         - handling of NULL pointers in allocation improved;
 *           no more ABORT but NULL-returns
 * v1.0.1
 *         - in check_rmalloc: global heap was not checked
 *         - consistency check of wilderness, full/non-empty containers
 *         - alignment of size for heap-allocation
 * v1.0
 *         - removed bug in r_realloc (wrong size for given block)
 * v1.0-pre5
 *         - NO_OF_CLASSES, NO_OF_SMALL must be even; changed
 * v1.0-pre4
 *         - reduced internal fragmentation for small chunks by
 *           using a special type (sblock_t)
 *         - minor bugfixes for non-pthread mode
 * v1.0-pre3
 *         - added USE_THR_HACKS to simulate thread-specific-data
 *         - restructured container-management; the empty-list not
 *           needed anymore
 * v1.0-pre2
 *         - removed needless data-segment-structure
 *         - removed needless empty-lists for containers
 * v1.0-pre1
 *         - another major rewrite: added coalescing similar to
 *           dmalloc but without caching (or deffered)
 *         - maintained container-managment for small blocks
 *           and coalescing strategy for middle sized blocks
 *         - removed rounding up of size-requests
 *         - no longer memory-transfers between heaps
 * v0.99   - major rewrite with containers which can be reused by
 *           other size-classes => led to massive fragmentation
 *         - code-cleanup
 * v0.98.2 - replaced "real" thread-private-data by simulated one
 * v0.98.1 - active fragmentation statistics with STAT_FRAGMENTATION
 *           (prints "real" consumption)
 *         - moved top_pad into data-segment-structure
 * v0.98   - used remaining blocks of each data-segment to serve future
 *           requests, thereby reducing fragmentation (left-over-blocks)
 * v0.97   - small bugfixes
 *         - top_pad is now heap-specific and adjusted dynamicly
 *           (1/TOP_PAD_FRACTION of mem_used)
 *         - changed size-field in nodes to slot; avoided to many calls
 *           to size_to_class
 *         - moved handling of sizes larger than biggest size-class to
 *           operating-system (via malloc or mmap)
 *         - fixed error-handling if no memory is available from the OS
 *           (sends kill to own process instead of exit(1))
 *         - added another trace-method: trace by allocation
 * v0.96   - rewritten chunk-handling: heap -> block -> chunk
 *         - exchange of blocks between thread-heaps and global-heaps
 *           to reduce overall memory-consumption
 * v0.95   - round sizes to next size-class to guarantee, that there
 *           is a free chunk in the list (avoid searching)
 *         - increased number of size-classes; changed size_to_class to
 *           use a binary search instead of linear
 *         - ok. maybe 10 MB for the DEFAULT_TOP_PAD is better ?
 *         - implemented trace-functionality for thread-heaps
 * v0.94   - rearranged code for system-allocation (mmap and malloc had
 *           too much in common to be separated)
 *         - removed most messages when mutices were locked
 *         - set MMAP to be default, even when malloc is not overwritten
 *         - changed "r_mallinfo" to "rmallinfo" and added declaration
 *           in rmalloc.h
 * v0.93   - changed heap handling, now each thread has a private
 *           heap; if the thread terminates, this heap is then marked
 *           as free and can be used by a new thread
 *         - removed creation of first heap (besides global-heap), this
 *           will now be done on demand (in r_malloc)
 *         - DEFAULT_TOP_PAD reduced to 5 MB (should also be enough)
 *         - heap-list is no longer a circular list (just a linear one)
 *         - heaps can now be aligned to the cache-line-size (ALIGN_TO_CACHE_SIZE)
 *         - added wrapper for mallinfo which calls internal r_mallinfo
 *           (r_mallinfo uses mallinfo-struct with ulong fields)
 * v0.92   - replaced simulated thread-specific-data by real one
 * v0.91   - fixed bug/feature when using malloc as system-alloc: only
 *           requested chunks of needed size were allocated, not plus
 *           DEFAULT_TOP_PAD
 * v0.90   - initial release
 *
 * ToDo
 * ----
 *
 * - check release of heaps during thread-cancelation, maybe reprogramm
 *   thread-private data
 * - add support for mremap
 */


/*
 * the version info
 */

#define RMALLOC_VERSION  "1.1"

/*****************************************************************
 *****************************************************************
 **
 ** some options
 **
 *****************************************************************
 *****************************************************************/

/*
 * use rmalloc as a malloc-replacment
 */
#if defined(__cplusplus) && !defined(OVERLOAD_MALLOC)
#  define OVERLOAD_MALLOC  0
#else
#  ifndef OVERLOAD_MALLOC
#  define OVERLOAD_MALLOC  1
#  endif
#endif

/*
 * use rmalloc as a new/delete replacment
 */
#if defined(__cplusplus) && !defined(OVERLOAD_NEW)
#  define OVERLOAD_NEW     1
#endif

/*
 * use pthreads
 */
#ifndef USE_PTHREAD
#  define USE_PTHREAD      1
#endif

/*
 * set to 1 to simulate thread-private
 * data
 */
#if OVERLOAD_MALLOC == 1 && ! defined(USE_THR_HACKS)
#  warning "rmalloc: using thread-hacks when overloading malloc\n"
#  define USE_THR_HACKS    1
#endif

#ifndef USE_THR_HACKS
#  define USE_THR_HACKS    0
#endif

/*
 * define system-allocation method: MALLOC, MMAP, SBRK
 */
#if OVERLOAD_MALLOC == 1
#  if !defined(USE_MMAP) && !defined(USE_SBRK)
#    if HAS_MMAP == 1
#      define USE_MMAP         1
#      define USE_SBRK         0
#    else
#      define USE_MMAP         0
#      define USE_SBRK         1
#    endif
#  endif
#else
#  if !defined(USE_MALLOC) && !defined(USE_MMAP) && !defined(USE_SBRK)
#    if HAS_MMAP == 1
#      define USE_MALLOC       0
#      define USE_MMAP         1
#      define USE_SBRK         0
#    else
#      define USE_MALLOC       1
#      define USE_MMAP         0
#      define USE_SBRK         0
#    endif
#  endif
#endif

#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>
#include <assert.h>
#if HAS_GETPAGESIZE == 1
#include <unistd.h>
#endif
#include <sys/types.h>
#include <sys/stat.h>
#if HAS_MMAP == 1
#include <sys/mman.h>
#endif
#include <fcntl.h>
#include <string.h>
#include <ctype.h>
#include <errno.h>

#if USE_PTHREAD == 1
#  include <pthread.h>
#endif

#include "rmalloc.h"

#ifdef __cplusplus
extern "C" {
#endif

/*****************************************************************
 *****************************************************************
 **
 ** defines and option setting
 **
 *****************************************************************
 *****************************************************************/

/* overload malloc functions */
#if OVERLOAD_MALLOC == 1
#  define  r_malloc          malloc
#  define  r_calloc          calloc
#  define  r_free            free
#  define  r_realloc         realloc
#  define  r_mallopt         mallopt
#endif

/* alignment of all memory (8 or 16 bytes) */
#define ALIGNMENT            8
#define IS_ALIGNED( s )      (((s) % ALIGNMENT) == 0)
#define IS_MISALIGNED( s )   (((s) % ALIGNMENT) != 0)
#define ALIGN( s )           ((s) += (ALIGNMENT - ((s) % ALIGNMENT)))
    
/*
 * total number of size classes,
 * number of small classes and small size
 * and first "middle" size
 */
#define NO_OF_CLASSES        222
#if ALIGNMENT == 8
#  define SMALL_SIZE         256
#  define MIN_BLOCK_SIZE     288
#else
#  define SMALL_SIZE         512
#  define MIN_BLOCK_SIZE     576
#endif
#define NO_OF_SMALL          (SMALL_SIZE / ALIGNMENT)

/* size of a container for small blocks and
 * maximal number of full cont. per heap */
#define CONTAINER_SIZE       8192
#define MAX_FULL_CONTAINERS  8

/* container-classes */
#define TYPE_FULL            0
#define TYPE_NONEMPTY        1
#define TYPE_EMPTY           2

/*
 * determine bitsize
 * (actually this is just a guess and hopefully covers all used defines)
 */
#if defined(__x86_64__)    || defined(__ia64__)   || defined(__powerpc64__) || \
    defined(__sparc_v9__)  || defined(__sparcv9)  || defined(__arch64__)    || \
    defined(__LP64__)      || defined(_LP64)
#  define BITSIZE 64
#else
#  define BITSIZE 32
#endif

/*
 * set to 1 to use one heap per thread which potentially increases
 * the parallel speedup, otherwise a new heap is only created if all
 * others are currently used
 * WARNING: ONE_HEAP_PER_THREAD == 0 is not working
 */
#define ONE_HEAP_PER_THREAD  1
    
/*
 * set to 1 to use binary search for size classes
 */
#define USE_BIN_SEARCH       1
    
/*
 * set to 1 to sort freed middle sized blocks
 *  - increases time per free but potentially decreases fragmentation
 */
#define SORT_MIDDLE_BLOCKS   0
    
/*
 * set to 1 to sort wilderness blocks in heap
 *  - might increases time per alloc/free but potentially decreases fragmentation
 */
#define SORT_WILDERNESS      0
    
/*
 * set to 1 to enable debugging
 */
#define DEBUG_BLOCK          0
#define DEBUG_CONT           0
#define DEBUG_RMALLOC        0

/* magic numbers for blocks and containers */
#define BLOCK_MAGIC          0xbabecafe
#define CONT_MAGIC           0xcafebabe
    
/* default trace file */
#define DEFAULT_TRACE_FILE   "rmalloc.trc"

/* trace-types: 1 : trace by steps, 2: trace by allocation */
#define TRACE_STEPS          1
#define TRACE_ALLOCATION     2
    
/* default value for top-pad in bytes */
#define DEFAULT_TOP_PAD      1*1024*1024

/* fraction of the heap, top_pad must not exceed */
#define TOP_PAD_FRACTION     50

/* exit command */
#define EXIT                 exit( 1 )

/* error value in sbrk */
#define SBRK_FAIL            -1

/* definitions for mmap-usage */
#if USE_MMAP == 1

#if !defined(MAP_ANONYMOUS) && defined(MAP_ANON)
#  define MAP_ANONYMOUS      MAP_ANON
#endif

#if !defined(MAP_FAILED)
#  define MAP_FAILED         ((char*)-1)
#endif

#if !defined(MAP_ANONYMOUS)

static int dev_zero_fd = -1;

#define MMAP(addr, size, prot, flags) ((dev_zero_fd < 0) ? \
 (dev_zero_fd = open("/dev/zero", O_RDWR), \
  mmap((addr), (size), (prot), (flags), dev_zero_fd, 0)) : \
   mmap((addr), (size), (prot), (flags), dev_zero_fd, 0))

#else

#define MMAP(addr, size, prot, flags) \
 (mmap((addr), (size), (prot), (flags)|MAP_ANONYMOUS, -1, 0))

#endif  /* !defined(MAP_ANONYMOUS) */

#endif  /* USE_MMAP == 1 */

/*****************************************************************
 *****************************************************************
 **
 ** data-types for rmalloc
 **
 *****************************************************************
 *****************************************************************/


/* internal size and pointer types */
typedef  unsigned long  internal_size_t;

/* bool looks better than int */
#ifndef __cplusplus
typedef enum { false, true } bool;
#endif

/*
 * mutex and thread-specific-data type
 */
#if USE_PTHREAD == 1
typedef pthread_mutex_t  mutex_t;
typedef pthread_key_t    tsd_key_t;
#else
typedef int              mutex_t;
typedef int              tsd_key_t;
#endif

/*
 * forward declaration of types
 */

struct sysconf_s;

struct sblock_s;
struct block_s;
struct container_s;
struct block_storage_s;
struct heap_s;
    
struct heaplist_s;
struct trace_s;

typedef struct sysconf_s   sysconf_t;
typedef struct stat_s      stat_t;
    
typedef struct sblock_s        sblock_t;
typedef struct block_s         block_t;
typedef struct container_s     container_t;
typedef struct block_storage_s block_storage_t;
typedef struct heap_s          heap_t;

typedef struct heap_list_s heap_list_t;
typedef struct trace_s     trace_t;

/*
 * for holding info about system
 */
struct sysconf_s
{
    internal_size_t   pagesize;  /* number of bytes per page                        */
    internal_size_t   top_pad;   /* number of extra bytes to allocated in each sbrk */
};

/*
 * holding statistical information
 */
struct stat_s
{
    internal_size_t  used[ NO_OF_CLASSES ];  /* number of used chunks */
    internal_size_t  free[ NO_OF_CLASSES ];  /* number of free chunks */
    
    internal_size_t  used_mem;     /* total number of bytes used in heap (alloc from system) */
    internal_size_t  max_mem;      /* maximal number of bytes allocated from system          */
    internal_size_t  mem_in_use;   /* current number of bytes, used by application           */
    internal_size_t  max_in_use;   /* maximal number of bytes, used by application           */
    internal_size_t  allocated;    /* sum of all allocated memory in heap                    */
};

/*
 * datatype for small memory chunks
 */
struct sblock_s
{
    internal_size_t  size;         /* size of block               */
    void           * hc_addr;      /* heap/container we belong to */

#if ALIGNMENT == 8
#  if DEBUG_BLOCK == 1
    unsigned         magic;        /* magic number for the block  */
    unsigned         __pad__;      /* for 8-byte-alignment        */
#  endif
#else /* ALIGNMENT == 16 */
#  if BITSIZE == 64
#    if DEBUG_BLOCK == 1
    unsigned         magic;        /* magic number for the block  */
    unsigned         __pad__[3];   /* for 16-byte-alignment       */
#    endif
#  else
#    if DEBUG_BLOCK == 1
    unsigned         magic;        /* magic number for the block  */
    unsigned         __pad__;      /* for 16-byte-alignment       */
#    else
    unsigned         __pad__[2];   /* for 16-byte-alignment       */
#    endif
#  endif
#endif
    
    sblock_t       * next;         /* next, prev pointer for list */
};

/*
 * datatype for memory chunks
 */
struct block_s
{
    internal_size_t  size;         /* size of block               */
    void           * hc_addr;      /* heap/container we belong to */

#if ALIGNMENT == 8
#  if DEBUG_BLOCK == 1
    unsigned         magic;        /* magic number for the block  */
    unsigned         __pad__;      /* for 8-byte-alignment        */
#  endif
#else /* ALIGNMENT == 16 */
#  if BITSIZE == 64
#    if DEBUG_BLOCK == 1
    unsigned         magic;        /* magic number for the block  */
    unsigned         __pad__[3];   /* for 16-byte-alignment       */
#    endif
#  else
#    if DEBUG_BLOCK == 1
    unsigned         magic;        /* magic number for the block  */
    unsigned         __pad__;      /* for 16-byte-alignment       */
#    else
    unsigned         __pad__[2];   /* for 16-byte-alignment       */
#    endif
#  endif /* ALIGNMENT */
#endif
    
    block_t        * next, * prev; /* next, prev pointer for list       */
};

/* size of extra-data to be stored in each chunk (pointers + size-field) */
#if ALIGNMENT == 8
#  if DEBUG_BLOCK == 1
#    define PREFIX_SIZE       (sizeof(internal_size_t) + sizeof(void*) + 2*sizeof(unsigned))
#  else
#    define PREFIX_SIZE       (sizeof(internal_size_t) + (sizeof(void*)))
#  endif
#else /* ALIGNMENT == 16 */
#  if BITSIZE == 64
#    if DEBUG_BLOCK == 1
#      define PREFIX_SIZE     (sizeof(internal_size_t) + sizeof(void*) + 4*sizeof(unsigned))
#    else
#      define PREFIX_SIZE     (sizeof(internal_size_t) + sizeof(void*))
#    endif
#  else
#    if DEBUG_BLOCK == 1
#      define PREFIX_SIZE     (sizeof(internal_size_t) + sizeof(void*) + 2*sizeof(unsigned))
#    else
#      define PREFIX_SIZE     (sizeof(internal_size_t) + sizeof(void*) + 2*sizeof(unsigned))
#    endif
#  endif
#endif /* ALIGNMENT */
#define POSTFIX_SIZE    (sizeof(internal_size_t))
#define EXTRA_DATA_SIZE (PREFIX_SIZE)

/*
 * container for a list of blocks
 */
struct container_s
{
    heap_t       * heap;               /* heap holding this container               */
    container_t  * next, * prev;       /* data for double-linked-list               */
    sblock_t     * blocks;             /* list of free memory-chunks in container   */
    unsigned       free, max;          /* number of free and total number of blocks */
    unsigned       sclass;             /* which size-class represents the container */
#if DEBUG_CONT == 1
    unsigned       magic;              /* magic number for the container            */
#else
    unsigned       __pad__;            /* for alignment                             */
#endif
    
};

/*
 * combines different container-types
 */
struct block_storage_s
{
    container_t  * empty[ NO_OF_SMALL ];    /* all blocks in container used         */
    container_t  * nonempty[ NO_OF_SMALL ]; /* at least one block used, but not all */
    container_t  * full, * last_full;       /* full container, no block used        */
    unsigned       no_of_full;              /* number of full containers            */
    block_t      * blocks[ NO_OF_CLASSES ]; /* size-classes for usual blocks        */
    block_t      * wilderness;              /* empty data-segments for future req.  */
};

/*
 * type for a heap
 */
struct heap_s
{
    int              id;                    /* unique id of heap                    */
    heap_t         * next;                  /* link for the heap-list               */
    block_storage_t  blocks;                /* holds all stored blocks              */
    stat_t           stat;                  /* status information                   */
    internal_size_t  top_pad;               /* number of bytes to alloc. in advance */
    mutex_t          mutex;                 /* mutex for locking heap               */
    bool             used;                  /* if false, heap is not used           */
};

/*
 * list of heaps for the threads
 */
struct heap_list_s
{
    heap_t   * heaps;                       /* linked list of heaps    */
    unsigned   nheaps;                      /* number of heaps in list */
    mutex_t    mutex;                       /* mutex for safety        */
};

/*
 * type for tracing memory-allocation
 */
struct trace_s
{
    FILE            * file;                 /* handle for trace-file                  */
    internal_size_t   old_rm_used;          /* old number of used bytes in Rmalloc    */
    internal_size_t   old_app_used;         /* old number of used bytes in app        */
    internal_size_t   old_free;             /* old number of free bytes               */
    internal_size_t   old_alloc;            /* old total allocation size (time in MB) */
    unsigned          step;                 /* allocation step                        */
    internal_size_t   size;                 /* size of chunk to trace                 */
    unsigned          type;                 /* type of trace                          */
};



/*****************************************************************
 *****************************************************************
 **
 ** thread definition
 **
 *****************************************************************
 *****************************************************************/

#if USE_PTHREAD == 1

/* maximal number of threads */
#define TSD_MAX   256

/*
 * handling of mutices
 */
#  define MUTEX_INITIALIZER  PTHREAD_MUTEX_INITIALIZER
#  define DEFINE_MUTEX( m )  static mutex_t m = PTHREAD_MUTEX_INITIALIZER;

#  define MUTEX_INIT(mutex,status)  if (( status = pthread_mutex_init( & mutex, NULL )) != 0) \
                                    { ERROR( __LINE__, "init_heap : error in mutex_init (%s)\n", \
                                             strerror( status ) );      \
                                      EXIT; }

#  define LOCK( m )         pthread_mutex_lock(    & (m) )
#  define UNLOCK( m )       pthread_mutex_unlock(  & (m) )
#  define TRYLOCK( m )      (pthread_mutex_trylock( & (m) ) != EBUSY)
#  define MUTEX_BUSY        EBUSY
#  define MUTEX_LOCKED( m ) (TRYLOCK( m ) == MUTEX_BUSY)

/*
 * handling thread-specific-data
 */
#if USE_THR_HACKS == 1
#  define TSD_KEY_INIT( key )       { int i; (*(key)) = 0; for ( i = 0; i < TSD_MAX; i++ ) tsd[i] = NULL; }
#  define TSD_GET_DATA( key )       (tsd[ ((unsigned) pthread_self()) % 256 ])
#  define TSD_SET_DATA( key, data ) tsd[ ((unsigned) pthread_self()) % 256 ] = ((heap_t*) data);
#else
#  define TSD_KEY_INIT( key )       pthread_key_create( key, heap_key_destructor )
#  define TSD_GET_DATA( key )       pthread_getspecific( key )
#  define TSD_SET_DATA( key, data ) pthread_setspecific( key, (void *) data )
#endif

/*
 * conversion thread <-> heap
 */

#  define GET_THR_HEAP()        ((heap_t*) TSD_GET_DATA( heap_key ))
#  define SET_THR_HEAP( heap )  TSD_SET_DATA( heap_key, (void*) heap )

#else

/* maximal number of threads */
#define TSD_MAX   256

/*
 * handling of mutices
 */
#  define MUTEX_INITIALIZER  0
#  define DEFINE_MUTEX( m )  
#  define MUTEX_INIT(mutex,status)

#  define LOCK( m )     
#  define UNLOCK( m )   
#  define TRYLOCK( m )       1
#  define MUTEX_BUSY         0
#  define MUTEX_LOCKED( m )  false

/*
 * handling thread-specific-data
 */
#  define TSD_KEY_INIT( key )  

/*
 * conversion thread <-> heap
 */
#  define THREAD_ID             0
#  define GET_THR_HEAP()        & global_heap
#  define SET_THR_HEAP( heap )  

#endif  /* USE_PTHREAD == 1 */


/*****************************************************************
 *****************************************************************
 **
 ** heap specific data
 **
 *****************************************************************
 *****************************************************************/

/* lookup table for size classes */
#if ALIGNMENT == 8
static internal_size_t size_classes[ NO_OF_CLASSES ] = { 
    8, 16, 24, 32, 40, 48, 56, 64, 72, 80, 88, 96, 104, 112, 120, 128, 
    136, 144, 152, 160, 168, 176, 184, 192, 200, 208, 216, 224, 232, 240, 
    248, 256, 288, 320, 352, 384, 416, 448, 480, 512, 544, 576, 608, 640, 
    672, 704, 736, 768, 800, 832, 864, 896, 928, 960, 992, 1024, 1056, 
    1088, 1120, 1152, 1184, 1216, 1248, 1280, 1408, 1536, 1664, 1792, 
    1920, 2048, 2176, 2304, 2432, 2560, 2688, 2816, 2944, 3072, 3200, 
    3328, 3456, 3584, 3712, 3840, 3968, 4096, 4224, 4352, 4864, 5376, 
    5888, 6400, 6912, 7424, 7936, 8448, 8960, 9472, 9984, 10496, 11008, 
    11520, 12032, 12544, 13056, 13568, 14080, 14592, 15104, 15616, 16128, 
    16640, 18688, 20736, 22784, 24832, 26880, 28928, 30976, 33024, 35072, 
    37120, 39168, 41216, 43264, 45312, 47360, 49408, 57600, 65792, 73984, 
    82176, 90368, 98560, 106752, 114944, 123136, 131328, 139520, 147712, 
    155904, 164096, 172288, 180480, 213248, 246016, 278784, 311552, 344320, 
    377088, 409856, 442624, 475392, 508160, 540928, 573696, 606464, 639232, 
    672000, 704768, 775248, 852776, 938048, 1031856, 1135040, 1248544, 
    1373400, 1510736, 1661808, 1827992, 2010792, 2211864, 2433056, 2676360, 
    2943992, 3238392, 3562232, 3918456, 4310304, 4741328, 5215464, 5737008, 
    6310712, 6941784, 7635960, 8399552, 9239512, 10163456, 11179808, 12297784, 
    13527560, 14880320, 16368352, 18005184, 19805704, 21786272, 23964896, 
    26361392, 28997528, 31897280, 35087008, 38595704, 42455280, 46700808, 
    51370888, 56507976, 62158768, 68374648, 75212112, 82733320, 91006656, 
    100107320, 110118048, 121129856, 133242840, 146567120, 161223832, 
    177346216, 195080840, 214588920, 236047808, 259652592
};
#else
static internal_size_t size_classes[ NO_OF_CLASSES ] = { 
    16, 32, 48, 64, 80, 96, 112, 128, 144, 160, 176, 192, 208, 224, 240, 
    256, 272, 288, 304, 320, 336, 352, 368, 384, 400, 416, 432, 448, 464, 
    480, 496, 512, 576, 640, 704, 768, 832, 896, 960, 1024, 1088, 1152, 
    1216, 1280, 1344, 1408, 1472, 1536, 1600, 1664, 1728, 1792, 1856, 
    1920, 1984, 2048, 2112, 2176, 2240, 2304, 2368, 2432, 2496, 2560, 
    2688, 2816, 2944, 3072, 3200, 3328, 3456, 3584, 3712, 3840, 3968, 
    4096, 4224, 4352, 4480, 4608, 4736, 4864, 4992, 5120, 5248, 5376, 
    5504, 5632, 6144, 6656, 7168, 7680, 8192, 8704, 9216, 9728, 10240, 
    10752, 11264, 11776, 12288, 12800, 13312, 13824, 14336, 14848, 15360, 
    15872, 16384, 16896, 17408, 17920, 19968, 22016, 24064, 26112, 28160, 
    30208, 32256, 34304, 36352, 38400, 40448, 42496, 44544, 46592, 48640, 
    50688, 58880, 67072, 75264, 83456, 91648, 99840, 108032, 116224, 124416, 
    132608, 140800, 148992, 157184, 165376, 173568, 181760, 214528, 247296, 
    280064, 312832, 345600, 378368, 411136, 443904, 476672, 509440, 542208, 
    574976, 607744, 640512, 673280, 706048, 776656, 854320, 939760, 1033728, 
    1137104, 1250816, 1375888, 1513488, 1664832, 1831312, 2014448, 2215888, 
    2437472, 2681232, 2949344, 3244272, 3568704, 3925584, 4318128, 4749952, 
    5224944, 5747440, 6322176, 6954400, 7649824, 8414816, 9256288, 10181920, 
    11200112, 12320128, 13552128, 14907344, 16398080, 18037888, 19841680, 
    21825840, 24008432, 26409264, 29050192, 31955216, 35150736, 38665808, 
    42532384, 46785632, 51464192, 56610608, 62271664, 68498832, 75348720, 
    82883584, 91171936, 100289136, 110318048, 121349856, 133484832, 146833328, 
    161516656, 177668320, 195435152, 214978656, 236476528, 260124176
};
#endif
    
/* holds info about system */
static sysconf_t  sys_conf = { 0, DEFAULT_TOP_PAD };

/* global heap */
static heap_t global_heap = {
    0,
    NULL,
    { { NULL }, { NULL }, NULL, NULL, 0, { NULL }, NULL },
    { { 0 }, { 0 }, 0, 0, 0, 0 },
    DEFAULT_TOP_PAD,
    MUTEX_INITIALIZER,
    true
};

/* counter for the allocated heaps */
static unsigned     heap_id = 0;

/* linked list of heaps for threads (global heap is not a part of it) */
static heap_list_t  heap_list = { NULL, 0, MUTEX_INITIALIZER };

/* thread specific data */
#if USE_PTHREAD == 1 && USE_THR_HACKS == 1
static heap_t     * tsd[ TSD_MAX ];
#endif

/* is the heap initialised */
static bool         is_initialised = false;

/* how much statistics should be printed */
static unsigned     heap_stat_lvl = 0;

/* for memory trace */
static bool         trace_mem = false;
static trace_t      trace     = { NULL, 0, 0, 0, 0 };

/* key for accessing thread-specific-data */
#if USE_PTHREAD == 1
static tsd_key_t    heap_key;
#endif

/*****************************************************************
 *****************************************************************
 **
 ** defines for easier access
 **
 *****************************************************************
 *****************************************************************/


#define BLOCK_SIZE( b )     ((b)->size & ~0x7)
#define PRED_SIZE( b )      (*(((internal_size_t*) (b)) - 1))
#define SET_SIZE( b, n )    ((b)->size = (n) + ((b)->size & 0x7))
#define SET_EOB_SIZE( b )   { char * p = ((char*) b) + BLOCK_SIZE( b ) - sizeof(internal_size_t); \
                              *((internal_size_t*) p) = BLOCK_SIZE( b ); }

#define GET_HEAP( b )       ((heap_t *)      (((internal_size_t) (b)->hc_addr) & ~0x7))
#define GET_CONT( b )       ((container_t *) (((internal_size_t) (b)->hc_addr) & ~0x7))

#define IS_PRED_USED( b )   ((b)->size &  0x1)
#define SET_PRED_USED( b )  ((b)->size |= 0x1)
#define SET_PRED_FREE( b )  ((b)->size &= ~0x1)

#define IS_USED( b )        ((b)->size &  0x2)
#define SET_USED( b )       ((b)->size |= 0x2)
#define SET_FREE( b )       ((b)->size &= ~0x2)

#define IS_WILD( b )        ((b)->size &  0x4)
#define SET_WILD( b )       ((b)->size |= 0x4)
#define SET_TAME( b )       ((b)->size &= ~0x4)

#define IS_DATASEG( b )     (((internal_size_t) (b)->hc_addr) &  0x1)
#define SET_DATASEG( b )    ((b)->hc_addr = (void *) (((internal_size_t) (b)->hc_addr) | 0x1))

#define BLOCK_TO_PTR( n )   ((void *)   (((char*) n) + PREFIX_SIZE))
#define PTR_TO_BLOCK( p )   ((block_t*) (((char*) p) - PREFIX_SIZE))

#ifndef MAX   
#define MAX(x,y) ((x) > (y) ? (x) : (y))
#endif

#if DEBUG_CONT == 1
#define SET_CONT_MAGIC( c )       (c)->magic = CONT_MAGIC
#define CHECK_CONT_MAGIC( c, m )  if ( (c)->magic != CONT_MAGIC ) \
                                  { ERROR( __LINE__, "%s : container corrupted\n", m ); return; }
#else
#define SET_CONT_MAGIC( c )       
#define CHECK_CONT_MAGIC( c, m )
#endif

#if DEBUG_BLOCK == 1
#  define SET_BLOCK_MAGIC( b )         { (b)->magic = BLOCK_MAGIC; (b)->__pad__ = 0; }
#  define CHECK_BLOCK_MAGIC( b, m, r ) if ( (b)->magic != BLOCK_MAGIC ) \
                                       { ERROR( __LINE__, "%s : block corrupted\n", m ); r; }
#  define CHECK_BLOCK( heap, block )   check_block( (heap), (block) )
#else
#  define SET_BLOCK_MAGIC( b )      
#  define CHECK_BLOCK_MAGIC( b, m, r )
#  define CHECK_BLOCK( heap, block )   
#endif

#if DEBUG_RMALLOC == 1
#  define CHECK_RMALLOC      check_rmalloc()
#else
#  define CHECK_RMALLOC      
#endif

/*****************************************************************
 *****************************************************************
 **
 ** forward declarations and inline functions
 **
 *****************************************************************
 *****************************************************************/


/*
 * initialise/finish malloc
 */
static void rmalloc_init     ();
static void rmalloc_finish   ();

/*
 * heap management
 */
static void     init_heap     ( heap_t * heap );
static void     insert_heap   ( heap_t * heap );
static heap_t * get_free_heap ();
static void     insert_wild   ( heap_t * heap, block_t * block );

/*
 * allocation/deallocation
 */
void * r_malloc  ( size_t size );
void * r_calloc  ( size_t nmemb, size_t size);
void   r_free    ( void * ptr  );
void * r_realloc ( void * ptr, size_t size );

/********************************************
 * methods for small sizes
 ********************************************/

/*
 * allocation and deallocation
 */

static block_t * small_alloc ( heap_t  * heap, unsigned sclass );
static void      small_free  ( block_t * block );

/*
 * container management
 */

static void          remove_container ( container_t * container,
                                        unsigned      cfrom );
static void          insert_container ( container_t * container,
                                        unsigned      cto );
static container_t * alloc_container  ( heap_t      * heap,
                                        unsigned      sclass );
static void          init_container_blocks ( container_t * container );

/********************************************
 * methods for middle (usual) sizes
 ********************************************/

/*
 * allocation and deallocation
 */

static block_t * middle_alloc  ( heap_t * heap, internal_size_t size, unsigned sclass );
static void      middle_free   ( heap_t * heap, block_t * block );

/* insert free middle-sized block into given heap */
static void      middle_insert ( heap_t * heap, block_t * block );

/********************************************
 * methods for large (unmanaged) sizes
 ********************************************/

static block_t * vm_alloc   ( internal_size_t size );
static void      vm_dealloc ( block_t * block );

/*
 * trace memory allocation
 */
static void rmalloc_init_trace   ( const char * name );
static void rmalloc_trace        ();
static void rmalloc_finish_trace ();

/*
 * print statistics
 */
static void rmalloc_stat ();

/*
 * translate size into size-class
 */
static unsigned size_to_class  ( internal_size_t  size );

/*
 * check given block
 */
#if DEBUG_RMALLOC == 1
static int check_block ( heap_t * heap, block_t * block );
#endif

/*
 * check status of rmalloc
 */
#if DEBUG_RMALLOC == 1
static void check_rmalloc ();
#endif

/*
 * destructor for thread-private data
 */
#if USE_PTHREAD == 1 && USE_THR_HACKS == 0
static void heap_key_destructor ( void * arg  );
#endif

/*
 * message output
 */
static void OUTPUT ( const char * msg, ... );
static void ERROR  ( int lineno, const char * msg, ... );

/*
 * local debug hook called on error detection
 */
void __RMALLOC_debug ();

/*****************************************************************
 *****************************************************************
 **
 ** allocation and deallocation of memory from system
 **
 *****************************************************************
 *****************************************************************/

DEFINE_MUTEX( system_mutex )
    
/*
 * allocation from system via sbrk / mmap etc.
 */
static void *
system_alloc ( heap_t * heap, internal_size_t size )
{
    char            * p;
    block_t         * block, * remainder;
    internal_size_t   n;
    
    if ( size == 0 )
        return NULL;

    assert( heap != NULL );

    /* adjust n to be multiple of ALIGNMENT */
    if ( IS_MISALIGNED( size ) )
    {
        ERROR( __LINE__, "system_alloc : warning, size is not a multiple of %d\n", ALIGNMENT );
        ALIGN( size );
    }

    /*
     * use wilderness of heap to serve malloc request
     */

    block = heap->blocks.wilderness;

    while ( block != NULL )
    {
        if ( BLOCK_SIZE( block ) >= size )
        {
            /* remove block from wilderness */
            if ( block->next != NULL ) block->next->prev = block->prev;
            if ( block->prev != NULL ) block->prev->next = block->next;
            else                       heap->blocks.wilderness = block->next;
                
            /* if block is large enough, it is split and the remainder       *
             * is moved into standard middle sized data managment, otherwise *
             * the block is not split and the remainder is regarded as       *
             * internal fragmentation                                        */

            n = BLOCK_SIZE( block ) - size;

            if ( n >= MIN_BLOCK_SIZE )
            {
                /* split block */
                remainder = (block_t*) (((char*) block) + size);

                SET_SIZE( block, size );
                SET_TAME( block );
                SET_BLOCK_MAGIC( block );

                /* overwrite data in size and hc_addr fields !!! */
                remainder->size = n;
                SET_WILD( remainder );
                SET_FREE( remainder );
                SET_BLOCK_MAGIC( remainder );
                remainder->hc_addr = GET_HEAP( block );
                remainder->prev    = NULL;
                remainder->next    = NULL;

                middle_insert( heap, remainder );
            }
            
            return block;
        }
        else
            block = block->next;
    }

    /*
     * no suitable block was found in the wilderness, create a
     * new block (plus padding) and serve the request
     */

#if USE_MMAP == 1 || USE_MALLOC == 1
    /*
     * use system-malloc to request a large memory-chunk
     */
        
    /* adjust top-pad of heap */
    if ( heap->stat.used_mem / TOP_PAD_FRACTION > heap->top_pad )
        heap->top_pad = heap->stat.used_mem / TOP_PAD_FRACTION;

    n = size + heap->top_pad;
    
    /* adjust size to be multiple of pagesize */
    if ( n % sys_conf.pagesize != 0 )
        n += sys_conf.pagesize - (n % sys_conf.pagesize);

#if USE_MALLOC == 1
    /* call malloc in a thread-safe way */
    LOCK( system_mutex );
    p = (char*) malloc( n );
    UNLOCK( system_mutex );
    
    if ( p == (char*) NULL )
        return NULL;
#else /* USE_MMAP == 1 */
    /* map new heap with mmap */
    p = (char*) MMAP( 0, n, PROT_READ | PROT_WRITE, MAP_PRIVATE );
    
    if ( p == (char*) MAP_FAILED )
        return NULL;
#endif

#elif USE_SBRK == 1
    /* to be done */
#else
    ERROR( __LINE__, "system_alloc : no heap-allocation method defined (MALLOC,MMAP,SBRK)\n" );
    EXIT;
#endif

    /*
     * update statistics
     */

    heap->stat.used_mem += n;
    heap->stat.max_mem   = MAX( heap->stat.used_mem, heap->stat.max_mem );

    /*
     * form block and insert remainder into middle sized management
     */

    block          = (block_t *) p;
    block->size    = size;
    block->hc_addr = heap;
    SET_BLOCK_MAGIC( block );

    /* this is the first block in the data-segment and therefore has no
     * predecessor. by setting pred_used, we avoid coalescing with it.
     * it is also marked as a data-segment block */
    SET_PRED_USED( block );
    SET_DATASEG( block );

    remainder          = (block_t*) (p + size);
    remainder->size    = n - size;
    remainder->hc_addr = heap;
    remainder->prev    = NULL;
    remainder->next    = NULL;
    SET_WILD( remainder );
    SET_BLOCK_MAGIC( remainder );

    middle_insert( heap, remainder );
    
    return block;
}

/*
 * allocation from system via sbrk / mmap etc.
 */
static void
system_dealloc ( heap_t * heap, block_t * block )
{
    assert(( heap != NULL ) && ( block != NULL ));

    /*
     * free memory according the previously used method
     */

    internal_size_t  n = BLOCK_SIZE( block );

#if USE_MMAP == 1 || USE_MALLOC == 1
    /*
     * use system-malloc to request a large memory-chunk
     */
        
#if USE_MALLOC == 1
    
    /* call free in a thread-safe way */
    LOCK( system_mutex );
    free( block );
    UNLOCK( system_mutex );
    
#else /* USE_MMAP == 1 */
    
    /* unmap block (munmap does not seem to be multi-thread safe) */
    LOCK( system_mutex );
    munmap( (char *) block, n );
    UNLOCK( system_mutex );
    
#endif

#elif USE_SBRK == 1
    /* to be done */
#else
    ERROR( __LINE__, "system_dealloc : no heap-deallocation method defined (MALLOC,MMAP,SBRK)\n" );
    EXIT;
#endif

    /*
     * update statistics
     */

    heap->stat.used_mem -= n;
}



/********************************************************************
 ********************************************************************
 **
 ** malloc/free interface
 **
 ********************************************************************
 ********************************************************************/

void *
r_malloc ( size_t size )
{
    block_t  * block  = NULL;
    heap_t   * heap   = NULL;
    unsigned   sclass = 0;
    

    if ( size == 0 )
        return NULL;

    if ( ! is_initialised )
        rmalloc_init();

    CHECK_RMALLOC;

    /* add management data for all sizes */
    size = size + PREFIX_SIZE;
        
    /* first check if size exceeds maximal managed size */
    if ( size > size_classes[ NO_OF_CLASSES - 1 ] )
    {
        block_t  * lblock;
        
        /* request is managed by the operating system */
        lblock = vm_alloc( size );

        if ( lblock == NULL )
            return NULL;
        
        /* block has no heap/container */
        lblock->hc_addr = NULL;

        /* update size because of possible pagesize alignment */
        size = BLOCK_SIZE( lblock );
        
        /* update statistics (of global heap) */
        LOCK( global_heap.mutex );
        global_heap.stat.used_mem   += size;
        global_heap.stat.max_mem     = MAX( heap->stat.max_mem, heap->stat.used_mem );
        global_heap.stat.mem_in_use += size;
        global_heap.stat.max_in_use  = MAX( heap->stat.max_in_use, heap->stat.mem_in_use );
        global_heap.stat.allocated  += size;
        UNLOCK( global_heap.mutex );
        
        return BLOCK_TO_PTR( lblock );
    }
    else
    {
        /* depending on size, adjust to hold extra data */
        if ( size <= SMALL_SIZE )
            size = MAX( size, sizeof(sblock_t) );
        else
            size = MAX( size, sizeof(block_t) + POSTFIX_SIZE );
        
        /* trim size to be multiple of ALIGNMENT */
        if ( IS_MISALIGNED( size ) )
            ALIGN( size );

        sclass = size_to_class( size );
    }
    
    /*
     * get heap for this thread
     */
    
    /* get private heap of thread or set new heap */
    heap = GET_THR_HEAP();
#if ONE_HEAP_PER_THREAD == 1
    if ( heap == NULL )
#else
    if (( heap == NULL ) || ( TRYLOCK( heap->mutex ) == false ))
#endif
    {
        /* look for free heap or create new one */
        if ((heap = get_free_heap()) == NULL)
        {
            internal_size_t  soh = EXTRA_DATA_SIZE + sizeof(heap_t);
            
            if ( IS_MISALIGNED( soh ) )
                ALIGN( soh );
            
            LOCK( global_heap.mutex );
            heap = (heap_t *) BLOCK_TO_PTR( (block_t *) middle_alloc( (heap_t *) & global_heap,
                                                                      soh, size_to_class( soh ) ) );
            UNLOCK( global_heap.mutex );

            if ( heap == NULL )
                return NULL;
            
            init_heap( heap );

            /* lock heap BEFORE inserting so we can be sure, no other grabs it */
            LOCK( heap->mutex );
            insert_heap( heap );
        }
#if ONE_HEAP_PER_THREAD == 1
        else
            LOCK( heap->mutex );
#endif

        /* set this heap to be the thread-private heap */
        SET_THR_HEAP( heap );
    }
#if ONE_HEAP_PER_THREAD == 1
    else
        LOCK( heap->mutex );
#endif

    /*
     * get free block, depending on size
     */

    if ( size <= SMALL_SIZE )
        block = small_alloc( heap, sclass );
    else
    {
        block = middle_alloc( heap, size, sclass );
        CHECK_BLOCK( heap, block );
    }
    
    if ( block == NULL )
    {
        UNLOCK( heap->mutex );
        return NULL;
    }

    /* update statistics */
    heap->stat.mem_in_use += BLOCK_SIZE( block );
    heap->stat.max_in_use  = MAX( heap->stat.max_in_use, heap->stat.mem_in_use );
    heap->stat.allocated  += BLOCK_SIZE( block );
    
    UNLOCK( heap->mutex );

    if ( trace_mem )
        rmalloc_trace();

    CHECK_RMALLOC;

    return BLOCK_TO_PTR( block );
}

void *
r_calloc ( size_t nmemb, size_t size )
{
    void * p;
    
    if ((nmemb == 0) || (size == 0))
        return NULL;
    
    p = r_malloc( nmemb * size );

    /* memset to zero */
    if ( p != NULL )
        memset( p, 0, nmemb * size );

    return p;
}

void
r_free ( void * ptr )
{
    block_t  * block;

    if ( ptr == NULL )
        return;

    /* should at least have called malloc for this ptr, if not : ERROR */
    if ( ! is_initialised )
    {
        ERROR( __LINE__, "free : rmalloc not initialised\n" );
        return;
    }

    CHECK_RMALLOC;

    block = PTR_TO_BLOCK( ptr );

    CHECK_BLOCK_MAGIC( block, "free", return );
    
    /* test whether block is managed by operating system */
    if ( block->hc_addr == NULL )
    {
        /* update statistics (of global heap) */
        LOCK( global_heap.mutex );
        global_heap.stat.used_mem   -= BLOCK_SIZE( block );
        global_heap.stat.mem_in_use -= BLOCK_SIZE( block );
        UNLOCK( global_heap.mutex );
        
        vm_dealloc( block );
    }
    else
    {
        /* now free block, depending on size */
        if ( BLOCK_SIZE( block ) <= SMALL_SIZE )
            small_free( block );
        else
        {
            heap_t * heap = GET_HEAP( block );

            assert( heap != NULL );

            /* lock mutex here, cause middle_free is also called
             * from small_free and mutex is already locked there
             */
            
            CHECK_BLOCK( heap, block );
            
            LOCK( heap->mutex );

            /* update statistics */
            heap->stat.mem_in_use -= BLOCK_SIZE( block );
    
            middle_free( heap, block );

            UNLOCK( heap->mutex );
        }
    }

    if ( trace_mem )
        rmalloc_trace();

    CHECK_RMALLOC;
}

void *
r_realloc ( void * ptr, size_t size )
{
    if ( ptr == NULL )
        return r_malloc( size );
    else
    {
        block_t        * block    = PTR_TO_BLOCK( ptr );
        void           * newptr   = NULL;
        internal_size_t  new_size = size + PREFIX_SIZE;
        internal_size_t  old_size = BLOCK_SIZE( block );

        CHECK_BLOCK_MAGIC( block, "realloc", return NULL );

        /*
         * adjust new size to get correct size of block
         * corresponding to this new size
         */

        if ( new_size <= SMALL_SIZE )
            new_size = MAX( new_size, sizeof(sblock_t) );
        else if ( new_size <= size_classes[ NO_OF_CLASSES - 1 ] )
            new_size  = MAX( new_size, sizeof(block_t) + POSTFIX_SIZE );

        if ( IS_MISALIGNED( size ) )
            ALIGN( size );

        /* if block has same size, just return */
        if ( new_size == old_size )
            return ptr;

        /* create new block and copy data */
        newptr = r_malloc( size );

        if ( newptr == NULL )
        {
            r_free( ptr );
            return NULL;
        }
        
        old_size -= PREFIX_SIZE;

        if ( size > old_size )
            memcpy( newptr, ptr, old_size );
        else
            memcpy( newptr, ptr, size );

        r_free( ptr );

        return newptr;
    }
}

/********************************************************************
 ********************************************************************
 **
 ** misc. functions for heap-management
 **
 ********************************************************************
 ********************************************************************/

/*
 * init/finish heap-manager
 */

static void
rmalloc_init ()
{
    DEFINE_MUTEX( mutex )

    if ( MUTEX_LOCKED( mutex ) )
    {
        ERROR( __LINE__, "rmalloc_init : mutex locked\n" );
        LOCK( mutex );
    }

    if ( ! is_initialised )
    {
        char * value;

#if HAS_GETPAGESIZE == 1
        sys_conf.pagesize = getpagesize();
#else
        sys_conf.pagesize = 8*1024;  /* some reasonable default value */
#endif

        /* init tsd-key */
        TSD_KEY_INIT( & heap_key );
        
        /* re-initialise global heap */
        init_heap( & global_heap );

        /*
         * setup memory trace
         */

        value = getenv( "RMALLOC_TRACE" );
        
        if (( value != NULL ) && ((value[0] == '1') || (value[0] == '2')))
        {
            /* set trace-type */
            if ( value[0] == '1' )
                trace.type = TRACE_STEPS;
            else
                trace.type = TRACE_ALLOCATION;
            
            /* get size of trace */

            if ((value = getenv( "RMALLOC_TRACE_SIZE" )) != NULL )
            {
                if ( strcmp( value, "all" ) == 0 )
                    trace.size = 0;
                else
                {
                    trace.size = atoi( value );

                    /* round to next multiple of ALIGNMENT */
                    trace.size += ALIGNMENT - (trace.size % ALIGNMENT);
                }
            }

            /* get name of tracefile and initialise */

            if ((value = getenv( "RMALLOC_TRACE_FILE" )) != NULL )
                rmalloc_init_trace( value );
            else
                rmalloc_init_trace( DEFAULT_TRACE_FILE );
        }
        
        /*
         * register cleanup-functions
         */

        if ( atexit( rmalloc_finish ) != 0 )
            ERROR( __LINE__, "init : error in atexit (%s)\n", strerror( errno ) );
        
        if ( atexit( rmalloc_finish_trace ) != 0 )
            ERROR( __LINE__, "init : error in atexit (%s)\n", strerror( errno ) );

        if ((value = getenv( "RMALLOC_STAT" )) != NULL )
        {
            heap_stat_lvl = atoi( value );

            if ( heap_stat_lvl > 0 )
            {
                if ( atexit( rmalloc_stat ) != 0 )
                    ERROR( __LINE__, "init : error in atexit (%s)\n",
                           strerror( errno ) );
            }
        }

        is_initialised = true;
    }
    
    UNLOCK( mutex );
}

/*
 * end it all
 */
static void
rmalloc_finish ()
{
    /*
     * clean will be called at the end and then
     * the memory will be deallocated by the system
     * so why bother ???
     */
}

/*
 * initialise a heap
 */
static void
init_heap ( heap_t * heap )
{
    unsigned  i;
#if USE_PTHREAD == 1
    int       status;
#endif
    
    heap->next    = NULL;
    heap->id      = heap_id++;
    heap->top_pad = sys_conf.top_pad;

    for ( i = 0; i < NO_OF_SMALL; i++ )
        heap->blocks.nonempty[i] = NULL;

    heap->blocks.full       = NULL;
    heap->blocks.last_full  = NULL;
    heap->blocks.no_of_full = 0;
    
    for ( i = 0; i < NO_OF_CLASSES; i++ )
    {
        heap->blocks.blocks[i] = NULL;
        
        heap->stat.free[i]  = 0;
        heap->stat.used[i]  = 0;
    }

    heap->blocks.wilderness = NULL;

    heap->stat.used_mem     = 0;
    heap->stat.max_mem      = 0;
    heap->stat.mem_in_use   = 0;
    heap->stat.max_in_use   = 0;
    heap->stat.allocated    = 0;

    heap->used = true;

    if ( heap != & global_heap )
    {
        MUTEX_INIT( heap->mutex, status );
    }
}

/*
 * inserts new heap into heap-list
 */
static void
insert_heap ( heap_t * heap )
{
    LOCK( heap_list.mutex );
    
    /* prepend into list */
    heap->next      = heap_list.heaps;
    heap_list.heaps = heap;

    heap_list.nheaps++;
    
    UNLOCK( heap_list.mutex );
}

/*
 * return unused/free heap from heap-list
 */
static heap_t *
get_free_heap ()
{
#if ONE_HEAP_PER_THREAD == 1
    heap_t  * heap;
    
    LOCK( heap_list.mutex );
    
    heap = heap_list.heaps;

    while ( heap != NULL )
    {
        if ( ! heap->used )
        {
            heap->used = true;
            UNLOCK( heap_list.mutex );
            return heap;
        }
        else
            heap = heap->next;
    }

    UNLOCK( heap_list.mutex );
    
    /* no free heap found */
    return NULL;
#else
    heap_t  * thr_heap = GET_THR_HEAP();
    heap_t  * heap     = NULL;

    if (( thr_heap != NULL ) && TRYLOCK( thr_heap->mutex ) )
        return thr_heap;
    
    LOCK( heap_list.mutex );
    
    heap = heap_list.heaps;

    while ( heap != NULL )
    {
        if (( heap != thr_heap ) && TRYLOCK( heap->mutex ) )
        {
            heap->used = true;
            UNLOCK( heap_list.mutex );
            return heap;
        }

        heap = heap->next;
    }

    UNLOCK( heap_list.mutex );
    
    /* no free heap found */
    return NULL;
#endif
}

/*
 * insert given block into wilderness of heap
 */
static void
insert_wild ( heap_t * heap, block_t * block )
{
    /*
     * before we insert the block, we compute the amount of memory
     * in the wilderness of the heap, if this is large enough, the
     * block is given back to the operating system
     */

    internal_size_t  mem  = 0;
    block_t        * wild = heap->blocks.wilderness;
            
    while ( wild != NULL )
    {
        mem  += BLOCK_SIZE( wild );
        wild  = wild->next;
    }

    /* if enough memory is available, free the block */
    if ( mem > 2 * BLOCK_SIZE( block ) )
    {
        system_dealloc( heap, block );
        return;
    }
    
#if SORT_WILDERNESS == 1
    /*
     * the block is put into the wilderness according to his size
     */
    
    block_t * next = NULL;
    block_t * prev = NULL;

    next = heap->blocks.wilderness;
    while ( next != NULL )
    {
        if ( BLOCK_SIZE( next ) < BLOCK_SIZE( block ) )
        {
            prev = next;
            next = next->next;
        }
        else
            break;
    }
    
    /* insert before nextessor */
    block->next = next;
    block->prev = prev;
    
    if ( next != NULL ) next->prev = block;
    if ( prev != NULL ) prev->next = block;
    else                heap->blocks.wilderness = block;
#else
    /*
     * simply put the block in front of the wilderness
     */

    block->next = heap->blocks.wilderness;
    if ( block->next != NULL )
        block->next->prev = block;
    block->prev = NULL;
    heap->blocks.wilderness = block;
#endif
}

/*
 * translate size into size-class
 */
static unsigned
size_to_class ( internal_size_t  size )
{
    if ( size <= SMALL_SIZE )
        return ((unsigned) ((size / ALIGNMENT) - 1));
    else
    {
#if USE_BIN_SEARCH == 1
        /*
         * use binary search for size-class
         */
        unsigned  lb = (SMALL_SIZE / ALIGNMENT)-1;
        unsigned  ub = NO_OF_CLASSES;
        unsigned  split;
        
        while ( ub - lb > 1 )
        {
            split = (ub + lb) / 2;
            
            if (( size_classes[ split-1 ] < size ) && (size <= size_classes[ split ]))
                return split;
            
            if ( size_classes[ split ] < size )
                lb = split;
            else
                ub = split;
        }
        
        return ub;
#else
        /*
         * use linear search for size-class
         */
        unsigned  sclass = SMALL_SIZE / ALIGNMENT;

        while ((sclass < NO_OF_CLASSES) && (size_classes[ sclass ] < size))
            sclass++;

        return sclass;
#endif
    }
}

/*
 * called when a thread is finished and releases
 * local heap into list of unused heaps
 */
#if USE_PTHREAD == 1 && USE_THR_HACKS == 0
static void
heap_key_destructor ( void * arg )
{
    heap_t  * heap = (heap_t *) arg;

    if ( heap != NULL )
        heap->used = false;
}
#endif

/********************************************************************
 ********************************************************************
 **
 ** allocation and deallocation of small sized blocks
 **
 ********************************************************************
 ********************************************************************/

/*
 * malloc/free for small sizes
 */
static block_t *
small_alloc ( heap_t * heap, unsigned sclass )
{
    container_t  * container;
    sblock_t     * block;

    assert( heap != NULL );

    /*
     * try to use block from non-empty container, if no such container
     * exists, try full ones or transfer a container from the global
     * heap. Only of no container is available at all, allocate a new one.
     */
    
    container = heap->blocks.nonempty[ sclass ];

    if ( container == NULL )
    {
        container = alloc_container( heap, sclass );

        if ( container == NULL )
            return NULL;
        
        insert_container( container, TYPE_NONEMPTY );
    }

    block = container->blocks;
    container->blocks = block->next;

    container->free--;
    
    heap->stat.free[ container->sclass ]--;
    heap->stat.used[ container->sclass ]++;
    
    if ( container->free == 0 )
        remove_container( container, TYPE_NONEMPTY );

    return (block_t*) block;
}

static void
small_free ( block_t * block )
{
    container_t * container;
    heap_t      * heap;
    sblock_t    * sblock = (sblock_t*) block;

    assert( sblock != NULL );
    
    container = GET_CONT( sblock );
    CHECK_CONT_MAGIC( container, "small_free" );
    
    heap = container->heap;
    
    LOCK( heap->mutex );
    
    /* update statistics */
    heap->stat.mem_in_use -= BLOCK_SIZE( sblock );
    
    /*
     * Insert block back into container, if container is full,
     * move it to the "full"-list and transfer last full
     * container to global heap if too many full containers exist
     * (see insert_container)
     */

    sblock->next      = container->blocks;
    container->blocks = sblock;
    container->free++;
    
    heap->stat.free[ container->sclass ]++;
    heap->stat.used[ container->sclass ]--;

    /* if free == 1, container was empty before */
    if ( container->free == 1 )
        insert_container( container, TYPE_NONEMPTY );

    if ( container->free == container->max )
    {
        remove_container( container, TYPE_NONEMPTY );
        insert_container( container, TYPE_FULL );
    }

    UNLOCK( heap->mutex );
}

/********************************************************************
 ********************************************************************
 **
 ** container managment
 **
 ********************************************************************
 ********************************************************************/

/*
 * remove given container from class "cfrom"
 */
static void
remove_container ( container_t * container, unsigned cfrom )
{
    heap_t   * heap;
    unsigned   sclass;

    assert( container != NULL );

    heap   = container->heap;
    sclass = container->sclass;

    /*
     * release container from old class
     */
    
    if ( container->prev != NULL )
        container->prev->next = container->next;
    
    if ( container->next != NULL )
        container->next->prev = container->prev;
            
    switch ( cfrom )
    {
    case TYPE_NONEMPTY :
        if ( heap->blocks.nonempty[ sclass ] == container )
            heap->blocks.nonempty[ sclass ] = container->next;
        break;
        
    case TYPE_FULL :
        if ( heap->blocks.full == container )
            heap->blocks.full = container->next;

        if ( heap->blocks.last_full == container )
            heap->blocks.last_full = container->prev;
        
        heap->blocks.no_of_full--;
        break;
        
    default:
        break;
    }
}

/*
 * insert (prepend) given container into class "cto"
 */
static void
insert_container ( container_t * container, unsigned cto )
{
    heap_t   * heap;
    unsigned   sclass;

    assert( container != NULL );

    heap   = container->heap;
    sclass = container->sclass;

    /*
     * prepend container to new class
     */
    
    container->prev = NULL;
    
    switch ( cto )
    {
    case TYPE_NONEMPTY :
        container->next = heap->blocks.nonempty[ sclass ];
        if ( heap->blocks.nonempty[ sclass ] != NULL )
            heap->blocks.nonempty[ sclass ]->prev = container;
        heap->blocks.nonempty[ sclass ] = container;
        break;
        
    case TYPE_FULL     :
        container->next = heap->blocks.full;
        if ( heap->blocks.full != NULL )
            heap->blocks.full->prev = container;
        else
            heap->blocks.last_full = container;
        heap->blocks.full = container;
        heap->blocks.no_of_full++;
        break;
        
    default:
        break;
    }

    /*
     * check if there are too many full containers in heap
     * and release last container to middle-malloc
     */
    
    if (( cto == TYPE_FULL ) && ( heap->blocks.no_of_full > MAX_FULL_CONTAINERS ))
    {
        container_t * last = heap->blocks.last_full;

        heap->blocks.last_full = last->prev;
        last->prev->next = NULL;
        heap->blocks.no_of_full--;

        middle_free( heap, PTR_TO_BLOCK( last ) );
    }
}

/*
 * return container
 */
static container_t *
alloc_container ( heap_t * heap, unsigned sclass )
{
    container_t * container;

    assert( heap != NULL );

    if ( heap->blocks.full != NULL )
    {
        /*
         * use old unused containers
         */
    
        container = heap->blocks.full;

        if ( container->next != NULL )
            container->next->prev = NULL;
        else
            heap->blocks.last_full = NULL;
        
        heap->blocks.full = container->next;
        heap->blocks.no_of_full--;

        container->next = NULL;
        container->prev = NULL;
        
        if ( container->sclass != sclass )
        {
            container->sclass = sclass;
            init_container_blocks( container );
        }
    }
    else
    {
        /*
         * allocate container from middle-malloc
         */

        const internal_size_t  csize = EXTRA_DATA_SIZE + sizeof(container_t) + CONTAINER_SIZE;
        void                 * p;

        p = middle_alloc( heap, csize, size_to_class( csize ) );

        if ( p == NULL )
            return NULL;
        
        container = (container_t *) BLOCK_TO_PTR( p );
        
        container->heap   = heap;
        container->next   = NULL;
        container->prev   = NULL;
        container->sclass = sclass;

        SET_CONT_MAGIC( container );
        
        /* build blocks in container */
        init_container_blocks( container );
    }
    
    return container;
}

static void
init_container_blocks ( container_t * container )
{
    internal_size_t  block_size;
    unsigned         count, i;
    char           * addr;
    block_t        * old_block = NULL;

    assert( container != NULL );

    block_size = size_classes[ container->sclass ];
    count      = (unsigned) (CONTAINER_SIZE / block_size);
    addr       = ((char*) container) + sizeof(container_t);

    container->max    = count;
    container->free   = count;
    container->blocks = (sblock_t*) addr;
    
    for ( i = 0; i < count; i++ )
    {
        block_t * block = (block_t*) addr;

        block->size    = block_size;
        block->hc_addr = container;
        if ( old_block != NULL )
            old_block->next = block;
        old_block = block;

        SET_BLOCK_MAGIC( block );
        
        addr += block_size;
    }

    old_block->next = NULL;

    container->heap->stat.free[ container->sclass ] += container->max;
}

/********************************************************************
 ********************************************************************
 **
 ** allocation and deallocation of middle-sized blocks
 **
 ********************************************************************
 ********************************************************************/

static block_t *
middle_alloc ( heap_t * heap, internal_size_t size, unsigned sclass )
{
    block_t  * block = NULL;
    unsigned   i;
    
    assert( heap != NULL );

    /*
     * look for an old block, and begin search with
     * given size-class
     */

    for ( i = sclass; i < NO_OF_CLASSES; i++ )
    {
        block = heap->blocks.blocks[ i ];
        
        while ( block != NULL )
        {
            if ( BLOCK_SIZE( block ) >= size )
            {
                /* remove block from list */
                if ( block->next != NULL ) block->next->prev = block->prev;
                if ( block->prev != NULL ) block->prev->next = block->next;
                else                       heap->blocks.blocks[ i ] = block->next;
            
                block->next = NULL;
                block->prev = NULL;

                /* update statistics */
                heap->stat.free[ i ]--;
            
                /*
                 * check if block is from a different size-class
                 * and if the remainder of the block can be used
                 * by another sclass
                 */
                
                if (( i != sclass ) && ( BLOCK_SIZE( block ) - size > MIN_BLOCK_SIZE))
                {
                    /*
                     * split block and put remainder into another size-class
                     */

                    block_t * remainder = (block_t*) (((char*) block) + size);

                    /* overwrite stored data in size field */
                    remainder->size = BLOCK_SIZE( block ) - size;
                    SET_FREE( remainder );
                    SET_BLOCK_MAGIC( remainder );

                    if ( IS_WILD( block ) )
                        SET_WILD( remainder );

                    remainder->hc_addr = heap;
                    
                    /* put block into proper size-class */
                    middle_insert( heap, remainder );

                    /* adjust blocksize of old block */
                    SET_SIZE( block, size );
                    SET_TAME( block );
                }
                else
                    sclass = i;

                /* finish search */
                i = NO_OF_CLASSES;
                
                break;
            }
            else
                block = block->next;
        }
    }

    /*
     * allocate block from system
     */

    if ( block == NULL )
    {
        block = (block_t*) system_alloc( heap, size );

        if ( block == NULL )
            return NULL;
    }

    /* update statistics */
    heap->stat.used[ sclass ]++;

    SET_USED( block );

    /* set boundary tag in successor */
    if ( ! IS_WILD( block ) )
    {
        block_t * succ = (block_t*) (((char*) block) + BLOCK_SIZE( block ));

        SET_PRED_USED( succ );
    }

    return block;
}
               
static void
middle_free ( heap_t * heap, block_t * block )
{
    unsigned  sclass;

    assert((heap != NULL) && ( block != NULL ));

    sclass = size_to_class( BLOCK_SIZE( block ) );
    heap->stat.used[ sclass ]--;

    /*
     * if block is not wild, we try to coalesce it with succ
     */

    if ( ! IS_WILD( block ) )
    {
        block_t * succ = (block_t*) (((char*) block) + BLOCK_SIZE( block ));

        if ( ! IS_USED( succ ) )
        {
            /* remove successor from his list */
            sclass = size_to_class( BLOCK_SIZE( succ ) );
            
            if ( succ->next != NULL ) succ->next->prev = succ->prev;
            if ( succ->prev != NULL ) succ->prev->next = succ->next;
            else                      heap->blocks.blocks[ sclass ] = succ->next;
            
            heap->stat.free[ sclass ]--;
            
            /* if succ is wild, coalesced block is also wild */
            if ( IS_WILD( succ ) )
                SET_WILD( block );
            
            /* coalesce, e.g. adjust size of block */
            SET_SIZE( block, BLOCK_SIZE( block ) + BLOCK_SIZE( succ ) );
        }
    }

    /*
     * try to coalesce with predecessor
     */

    if ( ! IS_PRED_USED( block ) )
    {
        internal_size_t  psize = PRED_SIZE( block );
        block_t        * pred  = (block_t*) (((char*) block) - psize);
        
        /* remove predecessor from his list */
        sclass = size_to_class( psize );
        
        if ( pred->next != NULL ) pred->next->prev = pred->prev;
        if ( pred->prev != NULL ) pred->prev->next = pred->next;
        else                      heap->blocks.blocks[ sclass ] = pred->next;

        heap->stat.free[ sclass ]--;
        
        /* if block is wild, so is the coalesced block */
        if ( IS_WILD( block ) )
            SET_WILD( pred );

        /* coalesce, e.g. adjust size */
        SET_SIZE( pred, BLOCK_SIZE( block ) + psize );

        /* replace block */
        block = pred;
    }

    /*
     * put empty data-segment into wilderness or standard block
     * corresponding size-class
     */

    if ( IS_WILD( block ) && IS_DATASEG( block ) )
        insert_wild( heap, block );
    else
        middle_insert( heap, block );
}

/*
 * insert free middle-sized block into given heap
 */
static void
middle_insert ( heap_t * heap, block_t * block )
{
    assert(( heap != NULL ) && ( block != NULL ));

#if SORT_MIDDLE_BLOCKS == 1
    
    block_t  * next;
    block_t  * prev   = NULL;
    unsigned   sclass = size_to_class( BLOCK_SIZE( block ) );

    /* sort block into list, according to size */
    next = heap->blocks.blocks[ sclass ];
        
    while ( next != NULL )
    {
        if ( BLOCK_SIZE( next ) < BLOCK_SIZE( block ) )
        {
            prev = next;
            next = next->next;
        }
        else
            break;
    }

    /* insert before successor */
    block->next = next;
    block->prev = prev;

    if ( next != NULL ) next->prev = block;
    if ( prev != NULL ) prev->next = block;
    else                heap->blocks.blocks[ sclass ] = block;
    
#else
    
    unsigned sclass = size_to_class( BLOCK_SIZE( block ) );

    /* put block at begining of list */
    block->next = heap->blocks.blocks[ sclass ];

    if ( block->next != NULL )
        block->next->prev = block;

    block->prev = NULL;
    heap->blocks.blocks[ sclass ] = block;
    
#endif
        
    heap->stat.free[ sclass ]++;

    /* copy size-information to end of block */
    SET_EOB_SIZE( block );

    /* adjust boundary tag */
    if ( ! IS_WILD( block ) )
    {
        block_t * succ = (block_t*) (((char*) block) + BLOCK_SIZE( block ));
        SET_PRED_FREE( succ );
    }

    SET_FREE( block );
}

/********************************************************************
 ********************************************************************
 **
 ** handling of unmanaged chunks
 **
 ********************************************************************
 ********************************************************************/

/*
 * allocate block directly
 */
static block_t *
vm_alloc ( internal_size_t size )
{
    block_t  * block = NULL;
    
#if USE_MMAP == 1
    char     * p     = NULL;
#endif
    
    /* round size up to next multiple of pagesize */
    if ( size % sys_conf.pagesize != 0 )
        size += sys_conf.pagesize - (size % sys_conf.pagesize);
    
#if USE_MALLOC == 1
    block = (block_t*) malloc( size );
#elif USE_MMAP == 1
    p     = (char*) MMAP( 0, size, PROT_READ | PROT_WRITE, MAP_PRIVATE );

    if ( p == (char*) MAP_FAILED )
        return NULL;

    block = (block_t *) p;
#else
    ERROR( __LINE__, "vm_alloc : no supported allocation method specified\n" );
    return NULL;
#endif

    if ( block != NULL )
        block->size = size;
    
    return block;
}

static void
vm_dealloc ( block_t * block )
{
#if USE_MALLOC == 1
    return free( block );
#else
#if USE_MMAP == 1
    if ( BLOCK_SIZE( block ) <= size_classes[ NO_OF_CLASSES - 1 ] )
    {
        ERROR( __LINE__, "vm_dealloc : size of block too small for vm_dealloc\n" );
        return;
    }
    
    if ( munmap( (char*) block, BLOCK_SIZE( block ) ) != 0 )
    {
        ERROR( __LINE__, "vm_dealloc : munmap failed (%s)\n",
               strerror( errno ) );
        return;
    }
#endif
#endif
}

/********************************************************************
 ********************************************************************
 **
 ** get info and set options from/for malloc
 **
 ********************************************************************
 ********************************************************************/

/*
 * set malloc options
 */

int
r_mallopt ( int param, int val )
{
    if ( param == M_TOP_PAD )
    {
        if ( val < 0 )
            val = 0;
        
        sys_conf.top_pad = val;

        return 0;
    }

    return 0;
}
    
/*
 * report memory usage
 */

struct ul_mallinfo
rmallinfo ( void )
{
    struct ul_mallinfo  mi;
    internal_size_t     total_size;
    internal_size_t     used_size;
    internal_size_t     free_size;
    heap_t            * heap;

    LOCK( global_heap.mutex );
    
    total_size = global_heap.stat.used_mem;
    used_size  = global_heap.stat.mem_in_use;
    free_size  = global_heap.stat.used_mem - global_heap.stat.mem_in_use;

    UNLOCK( global_heap.mutex );

    LOCK( heap_list.mutex );

    heap = heap_list.heaps;
    while ( heap )
    {
        LOCK( heap->mutex );
        
        total_size += heap->stat.used_mem;
        used_size  += heap->stat.mem_in_use;
        free_size  += heap->stat.used_mem - heap->stat.mem_in_use;
            
        UNLOCK( heap->mutex );
            
        heap = heap->next;
    }
    
    UNLOCK( heap_list.mutex );
    
    mi.arena    = total_size;      /* total space allocated from system */
    mi.ordblks  = 0;               /* number of non-inuse chunks */
    mi.smblks   = 0;               /* unused -- always zero */
    mi.hblks    = 0;               /* number of mmapped regions */
    mi.hblkhd   = 0;               /* total space in mmapped regions */
    mi.usmblks  = 0;               /* unused -- always zero */
    mi.fsmblks  = 0;               /* unused -- always zero */
    mi.uordblks = used_size;       /* total allocated space */
    mi.fordblks = free_size;       /* total non-inuse space */
    mi.keepcost = 0;               /* top-most, releasable (via malloc_trim) space */

    return mi;
}

/*
 * wrapper for mallinfo
 */

#if OVERLOAD_MALLOC == 1
struct mallinfo
mallinfo ( void )
{
    struct mallinfo     mi;
    struct ul_mallinfo  ul_mi = rmallinfo();

    mi.arena    = ul_mi.arena;
    mi.ordblks  = ul_mi.ordblks;
    mi.smblks   = ul_mi.smblks; 
    mi.hblks    = ul_mi.hblks;  
    mi.hblkhd   = ul_mi.hblkhd; 
    mi.usmblks  = ul_mi.usmblks;
    mi.fsmblks  = ul_mi.fsmblks;
    mi.uordblks = ul_mi.uordblks;
    mi.fordblks = ul_mi.fordblks;
    mi.keepcost = ul_mi.keepcost;

    return mi;
}
#endif


/********************************************************************
 ********************************************************************
 **
 ** tracing
 **
 ********************************************************************
 ********************************************************************/

/*
 * allocation trace
 */

static void
rmalloc_init_trace ( const char * name )
{
    assert( name != NULL );

    if ((trace.file = fopen( name, "w" )) == NULL )
    {
        ERROR( __LINE__, "init_trace : could not open \"%s\" (%s)\n",
               name, strerror( errno ) );
        return;
    }
    
    trace.old_rm_used  = 0;
    trace.old_app_used = 0;
    trace.old_free     = 0;
    trace.old_alloc    = 0;
    trace.step         = 0;

    trace_mem = true;
}

/*
 * get status and write trace-entry
 */

static void
rmalloc_trace ()
{
    DEFINE_MUTEX( trace_mutex )
    
    internal_size_t  rm_used     = 0;  /* used by Rmalloc              */
    internal_size_t  app_used    = 0;  /* used by application          */
    internal_size_t  free        = 0;  /* total free memory            */
    internal_size_t  blocks_free = 0;  /* free memory in unused blocks */
    internal_size_t  wild_free   = 0;  /* free memory in wilderness    */
    internal_size_t  alloc       = 0;  /* total allocated memory       */
    
    if ( ! trace_mem )
        return;

    if ( trace.size == 0 )
    {
        heap_t  * heap  = & global_heap;
        block_t * block = NULL;

        while ( heap != NULL )
        {
            LOCK( heap->mutex );
            
            rm_used  += heap->stat.used_mem;
            app_used += heap->stat.mem_in_use;
            free     += heap->stat.used_mem - heap->stat.mem_in_use;
            alloc    += heap->stat.allocated;

            block = heap->blocks.wilderness;

            while ( block != NULL )
            {
                wild_free += BLOCK_SIZE( block );
                block      = block->next;
            }

            blocks_free = free - wild_free;

            UNLOCK( heap->mutex ); /* change of heap_list is improbable but possible */
            
            if ( heap == & global_heap ) heap = heap_list.heaps;
            else                         heap = heap->next;
        }
    }
    else
    {
        unsigned   sclass = size_to_class( trace.size );
        heap_t   * heap   = & global_heap;
        
        while ( heap != NULL )
        {
            LOCK( heap->mutex );
            
            app_used  += heap->stat.used[ sclass ] * size_classes[ sclass ];
            free      += heap->stat.free[ sclass ] * size_classes[ sclass ];
            alloc     += heap->stat.allocated;

            UNLOCK( heap->mutex ); /* change of heap_list is improbable but possible */
            
            if ( heap == & global_heap ) heap = heap_list.heaps;
            else                         heap = heap->next;
        }
    } 

    rm_used     /= 1024;
    app_used    /= 1024;
    blocks_free /= 1024;
    wild_free   /= 1024;
    free        /= 1024;

    LOCK( trace_mutex );
    
    if ( trace.type == TRACE_STEPS )
    {
        if ((( rm_used  > 0 ) ||
             ( app_used > 0 ) ||
             ( free     > 0 ))
            &&
            (( rm_used  != trace.old_rm_used  ) ||
             ( app_used != trace.old_app_used ) ||
             ( free     != trace.old_free     )))
        {
            fprintf( trace.file, "%ld %.3f %.3f %.3f %.3f\n",
                     (long) trace.step,
                     ((double) rm_used)     / 1024.0,
                     ((double) app_used)    / 1024.0,
                     ((double) blocks_free) / 1024.0,
                     ((double) wild_free)   / 1024.0 );
            trace.old_rm_used   = rm_used;
            trace.old_app_used  = app_used;
            trace.old_free      = free;
        }

        trace.step++;
    }
    else if ( trace.type == TRACE_ALLOCATION )
    {
        alloc /= 1024;
        
        if (( alloc != trace.old_alloc )
            &&
            ( (( app_used > 0) ||
               ( free     > 0))
              &&
              (( app_used != trace.old_app_used ) ||
               ( free     != trace.old_free ))))
        {
            fprintf( trace.file, "%.3f %ld %ld\n", ((float)alloc)/1024.0, app_used, free );
            trace.old_app_used = app_used;
            trace.old_free     = free;
            trace.old_alloc    = alloc;
        }
    }

    UNLOCK( trace_mutex );
}

/*
 * finish tracing
 */

static void
rmalloc_finish_trace ()
{
    if ( trace_mem )
    {
        fclose( trace.file );
        trace.file = NULL;
        trace_mem  = false;
    }
}


/********************************************************************
 ********************************************************************
 **
 ** statistical methods
 **
 ********************************************************************
 ********************************************************************/

#define BYTE_TO_MB( n ) (((double)(n)) / (1000.0 * 1000.0))

/*
 * print statistics
 */

static void
print_heap_stat ( heap_t  * heap )
{
    unsigned         i;
    unsigned         nwild;
    internal_size_t  swild;
    block_t        * block;
    
    if ( heap_stat_lvl <= 1 )
        return;
    
    /*
     * statistic for container usage
     */

    if ( heap_stat_lvl >= 3 )
    {
        container_t  * cont;
        
        OUTPUT( " sizeclass |       containers        |\n" );
        OUTPUT( "           |  # empty   | # nonempty |\n" );
        OUTPUT( "-----------+------------+------------+\n" );

        for ( i = 0; i <= NO_OF_SMALL; i++ )
        {
            unsigned  empty    = 0;
            unsigned  nonempty = 0;

            cont = heap->blocks.nonempty[i];
            while ( cont != NULL )
            {
                nonempty++;
                cont = cont->next;
            }

            if (( empty != 0 ) || (nonempty != 0))
                OUTPUT( " %9ld | %10d | %10d |\n", size_classes[i], empty, nonempty );
        }

        OUTPUT( "-----------+------------+------------+\n" );
        OUTPUT( " full blocks : %d\n", heap->blocks.no_of_full );
    }

    /*
     * statistic for small-usage
     */

    if ( heap_stat_lvl >= 3 )
    {
        OUTPUT( " sizeclass |  # used    |   # free   |     kB     |    comment    \n" );
        OUTPUT( "-----------+------------+------------+------------+---------------\n" );
        
        for ( i = NO_OF_SMALL; i < NO_OF_CLASSES; i++ )
        {
            block = heap->blocks.blocks[i];
            
            if ( block != NULL )
            {
                internal_size_t  mem = 0;

                while ( block != NULL )
                {
                    mem += BLOCK_SIZE( block );
                    block = block->next;
                }

                OUTPUT( "%10ld | %10ld | %10ld | %10ld |",
                         (long) size_classes[ i ],
                         (long) heap->stat.used[i],
                         (long) heap->stat.free[i],
                         (long) mem / 1024 );

                if ( heap->stat.used[i] != 0 )
                    OUTPUT( " %ld chunk(s) alive", heap->stat.used[i] );
                
                OUTPUT( "\n" );
            }
        }
    }

    /*
     * sum up wilderness
     */

    nwild = 0;
    swild = 0;
    block = heap->blocks.wilderness;

    while ( block != NULL )
    {
        nwild++;
        swild += BLOCK_SIZE( block );
        block = block->next;
    }
    
    /*
     * output
     */

    if ( heap->stat.max_mem != 0 )
        OUTPUT( "    max allocated            = %.2f MB\n", BYTE_TO_MB(heap->stat.max_mem) );
    
    if ((heap->stat.mem_in_use != 0) || (heap->stat.used_mem != 0))
        OUTPUT( "    mem in use / allocated   = %.2f MB / %.2f MB\n",
                 BYTE_TO_MB(heap->stat.mem_in_use), BYTE_TO_MB(heap->stat.used_mem) );
    
    if ( nwild != 0 )
        OUTPUT( "    wilderness blocks/size   = %d / %.2f MB\n",
                 nwild, BYTE_TO_MB(swild) );

    OUTPUT( "                     top-pad = %.2f MB\n",
             BYTE_TO_MB( heap->top_pad ) );
}

static void
rmalloc_stat ()
{
    heap_t * heap;

    if ( heap_stat_lvl == 0 )
        return;
    
    if ( ! is_initialised )
    {
        ERROR( __LINE__, "stat : heap not initialised\n" );
        return;
    }

#if USE_MMAP == 1
    OUTPUT( "Rmalloc version %s (mmap)\n", RMALLOC_VERSION );
#elif USE_MALLOC == 1
    OUTPUT( "Rmalloc version %s (malloc)\n", RMALLOC_VERSION );
#elif USE_SBRK == 1
    OUTPUT( "Rmalloc version %s (sbrk)\n", RMALLOC_VERSION );
#else
    OUTPUT( "Rmalloc version %s\n", RMALLOC_VERSION );
#endif
    
    /*
     * print statstics about each heap
     */

    if ( heap_stat_lvl >= 2 )
    {
        print_heap_stat( & global_heap );
        
        heap = heap_list.heaps;
        while ( heap != NULL )
        {
            OUTPUT( "heap %ld :\n", (long) heap->id );
            print_heap_stat( heap );
            heap = heap->next;
        }
    }
    
    if ( heap_stat_lvl >= 1 )
    {
        internal_size_t  max_alloc = global_heap.stat.max_mem;
        internal_size_t  max_use   = global_heap.stat.max_in_use;

        heap = heap_list.heaps;
        while ( heap != NULL )
        {
            max_alloc  += heap->stat.max_mem;
            max_use    += heap->stat.max_in_use;
            heap = heap->next;
        }
        
        OUTPUT( "global stat:\n" );

        if ( max_alloc > 0 )
            OUTPUT( "    max. mem used by rmalloc    = %.2f MB\n", BYTE_TO_MB(max_alloc) );
        
        if ( max_use > 0 )
            OUTPUT( "    max. mem used by app (frag) = %.2f MB (%.2f %%)\n",
                     BYTE_TO_MB(max_use), 100.0 * ((double) max_alloc) / ((double) max_use) );
    }
}

#if DEBUG_RMALLOC == 1
static int
check_block ( heap_t * heap, block_t * block )
{
#if DEBUG_BLOCK == 1
    if ( block->magic != BLOCK_MAGIC )
    {
        ERROR( __LINE__, "(check_block) middle block is corrupted\n" );
        return 1;
    }

    if ( block->__pad__ != 0 )
    {
        ERROR( __LINE__, "(check_block) middle block is corrupted\n" );
        return 1;
    }

    /* check EOB size information */
    if ( ! IS_USED( block ) )
    {
        internal_size_t  eob_size = * ((internal_size_t *) ( (((char *) block) + BLOCK_SIZE( block )) -
                                                             sizeof(internal_size_t) ));
            
        if ( BLOCK_SIZE( block ) != eob_size )
        {
            ERROR( __LINE__, "(check_block) end-of-block size information is invalid\n" );
            return 1;
        }
    }

    /* check predecessor information */
    if ( ! IS_PRED_USED( block ) )
    {
        internal_size_t  psize = PRED_SIZE( block );

        if ( psize == 0 )
        {
            ERROR( __LINE__, "(check_block) predecessor size is invalid\n" );
            return 1;
        }
    }
#endif

    return 0;
}

static void
check_rmalloc ()
{
    /*
     * check all heaps, all containers and all blocks
     */

    heap_t  * heap = & global_heap;

    while ( heap != NULL )
    {
        container_t  * container;
        unsigned       i;

        if (( heap->blocks.wilderness       != NULL ) &&
            ( heap->blocks.wilderness->prev != NULL ))
            ERROR( __LINE__, "(check_rmalloc) first wilderness block has predecessor\n" );

        for ( i = 0; i < NO_OF_SMALL; i++ )
        {
            container = heap->blocks.nonempty[i];

            while ( container != NULL )
            {
                sblock_t  * block;
            
#if DEBUG_CONT == 1
                if ( container->magic != CONT_MAGIC )
                {
                    ERROR( __LINE__, "(check_rmalloc) nonempty container is corrupted\n" );
                    break;
                }
#endif
            
                block = container->blocks;

#if DEBUG_CONT == 1
                if ( block == NULL )
                    ERROR( __LINE__, "(check_rmalloc) nonempty container has empty blocklist\n" );
#endif
                
                while ( block != NULL )
                {
#if DEBUG_BLOCK == 1
                    if ( block->magic != BLOCK_MAGIC )
                    {
                        ERROR( __LINE__, "(check_rmalloc) small nonempty block is corrupted\n" );
                        break;
                    }
#endif
                    block = block->next;
                }

                container = container->next;
            }
        }

        container = heap->blocks.full;
            
        while ( container != NULL )
        {
            sblock_t  * block;
            
#if DEBUG_CONT == 1
            if ( container->magic != CONT_MAGIC )
            {
                ERROR( __LINE__, "(check_rmalloc) full container is corrupted\n" );
                break;
            }
#endif
            
            block = container->blocks;

#if DEBUG_CONT == 1
            if ( block == NULL )
                ERROR( __LINE__, "(check_rmalloc) full container has empty blocklist\n" );
#endif
                
            while ( block != NULL )
            {
#if DEBUG_BLOCK == 1
                if ( block->magic != BLOCK_MAGIC )
                {
                    ERROR( __LINE__, "(check_rmalloc) small full block is corrupted\n" );
                    break;
                }
#endif
                block = block->next;
            }
            
            container = container->next;
        }
        
        for ( i = NO_OF_SMALL; i < NO_OF_CLASSES; i++ )
        {
            block_t * block = heap->blocks.blocks[i];

            while ( block != NULL )
            {
                if ( check_block( heap, block ) )
                    break;
                
                block = block->next;
            }
        }

        /* choose next heap */
        if ( heap == & global_heap )
            heap = heap_list.heaps;
        else
            heap = heap->next;
    }
}
#endif

#ifdef __cplusplus
}
#endif

/*
 * message output
 */
static void
OUTPUT ( const char * msg, ... )
{
    va_list ap;

    fprintf( stdout, "(rmalloc) " );
    
    va_start( ap, msg );
    vfprintf( stdout, msg, ap );
    va_end(ap);
}

static void
ERROR  ( int lineno, const char * msg, ... )
{
    va_list ap;

    __RMALLOC_debug();

    fprintf( stderr, "(rmalloc:%d) : ", lineno );
    
    va_start( ap, msg );
    vfprintf( stderr, msg, ap );
    va_end(ap);
}

/*
 * local debug hook called on error detection
 */
void __RMALLOC_debug () {}

/********************************************************************
 ********************************************************************
 **
 ** C++ memory management
 **
 ********************************************************************
 ********************************************************************/

/*
 * overload global new/delete
 */
  
#if defined(__cplusplus) && OVERLOAD_NEW == 1

#include <new>

void *
operator new ( size_t n ) throw (std::bad_alloc)
{
    void * p = r_malloc( n );

    if (( n != 0 ) && ( p == NULL ))
        throw std::bad_alloc();

    return p;
}

void *
operator new[] ( size_t n ) throw (std::bad_alloc)
{
    void * p = r_malloc( n );

    if (( n != 0 ) && ( p == NULL ))
        throw std::bad_alloc();

    return p;
}

void
operator delete ( void * p ) throw ()
{
    r_free( p );
}

void
operator delete[] ( void * p ) throw ()
{
    r_free( p );
}

#endif  /* defined(__cplusplus) && OVERLOAD_NEW == 1 */
