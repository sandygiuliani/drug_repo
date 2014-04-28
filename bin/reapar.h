/***********************************************************************

  reapar.h - Include file for reapar.c

***********************************************************************/

#include <string.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <ctype.h>

/* Flag-values */     
#define FALSE                0
#define TRUE                 1

/* Array parameters */
#define FILENAME_LEN       512
#define LINELEN           2000
#define MAX_PARAM         1000
#define MAX_PDB_DIR         27
#define STRING_LEN        1000

/* Error status */
#define OK                   0
#define NO_TOKENS            1
#define NAME_TOO_LONG        2
#define FILE_MISSING         2
#define NO_MORE_TOKENS       3
#define TOKEN_MISSING        4


/* Structure for tokens
   -------------------- */
struct token
{
  char               *name;
  char               *value;
  struct token       *next_token_ptr;
};

