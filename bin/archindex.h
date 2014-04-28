/***********************************************************************

  archindex.h - Include file for archindex.c

***********************************************************************/

#include "rmalloc.h"


/* Array parameters */
#define FILENAME_LEN     10000
#define LINELEN          50000
#define MAX_TOKENS       ( LINELEN / 8 )
#define EC_LEN              60
#define NAME_LEN           512
#define PFAM_LEN            20
#define SPECIES_LEN          8
#define UNIPROT_LEN         15


/* Output file pointer */
FILE *outfile_ptr;

/* Run types */     
#define COMPILE_INDEXES      0
#define SEARCH               1

/* Search types */ 
#define NONE                -1
#define PFAM_ID              0
#define SEQUENCE_ID          1

/* Parameters defining how much data to return from the searches */
#define MAX_ARCHITECTURES  150
#define MAX_SEQUENCES      100

/* Maximum number of domains colours for PDBsum pfamstr.html */
#define MAX_DOMAIN_COLOURS 1264
#define PFAMB_COLOUR_START   10

/* Key types */
#define KEY_PFAM_ID          0
#define KEY_UNIPROT_ACC      1
#define KEY_UNIPROT_ID       2

#define NKEYS                3

static char *DAT_FILE_NAME[2][NKEYS] = {
  { "archpfam.dat", "archuacc.dat", "archuid.dat" },
  { "catharch.dat", "cathuacc.dat", "cathuid.dat" }
};

static char *IDX_FILE_NAME[2][NKEYS] = {
  { "archpfam.idx", "archuacc.idx", "archuid.idx" },
  { "archcath.idx", "cathuacc.idx", "cathuid.idx" }
};

#define NKEY_SPECIES         6

static char *KEY_SPECIES[NKEY_SPECIES] = {
  "_HUMAN",
  "_MOUSE",
  "_DROME",
  "_CAEEL",
  "_YEAST",
  "_ECOLI"
};

/* Indexing gap */
#define INDEX_GAP         5000

/* Parameters for scoring architectures */
#define MAX_SCORE          100
#define PFAMA_DIST           4
#define PFAMB_DIST           3
#define SPLIT_CATH_DIST      2

/* Architecture alignment parameters */
#define MAX_COMPARE        100
#define PARENT               0
#define CURRENT              1
#define FIRST                0
#define SECOND               1
#define GAP_PENALTY          5
#define PARENT_SUBSTRING_DIST 3

/* Pfam domain types */
#define PFAM_A                0
#define PFAM_B                1

/* Domain types */
#define PFAM                 0
#define CATH                 1

static char *DOM_SEP[2] = { ".", "_" };

/* File types */
#define DOMSEQS              0
#define DOMAINS              1

/* Pfam domain operators */
#define OR                   0
#define AND                  1

/* Output file types */
#define ARCHSCHEMA_OUTPUT    0
#define DAT_OUTPUT           1

/* Fields in archstruc.uniprot.joined.by.accession file */
#define F_UNIPROT_ACC        0
#define F_ARCHITECTURE       1
#define F_UNIPROT_ID         2
#define F_REVIEW_FLAG        3
#define F_PROTEIN_NAME       4
#define F_GENE               5
#define F_SPECIES            6
#define F_SEQ_LENGTH         7
#define F_PDB_CODES          8
#define F_COVERAGE           9

/* Run status */
#define RUN_ERROR           -1
#define RUN_SUCCESSFUL       0

/* Error types */
#define NO_PARAMETERS        0
#define MEMORY_ALLOCATION    1
#define FILE_NOT_FOUND       2
#define FILE_READ_ERROR      3
#define CODE_NOT_FOUND       4
#define NO_HITS              5
#define NO_PARENT            6

#define NERROR_TYPES         7

static char *ERROR_TYPE[NERROR_TYPES] = {
  "Parameter error",
  "Memory error",
  "File not found",
  "File read error",
  "Identifier not found",
  "No hits",
  "Index error"
};

/* Enzyme sort types */
#define EC_CODE_SORT         0
#define EC_COUNT_SORT        1

/* Structure for input parameters
   ------------------------------ */
struct parameters
{
  int                     run_type;
  BOOL                    cath_domains;
  BOOL                    old_cath;
  int                     npfam;
  int                     max_architectures;
  int                     max_sequences;
  char                    *architecture;
  char                    *uniprot_acc;
  char                    *uniprot_id;
  char                    *species;
  int                     operator;
  BOOL                    pdb_only;
  BOOL                    reviewed_only;
  int                     output_type;
  char                    out_name[FILENAME_LEN];
  int                     get_ssg;
  struct pfam             *first_pfam_ptr;
  char                    *filter_seqs_file;
  struct filterseq        *first_filterseq_ptr;
  int                     combine;
  int                     run_64;
};

/* Structure for file and directory names
   -------------------------------------- */
struct file_data
{
  char                    *archschema_dir;
  char                    *cath_data_dir;
  char                    *enzyme_dir;
};

/* Structure for each arch record
   ------------------------------ */
struct arch
{
  char                    *key;
  long int                offset;
  struct arch             *next_arch_ptr;
};

/* Structure for each pfam record
   ------------------------------ */
struct pfam
{
  char                    *pfam_id;
  char                    *short_name;
  char                    *pfam_name;
  //  long int                offset[2];
  long int                offset;
  BOOL                    copy;
  BOOL                    extra;
  int                     nseqs;
  int                     colour;
  int                     pstn;
  struct offset           *first_offset_ptr;
  struct pfam             *next_pfam_ptr;
};

/* Structure for each filter sequence record
   ----------------------------------------- */
struct filterseq
{
  char                    *seq_id;
  struct filterseq        *next_filterseq_ptr;
};

/* Structure for each offset record
   ------------------------------ */
struct offset
{
  long int                offset;
  struct offset           *next_offset_ptr;
};

/* Structure for each hit record
   ---------------------------------- */
struct hit
{
  char                    *uniprot_acc;
  char                    *uniprot_id;
  char                    *architecture;
  char                    *gene;
  char                    *protein_name;
  int                     seq_len;
  char                    *pdb_codes;
  char                    *coverage;
  struct hit              *next_hit_ptr;
};

/* Structure for each archit record
   -------------------------------- */
struct archit
{
  char                    *architecture;
  int                     score;
  int                     nhits;
  struct hitlist          *first_hitlist_ptr;
  struct enzlist          *first_enzlist_ptr;
  struct archit           *next_archit_ptr;
};

/* Structure for each hitlist record
   ------------------------------ */
struct hitlist
{
  struct hit              *hit_ptr;
  struct hitlist          *next_hitlist_ptr;
};

/* Structure for each enzlist record
   ------------------------------ */
struct enzlist
{
  struct enzyme           *enzyme_ptr;
  int                     nseqs;
  struct enzlist          *next_enzlist_ptr;
};

/* Structure for each species record
   --------------------------------- */
struct species
{
  char                     *species_id;
  char                     *species_name;
  int                      nseqs;
  struct species           *next_species_ptr;
};

/* Structure for each enzyme record
   -------------------------------- */
struct enzyme
{
  char                     *ec_number;
  char                     *sort_number;
  char                     *enzyme_name;
  int                      narch;
  int                      nseqs;
  int                      ntotal;
  struct enzyme            *next_enzyme_ptr;
};

