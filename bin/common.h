/***********************************************************************

  common.h - Include file for common global vraibles

***********************************************************************/

/* Standard system include files */
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <ctype.h>
#include <time.h>
#include <limits.h>
#include <sys/stat.h>

/* Include file for reading gzipped files */
#ifdef TEST
#else
#include <zlib.h>
#endif

/* Flag-values */     
#define FALSE                0
#define TRUE                 1

typedef short  BOOL;

/* Directory separator */
static char DIR_SEP[] = "/";

/* NULL character */
static char null = '\0';

/* Field for storing current date and time 
   Access by:  time(&current_time);
               printf("Current time: %s\n",ctime(&current_time)); */
time_t current_time;

/* Operating system we are running on */
#define LINUX_32              0
#define LINUX_64              1

/* PDB file types */
#define FIND_FIRST           -9
#define NOT_FOUND            -1
#define STANDARD_PDB          0
#define UPLOAD_PDB            1
#define MCSG_PDB              2

/* PDB code length */
#define PDB_CODELEN          4

/* Template characters (length must match PDB code length above) */
static char TEMPLATE_CHARS[] = "abcd";

/* PDB file locations */
#define PDB_CURRENT          0
#define PDB_OBSOLETE         1
#define PDB_MODELS           2
#define PDB_MODELS_OBSOLETE  3
#define PDB_MISC             4
#define PDB_UPLOAD           5
#define PDB_MCSG             6

#define NPDB_DIR             7

/* Names of the PDB directories in CATHPARAM */
static char *PDB_DIR_TEMPLATE[NPDB_DIR] = { 
  "FTP_PDB_TEMPLATE", "FTP_OBSOLETE_TEMPLATE",
  "FTP_MODELS_TEMPLATE", "FTP_OBSOLETE_MODELS_TEMPLATE",
  "MISC_PDB_TEMPLATE", "UPLOAD_PDB_TEMPLATE", "MCSG_PDB_TEMPLATE"
};

/* PDBsum file locations */
#define PDBSUM_CURRENT       0
#define PDBSUM_UPLOAD        1

#define NPDBSUM_DIR          2

/* Names of the PDBsum directory templates */
static char *PDBSUM_DIR_TEMPLATE[NPDBSUM_DIR] = { 
  "PDBSUM_REAL_TEMPLATE", "UPLOAD_PDBSUM_TEMPLATE"
};

/* Names of the PDBsum data directories */
static char *PDBSUM_DATA_DIR[NPDBSUM_DIR] = { 
  "PDBSUM_REAL_DATA_DIR", "UPLOAD_PDBSUM_DATA_DIR"
};

/* Executable file locations */
#define ACROREAD_EXE         0
#define BABEL_EXE            1
#define CONVERT_EXE          2
#define CRYPT_EXE            3
#define ESEARCH_PROG         4
#define FASTA                5
#define GHOSTSCRIPT_EXE      6
#define HBADD_EXE            7
#define HBPLUS_EXE           8
#define HSEARCH_PROG         9
#define JAVA_EXE            10
#define LIGPLOT_EXE         11
#define LIGSEARCH_PYMOL_EXE 12
#define NACCESS_EXE         13
#define NCBI_FETCH          14
#define NEW_PYMOL_EXE       15
#define OLD_CONVERT_EXE     16
#define PDFEXTRACT_EXE      17
#define PDFIMAGES_EXE       18
#define PPMQUANT_EXE        19
#define PPMTOGIF_EXE        20
#define PSTOIMG_EXE         21
#define PSTOTEXT_EXE        22
#define PYMOL_EXE           23
#define PYMOL_64_EXE        24
#define RASTER3D            25
#define RUNROMLAS_CMND      26
#define SPEEDFILL_EXE       27
#define SVMTEST_EXE         28
#define SVMTORCH_EXE        29
#define SW_JAVA             30
#define VIEWFIGS_PL         31

#define NEXE_DIR            32

/* Names of the PDBsum directories in CATHPARAM */
static char *EXE_NAME[NEXE_DIR] = { 
  "ACROREAD_EXE", "BABEL_EXE", "CONVERT_EXE", "CRYPT_EXE", "ESEARCH_PROG",
  "FASTA",
  "GHOSTSCRIPT_EXE", "HBADD_EXE", "HBPLUS_EXE", "HSEARCH_PROG", "JAVA_EXE",
  "LIGPLOT_EXE", "LIGSEARCH_PYMOL_EXE", "NACCESS_EXE", "NCBI_FETCH",
  "NEW_PYMOL_EXE", "OLD_CONVERT_EXE", "PDFEXTRACT_EXE",
  "PDFIMAGES_EXE", "PPMQUANT_EXE", "PPMTOGIF_EXE", "PSTOIMG_EXE",
  "PSTOTEXT_EXE", "PYMOL_EXE", "PYMOL_64_EXE", "RASTER3D",
  "RUNROMLAS_CMND", "SPEEDFILL_EXE", "SVMTEST_EXE", "SVMTORCH_EXE",
  "SW_JAVA", "VIEWFIGS_PL"
};

/* Other file locations */
#define ARCHCATH_OLD_REAL_DATA_DIR  0
#define ARCHCATH_REAL_DATA_DIR      1
#define ARCHSCHEMA_REAL_DATA_DIR    2
#define BINDDB_DIR                  3
#define CATH_DATA_DIR               4
#define CATH_LIST                   5
#define CITEXPLORE_CLASSPATH        6
#define CULLED_DIR                  7
#define DIR_SEPARATOR               8
#define DNA_REPORTS_DIR             9
#define DNA_TEMPLATES_DIR          10
#define DOMALL_FILE                11
#define DRAT_REAL_DATA_DIR         12
#define DRUGBANK_DIR               13
#define DRUGCARDS_FILE             14
#define DRUGPORT_REAL_DATA_DIR     15
#define DRUGPORT_REAL_DIR          16
#define EBI                        17
#define EBI_HOME_PAGE              18
#define EBI_INCLUDE_DIR            19
#define EHEADERS_FILE              20
#define ENTRIES_DIR                21
#define ENTRIES_REAL_DIR           22
#define ENZ_REPORTS_DIR            23
#define ENZ_TEMPLATES_DIR          24
#define ENZYME_GIF_DIR             25
#define ENZYME_HOME_PAGE           26
#define ENZYME_MOLGIF_DIR          27
#define ENZYME_REAL_DIR            28
#define ENZYME_REAL_DATA_DIR       29
#define ENZYME_TEMPLATES_DIR       30
#define EXPASY_ENZYME_DIR          31
#define FAKE_PDB_ENTRY             32
#define GENE3D_DATA_DIR            33

#define GENOMES_LOCAL_DIR          34
#define GENOMES_LOCAL_DAT_DIR      35
#define GET_ECPAGE_PL              36
#define GET_PDBSUM_ENTRY           37
#define GET_UNIPROT_DATA           38
#define GOA_DIR                    39
#define GOA_LINK                   40
#define GO_LOCAL_DIR               41
#define GO_NAMES_FILE              42
#define GO_ONTOLOGY_DIR            43
#define HETGIF_DIR                 44
#define HETGIF_REAL_DIR            45
#define HETGROUP_DICTIONARY        46
#define HETGROUPS_DIR              47
#define HETGROUPS_REAL_DIR         48
#define HTTP_BASE                  49
#define JAX_LIB                    50
#define KEGG_COMPOUND_FILE         51
#define KEGG_DATA_DIR              52
#define KEGG_DRUG_FILE             53
#define KEGG_GIF_DIR               54
#define KEGG_MOL_DIR               55
#define LIG_REPORTS_DIR            56
#define LIG_TEMPLATES_DIR          57
#define LIGANDS_DIR                58
#define LIGANDS_REAL_DIR           59
#define LIG_CLUSTERS_REAL_DIR      60
#define LIGPLOT_BIN_DIR            61
#define LIGSEARCH_REAL_RESULTS     62
#define LOCAL_PATTERNS_FILE        63
#define NEW_CONSURF_DIR            64
#define NISTHET_NAME               65
#define NISTLIST_NAME              66
#define PDB_XML_DIR                67
#define PDBLIB_FILE                68
#define PDBSUM_CGI_DIR             69
#define PDBSUM_DIR                 70
#define PDBSUM_DOCS_DIR            71
#define PDBSUM_GENDATA_DIR         72
#define PDBSUM_GIF_DIR             73
#define PDBSUM_HOME_PAGE           74
#define PDBSUM_PROFUNC_DIR         75
#define PDBSUM_REAL_DATA_DIR       76
#define PDBSUM_REAL_TEMPLATES_DIR  77
#define PDBSUM_TEMPLATES_DIR       78
#define PDBSUM_TEMPLATES_GIF_DIR   79
#define PFAM_DATA_DIR              80
#define PQS_DIR                    81
#define PQS_DATA_DIR               82
#define PROCHECK_EXE_64_DIR        83
#define PROFUNC_CGI_DIR            84
#define PROFUNC_DATA_DIR           85
#define PROFUNC_TEMPLATES_DIR      86
#define PROSITE_FILE               87
#define PUBMED_REAL_DIR            88
#define RESULTS_DIR                89
#define S20REP_DOMAINS             90
#define SAS_CGI_DIR                91
#define SAS_CSS_DIR                92
#define SAS_EXE_DIR                93
#define SAS_GIF_DIR                94
#define SAS_HOME_PAGE              95
#define SAS_REAL_TEMPLATES_DIR     96
#define SPECIES_DIR                97
#define SPECIES_HTML_FILE          98
#define SWISS_CODE                 99
#define SWISS_NAME                100
#define SWISSPROT_DIR             101
#define SWISSPROT_LOCAL_DIR       102
#define TAXONOMY_DIR              103
#define TMP_DIR                   104
#define TOF_MATCHES_FILE          105
#define UPLOAD_PDBLIB_FILE        106
#define UPLOAD_PROSITE_FILE       107
#define UPLOAD_TOF_MATCHES_FILE   108
#define URL_NIST_DIRECT           109
#define VARIANTS_REAL_TEMPLATE    110
#define VARIANTS_TEMPLATE         111
#define ZPDB_TO_GO_FILE           112
#define ZSPROT_FASTA              113
#define ZTREMBL_FASTA             114
#define ZUNIPROT_TO_GO_FILE       115

#define NOTHER_DIR                116

/* Names of the PDBsum directories in CATHPARAM */
static char *OTHER_DIR_NAME[NOTHER_DIR] = { 
  "ARCHCATH_OLD_REAL_DATA_DIR", "ARCHCATH_REAL_DATA_DIR",
  "ARCHSCHEMA_REAL_DATA_DIR", "BINDDB_DIR", "CATH_DATA_DIR",
  "CATH_LIST", "CITEXPLORE_CLASSPATH", "CULLED_DIR", "DIR_SEP",
  "DNA_REPORTS_DIR", "DNA_TEMPLATES_DIR", "DOMALL_FILE",
  "DRAT_REAL_DATA_DIR", "DRUGBANK_DIR", "DRUGCARDS_FILE",
  "DRUGPORT_REAL_DATA_DIR", "DRUGPORT_REAL_DIR", "EBI", "EBI_HOME_PAGE",
  "EBI_INCLUDE_DIR", "EHEADERS_FILE", "ENTRIES_DIR", "ENTRIES_REAL_DIR",
  "ENZ_REPORTS_DIR", "ENZ_TEMPLATES_DIR", "ENZYME_GIF_DIR",
  "ENZYME_HOME_PAGE", "ENZYME_MOLGIF_DIR", "ENZYME_REAL_DIR",
  "ENZYME_REAL_DATA_DIR", "ENZYME_TEMPLATES_DIR", "EXPASY_ENZYME_DIR",
  "FAKE_PDB_ENTRY", "GENE3D_DATA_DIR", "GENOMES_LOCAL_DIR",
  "GENOMES_LOCAL_DAT_DIR", "GET_ECPAGE_PL",
  "GET_PDBSUM_ENTRY", "GET_UNIPROT_DATA", "GOA_DIR", "GOA_LINK",
  "GO_LOCAL_DIR", "GO_NAMES_FILE", "GO_ONTOLOGY_DIR", "HETGIF_DIR",
  "HETGIF_REAL_DIR", "HETGROUP_DICTIONARY", "HETGROUPS_DIR",
  "HETGROUPS_REAL_DIR", "HTTP_BASE", "JAX_LIB", "KEGG_COMPOUND_FILE",
  "KEGG_DATA_DIR", "KEGG_DRUG_FILE", "KEGG_GIF_DIR", "KEGG_MOL_DIR",
  "LIG_REPORTS_DIR", "LIG_TEMPLATES_DIR",
  "LIGANDS_DIR", "LIGANDS_REAL_DIR", "LIG_CLUSTERS_REAL_DIR",
  "LIGPLOT_BIN_DIR", "LIGSEARCH_REAL_RESULTS", "LOCAL_PATTERNS_FILE",
  "NEW_CONSURF_DIR", "NISTHET_NAME", "NISTLIST_NAME", "PDB_XML_DIR",
  "PDBLIB_FILE", "PDBSUM_CGI_DIR", "PDBSUM_DIR",
  "PDBSUM_DOCS_DIR", "PDBSUM_GENDATA_DIR", "PDBSUM_GIF_DIR",
  "PDBSUM_HOME_PAGE", "PDBSUM_PROFUNC_DIR", "PDBSUM_REAL_DATA_DIR",
  "PDBSUM_REAL_TEMPLATES_DIR",
  "PDBSUM_TEMPLATES_DIR", "PDBSUM_TEMPLATES_GIF_DIR", "PFAM_DATA_DIR",
  "PQS_DIR", "PQS_DATA_DIR", "PROCHECK_EXE_64_DIR", "PROFUNC_CGI_DIR",
  "PROFUNC_DATA_DIR", "PROFUNC_TEMPLATES_DIR", "PROSITE_FILE",
  "PUBMED_REAL_DIR", "RESULTS_DIR", "S20REP_DOMAINS",
  "SAS_CGI_DIR", "SAS_CSS_DIR", "SAS_EXE_DIR", "SAS_GIF_DIR",
  "SAS_HOME_PAGE", "SAS_REAL_TEMPLATES_DIR", "SPECIES_DIR",
  "SPECIES_HTML_FILE", "SWISS_CODE", "SWISS_NAME", "SWISSPROT_DIR",
  "SWISSPROT_LOCAL_DIR", "TAXONOMY_DIR", "TMP_DIR", "TOF_MATCHES_FILE",
  "UPLOAD_PDBLIB_FILE", "UPLOAD_PROSITE_FILE",
  "UPLOAD_TOF_MATCHES_FILE", "URL_NIST_DIRECT",
  "VARIANTS_REAL_TEMPLATE", "VARIANTS_TEMPLATE", "ZPDB_TO_GO_FILE",
  "ZSPROT_FASTA", "ZTREMBL_FASTA", "ZUNIPROT_TO_GO_FILE"
};

/* Array parameters */
#define C_FILENAME_LEN    1024

/* Coordinate subscripts */
#define X                    0
#define Y                    1
#define Z                    2

/* Large and small constants */
#ifndef INFINITY
#define INFINITY   (9999999999.9)
#endif
#ifndef INT_INFINITY
#define INT_INFINITY   INT_MAX
#endif
#define TINY       (0.0000000001)

/* Angle constants */
#define PI             ((double)  3.141592654)
#define RADDEG         ((double)  (180.0 / PI))
#define PI_BY_2        ((double)  (PI / 2.0))
#define TWO_PI         ((double)  (2.0 * PI))
#define SQRT_TWO_PI    ((double)  (sqrt(TWO_PI)))
#define ROOT2          ((double)  sqrt(2.0))
#define ONE_OVER_ROOT2 ((double)  (1.0 / ROOT2))

/* Other constants */
#define E_VALUE        ( exp(1) )
#define E_VALUE_SQRD   ( E_VALUE * E_VALUE )

/* Open output file flags */
#define TERMINATE      TRUE
#define WARN_ONLY      FALSE

/* RGB colour codes */
#define BLACK                 0
#define WHITE             ( 256 * 256 * 256 - 1 )
#define OFF_WHITE         ( WHITE - 1 )
#define RED               ( 256 * 256 * 255 )
#define GREEN             ( 256 * 255 )
#define BLUE              ( 255 )
#define BROWN             ( 256 * 256 * 204 + 256 * 153 +   0 )
#define CREAM             ( 256 * 256 * 255 + 256 * 255 + 179 )
#define CYAN              ( 256 * 256 *   0 + 256 * 255 + 255 )
#define DEEP_PINK         ( 256 * 256 * 255 + 256 *  20 + 147 )
#define GOLD3             ( 256 * 256 * 205 + 256 * 173 +   0 )
#define GREY1             ( 256 * 256 * 128 + 256 * 128 + 128 )
#define GREY2             ( 256 * 256 * 168 + 256 * 168 + 168 )
#define GREY3             ( 256 * 256 * 188 + 256 * 188 + 188 )
#define GREY4             ( 256 * 256 * 208 + 256 * 208 + 208 )
#define GREY5             ( 256 * 256 * 228 + 256 * 228 + 228 )
#define LAWN_GREEN        ( 256 * 256 * 124 + 256 * 252 +   0 )
#define PINK              ( 256 * 256 * 255 + 256 * 105 + 180 )
#define PURPLE            ( 256 * 256 * 204 + 256 * 102 + 255 )
#define YELLOW            ( 256 * 256 * 255 + 256 * 255 +   0 )

/* Cutoff for deciding whether a colour is a shade of grey */
#define GREY_CUTOFF                 10

/* Shadow parameters */
#define NSHADOW_COLOURS              5
#define NPLOT_SIZES                  2

/* Plot sizes */
#define STANDARD_PLOT                0
#define GO_LARGE                     1

static const int SHADOW_RANGE[NPLOT_SIZES] = { 4, 8 };

static const int SHADOW_COLOUR[NSHADOW_COLOURS]
= { GREY1, GREY2, GREY3, GREY4, GREY5 };

static const int NEMPTY[NPLOT_SIZES][NSHADOW_COLOURS] =
  { 0, 1, 2, 3, 4,
    0, 3, 5, 6, 7 };

/* Chain colours */
#define NCHAIN_COLOURS        8
static char *CHAIN_LIST
= "ABCDEFGHIJKLMNOPQRSTUVWXYZ123456 7890!#_-=:<>|";

static char *CHAIN_COLOUR[NCHAIN_COLOURS]
= {
  "purple", "red", "brown", "pink",
  "blue", "green", "cyan", "yellow" 
};

/* Defined colours */
#define NCOLOURS             21
static char *COLOUR_LIST[] = {
  "black", "red", "green", "blue", "yellow", "purple", "brown",
  "sky_blue", "orange", "cyan", "grey", "white", "magenta", "pink",
  "dark_grey", "cream", "paleblue","greenblue", "redorange", "violet",
  "peach"
};
static double RGB_TRIPLET[NCOLOURS][3] = {
  {0.000,  0.000,  0.000},  /* black */
  {1.000,  0.000,  0.000},  /* red */
  {0.000,  1.000,  0.000},  /* green */
  {0.000,  0.000,  1.000},  /* blue */
  {1.000,  1.000,  0.000},  /* yellow */
  {0.800,  0.400,  1.000},  /* purple */
  {0.800,  0.600,  0.000},  /* brown */
  {0.000,  0.600,  1.000},  /* sky_blue */
  {1.000,  0.400,  0.000},  /* orange */
  {0.200,  0.800,  1.000},  /* cyan */
  {0.500,  0.500,  0.500},  /* grey */
  {1.000,  1.000,  1.000},  /* white */
  {0.938,  0.000,  0.500},  /* magenta */
  {1.000,  0.412,  0.706},  /* pink */
  {0.200,  0.200,  0.200},  /* dark_grey */
  {1.000,  1.000,  0.700},  /* cream */
  {0.680,  0.850,  0.900},  /* paleblue */
  {0.200,  0.600,  0.400},  /* greenblue */
  {1.000,  0.400,  0.000},  /* redorange */
  {0.930,  0.510,  0.930},  /* violet */
  {1.000,  0.800,  0.600}   /* peach */
};

/* Colours for CATH domains */
#define NCATH_COLOURS              8

static int CATH_COLOUR[NCATH_COLOURS] = {
  45,  // Black
  1,   // Red
  2,   // Blue
  0,   // Green
  37,  // SaddleBrown
  23,  // Orange
  5,   // Cyan
  13  // HotPink
};

/* Colours for Pfam domains */
#define NPFAM_COLOURS             46
static char *COLOUR_DEFN[NPFAM_COLOURS + 1][3] = {
  {"Green",           "#00FF00", "0 255 0"},     //  0
  {"Red",             "#FF0000", "255 0 0"},     //  1
  {"Blue",            "#0000FF", "0 0 255"},     //  2
  {"Gold",            "#FFD700", "255 215 0"},   //  3
  {"Purple",          "#A020F0", "160 32 240"},  //  4
  {"Cyan",            "#00FFFF", "0 255 255"},   //  5
  {"Firebrick",       "#B22222", "178 34 34"},   //  6
  {"LimeGreen",       "#32CD32", "50 205 50"},   //  7
  {"DeepSkyBlue",     "#00BFFF", "0 191 255"},   //  8
  {"LightGoldenrod1", "#FFEC8B", "255 236 139"}, //  9
  {"OrangeRed",       "#FF4500", "255 69 0"},    // 10
  {"LightGreen",      "#90EE90", "144 238 144"}, // 11
  {"SkyBlue",         "#87CEEB", "135 206 235"}, // 12
  {"HotPink",         "#FF69B4", "255 105 180"}, // 13
  {"SeaGreen",        "#2E8B57", "46 139 87"},   // 14
  {"MidnightBlue",    "#191970", "25 25 112"},   // 15
  {"DarkRed",         "#8B0000", "139 0 0"},     // 16
  {"Peru",            "#CD853F", "205 133 63"},  // 17
  {"SpringGreen",     "#00FF7F", "0 255 127"},   // 18
  {"Magenta",         "#FF00FF", "255 0 255"},   // 19
  {"DimGrey",         "#696969", "105 105 105"}, // 20
  {"DarkGoldenrod",   "#B8860B", "184 134 11"},  // 21
  {"Maroon",          "#B03060", "176 48 96"},   // 22
  {"Orange",          "#FFA500", "255 165 0"},   // 23
  {"OliveDrab",       "#6B8E23", "107 142 35"},  // 24
  {"CornflowerBlue",  "#6495ED", "100 149 237"}, // 25
  {"DeepPink",        "#FF1493", "255 20 147"},  // 26
  {"DarkSlateBlue",   "#483D8B", "72 61 139"},   // 27
  {"SandyBrown",      "#F4A460", "244 164 96"},  // 28
  {"DarkKhaki",       "#BDB76B", "189 183 107"}, // 29
  {"Tomato",          "#FF6347", "255 99 71"},   // 30
  {"SteelBlue",       "#4682B4", "70 130 180"},  // 31
  {"Chocolate",       "#D2691E", "210 105 30"},  // 32
  {"PaleGreen",       "#98FB98", "152 251 152"}, // 33
  {"MediumOrchid",    "#BA55D3", "186 85 211"},  // 34
  {"DarkSlateGray",   "#2F4F4F", "47 79 79"},    // 35
  {"Salmon",          "#FA8072", "250 128 114"}, // 36
  {"SaddleBrown",     "#8B4513", "139 69 19"},   // 37
  {"SlateGrey",       "#708090", "112 128 144"}, // 38
  {"LawnGreen",       "#7CFC00", "124 252 0"},   // 39
  {"Burlywood",       "#DEB887", "222 184 135"}, // 40
  {"DarkGreen",       "#006400", "0 100 0"},     // 41
  {"SlateBlue",       "#6A5ACD", "106 90 205"},  // 42
  {"DarkOliveGreen",  "#556B2F", "85 107 47"},   // 43
  {"Grey",            "#BEBEBE", "190 190 190"}, // 44
  {"Black",           "#000000", "0 0 0"},       // 45
  {"White",           "#FFFFFF", "255 255 255"}  // 46
};
static char COLOUR_CODE[]
= "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrst";
#define BLACK_COL            45
#define WHITE_COL            46
#define PFAMA_STARTCOL        0
#define PFAMB_STARTCOL       10
#define HALF_TEXT           330.0

/* Colour of dots indicating break in sequence */
static int RGB_BREAK[3] = { 64, 64, 64 };

/* Colour of labels showing names of Pfam domains */
static int RGB_TEXT[3] = { 170, 170, 170 };

/* Structure for common directory names
   ------------------------------------ */
struct c_file_data
{
  char                    *pdb_template[NPDB_DIR];
  char                    *pdbsum_exe_dir;
  char                    *pdbsum_gen_dir;
  char                    *pdbsum_template[NPDBSUM_DIR];
  char                    *pdbsum_data_dir[NPDBSUM_DIR];
  char                    *exe_name[NEXE_DIR];
  char                    *other_dir[NOTHER_DIR];
};

