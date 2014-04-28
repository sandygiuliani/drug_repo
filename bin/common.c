/***********************************************************************

  common.c - Common routines used by all programs

***********************************************************************

Date:-         22 Oct 2009

Written by:-   Roman Laskowski

----------------------------------------------------------------------*/

/* Compilation
   -----------

cc -c common.c

*/

/* Routines
   --------

   -> store_string
   -> c_get_paths
         -> get_param (in reapar.c)
   -> check_file
   -> fix_upper
   -> form_filename
         -> store_string
   -> find_pdbsum_dir
         -> form_filename
   -> fix_token
   -> free_tokens
   -> get_random_number
   -> get_tokens
   -> open_output_file
   -> parse_uniprot_name
   -> perform_splice
   -> remove_blanks
   -> remove_leading_blanks
   -> string_chop
   -> string_truncate
   -> to_lower
   -> to_upper
   -> tokenize
   -> type_set
   -> write_ppm_file
   -> convert_to_gif

*/


#include "common.h"

/* Prototypes
   ---------- */

/* In reapar.c */
int get_param(char *token_name,char *file_name,int filename_len,
              int exit_on_error);

/***********************************************************************

store_string  -  Allocate memory for given string and store pointer to
                 it

***********************************************************************/

void store_string(char **string_ptr,char *string)
{
  int len;

  /* Get the length of the string */
  len = strlen(string);

  /* Allocate memory for the string */
  *string_ptr = (char *) malloc(sizeof(char)*(len + 1));

  /* Store the string */
  strcpy(*string_ptr,string);
}
/***********************************************************************

c_get_paths  -  Read in the paths for the common PDBsum and PDB files

***********************************************************************/

void c_get_paths(struct c_file_data *file_name_ptr,int run_64)
{
  char file_name[C_FILENAME_LEN];

  int error_status, idir;
  static int exit_on_error = FALSE;

  /* PDB directory pathnames */
  for (idir = 0; idir < NPDB_DIR; idir++)
    {
      error_status = get_param(PDB_DIR_TEMPLATE[idir],file_name,
                               C_FILENAME_LEN,exit_on_error);
      store_string(&(file_name_ptr->pdb_template[idir]),file_name);
    }

  /* PDBsum directory pathnames */
  for (idir = 0; idir < NPDBSUM_DIR; idir++)
    {
      /* Get the PDBsum directory template */
      error_status = get_param(PDBSUM_DIR_TEMPLATE[idir],file_name,
                               C_FILENAME_LEN,exit_on_error);
      store_string(&(file_name_ptr->pdbsum_template[idir]),file_name);

      /* Get the name of the corresponding data directory */
      error_status = get_param(PDBSUM_DATA_DIR[idir],file_name,
                               C_FILENAME_LEN,exit_on_error);
      store_string(&(file_name_ptr->pdbsum_data_dir[idir]),file_name);
    }

  /* PDBsum generate executable directories */
  if (run_64 == TRUE)
    error_status = get_param("PDBSUM_EXE_64_DIR",file_name,
                             C_FILENAME_LEN,exit_on_error);
  else
    error_status = get_param("PDBSUM_EXE_DIR",file_name,
                             C_FILENAME_LEN,exit_on_error);
  store_string(&(file_name_ptr->pdbsum_exe_dir),file_name);

  /* PDBsum generate executable directories */
  if (run_64 == TRUE)
    error_status = get_param("PDBSUM_GEN_64_DIR",file_name,
                             C_FILENAME_LEN,exit_on_error);
  else
    error_status = get_param("PDBSUM_GEN_DIR",file_name,
                             C_FILENAME_LEN,exit_on_error);
  store_string(&(file_name_ptr->pdbsum_gen_dir),file_name);

  /* Executables */
  for (idir = 0; idir < NEXE_DIR; idir++)
    {
      error_status = get_param(EXE_NAME[idir],file_name,
                               C_FILENAME_LEN,exit_on_error);
      store_string(&(file_name_ptr->exe_name[idir]),file_name);
    }

  /* Other file and pathnames */
  for (idir = 0; idir < NOTHER_DIR; idir++)
    {
      error_status = get_param(OTHER_DIR_NAME[idir],file_name,
                               C_FILENAME_LEN,exit_on_error);
      store_string(&(file_name_ptr->other_dir[idir]),file_name);
    }
}
/***********************************************************************

check_file  -  Check whether the given file exists

***********************************************************************/

int check_file(char *file_name)
{
  int have_file;

  FILE *file_ptr;

  /* Initialise */
  have_file = FALSE;

  /* Open the file to see if it exists */
  if ((file_ptr = fopen(file_name,"r")) !=  NULL)
    {
      /* Have the file */
      have_file = TRUE;

      /* Close it */
      fclose(file_ptr);
    }

  /* Return whether file exists */
  return(have_file);
}
/***********************************************************************

fix_upper  -  Search for certain strings that must be set in upper-case

***********************************************************************/

void fix_upper(char *name)
{
#define NCAPWORDS 25

  char ch;

  char *string_ptr;

  int ipos, jpos, len, nchars, word;
  int check_word[NCAPWORDS], done, is_word;

  static char *lowword[] = {
    "pdb", "dna", "rna", "hiv", "nag", "nad", "nadp", "nadph", "trna",
    "mrna", "e.c", "ec", "fab",
    "i", "iii", "iii", "iv", "v", "vi", "vii", "viii", "ix", "x",
    "xi", "xii"
  };

  static char *capword[] = {
    "PDB", "DNA", "RNA", "HIV", "NAG", "NAD", "NADP", "NADPH", "tRNA",
    "mRNA", "E.C", "EC", "FAB",
    "I", "III", "III", "IV", "V", "VI", "VII", "VIII", "IX", "X",
    "XI", "XII"
  };

  /* Initialise variables */
  done = FALSE;
  nchars = strlen(name);
  for (word = 0; word < NCAPWORDS; word++)
    check_word[word] = TRUE;

  /* Loop until all cap-words have been replaced */
  while (done == FALSE)
    {
      /* Assume that we are done */
      done = TRUE;

      /* Loop through all the key words */
      for (word = 0; word < NCAPWORDS; word++)
        {
          /* If this word to be checked, then do so */
          if (check_word[word] == TRUE)
            {
              /* Search for presence of this word in the string */
              string_ptr = strstr(name,lowword[word]);
              if (string_ptr != NULL)
                {
                  /* Assume this is a complete word */
                  is_word = TRUE;

                  /* Get the string's start-position */
                  ipos = string_ptr - name;
                  len = strlen(lowword[word]);
                  jpos = len + ipos;

                  /* Check character to left of match to test for word
                     separator */
                  if (ipos > 0)
                    {
                      /* Get character to left */
                      ch = name[ipos - 1];

                      /* If character is a letter, then match does
                         not represent a word */
                      if ((ch >= 'a' && ch <= 'z') ||
                          (ch >= 'A' && ch <= 'Z'))
                        is_word = FALSE;
                    }

                  /* If still possibly a word, check the character
                     to the right of the word */
                  if (is_word == TRUE && jpos < nchars)
                    {
                      /* Get character to right */
                      ch = name[jpos];

                      /* If character is a letter, then match does
                         not represent a word */
                      if ((ch >= 'a' && ch <= 'z') ||
                          (ch >= 'A' && ch <= 'Z'))
                        is_word = FALSE;
                    }

                  /* If have a capital-word, then splice into string */
                  if (is_word == TRUE)
                    {
                      /* Splice into string */
                      strncpy(string_ptr,capword[word],len);

                      /* Set flag that we're not done yet */
                      done = FALSE;
                    }
                }

              /* If word not present, then don't need to check for
                 it again */
              else
                check_word[word] = FALSE;
            }
        }
    }
}
/***********************************************************************

form_filename  -  Use the filename template and entered PDB code to
                  determine full name of file

***********************************************************************/

char *form_filename(char *pdb_code,char *filename_template)
{
  char ch, copy_char[2], file_name[C_FILENAME_LEN];

  char *final_name = &null;
  char *string_ptr;

  int i, ichar, len, len_code;
  int incode;

  /* Initialise variables */
  file_name[0] = '\0';
  incode = FALSE;

  /* Get length of template and PDB code */
  len = strlen(filename_template);
  len_code = strlen(pdb_code);

  /* If not too long, form file name */
  if (len < C_FILENAME_LEN)
    {
      /* Find any strings between square brackets in this template */
      for (i = 0; i < len; i++)
        {
          /* Get this character */
          ch = filename_template[i];
          copy_char[0]= '\0';

          /* If this an open-bracket, set flag */
          if (ch == '[')
            incode = TRUE;

          /* If close-bracket, unset flag */
          else if (ch == ']')
            incode = FALSE;
          else
            {
              /* If in code, transfer appropriate part of PDB code to
                 file name */
              if (incode == TRUE)
                {
                  /* See if this character is a template character */
                  string_ptr = strchr(TEMPLATE_CHARS,ch);

                  /* If present, transfer corresponding character from
                     PDB code */
                  if (string_ptr != NULL)
                    {
                      /* Get this character's array position */
                      ichar = string_ptr - TEMPLATE_CHARS;

                      /* If valid, then save character as the copy
                         character */
                      if (ichar > -1 && ichar < len_code)
                        copy_char[0] = pdb_code[ichar];
                    }
                }

              /* Otherwise, just transfer character from template */
              else
                copy_char[0] = ch;

              /* If have a character to copy across, perform the transfer */
              if (copy_char[0] != '\0')
                {
                  /* Append null and tack onto end of file name */
                  copy_char[1] = '\0';
                  strcat(file_name,copy_char);
                }
            }
        }

      /* Store the final file name */
      store_string(&final_name,file_name);
    }

  /* Otherwise, show error message */
  else
    printf("*** ERROR. Filename template too long!\n"
           "           Template: %s\n",filename_template);

  /* Return name of file */
  return(final_name);
}
/***********************************************************************

find_pdbsum_dir  -  Locate the PDBsum directory for this PDB entry

***********************************************************************/

int find_pdbsum_dir(char *pdb_code,char *pdbsum_template[NPDBSUM_DIR],
                    int pdb_type,char *dir_name)
{
  char file_name[C_FILENAME_LEN];
  char *fname;

  int dir_type, end, len, start, type;
  int error;

  struct stat st;

  FILE *file_ptr;

  /* Initialise */
  strcpy(dir_name,"./");
  dir_type = NOT_FOUND;

  /* Determine the directory type */
  if (pdb_type == STANDARD_PDB)
    {
      start = PDBSUM_CURRENT;
      end = PDBSUM_CURRENT + 1;
    }
  else if (pdb_type == FIND_FIRST)
    {
      start = PDBSUM_CURRENT;
      end = NPDBSUM_DIR;
    }
  else
    {
      start = PDBSUM_UPLOAD;
      end = NPDBSUM_DIR;
    }

  /* Loop over the PDBsum directories to try */
  for (type = start; type < end && dir_type == NOT_FOUND; type++)
    {
      /* Form the name of the PDBsum directory */
      fname = form_filename(pdb_code,pdbsum_template[type]);
      strcpy(file_name,fname);

      /* If searching for PDBsum dir, see if this one exists */
      if (pdb_type == FIND_FIRST)
        error = stat(file_name,&st);
      else
        error = FALSE;

      /* If directory exists, save its name */
      if (error == FALSE)
        {
          /* Save the directory name */
          strcpy(dir_name,fname);

          /* Set flag */
          dir_type = type;
        }
    }

  /* Append directory separator if required */
  if (dir_type != NOT_FOUND)
    {
      len = strlen(dir_name);
      if (dir_name[len - 1] != '/')
        strcat(dir_name,"/");
    }

  /* Return the file type */
  return(dir_type);
}
/***********************************************************************

free_tokens - Free up the memory used for storing token strings

**********************************************************************/  

void free_tokens(char **token,int ntokens)
{
  int itoken;

  /* Loop over the tokens to free up the memory */
  for (itoken = 0; itoken < ntokens; itoken++)
    if (token[itoken][0] != '\0')
      free(token[itoken]);
}
/***********************************************************************

get_random_number  -  Routine for generating a random number from the
                      supplied seed. Returns a uniform random deviate
                      between 0.0 and 1.0. On the first call, the
                      seed must be set to any negative value to
                      initialise or reinitialise the sequence.
                      (Algorithm is that described for routine ran2
                      in Numerical Recipes).

***********************************************************************/

float get_random_number(int *iseed,int random_start)
{
  static int ir[97], iy, j;

  static float random_number, seed;

  static int first_pass = TRUE;
  static int m = 714025;
  static int ia = 1366;
  static int ic = 150889;
  static float rm = (float) (1.4005112E-6);

  FILE *file_ptr;

  /* Initialise the random number */
  random_number = 0.0;

  /* If this is the first pass, or sequence to be re-initialised, then
     pick up seed from file */
  if (*iseed < 0 || first_pass == TRUE)
    {
      /* Initialise the seed */
      seed = 0.0;

      /* If a completely random start is required, then read in the
         number written out to the file ranseed.dat last time this
         routine was run */
      if (random_start == TRUE)
        {
          /* Open the random seed file and read in seed */
          if ((file_ptr = fopen("ranseed.dat","r")) !=  NULL)
            {
              /* If file exists, then read in the seed */
              fscanf(file_ptr,"%f",&seed);

              /* Close the file */
              fclose(file_ptr);
            }

          /* Store integer version of the seed */
          *iseed = (int) seed;
        }

      /* Otherwise, use a standard starting point */
      else
        *iseed = 54813;

      /* Check that the seed is non-zero */
      if (*iseed == 0) 
        *iseed = 1;

      /* Generate random number from seed */
      *iseed = (ic - *iseed) % m;
      for (j = 0; j < 97; j++)
        {
          *iseed = (ia * (*iseed) + ic) % m;
          ir[j] = *iseed;
        }
      *iseed = (ia * (*iseed) + ic) % m;
      iy = *iseed;
    }

  /* Check for error */
  j = (97 * iy) / m;
  if (j >= 97 || j < 0)
    {
      printf("*** Random number error\n");
      exit(1);
    }

  /* Get the random number */
  iy = ir[j];
  random_number = iy * rm;
  *iseed = (ia * (*iseed) + ic) % m;
  ir[j] = *iseed;

  /* If this is the first call, use the generated random number to
     create a new seed and write to disk */
  if (first_pass == TRUE)
    {
      first_pass = FALSE;
      seed = iy * rm * 100000;

      /* Open the random seed file */
      if ((file_ptr = fopen("ranseed.dat","w")) !=  NULL)
        {
          /* If no file error, then write out the seed */
          fprintf(file_ptr,"%f\n",seed);

          /* Close the file */
          fclose(file_ptr);
        }
    }

  /* Return the random number generated */
  return(random_number);
}
/***********************************************************************

get_tokens - Get the space- or tab-delimited token from the given
             input line

**********************************************************************/  

int get_tokens(char *line,char **token,int max_token,char token_sep)
{
  char ch, last_ch, token_sep2;

  char *tmp, *string, *tmp_cpy;

  int intok, ipos, itoken, len, ntoken, tokenlen;
  int done, keep_token;

  /* Initialise variables */
  done = FALSE;
  keep_token = FALSE;
  last_ch = '\0';
  ntoken = 0;
  intok = 0;
  ipos = 0;
  if (token_sep == '[')
    {
      token_sep2 = ']';
      keep_token = TRUE;
    }
  else if (token_sep == ':')
    token_sep2 = ':';
  else
    token_sep2 = '\0';

  /* Get the length of the line */
  len = strlen(line);
  if (len == 0)
    return(ntoken);

  /* Allocate space for temporary string */
  tmp = (char *) malloc(sizeof(char)*(len + 1));
  tmp[0] = '\0';
  tmp_cpy = (char *) malloc(sizeof(char)*(len + 1));
  tmp_cpy[0] = '\0';

  /* Loop until end of line reached */
  while (done == FALSE && ipos < len + 1)
    {
      /* Get the next character in the string */
      ch = line[ipos];

      /* Check for end of line */
      if (ch == '\0' || ch == '\n')
        {
          done = TRUE;
          ch = token_sep;
        }

      /* Check for double-colon separator */
      if (token_sep == ':' && ch == token_sep && last_ch == token_sep)
        ch = '\0';

      /* If have a token separator, then take to be the end of the
         token if have already got something in it */
      if (ch == token_sep || ch == token_sep2)
        {
          /* If have a token, increment token count */
          if (tmp[0] != '\0' || token_sep != ' ')
            {
              /* If keeping tokens and this is a token-end, add it
                 to the current string */
              if (keep_token == TRUE && ch == token_sep2)
                {
                  tmp[intok] = ch;
                  intok++;
                }

              /* Terminate current token */
              tmp[intok] = '\0';

              /* Check for bounding quotes */
              if (tmp[0] == '"' &&
                  tmp[intok - 1] == '"')
                {
                  /* Strip out the quotes */
                  strcpy(tmp_cpy,tmp+1);
                  strcpy(tmp,tmp_cpy);
                  tmp[intok - 2] = '\0';
                }

              /* Save the current token */
              tokenlen = strlen(tmp);
              store_string(&(token[ntoken]),tmp);

              /* Increment token count */
              ntoken++;

              /* Reinitialise position within current token */
              intok = 0;
              tmp[intok] = '\0';

              /* If keeping tokens and this is a token-start, add it
                 at the start of the current string */
              if (keep_token == TRUE && ch == token_sep)
                {
                  tmp[intok] = ch;
                  intok++;
                }
            }

          /* If keeping separator and this is a start-separator,
             start next token with it */
          else if (keep_token == TRUE && ch == token_sep)
            {
              tmp[intok] = ch;
              intok++;
            }
        }

      /* Otherwise, add character to token-string provided not
         too long */
      else if (ch != '\0')
        {
          if (ntoken < max_token && intok < len)
            {
              tmp[intok] = ch;
              intok++;
            }
        }

      /* Save the current character */
      last_ch = ch;

      /* Go to the next character in the line */
      ipos++;
    }

  /* Free up memory taken by temporary line */
  free(tmp);
  free(tmp_cpy);

  /* Return the number of tokens picked up */
  return(ntoken);
} 
/***********************************************************************

open_output_file  -  Open the output file

***********************************************************************/

FILE *open_output_file(char *out_name,int terminate)
{
  FILE *outfile_ptr;

  /* Open the output file, or write to stdout */
  if (out_name[0] == '\0')
    outfile_ptr = stdout;
  else
    outfile_ptr = fopen(out_name,"w");

  /* Warn or terminate */
  if (outfile_ptr == NULL)
    {
      if (terminate == TRUE)
        {
          printf("*** ERROR. Unable to open output file: %s\n",
                 out_name);
          printf("***        Program aborted\n");
          exit(-1);
        }
      else
        printf("*** Warning. Unable to open output file: %s\n",
               out_name);
    }

  /* Return the pointer to the output file */
  return(outfile_ptr);
}
/***********************************************************************

parse_uniprot_name  -  Parse the UniProt title to split up into its
                       components: protein name, species, etc

***********************************************************************/

void parse_uniprot_name(char *string,char *protein,char *species,
                        char *uniprot_acc,char *uniprot_id)
{
  char *tmp, *string_ptr, *start_string_ptr;

  int i, last_space, len;
  int done;

  /* Initialise variables */
  strcpy(protein,string);
  species[0] = '\0';
  uniprot_acc[0] = '\0';
  uniprot_id[0] = '\0';

  /* Get the length of the string */
  len = strlen(string);

  /* Allocate memory for temporary string */
  tmp = (char *) malloc(sizeof(char)*(len + 1));

  /* Identify the start of the codes */
  start_string_ptr = string_ptr = strstr(string,"sp|");
  if (string_ptr == NULL)
    start_string_ptr = string_ptr = strstr(string,"tr|");

  /* If found, then extract them */
  if (string_ptr != NULL)
    {
      /* Copy to temporary string */
      strcpy(tmp,string_ptr + 3);

      /* Get the position of the divider */
      string_ptr = strchr(tmp,'|');

      /* If found, we have the two codes, so can store */
      if (string_ptr != NULL)
        {
          /* Get the length of the accession code */
          len = string_ptr - tmp;

          /* Copy across */
          strncpy(uniprot_acc,tmp,len);
          uniprot_acc[len] = '\0';

          /* Get the UniProt id */
          string_ptr = strchr(tmp,' ');

          /* If end of id found, store the id */
          if (string_ptr != NULL)
            {
              /* Terminate the string */
              string_ptr[0] = '\0';

              /* Copy across */
              strcpy(uniprot_id,tmp + len + 1);
            }
        }

      /* Look for first blank after codes block */
      string_ptr = strchr(start_string_ptr,' ');

      /* If found, then transfer */
      if (string_ptr != NULL)
        {
          /* Save the start of the name string */
          start_string_ptr = string_ptr + 1;

          /* Copy protein name to temporary string */
          strcpy(tmp,start_string_ptr);

          /* Find the start of the species */
          string_ptr = strstr(tmp," OS=");

          /* If found, can save the protein name */
          if (string_ptr != NULL)
            string_ptr[0] = '\0';

          /* Save the protein name */
          strcpy(protein,tmp);
        }
    }

  /* Extract the species name */
  string_ptr = strstr(string," OS=");

  /* If found, then extract species name */
  if (string_ptr != NULL)
    {
      /* Copy to temporary area */
      strcpy(tmp,string_ptr + 4);

      /* Get the length of the current string */
      len = strlen(tmp);

      /* Loop through the characters to find the last space before the
         next item */
      last_space = 0;
      done = FALSE;
      for (i = 0; i < len && done == FALSE; i++)
        {
          /* Check if this is a space or a new field */
          if (tmp[i] == ' ')
            last_space = i;
          else if (tmp[i] == '=')
            done = TRUE;
        }

      /* Terminate the species string and copy */
      tmp[last_space] = '\0';
      strcpy(species,tmp);
    }

  /* Free up memory taken by the temporary string */
  free(tmp);
}
/***********************************************************************

perform_splice  -  Splice out given string from line and replace with
                   new string

***********************************************************************/

int perform_splice(char *line,char *string,char *new_string)
{
  char *tmp;

  int end_pos, ipos, len, len_new, len_string;
  int done;

  /* Initialise variables */
  done = FALSE;
  end_pos = 0;

  /* Get the length of the input line */
  len = strlen(line);
  len_string = strlen(string);
  len_new = strlen(new_string);
/* debug 
  printf("Line: %s\n",line);
  printf("Old string: %s\n",string);
  printf("New string: %s\n",new_string);
  fflush(stdout);
 debug */

  /* If any string empty, then return */
  if (len == 0 || len_string == 0)
    return(len);

  /* Loop through the line searching for the replace-string */
  for (ipos = 0; ipos < len - len_string + 1 && done == FALSE; ipos++)
    {
      /* Check the string at the current position */
      if (!strncmp(line+ipos,string,len_string))
        {
          /* Create a temporary string to hold the spliced string */
          tmp = (char *) malloc(sizeof(char)*(len + len_new + 1));

          /* Initialise the temporary string */
          tmp[0] = '\0';

          /* Transfer first part of line */
          strncpy(tmp,line,ipos);
          tmp[ipos] = '\0';

          /* Splice in the replacement path */
          strcat(tmp,new_string);

          /* Store end position of splice */
          end_pos = ipos + len_new;

          /* Add on the remainder of the line */
          strcat(tmp,line+ipos+len_string);

          /* Transfer back across to the input line */
          strcpy(line,tmp);

          /* Free up memory taken by the temporary string */
          free(tmp);

          /* Set flag that splice done */
          done = TRUE;
        }
    }

  /* Return splice end position */
  return(end_pos);
}
/***********************************************************************

remove_blanks  -  Remove all blanks from the given string

***********************************************************************/

void remove_blanks(char *string)
{
  char *tmp;

  int i, len;
  int done;

  /* Initialise variables */
  done = FALSE;
  i = 0;

  /* Get the length of the string */
  len = strlen(string);

  /* Allocate memory for temporary string */
  tmp = (char *) malloc(sizeof(char)*(len + 1));

  /* Loop through the characters squeezing out the blanks */
  while (done == FALSE)
    {
      /* If current character is NULL, then we're done */
      if (string[i] == '\0')
        done = TRUE;

      /* If this is a blank, then shift string along */
      else if (string[i] == ' ')
        {
          strcpy(tmp,string + i + 1);
          strcpy(string + i,tmp);
        }

      /* Otherwise, increment position */
      else
        i++;
    }

  /* Free up memory taken by the temporary string */
  free(tmp);
}
/***********************************************************************

remove_leading_blanks  -  Remove leading blank spaces from the given
                          string

***********************************************************************/

void remove_leading_blanks(char *string)
{
  char *tmp;

  int ipos, fpos, len;

  /* If no leading blanks, then return */
  if (string[0] != ' ' && string[0] != '_')
    return;

  /* Get the length of the string */
  len = strlen(string);

  /* Allocate memory for temporary string */
  tmp = (char *) malloc(sizeof(char)*(len + 1));

  /* Initialise position of first non-blank character */
  fpos = 0;

  /* Loop through to first non-blank character */
  for (ipos = 0; ipos < len && fpos == 0; ipos++)
    if (string[ipos] != ' ' && string[ipos] != '_')
      fpos = ipos;

  /* Shift whole string across */
  if (fpos > 0)
    {
      strcpy(tmp,string+fpos);
      strcpy(string,tmp);
    }
  else
    string[0] = '\0';

  /* Free up memory taken by the temporary string */
  free(tmp);
}
/***********************************************************************

remove_marks  -  Remove all type-setting characters from the given
                 string

***********************************************************************/

void remove_marks(char *string)
{
  char *tmp;

  int ipos, len;

  /* Get the length of the string */
  len = strlen(string);
  if (len == 0)
    return;

  /* Allocate memory for temporary string */
  tmp = (char *) malloc(sizeof(char)*(len + 1));

  /* Initialise variables */
  ipos = 0;

  /* Loop until string has been processed */
  while (ipos < len)
    {
      /* If have a blank space at this position, remove */
      if (string[ipos] == '$' || string[ipos] == '/' ||
          string[ipos] == '\\')
        {
          /* Shift string along */
          strcpy(tmp,string+ipos+1);
          strcpy(string+ipos,tmp);
          len--;
        }

      /* Otherwise, move on to next position */
      else
        ipos++;
    }

  /* Free up memory taken by the temporary string */
  free(tmp);
}
/***********************************************************************

string_chop  -  Truncate the given string at the last non-blank
                character

***********************************************************************/

int string_chop(char *string)
{
  int iend, ipos, len;

  /* Get length of string */
  len = strlen(string);

  /* If empty, then return */
  if (len == 0)
    return(len);

  /* Remove endline, if there is one */
  if (string[len - 1] == 13 || string[len - 1] == '\n')
    {
      string[len - 1] = '\0';
      len--;
    }

  /* Initialise end character */
  iend = 0;

  /* Loop back through the string until we hit the first non-blank
     character */
  for (ipos = len - 1; ipos > -1 && iend == 0; ipos--)
    {
      /* If not a space at this position, then have end of string */
      if (string[ipos] != ' ')
        iend = ipos + 1;
    }

  /* Truncate the string and set its length */
  string[iend] = '\0';
  len = iend;

  /* Return the string length */
  return(len);
}
/***********************************************************************

string_truncate  -  Truncate the given string at the last non-blank
                    character, or at its maximum length

***********************************************************************/

int string_truncate(char *string,int max_length)
{
  char ch;

  int iend, ipos, len;

  /* Initialise variables */
  iend = -1;
  ch = string[0];
  for (ipos = 0; ipos < max_length && ch != '\0'; ipos++)
    {
      /* Get the current character and check for end of string */
      ch = string[ipos];
      if (ch == '\n' || ch == 13)
        ch = '\0';
      if (ch != '\0')
        {
          /* If not a space, then mark end of string */
          if (ch != ' ')
            iend = ipos;
        }
    }

  /* End the string after the last non-blank character */
  if (iend + 1 < max_length)
    {
      string[iend + 1] = '\0';
      len = iend + 1;
    }
  else
    {
      string[max_length - 1] = '\0';
      len = max_length - 1;
    }

  /* Return the string length */
  return(len);
}
/***********************************************************************

to_lower  -  Convert given character string to lower case

***********************************************************************/

void to_lower(char *search_string,int length)
{
  int ichar;

  int ipos;

  /* Loop through all the character positions, converting any letters
     to upper case */
  for (ipos = 0; ipos < length; ipos++)
    {
      /* Check the current character */
      ichar = search_string[ipos] - 'A';
      if (ichar >= 0 && ichar < 26)
        {
          ichar = ichar + 'a';
          search_string[ipos] = ichar;
        }
    }
}
/***********************************************************************

to_upper  -  Convert given character string to upper case

***********************************************************************/

void to_upper(char *search_string,int length)
{
  int ichar;

  int ipos;

  /* Loop through all the character positions, converting any letters
     to upper case */
  for (ipos = 0; ipos < length; ipos++)
    {
      /* Check the current character */
      ichar = search_string[ipos] - 'a';
      if (ichar >= 0 && ichar < 26)
        {
          ichar = ichar + 'A';
          search_string[ipos] = ichar;
        }
    }
}
/***********************************************************************

tokenize - Split the given string into tokens, replacing each instance
           of the token separator with a null and returning an array of
           string start-points for each token

           Note: this corrupts the original string

**********************************************************************/  

int tokenize(char *string,char **token,int max_tokens,char token_sep)
{
  char ch;

  char *string_ptr;

  int ipos, ntokens;
  int done;

  /* Initialise variables */
  done = FALSE;

  /* Set first token at rhe start of the string */
  token[0] = string;
  ntokens = 1;

  /* Loop until all tokens have been extracted */
  while (done == FALSE && ntokens < max_tokens)
    {
      /* Get location of the next token separator */
      string_ptr = strchr(string,token_sep);

      /* If have a token, null terminate and keep looking */
      if (string_ptr != NULL)
        {
          /* Shift start-pos of string and next token */
          token[ntokens] = string = string_ptr + 1;
          ntokens++;

          /* Null terminate current token */
          string_ptr[0] = '\0';
        }

      /* Otherwise, we're done */
      else
        done = TRUE;
    }

  /* Return number of tokens */
  return(ntokens);
}
/***********************************************************************

type_set  -  Remove typographical characters and convert into upper and
             lower case.

***********************************************************************/

void type_set(char *name,int cap_start)
{
  char ch, last_ch;

  int ifree, ipos, nchars;
  int capital, dollar, skip, slash, upcase, up_one;

  /* Initialise variables */
  capital = TRUE;
  last_ch = ' ';
  upcase = FALSE;
  up_one = FALSE;
  ifree = 0;
  nchars = strlen(name);

  /* Check whether the string contains a dollar after a slash */
  dollar = FALSE;
  slash = FALSE;
  for (ipos = 0; ipos < nchars; ipos++)
    {
      ch = name[ipos];
      if (ch == '/')
        slash = TRUE;
      if (ch == '$' && slash == TRUE)
        dollar = TRUE;
    }

  /* Loop over all the characters */
  for (ipos = 0; ipos < nchars; ipos++)
    {
      /* Retrieve this character */
      ch = name[ipos];
      skip = FALSE;

      /* If in lower case, convert to upper */
      if (ch >= 'a' && ch <= 'z')
        ch = ch - 'a' + 'A';

      /* Check if this is a letter */
      if (ch >= 'A' && ch <= 'Z')
        {
          /* If first character then leave it in upper-case */
          if (capital == TRUE || up_one == TRUE)
            {
              capital = FALSE;
              up_one = FALSE;
            }

          /* Otherwise, set in lower-case, if required */
          else if (upcase == FALSE && up_one == FALSE)
            ch = ch - 'A' + 'a';
        }

      /* Check if this is a number */
      else if (ch >= '0' && ch <= '9')
        {
          /* If first character then set whole word in lower-case (as
             possibly a PDB code) */
          if (capital == TRUE)
            {
              capital = FALSE;
              up_one = FALSE;
            }
        }

      /* Check if this is one of the typesetting characters */

      /* Upper-case the next letter only */
      else if (ch == '*')
        {
          up_one = TRUE;
          skip = TRUE;
        }

      /* Upper-case until dollar */
      else if (ch == '/' && dollar == TRUE)
        {
          upcase = TRUE;
          skip = TRUE;
        }

      /* End upper-case */
      else if (ch == '$')
        {
          upcase = FALSE;
          skip = TRUE;
        }

      /* After a full-stop, need to capitalise next character */
      else if (ch == '.')
        capital = TRUE;

      /* If in author-list, capitalise next letter after a comma, dash
         or apostrophe too */
      else if (cap_start == TRUE && (ch == ',' || ch == '-' ||
                                     ch == '\''))
        capital = TRUE;

      /* Save the current character */
      if (skip == FALSE)
        {
          name[ifree] = ch;
          last_ch = ch;
          ifree = ifree + 1;
        }
    }

  /* Terminate string */
  name[ifree] = '\0';
}
/***********************************************************************

write_ppm_file  -  Write out the given image array as a .ppm file

***********************************************************************/

void write_ppm_file(int *array,int ncols,int nrows,char *file_name,
                    char *plot_name)
{
  int colour, icol, iloc, irow, ntrip, rgb[3];

  FILE *file_ptr;

  /* Open the output file */
  file_ptr = open_output_file(file_name,FALSE);
  if (file_ptr == NULL)
    return;

  /* Write out header records */
  fprintf(file_ptr,"P3\n");
  fprintf(file_ptr,"# %s %s\n",plot_name,file_name);
  fprintf(file_ptr,"%d %d\n",ncols,nrows);
  fprintf(file_ptr,"255\n");

  /* Initialise count of RGB triplets on line */
  ntrip = 0;

  /* Loop through the rows making up the image */
  for (irow = nrows - 1; irow > -1; irow--)
    {
      /* Loop over the columns to be plotted */
      for (icol = 0; icol < ncols; icol++)
        {
          /* Get position of this pixel in the image array */
          iloc = irow * ncols + icol;

          /* Get the colour at this pixel position */
          colour = array[iloc];

          /* Convert colour into an RGB triplet */
          rgb[0] = colour / (256 * 256);
          colour = colour % (256 * 256);
          rgb[1] = colour / 256;
          rgb[2] = colour % 256;

          /* Write out this pixel */
          fprintf(file_ptr,"%d %d %d ",rgb[0],rgb[1],rgb[2]);

          /* Increment count of triplets on this line */
          ntrip++;

          /* If have too many triplets on this line, throw
             new line */
          if (ntrip > 4)
            {
              fprintf(file_ptr,"\n");
              ntrip = 0;
            }
        }
    }

  /* Close the .ppm image file */
  fclose(file_ptr);
}
/***********************************************************************

convert_to_gif  -  Convert .ppm file to a GIF image

***********************************************************************/

void convert_to_gif(char *convert_exe,char *file_name,char *out_name)
{
  char *command;

  int len;

  /* Get the length of the command line */
  len = strlen(convert_exe) + strlen(file_name) + strlen(out_name) + 50;

  /* Show name of image file to be converted */
  printf("Converting .ppm image to %s ...\n",out_name);

  /* Allocate memory for the command */
  command = (char *) malloc(sizeof(char)*(len + 1));

  /* Form the parameters to run convert */
  sprintf(command,"%s -trim %s %s  > convert.log",convert_exe,file_name,
          out_name);

  /* Run convert */
#ifdef TEST
/* debug */
  printf("COMMAND: %s\n",command);
/* debug */
#else
  system(command);
#endif

  /* Free up memory taken by the command */
  free(command);
}
