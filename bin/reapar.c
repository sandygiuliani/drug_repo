/***********************************************************************

  reapar.c - Program to pick up all the program parameters (paths and
             filenames) as given in the Perl-script format CATHPARAM file

************************************************************************

Date:-         18 Jan 1999

Written by:-   Roman Laskowski

----------------------------------------------------------------------*/

/* Datafiles
   ---------

Actual name      File name        Description
-----------      ---------        -----------
CATHPARAM        file_name        Input CATH parameter file giving all the
                                  directory paths for the files used by
                                  the CATH, PDBsum, etc. programs

*/

/* Function calling-tree
   ---------------------

main
   -> get_parameter  or  get_param  or  get_next_param
         -> get_file_parameters
               -> init_file_parameters
               -> interpret_tokens
               -> expand_token_value
                     -> locate_token
               -> store_token
         -> get_token_value

*/


#include "reapar.h"





/***********************************************************************

interpret_tokens  -  Interpret the tokens on the current input line

***********************************************************************/

void interpret_tokens(char input_line[LINELEN + 1],int len,
                      char token_name[STRING_LEN],
                      char token_value[STRING_LEN])
{
  char ch;
  int ipos, outpos, state;
  int done, in_name, savechar;

  /* Initialise variables */
  done = FALSE;
  in_name = TRUE;
  token_name[0] = '\0';
  token_value[0] = '\0';

  /* Initialise state for interpreting the token in the line */
  outpos = 0;
  state = 0;

  /* Loop through all the characters in the line */
  for (ipos = 0; ipos < len + 1 && done == FALSE; ipos++)
    {
      /* Get the character at this position */
      ch = input_line[ipos];
      savechar = FALSE;

      /* Check that haven't reached the end of the line */
      if (ch == '\n' || ch == '\0')
        done = TRUE;

      /* Otherwise, process according to current state */
      else
        {
          switch (state)
            {
              /* If looking for first non-blank character, check if
                 we have one */
            case 0 :

              /* If non-blank, then initialise and switch state */
              if (ch != ' ' && ch != '\t')
                {
                  /* Initialise for token or token-string */
                  outpos = 0;
                  savechar = TRUE;
                  state = 1;
                }

              break;

              /* Adding characters to current string */
            case 1 :

              /* If non-blank, then add to current string */
              if (ch != ' ' && ch != '\t')
                savechar = TRUE;

              /* Otherwise, have reached end of current token */
              else
                {
                  /* Terminate string and switch state */
                  ch = '\0';
                  if (token_value[0] != '\0')
                    done = TRUE;
                  else
                    state = 2;
                  savechar = TRUE;
                }

              break;

              /* Looking for equals sign */
            case 2 :

              /* If this is an equals sign, the switch state */
              if (ch == '=')
                {
                  state = 0;
                  in_name = FALSE;
                }

              break;

              /* Default */
            default :
              break;

            }

          /* If saving character, then do so */
          if (savechar == TRUE)
            {
              /* Check that not off end of string */
              if (outpos >= STRING_LEN - 1)
                {
                  outpos = STRING_LEN - 1;
                  ch = '\0';
                }

              /* Save in token name or token value */
              if (in_name == TRUE)
                token_name[outpos] = ch;
              else
                token_value[outpos] = ch;
              outpos++;
            }
        }
    }

  /* Close off the token-name or token value */
  if (in_name == TRUE)
    token_name[outpos] = '\0';
  else
    token_value[outpos] ='\0';
}
/***********************************************************************

locate_token  -  Locate the token and replace its name by its value

***********************************************************************/

int locate_token(char *token_name,char *token_value,
		 struct token *first_token_ptr)
{
  int found;

  struct token *token_ptr;

  /* Initialise variables */
  found = FALSE;
  token_value[0] = '\0';

  /* Start at the first token */
  token_ptr = first_token_ptr;

  /* Loop through all the tokens until get a match */
  while (token_ptr != NULL && found == FALSE)
    {
      /* Check if the token name matches */
      if (!strcmp(token_name,token_ptr->name))
        {
          /* Have a match, so retrieve the token value */
          strcpy(token_value,token_ptr->value);

          /* Set flag to end search here */
          found = TRUE;
        }

      /* Get pointer to next token */
      token_ptr = token_ptr->next_token_ptr;
    }

  /* Return whether token found */
  return(found);
}
/***********************************************************************

expand_token_value  -  Expand token value if it is a compound name
                       involving a prior token value

***********************************************************************/

void expand_token_value(char token_value[STRING_LEN],
                        struct token *first_token_ptr)
{
  char *string_ptr;
  char ch, token_name[STRING_LEN], temp_string[STRING_LEN];
  char converted_value[STRING_LEN];

  int start_pos, end_pos, ipos, len;

  /* Check whether expansion required */
  string_ptr = strstr(token_value,"$(");

  /* Expand each token, one by one */
  while (string_ptr != NULL)
    {
      /* Initialise start- and end-markers */
      start_pos = -1;
      end_pos = -1;
      string_ptr = NULL;

      /* Get current length of the value-string */
      len = strlen(token_value);

      /* Find the token-string to be substituted */
      for (ipos = 0; ipos < len && end_pos == -1; ipos++)
        {
          /* Get this character */
          ch = token_value[ipos];

          /* Check for start of substitution */
          if (ch == '$')
            start_pos = ipos;

          /* Check for end */
          else if (ch == ')' && start_pos != -1)
            end_pos = ipos;
        }

      /* If have valid start- and end-values, perform the substitution */
      if (start_pos > -1 && end_pos > -1)
        {
          len = end_pos - start_pos - 2;

          /* Splice out the token-name to be substituted */
          strncpy(token_name,token_value+start_pos + 2,len);
          token_name[len] = '\0';

          /* Locate the token and replace its name by its value */
          locate_token(token_name,converted_value,first_token_ptr);

          /* Splice together the new token value */
          strncpy(temp_string,token_value,start_pos);
          temp_string[start_pos] = '\0';
          strcat(temp_string,converted_value);
          strcat(temp_string,token_value+(start_pos + len + 3));

          /* Store the updated token value */
          strcpy(token_value,temp_string);

          /* Check whether there are any more dollar tokens left */
          string_ptr = strstr(token_value,"$(");
        }
    }
}
/***********************************************************************

store_token  -  Store the given token in the linked list

***********************************************************************/

void store_token(char token_name[STRING_LEN],
                 char token_value[STRING_LEN],
                 struct token **fst_token_ptr,
                 struct token **lst_token_ptr)
{
  int len;

  struct token *token_ptr, *first_token_ptr, *last_token_ptr;

  /* Initialise pointers */
  first_token_ptr = *fst_token_ptr;
  last_token_ptr = *lst_token_ptr;

  /* Allocate memory for token */
  token_ptr = (struct token *) malloc(sizeof(struct token));

  /* Check that memory allocation successful */
  if (token_ptr == NULL)
    {
      printf("*** ERROR. Can't allocate memory for struct token\n");
      exit (1);
    }

  /* Add link from previous token to the current one */
  if (last_token_ptr != NULL)
    last_token_ptr->next_token_ptr = token_ptr;
  else
    first_token_ptr = token_ptr;
  last_token_ptr = token_ptr;

  /* Get length of name token and allocate memory for the string */
  len = strlen(token_name) + 1;
  token_ptr->name = (char *) malloc(sizeof(char)*(len));

  /* Store the string */
  strcpy(token_ptr->name,token_name);

  /* Repeat for value */
  len = strlen(token_value) + 1;
  token_ptr->value = (char *) malloc(sizeof(char)*(len));

  /* Store the string */
  strcpy(token_ptr->value,token_value);

  /* Initialise pointer to next token */
  token_ptr->next_token_ptr = NULL;

  /* Save token pointers */
  *fst_token_ptr = first_token_ptr;
  *lst_token_ptr = last_token_ptr;
}
/***********************************************************************

get_file_parameters  -  Read in the parameters from the CATHPARAM file

***********************************************************************/

int get_file_parameters(struct token **fst_token_ptr,int exit_on_error)
{
  char cathparam_name[FILENAME_LEN];
  char input_line[LINELEN + 1];
  char token_name[STRING_LEN], token_value[STRING_LEN];

  char *path;

  static char *param_name[] = {
    "CATHPARAM",
    "../../param/CATHPARAM",
    "/ebi/research/thornton/www/databases/cgi-bin/param/CATHPARAM",
    "./param/CATHPARAM",
    "../cath/param/CATHPARAM",
    "../../cath/param/CATHPARAM",
    "/homes/roman/data/pdbsum/CATHPARAM",
    "/nfs/httpd/cgi-bin/cath/param/CATHPARAM",
    "/usr/local/httpd/cgi-bin/cath/param/CATHPARAM",
    "\0"
  };

  int error_status;
  int done, iparam, len, line, ntoken;

  struct token *first_token_ptr, *last_token_ptr;

  FILE *file_ptr;

  /* Initialise variables */
  done = FALSE;
  error_status = OK;
  line = 0;
  ntoken = 0;
  file_ptr = NULL;

  /* Initialise the token pointers */
  *fst_token_ptr = NULL;
  first_token_ptr = NULL;
  last_token_ptr = NULL;

  /* Check if we have the name of the parameter file in the environment
     variable CATHPARAM */
  path = getenv("CATHPARAM");

  if (path != NULL)
    {
      /* Copy across as the file name */
      strcpy(cathparam_name,path);

      /* Try opening this file */
      file_ptr = fopen(cathparam_name,"r");

      /* If file-pointer not null, then have hit a valid file */
      if (file_ptr != NULL)
        done = TRUE;
    }

  /* If don't have the CATHPARAM file, search through list of default
     file names until we reach one that exists */
  for (iparam = 0; done == FALSE; iparam++)
    {
      /* Check that not the end of the filenames */
      if (param_name[iparam][0] != '\0')
        {
          /* Copy across as the file name */
          strcpy(cathparam_name,param_name[iparam]);

          /* Try opening this file */
          file_ptr = fopen(cathparam_name,"r");

          /* If file-pointer not null, then have hit a valid file */
          if (file_ptr != NULL)
            done = TRUE;
        }

      /* If have finished all possibles, then nothing more to try */
      else
        {
          /* If error to be considered fatal, then end */
          if (exit_on_error == TRUE)
            {
              printf("*** ERROR. CATHPARAM file not found\n");
              exit(1);
            }

          /* Otherwise, just show warning */
          else
            {
              printf("*** Warning. CATHPARAM file not found\n");
              error_status = FILE_MISSING;
            }

          /* Set flag that done */
          done = TRUE;
        }
    }

  /* If OK to continue, read through the file */
  if (error_status == OK)
    {
      /* Initialise flag */
      done = FALSE;

      /* Loop while reading in records from the file */
      while (fgets(input_line,LINELEN,file_ptr) != NULL && done == FALSE)
        {
          /* Get the length of the line read in */
          len = strlen(input_line);

          /* Ignore if line is blank or is a comment line */
          if (len > 1 && input_line[0] != '#')
            {
              /* Interpret and expand the tokens on this line */
              interpret_tokens(input_line,len,token_name,token_value);

              /* Expand token value if it is a compound name involving a
                 prior token value */
              expand_token_value(token_value,first_token_ptr);

              /* Store the current token */
              store_token(token_name,token_value,&first_token_ptr,
                          &last_token_ptr);
            }

          /* Check that not the end of data */
          else if (!strncmp(input_line,"#---------",10))
            done = TRUE;

          /* Increment line-count */
          line++;
        }
  
      /* Close the file */
      fclose(file_ptr);
    }

  /* Return the pointer to the first of the tokens read in */
  *fst_token_ptr = first_token_ptr;

  /* Return error status */
  return(error_status);
}
/***********************************************************************

get_parameter  -  Get the required file parameter

***********************************************************************/

void get_parameter(char *token_name,char *file_name,int filename_len)
{
  char token_value[STRING_LEN];

  int len;
  static int exit_on_error = TRUE;

  static int first_time = TRUE;

  static struct token *first_token_ptr;

  /* Initialise variables */
  file_name[0] = '\0';

  /* If this is the first call, then read in all the parameters from
     the CATHPARAM file */
  if (first_time == TRUE)
    {
      /* Interpret the command-line parameters */
      get_file_parameters(&first_token_ptr,exit_on_error);

      /* If no tokens found, then return error */
      if (first_token_ptr == NULL)
        {
          printf("*** ERROR. No tokens found in CATHPARAM file\n");
          return;
        }

      /* Reset the flag */
      first_time = FALSE;
    }

  /* Get the value corresponding to the given token */
  locate_token(token_name,token_value,first_token_ptr);

  /* Check that the length of the string not too long */
  len = strlen(token_value);
  if (len > filename_len - 1)
    {
      printf("*** ERROR. Filename from CATHPARAM file too long\n");
      exit(1);
    }

  /* Transfer token value across */
  strcpy(file_name,token_value);
}
/***********************************************************************

get_param  -  Get the required file parameter (without exit on error)

***********************************************************************/

int get_param(char *token_name,char *file_name,int filename_len,
              int exit_on_error)
{
  char token_value[STRING_LEN];

  int error_status, len;
  int found;

  static int first_time = TRUE;

  static struct token *first_token_ptr;

  /* Initialise variables */
  error_status = OK;
  file_name[0] = '\0';

  /* If this is the first call, then read in all the parameters from
     the CATHPARAM file */
  if (first_time == TRUE)
    {
      /* Interpret the command-line parameters */
      error_status = get_file_parameters(&first_token_ptr,exit_on_error);

      /* If no tokens found, then return error */
      if (error_status == OK && first_token_ptr == NULL)
        {
          /* If error to be considered fatal, then end */
          if (exit_on_error == TRUE)
            {
              printf("*** ERROR. No tokens found in CATHPARAM file\n");
              exit(1);
            }

          /* Otherwise, just show warning */
          else
            {
              printf("*** Warning. No tokens found in CATHPARAM file\n");
              error_status = NO_TOKENS;
            }

          /* Return */
          return(error_status);
        }

      /* Reset the flag */
      first_time = FALSE;
    }

  /* If OK to proceed, then do so */
  if (error_status == OK)
    {
      /* Get the value corresponding to the given token */
      found = locate_token(token_name,token_value,first_token_ptr);

      /* If token found, then save the value */
      if (found == TRUE)
	{
	  /* Check that the length of the string not too long */
	  len = strlen(token_value);
	  if (len > filename_len - 1)
	    {
	      if (exit_on_error == TRUE)
		{
		  printf("*** ERROR. Filename from CATHPARAM file "
			 "too long\n");
		  exit(1);
		}

	      /* Otherwise, just show warning */
	      else
		{
		  printf("*** Warning. Filename from CATHPARAM file "
			 "too long\n");
		  error_status = NAME_TOO_LONG;
		}
	    }

	  /* Transfer token value across */
	  else
	    strcpy(file_name,token_value);
	}

      /* If token not found, then show warning */
      else
	{
	  printf("*** Token not found in CATHPARAM file: %s\n",token_name);
	  error_status = TOKEN_MISSING;
	}
    }

  /* Return error status */
  return(error_status);
}
/***********************************************************************

get_next_param  -  Get the next file parameter from the CATHPARAM file

***********************************************************************/

int get_next_param(char *token_name,char *file_name,int filename_len,
                   int exit_on_error)
{
  int error_status, len1, len2;

  static int first_time = TRUE;

  struct token *token_ptr;
  static struct token *first_token_ptr, *next_token_ptr;

  /* Initialise variables */
  error_status = OK;
  file_name[0] = '\0';

  /* If this is the first call, then read in all the parameters from
     the CATHPARAM file */
  if (first_time == TRUE)
    {
      /* Interpret the command-line parameters */
      error_status = get_file_parameters(&first_token_ptr,exit_on_error);

      /* If no tokens found, then return error */
      if (error_status == OK && first_token_ptr == NULL)
        {
          /* If error to be considered fatal, then end */
          if (exit_on_error == TRUE)
            {
              printf("*** ERROR. No tokens found in CATHPARAM file\n");
              exit(1);
            }

          /* Otherwise, just show warning */
          else
            {
              printf("*** Warning. No tokens found in CATHPARAM file\n");
              error_status = NO_TOKENS;
            }

          /* Return */
          return(error_status);
        }

      /* Initialise next token pointer */
      next_token_ptr = first_token_ptr;

      /* Reset the flag */
      first_time = FALSE;
    }

  /* If OK to proceed, then do so */
  if (error_status == OK)
    {
      /* Get the next token */
      token_ptr = next_token_ptr;

      /* If valid, then retrieve token name and value */
      if (token_ptr != NULL)
        {
          /* Check that the lengths of the strings are not too long */
          len1 = strlen(token_ptr->name);
          len2 = strlen(token_ptr->value);
          if (len1 > filename_len - 1 || len2 > filename_len - 1)
            {
              /* Show error/warning message */
              if (exit_on_error == TRUE)
                printf("*** ERROR. Name from CATHPARAM file too long:\n");
              else
                printf("*** Warning. Name from CATHPARAM file too long\n");
              if (len1 > filename_len - 1)
                printf("    %s\n",token_ptr->name);
              if (len2 > filename_len - 1)
                printf("    %s\n",token_ptr->value);

              /* Exit/return */
              if (exit_on_error == TRUE)
                exit(1);
              else
                error_status = NAME_TOO_LONG;
            }

          /* Store token name and value */
          else
            {
              strcpy(token_name,token_ptr->name);
              strcpy(file_name,token_ptr->value);
            }

          /* Store next token */
          next_token_ptr = token_ptr->next_token_ptr;
        }

      /* Otherwise, return flag indicating that no more tokens */
      else
        error_status = NO_MORE_TOKENS;
    }

  /* Return error status */
  return(error_status);
}

/***********************************************************************

                             M   A   I   N

***********************************************************************/

void main_test(void)
{
  char file_name[FILENAME_LEN];

  printf("Looking for parameter PDB_TEMPLATE\n");

  /* Get the required file parameter */
  get_parameter("PDB_TEMPLATE",file_name,FILENAME_LEN);

  /* Write out answer */
  printf("Parameter PDB_TEMPLATE:       [%s]\n",file_name);

}

