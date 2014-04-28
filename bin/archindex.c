/***********************************************************************

  archindex.c - Program to generate ArchSchema index files

************************************************************************

Date:-         1 Dec 2009
Dir update:-   6 Dec 2012

Written by:-   Roman Laskowski

----------------------------------------------------------------------*/

/* Datafiles
   ---------

Actual name           Description
-----------           -----------
archindex.dat         Output .dat file of structures relating to a given
                      Pfam domain
archstruc.uniprot.joined.by.accession
                      UniProt id / architecture / name, etc / 3D coverage
cathstruc.uniprot.joined.by.accession
                      as above, but for CATH domain data
domseqs.struc         Pfam id / UniProt ids - No longer used
eheaders.txt          List of E.C. numbers and their descriptions
speclist.out          Input file giving species suffix/ species name
ssg.txt               Input file giving SSG annotations

Compilation
-----------

cc -c archindex.c
cc -c common.c
cc -c reapar.c
cc -o archindex archindex.o common.o reapar.o -lm -lz

Optimization flags:-  cc -O2 -funroll-loops -c archindex.c

For verbose mode, compile with -DVERBOSE


/* Function calling-tree
   ---------------------

main
   -> get_command_arguments
         -> write_error_message
         -> tokenize (in common.c)
         -> find_pfam
         -> store_string (in common.c)
         -> append_string
         -> create_pfam_record
               -> write_error_message
         -> create_filterseq_record
   -> c_get_paths (in common.c)
         -> get_param (in reapar.c)
         -> store_string (in common.c)

   C O M P I L E   D A T A

   -> get_arch_data
         -> tokenize (in common.c)
         -> create_arch_record
               -> write_error_message
   -> sort_arch_records
         -> arch_sort
   -> write_offsets_file
         -> open_output_file (in common.c)
   -> write_offsets_index
         -> open_output_file (in common.c)
   -> transfer_to_pfam
         -> create_pfam_record
               -> write_error_message
         -> create_offset_record
               -> write_error_message
   -> get_offsets
   -> write_pfam_index
         -> open_output_file (in common.c)
   -> index_to_index
         -> open_output_file

   S E A R C H   D A T A

   -> sort_pfam_records
         -> pfam_sort
   -> perform_search
         -> uniprot_search
               -> find_uniprot
                     -> read_uniprot_idx
                           -> string_truncate
                     -> read_uniprot_dat
                           -> fseek
                           -> string_truncate (in common.c)
                           -> create_offset_record
                                 -> write_error_message
                     -> get_seq_architecture
                           -> fseek
                           -> string_truncate (in common.c)
                           -> create_hit_record
                                 -> write_error_message
                     -> select_hit
                     -> split_architecture
                           -> create_pfam_record
                                 -> write_error_message
                     -> sort_pfam_records
                           -> pfam_sort
         -> pfam_search
               -> read_archpfam_idx
                     -> string_truncate (in common.c)
               -> read_archpfam_dat
                     -> fseek
                     -> string_truncate (in common.c)
                     -> create_offset_record
                           -> write_error_message
               -> get_architectures
                     -> fseek
                     -> string_truncate (in common.c)
                     -> tokenize (in common.c)
                     -> check_selection
                     -> create_hit_record
                           -> write_error_message
                           -> store_string (in common.c)
                     -> create_species_record
                           -> write_error_message
                           -> store_string (in common.c)
               -> sort_species_records
                     -> species_sort
               -> get_species_names
               -> choose_parent
               -> sort_hit_records
                     -> hit_sort
               -> extract_architectures
                     -> create_archit_record
                     -> split_architecture
                           -> tokenize (in common.c)
                           -> create_pfam_record
                                 -> write_error_message
               -> find_architecture
               -> collapse_pfamb_domains
                     -> replace_pfamb_domains
                           -> tokenize (in common.c)
                           -> store_string (in common.c)
                     -> sort_archit_records
                           -> archit_sort
                     -> combine_architectures
               -> compare_to_parent
                     -> score_architectures
                           -> extract_domains
                           -> score_extra_domains
                           -> create_alignment_arrays
                           -> fill_array
                                 -> calc_similarity
                           -> nw_sweep
                           -> find_best_path
                           -> calc_score
                                 -> calc_similarity
                           -> free_alignment_arrays
               -> delete_architectures
               -> read_ssg_data
               -> get_enzyme_data
                     -> remove_partials
                           -> tokenize (in common.c)
                     -> find_enzyme
                     -> create_enzyme_record
                           -> store_string (in common.c)
                           -> convert_ec
                     -> create_enzlist_record
               -> prune_sequences
               -> extract_pfam_domains
                     -> split_architecture
                           -> tokenize (in common.c)
                           -> create_pfam_record
                                 -> write_error_message
               -> sort_pfam_records
                     -> pfam_sort
               -> delete_pfam_records
               -> get_pfam_names
               -> assign_pfam_colours
               -> read_ssg_data
                     -> annotate_sequence
               -> write_hit_summary
               -> write_dat_summary
                     -> tokenize (in common.c)
                     -> find_pfam
               -> sort_by_score
                     -> score_sort
               -> write_hits
                     -> tokenize (in common.c)
                     -> find_pfam
               -> write_species
               -> write_pfam
               -> get_enzyme_names
                     -> sort_enzyme_records
                           -> enzyme_sort1
               -> sort_enzyme_records
                     -> enzyme_sort1
               -> write_enzymes
               -> calc_arch_dists
               -> get_connectivities

*/

#include "common.h"
#include "archindex.h"

/* Prototypes
   ---------- */

/* In common.c */
void c_get_paths(struct c_file_data *file_name_ptr,int run_64);
int find_pdbsum_dir(char *pdb_code,char *pdbsum_template[NPDBSUM_DIR],
                    int pdb_type,char *dir_name);
char *form_filename(char *pdb_code,char *filename_template);
FILE *open_output_file(char *out_name,int terminate);
void store_string(char **string_ptr,char *string);
int string_truncate(char *string,int max_length);
int tokenize(char *string,char **token,int max_tokens,char token_sep);
void to_upper(char *search_string,int length);

/* Function prototypes from reapar.c */
int get_param(char *token_name,char *file_name,int filename_len,
              int exit_on_error);

/***********************************************************************

append_string  -  Append given string to existing string, allocating
                  memory for the expanded string, storing pointer to it

***********************************************************************/

void append_string(char **string_ptr,char *string)
{
  char *new;

  int len;

  /* Get the total length of the two strings */
  len = strlen(*string_ptr) + strlen(string);

  /* Allocate memory for the new string */
  new = (char *) malloc(sizeof(char)*(len + 1));

  /* Append second string onto first */
  strcpy(new,*string_ptr);
  strcat(new,string);

  /* Free up memory taken by old string */
  free(*string_ptr);

  /* Return the new atring pointer */
  *string_ptr = new;
}
/***********************************************************************

find_pfam  -  Find the given pfam entry

***********************************************************************/

struct pfam *find_pfam(struct pfam *first_pfam_ptr,char *pfam_id)
{
  struct pfam *pfam_ptr, *found_pfam_ptr;

  /* Initialise variables */
  found_pfam_ptr = NULL;

  /* Get pointer to the first pfam */
  pfam_ptr = first_pfam_ptr;

  /* Loop over all the pfams to find the one of interest */
  while (pfam_ptr != NULL && found_pfam_ptr == NULL)
    {
      /* If this is the right pfam, store pointer */
      if (!strcmp(pfam_id,pfam_ptr->pfam_id))
        found_pfam_ptr = pfam_ptr;

      /* Get pointer to the next pfam */
      pfam_ptr = pfam_ptr->next_pfam_ptr;
    }

  /* Return the pfaming pfam pointer */
  return(found_pfam_ptr);
}
/***********************************************************************

write_error_message  -  Have error, so write error type and message to
                        the output

***********************************************************************/

void write_error_message(int type,char *message,int output_type)
{
  /* Write error message depending on the output file type */
  if (output_type == ARCHSCHEMA_OUTPUT)
    {
      fprintf(outfile_ptr,":ERROR\n");
      fprintf(outfile_ptr,"%s\t%s\n",ERROR_TYPE[type],message);
    }
  else
    {
      fprintf(outfile_ptr,"ERROR_TYPE %s\n",ERROR_TYPE[type]);
      fprintf(outfile_ptr,"ERROR_MESSAGE %s\n",message);
      fprintf(outfile_ptr,"NHITS 0\n");
    }
  fflush(outfile_ptr);
  fclose(outfile_ptr);
  /*
  fprintf(outfile_ptr,"<void method=\"add\">\n");
  fprintf(outfile_ptr,"<int>%d</int>\n",RUN_ERROR);
  fprintf(outfile_ptr,"</void>\n");
  fprintf(outfile_ptr,"<void method=\"add\">\n");
  fprintf(outfile_ptr,"<string>%s</string>\n",ERROR_TYPE[type]);
  fprintf(outfile_ptr,"</void>\n");
  fprintf(outfile_ptr,"<void method=\"add\">\n");
  fprintf(outfile_ptr,"<string>%s</string>\n",message);
  fprintf(outfile_ptr,"</void>\n");
  fprintf(outfile_ptr,"</object>\n");
  fprintf(outfile_ptr,"</java>\n");
  */
  exit(1);
}
/***********************************************************************

create_pfam_record  -  Create and initialise a new pfam record

***********************************************************************/

struct pfam *create_pfam_record(struct pfam **fst_pfam_ptr,
                                struct pfam **lst_pfam_ptr,
                                char *pfam_id,BOOL copy,
                                struct parameters *params)
{
  struct pfam *pfam_ptr, *first_pfam_ptr, *last_pfam_ptr;

  /* Initialise pfam pointers */
  pfam_ptr = NULL;
  first_pfam_ptr = *fst_pfam_ptr;
  last_pfam_ptr = *lst_pfam_ptr;

  /* Allocate memory for structure to hold pfam info */
  pfam_ptr = (struct pfam *) malloc(sizeof(struct pfam));
  if (pfam_ptr == NULL)
    {
      if (outfile_ptr != NULL)
        write_error_message(MEMORY_ALLOCATION,
                            "Cannot allocate memory for struct pfam",
                            params->output_type);
      else
        printf("*** Cannot allocate memory for struct pfam\n");
      exit(1);
    }

  /* If this is the very first, save its pointer */
  if (first_pfam_ptr == NULL)
    first_pfam_ptr = pfam_ptr;

  /* Add link from previous pfam to the current one */
  if (last_pfam_ptr != NULL)
    last_pfam_ptr->next_pfam_ptr = pfam_ptr;
  last_pfam_ptr = pfam_ptr;

  /* Store known elements of the structure */
  store_string(&(pfam_ptr->pfam_id),pfam_id);
  pfam_ptr->copy = copy;

  /* Initialise remaining fields */
  //  pfam_ptr->offset[DOMSEQS] = -1;
  //  pfam_ptr->offset[DOMAINS] = -1;
  pfam_ptr->extra = FALSE;
  pfam_ptr->offset = -1;
  pfam_ptr->pfam_name = &null;
  pfam_ptr->short_name = &null;
  pfam_ptr->nseqs = 0;
  pfam_ptr->colour = 0;
  pfam_ptr->pstn = 0;
  pfam_ptr->first_offset_ptr = NULL;

  /* Initialise pointer to next pfam record */
  pfam_ptr->next_pfam_ptr = NULL;

  /* Return current pointers */
  *fst_pfam_ptr = first_pfam_ptr;
  *lst_pfam_ptr = last_pfam_ptr;

  /* Return new pfam pointer */
  return(pfam_ptr);
}
/***********************************************************************

get_command_arguments  -  Interpret the command line arguments

***********************************************************************/

void get_command_arguments(char *string[],int ntoken,
                           struct parameters *params)
{
  BOOL copy;

  char architecture[LINELEN + 1], dom_sep, species[SPECIES_LEN + 3];

  char *string_ptr, *token[MAX_TOKENS];

  int i, itoken, n;
  
  struct pfam *pfam_ptr, *last_pfam_ptr;
  
  /* Initialise variables */
  architecture[0] = '\0';
  pfam_ptr = last_pfam_ptr = NULL;

  /* Initialise parameters */
  params->run_type = SEARCH;
  params->cath_domains = FALSE;
  params->old_cath = FALSE;
  params->npfam = 0;
  params->max_architectures = MAX_ARCHITECTURES;
  params->max_sequences = MAX_SEQUENCES;
  params->architecture = &null;
  params->uniprot_acc = &null;
  params->uniprot_id = &null;
  params->species = &null;
  params->filter_seqs_file = &null;
  params->operator = OR;
  params->pdb_only = FALSE;
  params->reviewed_only = FALSE;
  params->output_type = ARCHSCHEMA_OUTPUT;
  params->out_name[0] = '\0';
  params->get_ssg = FALSE;
  params->first_pfam_ptr = NULL;
  params->first_filterseq_ptr = NULL;
  params->combine = FALSE;
  params->run_64 = FALSE;

  /* If -help entered on the command line, show options available */
  if (ntoken < 1 || (ntoken == 1 && (!strcmp(string[1],"-help") ||
                                     !strcmp(string[1],"-h"))))
    {
      if (outfile_ptr != NULL)
        write_error_message(NO_PARAMETERS,
                            "No parameters entered for program archindex",
                            params->output_type);
      else
        {
          printf("\n");
          printf("Correct usage is ...\n");
          printf("\n");
          printf("    archindex [-c]  [-cath]  [-p pfam_id|cath_id]  "
                 "[-u uniprot_code]\n"
                 "              [-s species]  [-o operator]  [-pdb]  "
                 "[-rev]\n"
                 "              [-maxa max_no]  [-maxs max_no]  [-ssg]\n"
                 "              [-output file_name]  [-dat]  [-old]\n"
                 "              [-inclist file_name]  [-combine]  [-64]\n");
          printf("where\n");
          printf("     * -c          = compile ArchSchema indexes using the "
                 "relevant local\n"
                 "                     files\n");
          printf("     * -cath       = use CATH domain data rather than "
                 "Pfam\n");
          printf("     * -maxa max_no  = maximum number of architectures "
                 "to return\n");
          printf("     * -maxs max_no  = maximum number of sequences "
                 "to return (-1 to return all)\n");
          printf("     * -p pfam_id|cath_id = Pfam-id/architecture, or CATH "
                 "id/architecture \n"
                 "                     to search for. CATH ids should contain "
                 "underscores\n");
          printf("     * -u uniprot_code = UniProt accession or id\n");
          printf("     * -s species  = species identifier (eg BACSU)\n");
          printf("     * -o operator = Pfam domain operator (AND/OR)\n");
          printf("     * -pdb        = flag indicating that only architecture "
                 "with structural\n"
                 "                     info are required\n");
          printf("     * -rev        = flag indicating that only reviewed "
                 "UniProt sequences\n"
                 "                     to be included\n");
          printf("     * -ssg        = flag indicating that SSG clusters "
                 "to be applied to data\n");
          printf("     * -output file_name = write output to file\n");
          printf("     * -dat        = flag indicating that output to be "
                 "in .dat format\n");
          printf("     * -old        = flag indicating that old version of "
                 "CATH data to be used\n");
          printf("     * -inclist file_name = file listing sequence "
                 "identifiers to be included\n");
          printf("     * -combine    = create a dummy architecture from "
                 "the domain list entered to connect\n"
                 "unconnected domain maps\n");
          printf("     * -64         = run on 64-bit processor\n");
          printf("\n");
        }
      exit(1);
    }

  /* Loop through any command arguments */
  for (itoken = 1; itoken < ntoken + 1; itoken++)
    {
      /* Check for the -c option */
      if (!strcmp(string[itoken],"-c"))
        params->run_type = COMPILE_INDEXES;

      /* Check for the -cath option indicating that CATH, rather than Pfam,
         domains to be used */
      else if (!strcmp(string[itoken],"-cath"))
        params->cath_domains = TRUE;

      /* Check for the -p option giving the Pfam id or architecture to
         be searched */
      else if (!strcmp(string[itoken],"-p"))
        {
          /* Get the Pfam id from the next token, if there is one */
          if (itoken < ntoken)
            {
              strcpy(architecture,string[itoken + 1]);
              itoken++;
            }          
        }

      /* Check for the -u option giving the UniProt accession or identifier
         code */
      else if (!strcmp(string[itoken],"-u"))
        {
          /* Get the code from the next token, if there is one */
          if (itoken < ntoken)
            {
              if (strchr(string[itoken + 1],'_') != NULL)
                store_string(&(params->uniprot_id),string[itoken + 1]);
              else
                store_string(&(params->uniprot_acc),string[itoken + 1]);
              itoken++;
            }          
        }

      /* Check for the -s option giving the selected species identifier */
      else if (!strcmp(string[itoken],"-s"))
        {
          /* Get the identifier from the next token, if there is one */
          if (itoken < ntoken && strcmp(string[itoken + 1],"ALL") &&
              strlen(string[itoken + 1]) < SPECIES_LEN)
            {
              /* If species string is empty, then store this species */
              if (params->species == &null)
                {
                  /* Form species name, bounded by separators */
                  strcpy(species,"|");
                  strcat(species,string[itoken + 1]);
                  strcat(species,"|");

                  /* Store it */
                  store_string(&(params->species),species);
                }

              /* Otherwise, append to current species list */
              else
                {
                  /* Form species name, terminated by separator */
                  strcpy(species,string[itoken + 1]);
                  strcat(species,"|");

                  /* Store it */
                  append_string(&(params->species),species);
                }
              itoken++;
            }          
        }

      /* Check for the -maxa option giving the maximum number of
         architectures to output */
      else if (!strcmp(string[itoken],"-maxa"))
        {
          /* Get the number from the next token, if there is one */
          if (itoken < ntoken)
            {
              params->max_architectures = atoi(string[itoken + 1]);
              itoken++;
            }          
        }

      /* Check for the -maxs option giving the maximum number of
         sequences to output */
      else if (!strcmp(string[itoken],"-maxs"))
        {
          /* Get the number from the next token, if there is one */
          if (itoken < ntoken)
            {
              params->max_sequences = atoi(string[itoken + 1]);
              itoken++;
            }          
        }

      /* Check for the -o option giving the Pfam domain operator */
      else if (!strcmp(string[itoken],"-o"))
        {
          /* Get the operator from the next token, if there is one */
          if (itoken < ntoken)
            {
              if (!strcmp(string[itoken + 1],"AND"))
                params->operator = AND;
              else if (!strcmp(string[itoken + 1],"OR"))
                params->operator = OR;
              itoken++;
            }          
        }

      /* Check for the -pdb option indicating only hits with structural
         info to be returned */
      else if (!strcmp(string[itoken],"-pdb"))
        params->pdb_only = TRUE;

      /* Check for the -rev option indicating only reviewed UniProt
         sequences to be included */
      else if (!strcmp(string[itoken],"-rev"))
        params->reviewed_only = TRUE;

      /* Check for the -ssg option indicating that sequences to be
         annotated with SSG data */
      else if (!strcmp(string[itoken],"-ssg"))
        params->get_ssg = TRUE;

      /* Check for the -dat option indicating text output listing
         PDB structures for a given Pfam domain */
      else if (!strcmp(string[itoken],"-dat"))
        params->output_type = DAT_OUTPUT;

      /* Check for the -output option giving the name of the output file */
      else if (!strcmp(string[itoken],"-output"))
        {
          /* Get the file name from the next token, if there is one */
          if (itoken < ntoken)
            {
              strcpy(params->out_name,string[itoken + 1]);
              itoken++;
            }          
        }

      /* Check for the -inclist option giving the name of the file
       containing the list of sequences to be included */
      else if (!strcmp(string[itoken],"-inclist"))
        {
          /* Get the file name from the next token, if there is one */
          if (itoken < ntoken)
            {
              store_string(&(params->filter_seqs_file),string[itoken + 1]);
              itoken++;
            }          
        }

      /* Check for the -old option indicating that old CATH to be used */
      else if (!strcmp(string[itoken],"-old"))
        params->old_cath = TRUE;

      /* Check for the -combine option indicating that the entered 
         architecture may be a fake one intended to pull in unconnected
         domain architectures */
      else if (!strcmp(string[itoken],"-combine"))
        params->combine = TRUE;

      /* Check for the 64-bit option */
      else if (!strcmp(string[itoken],"-64"))
        params->run_64 = TRUE;
    }

  /* If we have an architecture, split it up into its tokens */
  if (architecture[0] != '\0')
    {
      /* Save the architecture */
      store_string(&(params->architecture),architecture);

      /* Determine the domain separator */
      if (params->cath_domains == TRUE)
        dom_sep = '_';
      else
        dom_sep = '.';

      /* Assume have architecture, being a list of dot-separated
         Pfam domains */
      n = tokenize(architecture,token,MAX_TOKENS,dom_sep);

      /* Loop over the Pfam ids, creating a Pfam record for
         each */
      for (i = 0; i < n; i++)
        {
          /* Check that we don't already have this Pfam domain */
          pfam_ptr = find_pfam(params->first_pfam_ptr,token[i]);

          /* Set copy flag depending on whether this is a
             copy of an earlier domain */
          if (pfam_ptr != NULL)
            copy = TRUE;
          else
            copy = FALSE;

          /* Create a new Pfam record */
          pfam_ptr
            = create_pfam_record(&(params->first_pfam_ptr),
                                 &last_pfam_ptr,token[i],copy,params);

          /* Increment count of Pfam domains to search on */
          params->npfam++;
        }
    }
  else
    params->combine = FALSE;

  /* If combine option selected, make sure there is no UniProt id */
  if (params->combine == TRUE)
    params->uniprot_acc = params->uniprot_id = &null;

  /* Check that have a sensible number of architectures */
  if (params->max_architectures < 1)
    params->max_architectures = MAX_ARCHITECTURES;
}
/***********************************************************************

create_arch_record  -  Create and initialise a new pfam record

***********************************************************************/

struct arch *create_arch_record(struct arch **fst_arch_ptr,
                                struct arch **lst_arch_ptr,
                                char *pfam_id,long offset,
                                struct parameters *params)
{
  struct arch *arch_ptr, *first_arch_ptr, *last_arch_ptr;

  /* Initialise pfam pointers */
  arch_ptr = NULL;
  first_arch_ptr = *fst_arch_ptr;
  last_arch_ptr = *lst_arch_ptr;

  /* Allocate memory for structure to hold pfam info */
  arch_ptr = (struct arch *) malloc(sizeof(struct arch));
  if (arch_ptr == NULL)
    {
      if (params->run_type == COMPILE_INDEXES)
        {
          printf("*** Can't allocate memory for struct arch\n");
          exit (1);
        }
      else
        write_error_message(MEMORY_ALLOCATION,
                            "Cannot allocate memory for struct arch",
                            params->output_type);
    }

  /* If this is the very first, save its pointer */
  if (first_arch_ptr == NULL)
    first_arch_ptr = arch_ptr;

  /* Add link from previous pfam to the current one */
  if (last_arch_ptr != NULL)
    last_arch_ptr->next_arch_ptr = arch_ptr;
  last_arch_ptr = arch_ptr;

  /* Store known elements of the structure */
  store_string(&(arch_ptr->key),pfam_id);
  arch_ptr->offset = offset;

  /* Initialise pointer to next pfam record */
  arch_ptr->next_arch_ptr = NULL;

  /* Return current pointers */
  *fst_arch_ptr = first_arch_ptr;
  *lst_arch_ptr = last_arch_ptr;

  /* Return new pfam pointer */
  return(arch_ptr);
}
/***********************************************************************

get_arch_data  -  Read through the gzipped swisspfam file to pick up
                  offsets of each Pfam archtecture record

***********************************************************************/

void get_arch_data(int type,struct arch **first_arch_ptr,int *nrecords,
                   struct parameters *params)
{
  char file_name[FILENAME_LEN], number_string[20], message[LINELEN + 1];
  char input_line[LINELEN + 1];
  char pfam_id[PFAM_LEN];
  char uniprot_acc[UNIPROT_LEN], uniprot_id[UNIPROT_LEN];
  char architecture[LINELEN + 1], empty[2];

  char *string_ptr, *token[MAX_TOKENS];

  int i, itoken, jtoken, len, longest, nlines, narch, ntokens;
  int wanted;

  long offset;

  struct arch *arch_ptr, *last_arch_ptr[NKEYS];

  FILE *file_ptr;

  /* Initialise variables */
  narch = nlines = longest = 0;

  /* Initialise pointers */
  for (i = 0; i < NKEYS; i++)
    {
      first_arch_ptr[i] = last_arch_ptr[i] = NULL;
      nrecords[i] = 0;
    }
  arch_ptr = NULL;

  /* Form name of input file */
  strcpy(file_name,"archstruc.uniprot.joined.by.accession");
  if (params->cath_domains == TRUE)
    strcpy(file_name,"cathstruc.uniprot.joined.by.accession");

  /* Open the input data file */
  if ((file_ptr = fopen(file_name,"r")) == NULL)
    {
      /* Otherwise, error must be due to missing file */
      sprintf(message,"Unable to open data file [%s]\n",file_name);
      write_error_message(FILE_NOT_FOUND,message,params->output_type);
    }

  /* Get current file position */
  offset = ftell(file_ptr);

  /* Read through the file */
  while (fgets(input_line,LINELEN,file_ptr) != NULL)
    {
      /* Retrieve the first three tokens: UniProt accession, Pfam
         domain architecture, and UniProt id */
      ntokens = tokenize(input_line,token,MAX_TOKENS,'\t');

      /* If have three, save them */
      if (ntokens > 3)
        {
          /* Save the data read in */
          strcpy(uniprot_acc,token[0]);
          strcpy(architecture,token[1]);
          strcpy(uniprot_id,token[2]);

          /* Create a record to store the location of this UniProt
             accession number */
          arch_ptr
            = create_arch_record(&(first_arch_ptr[KEY_UNIPROT_ACC]),
                                 &(last_arch_ptr[KEY_UNIPROT_ACC]),
                                 uniprot_acc,offset,params);

          /* Increment count of UniProt records created */
          nrecords[KEY_UNIPROT_ACC]++;

          /* Repeat for the UniProt id */
          arch_ptr
            = create_arch_record(&(first_arch_ptr[KEY_UNIPROT_ID]),
                                 &(last_arch_ptr[KEY_UNIPROT_ID]),
                                 uniprot_id,offset,params);

          /* Increment count of UniProt records created */
          nrecords[KEY_UNIPROT_ID]++;

          /* Split the architecture up into its constituent Pfam
             domains */
          ntokens
            = tokenize(architecture,token,MAX_TOKENS,DOM_SEP[type][0]);

          /* Create a Pfam record for each unique domain in this
             architecture */
          for (itoken = 0; itoken < ntokens; itoken++)
            {
              /* Initialise flag */
              wanted = TRUE;

              /* Check previous domains to ensure this one is a new one */
              for (jtoken = 0; jtoken < itoken && wanted == TRUE; jtoken++)
                {
                  /* If domain matches, then don't want it */
                  if (!strcmp(token[jtoken],token[itoken]))
                    wanted = FALSE;
                }

              /* If wanted, then create a Pfam record */
              if (wanted == TRUE)
                {
                  arch_ptr
                    = create_arch_record(&(first_arch_ptr[KEY_PFAM_ID]),
                                         &(last_arch_ptr[KEY_PFAM_ID]),
                                         token[itoken],offset,params);

                  /* Increment count of Pfam records created */
                  nrecords[KEY_PFAM_ID]++;
                }
            }
        }

      /* Increment line count */
      nlines++;

      /* Get current file position */
      offset = ftell(file_ptr);
    }

  /* Close the input data file */
  fclose(file_ptr);
}
/***********************************************************************

arch_sort  -  Shell-sort to arrange the Pfam records in ascending order

***********************************************************************/

void arch_sort(struct arch **arch_index_ptr,int narch)
{
  int match;
  long endloop, i, indi, indl, j, k, l, lognb2, m, n, nn;
  float aln2i = (float) 1.4426950, tiny = (float) 1.e-5;
  float calc;
  double dble, dble1;

  struct arch *arch1_ptr, *arch2_ptr, *swap_arch_ptr;

  /*--Initialise values */
  n = narch;
  dble = n;
  dble1 = log(dble);
  calc = (float) dble1;
  calc = calc * aln2i + tiny;
  lognb2 = (long) calc;

  /* Perform sort */
  m = n;
  for (nn = 1; nn <= lognb2; nn++)
    {
      m = m / 2;
      k = n - m;
      for (j = 1; j <= k; j++)
        {
          i = j;
          endloop = FALSE;
          while (endloop == FALSE)
            {
              l = i + m;
              indi = i - 1;
              indl = l - 1;

              /* Get the two code records */
              arch1_ptr = arch_index_ptr[indi];
              arch2_ptr = arch_index_ptr[indl];

              /* Determine which string is greater */
              match = strcmp(arch1_ptr->key,arch2_ptr->key);

              /* If in wrong order, then swap sort-pointers */
              if (match > 0 || (match == 0 &&
                                arch1_ptr->offset > arch2_ptr->offset))
                {
                  swap_arch_ptr = arch_index_ptr[indi];
                  arch_index_ptr[indi] = arch_index_ptr[indl];
                  arch_index_ptr[indl] = swap_arch_ptr;

                  i = i - m;
                  if (i < 1)
                    endloop = TRUE;
                }
              else
                endloop = TRUE;
            }
        }
    }
}
/***********************************************************************

sort_arch_records  -  Create index of architecture records

***********************************************************************/

void sort_arch_records(struct arch *first_arch_ptr,
                       struct arch ***arch_ind_ptr,int narch,
                       struct parameters *params)
{
  int nentry;

  struct arch *arch_ptr;

  struct arch **arch_index_ptr;

  /* Create the pfam array as one large array */
  arch_index_ptr = (struct arch **) malloc(narch * sizeof(struct arch *));

  /* Check that sufficient memory was allocated */
  if (arch_index_ptr == NULL)
    write_error_message(MEMORY_ALLOCATION,
                        "Cannot allocate memory for pfam array",
                        params->output_type);

  /* Get pointer to the first code entry */
  arch_ptr = first_arch_ptr;
  nentry = 0;

  /* Loop through all entries to place each one in the sorted pfam
     array */
  while (arch_ptr != NULL)
    {
      /* Store current pointer at end of the pfam array */
      if (nentry < narch)
        arch_index_ptr[nentry] = arch_ptr;
      else
        {
          printf("*** No. architectures > size of index!\n");
          exit(1);
        }

      /* Increment count of stored entries */
      nentry++;

      /* Get pointer to the next entry in the linked list */
      arch_ptr = arch_ptr->next_arch_ptr;
    }

  /* If number of records doesn't match size of index, then have error */
  if (nentry != narch)
    {
      printf("*** No. architectures (%d) != size of index (%d)!\n",
             nentry,narch);
      exit(1);
    }

  /* Sort the code pointers into order of Pfam id */
  if (narch > 1)
    arch_sort(arch_index_ptr,narch);

  /* Return the array pointer */
  *arch_ind_ptr = arch_index_ptr;
}
/***********************************************************************

write_offsets_file  -  Write the list of UniProt codes and their offsets

***********************************************************************/

void write_offsets_file(int type,int ikey,struct arch **arch_index_ptr,
                        int narch)
{
  int iarch;

  struct arch *arch_ptr;

  FILE *file_ptr;

  /* Open the output arch file */
  printf("Opening output file %s ...\n",DAT_FILE_NAME[type][ikey]);
  fflush(stdout);
  file_ptr = open_output_file(DAT_FILE_NAME[type][ikey],TRUE);

  /* Loop over all the records to be written out */
  for (iarch = 0; iarch < narch; iarch++)
    {
      /* Get pointer to this record */
      arch_ptr = arch_index_ptr[iarch];

      /* Write out the key and offset */
      fprintf(file_ptr,"%s %ld\n",arch_ptr->key,arch_ptr->offset);

      /* Free up the memory taken by this record */
      //OUT if (arch_ptr->key != &null)
        //OUT free(arch_ptr->key);
      free(arch_ptr);
    }

  /* Close the output arch file */
  fclose(file_ptr);

  /* Free up the memory taken by the index */
  free(arch_index_ptr);

  /* Show number of offset records written out */
  printf("Number of offset records written out = %8d\n",narch);
  fflush(stdout);

}
/***********************************************************************

write_offsets_index  -  Create an index to the offsets file

***********************************************************************/

void write_offsets_index(int type,int ikey,struct parameters *params)
{
  char file_name[FILENAME_LEN], input_line[LINELEN + 1];
  char code[UNIPROT_LEN], message[LINELEN + 1];

  char *string_ptr;

  int len, nlines, next_index;

  long nout, offset;

  FILE *file_ptr, *idx_file_ptr;

  /* Initialise variables */
  nlines = nout = 0;
  next_index = INDEX_GAP;

  /* Get the name of the input file */
  strcpy(file_name,DAT_FILE_NAME[type][ikey]);

  /* Open the input file */
  printf("Opening input file %s ...\n",file_name);
  fflush(stdout);
  if ((file_ptr = fopen(file_name,"r")) == NULL)
    {
      /* Otherwise, error must be due to missing file */
      sprintf(message,"Unable to open data file [%s]\n",file_name);
      write_error_message(FILE_NOT_FOUND,message,params->output_type);
    }

  /* Open the output index file */
  printf("Opening output file %s ...\n",IDX_FILE_NAME[type][ikey]);
  fflush(stdout);
  idx_file_ptr = open_output_file(IDX_FILE_NAME[type][ikey],TRUE);

  /* Get current file position */
  offset = ftell(file_ptr);

  /* Loop while reading the file */
  while (fgets(input_line,LINELEN,file_ptr) != NULL)
    {
      /* Increment line count */
      nlines++;

      /* If this is the next indexing position, write out the code
         and offset */
      if (nlines > next_index)
        {
          /* Get the end of the code */
          string_ptr = strchr(input_line,' ');

          /* If found, then extract code */
          if (string_ptr != NULL)
            {
              /* Get the code if not too long */
              string_ptr[0] = '\0';
              len = strlen(input_line);
              if (len < UNIPROT_LEN)
                {
                  /* Write out the offset */
                  strcpy(code,input_line);
                  fprintf(idx_file_ptr,"%s %ld\n",code,offset);

                  /* Increment count */
                  nout++;
                }
            }

          /* Get the next indexing position */
          next_index = nlines + INDEX_GAP;
        }

      /* Get current file position */
      offset = ftell(file_ptr);
    }

  /* Show number of index records written out */
  printf("Number of offset indexes written out = %16ld\n",nout);
  fflush(stdout);

  /* Close both files */
  fclose(file_ptr);
  fclose(idx_file_ptr);
}
/***********************************************************************

create_offset_record  -  Create and initialise a new offset record

***********************************************************************/

struct offset *create_offset_record(struct offset **fst_offset_ptr,
                                    struct offset **lst_offset_ptr,
                                    long offset,
                                    struct parameters *params)
{
  struct offset *offset_ptr, *first_offset_ptr, *last_offset_ptr;

  /* Initialise offset pointers */
  offset_ptr = NULL;
  first_offset_ptr = *fst_offset_ptr;
  last_offset_ptr = *lst_offset_ptr;

  /* Allocate memory for structure to hold offset info */
  offset_ptr = (struct offset *) malloc(sizeof(struct offset));
  if (offset_ptr == NULL)
    {
      if (params->run_type == COMPILE_INDEXES)
        {
          printf("*** Can't allocate memory for struct offset\n");
          exit (1);
        }
      else
        write_error_message(MEMORY_ALLOCATION,
                            "Cannot allocate memory for struct offset",
                            params->output_type);
    }

  /* If this is the very first, save its pointer */
  if (first_offset_ptr == NULL)
    first_offset_ptr = offset_ptr;

  /* Add link from previous offset to the current one */
  if (last_offset_ptr != NULL)
    last_offset_ptr->next_offset_ptr = offset_ptr;
  last_offset_ptr = offset_ptr;

  /* Store known elements of the structure */
  offset_ptr->offset = offset;

  /* Initialise pointer to next offset record */
  offset_ptr->next_offset_ptr = NULL;

  /* Return current pointers */
  *fst_offset_ptr = first_offset_ptr;
  *lst_offset_ptr = last_offset_ptr;

  /* Return new offset pointer */
  return(offset_ptr);
}
/***********************************************************************

transfer_to_pfam  -  Transfer the Pfam ids and offsets to a new set of
                     Pfam records

***********************************************************************/

int transfer_to_pfam(struct pfam **fst_pfam_ptr,
                     struct arch **arch_index_ptr,int narch,
                     struct parameters *params)
{
  BOOL copy;

  char last_pfam_id[PFAM_LEN];

  int iarch, npfam;

  struct arch *arch_ptr;
  struct offset *offset_ptr, *last_offset_ptr;
  struct pfam *pfam_ptr, *first_pfam_ptr, *last_pfam_ptr;

  /* Initialise variables */
  printf("In transfer_to_pfam ... narch = %d\n",narch);
  fflush(stdout);
  copy = FALSE;
  last_pfam_id[0] = '\0';
  npfam = 0;

  /* Initialise pointers */
  pfam_ptr = first_pfam_ptr = last_pfam_ptr = NULL;

  /* Loop over the sorted list of Pfam ids and their offset */
  for (iarch = 0; iarch < narch; iarch++)
    {
      /* Get this record */
      arch_ptr = arch_index_ptr[iarch];

      /* Process if point not null */
      if (arch_ptr != NULL)
        {
          /* If it is a new Pfam id, then create a new Pfam record to store
             it */
          if (strcmp(arch_ptr->key,last_pfam_id))
            {
              /* Create a new Pfam record */
              pfam_ptr
                = create_pfam_record(&first_pfam_ptr,&last_pfam_ptr,
                                     arch_ptr->key,copy,params);

              /* Increment count of Pfam records */
              npfam++;

              /* Save the Pfam id */
              strcpy(last_pfam_id,arch_ptr->key);

              /* Initialize the last offset pointer */
              last_offset_ptr = NULL;
            }

          /* Create an offset record to store the current offset in the
             archstruc.uniprot.joined.by.accession file */
          create_offset_record(&(pfam_ptr->first_offset_ptr),&last_offset_ptr,
                               arch_ptr->offset,params);

          /* Free the memory taken by the current arch record */
          free(arch_ptr);
        }
    }

  /* Free up the memory taken by the index */
  free(arch_index_ptr);

  /* Return first Pfam pointer */
  *fst_pfam_ptr = first_pfam_ptr;

  /* Show number of records created */
  printf("Number of Pfam records created = %8d\n",npfam);

  /* Return number of Pfam records created */
  return(npfam);
}
/***********************************************************************

get_offsets  -  Read in the offsets for the given file

***********************************************************************/

void get_offsets(struct pfam *first_pfam_ptr,char *file_name,
                 int type,struct parameters *params)
{
  char input_line[LINELEN + 1], last_pfam_id[PFAM_LEN + 1];
  char message[LINELEN + 1];

  char *string_ptr;

  int len, nlines, nmatched;
  int done, next_code, next_record;

  long offset;

  struct pfam *pfam_ptr;

  FILE *file_ptr;

  /* Initialise variables */
  done = FALSE;
  last_pfam_id[0] = '\0';
  next_code = TRUE;
  next_record = TRUE;
  nlines = nmatched = 0;

  /* Get pointer to the first pfam record */
  pfam_ptr = first_pfam_ptr;
  next_code = FALSE;
  if (pfam_ptr == NULL)
    done = TRUE;

  /* Open the input data file */
  if ((file_ptr = fopen(file_name,"r")) == NULL)
    {
      /* Otherwise, error must be due to missing file */
      sprintf(message,"Unable to open data file [%s]\n",file_name);
      write_error_message(FILE_NOT_FOUND,message,params->output_type);
    }

  /* Get current file position */
  offset = ftell(file_ptr);

  /* Loop until all codes have been assigned a protein name */
  while (done == FALSE)
    {
      /* If need to go to the next code, then do so */
      if (next_code == TRUE)
        {
          /* Get pointer to the next pfam record */
          if (pfam_ptr != NULL)
            pfam_ptr = pfam_ptr->next_pfam_ptr;
          else
            done = TRUE;

          /* Check for end */
          if (pfam_ptr == NULL)
            done = TRUE;

          /* Set flag */
          next_code = FALSE;
        }

      /* If need to read to next record, then loop until we have the
         next useful one */
      while (done == FALSE && next_record == TRUE)
        {
          /* Read in the next record from the file */
          if (fgets(input_line,LINELEN,file_ptr) != NULL)
            {
              /* Increment line count */
              nlines++;

              /* Get the Pfam id */
              string_ptr = strchr(input_line,'\t');

              /* If have a tab, terminate the string */
              if (string_ptr != NULL)
                {
                  string_ptr[0] = '\0';
                  len = strlen(input_line);

                  /* Check that Pfam id not too long */
                  if (len > PFAM_LEN - 1)
                    {
                      sprintf(message,"Invalid Pfam id: [%s]\n",input_line);
                      write_error_message(FILE_READ_ERROR,message,
                                          params->output_type);
                    }

                  /* Check that Pfam ids in the file are in sorted
                     order */
                  if (strcmp(input_line,last_pfam_id) < 0)
                    {
                      sprintf(message,"Codes not sorted current [%s] "
                              "lt previous [%s]\n",input_line,
                              last_pfam_id);
                      write_error_message(FILE_READ_ERROR,message,
                                          params->output_type);
                    }

                  /* Save the Pfam id */
                  strcpy(last_pfam_id,input_line);

                  /* If this is less than our current Pfam record,
                     need to keep reading */
                  if (strcmp(input_line,pfam_ptr->pfam_id) < 0)
                    next_record = TRUE;
                  else
                    next_record = FALSE;
                }



              /* Get current file position */
              offset = ftell(file_ptr);
            }

          /* Have reached the end of the file, so we're done */
          else
            done = TRUE;
        }

      /* If not yet done, check the match */
      if (done == FALSE)
        {
          /* If this matches our current Pfam record, store the offset */
          if (!strcmp(input_line,pfam_ptr->pfam_id))
            {
              /* Store the offset */
              //              pfam_ptr->offset[type] = offset;

              /* Increment count of matches */
              nmatched++;

              /* Set both flags to go to the next record */
              next_code = TRUE;
              next_record = TRUE;
            }

          /* Otherwise, if it is higher, need to get the next code */
          else if (strcmp(input_line,pfam_ptr->pfam_id) > 0)
            {
              /* Set flags */
              next_code = TRUE;
              next_record = FALSE;
            }

          /* Otherwise, must be lower, so need to read the next
             record */
          else
              next_record = TRUE;
        }
    }

  /* Close the input data file */
  fclose(file_ptr);
}
/***********************************************************************

write_pfam_index  -  Write out the Pfam offsets index

***********************************************************************/

void write_pfam_index(char *file_name,struct pfam *first_pfam_ptr)
{
  struct offset *offset_ptr;
  struct pfam *pfam_ptr;

  FILE *file_ptr;

  /* Open the output pfam file */
  printf("Opening output file %s ...\n",file_name);
  fflush(stdout);
  file_ptr = open_output_file(file_name,TRUE);

  /* Get pointer to the first Pfam record */
  pfam_ptr = first_pfam_ptr;

  /* Loop through all the Pfam records to write out */
  while (pfam_ptr != NULL)
    {
      /* Write out the Pfam identifier */
      fprintf(file_ptr,":%s\n",pfam_ptr->pfam_id);

      /* Write out the offsets to the pfam_domains.parsed.sorted and
         domseqs.struc files */
      //      fprintf(file_ptr,"%d\n",pfam_ptr->offset[DOMSEQS]);
      //      fprintf(file_ptr,"%d\n",pfam_ptr->offset[DOMAINS]);

      /* Get this Pfam entry's first offset record */
      offset_ptr = pfam_ptr->first_offset_ptr;

      /* Loop through all the offset records to write out */
      while (offset_ptr != NULL)
        {
          /* Write out the current offset to the
             archstruc.uniprot.joined.by.accession file */
          fprintf(file_ptr,"%ld\n",offset_ptr->offset);

          /* Get the next Offset record */
          offset_ptr = offset_ptr->next_offset_ptr;
        }

      /* get the next Pfam record */
      pfam_ptr = pfam_ptr->next_pfam_ptr;
    }

  /* Close the output pfam file */
  fclose(file_ptr);
}
/***********************************************************************

index_to_index  -  Create an index file, archpfam.idx, to archpfam.dat

***********************************************************************/

void index_to_index(char *file_name,char *idx_file_name)
{
  char input_line[LINELEN + 1];

  int len, nlines, npfam;

  long offset;

  FILE *file_ptr, *idx_file_ptr;

  /* Initialise variables */
  nlines = npfam = 0;

  /* Open the input file */
  if ((file_ptr = fopen(file_name,"r")) == NULL)
    {
      /* Otherwise, error must be due to missing file */
      printf("\n*** Unable to open data file [%s]\n",file_name);
      exit(1);
    }

  /* Open the output index file */
  idx_file_ptr = open_output_file(idx_file_name,TRUE);

  /* Show message */
  printf("Reading file %s ...\n",file_name);
  fflush(stdout);

  /* Get current file position */
  offset = ftell(file_ptr);

  /* Loop while reading the file */
  while (fgets(input_line,LINELEN,file_ptr) != NULL)
    {
      /* Increment line count */
      nlines++;

      /* If this is a Pfam id, then write out its offset */
      if (input_line[0] == ':')
        {
          /* Truncate the string */
          len = string_truncate(input_line,LINELEN);

          /* Write out */
          fprintf(idx_file_ptr,"%s %ld\n",input_line+1,offset);

          /* Increment count */
          npfam++;
        }

      /* Get current file position */
      offset = ftell(file_ptr);
    }

  /* Show number of Pfam domains indexed */
  printf("  Number of Pfam domains indexed: %8d\n",npfam);
  fflush(stdout);

  /* Close both files */
  fclose(file_ptr);
  fclose(idx_file_ptr);
}
/***********************************************************************

read_uniprot_idx  -  Read archuacc.idx or archuid.idx index file to
                     pick up the starting offset for the corresponding
                     .dat file

***********************************************************************/

long read_uniprot_idx(int type,int ikey,char *code,char *archschema_dir,
                      struct parameters *params)
{
  char file_name[FILENAME_LEN], input_line[LINELEN + 1];
  char number_string[20], message[LINELEN + 1];

  char *string_ptr;

  int len, nlines;
  int done;

  long offset;

  struct pfam *pfam_ptr;

  FILE *file_ptr;

  /* Initialise variables */
  done = FALSE;
  offset = 0;

  /* Form name of input file */
  strcpy(file_name,archschema_dir);
  strcat(file_name,"/");
  strcat(file_name,IDX_FILE_NAME[type][ikey]);

  /* Open the index file */
  if ((file_ptr = fopen(file_name,"r")) == NULL)
    {
      /* Otherwise, error must be due to missing file */
      sprintf(message,"Unable to open data file [%s]\n",file_name);
      write_error_message(FILE_NOT_FOUND,message,params->output_type);
    }

  /* Read the file until we hit the right start-point */
  while (fgets(input_line,LINELEN,file_ptr) != NULL && done == FALSE)
    {
      /* Find the space character */
      string_ptr = strchr(input_line,' ');

      /* If found, extract the UniProt code */
      if (string_ptr != NULL)
        {
          /* Truncate the string */
          len = string_chop(input_line);

          /* Terminate the code */
          string_ptr[0] = '\0';

          /* If this code is larger than ours, then we are done */
          if (strcmp(input_line,code) > 0)
            done = TRUE;

          /* If we're not done, extract the offset from the current
             record */
          if (done == FALSE)
            {
              /* Extract the offset */
              strcpy(number_string,string_ptr + 1);
              offset = atol(number_string);
            }
        }
    }

  /* Return the starting offset */
  return(offset);
}
/***********************************************************************

read_uniprot_dat  -  Read archuacc.dat or archuid.dat file, from the
                     given offset to find the offset(s) for our UniProt
                     code in the architectures file

***********************************************************************/

long read_uniprot_dat(int type,int ikey,char *code,long offset,
                      struct offset **fst_offset_ptr,char *archschema_dir,
                      struct parameters *params)
{
  char file_name[FILENAME_LEN], input_line[LINELEN + 1];
  char number_string[20], message[LINELEN + 1];

  char *string_ptr;

  int len, nlines;
  int done, error;

  long noffset;

  struct offset *offset_ptr, *first_offset_ptr, *last_offset_ptr;

  FILE *file_ptr;

  /* Initialise variables */
  done = FALSE;
  noffset = 0;

  /* Initialise pointers */
  offset_ptr = first_offset_ptr = last_offset_ptr = NULL;

  /* Form name of input file */
  strcpy(file_name,archschema_dir);
  strcat(file_name,"/");
  strcat(file_name,DAT_FILE_NAME[type][ikey]);

  /* Open the index file */
  if ((file_ptr = fopen(file_name,"r")) == NULL)
    {
      /* Otherwise, error must be due to missing file */
      sprintf(message,"Unable to open data file [%s]\n",file_name);
      write_error_message(FILE_NOT_FOUND,message,params->output_type);
    }

  /* Find start-point in data file using fseek */
  error = fseek(file_ptr,offset,SEEK_SET);

  /* Read the file until we hit the right start-point */
  while (fgets(input_line,LINELEN,file_ptr) != NULL && done == FALSE)
    {
      /* Find the space character */
      string_ptr = strchr(input_line,' ');

      /* If found, extract the UniProt code */
      if (string_ptr != NULL)
        {
          /* Terminate the code */
          string_ptr[0] = '\0';

          /* If this is our code, then extract the offset and save */
          if (!strcmp(input_line,code))
            {
              /* Truncate the string */
              len = string_chop(string_ptr + 1);

              /* Extract the offset */
              offset = atol(string_ptr + 1);

              /* Create an offset record */
              offset_ptr
                = create_offset_record(&first_offset_ptr,
                                       &last_offset_ptr,offset,params);

              /* Increment count of records */
              noffset++;
            }

          /* Otherwise, if this code is larger than ours, then we
             are done */
          else if (strcmp(input_line,code) > 0)
            done = TRUE;
        }
    }

  /* Return the first of the offset records */
  *fst_offset_ptr = first_offset_ptr;

  /* Return the number of offset records we've created */
  return(noffset);
}
/***********************************************************************

create_hit_record  -  Create and initialise a new hit record

***********************************************************************/

struct hit *create_hit_record(struct hit **fst_hit_ptr,
                              struct hit **lst_hit_ptr,
                              char *uniprot_acc,char *uniprot_id,
                              char *architecture,char *gene,
                              char *protein_name,int seq_len,
                              char *pdb_codes,char *coverage,
                              struct parameters *params)
{
  struct hit *hit_ptr, *first_hit_ptr, *last_hit_ptr;

  /* Initialise hit pointers */
  hit_ptr = NULL;
  first_hit_ptr = *fst_hit_ptr;
  last_hit_ptr = *lst_hit_ptr;

  /* Allocate memory for structure to hold hit info */
  hit_ptr = (struct hit *) malloc(sizeof(struct hit));
  if (hit_ptr == NULL)
    write_error_message(MEMORY_ALLOCATION,
                        "Cannot allocate memory for struct hit",
                        params->output_type);

  /* If this is the very first, save its pointer */
  if (first_hit_ptr == NULL)
    first_hit_ptr = hit_ptr;

  /* Add link from previous hit to the current one */
  if (last_hit_ptr != NULL)
    last_hit_ptr->next_hit_ptr = hit_ptr;
  last_hit_ptr = hit_ptr;

  /* Store known elements of the structure */
  store_string(&(hit_ptr->uniprot_acc),uniprot_acc);
  store_string(&(hit_ptr->uniprot_id),uniprot_id);
  store_string(&(hit_ptr->architecture),architecture);
  store_string(&(hit_ptr->gene),gene);
  store_string(&(hit_ptr->protein_name),protein_name);
  store_string(&(hit_ptr->pdb_codes),pdb_codes);
  store_string(&(hit_ptr->coverage),coverage);
  hit_ptr->seq_len = seq_len;

  /* Initialise pointer to next hit record */
  hit_ptr->next_hit_ptr = NULL;

  /* Return current pointers */
  *fst_hit_ptr = first_hit_ptr;
  *lst_hit_ptr = last_hit_ptr;

  /* Return new hit pointer */
  return(hit_ptr);
}
/***********************************************************************

get_seq_architecture  -  Get the architecture for the give UniProt seq
                         from the archstruc.uniprot.joined.by.accession
                         file

***********************************************************************/

int get_seq_architecture(int type,struct offset *first_offset_ptr,
                         struct hit **fst_hit_ptr,char *archschema_dir,
                         struct parameters *params)
{
  char file_name[FILENAME_LEN], input_line[LINELEN + 1];
  char message[LINELEN + 1];

  char *string_ptr, *token[MAX_TOKENS];

  int itoken, len, line, nhits, ntokens;
  int done, error, wanted;

  long offset;

  struct hit *hit_ptr, *first_hit_ptr, *last_hit_ptr;
  struct offset *offset_ptr;

  FILE *file_ptr;

  /* Initialise variables */
  nhits = 0;

  /* Initialise pointers */
  hit_ptr = first_hit_ptr = last_hit_ptr = NULL;

  /* Form name of input file */
  strcpy(file_name,archschema_dir);
  strcat(file_name,"/");
  if (type == PFAM)
    strcat(file_name,"archstruc.uniprot.joined.by.accession");
  else
    strcat(file_name,"cathstruc.uniprot.joined.by.accession");

  /* Open the input file */
  if ((file_ptr = fopen(file_name,"r")) == NULL)
    {
      /* Otherwise, error must be due to missing file */
      sprintf(message,"Unable to open data file [%s]\n",file_name);
      write_error_message(FILE_NOT_FOUND,message,params->output_type);
    }

  /* Get pointer to the first offset record */
  offset_ptr = first_offset_ptr;

  /* Loop over the offset records */
  while (offset_ptr != NULL)
    {
      /* Get the offset from this record */
      offset = offset_ptr->offset;

      /* If valid, then find this point in the file */
      if (offset > -1)
        {
          /* Find start-point in data file using fseek */
          error = fseek(file_ptr,offset,SEEK_SET);

          /* Initialise parameters for this read */
          done = FALSE;
          line = 0;

          /* Read in this record */
          if (fgets(input_line,LINELEN,file_ptr) != NULL)
            {
              /* Truncate the string */
              len = string_chop(input_line);

              /* Get the tokens on this line */
              ntokens = tokenize(input_line,token,MAX_TOKENS,'\t');

              /* Create a hit record to store the data */
              hit_ptr
                = create_hit_record(&first_hit_ptr,&last_hit_ptr,
                                    token[F_UNIPROT_ACC],
                                    token[F_UNIPROT_ID],
                                    token[F_ARCHITECTURE],
                                    token[F_GENE],
                                    token[F_PROTEIN_NAME],
                                    atoi(token[F_SEQ_LENGTH]),
                                    token[F_PDB_CODES],
                                    token[F_COVERAGE],params);

              /* Increment count of hits */
              nhits++;
            }

          /* Increment line counter */
          line++;
        }
      
      /* Get the next offset record */
      offset_ptr = offset_ptr->next_offset_ptr;
    }

  /* Close the file */
  fclose(file_ptr);

  /* Return pointer */
  *fst_hit_ptr = first_hit_ptr;

  /* Return the number of records matched */
  return(nhits);
}
/***********************************************************************

select_hit  -  If we have more than one sequence hit, determine which to
               use

***********************************************************************/

int select_hit(struct hit **fst_hit_ptr)
{
  struct hit *hit_ptr, *first_hit_ptr;

  /* Initialise pointers */
  first_hit_ptr = *fst_hit_ptr;

  /* Get the first hit record */
  hit_ptr = first_hit_ptr;

  /* Loop over the hit records to find the best to use */
  while (hit_ptr != NULL)
    {
/* @@@ */

      /* Get the next hit record */
      hit_ptr = hit_ptr->next_hit_ptr;
    }

  /* Return first pointer */
  *fst_hit_ptr = first_hit_ptr;
}
/***********************************************************************

split_architecture  -  Split up the architecture into its constituent
                       Pfam domains

***********************************************************************/

int split_architecture(int type,struct pfam **fst_pfam_ptr,
                       struct pfam **lst_pfam_ptr,char *architecture,
                       int npfam,struct parameters *params)
{
  BOOL copy;

  char *tmp, *token[MAX_TOKENS];

  int itoken, jtoken, len, ntokens;

  struct pfam *pfam_ptr, *first_pfam_ptr, *last_pfam_ptr;

  /* Initialise variables */
  copy = FALSE;

  /* Get the length of the architecture */
  len = strlen(architecture);

  /* Allocate temporary memory for the architecture */
  tmp = (char *) malloc (sizeof(char)*(len + 1));

  /* Store the architecture */
  strcpy(tmp,architecture);

  /* Initialise pointers */
  pfam_ptr = NULL;
  first_pfam_ptr = *fst_pfam_ptr;
  last_pfam_ptr = *lst_pfam_ptr;

  /* Extract the dot-separated Pfam domains from the architecture */
  ntokens = tokenize(tmp,token,MAX_TOKENS,DOM_SEP[type][0]);

  /* Loop over the Pfam ids, to create a Pfam record for each */
  for (itoken = 0; itoken < ntokens; itoken++)
    {
      /* Determine whether this a copy domain */
      copy = FALSE;
      for (jtoken = 0; jtoken < itoken; jtoken++)
        {
          if (!strcmp(token[itoken],token[jtoken]))
            copy = TRUE;
        }

      /* Create a new Pfam record */
      pfam_ptr
        = create_pfam_record(&first_pfam_ptr,&last_pfam_ptr,
                             token[itoken],copy,params);

      /* If this is the dummy Pfam-B record, save it name */
      if (!strcmp(token[itoken],"PB000000"))
        {
          store_string(&(pfam_ptr->pfam_name),"Pfan-B domain");
          store_string(&(pfam_ptr->short_name),"Pfam-B");
        }

      /* Increment count of records created */
      npfam++;
    }

  /* Return current pointers */
  *fst_pfam_ptr = first_pfam_ptr;
  *lst_pfam_ptr = last_pfam_ptr;

  /* Free up memory taken by the architecture */
  free(tmp);

  /* Return number of Pfam records */
  return(npfam);
}
/***********************************************************************

pfam_sort  -  Shell-sort to arrange the Pfam records in ascending order

***********************************************************************/

void pfam_sort(struct pfam **pfam_index_ptr,int npfam)
{
  int match;
  long endloop, i, indi, indl, j, k, l, lognb2, m, n, nn;
  float aln2i = (float) 1.4426950, tiny = (float) 1.e-5;
  float calc;
  double dble, dble1;

  struct pfam *pfam1_ptr, *pfam2_ptr, *swap_pfam_ptr;

  /*--Initialise values */
  n = npfam;
  dble = n;
  dble1 = log(dble);
  calc = (float) dble1;
  calc = calc * aln2i + tiny;
  lognb2 = (long) calc;

  /* Perform sort */
  m = n;
  for (nn = 1; nn <= lognb2; nn++)
    {
      m = m / 2;
      k = n - m;
      for (j = 1; j <= k; j++)
        {
          i = j;
          endloop = FALSE;
          while (endloop == FALSE)
            {
              l = i + m;
              indi = i - 1;
              indl = l - 1;

              /* Get the two code records */
              pfam1_ptr = pfam_index_ptr[indi];
              pfam2_ptr = pfam_index_ptr[indl];

              /* Determine which string is greater */
              match = strcmp(pfam1_ptr->pfam_id,pfam2_ptr->pfam_id);

              /* If in wrong order, then swap sort-pointers */
              if (match > 0 || (match == 0 &&
                                pfam1_ptr->copy == TRUE &&
                                pfam2_ptr->copy == FALSE))
                {
                  swap_pfam_ptr = pfam_index_ptr[indi];
                  pfam_index_ptr[indi] = pfam_index_ptr[indl];
                  pfam_index_ptr[indl] = swap_pfam_ptr;

                  i = i - m;
                  if (i < 1)
                    endloop = TRUE;
                }
              else
                endloop = TRUE;
            }
        }
    }
}
/***********************************************************************

sort_pfam_records  -  Sort the Pfam domains into search order

***********************************************************************/

void sort_pfam_records(struct pfam **fst_pfam_ptr,int npfam,
                       struct parameters *params)
{
  int iorder, nentry;

  struct pfam *pfam_ptr, *first_pfam_ptr, *last_pfam_ptr;

  struct pfam **pfam_index_ptr;

  /* Initialise variables */
  first_pfam_ptr = *fst_pfam_ptr;
  last_pfam_ptr = NULL;

  /* Create the pfam array as one large array */
  pfam_index_ptr = (struct pfam **) malloc(npfam * sizeof(struct pfam *));

  /* Check that sufficient memory was allocated */
  if (pfam_index_ptr == NULL)
    write_error_message(MEMORY_ALLOCATION,
                        "Cannot allocate memory for pfam array",
                        params->output_type);

  /* Get pointer to the first code entry */
  pfam_ptr = first_pfam_ptr;
  nentry = 0;

  /* Loop through all entries to place each one in the pfam array */
  while (pfam_ptr != NULL)
    {
      /* Store current pointer at end of the pfam array */
      if (nentry < npfam)
        pfam_index_ptr[nentry] = pfam_ptr;

      /* Increment count of stored entries */
      nentry++;

      /* Get pointer to the next entry in the linked list */
      pfam_ptr = pfam_ptr->next_pfam_ptr;
    }

  /* Sort the code pointers into order of primary pfam */
  pfam_sort(pfam_index_ptr,npfam);

  /* Having sorted the records, rearrange order of the pointers */
  for (iorder = 0; iorder < npfam; iorder++)
    {
      /* Get the current pfam record */
      pfam_ptr = pfam_index_ptr[iorder];

      /* If this is the first, then make the first record */
      if (iorder == 0)
        first_pfam_ptr = pfam_ptr;

      /* Otherwise, get the previous pfam record to point to this one */
      else
        last_pfam_ptr->next_pfam_ptr = pfam_ptr;

      /* Blank out current record's next pfam pointer */
      pfam_ptr->next_pfam_ptr = NULL;

      /* Store current record as last encountered */
      last_pfam_ptr = pfam_ptr;
    }

  /* Free up the index array created here */
  free(pfam_index_ptr);

  /* Return the first Pfam record */
  *fst_pfam_ptr = first_pfam_ptr;
}
/***********************************************************************

find_uniprot  -  Find the architecture record(s) for the given UniProt
                 code

***********************************************************************/

int find_uniprot(int type,int ikey,char *code,
                 char *archschema_dir,struct parameters *params)
{
  int nhits;

  long offset, noffset;

  struct offset *offset_ptr, *first_offset_ptr;
  struct hit *first_hit_ptr;
  struct pfam *last_pfam_ptr;

  /* Initialise variables */
  nhits = 0;
  first_hit_ptr = NULL;
  last_pfam_ptr = NULL;

  /* Read in the starting offset for the .dat file */
  offset = read_uniprot_idx(type,ikey,code,archschema_dir,params);

  /* Use the offset value to read the corresponding .dat file to pick
     up the offset(s) in the archstruc.uniprot.joined.by.accession file */
  noffset
    = read_uniprot_dat(type,ikey,code,offset,&first_offset_ptr,
                       archschema_dir,params);

  /* If no offsets found, then can't find the sequence */
  if (noffset == 0)
    return(nhits);

  /* Retrieve the architecture record(s) for the current UniProt seq */
  nhits = get_seq_architecture(type,first_offset_ptr,&first_hit_ptr,
                               archschema_dir,params);

  /* If no hits, then return */
  if (nhits == 0)
    return(nhits);

  /* If we have more than one hit, determine which to use */
  if (nhits > 1)
    select_hit(&first_hit_ptr);

  /* Get the architecture for this hit */
  if (params->npfam == 0)
    {
      /* Split up the architecture into its constituent
         Pfam domains */
      params->npfam
        = split_architecture(type,&(params->first_pfam_ptr),
                             &last_pfam_ptr,first_hit_ptr->architecture,
                             params->npfam,params);
    }

  /* Return the number of hits */
  return(nhits);
}
/***********************************************************************

uniprot_search  -  Find the architecture record of the given UniProt
                   sequence

***********************************************************************/

int uniprot_search(int type,char *archschema_dir,
                   struct parameters *params)
{
  char search_code[UNIPROT_LEN + 1], err_message[LINELEN];

  int ikey, nhits;

  /* Initialise variables */
  nhits = 0;

  /* Search according to whether we have the UniProt accession or the
     UniProt id */
  if (params->uniprot_id != &null)
    {
      ikey = KEY_UNIPROT_ID;
      strcpy(search_code,params->uniprot_id);
    }
  else
    {
      ikey = KEY_UNIPROT_ACC;
      strcpy(search_code,params->uniprot_acc);
    }

  /* Perform the search */
  nhits = find_uniprot(type,ikey,search_code,archschema_dir,params);

  /* If search code not found, show error message */
  if (nhits == 0)
    {
      sprintf(err_message,"UniProt code %s not found in database!",
              search_code);
      write_error_message(CODE_NOT_FOUND,err_message,params->output_type);
    }
    
  /* Return number of hits */
  return(nhits);
}
/***********************************************************************

add_split_domains  -  Add the split CATH domain to current list of
                      domains

***********************************************************************/

void add_split_domains(struct parameters *params)
{
  char pfam_id[PFAM_LEN + 1];

  int ipfam, npfam;
  int copy, found;

  struct pfam *pfam_ptr, *other_pfam_ptr, *last_pfam_ptr;

  /* Initialise variables */
  copy = FALSE;
  ipfam = 0;
  last_pfam_ptr = NULL;

  /* Get pointer to the first pfam */
  pfam_ptr = params->first_pfam_ptr;

  /* Loop over all the pfams to stroe pointer to last one */
  while (pfam_ptr != NULL && ipfam < params->npfam)
    {
      /* Save this pointer */
      last_pfam_ptr = pfam_ptr;

      /* Increment domain count */
      ipfam++;

      /* Get pointer to the next pfam */
      pfam_ptr = pfam_ptr->next_pfam_ptr;
    }

  /* Re-initialise counter */
  ipfam = 0;
  npfam = params->npfam;

  /* Get pointer to the first pfam */
  pfam_ptr = params->first_pfam_ptr;

  /* Loop over all the pfams to add split-domain versions */
  while (pfam_ptr != NULL && ipfam < npfam)
    {
      /* If a split-domain, the form unsplit version */
      if (pfam_ptr->pfam_id[0] == 'p')
        strcpy(pfam_id,pfam_ptr->pfam_id + 1);

      /* Otherwise, if not a split domain, form split-domain version */
      else
        {
          strcpy(pfam_id,"p");
          strcat(pfam_id,pfam_ptr->pfam_id);
        }

      /* Initialise flag in search to see if we already have this domain */
      found = FALSE;

      /* Get pointer to the first record pfam */
      other_pfam_ptr = params->first_pfam_ptr;

      /* Loop over all the pfams to see if we already ahve this one */
      while (other_pfam_ptr != NULL && found == FALSE)
        {
          /* If we already have this domain, set flag */
          if (!strcmp(other_pfam_ptr->pfam_id,pfam_id))
            found = TRUE;

          /* Get pointer to the next pfam */
          other_pfam_ptr = other_pfam_ptr->next_pfam_ptr;
        }

      /* If domain not found, then add it */
      if (found == FALSE)
        {
          /* Create a new Pfam record */
          other_pfam_ptr
            = create_pfam_record(&(params->first_pfam_ptr),
                                 &last_pfam_ptr,pfam_id,copy,params);

          /* Mark it as an extra domain */
          other_pfam_ptr->extra = TRUE;

          /* Increment count of Pfam domains to search on */
          params->npfam++;
        }

      /* Increment domain count */
      ipfam++;

      /* Get pointer to the next pfam */
      pfam_ptr = pfam_ptr->next_pfam_ptr;
    }
}
/***********************************************************************

read_archpfam_idx  -  Read in the offsets from the archpfam.idx file for
                      all the Pfam domains to be searched

***********************************************************************/

int read_archpfam_idx(int type,char *archschema_dir,
                      struct parameters *params)
{
  char file_name[FILENAME_LEN], input_line[LINELEN + 1];
  char number_string[20], pfam_id[PFAM_LEN + 1], message[LINELEN + 1];

  char *string_ptr;

  int len, nhits, nlines;
  int done, next_code, next_record;

  long offset;

  struct pfam *pfam_ptr;

  FILE *file_ptr;

  /* Initialise variables */
  done = FALSE;
  next_code = TRUE;
  next_record = TRUE;
  nlines = nhits = 0;
  offset = -1;
  pfam_id[0] = '\0';

  /* Form name of input file */
  strcpy(file_name,archschema_dir);
  strcat(file_name,"/");
  if (type == PFAM)
    strcat(file_name,"archpfam.idx");
  else
    strcat(file_name,"archcath.idx");

  /* Get pointer to the first pfam record */
  pfam_ptr = params->first_pfam_ptr;
  next_code = FALSE;
  if (pfam_ptr == NULL)
    return(nhits);

  /* Open the index file, archpfam.idx */
  if ((file_ptr = fopen(file_name,"r")) == NULL)
    {
      /* Otherwise, error must be due to missing file */
      sprintf(message,"Unable to open data file [%s]\n",file_name);
      write_error_message(FILE_NOT_FOUND,message,params->output_type);
    }

  /* Loop until all codes have been assigned a protein name */
  while (done == FALSE)
    {
      /* If need to go to the next code, then do so */
      if (next_code == TRUE)
        {
          /* Get pointer to the next pfam record */
          if (pfam_ptr != NULL)
            pfam_ptr = pfam_ptr->next_pfam_ptr;
          else
            done = TRUE;

          /* Check for end */
          if (pfam_ptr == NULL)
            done = TRUE;

          /* Set flag */
          next_code = FALSE;
        }

      /* If need to read to next record, then loop until we have the
         next useful one */
      while (done == FALSE && next_record == TRUE)
        {
          /* Read in the next record from the file */
          if (fgets(input_line,LINELEN,file_ptr) != NULL)
            {
              /* Increment line count */
              nlines++;

              /* Get the Pfam id */
              pfam_id[0] = '\0';
              string_ptr = strchr(input_line,' ');

              /* If have a space, extract the Pfam id */
              if (string_ptr != NULL)
                {
                  /* Get the length of the code */
                  len = string_ptr - input_line;

                  /* Check that Pfam id not too long */
                  if (len < PFAM_LEN)
                    {
                      strncpy(pfam_id,input_line,len);
                      pfam_id[len] = '\0';

                      /* Truncate the string */
                      len = string_chop(string_ptr + 1);

                      /* Get the offset */
                      strcpy(number_string,string_ptr + 1);
                      offset = atol(number_string);
                    }
                }

              /* If this is less than our current Pfam record,
                 need to keep reading */
              if (strcmp(pfam_id,pfam_ptr->pfam_id) < 0)
                next_record = TRUE;
              else
                next_record = FALSE;
            }

          /* Have reached the end of the file, so we're done */
          else
            done = TRUE;
        }

      /* If not yet done, check the match */
      if (done == FALSE)
        {
          /* If this matches our current Pfam record, store the offset */
          if (!strcmp(pfam_id,pfam_ptr->pfam_id))
            {
              /* Store the offset */
              pfam_ptr->offset = offset;

              /* Increment count of matches */
              nhits++;

              /* Set both flags to go to the next record */
              next_code = TRUE;
              next_record = TRUE;
            }

          /* Otherwise, if it is higher, need to get the next code */
          else if (strcmp(input_line,pfam_ptr->pfam_id) > 0)
            {
              /* Set flags */
              next_code = TRUE;
              next_record = FALSE;
            }

          /* Otherwise, must be lower, so need to read the next
             record */
          else
            next_record = TRUE;
        }
    }

  /* Return number of Pfam ids matched */
  return(nhits);
}
/***********************************************************************

read_archpfam_dat  -  Read in the archpfam.dat offsets from the
                      archpfam.idx file for all the Pfam domains to be
                      searched

***********************************************************************/

int read_archpfam_dat(int type,char *archschema_dir,
                      struct parameters *params)
{
  char file_name[FILENAME_LEN], input_line[LINELEN + 1];
  char pfam_id[PFAM_LEN + 1], message[LINELEN + 1];

  char *string_ptr;

  int len, line, nhits;
  int done, error;

  long offset;

  struct offset *offset_ptr, *last_offset_ptr;
  struct pfam *pfam_ptr;

  FILE *file_ptr;

  /* Initialise variables */
  nhits = 0;
  pfam_id[0] = '\0';

  /* Form name of input file */
  strcpy(file_name,archschema_dir);
  strcat(file_name,"/");
  if (type == PFAM)
    strcat(file_name,"archpfam.dat");
  else
    strcat(file_name,"archcath.dat");

  /* Open the index file, archpfam.idx */
  if ((file_ptr = fopen(file_name,"r")) == NULL)
    {
      /* Otherwise, error must be due to missing file */
      sprintf(message,"Unable to open data file [%s]\n",file_name);
      write_error_message(FILE_NOT_FOUND,message,params->output_type);
    }

  /* Get pointer to the first pfam record */
  pfam_ptr = params->first_pfam_ptr;

  /* Loop over the Pfam records */
  while (pfam_ptr != NULL)
    {
      /* Initialise variables */
      last_offset_ptr = NULL;
      
      /* Get the offset from this record */
      offset = pfam_ptr->offset;

      /* If valid, then find this point in the file */
      if (offset > -1)
        {
          /* Find start-point in data file using fseek */
          error = fseek(file_ptr,offset,SEEK_SET);

          /* Initialise parameters for this read */
          done = FALSE;
          line = 0;

          /* Loop while reading in records from the file */
          while (fgets(input_line,LINELEN,file_ptr) != NULL &&
                 done == FALSE)
            {
              /* Truncate the string */
              len = string_chop(input_line);

              /* If this is the first line, check that the Pfam id
                 matches */
              if (line == 0)
                {
                  /* If have a match, then increment match count */
                  if (!strcmp(input_line + 1,pfam_ptr->pfam_id))
                    nhits++;
                  else
                    done = TRUE;
                }

              /* Otherwise, check for start of next Pfam record */
              else if (input_line[0] == ':')
                done = TRUE;

              /* Otherwise, get this offset */
              else
                {
                  /* Get offset value */
                  offset = atol(input_line);

                  /* If one of the first two, store in Pfam record */
                  //                  if (line < 3)
                  //                    pfam_ptr->offset[line - 1] = offset;

                  /* Otherwise, need to create an offset record */
                  //                  else
                  //                    {
                  /* Create an offset record to store the current
                     offset in the
                     archstruc.uniprot.joined.by.accession file */
                  create_offset_record(&(pfam_ptr->first_offset_ptr),
                                       &last_offset_ptr,offset,params);
                      //                    }
                }

              /* Increment line counter */
              line++;
            }
        }

      /* Get the next Pfam record */
      pfam_ptr = pfam_ptr->next_pfam_ptr;
    }

  /* Close the file */
  fclose(file_ptr);

  /* Return the number of records matched */
  return(nhits);
}
/***********************************************************************

check_selection  -  Check whether the given hit satisfies the selection
                    criteria

***********************************************************************/

int check_selection(int type,char *uniprot_id,char *uniprot_acc,
                    char *review_flag,char *architecture,char *pdb_codes,
                    int *is_parent,struct parameters *params)
{
  char species[SPECIES_LEN + 3];

  char *tmp, *string_ptr;

  int len, nfound, npfam;
  int wanted;

  struct filterseq *filterseq_ptr;
  struct pfam *pfam_ptr;

  /* Initialise variables */
  *is_parent = FALSE;
  wanted = TRUE;

  /* If this is the parent sequence, then keep irrespective of the other
     selection criteria */
  if (!strcmp(uniprot_id,params->uniprot_id) ||
      !strcmp(uniprot_acc,params->uniprot_acc))
    {
      /* Set flag that this is the parent sequence */
      *is_parent = TRUE;
      return(TRUE);
    }

  /* If require reviewed UniProt entries only, then check the review flag */
  if (params->reviewed_only == TRUE && strcmp(review_flag,"Reviewed"))
    return(FALSE);

  /* If require a specific species, then check the UniProt id matches */
  if (params->species != &null)
    {
      /* Get the species from this Uniprot id */
      string_ptr = strchr(uniprot_id,'_');

      /* If have a species name, then prepare to search against our list */
      if (string_ptr != NULL)
        {
          /* Form the species search string */
          strcpy(species,"|");
          strcat(species,string_ptr + 1);
          strcat(species,"|");

          /* Search to the species in the list */
          string_ptr = strstr(params->species,species);

          /* If not found then don't want this sequence */
          if (string_ptr == NULL)
            return(FALSE);
        }

      /* If we don't have a species name, then reject */
      else
        return(FALSE);
    }

  /* If structure-only hit resquired, check that this one has PDB codes
     assigned */
  if (params->pdb_only == TRUE && !strcmp(pdb_codes,"NONE"))
    return(FALSE);

  /* If domain operator is AND, need all the domains to match */
  if (params->operator == AND)
    {
      /* Take a copy of the Pfam domain architecture */
      len = strlen(architecture);

      /* Create a temporary string */
      tmp = (char *) malloc(sizeof(char)*(len + 1));

      /* Get pointer to the first pfam record */
      pfam_ptr = params->first_pfam_ptr;

      /* Loop over the Pfam records to initialise selection flag */
      while (pfam_ptr != NULL && wanted == TRUE)
        {
          /* Check to see if this Pfam domain is present in the
             architecture */
          string_ptr = strstr(architecture,pfam_ptr->pfam_id);

          /* If found, knock out the domain so it isn't found a 
             second time */
          if (string_ptr != NULL)
            string_ptr[0] = DOM_SEP[type][0];

          /* Otherwise, hit is not wanted as it doesn't contain all
             our Pfam domains */
          else
            wanted = FALSE;

          /* Get the next Pfam record */
          pfam_ptr = pfam_ptr->next_pfam_ptr;
        }

      /* Free up the memory */
      free(tmp);
    }

  /* Otherwise, for the OR operator, only need one domain to match */
  else
    {
      /* Assume hit is not wanted */
      wanted = FALSE;

      /* Get pointer to the first pfam record */
      pfam_ptr = params->first_pfam_ptr;

      /* Loop over the Pfam records to initialise selection flag */
      while (pfam_ptr != NULL && wanted == FALSE)
        {
          /* If this domain is in our architecture, then architecture
             is wanted */
          if (strstr(architecture,pfam_ptr->pfam_id) != NULL)
            wanted = TRUE;

          /* Get the next Pfam record */
          pfam_ptr = pfam_ptr->next_pfam_ptr;
        }
    }

  /* If this hit is wanted, but have a list of filter sequences, check
     that this sequence is on that list */
  if (wanted == TRUE && params->first_filterseq_ptr != NULL)
    {
      /* Assume that sequence is not on the list */
      wanted = FALSE;

        /* Get pointer to the first sequence record */
      filterseq_ptr = params->first_filterseq_ptr;

      /* Loop over the sequence records to see if this one is amongst
         them */
      while (filterseq_ptr != NULL && wanted == FALSE)
        {
          /* Check if the current hit's UniProt code or accession code
             matches the one from the list */
          if (!strcmp(uniprot_id,filterseq_ptr->seq_id) ||
              !strcmp(uniprot_acc,filterseq_ptr->seq_id))
            wanted = TRUE;

          /* Get the next record */
          filterseq_ptr = filterseq_ptr->next_filterseq_ptr;
        }
    }

  /* Return whether the hit is wanted */
  return(wanted);
}
/***********************************************************************

check_previous  -  Check previous Pfam domains to see if architecture
                   contains one of them (and hence should already have
                   been picked up)

***********************************************************************/

int check_previous(int type,struct pfam *first_pfam_ptr,
                   struct pfam *pfam_ptr,char *architecture)
{
  char *tmp, *token[MAX_TOKENS];

  int itoken, len, ntokens;
  int wanted;

  struct pfam *other_pfam_ptr;

  /* Initialise variables */
  wanted = TRUE;

  /* Get the length of the architecture */
  len = strlen(architecture);

  /* Allocate temporary memory for the architecture */
  tmp = (char *) malloc (sizeof(char)*(len + 1));

  /* Store the architecture */
  strcpy(tmp,architecture);

  /* Extract the dot-separated Pfam domains from the architecture */
  ntokens = tokenize(tmp,token,MAX_TOKENS,DOM_SEP[type][0]);

  /* Get pointer to the first pfam record */
  other_pfam_ptr = first_pfam_ptr;

  /* Loop over the Pfam records to see if any are in this architecture */
  while (other_pfam_ptr != NULL && other_pfam_ptr != pfam_ptr &&
         wanted == TRUE)
    {
      /* Loop over the Pfam ids, to see if we have current one */
      for (itoken = 0; itoken < ntokens; itoken++)
        {
          /* If have a match, then current architecture not wanted */
          if (!strcmp(other_pfam_ptr->pfam_id,token[itoken]))
            wanted = FALSE;
        }

      /* Get the next Pfam record */
      other_pfam_ptr = other_pfam_ptr->next_pfam_ptr;
    }

  /* Free up memory taken by the architecture */
  free(tmp);

  /* Return whether the hit is wanted */
  return(wanted);
}
/***********************************************************************

create_species_record  -  Create and initialise a new species record

***********************************************************************/

struct species *create_species_record(struct species **fst_species_ptr,
                                      struct species **lst_species_ptr,
                                      char *species_id,
                                      struct parameters *params)
{
  struct species *species_ptr, *first_species_ptr, *last_species_ptr;

  /* Initialise species pointers */
  species_ptr = NULL;
  first_species_ptr = *fst_species_ptr;
  last_species_ptr = *lst_species_ptr;

  /* Allocate memory for structure to hold species info */
  species_ptr = (struct species *) malloc(sizeof(struct species));
  if (species_ptr == NULL)
    write_error_message(MEMORY_ALLOCATION,
                        "Cannot allocate memory for struct species",
                        params->output_type);

  /* If this is the very first, save its pointer */
  if (first_species_ptr == NULL)
    first_species_ptr = species_ptr;

  /* Add link from previous species to the current one */
  if (last_species_ptr != NULL)
    last_species_ptr->next_species_ptr = species_ptr;
  last_species_ptr = species_ptr;

  /* Store known elements of the structure */
  store_string(&(species_ptr->species_id),species_id);

  /* Initialise remaining fields */
  species_ptr->species_name = &null;
  species_ptr->nseqs = 1;

  /* Initialise pointer to next species record */
  species_ptr->next_species_ptr = NULL;

  /* Return current pointers */
  *fst_species_ptr = first_species_ptr;
  *lst_species_ptr = last_species_ptr;

  /* Return new species pointer */
  return(species_ptr);
}
/***********************************************************************

get_architectures  -  If have any offsets, pick up the data from the
                      archstruc.uniprot.joined.by.accession file

***********************************************************************/

int get_architectures(int type,struct hit **fst_hit_ptr,
                      struct hit **prnt_hit_ptr,
                      struct species **fst_species_ptr,int *nsp,
                      char *archschema_dir,
                      struct parameters *params)
{
  char file_name[FILENAME_LEN], input_line[LINELEN + 1];
  char message[LINELEN + 1];

  char *string_ptr, *token[MAX_TOKENS];

  int itoken, len, line, nhits, nspecies, ntokens;
  int done, error, is_parent, wanted;

  long offset;

  struct hit *hit_ptr, *first_hit_ptr, *last_hit_ptr;
  struct hit *parent_hit_ptr;
  struct offset *offset_ptr;
  struct pfam *pfam_ptr;
  struct species *species_ptr, *first_species_ptr, *last_species_ptr;

  FILE *file_ptr;

  /* Initialise variables */
  nhits = 0;
  *nsp = nspecies = 0;

  /* Initialise pointers */
  hit_ptr = first_hit_ptr = last_hit_ptr = NULL;
  *prnt_hit_ptr = parent_hit_ptr = NULL;
  species_ptr = first_species_ptr = last_species_ptr = NULL;

  /* Form name of input file */
  strcpy(file_name,archschema_dir);
  strcat(file_name,"/");
  if (type == PFAM)
    strcat(file_name,"archstruc.uniprot.joined.by.accession");
  else
    strcat(file_name,"cathstruc.uniprot.joined.by.accession");

  /* Open the input file */
  if ((file_ptr = fopen(file_name,"r")) == NULL)
    {
      /* Otherwise, error must be due to missing file */
      sprintf(message,"Unable to open data file [%s]\n",file_name);
      write_error_message(FILE_NOT_FOUND,message,params->output_type);
    }

  /* Get pointer to the first pfam record */
  pfam_ptr = params->first_pfam_ptr;

  /* Loop over the Pfam records */
  while (pfam_ptr != NULL)
    {
      /* Get pointer to the first offset record */
      offset_ptr = pfam_ptr->first_offset_ptr;

      /* Loop over the offset records */
      while (pfam_ptr->copy == FALSE && offset_ptr != NULL)
        {
          /* Get the offset from this record */
          offset = offset_ptr->offset;

          /* If valid, then find this point in the file */
          if (offset > -1)
            {
              /* Find start-point in data file using fseek */
              error = fseek(file_ptr,offset,SEEK_SET);

              /* Initialise parameters for this read */
              done = FALSE;
              line = 0;

              /* Read in this record */
              if (fgets(input_line,LINELEN,file_ptr) != NULL)
                {
                  /* Truncate the string */
                  len = string_chop(input_line);

                  /* Get the tokens on this line */
                  ntokens = tokenize(input_line,token,MAX_TOKENS,'\t');

                  /* Check whether this record satisfies the selection
                     criteria */
                  wanted
                    = check_selection(type,token[F_UNIPROT_ID],
                                      token[F_UNIPROT_ACC],
                                      token[F_REVIEW_FLAG],
                                      token[F_ARCHITECTURE],
                                      token[F_PDB_CODES],&is_parent,
                                      params);

                  /* If this is not the parent, and is not the first
                     Pfam domain, check whether architecture contains
                     a domain we have already searched on (and so is
                     one we already have) */
                  /* if (wanted == TRUE && is_parent == FALSE &&
                     pfam_ptr != params->first_pfam_ptr) */
                  if (wanted == TRUE)
                    wanted
                      = check_previous(type,params->first_pfam_ptr,
                                       pfam_ptr,token[F_ARCHITECTURE]);

                  /* If this hit is wanted, save the data */
                  if (wanted == TRUE)
                    {
                      hit_ptr
                        = create_hit_record(&first_hit_ptr,
                                            &last_hit_ptr,
                                            token[F_UNIPROT_ACC],
                                            token[F_UNIPROT_ID],
                                            token[F_ARCHITECTURE],
                                            token[F_GENE],
                                            token[F_PROTEIN_NAME],
                                            atoi(token[F_SEQ_LENGTH]),
                                            token[F_PDB_CODES],
                                            token[F_COVERAGE],params);

                      /* Increment count of hits */
                      nhits++;

                      /* If this is the parent sequence, save its
                         pointer */
                      if (is_parent == TRUE)
                        parent_hit_ptr = hit_ptr;

                      /* Get the species identifier for this sequence */
                      string_ptr = strchr(hit_ptr->uniprot_id,'_');
                      if (string_ptr != NULL)
                        {
                          /* Create a species record */
                          species_ptr
                            = create_species_record(&first_species_ptr,
                                                    &last_species_ptr,
                                                    string_ptr + 1,params);

                          /* Increment count of species */
                          nspecies++;
                        }
                    }

                  /* Increment line counter */
                  line++;
                }
            }

          /* Get the next offset record */
          offset_ptr = offset_ptr->next_offset_ptr;
        }

      /* Get the next Pfam record */
      pfam_ptr = pfam_ptr->next_pfam_ptr;
    }

  /* Close the file */
  fclose(file_ptr);

  /* If combining architectures, create a dummy species record */
  if (params->combine == TRUE)
    {
      /* Create the dummy species record */
      species_ptr
        = create_species_record(&first_species_ptr,&last_species_ptr,
                                "NONE",params);

      /* Increment count of species */
      nspecies++;

      /* Store the name */
      store_string(&(species_ptr->species_name),"No species");
    }

  /* Return pointers */
  *fst_hit_ptr = first_hit_ptr;
  *fst_species_ptr = first_species_ptr;
  *prnt_hit_ptr = parent_hit_ptr;

  /* Return count of species records created */
  *nsp = nspecies;

  /* Return the number of records matched */
  return(nhits);
}
/***********************************************************************

choose_parent  -  Select one of returned sequences to act as the parent
                  sequence

***********************************************************************/

struct hit *choose_parent(int type,struct hit *first_hit_ptr,int *nh,
                          struct parameters *params)
{
  char *architecture;
  
  int best_score, ikey, ipfam, len, ndomains, nhits, score, sscore;
  int have_full;

  struct hit *hit_ptr, *last_hit_ptr, *parent_hit_ptr;
  struct pfam *pfam_ptr;

  /* Initialise variables */
  best_score = ipfam = ndomains = score = 0;
  have_full = FALSE;
  last_hit_ptr = NULL;
  nhits = *nh;

  /* Determine length of parent architecture */
  len = params->npfam * (PFAM_LEN + 1);

  /* Allocate temporary memory for the architecture */
  architecture = (char *) malloc (sizeof(char)*(len + 1));
  architecture[0] = '\0';

  /* Get pointer to the first pfam */
  pfam_ptr = params->first_pfam_ptr;

  /* Loop over all the pfams to form the architecture */
  while (pfam_ptr != NULL && ipfam < params->npfam)
    {
      /* Add this domain to the architecture */
      if (pfam_ptr->extra == FALSE)
        {
          if (ipfam > 0)
            strcat(architecture,DOM_SEP[type]);
          strcat(architecture,pfam_ptr->pfam_id);

          /* Increment count of domains */
          ndomains++;
        }

      /* Increment domain count */
      ipfam++;

      /* Get pointer to the next pfam */
      pfam_ptr = pfam_ptr->next_pfam_ptr;
    }

  /* Get the first hit record */
  parent_hit_ptr = hit_ptr = first_hit_ptr;

  /* Loop over the hit records to write out */
  while (hit_ptr != NULL)
    {
      /* Save this hit as the last so far */
      last_hit_ptr = hit_ptr;

      /* If architecture matches the parent architecture exactly,
         check the species */
      if (!strcmp(architecture,hit_ptr->architecture))
        {
          score = 100;
          have_full = TRUE;
        }

      /* If hit's architecture contains the parent architecture, then
         give lower score */
      else if (strstr(hit_ptr->architecture,architecture) != NULL )
        {
          score = 20;
          have_full = TRUE;
        }

      /* If searching for a combined architecure, search for each domain
         separately */
      else if (params->combine == TRUE)
        {
          /* Initialise score */
          ipfam = score = 0;

          /* Get pointer to the first pfam */
          pfam_ptr = params->first_pfam_ptr;

          /* Loop over all the pfams to form the architecture */
          while (pfam_ptr != NULL && ipfam < params->npfam)
            {
              /* If not an added split domain, check if it is in the
                 current architecture */
              if (pfam_ptr->extra == FALSE)
                {
                  if (strstr(hit_ptr->architecture,pfam_ptr->pfam_id) != NULL)
                    score++;
                }

              /* Increment domain count */
              ipfam++;

              /* Get pointer to the next pfam */
              pfam_ptr = pfam_ptr->next_pfam_ptr;
            }

          /* If have all the domains, then set full flag */
          if (score == ndomains)
            have_full = TRUE;
        }
      
      /* Otherwise, score zero */
      else
        score = 0;

      /* If have a non-zero score, assess according to how "interesting"
         the species is */
      if (score > 0)
        {
          /* Set default species score */
          sscore = NKEY_SPECIES;

          /* Loop over the key species and score accordingly */
          for (ikey = 0; ikey < NKEY_SPECIES; ikey++)
            {
              /* If have a match, then store key position as the match
                 score */
              if (strstr(hit_ptr->uniprot_id,KEY_SPECIES[ikey]) != NULL)
                sscore = ikey;
            }

          /* Subtract the species score from the overall score */
          score = score - sscore;
        }

      /* If this is the best score so far, save this hit as a potential
         parent */
      if (score > best_score)
        {
          parent_hit_ptr = hit_ptr;
          best_score = score;
        }

      /* Get the next hit record */
      hit_ptr = hit_ptr->next_hit_ptr;
    }

  /* If combining and all domains weren't found, create a dummy
     architecture and make it the parent */
  if (params->combine == TRUE && have_full == FALSE)
    {
      /* Create a dummy hit record */
      hit_ptr = create_hit_record(&first_hit_ptr,&last_hit_ptr,
                                  "NONE","NONE_NONE",architecture,"\0",
                                  "\0",0,"NONE","NONE",params);

      /* Increment hit count */
      nhits++;

      /* Save this as the parent architecture */
      parent_hit_ptr = hit_ptr;
    }

  /* Free up the memory allocated in routine */
  free(architecture);

  /* Return number of hits */
  *nh = nhits;

  /* Return the parent hit, if we have one */
  return(parent_hit_ptr);
}
/***********************************************************************

hit_sort  -  Shell-sort to arrange the hits by architecture

***********************************************************************/

void hit_sort(struct hit **hit_index_ptr,int nhits)
{
  long endloop, i, indi, indl, j, k, l, lognb2, m, n, nn;
  float aln2i = (float) 1.4426950, tiny = (float) 1.e-5;
  float calc;
  double dble, dble1;

  struct hit *hit1_ptr, *hit2_ptr, *swap_hit_ptr;

  /*--Initialise values */
  n = nhits;
  dble = n;
  dble1 = log(dble);
  calc = (float) dble1;
  calc = calc * aln2i + tiny;
  lognb2 = (long) calc;

  /* Perform sort */
  m = n;
  for (nn = 1; nn <= lognb2; nn++)
    {
      m = m / 2;
      k = n - m;
      for (j = 1; j <= k; j++)
        {
          i = j;
          endloop = FALSE;
          while (endloop == FALSE)
            {
              l = i + m;
              indi = i - 1;
              indl = l - 1;

              /* Get the two code records */
              hit1_ptr = hit_index_ptr[indi];
              hit2_ptr = hit_index_ptr[indl];

              /* If in wrong order, then swap sort-pointers */
              if (strcmp(hit1_ptr->architecture,hit2_ptr->architecture) > 0)
                {
                  swap_hit_ptr = hit_index_ptr[indi];
                  hit_index_ptr[indi] = hit_index_ptr[indl];
                  hit_index_ptr[indl] = swap_hit_ptr;

                  i = i - m;
                  if (i < 1)
                    endloop = TRUE;
                }
              else
                endloop = TRUE;
            }
        }
    }
}
/***********************************************************************

sort_hit_records  -  Sort the hit records by architecture

***********************************************************************/

void sort_hit_records(struct hit *first_hit_ptr,
                      struct hit ***hit_ind_ptr,int nhits,
                      struct parameters *params)
{
  int nentry;

  struct hit *hit_ptr;

  struct hit **hit_index_ptr;

  /* Create the pfam array as one large array */
  hit_index_ptr = (struct hit **) malloc(nhits * sizeof(struct hit *));

  /* Check that sufficient memory was allocated */
  if (hit_index_ptr == NULL)
    write_error_message(MEMORY_ALLOCATION,
                        "Cannot allocate memory for hits array",
                        params->output_type);

  /* Get pointer to the first code entry */
  hit_ptr = first_hit_ptr;
  nentry = 0;

  /* Loop through all entries to place each one in the array */
  while (hit_ptr != NULL)
    {
      /* Store current pointer at end of array */
      if (nentry < nhits)
        hit_index_ptr[nentry] = hit_ptr;

      /* Increment count of stored entries */
      nentry++;

      /* Get pointer to the next entry in the linked list */
      hit_ptr = hit_ptr->next_hit_ptr;
    }

  /* Sort the hit pointers into order of architecture */
  if (nhits > 1)
    hit_sort(hit_index_ptr,nhits);

  /* Return the array pointer */
  *hit_ind_ptr = hit_index_ptr;
}
/***********************************************************************

create_archit_record  -  Create and initialise a new architecture record

***********************************************************************/

struct archit *create_archit_record(struct archit **fst_archit_ptr,
                                    struct archit **lst_archit_ptr,
                                    char *architecture,
                                    struct parameters *params)
{
  struct archit *archit_ptr, *first_archit_ptr, *last_archit_ptr;

  /* Initialise pointers */
  archit_ptr = NULL;
  first_archit_ptr = *fst_archit_ptr;
  last_archit_ptr = *lst_archit_ptr;

  /* Allocate memory for structure */
  archit_ptr = (struct archit *) malloc(sizeof(struct archit));
  if (archit_ptr == NULL)
    write_error_message(MEMORY_ALLOCATION,
                        "Cannot allocate memory for struct archit",
                        params->output_type);

  /* If this is the very first, save its pointer */
  if (first_archit_ptr == NULL)
    first_archit_ptr = archit_ptr;

  /* Add link from previous pfam to the current one */
  if (last_archit_ptr != NULL)
    last_archit_ptr->next_archit_ptr = archit_ptr;
  last_archit_ptr = archit_ptr;

  /* Store known elements of the structure */
  store_string(&(archit_ptr->architecture),architecture);

  /* Initialise remaining fields */
  archit_ptr->nhits = 0;
  archit_ptr->score = 0;
  archit_ptr->first_enzlist_ptr = NULL;
  archit_ptr->first_hitlist_ptr = NULL;

  /* Initialise pointer to next record */
  archit_ptr->next_archit_ptr = NULL;

  /* Return current pointers */
  *fst_archit_ptr = first_archit_ptr;
  *lst_archit_ptr = last_archit_ptr;

  /* Return new pointer */
  return(archit_ptr);
}
/***********************************************************************

create_hitlist_record  -  Create and initialise a new hitlist record

***********************************************************************/

struct hitlist *create_hitlist_record(struct hitlist **fst_hitlist_ptr,
                                      struct hitlist **lst_hitlist_ptr,
                                      struct hit *hit_ptr,
                                      struct parameters *params)
{
  struct hitlist *hitlist_ptr, *first_hitlist_ptr, *last_hitlist_ptr;

  /* Initialise pointers */
  hitlist_ptr = NULL;
  first_hitlist_ptr = *fst_hitlist_ptr;
  last_hitlist_ptr = *lst_hitlist_ptr;

  /* Allocate memory for structure */
  hitlist_ptr = (struct hitlist *) malloc(sizeof(struct hitlist));
  if (hitlist_ptr == NULL)
    write_error_message(MEMORY_ALLOCATION,
                        "Cannot allocate memory for struct hitlist",
                        params->output_type);

  /* If this is the very first, save its pointer */
  if (first_hitlist_ptr == NULL)
    first_hitlist_ptr = hitlist_ptr;

  /* Add link from previous pfam to the current one */
  if (last_hitlist_ptr != NULL)
    last_hitlist_ptr->next_hitlist_ptr = hitlist_ptr;
  last_hitlist_ptr = hitlist_ptr;

  /* Store known elements of the structure */
  hitlist_ptr->hit_ptr = hit_ptr;

  /* Initialise pointer to next record */
  hitlist_ptr->next_hitlist_ptr = NULL;

  /* Return current pointers */
  *fst_hitlist_ptr = first_hitlist_ptr;
  *lst_hitlist_ptr = last_hitlist_ptr;

  /* Return new pointer */
  return(hitlist_ptr);
}
/***********************************************************************

extract_architectures  -  Create an architecture record for each unique
                          architecture

***********************************************************************/

int extract_architectures(struct hit **hit_index_ptr,int nhits,
                          struct hit *parent_hit_ptr,
                          struct archit **fst_archit_ptr,
                          struct parameters *params)
{
  int ihit, narchit;

  struct archit *archit_ptr, *first_archit_ptr, *last_archit_ptr;
  struct hit *hit_ptr;
  struct hitlist *hitlist_ptr, *last_hitlist_ptr;

  /* Initialise variables */
  narchit = 0;

  /* Initialise pointers */
  archit_ptr = first_archit_ptr = last_archit_ptr = NULL;
  last_hitlist_ptr = NULL;

  /* Loop over the sorted hit records */
  for (ihit = 0; ihit < nhits; ihit++)
    {
      /* Get this hit record */
      hit_ptr = hit_index_ptr[ihit];

      /* If this is the first hit, or the architecture is a new one,
         create a new architecture record */
      if (narchit == 0 ||
          strcmp(hit_ptr->architecture,last_archit_ptr->architecture))
        {
          /* Create a new architecture record */
          archit_ptr
            = create_archit_record(&first_archit_ptr,
                                   &last_archit_ptr,
                                   hit_ptr->architecture,params);

          /* Increment count of architecture records created */
          narchit++;

          /* Initialise last hitlist pointer */
          last_hitlist_ptr = NULL;
        }

      /* Attach the current hit to its corresponding architecture
         record */
      hitlist_ptr
        = create_hitlist_record(&(archit_ptr->first_hitlist_ptr),
                                &last_hitlist_ptr,hit_ptr,params);

      /* Increment count of this architecture's hits */
      archit_ptr->nhits++;
    }

  /* Return pointers */
  *fst_archit_ptr = first_archit_ptr;

  /* Return number of unique architectures */
  return(narchit);
}
/***********************************************************************

find_architecture  -  Find the given architecture entry

***********************************************************************/

struct archit *find_architecture(struct archit *first_archit_ptr,
                                 char *architecture)
{
  struct archit *archit_ptr, *found_archit_ptr;

  /* Initialise variables */
  found_archit_ptr = NULL;

  /* Get pointer to the first architecture */
  archit_ptr = first_archit_ptr;

  /* Loop over all the architectures to find the one of interest */
  while (archit_ptr != NULL && found_archit_ptr == NULL)
    {
      /* If this is the right architecture, store pointer */
      if (!strcmp(architecture,archit_ptr->architecture))
        found_archit_ptr = archit_ptr;

      /* Get pointer to the next architecture */
      archit_ptr = archit_ptr->next_archit_ptr;
    }

  /* Return the architecture pointer */
  return(found_archit_ptr);
}
/***********************************************************************

replace_pfamb_domains  -  Replace all the Pfam-B domnains in the
                          architecture records with PB000000

***********************************************************************/

int replace_pfamb_domains(int type,struct archit *first_archit_ptr,
                          char *parent_architecture)
{
  char *string_ptr;
  char *newtmp, *tmp, *token[MAX_TOKENS];

  int ipos, itoken, len, narchit, ndeleted, nseqs, ntokens;

  struct archit *archit_ptr;

  /* Initialise variables */
  narchit = 0;

  /* Get the first archit record */
  archit_ptr = first_archit_ptr;

  /* Loop over the archit records perform the collapse */
  while (archit_ptr != NULL)
    {
      /* Increment architecture count */
      narchit++;

      /* If architecture contains any Pfam-B domains, and is not the
         parent, then process */
      if (strstr(archit_ptr->architecture,"PB") != NULL &&
          strcmp(archit_ptr->architecture,parent_architecture))
        {
          /* Get this architecture's length */
          len = strlen(archit_ptr->architecture);

          /* Allocate temporary memory for temporary strings */
          tmp = (char *) malloc (sizeof(char)*(len + 1));
          newtmp = (char *) malloc (sizeof(char)*(2 * len + 1));
          newtmp[0] = '\0';

          /* Store the architecture */
          strcpy(tmp,archit_ptr->architecture);

          /* Free up the memory taken by old string */
          //OUT free(archit_ptr->architecture);

          /* Extract the dot-separated Pfam domains from the architecture */
          ntokens = tokenize(tmp,token,MAX_TOKENS,DOM_SEP[type][0]);

          /* Loop over the tokens to rebuild the architecture */
          for (itoken = 0; itoken < ntokens; itoken++)
            {
              /* If this is not the first token, then append a dot */
              if (itoken != 0)
                strcat(newtmp,DOM_SEP[type]);

              /* If a Pfam-A domain, then append */
              if (!strncmp(token[itoken],"PF",2))
                strcat(newtmp,token[itoken]);

              /* If domain is in the parent architecture, then leave
                 unchanged */
              else if (strstr(parent_architecture,token[itoken]) != NULL)
                strcat(newtmp,token[itoken]);

              /* Otherwise, add dummy Pfam-B domain */
              else
                strcat(newtmp,"PB000000");
            }

          /* Store the new architecture */
          store_string(&(archit_ptr->architecture),newtmp);

          /* Free up memory taken by the temporary strings */
          free(tmp);
          free(newtmp);
        }

      /* Get the next archit record */
      archit_ptr = archit_ptr->next_archit_ptr;
    }

  /* Return number of architectures */
  return(narchit);
}
/***********************************************************************

archit_sort  -  Shell-sort to arrange the archits by arcarchitecture

***********************************************************************/

void archit_sort(struct archit **archit_index_ptr,int narchits)
{
  long endloop, i, indi, indl, j, k, l, lognb2, m, n, nn;
  float aln2i = (float) 1.4426950, tiny = (float) 1.e-5;
  float calc;
  double dble, dble1;

  struct archit *archit1_ptr, *archit2_ptr, *swap_archit_ptr;

  /*--Initialise values */
  n = narchits;
  dble = n;
  dble1 = log(dble);
  calc = (float) dble1;
  calc = calc * aln2i + tiny;
  lognb2 = (long) calc;

  /* Perform sort */
  m = n;
  for (nn = 1; nn <= lognb2; nn++)
    {
      m = m / 2;
      k = n - m;
      for (j = 1; j <= k; j++)
        {
          i = j;
          endloop = FALSE;
          while (endloop == FALSE)
            {
              l = i + m;
              indi = i - 1;
              indl = l - 1;

              /* Get the two code records */
              archit1_ptr = archit_index_ptr[indi];
              archit2_ptr = archit_index_ptr[indl];

              /* If in wrong order, then swap sort-pointers */
              if (strcmp(archit1_ptr->architecture,archit2_ptr->architecture) > 0)
                {
                  swap_archit_ptr = archit_index_ptr[indi];
                  archit_index_ptr[indi] = archit_index_ptr[indl];
                  archit_index_ptr[indl] = swap_archit_ptr;

                  i = i - m;
                  if (i < 1)
                    endloop = TRUE;
                }
              else
                endloop = TRUE;
            }
        }
    }
}
/***********************************************************************

sort_archit_records  -  Sort the archit records by architecture

***********************************************************************/

void sort_archit_records(struct archit *first_archit_ptr,
                         struct archit ***archit_ind_ptr,int narchit,
                         struct parameters *params)
{
  int nentry;

  struct archit *archit_ptr;

  struct archit **archit_index_ptr;

  /* Create the pfam array as one large array */
  archit_index_ptr
    = (struct archit **) malloc(narchit * sizeof(struct archit *));

  /* Check that sufficient memory was allocated */
  if (archit_index_ptr == NULL)
    write_error_message(MEMORY_ALLOCATION,
                        "Cannot allocate memory for archits array",
                        params->output_type);

  /* Get pointer to the first code entry */
  archit_ptr = first_archit_ptr;
  nentry = 0;

  /* Loop through all entries to place each one in the array */
  while (archit_ptr != NULL)
    {
      /* Store current pointer at end of array */
      if (nentry < narchit)
        archit_index_ptr[nentry] = archit_ptr;

      /* Increment count of stored entries */
      nentry++;

      /* Get pointer to the next entry in the linked list */
      archit_ptr = archit_ptr->next_archit_ptr;
    }

  /* Sort the archit pointers into order of architecture */
  if (narchit > 1)
    archit_sort(archit_index_ptr,narchit);

  /* Return the array pointer */
  *archit_ind_ptr = archit_index_ptr;
}
/***********************************************************************

combine_architectures  -  Combine identical architectures into one

***********************************************************************/

int combine_architectures(struct archit **fst_archit_ptr,
                          struct archit **archit_index_ptr,int narchit)
{
  char last_archit_id[20];

  int ientry, nentry;

  struct archit *archit_ptr, *first_archit_ptr, *last_archit_ptr;
  struct hitlist *hitlist_ptr, *last_hitlist_ptr;

  /* Initialise variables */
  first_archit_ptr = *fst_archit_ptr;
  last_archit_ptr = NULL;
  last_hitlist_ptr = NULL;

  /* Get number of architecture entries */
  nentry = narchit;

  /* Re-initialise count of unique architectures */
  narchit = 0;

  /* Set the first pointer */
  first_archit_ptr = last_archit_ptr = archit_index_ptr[0];

  /* Initialise its pointer to the next record */
  first_archit_ptr->next_archit_ptr = NULL;

  /* Increment count of archit records */
  narchit++;

  /* Loop through the sorted list of architectures */
  for (ientry = 1; ientry < nentry; ientry++)
    {
      /* Get the current archit record */
      archit_ptr = archit_index_ptr[ientry];

      /* Get to point to NULL */
      archit_ptr->next_archit_ptr = NULL;

      /* If this is the same architecture as the last record, then
         combine the two */
      if (!strcmp(archit_ptr->architecture,last_archit_ptr->architecture))
        {
          /* If have a last hit for the previous architecture, start
             search for end at that one */
          if (last_hitlist_ptr != NULL)
            hitlist_ptr = last_hitlist_ptr;

          /* Otherwise, start search at its first hit */
          else
            hitlist_ptr = last_archit_ptr->first_hitlist_ptr;

          /* Re-initialise lst hit */
          last_hitlist_ptr = NULL;

          /* Loop through the previous record's hits until we reach
             the end */
          while (hitlist_ptr != NULL)
            {
              /* Save as the current last hit */
              last_hitlist_ptr = hitlist_ptr;

              /* Get the next hit record */
              hitlist_ptr = hitlist_ptr->next_hitlist_ptr;
            }

          /* Get the last hit to point to the first hit of the current
             architecture */
          if (last_hitlist_ptr != NULL)
            last_hitlist_ptr->next_hitlist_ptr
              = archit_ptr->first_hitlist_ptr;
          else
            last_archit_ptr->first_hitlist_ptr
              = archit_ptr->first_hitlist_ptr;

          /* Add number of hits */
          last_archit_ptr->nhits
            = last_archit_ptr->nhits + archit_ptr->nhits;

          /* If this is the parent architecture, make the last */

          /* Free up memory taken by current architecture record */
          //OUT free(archit_ptr->architecture);
          free(archit_ptr);
        }

      /* Otherwise, keep this architecture */
      else
        {
          /* Get last architecture to point to this one */
          last_archit_ptr->next_archit_ptr = archit_ptr;

          /* Null terminate list */
          archit_ptr->next_archit_ptr = NULL;

          /* Save this one as the last */
          last_archit_ptr = archit_ptr;

          /* Re-initialise the last hit-list pointer */
          last_hitlist_ptr = NULL;

          /* Increment count of archit records */
          narchit++;
        }
    }

  /* Free up the index array */
  free(archit_index_ptr);

  /* Return the first archit record */
  *fst_archit_ptr = first_archit_ptr;

  /* Return new count of archit records */
  return(narchit);
}
/***********************************************************************

test_hit  -  Test architecture of selected hit

***********************************************************************/

void test_hit(struct archit *first_archit_ptr,char *uniprot_acc)
{
  struct archit *archit_ptr;
  struct hit *hit_ptr;
  struct hitlist *hitlist_ptr, *next_hitlist_ptr, *prev_hitlist_ptr;

  /* Get the first architecture record */
  archit_ptr = first_archit_ptr;

  /* Loop over the architecture records to write out */
  while (archit_ptr != NULL)
    {
      /* Get the pointer to the first of the hits */
      hitlist_ptr = archit_ptr->first_hitlist_ptr;

      /* Loop through the hits to write out */
      while (hitlist_ptr != NULL)
        {
          /* Get the corresponding hit record */
          hit_ptr = hitlist_ptr->hit_ptr;

          /* If this is our record, then show architecture */
          if (!strcmp(uniprot_acc,hit_ptr->uniprot_acc))
            printf("Have hit record for %s -> Architecture %s\n",
                   hit_ptr->uniprot_acc,archit_ptr->architecture);

          /* Get the next hit record */
          hitlist_ptr = hitlist_ptr->next_hitlist_ptr;
        }

      /* Get the next architecture record */
      archit_ptr = archit_ptr->next_archit_ptr;
    }
}
/***********************************************************************

collapse_pfamb_domains  -  If have too many hits, then replace Pfam-B
                           domains by a single catch-all domain

***********************************************************************/

int collapse_pfamb_domains(int type,struct archit **fst_archit_ptr,
                           char *parent_architecture,
                           struct parameters *params)
{
  int narchit;

  struct archit *first_archit_ptr;

  struct archit **archit_index_ptr;

  /* Save the starting pointers */
  first_archit_ptr = *fst_archit_ptr;

  /* Replace all the Pfam-B domnains in the architecture records
     with PB000000 */
  narchit = replace_pfamb_domains(type,first_archit_ptr,parent_architecture);

  /* Sort the new architectures so that we can merge the new duplicates */
  sort_archit_records(first_archit_ptr,&archit_index_ptr,narchit,params);

  /* Combine identical architectures into one */
  narchit
    = combine_architectures(&first_archit_ptr,archit_index_ptr,narchit);

  /* Return new pointers */
  *fst_archit_ptr = first_archit_ptr;

  /* Return new number of architectures */
  return(narchit);
}
/***********************************************************************

strccnt  -  Count the number of occurrences of a given character in the
            given string

***********************************************************************/

int strccnt(const char *s,char c)
{
  int n = 0;

  while(*s != '\0') {
    if(*s++ == c)
      n++;
  }

  return n;
}
/***********************************************************************

modify_split_domains  -  For any split-CATH domains, replace dot by
                         colons

***********************************************************************/

void modify_split_domains(char *architecture)
{
  char ch;

  int i, len;
  int havep;

  /* Initialise variables */
  havep = FALSE;
  len = strlen(architecture);
  
  /* Loop through all the characters in this architecture */
  for (i = 0; i < len; i++)
    {
      /* Get this character */
      ch = architecture[i];

      /* If replacing, then perform replacament */
      if (ch == '.' && havep == TRUE)
        architecture[i] = ':';

      /* If have reached domain separator thenunset flag */
      else if (ch == '_')
        havep = FALSE;

      /* If have hist a 'p', set flag */
      else if (ch == 'p')
        havep = TRUE;
    }
}
/***********************************************************************

extract_domains  -  Extract the Pfam domains from the given architecture

***********************************************************************/

int extract_domains(int type,char *architecture,
                    char pfam_id[MAX_COMPARE][PFAM_LEN],int *domain_no,
                    int *ndom,int istart)
{
  char *tmp, *token[MAX_TOKENS];

  int copy, itoken, idomain, len, ndomains, ntokens, nstored, nunique;

  /* Initialise variables */
  ntokens = nunique = 0;
  nstored = istart;

  /* Get the length of the architecture */
  len = strlen(architecture);

  /* Allocate temporary memory for the architecture */
  tmp = (char *) malloc (sizeof(char)*(len + 1));

  /* Store the architecture */
  strcpy(tmp,architecture);

  /* Extract the dot-separated Pfam domains from the architecture */
  ndomains = ntokens = tokenize(tmp,token,MAX_TOKENS,DOM_SEP[type][0]);

  /* Loop over the Pfam ids to store the unique ones */
  for (itoken = 0; itoken < ntokens; itoken++)
    {
      /* Determine whether this a copy domain */
      copy = -1;
      for (idomain = 0; idomain < nstored && copy == -1; idomain++)
        {
          if (!strcmp(token[itoken],pfam_id[idomain]))
            copy = idomain;
        }

      /* If this is a copy domain, store the number identified */
      if (copy != -1)
        {
          /* Store as negative, if this is a PfamB domain */
          if (!strncmp(token[itoken],"PF",2) ||
              (token[itoken][0] != 'P' && token[itoken][0] != 'p'))
            domain_no[itoken] = copy;
          else
            domain_no[itoken] = -copy;
        }

      /* Otherwise, store this domain */
      else if (nunique < MAX_COMPARE)
        {
          /* Save the domain id */
          strcpy(pfam_id[nstored],token[itoken]);

          /* Save the domain number */
          if (!strncmp(token[itoken],"PF",2) ||
              (token[itoken][0] != 'P' && token[itoken][0] != 'p'))
            domain_no[itoken] = nstored;
          else
            domain_no[itoken] = -nstored;

          /* Increment counts of unique domains stored */
          nstored++;
          nunique++;
        }
    }

  /* Free up memory taken by the architecture */
  free(tmp);

  /* Return number of domains in the architecture */
  *ndom = ndomains;

  /* Return number of unique domains stored */
  return(nunique);
}
/***********************************************************************

score_extra_domains  -  Score the domains outside the matched region

***********************************************************************/

int score_extra_domains(int type,char *architecture,int istart,int iend,
                        int to_parent)
{
  char *string_ptr, *tmp;

  int add, dist, ipos, len;

  /* Initialise variables */
  dist = 0;
  ipos = 0;
  len = strlen(architecture);

#ifdef VERBOSE
  /* Create string for whole architecture */
  store_string(&tmp,architecture);
#endif

  /* Count the domains prior to the start of the matched region */
  while (ipos < istart)
    {
      /* Determine whether we have a PFam-A, Pfam-B CATH or split-CATH
         domain */
      if (!strncmp(architecture + ipos,"PF",2) ||
          (architecture[ipos] != 'P' && architecture[ipos] != 'p'))
        add = PFAMA_DIST;
      else if (architecture[ipos] == 'p')
        add = SPLIT_CATH_DIST;
      else
        add = PFAMB_DIST;

      /* Add to distance */
      dist = dist + add;

      /* Show this domain and its score */

#ifdef VERBOSE
      /* Extract just this domain */
      strcpy(tmp,architecture + ipos);
      string_ptr = strchr(tmp,DOM_SEP[type][0]);
      if (string_ptr != NULL)
        string_ptr[0] = '\0';

      /* Show this domain and its score */
      printf("   Extra domain %s, dist = %d",tmp,dist);
#endif

      /* Get start of next domain */
      string_ptr = strchr(architecture + ipos,DOM_SEP[type][0]);
      if (string_ptr != NULL)
        ipos = string_ptr - architecture + 1;
      else
        ipos = istart;
    }

  /* Start at the end of the range and continue to end of architecture */
  ipos = iend + 1;

  /* Count the domains after the matched region */
  while (ipos < len)
    {
      /* Determine whether we have a PFam-A, Pfam-B CATH or split-CATH
         domain */
      if (!strncmp(architecture + ipos,"PF",2) ||
          (architecture[ipos] != 'P' && architecture[ipos] != 'p'))
        add = PFAMA_DIST;
      else if (architecture[ipos] == 'p')
        add = SPLIT_CATH_DIST;
      else
        add = PFAMB_DIST;

      /* Add to distance */
      dist = dist + add;

#ifdef VERBOSE
      /* Extract just this domain */
      strcpy(tmp,architecture + ipos);
      string_ptr = strchr(tmp,DOM_SEP[type][0]);
      if (string_ptr != NULL)
        string_ptr[0] = '\0';

      /* Show this domain and its score */
      printf("   Extra domain %s, dist = %d",tmp,dist);
#endif

      /* Get start of next domain */
      string_ptr = strchr(architecture + ipos,DOM_SEP[type][0]);
      if (string_ptr != NULL)
        ipos = string_ptr - architecture + 1;
      else
        ipos = len;
    }

  /* If architecture is subset of parent (or vice versa), reduce distance */
  if (to_parent == TRUE)
    dist--;

  /* Return the scored distance */
  return(dist);
}
/***********************************************************************

create_alignment_arrays  -  Create the arrays for the sequence alignment
                            of the two sequences of title words and
                            author names

***********************************************************************/

void create_alignment_arrays(int nterm[2],int **arr,int **tmparr,
                             int **pth)
{
  int iterm, iseq, max_path, narray;

  int *array, *path, *tmp_array;

  /* Initialise variables */
  narray = nterm[0] * nterm[1];
  max_path = nterm[0] + nterm[1];
  
  /* Create Needleman-Wunsch array for performing the alignment */
  array = (int *) malloc(narray * sizeof(int));
  if (array == NULL)
    {
      printf("*** ERROR. Unable to allocate memory for NW array\n");
      exit(1);
    }

  /* Create an array for storing the path through the matrix that gives
     the best alignment score */
  path = (int *) malloc(max_path * sizeof(int));
  if (path == NULL)
    {
      printf("*** ERROR. Unable to allocate memory for path\n");
      exit(1);
    }

  /* Create working array for use in Needleman-Wunsch process */
  tmp_array = (int *) malloc(narray * sizeof(int));
  if (tmp_array == NULL)
    {
      printf("*** ERROR. Unable to allocate memory for temporary array\n");
      exit(1);
    }

  /* Return the array pointers */
  *arr = array;
  *pth = path;
  *tmparr = tmp_array;
}
/***********************************************************************

fill_array  -  Initialise array with domain identities

***********************************************************************/

void fill_array(int nterm[2],int domain_no[2][MAX_COMPARE],int *array,
                int *tmp_array)
{
  int i, ipoint, iterm, jterm;

  float similarity;

  /* Loop over the terms in first sequence */
  for (iterm = 0; iterm < nterm[FIRST]; iterm++)
    {
      /* Loop over the terms in the second sequence */
      for (jterm = 0; jterm < nterm[SECOND]; jterm++)
        {
          /* Get the array location of this term pair */
          ipoint = iterm * nterm[SECOND] + jterm;

          /* Update array point depending on whether the terms
             match */
          if (domain_no[FIRST][iterm] == domain_no[SECOND][jterm])
            array[ipoint] = 10;
          else
            array[ipoint] = 0;

          /* Store value in temporary array, too */
          tmp_array[ipoint] = array[ipoint];
        }
    }
}
/***********************************************************************

nw_sweep  -  Sweep through the alignment matrix, updating the cells
             according to the best path from any cell on the previous
             row or column - slow version

***********************************************************************/

void nw_sweep(int *array,int nterm[2])
{
  int ipoint, iterm, jpoint, jterm, kres, lres;
  int max_score, score;

  /* Loop backwards through the array, starting at the bottom right
     cell and working leftwards and upwards */
  for (iterm = nterm[0] - 2; iterm > -1; iterm--)
    { 
      for (jterm = nterm[1] - 2; jterm > -1; jterm--)
        {
          /* Get the array location of this residue pair */
          ipoint = iterm * nterm[1] + jterm;

          /* Initialise maximum score */
          max_score = -100;

          /* Loop over the cells to the right in the row below the
             current one */
          for (kres = iterm + 1, lres = jterm + 1;
               kres < nterm[0] && lres < nterm[1]; lres++)
            {
              /* Get this cell's array location */
              jpoint = kres * nterm[1] + lres;

              /* Calculate score to get from here to our current cell */
              score = array[jpoint];
              if (lres > jterm + 1)
                score = score - GAP_PENALTY;

              /* If this is the largest score so far, then store */
              if (score > max_score)
                max_score = score;
            }

          /* Loop over the downward cells in the column to the right
             of the current one */
          for (kres = iterm + 1, lres = jterm + 1;
               kres < nterm[0] && lres < nterm[1]; kres++)
            {
              /* Get this cell's array location */
              jpoint = kres * nterm[1] + lres;

              /* Calculate score to get from here to our current cell */
              score = array[jpoint];
              if (kres > iterm + 1)
                score = score - GAP_PENALTY;

              /* If this is the largest score so far, then store */
              if (score > max_score)
                max_score = score;
            }

          /* Add maximum score to cell's current score */
          array[ipoint] = array[ipoint] + max_score;
        }
    }
}
/***********************************************************************

find_best_path  -  Find the best path through the matrix, corresponding
                   to the best alignment of the two sequences

***********************************************************************/

int find_best_path(int *array,int nterm[2],int *path,int *alnlen)
{
  int ipoint, iterm, jterm, start_iterm, start_jterm;
  int align_len, last_pos[2], max_path, max_pos[2], path_len;
  int best_dist, dist, gap_penalty, max_score, score;
  int done, first;

  /* Initialise variables */
  align_len = 0;
  best_dist = 999;
  path_len = 0;
  start_iterm = 0;
  start_jterm = 0;
  done = FALSE;
  first = TRUE;
  last_pos[0] = -1;
  last_pos[1] = -1;
  max_path = nterm[0] + nterm[1];

  /* Loop until entire array has been traversed */
  while (done == FALSE)
    {
      /* Initialise maximum score and its position */
      gap_penalty = 0;
      max_score = -9999999;
      max_pos[0] = -1;
      max_pos[1] = -1;

      /* Loop through row to get maximum score and position */
      iterm = start_iterm;
      for (jterm = start_jterm; jterm < nterm[1]; jterm++)
        {
          /* Get the array location of this cell */
          ipoint = iterm * nterm[1] + jterm;
          dist = jterm - start_jterm;
          score = array[ipoint] - gap_penalty;

          /* If this is the highest score so far, store score and
             position */
          if (score > max_score)
            {
              max_score = score;
              max_pos[0] = iterm;
              max_pos[1] = jterm;
              best_dist = dist;
            }

          /* Set the gap penalty */
          if (first == FALSE)
            gap_penalty = GAP_PENALTY;
        }

      /* Re-initialise gap penalty */
      gap_penalty = 0;

      /* Loop down the column to get maximum score and position */
      jterm = start_jterm;
      for (iterm = start_iterm; iterm < nterm[0]; iterm++)
        {
          /* Get the array location of this cell */
          ipoint = iterm * nterm[1] + jterm;
          dist = iterm - start_iterm;
          score = array[ipoint] - gap_penalty;

          /* If this is the highest score so far, store score and
             position */
          if (score > max_score ||
              (score == max_score && dist < best_dist))
            {
              max_score = score;
              max_pos[0] = iterm;
              max_pos[1] = jterm;
              best_dist = dist;
            }

          /* Set the gap penalty */
          if (first == FALSE)
            gap_penalty = GAP_PENALTY;
        }

      /* Store best position in the path array */
      if (path_len < max_path - 1)
        {
          /* Get array position */
          ipoint = max_pos[0] * nterm[1] + max_pos[1];

          /* Check for gap in first sequence */
          if (max_pos[0] - last_pos[0] > 1)
            {
              /* Add insertion in first sequence to overall alignment
                 length */
              align_len = align_len + max_pos[0] - last_pos[0] - 1;
            }

          /* Check for gap in second sequence */
          else if (max_pos[1] - last_pos[1] > 1)
            {
              /* Add insertion in second sequence to overall alignment
                 length */
              align_len = align_len + max_pos[1] - last_pos[1] - 1;
            }

          /* Increment total alignment length */
          align_len++;

          /* Store current position */
          last_pos[0] = max_pos[0];
          last_pos[1] = max_pos[1];

          /* Store in path */
          path[path_len] = ipoint;
          path_len++;
        }

      /* If have exceeded maximum allowed path length, then stop
         here */
      else
        {
          printf("*** ERROR. Alignment path length exceeded!!\n");
          printf("***        Serious program error ...\n");
          exit(-1);
        }

      /* Set the new start-position as diagonally right and below
         the point we have just reached */
      start_iterm = max_pos[0] + 1;
      start_jterm = max_pos[1] + 1;

      /* Check if we have reached the end of the array */
      if (start_iterm > nterm[0] - 1 || start_jterm > nterm[1] - 1)
        {
          /* If haven't finished on the very last cell of the array,
             make path end one cell beyond it */
          if (start_iterm != nterm[0] || start_jterm != nterm[1])
            {
              /* Check for gap in first sequence */
              if (nterm[0] - last_pos[0] > 1)
                {
                  /* Add insertion in first sequence to overall alignment
                     length */
                  align_len = align_len + nterm[0] - last_pos[0] - 1;
                }

              /* Check for gap in second sequence */
              else if (nterm[1] - last_pos[1] > 1)
                {
                  /* Add insertion in second sequence to overall alignment
                     length */
                  align_len = align_len + nterm[1] - last_pos[1] - 1;
                }

              /* Store in path */
              path[path_len] = -1;
              path_len++;
            }

          /* Set flag that we're done */
          done = TRUE;
        }

      /* Unset flag */
      first = FALSE;
    }

  /* Return total length of the alignment */
  *alnlen = align_len;

  /* Return the alignment path length */
  return(path_len);
}
/***********************************************************************

calc_score  -  Calculate the distance between the two aligned
               architectures

***********************************************************************/

int calc_score(int *array,int *path,int path_len,int align_len,
               int nterm[2],int domain_no[2][MAX_COMPARE])
{
#ifdef VERBOSE
  char seq1[LINELEN], seq2[LINELEN], number_string[20];
#endif

  int align_pos, array_pos, ipoint, ipos;
  int i, iterm, jterm, kterm, last_iterm, last_jterm;
  int nidentical, noverlap, nlong, nmax, nsame;
  int score, similarity;
  int in_overlap;

  /* Initialise variables */
  align_pos = 0;
  last_iterm = -1;
  last_jterm = -1;
  in_overlap = FALSE;
  noverlap = 0;
  nidentical = 0;
  nlong = nterm[0];
  if (nterm[1] > nlong)
    nlong = nterm[1];
  if (nlong == 0)
    nlong = 1;
  score = 0;
#ifdef VERBOSE
  seq1[0] = seq2[0] = '\0';
#endif

  /* Loop over the points along the path */
  for (ipos = 0; ipos < path_len && align_pos < align_len; ipos++)
    {
      /* Get this point in the path */
      ipoint = path[ipos];

      /* Extract row and column number of this cell */
      if (ipoint == -1)
        {
          iterm = nterm[0];
          jterm = nterm[1];
          in_overlap = FALSE;
        }
      else
        {
          iterm = ipoint / nterm[1];
          jterm = ipoint - iterm * nterm[1];
        }

      /* Check for gap in second sequence */
      if (iterm - last_iterm > 1)
        {
          /* Store insertion from first sequence */
          for (kterm = last_iterm + 1; kterm < iterm; kterm++)
            {
              /* Increment alignment position */
              align_pos++;

              /* Increment overlap length */
              if (in_overlap == TRUE)
                noverlap++;

              /* Add this mismatched domain to overall distance */
              if (domain_no[PARENT][kterm] < 0)
                score = score + PFAMB_DIST;
              else
                score = score + PFAMA_DIST;
              if (in_overlap == TRUE)
                score = score + 1;
#ifdef VERBOSE
              sprintf(number_string,":%3d",domain_no[PARENT][kterm]);
              strcat(seq1,number_string);
              strcat(seq2,"    ");
#endif
            }
        }

      /* Check for gap in first sequence */
      else if (jterm - last_jterm > 1)
        {
          /* Store insertion from second sequence */
          for (kterm = last_jterm + 1; kterm < jterm; kterm++)
            {
              /* Increment alignment position */
              align_pos++;

              /* Increment overlap length */
              if (in_overlap == TRUE)
                noverlap++;

              /* Add this mismatched domain to overall distance */
              if (domain_no[CURRENT][kterm] < 0)
                score = score + PFAMB_DIST;
              else
                score = score + PFAMA_DIST;
              if (in_overlap == TRUE)
                score = score + 1;
#ifdef VERBOSE
              sprintf(number_string,":%3d",domain_no[CURRENT][kterm]);
              strcat(seq2,number_string);
              strcat(seq1,"    ");
#endif
            }
        }

      /* If valid, store the current alignment position */
      if (iterm < nterm[0] && jterm < nterm[1])
        {
          /* Increment alignment position */
          align_pos++;

          /* Calculate the similarity score of the two strings at this
             alignment position */
          if (domain_no[PARENT][iterm] == domain_no[CURRENT][jterm])
            nidentical++;

          /* If they don't match, then increment distance score */
          else if (domain_no[PARENT][iterm] < 0 &&
                   domain_no[CURRENT][jterm] < 0)
            score = score + 2 * PFAMB_DIST;
          else
            score = score + 2 * PFAMA_DIST;

          /* Increment overlap length */
          in_overlap = TRUE;
          noverlap++;
#ifdef VERBOSE
          sprintf(number_string,":%3d",domain_no[PARENT][iterm]);
          strcat(seq1,number_string);
          sprintf(number_string,":%3d",domain_no[CURRENT][jterm]);
          strcat(seq2,number_string);
#endif
        }

      /* Store current position */
      last_iterm = iterm;
      last_jterm = jterm;
    }

#ifdef VERBOSE
   printf("Alignment:\n");
   printf("%s\n",seq1);
   printf("%s\n",seq2);
#endif

  /* Return the calculated match score */
  return(score);
}
/***********************************************************************

free_alignment_arrays  -  Free up the memory taken up by the alignment
                          arrays

***********************************************************************/

void free_alignment_arrays(int *array,int *path,int *tmp_array)
{
  /* Clear the NW arrays */
  free(array);
  free(path);
  free(tmp_array);
}
/***********************************************************************

score_architectures  -  Score the given pair of architectures

***********************************************************************/

int score_architectures(int type,char *architecture1,char *architecture2,
                        int *ndomains,int *nunique,
                        char pfam_id[MAX_COMPARE][PFAM_LEN],
                        int domain_no[2][MAX_COMPARE],
                        char *parent_architecture)
{
  char *string_ptr;

  int align_len, dist, iend, istart, path_len;
  int to_parent;

  int *array, *path, *tmp_array;

  /* Initialise variables */
  dist = 0;
  to_parent = FALSE;

#ifdef VERBOSE
   printf("Comparing = %s ",architecture1);
   if (!strcmp(architecture1,parent_architecture))
     printf("[P]");
   printf(" vs %s ",architecture2);
   if (!strcmp(architecture2,parent_architecture))
     printf("[P]");
   printf("\n");
#endif

  /* If architecture1 matches architecture2 exactly, then set the
     distance as zero */
  if (!strcmp(architecture1,architecture2))
    dist = 0;

  /* If architecture2 contains architecture1, then calculate
     distance according to number of additional domains
     architecture2 contains */
  else if ((string_ptr = strstr(architecture2,architecture1))
           != NULL )
    {
      /* If architecture1 is the parent architecture, set flag
         to halve score of additional domains */
      if (!strcmp(architecture1,parent_architecture))
        {
          /* Set flag to reduce extra-domains score */
          to_parent = TRUE;
          // dist = PARENT_SUBSTRING_DIST;
        }

      /* Get the start-position of the match */
      istart = string_ptr - architecture2;
      iend = istart + strlen(architecture1);

      /* Score the extra domains outside the matched region */
      dist = score_extra_domains(type,architecture2,istart,iend,
                                 to_parent);
    }

  /* Try the match the other way */
  else if ((string_ptr = strstr(architecture1,architecture2))
           != NULL )
    {
      /* If architecture2 is the parent architecture, set flag
         to halve score of additional domains */
      if (!strcmp(architecture2,parent_architecture))
        {
          /* Set flag to reduce extra-domains score */
          to_parent = TRUE;
        }

      /* Get the start-position of the match */
      istart = string_ptr - architecture1;
      iend = istart + strlen(architecture2);
      
      /* Score the extra domains outside the matched region */
      dist = score_extra_domains(type,architecture1,istart,iend,
                                 to_parent);
    }

  /* Otherwise, perform N & W alignment to determine distance */
  else
    {
      /* Determine the number of domains in architecture2 */
      //ndomains = strccnt(architecture2,'.') + 1;

      /* Create the arrays to be used in computing the
         'sequence alignment' between the two architectures */
      create_alignment_arrays(ndomains,&array,&tmp_array,&path);

      /* Initialise array with identities */
      fill_array(ndomains,domain_no,array,tmp_array);

      /* Apply Needleman and Wunsch algorithm to adjust scores
         matrix */
      nw_sweep(array,ndomains);

      /* Find best route through matrix to give best alignment
         of the two sequences */
      path_len = find_best_path(array,ndomains,path,&align_len);

      /* Calculate the match score between the aligned
         architectures */
      dist = calc_score(array,path,path_len,align_len,ndomains,
                        domain_no);
          
      /* Free up memory taken by the alignment arrays */
      free_alignment_arrays(array,tmp_array,path);
    }
#ifdef VERBOSE
   printf("   Distance score = %d\n",dist);
#endif

  /* Return the distance between the to architectures */
  return(dist);
}
/***********************************************************************

compare_to_parent  -  Score all architecture distances from the parent
                      architecture to calculate a cut-off score for
                      which architectures to keep

***********************************************************************/

int compare_to_parent(int type,struct archit *first_archit_ptr,
                      char *parent_architecture,struct parameters *params)
{
  char pfam_id[MAX_COMPARE][PFAM_LEN];
  char p_architecture[LINELEN + 1], n_architecture[LINELEN + 1];

  char *string_ptr;

  int count[MAX_SCORE + 1], dist, i, istart, iend, len, narchit;
  int domain_no[2][MAX_COMPARE], ndomains[2], nunique[2];
  int cut_off;

  struct archit *archit_ptr;

  /* Initialise variables */
  cut_off = -1;
  for (i = 0; i < MAX_SCORE + 1; i++)
    count[i] = 0;

  /* Save the parent architecture */
  strcpy(p_architecture,parent_architecture);

  /* Modify any split-CATH domains */
  if (params->cath_domains == TRUE)
    modify_split_domains(p_architecture);

  /* Extract the Pfam domains */
  nunique[PARENT]
    = extract_domains(type,p_architecture,pfam_id,domain_no[PARENT],
                      &(ndomains[PARENT]),0);

  /* Get the first hit record */
  archit_ptr = first_archit_ptr;

  /* Loop over the hit records to write out */
  while (archit_ptr != NULL)
    {
      /* Save this architecture */
      strcpy(n_architecture,archit_ptr->architecture);

      /* Modify any split-CATH domains */
      if (params->cath_domains == TRUE)
        modify_split_domains(n_architecture);

      /* Extract the domains from this architecture */
      nunique[CURRENT]
        = extract_domains(type,n_architecture,pfam_id,
                          domain_no[CURRENT],&(ndomains[CURRENT]),
                          nunique[PARENT]);

      /* Calculate the distance score between this architecture and the
         parent architecture */
      dist = score_architectures(type,p_architecture,n_architecture,
                                 ndomains,nunique,pfam_id,domain_no,
                                 p_architecture);

      /* Store the score count */
      if (dist > MAX_SCORE)
        dist = MAX_SCORE;
      count[dist] = count[dist] + 1;

      /* Attach the score to this architecture */
      archit_ptr->score = dist;

      /* Get the next hit record */
      archit_ptr = archit_ptr->next_archit_ptr;
    }

  /* Determine the cut-off score which will give the set of
     architectures that are the most similar to the parent architecture */
  narchit = 0;
  for (i = 0; i < MAX_SCORE + 1 && cut_off == -1; i++)
    {
      /* Add the count for the current score to the total */
      narchit = narchit + count[i];

      /* If count has exceeded the maximum, then set the cut-off at this
         score */
      if (narchit > params->max_architectures)
        cut_off = i;
    }

  /* Return the cut-off value */
  return(cut_off);
}
/***********************************************************************

delete_architectures  -  Delete the duplicate Pfam records

***********************************************************************/

int delete_architectures(struct archit **fst_archit_ptr,int cut_off)
{
  int ndeleted, narchit;

  struct archit *archit_ptr, *next_archit_ptr, *prev_archit_ptr;
  struct archit *first_archit_ptr;

  /* Initialise variables */
  first_archit_ptr = *fst_archit_ptr;
  ndeleted = 0;
  narchit = 0;

  /* Initialise pointers */
  prev_archit_ptr = NULL;

  /* Get the first archit */
  archit_ptr = first_archit_ptr;

  /* Loop through the archit records to delete those above the cut-off */
  while (archit_ptr != NULL)
    {
      /* Get the next archit */
      next_archit_ptr = archit_ptr->next_archit_ptr;

      /* If score is above the cut-off, then delete */
      if (archit_ptr->score > cut_off)
        {
          /* If have a previous archit record, then get it to point
             around this one */
          if (prev_archit_ptr != NULL)
            prev_archit_ptr->next_archit_ptr = next_archit_ptr;

          /* Otherwise, this must be the first archit, so make
             the next one the new first */
          else
            first_archit_ptr = next_archit_ptr;

          /* Free up the memory taken by the deleted record */
          //OUT free(archit_ptr->architecture);
          free(archit_ptr);

          /* Increment count of deleted archits */
          ndeleted++;
        }

      /* Otherwise, save as the new previous record */
      else
        {
          /* Save current archit */
          prev_archit_ptr = archit_ptr;

          /* Increment count of retained archit records */
          narchit++;
        }

      /* Get the next archit record */
      archit_ptr = next_archit_ptr;
    }

  /* Return new starting archit */
  *fst_archit_ptr = first_archit_ptr;

  /* Return number of retained archits */
  return(narchit);
}
/***********************************************************************

extract_pfam_domains  -  Extract the Pfam domains from the remaining
                         architectures

***********************************************************************/

int extract_pfam_domains(int type,struct archit *first_archit_ptr,
                         struct pfam **fst_pfam_ptr,
                         struct parameters *params)
{
  char *string_ptr;

  int len, npfam;

  struct archit *archit_ptr;
  struct pfam *pfam_ptr, *first_pfam_ptr, *last_pfam_ptr;

  /* Initialise variables */
  npfam = 0;

  /* Initialise pointers */
  first_pfam_ptr = last_pfam_ptr = NULL;

  /* Get the first hit record */
  archit_ptr = first_archit_ptr;

  /* Loop over the hit records to write out */
  while (archit_ptr != NULL)
    {
      /* Split up the architecture into its constituent Pfam
         domains */
      npfam
        = split_architecture(type,&first_pfam_ptr,&last_pfam_ptr,
                             archit_ptr->architecture,npfam,params);

      /* Get the next hit record */
      archit_ptr = archit_ptr->next_archit_ptr;
    }

  /* Return pointer to first Pfam record */
  *fst_pfam_ptr = first_pfam_ptr;

  /* Return the number of Pfam domains created */
  return(npfam);
}
/***********************************************************************

prune_sequences  -  Prune sequences for any architectures that have a
                    huge number

***********************************************************************/

int prune_sequences(struct archit *first_archit_ptr,
                    struct hit *parent_hit_ptr,struct parameters *params)
{
  char *string_ptr;

  int loop, nhits, nseqs;
  int done, have_parent;

  struct archit *archit_ptr;
  struct hit *hit_ptr;
  struct hitlist *hitlist_ptr, *next_hitlist_ptr, *prev_hitlist_ptr;

  /* Initialise variables */
  nseqs = 0;

  /* Get the first hit record */
  archit_ptr = first_archit_ptr;

  /* Loop over the hit records to write out */
  while (archit_ptr != NULL)
    {
      /* If this architecture's number of hits exceeds the limit, then
         remove the excess */
      if (archit_ptr->nhits > params->max_sequences &&
          params->max_sequences != -1)
        {
          /* Initialise flags */
          done = FALSE;
          have_parent = FALSE;

          /* Loop twice: first to shoft sequences with structural info
             to the front of the linked list, and then to shift the
             parent sequence, if present */
          for (loop = 0; loop < 2 && done == FALSE; loop++)
            {
              /* Get the pointer to the first hit list record */
              hitlist_ptr = archit_ptr->first_hitlist_ptr;

              /* Initialize previous hit */
              prev_hitlist_ptr = NULL;

              /* Loop through the hits to shift any with 3D structure to
                 the top of the list */
              while (hitlist_ptr != NULL)
                {
                  /* Get the corresponding hit record */
                  hit_ptr = hitlist_ptr->hit_ptr;

                  /* If this is the parent hit, set flag */
                  if (hit_ptr == parent_hit_ptr)
                    have_parent = TRUE;

                  /* Save the next hitlist record */
                  next_hitlist_ptr = hitlist_ptr->next_hitlist_ptr;

                  /* If this hit has some structural information, or is
                     the parent sequence, send it to the top of the list */
                  if (((loop == 0 && strcmp(hit_ptr->coverage,"NONE")) ||
                       (loop == 1 && hit_ptr == parent_hit_ptr)) &&
                      hitlist_ptr != archit_ptr->first_hitlist_ptr)
                    {
                      /* Get the previous record to point round this one */
                      prev_hitlist_ptr->next_hitlist_ptr
                        = hitlist_ptr->next_hitlist_ptr;

                      /* Get this record to point to the first record */
                      hitlist_ptr->next_hitlist_ptr
                        = archit_ptr->first_hitlist_ptr;

                      /* Make it the first record */
                      archit_ptr->first_hitlist_ptr = hitlist_ptr;
                    }

                  /* Otherwise, set this as the previous record */
                  else
                    prev_hitlist_ptr = hitlist_ptr;

                  /* Get the next hit record */
                  hitlist_ptr = next_hitlist_ptr;
                }

              /* If not parent hit found, then don't need to loop a 
                 second time */
              if (have_parent == FALSE)
                done = TRUE;
            }

          /* Initialise count of hits */
          nhits = 0;

          /* Get the pointer to the first hit list record again */
          hitlist_ptr = archit_ptr->first_hitlist_ptr;

          /* Loop through a second time to truncate at the required
             number of hits */
          while (hitlist_ptr != NULL)
            {
              /* Increment hit counter */
              nhits++;

              /* If this is the maximum number of sequences, then
                 truncate at this hit */
              if (nhits >= params->max_sequences)
                hitlist_ptr->next_hitlist_ptr = NULL;

              /* Get the next hit record */
              hitlist_ptr = hitlist_ptr->next_hitlist_ptr;
            }

          /* Add to total number of sequences */
          nseqs = nseqs + params->max_sequences;
        }

      /* Accumulate total number of sequences */
      else
        nseqs = nseqs + archit_ptr->nhits;

      /* Get the next hit record */
      archit_ptr = archit_ptr->next_archit_ptr;
    }

  /* Return total number of sequences retained */
  return(nseqs);
}
/***********************************************************************

write_hit_summary  -  Write out the number of hits returned by the
                      search, plus details of the parent sequence

***********************************************************************/

void write_hit_summary(int nhits,int ntotal,int ncollapsed,int nfiltered,
                       struct hit *parent_hit_ptr,int npfam,int nspecies,
                       struct parameters *params)
{
  fprintf(outfile_ptr,":SUMMARY\n");

  /* Write the number of hits */
  fprintf(outfile_ptr,"%d\n",nhits);
  fprintf(outfile_ptr,"%d\n",ntotal);
  fprintf(outfile_ptr,"%d\n",ncollapsed);
  fprintf(outfile_ptr,"%d\n",nfiltered);
  fprintf(outfile_ptr,"%d\n",params->npfam);
  fprintf(outfile_ptr,"%d\n",npfam);
  fprintf(outfile_ptr,"%d\n",nspecies);

  /* Write out the details of the parent sequence */
  fprintf(outfile_ptr,":PARENT\n");
  fprintf(outfile_ptr,"%s\t%s\t%s\t%s\t%d\n",
          parent_hit_ptr->uniprot_acc,
          parent_hit_ptr->uniprot_id,
          parent_hit_ptr->architecture,
          parent_hit_ptr->protein_name,
          parent_hit_ptr->seq_len);
}
/***********************************************************************

write_dat_summary  -  Write out hit stats to the .dat file

***********************************************************************/

void write_dat_summary(int nhits,int ntotal,int ncollapsed,int nfiltered,
                       struct hit *parent_hit_ptr,int npfam,int nspecies,
                       struct pfam *first_pfam_ptr,
                       struct parameters *params)
{
  int itoken, ntokens;

  char *token[MAX_TOKENS];

  struct pfam *pfam_ptr;

  /* Write the number of hits */
  fprintf(outfile_ptr,"NHITS %d\n",nhits);
  fprintf(outfile_ptr,"NTOTAL %d\n",ntotal);
  fprintf(outfile_ptr,"NCOLLAPSED %d\n",ncollapsed);
  fprintf(outfile_ptr,"NFILTERED %d\n",nfiltered);
  fprintf(outfile_ptr,"PARAM_NPFAM %d\n",params->npfam);
  fprintf(outfile_ptr,"NPFAM %d\n",npfam);
  fprintf(outfile_ptr,"NSPECIES %d\n",nspecies);

  /* Write out the details of the domain being searched */
  fprintf(outfile_ptr,"SEARCH_ARCHITECTURE %s\n",
          params->architecture);

  /* Break up architecture into individual Pfam domains */
  ntokens
    = tokenize(params->architecture,token,MAX_TOKENS,'.');
  fprintf(outfile_ptr,"NSEARCH_ID %d\n",ntokens);

  /* Loop over the domains to write each one out */
  for (itoken = 0; itoken < ntokens; itoken++)
    {
      /* Write the domain id */
      fprintf(outfile_ptr,"SEARCH_ID[%d] %s\n",itoken + 1,token[itoken]);

      /* Find the corresponding Pfam record to get domain name */
      pfam_ptr = find_pfam(first_pfam_ptr,token[itoken]);
      if (pfam_ptr != NULL)
        {
          fprintf(outfile_ptr,"SEARCH_NAME[%d] %s\n",itoken + 1,
                  pfam_ptr->pfam_name);
          fprintf(outfile_ptr,"SEARCH_SHORT[%d] %s\n",itoken + 1,
                  pfam_ptr->short_name);

          /* Write out Pfam-B flag */
          if (!strncmp(pfam_ptr->pfam_id,"PF",2))
            {
              fprintf(outfile_ptr,"SEARCH_PFAMB[%d] FALSE\n",itoken + 1);
              fprintf(outfile_ptr,"SEARCH_PFAM_CAT[%d] A\n",itoken + 1);
            }
          else
            {
              fprintf(outfile_ptr,"SEARCH_PFAMB[%d] TRUE\n",itoken + 1);
              fprintf(outfile_ptr,"SEARCH_PFAM_CAT[%d] B\n",itoken + 1);
            }

          fprintf(outfile_ptr,"SEARCH_PFAM_COLOUR[%d] %d\n",itoken + 1,
                  pfam_ptr->colour);
        }
    }

  /* Write out the remaining details */
  fprintf(outfile_ptr,"PARENT_UNIPROT_ACC %s\n",
          parent_hit_ptr->uniprot_acc);
  fprintf(outfile_ptr,"PARENT_UNIPROT_ID %s\n",
          parent_hit_ptr->uniprot_id);
  fprintf(outfile_ptr,"PARENT_ARCHITECTURE %s\n",
          parent_hit_ptr->architecture);
  fprintf(outfile_ptr,"PARENT_PROTEIN_NAME %s\n",
          parent_hit_ptr->protein_name);
  fprintf(outfile_ptr,"PARENT_SEQ_LEN %d\n",parent_hit_ptr->seq_len);
}
/***********************************************************************

score_sort  -  Shell-sort to sort architectures by similarity to
               parent architecture

***********************************************************************/

void score_sort(struct archit **archit_index_ptr,int narchit)
{
  int match;
  long endloop, i, indi, indl, j, k, l, lognb2, m, n, nn;
  float aln2i = (float) 1.4426950, tiny = (float) 1.e-5;
  float calc;
  double dble, dble1;

  struct archit *archit1_ptr, *archit2_ptr, *swap_archit_ptr;

  /*--Initialise values */
  n = narchit;
  dble = n;
  dble1 = log(dble);
  calc = (float) dble1;
  calc = calc * aln2i + tiny;
  lognb2 = (long) calc;

  /* Perform sort */
  m = n;
  for (nn = 1; nn <= lognb2; nn++)
    {
      m = m / 2;
      k = n - m;
      for (j = 1; j <= k; j++)
        {
          i = j;
          endloop = FALSE;
          while (endloop == FALSE)
            {
              l = i + m;
              indi = i - 1;
              indl = l - 1;

              /* Get the two code records */
              archit1_ptr = archit_index_ptr[indi];
              archit2_ptr = archit_index_ptr[indl];

              /* If in wrong order, then swap sort-pointers */
              if (archit1_ptr->score > archit2_ptr->score ||
                  (archit1_ptr->score == archit2_ptr->score &&
                   strcmp(archit1_ptr->architecture,
                          archit2_ptr->architecture) > 0))
                {
                  swap_archit_ptr = archit_index_ptr[indi];
                  archit_index_ptr[indi] = archit_index_ptr[indl];
                  archit_index_ptr[indl] = swap_archit_ptr;

                  i = i - m;
                  if (i < 1)
                    endloop = TRUE;
                }
              else
                endloop = TRUE;
            }
        }
    }
}
/***********************************************************************

sort_by_score  -  Sort the architectures by score

***********************************************************************/

void sort_by_score(struct archit **fst_archit_ptr,int narchit,
                   struct parameters *params)
{
  int iorder, nentry;

  struct archit *archit_ptr, *first_archit_ptr, *last_archit_ptr;

  struct archit **archit_index_ptr;

  /* Initialise variables */
  first_archit_ptr = *fst_archit_ptr;
  last_archit_ptr = NULL;

  /* Create the archit array as one large array */
  archit_index_ptr
    = (struct archit **) malloc(narchit * sizeof(struct archit *));

  /* Check that sufficient memory was allocated */
  if (archit_index_ptr == NULL)
    write_error_message(MEMORY_ALLOCATION,
                        "Cannot allocate memory for archit array",
                        params->output_type);

  /* Get pointer to the first code entry */
  archit_ptr = first_archit_ptr;
  nentry = 0;

  /* Loop through all entries to place each one in the archit array */
  while (archit_ptr != NULL)
    {
      /* Store current pointer at end of the archit array */
      if (nentry < narchit)
        archit_index_ptr[nentry] = archit_ptr;

      /* Increment count of stored entries */
      nentry++;

      /* Get pointer to the next entry in the linked list */
      archit_ptr = archit_ptr->next_archit_ptr;
    }

  /* Sort the code pointers into order of primary archit */
  score_sort(archit_index_ptr,narchit);

  /* Having sorted the records, rearrange order of the pointers */
  for (iorder = 0; iorder < narchit; iorder++)
    {
      /* Get the current archit record */
      archit_ptr = archit_index_ptr[iorder];

      /* If this is the first, then make the first record */
      if (iorder == 0)
        first_archit_ptr = archit_ptr;

      /* Otherwise, get the previous archit record to point to this one */
      else
        last_archit_ptr->next_archit_ptr = archit_ptr;

      /* Blank out current record's next archit pointer */
      archit_ptr->next_archit_ptr = NULL;

      /* Store current record as last encountered */
      last_archit_ptr = archit_ptr;
    }

  /* Free up the index array created here */
  free(archit_index_ptr);

  /* Return the first Archit record */
  *fst_archit_ptr = first_archit_ptr;
}
/***********************************************************************

find_enzyme  -  Find the given enzyme entry

***********************************************************************/

struct enzyme *find_enzyme(struct enzyme *first_enzyme_ptr,
                           char *ec_number)
{
  struct enzyme *enzyme_ptr, *found_enzyme_ptr;

  /* Initialise variables */
  found_enzyme_ptr = NULL;

  /* Get pointer to the first enzyme */
  enzyme_ptr = first_enzyme_ptr;

  /* Loop over all the enzymes to find the one of interest */
  while (enzyme_ptr != NULL && found_enzyme_ptr == NULL)
    {
      /* If this is the right enzyme, store pointer */
      if (!strcmp(ec_number,enzyme_ptr->ec_number))
        found_enzyme_ptr = enzyme_ptr;

      /* Get pointer to the next enzyme */
      enzyme_ptr = enzyme_ptr->next_enzyme_ptr;
    }

  /* Return the enzymeing enzyme pointer */
  return(found_enzyme_ptr);
}
/***********************************************************************

convert_ec  -  Convert the E.C. number into form suitable for sorting

***********************************************************************/

void convert_ec(char *ec_num,char *sort_number)
{
  char ec_number[EC_LEN], ch, number_string[20];

  char *string_ptr, *token[MAX_TOKENS];

  int inum, ipos, itoken, ntokens;

  /* initialise variables */
  sort_number[0] = '\0';

  /* Check for comma, indicating a list */
  string_ptr = strchr(ec_num,',');

  /* If have a list of E.C. numbers, consider only the first */
  if (string_ptr != NULL)
    {
      /* Get the position of the comma */
      ipos = string_ptr - ec_num;
      
      /* Strip off all but the first E.C. number */
      if (ipos < EC_LEN / 2)
        {
          strncpy(ec_number,ec_num,ipos);
          ec_number[ipos] = '\0';
        }
      else
        strcpy(ec_number,"XXXX");
    }
  else
    strcpy(ec_number,ec_num);

  /* Break up the number into dot-separated tokens */
  ntokens = tokenize(ec_number,token,MAX_TOKENS,'.');
  for (itoken = 0; itoken < ntokens; itoken++)
    {
      /* If not the first number, append a dot */
      if (itoken > 0)
        strcat(sort_number,".");

      /* Format the sort number */
      ch = token[itoken][0];
      if (ch > '0' && ch <= '9')
        {
          inum = atoi(token[itoken]);
          sprintf(number_string,"%4.4d",inum);
          strcat(sort_number,number_string);
        }
      else if (ch == '-')
        strcat(sort_number,"0000");
      else
        {
          inum = atoi(token[itoken]+1);
          sprintf(number_string,"%c%3.3d",ch,inum);
          strcat(sort_number,number_string);
        }
    }
}
/***********************************************************************

create_enzyme_record  -  Create and initialise a new enzyme record

***********************************************************************/

struct enzyme *create_enzyme_record(struct enzyme **fst_enzyme_ptr,
                                    struct enzyme **lst_enzyme_ptr,
                                    char *ec_number,
                                    struct parameters *params)
{
  char sort_number[EC_LEN + 1];

  struct enzyme *enzyme_ptr, *first_enzyme_ptr, *last_enzyme_ptr;

  /* Initialise enzyme pointers */
  enzyme_ptr = NULL;
  first_enzyme_ptr = *fst_enzyme_ptr;
  last_enzyme_ptr = *lst_enzyme_ptr;

  /* Allocate memory for structure to hold enzyme info */
  enzyme_ptr = (struct enzyme *) malloc(sizeof(struct enzyme));
  if (enzyme_ptr == NULL)
    write_error_message(MEMORY_ALLOCATION,
                        "Cannot allocate memory for struct enzyme",
                        params->output_type);

  /* Store known elements of the structure */
  store_string(&(enzyme_ptr->ec_number),ec_number);

  /* Convert E.C. number into format suitable for sorting */
  convert_ec(ec_number,sort_number);

  /* Save the sort-number */
  store_string(&(enzyme_ptr->sort_number),sort_number);

  /* Initialise remaining fields */
  enzyme_ptr->enzyme_name = &null;
  enzyme_ptr->narch = 0;
  enzyme_ptr->nseqs = 0;
  enzyme_ptr->ntotal = 0;

  /* If this is the very first, save its pointer */
  if (first_enzyme_ptr == NULL)
    first_enzyme_ptr = enzyme_ptr;

  /* Add link from previous enzyme to the current one */
  if (last_enzyme_ptr != NULL)
    last_enzyme_ptr->next_enzyme_ptr = enzyme_ptr;
  last_enzyme_ptr = enzyme_ptr;

  /* Initialise pointer to next enzyme record */
  enzyme_ptr->next_enzyme_ptr = NULL;

  /* Return current pointers */
  *fst_enzyme_ptr = first_enzyme_ptr;
  *lst_enzyme_ptr = last_enzyme_ptr;

  /* Return new enzyme pointer */
  return(enzyme_ptr);
}
/***********************************************************************

write_hits  -  Write out all the hit records

***********************************************************************/

int write_hits(struct archit *first_archit_ptr,char *parent_architecture,
               struct pfam *first_pfam_ptr,struct parameters *params)
{
  char chain, ec_number[LINELEN + 1], number_string[20], pdb_code[5];
  char cov_code, coverage[3];

  char *architecture;
  char *domain[MAX_TOKENS], *string_ptr, *token[MAX_TOKENS];

  int iarch, ipos, istart, len, narchit, nseqs, numpos, parent_id;
  int itoken, ndomains, nhits, npdb, nstart, ntokens, total_hits;
  int join_code;

  int done, have_ec;

  struct archit *archit_ptr;
  struct enzlist *enzlist_ptr;
  struct enzyme *enzyme_ptr;
  struct hit *hit_ptr;
  struct hitlist *hitlist_ptr, *next_hitlist_ptr, *prev_hitlist_ptr;
  struct pfam *pfam_ptr;

  /* Initialise variables */
  architecture = &null;
  iarch = total_hits = parent_id = 0;
  if (params->get_ssg == TRUE)
    numpos = 6;
  else
    numpos = 7;
  narchit = 0;
  nhits = 0;
  npdb = 0;

  /* Write out numbr of hits */
  if (params->output_type == ARCHSCHEMA_OUTPUT)
    {
      fprintf(outfile_ptr,":HITS\n");
      fprintf(outfile_ptr,":HIT_LIST\n");
    }

  /* Get the first architecture record */
  archit_ptr = first_archit_ptr;

  /* Loop over the architecture records to write out */
  while (archit_ptr != NULL)
    {
      /* Increment architecture count */
      narchit++;

      /* If this is the parent architecture, store its identifier */
      if (!strcmp(archit_ptr->architecture,parent_architecture))
        parent_id = iarch;

      /* Write out this architecture and its total number of sequences */
      if (params->output_type == ARCHSCHEMA_OUTPUT)
        {
          fprintf(outfile_ptr,":A\n");
          fprintf(outfile_ptr,"%s\t%d\t%d\n",archit_ptr->architecture,
                  archit_ptr->nhits,archit_ptr->score);
        }
      else
        {
          /* Save the architecture */
          //OUT if (architecture != &null)
            //OUT free(architecture);
          store_string(&architecture,archit_ptr->architecture);

          /* Break up the architecture into its constituent domains */ 
          ndomains = tokenize(architecture,domain,MAX_TOKENS,'.');
        }

      /* Initialise flag */
      have_ec = FALSE;
      nseqs = 0;
      nstart = nhits;

      /* Get the pointer to the first of the hits */
      hitlist_ptr = archit_ptr->first_hitlist_ptr;

      /* Loop through the hits to write out */
      while (hitlist_ptr != NULL)
        {
          /* Get the corresponding hit record */
          hit_ptr = hitlist_ptr->hit_ptr;

          /* Increment hit count */
          total_hits++;
          nseqs++;
          nhits++;

          /* Write out the hit details */
          if (params->output_type == ARCHSCHEMA_OUTPUT)
            fprintf(outfile_ptr,"%s\t%s\t%s\t%d\t%s\t%s\n",
                    hit_ptr->uniprot_acc,
                    hit_ptr->uniprot_id,
                    hit_ptr->protein_name,
                    hit_ptr->seq_len,
                    hit_ptr->pdb_codes,
                    hit_ptr->coverage);
          else
            {
              /* Write out the architecture */
              fprintf(outfile_ptr,"ARCHITECTURE[%d] %s\n",nhits,
                      archit_ptr->architecture);
              fprintf(outfile_ptr,"ARCHITECTURE_SCORE[%d] %d\n",nhits,
                      archit_ptr->score);

              /* If not the first sequence, write out a ditto flag */
              if (nseqs == 1)
                fprintf(outfile_ptr,"NEW[%d] TRUE\n",nhits);

              /* Write out this sequence */
              fprintf(outfile_ptr,"UNIPROT_ACC[%d] %s\n",nhits,
                      hit_ptr->uniprot_acc);
              fprintf(outfile_ptr,"UNIPROT_ID[%d] %s\n",nhits,
                      hit_ptr->uniprot_id);
              fprintf(outfile_ptr,"PROTEIN_NAME[%d] %s\n",nhits,
                      hit_ptr->protein_name);
              fprintf(outfile_ptr,"SEQ_LEN[%d] %d\n",nhits,
                      hit_ptr->seq_len);
              fprintf(outfile_ptr,"COVERAGE[%d] %s\n",nhits,
                      hit_ptr->coverage);

              /* Write out the number of domains in current architecture */
              fprintf(outfile_ptr,"NDOMAINS[%d] %d\n",nhits,ndomains);

              /* Loop over the domains to write each one out */
              for (itoken = 0; itoken < ndomains; itoken++)
                {
                  /* Write the domain id */
                  fprintf(outfile_ptr,"ARCH_PFAM_ID[%d][%d] %s\n",nhits,
                          itoken + 1,domain[itoken]);

                  /* Write out Pfam-B flag */
                  if (!strncmp(domain[itoken],"PF",2))
                    {
                      fprintf(outfile_ptr,"ARCH_PFAMB[%d][%d] FALSE\n",nhits,
                              itoken + 1);
                      fprintf(outfile_ptr,"ARCH_PFAM_CAT[%d][%d] A\n",nhits,
                              itoken + 1);
                    }
                  else
                    {
                      fprintf(outfile_ptr,"ARCH_PFAMB[%d][%d] TRUE\n",nhits,
                              itoken + 1);
                      fprintf(outfile_ptr,"ARCH_PFAM_CAT[%d][%d] B\n",nhits,
                              itoken + 1);
                    }

                  /* Find the corresponding Pfam record to get domain name */
                  pfam_ptr = find_pfam(first_pfam_ptr,domain[itoken]);
                  if (pfam_ptr != NULL)
                    {
                      /* Write out this domain's details */
                      fprintf(outfile_ptr,"ARCH_PFAM_NAME[%d][%d] %s\n",
                              nhits,itoken + 1,pfam_ptr->pfam_name);
                      fprintf(outfile_ptr,"ARCH_PFAM_SHORT[%d][%d] %s\n",
                              nhits,itoken + 1,pfam_ptr->short_name);
                      fprintf(outfile_ptr,"ARCH_PFAM_COLOUR[%d][%d] %d\n",
                              nhits,itoken + 1,pfam_ptr->colour);
                    }

                  /* Get this domain's coverage */
                  strncpy(coverage,hit_ptr->coverage + 2 * itoken,2);
                  coverage[2] = '\0';

                  /* Get the coverage identifier */
                  cov_code = coverage[0];
                  if (cov_code != 'A' && cov_code != 'P' &&
                      cov_code != 'F')
                    cov_code = 'X';

                  /* Determine what the join to next domain is */
                  if (coverage[1] == '-')
                    join_code = 1;
                  else
                    join_code = 0;

                  /* Write out the coverage codes */
                  fprintf(outfile_ptr,"ARCH_COVERAGE[%d][%d] %c\n",
                          nhits,itoken + 1,cov_code);
                  fprintf(outfile_ptr,"ARCH_JOIN[%d][%d] %d\n",
                          nhits,itoken + 1,join_code);
                }

              /* List the individual PDB codes */
              fprintf(outfile_ptr,"PDB_CODES[%d] %s\n",nhits,
                      hit_ptr->pdb_codes);
              ntokens
                = tokenize(hit_ptr->pdb_codes,token,MAX_TOKENS,' ');
              fprintf(outfile_ptr,"ARCH_NPDB[%d] %d\n",nhits,ntokens);
              npdb = npdb + ntokens;

              /* Loop over the PDB codes to write each one out */
              for (itoken = 0; itoken < ntokens; itoken++)
                {
                  /* Get the PDB code and chain id */
                  len = strlen(token[itoken]);
                  strncpy(pdb_code,token[itoken],4);
                  pdb_code[4] = '\0';
                  if (len > 4)
                    chain = token[itoken][4];
                  else
                    chain = ' ';

                  /* Write the PDB code */
                  fprintf(outfile_ptr,"ARCH_PDB[%d][%d] %s\n",nhits,
                          itoken + 1,pdb_code);

                  /* Write the chain id */
                  fprintf(outfile_ptr,"ARCH_CHAIN[%d][%d] %c\n",nhits,
                          itoken + 1,chain);
                }
            }

          /* Get the next hit record */
          hitlist_ptr = hitlist_ptr->next_hitlist_ptr;
        }

      /* Write out the number of hits */
      if (params->output_type == DAT_OUTPUT)
        fprintf(outfile_ptr,"ARCHITECTURE_NSEQS[%d] %d\n",nstart,
                nseqs);

      /* If any of the sequences were enzymes, then write out list
         of E.C. numbers encountered */
      if (archit_ptr->first_enzlist_ptr != NULL &&
          params->output_type == ARCHSCHEMA_OUTPUT)
        {
          /* Write out the enzyme header record */
          fprintf(outfile_ptr,":E\n");

          /* Get the pointer to the first enzlist record */
          enzlist_ptr = archit_ptr->first_enzlist_ptr;

          /* Loop over all the enzymes to list this that fund for current
             architecture */
          while (enzlist_ptr != NULL)
            {
              /* If any sequences were annotated for this code, then
                 write out */
              if (enzlist_ptr->nseqs > 0)
                {
                  /* Get the corresponding enzyme record */
                  enzyme_ptr = enzlist_ptr->enzyme_ptr;

                  /* Write out */
                  fprintf(outfile_ptr,"%s\t%d\n",enzyme_ptr->ec_number,
                          enzlist_ptr->nseqs);

                  /* Re-initialise the count */
                  enzlist_ptr->nseqs = 0;
                }

              /* Get pointer to the next enzyme */
              enzlist_ptr = enzlist_ptr->next_enzlist_ptr;
            }
        }

      /* Get the next architecture record */
      archit_ptr = archit_ptr->next_archit_ptr;
      iarch++;
    }

  /* Write out total number of architectures */
  if (params->output_type == DAT_OUTPUT)
    {
      fprintf(outfile_ptr,"NARCHIT %d\n",narchit);
      fprintf(outfile_ptr,"NSEQS %d\n",nhits);
    }

  /* Write out the total number of PDB structures */
  if (params->output_type == DAT_OUTPUT)
    fprintf(outfile_ptr,"NPDB %d\n",npdb);

  /* Return the architecture position of the parent architecture */
  return(parent_id);
}
/***********************************************************************

species_sort  -  Shell-sort to arrange the species records in ascending
                 order

***********************************************************************/

void species_sort(struct species **species_index_ptr,int nspecies)
{
  int match;
  long endloop, i, indi, indl, j, k, l, lognb2, m, n, nn;
  float aln2i = (float) 1.4426950, tiny = (float) 1.e-5;
  float calc;
  double dble, dble1;

  struct species *species1_ptr, *species2_ptr, *swap_species_ptr;

  /*--Initialise values */
  n = nspecies;
  dble = n;
  dble1 = log(dble);
  calc = (float) dble1;
  calc = calc * aln2i + tiny;
  lognb2 = (long) calc;

  /* Perform sort */
  m = n;
  for (nn = 1; nn <= lognb2; nn++)
    {
      m = m / 2;
      k = n - m;
      for (j = 1; j <= k; j++)
        {
          i = j;
          endloop = FALSE;
          while (endloop == FALSE)
            {
              l = i + m;
              indi = i - 1;
              indl = l - 1;

              /* Get the two code records */
              species1_ptr = species_index_ptr[indi];
              species2_ptr = species_index_ptr[indl];

              /* Determine which string is greater */
              match
                = strcmp(species1_ptr->species_id,species2_ptr->species_id);

              /* If in wrong order, then swap sort-pointers */
              if (match > 0)
                {
                  swap_species_ptr = species_index_ptr[indi];
                  species_index_ptr[indi] = species_index_ptr[indl];
                  species_index_ptr[indl] = swap_species_ptr;

                  i = i - m;
                  if (i < 1)
                    endloop = TRUE;
                }
              else
                endloop = TRUE;
            }
        }
    }
}
/***********************************************************************

sort_species_records  -  Sort the Species domains into search order

***********************************************************************/

int sort_species_records(struct species **fst_species_ptr,int nspecies,
                         struct parameters *params)
{
  char last_species_id[20];

  int iorder, nentry;

  struct species *species_ptr, *first_species_ptr, *last_species_ptr;

  struct species **species_index_ptr;

  /* Initialise variables */
  first_species_ptr = *fst_species_ptr;
  last_species_ptr = NULL;
  last_species_id[0] = '\0';

  /* Create the species array as one large array */
  species_index_ptr
    = (struct species **) malloc(nspecies * sizeof(struct species *));

  /* Check that sufficient memory was allocated */
  if (species_index_ptr == NULL)
    write_error_message(MEMORY_ALLOCATION,
                        "Cannot allocate memory for species array",
                        params->output_type);

  /* Get pointer to the first code entry */
  species_ptr = first_species_ptr;
  nentry = 0;

  /* Loop through all entries to place each one in the species array */
  while (species_ptr != NULL)
    {
      /* Store current pointer at end of the species array */
      if (nentry < nspecies)
        species_index_ptr[nentry] = species_ptr;

      /* Increment count of stored entries */
      nentry++;

      /* Get pointer to the next entry in the linked list */
      species_ptr = species_ptr->next_species_ptr;
    }

  /* Sort the code pointers into order of primary species */
  species_sort(species_index_ptr,nentry);

  /* Re-intialise species count */
  nspecies = 0;

  /* Having sorted the records, rearrange order of the pointers */
  for (iorder = 0; iorder < nentry; iorder++)
    {
      /* Get the current species record */
      species_ptr = species_index_ptr[iorder];

      /* Initialise pointer to next record */
      species_ptr->next_species_ptr = NULL;

      /* Increment count of species records */
      nspecies++;

      /* If this is the first, then make the first record */
      if (iorder == 0)
        first_species_ptr = last_species_ptr = species_ptr;

      /* Otherwise, check if we even want this one */
      else
        {
          /* If this is the same species-id as the previous record,
             then don't need this duplicate */
          if (!strcmp(last_species_ptr->species_id,species_ptr->species_id))
            {
              /* Increment previous record's sequence count */
              last_species_ptr->nseqs++;

              /* Free up the memory taken by the current record */
              //OUT free(species_ptr->species_id);
              free(species_ptr);

              /* Reduce species count */
              nspecies--;
            }

          /* Otherwise, we want this record */
          else
            {
              /* Blank out current record's next species pointer */
              species_ptr->next_species_ptr = NULL;

              /* Get the last species record to point to this one */
              last_species_ptr->next_species_ptr = species_ptr;

              /* Store current record as last encountered */
              last_species_ptr = species_ptr;
            }
        }
    }

  /* Free up the index array created here */
  free(species_index_ptr);

  /* Return the first species record */
  *fst_species_ptr = first_species_ptr;

  /* Return new count of species */
  return(nspecies);
}
/***********************************************************************

get_species_names  -  Read in the speices names from speclist.out

***********************************************************************/

void get_species_names(struct species *first_species_ptr,
                       char *archschema_dir,struct parameters *params)
{
  char file_name[FILENAME_LEN], input_line[LINELEN + 1];
  char last_species_id[10], species_name[NAME_LEN], message[LINELEN + 1];

  char *string_ptr;

  int len, nlines, nmatched;
  int done, next_code, next_record;

  struct species *species_ptr;

  FILE *file_ptr;

  /* Initialise variables */
  done = FALSE;
  last_species_id[0] = '\0';
  next_code = TRUE;
  next_record = TRUE;
  nlines = nmatched = 0;

  /* Get pointer to the first species record */
  species_ptr = first_species_ptr;
  next_code = FALSE;
  if (species_ptr == NULL)
    return;

  /* Form name of input file */
  strcpy(file_name,archschema_dir);
  strcat(file_name,"/");
  strcat(file_name,"speclist.out");

  /* Open the input data file */
  if ((file_ptr = fopen(file_name,"r")) == NULL)
    {
      /* Otherwise, error must be due to missing file */
      sprintf(message,"Unable to open data file [%s]\n",file_name);
      write_error_message(FILE_NOT_FOUND,message,params->output_type);
    }

  /* Loop until all codes have been assigned a protein name */
  while (done == FALSE)
    {
      /* If need to go to the next code, then do so */
      if (next_code == TRUE)
        {
          /* Get pointer to the next species record */
          if (species_ptr != NULL)
            species_ptr = species_ptr->next_species_ptr;
          else
            done = TRUE;

          /* Check for end */
          if (species_ptr == NULL)
            done = TRUE;

          /* Set flag */
          next_code = FALSE;
        }

      /* If need to read to next record, then loop until we have the
         next useful one */
      while (done == FALSE && next_record == TRUE)
        {
          /* Read in the next record from the file */
          if (fgets(input_line,LINELEN,file_ptr) != NULL)
            {
              /* Increment line count */
              nlines++;

              /* Get the species id */
              string_ptr = strchr(input_line,'\t');

              /* If have a tab, terminate the string */
              if (string_ptr != NULL)
                {
                  string_ptr[0] = '\0';
                  len = strlen(input_line);

                  /* Check that species id not too long */
                  if (len > SPECIES_LEN - 1)
                    {
                      sprintf(message,"Invalid Species id: [%s]\n",
                              input_line);
                      write_error_message(FILE_READ_ERROR,message,
                                          params->output_type);
                    }

                  /* Check that species ids in the file are in sorted
                     order */
                  if (strcmp(input_line,last_species_id) < 0)
                    {
                      sprintf(message,"Codes not sorted. "
                              "Current [%s] lt previous [%s]\n",
                              input_line,last_species_id);
                      write_error_message(FILE_READ_ERROR,message,
                                          params->output_type);
                    }

                  /* Save the species name */
                  strcpy(species_name,string_ptr + 1);

                  /* Truncate the name */
                  string_chop(species_name);

                  /* If this is less than our current species record,
                     need to keep reading */
                  if (strcmp(input_line,species_ptr->species_id) < 0)
                    next_record = TRUE;
                  else
                    next_record = FALSE;
                }

              /* Save the current species id */
              strcpy(last_species_id,input_line);
            }

          /* Have reached the end of the file, so we're done */
          else
            done = TRUE;
        }

      /* If not yet done, check the match */
      if (done == FALSE)
        {
          /* If this matches our current species record, store the name */
          if (!strcmp(input_line,species_ptr->species_id))
            {
              /* Store the name */
              store_string(&(species_ptr->species_name),species_name);

              /* Increment count of matches */
              nmatched++;

              /* Set both flags to go to the next record */
              next_code = TRUE;
              next_record = TRUE;
            }

          /* Otherwise, if it is higher, need to get the next code */
          else if (strcmp(input_line,species_ptr->species_id) > 0)
            {
              /* Set flags */
              next_code = TRUE;
              next_record = FALSE;
            }

          /* Otherwise, must be lower, so need to read the next
             record */
          else
              next_record = TRUE;
        }
    }

  /* Close the input data file */
  fclose(file_ptr);
}
/***********************************************************************

write_species  -  Write out all the species records

***********************************************************************/

void write_species(struct species *first_species_ptr,int nspecies)
{
  char species_name[LINELEN];

  int ispecies;

  struct species *species_ptr;

  /* Initialise variables */
  ispecies = 0;

  /* Start the species listing */
  fprintf(outfile_ptr,":SPECIES\n");
  /*
  fprintf(outfile_ptr,"<void method=\"add\">\n");
  fprintf(outfile_ptr,"<object class=\"java.util.HashMap\">\n");
  */

  /* Get the first species record */
  species_ptr = first_species_ptr;

  /* Loop over the species records to write out */
  while (species_ptr != NULL)
    {
      /* Increment species count */
      ispecies++;

      /* Check for missing species name */
      strcpy(species_name,species_ptr->species_name);
      if (species_name[0] == '\0')
        strcpy(species_name,"Unknown");

      /* Write out the species details */
      fprintf(outfile_ptr,"%s\t%s\t%d\n",
              species_ptr->species_id,species_name,species_ptr->nseqs);

      /*
      fprintf(outfile_ptr,"<void method=\"put\">\n");
      fprintf(outfile_ptr,"<string>%s</string>\n",
              species_ptr->species_id);
      fprintf(outfile_ptr,"<object class=\"java.util.HashMap\">\n");
      */

      /* Species name */
      /*
      fprintf(outfile_ptr,"<void method=\"put\">\n");
      fprintf(outfile_ptr,"<string>name</string>\n");
      fprintf(outfile_ptr,"<string>%s</string>\n",
              species_ptr->species_name);
      fprintf(outfile_ptr,"</void>\n");
      */

      /* Number of sequences */
      /*
      fprintf(outfile_ptr,"<void method=\"put\">\n");
      fprintf(outfile_ptr,"<string>nseqs</string>\n");
      fprintf(outfile_ptr,"<string>%d</string>\n",species_ptr->nseqs);
      fprintf(outfile_ptr,"</void>\n");
      */

      /* Close off this species entry */
      /*
      fprintf(outfile_ptr,"</object>\n");
      fprintf(outfile_ptr,"</void>\n");
      */

      /* Get the next species record */
      species_ptr = species_ptr->next_species_ptr;
    }

  /* Close off species list */
  /*
  fprintf(outfile_ptr,"</object>\n");
  fprintf(outfile_ptr,"</void>\n");
  */
}
/***********************************************************************

delete_pfam_records  -  Delete the duplicate Pfam records

***********************************************************************/

int delete_pfam_records(struct pfam **fst_pfam_ptr)
{
  char last_pfam_id[PFAM_LEN];

  int ndeleted, npfam;

  struct pfam *pfam_ptr, *next_pfam_ptr, *prev_pfam_ptr;
  struct pfam *first_pfam_ptr;

  /* Initialise variables */
  first_pfam_ptr = *fst_pfam_ptr;
  last_pfam_id[0] = '\0';
  ndeleted = 0;
  npfam = 0;

  /* Initialise pointers */
  prev_pfam_ptr = NULL;

  /* Get the first pfam */
  pfam_ptr = first_pfam_ptr;

  /* Loop through the pfams to see which are to be deleted */
  while (pfam_ptr != NULL)
    {
      /* Get the next pfam */
      next_pfam_ptr = pfam_ptr->next_pfam_ptr;

      /* Set its number of sequences to 1 */
      pfam_ptr->nseqs = 1;

      /* If Pfam id is the same as previous record, then delete */
      if (!strcmp(pfam_ptr->pfam_id,last_pfam_id))
        {
          /* If have a previous pfam, then get it to point
             around this one */
          if (prev_pfam_ptr != NULL)
            {
              prev_pfam_ptr->next_pfam_ptr = next_pfam_ptr;

              /* Increment its count of sequences */
              if (pfam_ptr->copy == FALSE)
                prev_pfam_ptr->nseqs++;
            }

          /* Otherwise, this must be the first pfam, so make
             the next one the new first */
          else
            first_pfam_ptr = next_pfam_ptr;

          /* Free up the memory taken by the deleted pfam */
          free(pfam_ptr);

          /* Increment count of deleted pfams */
          ndeleted++;
        }

      /* Otherwise, save as the new previous pfam */
      else
        {
          /* Save current pfam */
          prev_pfam_ptr = pfam_ptr;

          /* Save the Pfam id */
          strcpy(last_pfam_id,pfam_ptr->pfam_id);

          /* Increment count of retained pfams */
          npfam++;
        }

      /* Get the next pfam */
      pfam_ptr = next_pfam_ptr;
    }

  /* Return new starting pfam */
  *fst_pfam_ptr = first_pfam_ptr;

  /* Return number of retained pfams */
  return(npfam);
}
/***********************************************************************

annotate_sequence  -  Find the given UniProt sequence and annotate its
                      name to add SSG info

***********************************************************************/

int annotate_sequence(struct archit *first_archit_ptr,char *uniprot_acc,
                      char *ssg_class)
{
  char protein_name[LINELEN];

  int done;

  struct archit *archit_ptr;
  struct hit *hit_ptr;
  struct hitlist *hitlist_ptr;

  /* Initialise variables */
  done = FALSE;

  /* Get the first architecture record */
  archit_ptr = first_archit_ptr;

  /* Loop over the architecture records */
  while (archit_ptr != NULL && done == FALSE)
    {
      /* Get the pointer to the first of this architecture's hits */
      hitlist_ptr = archit_ptr->first_hitlist_ptr;

      /* Loop through the hits to find our sequence */
      while (hitlist_ptr != NULL && done == FALSE)
        {
          /* Get the corresponding hit record */
          hit_ptr = hitlist_ptr->hit_ptr;

          /* If this is the sequence we want, annotate its name */
          if (!strcmp(hit_ptr->uniprot_acc,uniprot_acc))
            {
              /* Get the protein name */
              strcpy(protein_name,hit_ptr->protein_name);

              /* Free up memory */
              //OUT if (hit_ptr->protein_name != &null)
                //OUT free(hit_ptr->protein_name);

              /* Append the SSG class */
              strcat(protein_name," [SSG ");
              strcat(protein_name,ssg_class);
              strcat(protein_name,"]");

              /* Store in the hit record */
              store_string(&(hit_ptr->protein_name),protein_name);

              /* Set flag that we're done */
              done = TRUE;
            }

          /* Get the next hit record */
          hitlist_ptr = hitlist_ptr->next_hitlist_ptr;
        }

      /* Get the next architecture record */
      archit_ptr = archit_ptr->next_archit_ptr;
    }

  /* Return whether done */
  return(done);
}
/***********************************************************************

read_ssg_data  -  Read in the SSG assignments for each UniProt accession
                  code

***********************************************************************/

int read_ssg_data(struct archit *first_archit_ptr,
                  struct c_file_data file_name_ptr)
{
  char file_name[FILENAME_LEN], input_line[LINELEN + 1];
  char ssg_class[20], uniprot_acc[UNIPROT_LEN];

  char *token[MAX_TOKENS];

  int len, ntokens;
  int ok;

  FILE *file_ptr;

  /* Initialise variables */
  ok = TRUE;

  /* Form name of input file */
  strcpy(file_name,file_name_ptr.other_dir[CATH_DATA_DIR]);
  strcat(file_name,"/ssg.txt");

  /* Open the input data file */
  if ((file_ptr = fopen(file_name,"r")) != NULL)
    {
      /* Read in the next record from the file */
      while (fgets(input_line,LINELEN,file_ptr) != NULL)
        {
          /* Truncate the string */
          len = string_truncate(input_line,LINELEN);
          
          /* Extract the UniProt accession and SSG cluster */
          ntokens = tokenize(input_line,token,MAX_TOKENS,' ');

          /* If have 2 tokens, process them */
          if (ntokens == 2)
            {
              /* First token is the UniProt accession code */
              strcpy(uniprot_acc,token[0]);

              /* Second token is the SSG class */
              strcpy(ssg_class,token[1]);

              /* Find the sequence and append the SSG class to it */
              annotate_sequence(first_archit_ptr,uniprot_acc,ssg_class);
            }
        }

      /* Close the input file */
      fclose(file_ptr);
    }
  else
    ok = FALSE;

  /* Return whether file read OK */
  return(ok);
}
/***********************************************************************

remove_partials  -  Remove any redundant partial E.C. assignments

***********************************************************************/

void remove_partials(char *ec_number)
{
  char tmp[LINELEN + 1];

  char *string_ptr, *token[MAX_TOKENS];

  int itoken, jtoken, len, ntokens;
  int blank[MAX_TOKENS], blanked, done;

  /* Initialise variables */
  blanked = FALSE;
  for (itoken = 0; itoken < MAX_TOKENS; itoken++)
    blank[itoken] = FALSE;

  /* Get the space-delimited tokens on this line */
  strcpy(tmp,ec_number);
  ntokens = tokenize(tmp,token,MAX_TOKENS,' ');

  /* Loop over the tokens to see which ones can be removed */
  for (itoken = 0; itoken < ntokens; itoken++)
    {
      /* If this is a partial, compare against all the others */
      string_ptr = strchr(token[itoken],'-');

      /* If a partial, then get part to be compared */
      if (string_ptr != NULL)
        {
          /* Get length to be compared */
          len = string_ptr - token[itoken];
          done = FALSE;

          /* Loop over all the other E.C. tokens */
          for (jtoken = 0; jtoken < ntokens && done == FALSE; jtoken++)
            {
              /* If not the current E.C. number, but matches, then
                 blank out current */
              if (jtoken != itoken &&
                  !strncmp(token[jtoken],token[itoken],len) &&
                  blank[jtoken] == FALSE)
                {
                  /* Blank out the current token */
                  blank[itoken] = TRUE;
                  blanked = TRUE;

                  /* Set flag that we're done */
                  done = TRUE;
                }

            }
        }
    }

  /* If any tokens blanked out, then reform the E.C number without them */
  if (blanked == TRUE)
    {
      /* Initialise E.C. number string */
      ec_number[0] = '\0';

      /* Loop over the tokens to place back any that weren't blanked out */
      for (itoken = 0; itoken < ntokens; itoken++)
        {
          /* If token not blank, then place back */
          if (blank[itoken] == FALSE)
            {
              /* If not the first, then append a space */
              if (ec_number[0] != '\0')
                strcat(ec_number," ");
              strcat(ec_number,token[itoken]);
            }              
        }

      /* Get the length of the E.C. number */
      len = strlen(ec_number);

      /* If last character is a comma, then remove it */
      if (ec_number[len - 1] == ',')
        ec_number[len - 1] = '\0';
    }
}
/***********************************************************************

find_enzlist  -  Find the given enzlist entry

***********************************************************************/

struct enzlist *find_enzlist(struct enzlist *first_enzlist_ptr,
                             struct enzyme *enzyme_ptr)
{
  struct enzlist *enzlist_ptr, *found_enzlist_ptr;

  /* Initialise variables */
  found_enzlist_ptr = NULL;

  /* Get pointer to the first enzlist */
  enzlist_ptr = first_enzlist_ptr;

  /* Loop over all the enzlists to find the one of interest */
  while (enzlist_ptr != NULL && found_enzlist_ptr == NULL)
    {
      /* If this is the right enzlist, store pointer */
      if (enzyme_ptr == enzlist_ptr->enzyme_ptr)
        found_enzlist_ptr = enzlist_ptr;

      /* Get pointer to the next enzlist */
      enzlist_ptr = enzlist_ptr->next_enzlist_ptr;
    }

  /* Return the enzlisting enzlist pointer */
  return(found_enzlist_ptr);
}
/***********************************************************************

create_enzlist_record  -  Create and initialise a new enzlist record

***********************************************************************/

struct enzlist *create_enzlist_record(struct enzlist **fst_enzlist_ptr,
                                      struct enzlist **lst_enzlist_ptr,
                                      struct enzyme *enzyme_ptr,
                                      struct parameters *params)
{
  struct enzlist *enzlist_ptr, *first_enzlist_ptr, *last_enzlist_ptr;

  /* Initialise pointers */
  enzlist_ptr = NULL;
  first_enzlist_ptr = *fst_enzlist_ptr;
  last_enzlist_ptr = *lst_enzlist_ptr;

  /* Allocate memory for structure */
  enzlist_ptr = (struct enzlist *) malloc(sizeof(struct enzlist));
  if (enzlist_ptr == NULL)
    write_error_message(MEMORY_ALLOCATION,
                        "Cannot allocate memory for struct enzlist",
                        params->output_type);

  /* If this is the very first, save its pointer */
  if (first_enzlist_ptr == NULL)
    first_enzlist_ptr = enzlist_ptr;

  /* Add link from previous record to the current one */
  if (last_enzlist_ptr != NULL)
    last_enzlist_ptr->next_enzlist_ptr = enzlist_ptr;
  last_enzlist_ptr = enzlist_ptr;

  /* Store known elements of the structure */
  enzlist_ptr->enzyme_ptr = enzyme_ptr;

  /* Initialise remaining elements */
  enzlist_ptr->nseqs = 0;

  /* Initialise pointer to next record */
  enzlist_ptr->next_enzlist_ptr = NULL;

  /* Return current pointers */
  *fst_enzlist_ptr = first_enzlist_ptr;
  *lst_enzlist_ptr = last_enzlist_ptr;

  /* Return new pointer */
  return(enzlist_ptr);
}
/***********************************************************************

get_enzyme_data  -  Pick up all gthe E.C. numbers attached to each
                    architecture from the sequences belonging to it

***********************************************************************/

int get_enzyme_data(struct archit *first_archit_ptr,
                    struct enzyme **fst_enzyme_ptr,
                    struct parameters *params)
{
  char ec_number[LINELEN + 1], number_string[20];

  char *string_ptr;

  int iarch, ihit, ipos, istart, len, nec, numpos, parent_id;

  int done, have_ec;

  struct archit *archit_ptr;
  struct enzyme *enzyme_ptr, *first_enzyme_ptr, *last_enzyme_ptr;
  struct enzlist *enzlist_ptr, *last_enzlist_ptr;
  struct hit *hit_ptr;
  struct hitlist *hitlist_ptr, *next_hitlist_ptr, *prev_hitlist_ptr;

  /* Initialise variables */
  iarch = ihit = nec = parent_id = 0;
  if (params->get_ssg == TRUE)
    numpos = 6;
  else
    numpos = 7;

  /* Initialise pointers */
  enzyme_ptr = first_enzyme_ptr = last_enzyme_ptr = NULL;

  /* Get the first architecture record */
  archit_ptr = first_archit_ptr;

  /* Loop over the architecture records to process their sequences */
  while (archit_ptr != NULL)
    {
      /* Initialise flag */
      have_ec = FALSE;

      /* Initialise enzyme list pointer */
      last_enzlist_ptr = NULL;

      /* Get the pointer to the first of the hits */
      hitlist_ptr = archit_ptr->first_hitlist_ptr;

      /* Loop through the hits to get E.C. assignments */
      while (hitlist_ptr != NULL)
        {
          /* Get the corresponding hit record */
          hit_ptr = hitlist_ptr->hit_ptr;

          /* Increment hit count */
          ihit++;

          /* Check for E.C./SSG data, identifier in protein name */
          string_ptr = NULL;
          if (params->get_ssg == TRUE)
            string_ptr = strstr(hit_ptr->protein_name," [SSG ");
          else
            string_ptr = strstr(hit_ptr->protein_name," (E.C. ");

          /* Extract any E.C. numbers from protein name and update
             enzyme records */
          if (string_ptr != NULL)
            {
              /* Set start position */
              ipos = string_ptr - hit_ptr->protein_name + numpos;

              /* Initialise flag */
              done = FALSE;

              /* Loop until all E.C numbers have been extracted */
              while (done == FALSE)
                {
                  /* Extract this E.C. number */
                  strcpy(ec_number,hit_ptr->protein_name + ipos);
                  len = strlen(ec_number);

                  /* Find the end of the code */
                  string_ptr = strchr(ec_number,')');
                  if (string_ptr == NULL)
                    string_ptr = strchr(ec_number,']');

                  /* If have end, then process */
                  if (len > 1 && string_ptr != NULL)
                    {
                      /* Terminate code */
                      string_ptr[0] = '\0';

                      /* Set flag to indicate we have enzymes */
                      have_ec = TRUE;

                      /* Get the length of the code */
                      len = strlen(ec_number);

                      /* Increment position */
                      ipos = ipos + len + 2;

                      /* If have a list of E.C. numbers and any are
                         partial, remove any redundant partials */
                      if (strchr(ec_number,',') != NULL &&
                          strchr(ec_number,'-') != NULL)
                        remove_partials(ec_number);

                      /* See if we already have this E.C. number */
                      enzyme_ptr
                        = find_enzyme(first_enzyme_ptr,ec_number);

                      /* If this is a new code, create a record */
                      if (enzyme_ptr == NULL)
                        {
                          /* Create a enzyme record */
                          enzyme_ptr
                            = create_enzyme_record(&first_enzyme_ptr,
                                                   &last_enzyme_ptr,
                                                   ec_number,params);

                          /* Increment count of enzyme */
                          nec++;

                          /* If this is an SSG record, then create name
                             for it */
                          if (params->get_ssg == TRUE)
                            {
                              /* Form name and store */
                              sprintf(number_string,"SSG class %d",nec);
                              store_string(&(enzyme_ptr->enzyme_name),
                                           number_string);

                              /* Set flag that we're done */
                              done = TRUE;
                            }
                        }

                      /* See if we have this enzyme assigned to this
                         architecture */
                      /* See if we already have this E.C. number */
                      enzlist_ptr
                        = find_enzlist(archit_ptr->first_enzlist_ptr,
                                       enzyme_ptr);

                      /* If this is a new code, create a record */
                      if (enzlist_ptr == NULL)
                        {
                          /* Create an enzlist record */
                          enzlist_ptr
                            = create_enzlist_record(&(archit_ptr->first_enzlist_ptr),
                                                    &last_enzlist_ptr,
                                                    enzyme_ptr,params);

                          /* Increment count of architectures associated
                             with the current E.C. code */
                          enzyme_ptr->narch++;
                        }

                      /* Update sequence count */
                      enzlist_ptr->nseqs++;

                      /* Update count of architectures assigned to this
                         E.C. class */
                      enzyme_ptr->nseqs++;
                      enzyme_ptr->ntotal++; 
                   }

                  /* Otherwise, we are done */
                  else
                    done = TRUE;
                }
            }

          /* Get the next hit record */
          hitlist_ptr = hitlist_ptr->next_hitlist_ptr;
        }

      /* Get the next architecture record */
      archit_ptr = archit_ptr->next_archit_ptr;
      iarch++;
    }

  /* Return first enzyme pointer and number of EC codes */
  *fst_enzyme_ptr = first_enzyme_ptr;

  /* Return number of enzyme records created */
  return(nec);
}
/***********************************************************************

get_pfam_names  -  Read in the pfam names from the
                   pfam_domains.parsed.sorted file

***********************************************************************/

void get_pfam_names(int type,struct pfam *first_pfam_ptr,
                    char *archschema_dir,struct parameters *params)
{
  char file_name[FILENAME_LEN], input_line[LINELEN + 1];
  char dom_type[6], last_pfam_id[PFAM_LEN], pfam_id[PFAM_LEN];
  char pfam_name[NAME_LEN], short_name[NAME_LEN], message[LINELEN + 1];

  char *string_ptr, *token[MAX_TOKENS];

  int itoken, len, nlines, ntokens, nmatched;
  int done, next_code, next_record;

  struct pfam *pfam_ptr;

  FILE *file_ptr;

  /* Initialise variables */
  if (params->cath_domains == FALSE)
    strcpy(dom_type,"Pfam");
  else
    strcpy(dom_type,"CATH");
  done = FALSE;
  last_pfam_id[0] = '\0';
  next_code = TRUE;
  next_record = TRUE;
  nlines = nmatched = 0;

  /* Get pointer to the first pfam record */
  pfam_ptr = first_pfam_ptr;
  next_code = FALSE;
  if (pfam_ptr == NULL)
    return;

  /* Form name of input file */
  strcpy(file_name,archschema_dir);
  strcat(file_name,"/");
  if (type == PFAM)
    strcat(file_name,"pfam_domains.parsed.sorted");
  else
    strcat(file_name,"cath_domains.parsed.sorted");

  /* Open the input data file */
  if ((file_ptr = fopen(file_name,"r")) == NULL)
    {
      /* Otherwise, error must be due to missing file */
      sprintf(message,"Unable to open data file [%s]\n",file_name);
      write_error_message(FILE_NOT_FOUND,message,params->output_type);
    }

  /* Loop until all codes have been assigned a protein name */
  while (done == FALSE)
    {
      /* If need to go to the next code, then do so */
      if (next_code == TRUE)
        {
          /* Get pointer to the next pfam record */
          if (pfam_ptr != NULL)
            pfam_ptr = pfam_ptr->next_pfam_ptr;
          else
            done = TRUE;

          /* Check for end */
          if (pfam_ptr == NULL)
            done = TRUE;

          /* Set flag */
          next_code = FALSE;
        }

      /* If need to read to next record, then loop until we have the
         next useful one */
      while (done == FALSE && next_record == TRUE)
        {
          /* Read in the next record from the file */
          if (fgets(input_line,LINELEN,file_ptr) != NULL)
            {
              /* Increment line count */
              nlines++;

              /* Get the tab-delimited tokens on this line */
              ntokens = tokenize(input_line,token,MAX_TOKENS,'\t');

              /* If have 3 tokens, process them */
              if (ntokens == 3)
                {
                  /* Check that pfam id not too long */
                  len = strlen(token[0]);
                  if (len > PFAM_LEN - 1)
                    {
                      sprintf(message,"Invalid %s id: [%s]\n",
                              dom_type,input_line);
                      write_error_message(FILE_READ_ERROR,message,
                                          params->output_type);
                    }
                  strcpy(pfam_id,token[0]);

                  /* Check that pfam ids in the file are in sorted
                     order */
                  if (strcmp(pfam_id,last_pfam_id) < 0)
                    {
                      sprintf(message,"Codes not sorted. "
                              "Current [%s] lt previous [%s]\n",
                              pfam_id,last_pfam_id);
                    }

                  /* Save the long and short Pfam names */
                  strcpy(pfam_name,token[1]);
                  strcpy(short_name,token[2]);
                  string_chop(short_name);

                  /* If this is less than our current pfam record,
                     need to keep reading */
                  if (strcmp(pfam_id,pfam_ptr->pfam_id) < 0)
                    next_record = TRUE;
                  else
                    next_record = FALSE;
                }

              /* Save the current pfam id */
              strcpy(last_pfam_id,pfam_id);
            }

          /* Have reached the end of the file, so we're done */
          else
            done = TRUE;
        }

      /* If not yet done, check the match */
      if (done == FALSE)
        {
          /* If this matches our current pfam record, store the name */
          if (!strcmp(pfam_id,pfam_ptr->pfam_id))
            {
              /* Store the name */
              store_string(&(pfam_ptr->pfam_name),pfam_name);
              store_string(&(pfam_ptr->short_name),short_name);

              /* Increment count of matches */
              nmatched++;

              /* Set both flags to go to the next record */
              next_code = TRUE;
              next_record = TRUE;
            }

          /* Otherwise, if it is higher, need to get the next code */
          else if (strcmp(pfam_id,pfam_ptr->pfam_id) > 0)
            {
              /* Set flags */
              next_code = TRUE;
              next_record = FALSE;
            }

          /* Otherwise, must be lower, so need to read the next
             record */
          else
              next_record = TRUE;
        }
    }

  /* Close the input data file */
  fclose(file_ptr);
}
/***********************************************************************

assign_pfam_colours  -  Assign colours to each of the Pfam domains

***********************************************************************/

void assign_pfam_colours(struct pfam *first_pfam_ptr,
                         struct parameters *params)
{
  int colour[2], itoken, loop, ntokens, pstn, type;

  char *token[MAX_TOKENS];

  struct pfam *pfam_ptr;

  /* Initialise variables */
  colour[0] = 0;
  colour[1] = PFAMB_COLOUR_START;
  pstn = 0;

  /* Get the Pfam domains in the search architecture */
  ntokens = tokenize(params->architecture,token,MAX_TOKENS,'.');

  /* Loop over the domains to write each one out */
  for (itoken = 0; itoken < ntokens; itoken++)
    {
      /* Find the corresponding Pfam record */
      pfam_ptr = find_pfam(first_pfam_ptr,token[itoken]);

      /* If found and not yet assigned a colour, give it the next
         colour */
      if (pfam_ptr != NULL && pfam_ptr->colour == 0)
        {
          /* Determine whether this is a Pfam-A of Pfam-B domain */
          if (!strncmp(pfam_ptr->pfam_id,"PF",2))
            type = PFAM_A;
          else
            type = PFAM_B;

          /* Increment the corresponding colour number and store */
          colour[type]++;
          pfam_ptr->colour = colour[type];
          if (pfam_ptr->colour > MAX_DOMAIN_COLOURS)
            pfam_ptr->colour = 0;

          /* Increment position in domain list and store */
          pstn++;
          pfam_ptr->pstn = pstn;
        }
    }

  /* Loop twice to assign colours and positions (first loop assigns
     positions to Pfam-A domains, and second loop does Pfam-B domains) */
  for (loop = 0; loop < 2; loop++)
    {
      /* Get pointer to the first pfam */
      pfam_ptr = first_pfam_ptr;

      /* Loop over all the pfams to assign colour to those that still
         don't have them */
      while (pfam_ptr != NULL)
        {
          /* Determine whether this is a Pfam-A of Pfam-B domain */
          if (!strncmp(pfam_ptr->pfam_id,"PF",2))
            type = PFAM_A;
          else
            type = PFAM_B;

          /* If no colour assigned, then assign colour */
          if (pfam_ptr->colour == 0)
            {
              /* Increment the corresponding colour number and store */
              colour[type]++;
              pfam_ptr->colour = colour[type];
              if (pfam_ptr->colour > MAX_DOMAIN_COLOURS)
                pfam_ptr->colour = 0;
            }

          /* If no position assigned, then assign position if this
             is the right loop */
          if (pfam_ptr->pstn == 0 && type == loop)
            {
              /* Increment position in domain list and store */
              pstn++;
              pfam_ptr->pstn = pstn;
            }
  
          /* Get pointer to the next pfam */
          pfam_ptr = pfam_ptr->next_pfam_ptr;
        }
    }
}
/***********************************************************************

write_pfam  -  Write out all the pfam records

***********************************************************************/

void write_pfam(struct pfam *first_pfam_ptr,struct parameters *params)
{
  int ipfam, npfam;

  struct pfam *pfam_ptr;

  /* initialise variables */
  npfam = 0;

  /* Start the Pfam descriptions listing */
  if (params->output_type == ARCHSCHEMA_OUTPUT)
    fprintf(outfile_ptr,":PFAM\n");

  /*
  fprintf(outfile_ptr,"<void method=\"add\">\n");
  fprintf(outfile_ptr,"<object class=\"java.util.HashMap\">\n");
  */

  /* Get the first pfam record */
  pfam_ptr = first_pfam_ptr;

  /* Loop over the pfam records to write out */
  while (pfam_ptr != NULL)
    {
      /* Write out the pfam details */
      if (params->output_type == ARCHSCHEMA_OUTPUT)
        fprintf(outfile_ptr,"%s\t%s\t%s\t%d\n",
                pfam_ptr->pfam_id,
                pfam_ptr->short_name,
                pfam_ptr->pfam_name,
                pfam_ptr->nseqs);

      else
        {
          /* If this is the search domain, set subscript to 1 */
          ipfam = pfam_ptr->pstn;
          if (ipfam > npfam)
            npfam = ipfam;

          /* Write out the Pfam details */
          fprintf(outfile_ptr,"PFAM_ID[%d] %s\n",ipfam,
                  pfam_ptr->pfam_id);
          if (!strncmp(pfam_ptr->pfam_id,"PF",2))
            {
              fprintf(outfile_ptr,"PFAMB[%d] FALSE\n",ipfam);
              fprintf(outfile_ptr,"PFAM_CAT[%d] A\n",ipfam);
            }
          else
            {
              fprintf(outfile_ptr,"PFAMB[%d] TRUE\n",ipfam);
              fprintf(outfile_ptr,"PFAM_CAT[%d] B\n",ipfam);
            }
          fprintf(outfile_ptr,"PFAM_SHORT_NAME[%d] %s\n",ipfam,
                  pfam_ptr->short_name);
          fprintf(outfile_ptr,"PFAM_PFAM_NAME[%d] %s\n",ipfam,
                  pfam_ptr->pfam_name);
          fprintf(outfile_ptr,"PFAM_NSEQS[%d] %d\n",ipfam,
                  pfam_ptr->nseqs);
          fprintf(outfile_ptr,"PFAM_COLOUR[%d] %d\n",ipfam,
                  pfam_ptr->colour);
        }

      /*
      fprintf(outfile_ptr,"<void method=\"put\">\n");
      fprintf(outfile_ptr,"<string>%s</string>\n",pfam_ptr->pfam_id);
      fprintf(outfile_ptr,"<object class=\"java.util.HashMap\">\n");
      */

      /* Short Pfam name */
      /*
      fprintf(outfile_ptr,"<void method=\"put\">\n");
      fprintf(outfile_ptr,"<string>id</string>\n");
      fprintf(outfile_ptr,"<string>%s</string>\n",pfam_ptr->short_name);
      fprintf(outfile_ptr,"</void>\n");
      */

      /* Pfam description */
      /*
      fprintf(outfile_ptr,"<void method=\"put\">\n");
      fprintf(outfile_ptr,"<string>description</string>\n");
      fprintf(outfile_ptr,"<string>%s</string>\n",pfam_ptr->pfam_name);
      fprintf(outfile_ptr,"</void>\n");
      */

      /* Number of sequences */
      /*
      fprintf(outfile_ptr,"<void method=\"put\">\n");
      fprintf(outfile_ptr,"<string>nseqs</string>\n");
      fprintf(outfile_ptr,"<string>%d</string>\n",pfam_ptr->nseqs);
      fprintf(outfile_ptr,"</void>\n");
      */

      /* Close off this Pfam entry */
      /*
      fprintf(outfile_ptr,"</object>\n");
      fprintf(outfile_ptr,"</void>\n");
      */

      /* Get the next pfam record */
      pfam_ptr = pfam_ptr->next_pfam_ptr;
    }

  /* Write out number of Pfam domains */
  if (params->output_type == DAT_OUTPUT)
    fprintf(outfile_ptr,"NPFAM %d\n",npfam);
 
  /* Close off pfam list */
  /*
  fprintf(outfile_ptr,"</object>\n");
  fprintf(outfile_ptr,"</void>\n");
  */
}
/***********************************************************************

enzyme_sort1  -  Shell-sort to arrange the Enzyme records in ascending
                 order of E.C. number

***********************************************************************/

void enzyme_sort1(struct enzyme **enzyme_index_ptr,int nec)
{
  int match;
  long endloop, i, indi, indl, j, k, l, lognb2, m, n, nn;
  float aln2i = (float) 1.4426950, tiny = (float) 1.e-5;
  float calc;
  double dble, dble1;

  struct enzyme *enzyme1_ptr, *enzyme2_ptr, *swap_enzyme_ptr;

  /*--Initialise values */
  n = nec;
  dble = n;
  dble1 = log(dble);
  calc = (float) dble1;
  calc = calc * aln2i + tiny;
  lognb2 = (long) calc;

  /* Perform sort */
  m = n;
  for (nn = 1; nn <= lognb2; nn++)
    {
      m = m / 2;
      k = n - m;
      for (j = 1; j <= k; j++)
        {
          i = j;
          endloop = FALSE;
          while (endloop == FALSE)
            {
              l = i + m;
              indi = i - 1;
              indl = l - 1;

              /* Get the two code records */
              enzyme1_ptr = enzyme_index_ptr[indi];
              enzyme2_ptr = enzyme_index_ptr[indl];

              /* Determine which string is greater */
              match
                = strcmp(enzyme1_ptr->sort_number,
                         enzyme2_ptr->sort_number);

              /* If in wrong order, then swap sort-pointers */
              if (match > 0)
                {
                  swap_enzyme_ptr = enzyme_index_ptr[indi];
                  enzyme_index_ptr[indi] = enzyme_index_ptr[indl];
                  enzyme_index_ptr[indl] = swap_enzyme_ptr;

                  i = i - m;
                  if (i < 1)
                    endloop = TRUE;
                }
              else
                endloop = TRUE;
            }
        }
    }
}
/***********************************************************************

enzyme_sort2  -  Shell-sort to arrange the Enzyme records in descending
                 order of occurrence on the plot

***********************************************************************/

void enzyme_sort2(struct enzyme **enzyme_index_ptr,int nec)
{
  int swap;
  long endloop, i, indi, indl, j, k, l, lognb2, m, n, nn;
  float aln2i = (float) 1.4426950, tiny = (float) 1.e-5;
  float calc;
  double dble, dble1;

  struct enzyme *enzyme1_ptr, *enzyme2_ptr, *swap_enzyme_ptr;

  /*--Initialise values */
  n = nec;
  dble = n;
  dble1 = log(dble);
  calc = (float) dble1;
  calc = calc * aln2i + tiny;
  lognb2 = (long) calc;

  /* Perform sort */
  m = n;
  for (nn = 1; nn <= lognb2; nn++)
    {
      m = m / 2;
      k = n - m;
      for (j = 1; j <= k; j++)
        {
          i = j;
          endloop = FALSE;
          while (endloop == FALSE)
            {
              l = i + m;
              indi = i - 1;
              indl = l - 1;

              /* Get the two code records */
              enzyme1_ptr = enzyme_index_ptr[indi];
              enzyme2_ptr = enzyme_index_ptr[indl];

              /* Determine whether to swap these two records */
              swap = FALSE;
              if (enzyme1_ptr->narch < enzyme2_ptr->narch)
                swap = TRUE;
              else if (enzyme1_ptr->narch == enzyme2_ptr->narch)
                {
                  /* Check number of sequences */
                  if (enzyme1_ptr->ntotal < enzyme2_ptr->ntotal)
                    swap = TRUE;

                  /* If equal, sort on the E.C. number */
                  else if (enzyme1_ptr->ntotal == enzyme2_ptr->ntotal)
                    {
                      if (strcmp(enzyme1_ptr->sort_number,
                                 enzyme2_ptr->sort_number) > 0)
                        swap = TRUE;
                    }
                }

              /* If in wrong order, then swap sort-pointers */
              if (swap == TRUE)
                {
                  swap_enzyme_ptr = enzyme_index_ptr[indi];
                  enzyme_index_ptr[indi] = enzyme_index_ptr[indl];
                  enzyme_index_ptr[indl] = swap_enzyme_ptr;

                  i = i - m;
                  if (i < 1)
                    endloop = TRUE;
                }
              else
                endloop = TRUE;
            }
        }
    }
}
/***********************************************************************

sort_enzyme_records  -  Sort the Enzyme domains into search order

***********************************************************************/

void sort_enzyme_records(struct enzyme **fst_enzyme_ptr,int nec,
                         int type,struct parameters *params)
{
  char ch, chstring[2], number_string[20], sort_number[40];

  int inum, iorder, ipos, nentry;

  struct enzyme *enzyme_ptr, *first_enzyme_ptr, *last_enzyme_ptr;

  struct enzyme **enzyme_index_ptr;

  /* Initialise variables */
  first_enzyme_ptr = *fst_enzyme_ptr;
  last_enzyme_ptr = NULL;

  /* Create the enzyme array as one large array */
  enzyme_index_ptr = (struct enzyme **) malloc(nec * sizeof(struct enzyme *));

  /* Check that sufficient memory was allocated */
  if (enzyme_index_ptr == NULL)
    write_error_message(MEMORY_ALLOCATION,
                        "Cannot allocate memory for enzyme array",
                        params->output_type);

  /* Get pointer to the first code entry */
  enzyme_ptr = first_enzyme_ptr;
  nentry = 0;

  /* Loop through all entries to place each one in the enzyme array */
  while (enzyme_ptr != NULL)
    {
      /* Store current pointer at end of the enzyme array */
      if (nentry < nec)
        enzyme_index_ptr[nentry] = enzyme_ptr;

      /* Increment count of stored entries */
      nentry++;

      /* Get pointer to the next entry in the linked list */
      enzyme_ptr = enzyme_ptr->next_enzyme_ptr;
    }

  /* Sort the code pointers into order of primary enzyme */
  if (type == EC_CODE_SORT)
    enzyme_sort1(enzyme_index_ptr,nec);
  else
    enzyme_sort2(enzyme_index_ptr,nec);

  /* Having sorted the records, rearrange order of the pointers */
  for (iorder = 0; iorder < nec; iorder++)
    {
      /* Get the current enzyme record */
      enzyme_ptr = enzyme_index_ptr[iorder];

      /* Convert the sort number back into an E.C. number */
      if (type == EC_COUNT_SORT)
        {
          sort_number[0] = '\0';
          for (ipos = 0; ipos < strlen(enzyme_ptr->sort_number);
               ipos = ipos + 5)
            {
              /* If not first number, append a dot */
              if (ipos > 0)
                strcat(sort_number,".");

              /* Convert back into a number */
              ch = enzyme_ptr->sort_number[ipos];
              if (ch >= '0' && ch <= '9')
                {
                  strncpy(number_string,enzyme_ptr->sort_number + ipos,4);
                  number_string[4] = '\0';
                }
              else
                {
                  strncpy(number_string,
                          enzyme_ptr->sort_number + ipos + 1,3);
                  number_string[3] = '\0';
                  chstring[0] = ch;
                  chstring[1] = '\0';
                  strcat(sort_number,chstring);
                }
              inum = atoi(number_string);

              /* Append the number, or dash */
              if (inum > 0)
                sprintf(number_string,"%d",inum);
              else
                strcpy(number_string,"-");
              strcat(sort_number,number_string);
            }

          /* Store the new sort number */
          //OUT if (enzyme_ptr->sort_number != &null)
          //OUT free(enzyme_ptr->sort_number);
          store_string(&(enzyme_ptr->sort_number),sort_number);
        }

      /* If this is the first, then make the first record */
      if (iorder == 0)
        first_enzyme_ptr = enzyme_ptr;

      /* Otherwise, get the previous enzyme record to point to this one */
      else
        last_enzyme_ptr->next_enzyme_ptr = enzyme_ptr;

      /* Blank out current record's next enzyme pointer */
      enzyme_ptr->next_enzyme_ptr = NULL;

      /* Store current record as last encountered */
      last_enzyme_ptr = enzyme_ptr;
    }

  /* Free up the index array created here */
  free(enzyme_index_ptr);

  /* Return the first Enzyme record */
  *fst_enzyme_ptr = first_enzyme_ptr;
}
/***********************************************************************

get_enzyme_names  -  Pick up the enzyme names

***********************************************************************/

int get_enzyme_names(struct enzyme *first_enzyme_ptr,int nec,
                     struct c_file_data file_name_ptr)
{
  char ec_number[20], file_ec_number[20], sort_number[60];
  char file_name[FILENAME_LEN], input_line[LINELEN + 1];
  char enzyme_name[LINELEN + 1], last_line[LINELEN + 1];

  char *string_ptr;

  int len;
  int done, got_record;

  struct enzyme *enzyme_ptr, *next_enzyme_ptr;

  FILE *file_ptr;

  /* Initialise variables */
  ec_number[0] = file_ec_number[0] = enzyme_name[0] = '\0';
  sort_number[0] = '\0';
  done = FALSE;
  input_line[0] = last_line[0] = '\0';

  /* Form name of input file */
  strcpy(file_name,file_name_ptr.other_dir[ENZYME_REAL_DATA_DIR]);
  strcat(file_name,"/eheaders.txt");

  /* Open the input data file */
  if ((file_ptr = fopen(file_name,"r")) == NULL)
    {
      /* Otherwise, error must be due to missing file */
      //printf("*** Warning. Unable to open file: %s\n",file_name);
      done = TRUE;
    }

  /* Get the first enzyme record */
  enzyme_ptr = next_enzyme_ptr = first_enzyme_ptr;

  /* Loop until file has been processed */
  while (done == FALSE)
    {
      /* If we don't have an E.C. record, get the next one */
      if (ec_number[0] == '\0')
        {
          /* If have a next enzyme record, get its E.C. number */
          if (next_enzyme_ptr != NULL)
            {
              /* Get the next record */
              enzyme_ptr = next_enzyme_ptr;

              /* Check if E.C. number is a list */
              string_ptr = strchr(enzyme_ptr->ec_number,',');

              /* Get the E.C. number */
              if (string_ptr == NULL)
                strcpy(ec_number,enzyme_ptr->ec_number);
              else
                {
                  ec_number[0] = '\0';
                  store_string(&(enzyme_ptr->enzyme_name),
                               "Multiple E.C. codes");
                }

              /* Get pointer to the next enzyme */
              next_enzyme_ptr = enzyme_ptr->next_enzyme_ptr;
            }

          /* Otherwise, we are done */
          else
            done = TRUE;
        }

      /* If we don't have an enzyme record from the file, get one now */
      if (file_ec_number[0] == '\0')
        {
          /* Initialise flag */
          got_record = FALSE;

          /* Read until we have an E.C. number and associated
             description */
          while (got_record == FALSE &&
                 fgets(input_line,LINELEN,file_ptr) != NULL)
            {
              /* Truncate the string */
              string_chop(input_line);

              /* If this is the E.C. number, then store */
              if (!strcmp(last_line,"EC_NO"))
                {
                  /* Store the E.C. number */
                  strcpy(file_ec_number,input_line + 1);
                  sort_number[0] = '\0';
                }

              /* If this is the description, then store and set flag */
              else if (!strcmp(last_line,"DESC"))
                {
                  /* Save the enzyme name */
                  strcpy(enzyme_name,input_line + 1);

                  /* Set flag */
                  got_record = TRUE;
                }

              /* Save the current record */
              strcpy(last_line,input_line);
            }

          /* If end of file encountered, then we're done */
          if (got_record == FALSE)
            {
              /* Set flags */
              file_ec_number[0] = sort_number[0] = '\0';
              done = TRUE;
            }
        }

      /* If have E.C. numbers to compare, then do so */
      if (done == FALSE && ec_number[0] != '\0' &&
          file_ec_number[0] != '\0')
        {
          /* If both match, then store the enzyme name and reinitialise */
          if (!strcmp(ec_number,file_ec_number))
            {
              /* Save the enzyme name */
              store_string(&(enzyme_ptr->enzyme_name),enzyme_name);

              /* Re-initialise both E.C. numbers */
              ec_number[0] = file_ec_number[0] = sort_number[0] = '\0';
            }

          /* Otherwise, need to compare which record is higher */
          else
            {
              /* Convert the E.C. number from file unless already
                 converted */
              if (sort_number[0] == '\0')
                convert_ec(file_ec_number,sort_number);

              /* If enzyme from record too high, then need to read next
                 record from file */
              if (strcmp(enzyme_ptr->sort_number,sort_number) > 0)
                file_ec_number[0] = '\0';
 
              /* Otherwise, get next record */
              else
                ec_number[0] = '\0';
            }
        }
    }

  /* Close the input file */
  fclose(file_ptr);
}
/***********************************************************************

write_enzymes  -  Write out all the enzymes names

***********************************************************************/

void write_enzymes(struct enzyme *first_enzyme_ptr,int nec)
{
  int ienzymes;

  struct enzyme *enzyme_ptr;

  /* Initialise variables */
  ienzymes = 0;

  /* Start the Enzymes descriptions listing */
  fprintf(outfile_ptr,":EC_NAME\n");

  /* Get the first enzymes record */
  enzyme_ptr = first_enzyme_ptr;

  /* Loop over the enzymes records to write out */
  while (enzyme_ptr != NULL)
    {
      /* Increment enzymes count */
      ienzymes++;

      /* Write out the enzymes details */
      fprintf(outfile_ptr,"%s\t%s\t%d\t%d\n",enzyme_ptr->ec_number,
              enzyme_ptr->enzyme_name,enzyme_ptr->narch,
              enzyme_ptr->ntotal);

      /* Get the next enzymes record */
      enzyme_ptr = enzyme_ptr->next_enzyme_ptr;
    }
}
/***********************************************************************

calc_arch_dists  -  Write out the architecture distances matrix

***********************************************************************/

void calc_arch_dists(int type,struct archit *first_archit_ptr,int narchit,
                     int **dist_matrix,struct archit ***archit_ind_ptr,
                     char *parent_architecture,struct parameters *params)
{
  char pfam_id[MAX_COMPARE][PFAM_LEN];
  char p_architecture[LINELEN + 1], n1_architecture[LINELEN + 1];
  char n2_architecture[LINELEN + 1];

  int dist, iarch, ipoint, jarch, narray;
  int domain_no[2][MAX_COMPARE], ndomains[2], nunique[2];
  int verbose;

  int *array;

  struct archit *archit1_ptr, *archit2_ptr;

  struct archit **archit_index_ptr;

  /* Initialise variables */
  iarch = jarch = 0;
  verbose = FALSE;

  /* Save the parent architecture */
  strcpy(p_architecture,parent_architecture);

  /* Modify any split-CATH domains */
  if (params->cath_domains == TRUE)
    modify_split_domains(p_architecture);

  /* Create an array to store all the architecture-architecture distances */
  narray = narchit * narchit;

  /* Create a distances matrix for holding the architecture-architecture
     distances */
  array = (int *) malloc(narray * sizeof(int));
  if (array == NULL)
    {
      printf("*** ERROR. Unable to allocate memory for distances array\n");
      exit(1);
    }

  /* Initialise array */
  for (ipoint = 0; ipoint < narray; ipoint++)
    array[ipoint] = 0;

  /* Create an array for storing the architectures */
  archit_index_ptr
    = (struct archit **) malloc(narchit * sizeof(struct archit *));
 
  /* Get the first architecture record */
  archit1_ptr = first_archit_ptr;

  /* Loop over the architecture records */
  while (archit1_ptr != NULL)
    {
      /* Store pointer in index array */
      archit_index_ptr[iarch] = archit1_ptr;

      /* Save this architecture */
      strcpy(n1_architecture,archit1_ptr->architecture);

      /* Modify any split-CATH domains */
      if (params->cath_domains == TRUE)
        modify_split_domains(n1_architecture);

      /* Extract this architecture's Pfam domains */
      nunique[FIRST]
        = extract_domains(type,n1_architecture,pfam_id,
                          domain_no[FIRST],&(ndomains[FIRST]),0);

      /* Get the next architecture record */
      archit2_ptr = archit1_ptr->next_archit_ptr;
      jarch = iarch + 1;

      /* Loop over all other architecture records */
      while (archit2_ptr != NULL)
        {
          /* Save this architecture */
          strcpy(n2_architecture,archit2_ptr->architecture);

          /* Modify any split-CATH domains */
          if (params->cath_domains == TRUE)
            modify_split_domains(n2_architecture);

          /* Extract the domains from this architecture */
          nunique[SECOND]
            = extract_domains(type,n2_architecture,pfam_id,
                              domain_no[SECOND],&(ndomains[SECOND]),
                              nunique[FIRST]);

          /* Calculate the distance score between these two architectures */
          dist = score_architectures(type,n1_architecture,
                                     n2_architecture,ndomains,
                                     nunique,pfam_id,domain_no,
                                     p_architecture);
          
          /* Store the score in the distances array */
          ipoint = iarch + narchit * jarch;
          array[ipoint] = dist;

          /* Store the other way round, too */
          ipoint = jarch + narchit * iarch;
          array[ipoint] = dist;

          /* Get the next architecture record */
          archit2_ptr = archit2_ptr->next_archit_ptr;
          jarch++;
        }

      /* Get the next architecture record */
      archit1_ptr = archit1_ptr->next_archit_ptr;
      iarch++;
    }

  /* For debugging, list the architectures and then show the distances
     matrix */
  if (verbose == TRUE)
    {
      /* List the architecture nodes */
      printf("\n");
      printf("\n");
      printf("Architectures:-\n");
      printf("-------------\n");
      for (iarch = 0; iarch < narchit; iarch++)
        printf("%3d. %s Score = %d\n",iarch,
               archit_index_ptr[iarch]->architecture,
               archit_index_ptr[iarch]->score);

      /* Print out the distances array */
      printf("\n");
      printf("Distance matrix\n");
      printf("\n");
      printf("    ");
      for (iarch = 0; iarch < narchit; iarch++)
        printf("%3d",iarch);
      printf("\n");
      printf("    ");
      for (iarch = 0; iarch < narchit; iarch++)
        printf("---");
      printf("\n");

      /* Loop down the page */
      for (iarch = 0; iarch < narchit; iarch++)
        {
          /* Print current architecture */
          printf("%3d.",iarch);

          /* Loop across the page */
          for (jarch = 0; jarch < narchit; jarch++)
            {
              /* Get this point in the distance matrix */
              ipoint = iarch + narchit * jarch;

              /* Print this distance */
              printf("%3d",array[ipoint]);
            }
          printf("\n");
        }

      /* Close off print */
      printf("\n");
      printf("\n");
    }

  /* Return the distance matrix */
  *dist_matrix = array;

  /* Return the architectures index array */
  *archit_ind_ptr = archit_index_ptr;
}
/***********************************************************************

get_connectivities  -  Determine the graph connectivities 

***********************************************************************/

void get_connectivities(struct archit **archit_index_ptr,int narchit,
                        int *dist_matrix,int parent_id)
{
  int dist, iarch, iconnect, inode, ipoint, jarch, jnode, loop;
  int best_node, next_free;
  int mindist, minparent_dist, nconn, nunconn, ntoconnect;

  int *connected, *unconnected, *toconnect;

  struct archit *archit_ptr;

  /* Initialise variables */
  loop = 0;
  nconn = nunconn = 0;

  /* Create arrays holding the connected and unconnected architectures */
  connected = (int *) malloc(narchit * sizeof(int));
  toconnect = (int *) malloc(narchit * narchit * sizeof(int));
  unconnected = (int *) malloc(narchit * sizeof(int));
  if (connected == NULL || toconnect == NULL || unconnected == NULL)
    {
      printf("*** ERROR. Unable to allocate memory for "
             "connected/unconnected arrays\n");
      exit(1);
    }

  /* Initialise arrays */
  for (iarch = 0; iarch < narchit; iarch++)
    connected[iarch] = unconnected[iarch] = -1;
  for (iarch = 0; iarch < narchit * narchit; iarch++)
    toconnect[iarch] = -1;

  /* Place the parent architecture into the connected list */
  connected[nconn] = parent_id;
  nconn++;

  /* Place all other architectures into the unconnected list */
  for (iarch = 0; iarch < narchit; iarch++)
    {
      if (iarch != parent_id)
        {
          unconnected[nunconn] = iarch;
          nunconn++;
        }
    }

  /* Write out the parent id */
  fprintf(outfile_ptr,":ID_PARENT\n");
  fprintf(outfile_ptr,"%d\n",parent_id);

  /* Write out heading record */
  fprintf(outfile_ptr,":CONNECTIONS\n");

  /* Loop while there are entries in the unconnected list */
  while (nunconn > 0)
    {
      /* Increment loop counter */
      loop++;

      /* Initialise minimum distance */
      mindist = 10000;
      ntoconnect = 0;

      /* Loop over the connected nodes to find the minimum distance
         between connected and unconnected nodes */
      for (inode = 0; inode < nconn; inode++)
        {
          /* Get the corresponding architecture identifier */
          iarch = connected[inode];

          /* Loop over the unconnected nodes */
          for (jnode = 0; jnode < nunconn; jnode++)
            {
              /* Get the corresponding architecture identifier */
              jarch = unconnected[jnode];

              /* Get the corresponding cell in the distance matrix */
              ipoint = iarch + narchit * jarch;

              /* If this is the lowest distance so far, then store */
              if (dist_matrix[ipoint] < mindist)
                {
                  /* Initialise the to-connect list */
                  ntoconnect = 0;

                  /* Store the minimum distance and add this node to
                     the to-connect list */
                  mindist = dist_matrix[ipoint];
                  toconnect[ntoconnect] = jnode;

                  /* Increment count of nodes to connect */
                  ntoconnect++;
                }

              /* If distance matches the lowest distance so far, then
                 add to list of nodes to connect */
              else if (dist_matrix[ipoint] == mindist)
                {
                  /* Add this node to the to-connect list */
                  toconnect[ntoconnect] = jnode;

                  /* Increment count of nodes to connect */
                  ntoconnect++;
                }
            }
        }

      /* Loop over all the nodes to be connected to connect them to the
         best already-connected nodes */
      for (iconnect = 0; iconnect < ntoconnect; iconnect++)
        {
          /* Get this node */
          jnode = toconnect[iconnect];

          /* Get the corresponding architecture id */
          jarch = unconnected[jnode];

          /* Check that this node hasn't already been processed */
          if (jarch > -1)
            {
              /* Re-initialise minimum distance */
              minparent_dist = 10000;
              best_node = 0;

              /* Loop over all the connected nodes to determine which one
                 to attach this node to */
              for (inode = 0; inode < nconn; inode++)
                {
                  /* Get the corresponding architecture identifier */
                  iarch = connected[inode];

                  /* Get the corresponding cell in the distance matrix */
                  ipoint = iarch + narchit * jarch;

                  /* If distance matches the minimum distance, then see
                     whether to connect to this node */
                  if (dist_matrix[ipoint] == mindist)
                    {
                      /* Get its architecture record and its distance from the
                         parent architecture */
                      archit_ptr = archit_index_ptr[iarch];
                      dist = archit_ptr->score;

                      /* If this is the minimum distance so far, store as the
                         best node to attach to */
                      if (dist < minparent_dist)
                        {
                          minparent_dist = dist;
                          best_node = inode;
                        }
                    }
                }

              /* Get the best connection */
              iarch = connected[best_node];

              /* Write out the connection */
              fprintf(outfile_ptr,"%d\t%d\t%d\n",iarch,jarch,mindist);

              /* Add current node to the connected list */
              connected[nconn] = jarch;
              nconn++;

              /* Delete the current node from the unconnected list */
              unconnected[jnode] = -1;
            }
        }

      /* Initialise location of next free position in the unconnected
         list */
      next_free = 0;

      /* Loop over the unconnected nodes to remove all the deleted ones */
      for (jnode = 0; jnode < nunconn; jnode++)
        {
          /* If not deleted, then move to next free position */
          if (unconnected[jnode] > -1)
            {
              /* Shift to free space */
              unconnected[next_free] = unconnected[jnode];

              /* Shift free space by one */
              next_free++;
            }
        }

      /* Save number of nodes now in the unconnected list */
      nunconn = next_free;
    }
}
/***********************************************************************

write_extra_data  -   Write out additional data

***********************************************************************/

void write_extra_data(int max_architectures)
{
  /* Write out the Maximum architecture */
  fprintf(outfile_ptr,":MAX_ARCH\n");
  fprintf(outfile_ptr,"%d\n",max_architectures);
}
/***********************************************************************

pfam_search  -  Search for all the given Pfam ids

***********************************************************************/

int pfam_search(int type,char *archschema_dir,
                struct c_file_data file_name_ptr,
                struct parameters *params)
{
  char dom_type[6], err_message[LINELEN];

  int cut_off, narchit, nec, nhits, npfam, nseqs, nspecies;
  int ntotal, ncollapsed, nfiltered, parent_id;

  int *dist_matrix;

  struct archit *first_archit_ptr, *parent_archit_ptr;
  struct enzyme *first_enzyme_ptr;
  struct hit *first_hit_ptr, *parent_hit_ptr;
  struct pfam *first_pfam_ptr;


  struct species *first_species_ptr, *last_species_ptr;

  struct archit **archit_index_ptr;
  struct hit **hit_index_ptr;

  /* Initialise variables */
  cut_off = -1;
  if (params->cath_domains == FALSE)
    strcpy(dom_type,"Pfam");
  else
    strcpy(dom_type,"CATH");
  nec = nhits = nspecies = 0;
  ntotal = ncollapsed = nfiltered = 0;

  /* Initialise pointers */
  first_archit_ptr = NULL;
  first_enzyme_ptr = NULL;
  first_hit_ptr = parent_hit_ptr = NULL;
  first_pfam_ptr = NULL;
  first_species_ptr = NULL;

  /* Read in the archpfam.dat offsets from the archpfam.idx file for all
     the Pfam domains to be searched */
  nhits = read_archpfam_idx(type,archschema_dir,params);

  /* If no hits found, then exit */
  if (nhits == 0)
    {
      if (params->first_pfam_ptr == NULL)
        write_error_message(NO_HITS,"No hits found",params->output_type);
      else
        {
          sprintf(err_message,"%s code %s not found in database!",
                  dom_type,params->first_pfam_ptr->pfam_id);
          write_error_message(CODE_NOT_FOUND,err_message,
                              params->output_type);
        }
    }

  /* Pick up the offsets from the archpfam.dat file */
  nhits = read_archpfam_dat(type,archschema_dir,params);

  /* If have any offsets, pick up the data from the
     archstruc.uniprot.joined.by.accession file */
  nhits = get_architectures(type,&first_hit_ptr,&parent_hit_ptr,
                            &first_species_ptr,&nspecies,
                            archschema_dir,params);

  /* If no hits returned, write error message */
  if (nhits == 0)
    write_error_message(NO_HITS,"No hits found",params->output_type);

  /* Sort the species and filter out the duplicate records */
  if (nspecies > 1)
    nspecies = sort_species_records(&first_species_ptr,nspecies,params);

  /* Pick up the species names */
  get_species_names(first_species_ptr,archschema_dir,params);

  /* If initial search was by Pfam domain, rather than by UniProt code,
     then need to identify one sequence as the "parent sequence" */
  if (parent_hit_ptr == NULL)
    parent_hit_ptr = choose_parent(type,first_hit_ptr,&nhits,params);

  /* If parent sequence not found, write error message */
  if (parent_hit_ptr == NULL)
    write_error_message(NO_PARENT,"Indexing error: parent sequence "
                        "not found",params->output_type);

  /* Sort the hit records by architecture */
  sort_hit_records(first_hit_ptr,&hit_index_ptr,nhits,params);

  /* Create an architecture record for each unique architecture */
  narchit
    = extract_architectures(hit_index_ptr,nhits,parent_hit_ptr,
                            &first_archit_ptr,params);
  ntotal = ncollapsed = nfiltered = narchit;

  /* Identify the parent architecture */
  parent_archit_ptr
    = find_architecture(first_archit_ptr,parent_hit_ptr->architecture);

  /* If parent architecture not found, write error message */
  if (parent_archit_ptr == NULL)
    write_error_message(NO_PARENT,"Indexing error: parent architecture "
                        "not found",params->output_type);

  /* If have too many architectures, then replace Pfam-B domains by a
     single catch-all domain */
  if (narchit > params->max_architectures &&
      params->cath_domains == FALSE &&
      params->output_type == ARCHSCHEMA_OUTPUT)
    {
      /* Collapse all Pfam-B domains to PB000000 */
      narchit
        = collapse_pfamb_domains(type,&first_archit_ptr,
                                 parent_hit_ptr->architecture,params);
      ncollapsed = nfiltered = narchit;
    }

  /* Score all architecture distances from the parent architecture to
     calculate a cut-off score for which architectures to keep */
  cut_off
    = compare_to_parent(type,first_archit_ptr,
                        parent_archit_ptr->architecture,params);

  /* If we still have too many architectures, then need to score them
     so that we only keep those that are the closest to the parent
     architecture */
  if (narchit > params->max_architectures &&
      params->output_type == ARCHSCHEMA_OUTPUT)
    {
      /* Delete any unwanted architectures */
      if (cut_off > -1)
        narchit = delete_architectures(&first_archit_ptr,cut_off);
      nfiltered = narchit;
    }

  /* If replacing enzyme data with SSG data, then read in the assignments */
  if (params->get_ssg == TRUE)
    nec = read_ssg_data(first_archit_ptr,file_name_ptr);

  /* Get the enzyme or SSG data */
  nec = get_enzyme_data(first_archit_ptr,&first_enzyme_ptr,params);

  /* Prune sequences for any architectures that have a huge number */
  if (params->output_type == ARCHSCHEMA_OUTPUT)
    nseqs = prune_sequences(first_archit_ptr,parent_hit_ptr,params);

  /* Extract the Pfam domains from the remaining architectures */
  npfam
    = extract_pfam_domains(type,first_archit_ptr,&first_pfam_ptr,params);

  /* Sort the pfam ids and filter out the duplicate records */
  sort_pfam_records(&first_pfam_ptr,npfam,params);

  /* Delete all duplicate Pfam records */
  npfam = delete_pfam_records(&first_pfam_ptr);

  /* Pick up the pfam names */
  get_pfam_names(type,first_pfam_ptr,archschema_dir,params);

  /* Assign the colours to the Pfam domains */
  assign_pfam_colours(first_pfam_ptr,params);

  /* Write out the hit summary */
  if (params->output_type == ARCHSCHEMA_OUTPUT)
    write_hit_summary(nhits,ntotal,ncollapsed,nfiltered,parent_hit_ptr,
                      npfam,nspecies,params);
  else
    write_dat_summary(nhits,ntotal,ncollapsed,nfiltered,parent_hit_ptr,
                      npfam,nspecies,first_pfam_ptr,params);

  /* If writing out to .dat file, sort the architectures in order of
     similarity to parent */
  if (params->output_type == DAT_OUTPUT)
    sort_by_score(&first_archit_ptr,narchit,params);

  /* Write out the hits */
  parent_id
    = write_hits(first_archit_ptr,parent_archit_ptr->architecture,
                 first_pfam_ptr,params);

  /* Write out the species list */
  if (params->output_type == ARCHSCHEMA_OUTPUT)
    write_species(first_species_ptr,nspecies);

  /* Write out the Pfam ids and names */
  write_pfam(first_pfam_ptr,params);

  /* Write out the remaining data tables */
  if (params->output_type == ARCHSCHEMA_OUTPUT)
    {
      /* If not using enzyme record for storing SSG data, get the enzyme
         names from the eheaders.txt file */
      if (params->get_ssg == FALSE)
        {
          /* Sort the enzyme records by E.C. code */
          sort_enzyme_records(&first_enzyme_ptr,nec,EC_CODE_SORT,params);

          /* Get the enzyme names from the eheaders.txt file */
          get_enzyme_names(first_enzyme_ptr,nec,file_name_ptr);
        }

      /* Sort the enzyme records by architecture and sequence counts */
      sort_enzyme_records(&first_enzyme_ptr,nec,EC_COUNT_SORT,params);

      /* Write out the list of enzyme names and their counts */
      write_enzymes(first_enzyme_ptr,nec);

      /* Calculate the architecture distances matrix */
      calc_arch_dists(type,first_archit_ptr,narchit,&dist_matrix,
                      &archit_index_ptr,parent_archit_ptr->architecture,
                      params);

      /* Determine the graph connectivities */
      get_connectivities(archit_index_ptr,narchit,dist_matrix,parent_id);

      /* Write out additional data */
      write_extra_data(params->max_architectures);
    }

  /* Return number of hits */
  return(nhits);
}
/***********************************************************************

create_filterseq_record  -  Create and initialise a new filterseq record

***********************************************************************/

struct filterseq *create_filterseq_record(struct filterseq **fst_filterseq_ptr,
                                      struct filterseq **lst_filterseq_ptr,
                                      char *seq_id,
                                      struct parameters *params)
{
  struct filterseq *filterseq_ptr, *first_filterseq_ptr, *last_filterseq_ptr;

  /* Initialise filterseq pointers */
  filterseq_ptr = NULL;
  first_filterseq_ptr = *fst_filterseq_ptr;
  last_filterseq_ptr = *lst_filterseq_ptr;

  /* Allocate memory for structure to hold filterseq info */
  filterseq_ptr = (struct filterseq *) malloc(sizeof(struct filterseq));
  if (filterseq_ptr == NULL)
    {
      write_error_message(MEMORY_ALLOCATION,
                          "Cannot allocate memory for struct filterseq",
                          params->output_type);
    }

  /* If this is the very first, save its pointer */
  if (first_filterseq_ptr == NULL)
    first_filterseq_ptr = filterseq_ptr;

  /* Add link from previous filterseq to the current one */
  if (last_filterseq_ptr != NULL)
    last_filterseq_ptr->next_filterseq_ptr = filterseq_ptr;
  last_filterseq_ptr = filterseq_ptr;

  /* Store known elements of the structure */
  store_string(&(filterseq_ptr->seq_id),seq_id);

  /* Initialise pointer to next filterseq record */
  filterseq_ptr->next_filterseq_ptr = NULL;

  /* Return current pointers */
  *fst_filterseq_ptr = first_filterseq_ptr;
  *lst_filterseq_ptr = last_filterseq_ptr;

  /* Return new filterseq pointer */
  return(filterseq_ptr);
}
/***********************************************************************

get_filter_list  -  Read in the list of filter sequences to be included

***********************************************************************/

int get_filter_list(struct parameters *params)
{
  char file_name[FILENAME_LEN], input_line[LINELEN + 1];
  char error_message[FILENAME_LEN + 50];

  int len, nseqs;

  struct filterseq *filterseq_ptr;
  struct filterseq *first_filterseq_ptr, *last_filterseq_ptr;

  FILE *file_ptr;

  /* Initialise variables */
  filterseq_ptr = first_filterseq_ptr = last_filterseq_ptr = NULL;
  nseqs = 0;

  /* Open the input data file */
  if ((file_ptr = fopen(params->filter_seqs_file,"r")) != NULL)
    {
      /* Read in the list of sequences to be retained */
      while (fgets(input_line,LINELEN,file_ptr) != NULL)
        {
          /* Truncate the string */
          len = string_truncate(input_line,LINELEN);

          /* Create a record to store this sequence */
          filterseq_ptr
            = create_filterseq_record(&first_filterseq_ptr,
                                      &last_filterseq_ptr,
                                      input_line,params);
          
          /* Increment count of filter sequences */
          nseqs++;
        }

      /* Close the input file */
      fclose(file_ptr);

    }

  /* If file not found, show error message */
  else
    {
      strcpy(error_message,"*** Warning. File not found: ");
      strcat(error_message,params->filter_seqs_file);
      write_error_message(FILE_NOT_FOUND,error_message,
                          params->output_type);
    }

  /* Save the first filter sequence pointer in the parameters */
  params->first_filterseq_ptr = first_filterseq_ptr;

  /* Return number of sequences stored */
  return(nseqs);
}
/***********************************************************************

perform_search  -  Perform the required search

***********************************************************************/

int perform_search(int type,struct c_file_data file_name_ptr,
                   struct parameters *params)
{
  char archschema_dir[FILENAME_LEN];

  int narchit, nhits;

  /* Determine which Archschema directory the searches are to be
     performed in */
  if (params->cath_domains == FALSE)
    strcpy(archschema_dir,
           file_name_ptr.other_dir[ARCHSCHEMA_REAL_DATA_DIR]);
  else
    {
      if (params->old_cath == FALSE)
        strcpy(archschema_dir,
               file_name_ptr.other_dir[ARCHCATH_REAL_DATA_DIR]);
      else
        strcpy(archschema_dir,
               file_name_ptr.other_dir[ARCHCATH_OLD_REAL_DATA_DIR]);
    }

  /* If have a UniProt code, find its architecture record */
  if (params->uniprot_acc != &null || params->uniprot_id != &null)
    {
      /* Perform the search */
      nhits = uniprot_search(type,archschema_dir,params);
    }

  /* If searching on CATH domains, add split domains to the mix */
  if (params->cath_domains == TRUE)
    add_split_domains(params);

  /* Sort the Pfam domains if we have more than one */
  if (params->npfam > 1)
    sort_pfam_records(&(params->first_pfam_ptr),params->npfam,params);

  /* Search for all the given Pfam ids */
  nhits = pfam_search(type,archschema_dir,file_name_ptr,params);

  /* Return number of hits */
  return(nhits);
}
/***********************************************************************

test_alignment  -  Test alignment

***********************************************************************/

void test_alignment(int type,char *architecture1,char *architecture2,
                    char *parent_architecture)
{
  char pfam_id[MAX_COMPARE][PFAM_LEN];

  int dist, iarch;
  int domain_no[2][MAX_COMPARE], ndomains[2], nunique[2];

  /* Extract this architecture's Pfam domains */
  nunique[FIRST]
    = extract_domains(type,architecture1,pfam_id,domain_no[FIRST],
                      &(ndomains[FIRST]),0);

  /* Extract the domains from this architecture */
  nunique[SECOND]
    = extract_domains(type,architecture2,pfam_id,domain_no[SECOND],
                      &(ndomains[SECOND]),nunique[FIRST]);

  /* Calculate the distance score between these two architectures */
  dist = score_architectures(type,architecture1,architecture2,ndomains,
                             nunique,pfam_id,domain_no,
                             parent_architecture);
}
/***********************************************************************

                             M   A   I   N

***********************************************************************/

int main(int argc,char *argv[])
{
  int ikey, narch, nhits, npfam, nrecords[NKEYS], type;

  struct c_file_data file_name_ptr;
  struct parameters params;

  struct arch *first_arch_ptr[NKEYS];
  struct pfam *first_pfam_ptr;

  struct arch **arch_index_ptr[NKEYS];

  /* Initialise variables */
  outfile_ptr = NULL;

  /* Interpret the command-line parameters */
  get_command_arguments(argv,argc - 1,&params);

  /* Determine which domain type we are dealing with */
  if (params.cath_domains == FALSE)
    type = PFAM;
  else
    type = CATH;

  /* Read in all the directory paths from the parameter file, CATHPARAM */
  c_get_paths(&file_name_ptr,params.run_64);

  /* Open the output file, or write to stdout */
  if (params.out_name[0] == '\0')
    outfile_ptr = stdout;
  else
    {
      if ((outfile_ptr = fopen(params.out_name,"w")) ==  NULL)
        {
          printf("\n*** Unable to open output file [%s]\n",
                 params.out_name);
          exit(1);
        }
    }

  /* If compiling indexes, process the input files */
  if (params.run_type == COMPILE_INDEXES)
    {
      /* Read in the archstruc.uniprot.joined.by.accession file */
      get_arch_data(type,first_arch_ptr,nrecords,&params);

      /* Sort the records just created */
      for (ikey = 0; ikey < NKEYS; ikey++)
        {
          if (nrecords[ikey] > 0)
            {
              sort_arch_records(first_arch_ptr[ikey],&(arch_index_ptr[ikey]),
                                nrecords[ikey],&params);
            }
        }

      /* Write out the UniProt index files */
      for (ikey = KEY_UNIPROT_ACC; ikey < NKEYS; ikey++)
        {
         /* Write the list of UniProt codes and their offsets */
          write_offsets_file(type,ikey,arch_index_ptr[ikey],nrecords[ikey]);

          /* Create an index to the offsets file */
          write_offsets_index(type,ikey,&params);
        }

      /* Transfer the Pfam ids and offsets to a new set of Pfam records */
      npfam = transfer_to_pfam(&first_pfam_ptr,arch_index_ptr[KEY_PFAM_ID],
                               nrecords[KEY_PFAM_ID],&params);

      /* Read in the offsets for the pfam_domains.parsed.sorted file */
      // OUT  get_offsets(first_pfam_ptr,"pfam_domains.parsed.sorted",DOMAINS);

      /* Repeat for the domseqs.struc file */
      // OUT  get_offsets(first_pfam_ptr,"domseqs.struc",DOMSEQS);

      /* Write out the offsets index */
      if (type == PFAM)
        {
          /* Write out the index and the index to the index */
          write_pfam_index("archpfam.dat",first_pfam_ptr);
          index_to_index("archpfam.dat","archpfam.idx");
        }
      else
        {
          /* Write out the index and the index to the index */
          write_pfam_index("archcath.dat",first_pfam_ptr);
          index_to_index("archcath.dat","archcath.idx");
        }
    }

  /* Otherwise, perform the required search */
  else
    {
      /* If we have a file listing the sequences to be included, read
         in the sequences */
      if (params.filter_seqs_file != &null)
        get_filter_list(&params);

      /* Write XML header */
      /*
      fprintf(outfile_ptr,"<?xml version=\"1.0\" encoding=\"UTF-8\" ?>\n");
      fprintf(outfile_ptr,"<java>\n");
      fprintf(outfile_ptr,"<object class=\"java.util.Vector\">\n");
      */

      /* Perform search */
      nhits = perform_search(type,file_name_ptr,&params);

      /* Close the XML file */
      /*
      fprintf(outfile_ptr,"</object>\n");
      fprintf(outfile_ptr,"</java>\n");
      */
      fclose(outfile_ptr);
    }

  return(0);
}
