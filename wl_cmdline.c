/*
  File autogenerated by gengetopt version 2.22.5
  generated with the following command:
  gengetopt --file-name=wl_cmdline --unamed-opts

  The developers of gengetopt consider the fixed text that goes in all
  gengetopt output files to be in the public domain:
  we make no copyright claims on it.
*/

/* If we use autoconf.  */
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifndef FIX_UNUSED
#define FIX_UNUSED(X) (void) (X) /* avoid warnings for unused params */
#endif

#include <getopt.h>

#include "wl_cmdline.h"

const char *gengetopt_args_info_purpose = "Sample the Density of States by a Wang-Landau MC simulation";

const char *gengetopt_args_info_usage = "Usage: " CMDLINE_PARSER_PACKAGE " [OPTIONS]... [FILES]...";

const char *gengetopt_args_info_description = "";

const char *gengetopt_args_info_help[] = {
  "  -h, --help                 Print help and exit",
  "  -V, --version              Print version and exit",
  "\nGeneral options:",
  "  -b, --bins=INT             Number of (equidistant) histogram bins  \n                               (default=`100')",
  "  -c, --checksteps=LONGLONG  Number of Wang-Landau steps before histogram is \n                               checked for flatness  (default=`1000000')",
  "      --elow=DOUBLE          Lower limit of sampling window (currently n/a)",
  "      --ehigh=DOUBLE         Upper limit of sampling window (currently n/a)",
  "      --flat=FLOAT           Flatness criterion for the histogram",
  "      --info                 Show settings  (default=off)",
  "  -m, --max=DOUBLE           Upper energy bound for sampling",
  "  -f, --mod=DOUBLE           Final value of Wang-Landau modification factor",
  "  -n, --norm=INT             Number of bins used for normalization",
  "  -r, --resolution=DOUBLE    Sampling resolution (histogram bin width)  \n                               (default=`0.5')",
  "  -l, --steplimit=LONGLONG   Maximum number of MC steps to perform  \n                               (default=`100000000')",
  "  -S, --seed=LONG            Seed for random number generation",
  "  -T, --Temp=FLOAT           Simulation temperature in Celsius (currently n/a)",
  "  -v, --verbose              Verbose output  (default=off)",
  "  -d, --debug                Debugging output  (default=off)",
    0
};

typedef enum {ARG_NO
  , ARG_FLAG
  , ARG_INT
  , ARG_LONG
  , ARG_FLOAT
  , ARG_DOUBLE
  , ARG_LONGLONG
} cmdline_parser_arg_type;

static
void clear_given (struct gengetopt_args_info *args_info);
static
void clear_args (struct gengetopt_args_info *args_info);

static int
cmdline_parser_internal (int argc, char **argv, struct gengetopt_args_info *args_info,
                        struct cmdline_parser_params *params, const char *additional_error);


static char *
gengetopt_strdup (const char *s);

static
void clear_given (struct gengetopt_args_info *args_info)
{
  args_info->help_given = 0 ;
  args_info->version_given = 0 ;
  args_info->bins_given = 0 ;
  args_info->checksteps_given = 0 ;
  args_info->elow_given = 0 ;
  args_info->ehigh_given = 0 ;
  args_info->flat_given = 0 ;
  args_info->info_given = 0 ;
  args_info->max_given = 0 ;
  args_info->mod_given = 0 ;
  args_info->norm_given = 0 ;
  args_info->resolution_given = 0 ;
  args_info->steplimit_given = 0 ;
  args_info->seed_given = 0 ;
  args_info->Temp_given = 0 ;
  args_info->verbose_given = 0 ;
  args_info->debug_given = 0 ;
}

static
void clear_args (struct gengetopt_args_info *args_info)
{
  FIX_UNUSED (args_info);
  args_info->bins_arg = 100;
  args_info->bins_orig = NULL;
  args_info->checksteps_arg = 1000000;
  args_info->checksteps_orig = NULL;
  args_info->elow_orig = NULL;
  args_info->ehigh_orig = NULL;
  args_info->flat_orig = NULL;
  args_info->info_flag = 0;
  args_info->max_orig = NULL;
  args_info->mod_orig = NULL;
  args_info->norm_orig = NULL;
  args_info->resolution_arg = 0.5;
  args_info->resolution_orig = NULL;
  args_info->steplimit_arg = 100000000;
  args_info->steplimit_orig = NULL;
  args_info->seed_orig = NULL;
  args_info->Temp_orig = NULL;
  args_info->verbose_flag = 0;
  args_info->debug_flag = 0;
  
}

static
void init_args_info(struct gengetopt_args_info *args_info)
{


  args_info->help_help = gengetopt_args_info_help[0] ;
  args_info->version_help = gengetopt_args_info_help[1] ;
  args_info->bins_help = gengetopt_args_info_help[3] ;
  args_info->checksteps_help = gengetopt_args_info_help[4] ;
  args_info->elow_help = gengetopt_args_info_help[5] ;
  args_info->ehigh_help = gengetopt_args_info_help[6] ;
  args_info->flat_help = gengetopt_args_info_help[7] ;
  args_info->info_help = gengetopt_args_info_help[8] ;
  args_info->max_help = gengetopt_args_info_help[9] ;
  args_info->mod_help = gengetopt_args_info_help[10] ;
  args_info->norm_help = gengetopt_args_info_help[11] ;
  args_info->resolution_help = gengetopt_args_info_help[12] ;
  args_info->steplimit_help = gengetopt_args_info_help[13] ;
  args_info->seed_help = gengetopt_args_info_help[14] ;
  args_info->Temp_help = gengetopt_args_info_help[15] ;
  args_info->verbose_help = gengetopt_args_info_help[16] ;
  args_info->debug_help = gengetopt_args_info_help[17] ;
  
}

void
cmdline_parser_print_version (void)
{
  printf ("%s %s\n",
     (strlen(CMDLINE_PARSER_PACKAGE_NAME) ? CMDLINE_PARSER_PACKAGE_NAME : CMDLINE_PARSER_PACKAGE),
     CMDLINE_PARSER_VERSION);
}

static void print_help_common(void) {
  cmdline_parser_print_version ();

  if (strlen(gengetopt_args_info_purpose) > 0)
    printf("\n%s\n", gengetopt_args_info_purpose);

  if (strlen(gengetopt_args_info_usage) > 0)
    printf("\n%s\n", gengetopt_args_info_usage);

  printf("\n");

  if (strlen(gengetopt_args_info_description) > 0)
    printf("%s\n\n", gengetopt_args_info_description);
}

void
cmdline_parser_print_help (void)
{
  int i = 0;
  print_help_common();
  while (gengetopt_args_info_help[i])
    printf("%s\n", gengetopt_args_info_help[i++]);
}

void
cmdline_parser_init (struct gengetopt_args_info *args_info)
{
  clear_given (args_info);
  clear_args (args_info);
  init_args_info (args_info);

  args_info->inputs = 0;
  args_info->inputs_num = 0;
}

void
cmdline_parser_params_init(struct cmdline_parser_params *params)
{
  if (params)
    { 
      params->override = 0;
      params->initialize = 1;
      params->check_required = 1;
      params->check_ambiguity = 0;
      params->print_errors = 1;
    }
}

struct cmdline_parser_params *
cmdline_parser_params_create(void)
{
  struct cmdline_parser_params *params = 
    (struct cmdline_parser_params *)malloc(sizeof(struct cmdline_parser_params));
  cmdline_parser_params_init(params);  
  return params;
}

static void
free_string_field (char **s)
{
  if (*s)
    {
      free (*s);
      *s = 0;
    }
}


static void
cmdline_parser_release (struct gengetopt_args_info *args_info)
{
  unsigned int i;
  free_string_field (&(args_info->bins_orig));
  free_string_field (&(args_info->checksteps_orig));
  free_string_field (&(args_info->elow_orig));
  free_string_field (&(args_info->ehigh_orig));
  free_string_field (&(args_info->flat_orig));
  free_string_field (&(args_info->max_orig));
  free_string_field (&(args_info->mod_orig));
  free_string_field (&(args_info->norm_orig));
  free_string_field (&(args_info->resolution_orig));
  free_string_field (&(args_info->steplimit_orig));
  free_string_field (&(args_info->seed_orig));
  free_string_field (&(args_info->Temp_orig));
  
  
  for (i = 0; i < args_info->inputs_num; ++i)
    free (args_info->inputs [i]);

  if (args_info->inputs_num)
    free (args_info->inputs);

  clear_given (args_info);
}


static void
write_into_file(FILE *outfile, const char *opt, const char *arg, const char *values[])
{
  FIX_UNUSED (values);
  if (arg) {
    fprintf(outfile, "%s=\"%s\"\n", opt, arg);
  } else {
    fprintf(outfile, "%s\n", opt);
  }
}


int
cmdline_parser_dump(FILE *outfile, struct gengetopt_args_info *args_info)
{
  int i = 0;

  if (!outfile)
    {
      fprintf (stderr, "%s: cannot dump options to stream\n", CMDLINE_PARSER_PACKAGE);
      return EXIT_FAILURE;
    }

  if (args_info->help_given)
    write_into_file(outfile, "help", 0, 0 );
  if (args_info->version_given)
    write_into_file(outfile, "version", 0, 0 );
  if (args_info->bins_given)
    write_into_file(outfile, "bins", args_info->bins_orig, 0);
  if (args_info->checksteps_given)
    write_into_file(outfile, "checksteps", args_info->checksteps_orig, 0);
  if (args_info->elow_given)
    write_into_file(outfile, "elow", args_info->elow_orig, 0);
  if (args_info->ehigh_given)
    write_into_file(outfile, "ehigh", args_info->ehigh_orig, 0);
  if (args_info->flat_given)
    write_into_file(outfile, "flat", args_info->flat_orig, 0);
  if (args_info->info_given)
    write_into_file(outfile, "info", 0, 0 );
  if (args_info->max_given)
    write_into_file(outfile, "max", args_info->max_orig, 0);
  if (args_info->mod_given)
    write_into_file(outfile, "mod", args_info->mod_orig, 0);
  if (args_info->norm_given)
    write_into_file(outfile, "norm", args_info->norm_orig, 0);
  if (args_info->resolution_given)
    write_into_file(outfile, "resolution", args_info->resolution_orig, 0);
  if (args_info->steplimit_given)
    write_into_file(outfile, "steplimit", args_info->steplimit_orig, 0);
  if (args_info->seed_given)
    write_into_file(outfile, "seed", args_info->seed_orig, 0);
  if (args_info->Temp_given)
    write_into_file(outfile, "Temp", args_info->Temp_orig, 0);
  if (args_info->verbose_given)
    write_into_file(outfile, "verbose", 0, 0 );
  if (args_info->debug_given)
    write_into_file(outfile, "debug", 0, 0 );
  

  i = EXIT_SUCCESS;
  return i;
}

int
cmdline_parser_file_save(const char *filename, struct gengetopt_args_info *args_info)
{
  FILE *outfile;
  int i = 0;

  outfile = fopen(filename, "w");

  if (!outfile)
    {
      fprintf (stderr, "%s: cannot open file for writing: %s\n", CMDLINE_PARSER_PACKAGE, filename);
      return EXIT_FAILURE;
    }

  i = cmdline_parser_dump(outfile, args_info);
  fclose (outfile);

  return i;
}

void
cmdline_parser_free (struct gengetopt_args_info *args_info)
{
  cmdline_parser_release (args_info);
}

/** @brief replacement of strdup, which is not standard */
char *
gengetopt_strdup (const char *s)
{
  char *result = 0;
  if (!s)
    return result;

  result = (char*)malloc(strlen(s) + 1);
  if (result == (char*)0)
    return (char*)0;
  strcpy(result, s);
  return result;
}

int
cmdline_parser (int argc, char **argv, struct gengetopt_args_info *args_info)
{
  return cmdline_parser2 (argc, argv, args_info, 0, 1, 1);
}

int
cmdline_parser_ext (int argc, char **argv, struct gengetopt_args_info *args_info,
                   struct cmdline_parser_params *params)
{
  int result;
  result = cmdline_parser_internal (argc, argv, args_info, params, 0);

  if (result == EXIT_FAILURE)
    {
      cmdline_parser_free (args_info);
      exit (EXIT_FAILURE);
    }
  
  return result;
}

int
cmdline_parser2 (int argc, char **argv, struct gengetopt_args_info *args_info, int override, int initialize, int check_required)
{
  int result;
  struct cmdline_parser_params params;
  
  params.override = override;
  params.initialize = initialize;
  params.check_required = check_required;
  params.check_ambiguity = 0;
  params.print_errors = 1;

  result = cmdline_parser_internal (argc, argv, args_info, &params, 0);

  if (result == EXIT_FAILURE)
    {
      cmdline_parser_free (args_info);
      exit (EXIT_FAILURE);
    }
  
  return result;
}

int
cmdline_parser_required (struct gengetopt_args_info *args_info, const char *prog_name)
{
  FIX_UNUSED (args_info);
  FIX_UNUSED (prog_name);
  return EXIT_SUCCESS;
}


static char *package_name = 0;

/**
 * @brief updates an option
 * @param field the generic pointer to the field to update
 * @param orig_field the pointer to the orig field
 * @param field_given the pointer to the number of occurrence of this option
 * @param prev_given the pointer to the number of occurrence already seen
 * @param value the argument for this option (if null no arg was specified)
 * @param possible_values the possible values for this option (if specified)
 * @param default_value the default value (in case the option only accepts fixed values)
 * @param arg_type the type of this option
 * @param check_ambiguity @see cmdline_parser_params.check_ambiguity
 * @param override @see cmdline_parser_params.override
 * @param no_free whether to free a possible previous value
 * @param multiple_option whether this is a multiple option
 * @param long_opt the corresponding long option
 * @param short_opt the corresponding short option (or '-' if none)
 * @param additional_error possible further error specification
 */
static
int update_arg(void *field, char **orig_field,
               unsigned int *field_given, unsigned int *prev_given, 
               char *value, const char *possible_values[],
               const char *default_value,
               cmdline_parser_arg_type arg_type,
               int check_ambiguity, int override,
               int no_free, int multiple_option,
               const char *long_opt, char short_opt,
               const char *additional_error)
{
  char *stop_char = 0;
  const char *val = value;
  int found;
  FIX_UNUSED (field);

  stop_char = 0;
  found = 0;

  if (!multiple_option && prev_given && (*prev_given || (check_ambiguity && *field_given)))
    {
      if (short_opt != '-')
        fprintf (stderr, "%s: `--%s' (`-%c') option given more than once%s\n", 
               package_name, long_opt, short_opt,
               (additional_error ? additional_error : ""));
      else
        fprintf (stderr, "%s: `--%s' option given more than once%s\n", 
               package_name, long_opt,
               (additional_error ? additional_error : ""));
      return 1; /* failure */
    }

  FIX_UNUSED (default_value);
    
  if (field_given && *field_given && ! override)
    return 0;
  if (prev_given)
    (*prev_given)++;
  if (field_given)
    (*field_given)++;
  if (possible_values)
    val = possible_values[found];

  switch(arg_type) {
  case ARG_FLAG:
    *((int *)field) = !*((int *)field);
    break;
  case ARG_INT:
    if (val) *((int *)field) = strtol (val, &stop_char, 0);
    break;
  case ARG_LONG:
    if (val) *((long *)field) = (long)strtol (val, &stop_char, 0);
    break;
  case ARG_FLOAT:
    if (val) *((float *)field) = (float)strtod (val, &stop_char);
    break;
  case ARG_DOUBLE:
    if (val) *((double *)field) = strtod (val, &stop_char);
    break;
  case ARG_LONGLONG:
#ifdef HAVE_LONG_LONG
    if (val) *((long long int*)field) = (long long int) strtol (val, &stop_char, 0);
#else
    if (val) *((long *)field) = (long)strtol (val, &stop_char, 0);
#endif
    break;
  default:
    break;
  };

  /* check numeric conversion */
  switch(arg_type) {
  case ARG_INT:
  case ARG_LONG:
  case ARG_FLOAT:
  case ARG_DOUBLE:
  case ARG_LONGLONG:
    if (val && !(stop_char && *stop_char == '\0')) {
      fprintf(stderr, "%s: invalid numeric value: %s\n", package_name, val);
      return 1; /* failure */
    }
    break;
  default:
    ;
  };

  /* store the original value */
  switch(arg_type) {
  case ARG_NO:
  case ARG_FLAG:
    break;
  default:
    if (value && orig_field) {
      if (no_free) {
        *orig_field = value;
      } else {
        if (*orig_field)
          free (*orig_field); /* free previous string */
        *orig_field = gengetopt_strdup (value);
      }
    }
  };

  return 0; /* OK */
}


int
cmdline_parser_internal (
  int argc, char **argv, struct gengetopt_args_info *args_info,
                        struct cmdline_parser_params *params, const char *additional_error)
{
  int c;	/* Character of the parsed option.  */

  int error = 0;
  struct gengetopt_args_info local_args_info;
  
  int override;
  int initialize;
  int check_required;
  int check_ambiguity;
  
  package_name = argv[0];
  
  override = params->override;
  initialize = params->initialize;
  check_required = params->check_required;
  check_ambiguity = params->check_ambiguity;

  if (initialize)
    cmdline_parser_init (args_info);

  cmdline_parser_init (&local_args_info);

  optarg = 0;
  optind = 0;
  opterr = params->print_errors;
  optopt = '?';

  while (1)
    {
      int option_index = 0;

      static struct option long_options[] = {
        { "help",	0, NULL, 'h' },
        { "version",	0, NULL, 'V' },
        { "bins",	1, NULL, 'b' },
        { "checksteps",	1, NULL, 'c' },
        { "elow",	1, NULL, 0 },
        { "ehigh",	1, NULL, 0 },
        { "flat",	1, NULL, 0 },
        { "info",	0, NULL, 0 },
        { "max",	1, NULL, 'm' },
        { "mod",	1, NULL, 'f' },
        { "norm",	1, NULL, 'n' },
        { "resolution",	1, NULL, 'r' },
        { "steplimit",	1, NULL, 'l' },
        { "seed",	1, NULL, 'S' },
        { "Temp",	1, NULL, 'T' },
        { "verbose",	0, NULL, 'v' },
        { "debug",	0, NULL, 'd' },
        { 0,  0, 0, 0 }
      };

      c = getopt_long (argc, argv, "hVb:c:m:f:n:r:l:S:T:vd", long_options, &option_index);

      if (c == -1) break;	/* Exit from `while (1)' loop.  */

      switch (c)
        {
        case 'h':	/* Print help and exit.  */
          cmdline_parser_print_help ();
          cmdline_parser_free (&local_args_info);
          exit (EXIT_SUCCESS);

        case 'V':	/* Print version and exit.  */
          cmdline_parser_print_version ();
          cmdline_parser_free (&local_args_info);
          exit (EXIT_SUCCESS);

        case 'b':	/* Number of (equidistant) histogram bins.  */
        
        
          if (update_arg( (void *)&(args_info->bins_arg), 
               &(args_info->bins_orig), &(args_info->bins_given),
              &(local_args_info.bins_given), optarg, 0, "100", ARG_INT,
              check_ambiguity, override, 0, 0,
              "bins", 'b',
              additional_error))
            goto failure;
        
          break;
        case 'c':	/* Number of Wang-Landau steps before histogram is checked for flatness.  */
        
        
          if (update_arg( (void *)&(args_info->checksteps_arg), 
               &(args_info->checksteps_orig), &(args_info->checksteps_given),
              &(local_args_info.checksteps_given), optarg, 0, "1000000", ARG_LONGLONG,
              check_ambiguity, override, 0, 0,
              "checksteps", 'c',
              additional_error))
            goto failure;
        
          break;
        case 'm':	/* Upper energy bound for sampling.  */
        
        
          if (update_arg( (void *)&(args_info->max_arg), 
               &(args_info->max_orig), &(args_info->max_given),
              &(local_args_info.max_given), optarg, 0, 0, ARG_DOUBLE,
              check_ambiguity, override, 0, 0,
              "max", 'm',
              additional_error))
            goto failure;
        
          break;
        case 'f':	/* Final value of Wang-Landau modification factor.  */
        
        
          if (update_arg( (void *)&(args_info->mod_arg), 
               &(args_info->mod_orig), &(args_info->mod_given),
              &(local_args_info.mod_given), optarg, 0, 0, ARG_DOUBLE,
              check_ambiguity, override, 0, 0,
              "mod", 'f',
              additional_error))
            goto failure;
        
          break;
        case 'n':	/* Number of bins used for normalization.  */
        
        
          if (update_arg( (void *)&(args_info->norm_arg), 
               &(args_info->norm_orig), &(args_info->norm_given),
              &(local_args_info.norm_given), optarg, 0, 0, ARG_INT,
              check_ambiguity, override, 0, 0,
              "norm", 'n',
              additional_error))
            goto failure;
        
          break;
        case 'r':	/* Sampling resolution (histogram bin width).  */
        
        
          if (update_arg( (void *)&(args_info->resolution_arg), 
               &(args_info->resolution_orig), &(args_info->resolution_given),
              &(local_args_info.resolution_given), optarg, 0, "0.5", ARG_DOUBLE,
              check_ambiguity, override, 0, 0,
              "resolution", 'r',
              additional_error))
            goto failure;
        
          break;
        case 'l':	/* Maximum number of MC steps to perform.  */
        
        
          if (update_arg( (void *)&(args_info->steplimit_arg), 
               &(args_info->steplimit_orig), &(args_info->steplimit_given),
              &(local_args_info.steplimit_given), optarg, 0, "100000000", ARG_LONGLONG,
              check_ambiguity, override, 0, 0,
              "steplimit", 'l',
              additional_error))
            goto failure;
        
          break;
        case 'S':	/* Seed for random number generation.  */
        
        
          if (update_arg( (void *)&(args_info->seed_arg), 
               &(args_info->seed_orig), &(args_info->seed_given),
              &(local_args_info.seed_given), optarg, 0, 0, ARG_LONG,
              check_ambiguity, override, 0, 0,
              "seed", 'S',
              additional_error))
            goto failure;
        
          break;
        case 'T':	/* Simulation temperature in Celsius (currently n/a).  */
        
        
          if (update_arg( (void *)&(args_info->Temp_arg), 
               &(args_info->Temp_orig), &(args_info->Temp_given),
              &(local_args_info.Temp_given), optarg, 0, 0, ARG_FLOAT,
              check_ambiguity, override, 0, 0,
              "Temp", 'T',
              additional_error))
            goto failure;
        
          break;
        case 'v':	/* Verbose output.  */
        
        
          if (update_arg((void *)&(args_info->verbose_flag), 0, &(args_info->verbose_given),
              &(local_args_info.verbose_given), optarg, 0, 0, ARG_FLAG,
              check_ambiguity, override, 1, 0, "verbose", 'v',
              additional_error))
            goto failure;
        
          break;
        case 'd':	/* Debugging output.  */
        
        
          if (update_arg((void *)&(args_info->debug_flag), 0, &(args_info->debug_given),
              &(local_args_info.debug_given), optarg, 0, 0, ARG_FLAG,
              check_ambiguity, override, 1, 0, "debug", 'd',
              additional_error))
            goto failure;
        
          break;

        case 0:	/* Long option with no short option */
          /* Lower limit of sampling window (currently n/a).  */
          if (strcmp (long_options[option_index].name, "elow") == 0)
          {
          
          
            if (update_arg( (void *)&(args_info->elow_arg), 
                 &(args_info->elow_orig), &(args_info->elow_given),
                &(local_args_info.elow_given), optarg, 0, 0, ARG_DOUBLE,
                check_ambiguity, override, 0, 0,
                "elow", '-',
                additional_error))
              goto failure;
          
          }
          /* Upper limit of sampling window (currently n/a).  */
          else if (strcmp (long_options[option_index].name, "ehigh") == 0)
          {
          
          
            if (update_arg( (void *)&(args_info->ehigh_arg), 
                 &(args_info->ehigh_orig), &(args_info->ehigh_given),
                &(local_args_info.ehigh_given), optarg, 0, 0, ARG_DOUBLE,
                check_ambiguity, override, 0, 0,
                "ehigh", '-',
                additional_error))
              goto failure;
          
          }
          /* Flatness criterion for the histogram.  */
          else if (strcmp (long_options[option_index].name, "flat") == 0)
          {
          
          
            if (update_arg( (void *)&(args_info->flat_arg), 
                 &(args_info->flat_orig), &(args_info->flat_given),
                &(local_args_info.flat_given), optarg, 0, 0, ARG_FLOAT,
                check_ambiguity, override, 0, 0,
                "flat", '-',
                additional_error))
              goto failure;
          
          }
          /* Show settings.  */
          else if (strcmp (long_options[option_index].name, "info") == 0)
          {
          
          
            if (update_arg((void *)&(args_info->info_flag), 0, &(args_info->info_given),
                &(local_args_info.info_given), optarg, 0, 0, ARG_FLAG,
                check_ambiguity, override, 1, 0, "info", '-',
                additional_error))
              goto failure;
          
          }
          
          break;
        case '?':	/* Invalid option.  */
          /* `getopt_long' already printed an error message.  */
          goto failure;

        default:	/* bug: option not considered.  */
          fprintf (stderr, "%s: option unknown: %c%s\n", CMDLINE_PARSER_PACKAGE, c, (additional_error ? additional_error : ""));
          abort ();
        } /* switch */
    } /* while */




  cmdline_parser_release (&local_args_info);

  if ( error )
    return (EXIT_FAILURE);

  if (optind < argc)
    {
      int i = 0 ;
      int found_prog_name = 0;
      /* whether program name, i.e., argv[0], is in the remaining args
         (this may happen with some implementations of getopt,
          but surely not with the one included by gengetopt) */

      i = optind;
      while (i < argc)
        if (argv[i++] == argv[0]) {
          found_prog_name = 1;
          break;
        }
      i = 0;

      args_info->inputs_num = argc - optind - found_prog_name;
      args_info->inputs =
        (char **)(malloc ((args_info->inputs_num)*sizeof(char *))) ;
      while (optind < argc)
        if (argv[optind++] != argv[0])
          args_info->inputs[ i++ ] = gengetopt_strdup (argv[optind-1]) ;
    }

  return 0;

failure:
  
  cmdline_parser_release (&local_args_info);
  return (EXIT_FAILURE);
}
