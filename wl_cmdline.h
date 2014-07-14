/** @file wl_cmdline.h
 *  @brief The header file for the command line option parser
 *  generated by GNU Gengetopt version 2.22.5
 *  http://www.gnu.org/software/gengetopt.
 *  DO NOT modify this file, since it can be overwritten
 *  @author GNU Gengetopt by Lorenzo Bettini */

#ifndef WL_CMDLINE_H
#define WL_CMDLINE_H

/* If we use autoconf.  */
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdio.h> /* for FILE */

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

#ifndef CMDLINE_PARSER_PACKAGE
/** @brief the program name (used for printing errors) */
#define CMDLINE_PARSER_PACKAGE PACKAGE
#endif

#ifndef CMDLINE_PARSER_PACKAGE_NAME
/** @brief the complete program name (used for help and version) */
#ifdef PACKAGE_NAME
#define CMDLINE_PARSER_PACKAGE_NAME PACKAGE_NAME
#else
#define CMDLINE_PARSER_PACKAGE_NAME PACKAGE
#endif
#endif

#ifndef CMDLINE_PARSER_VERSION
/** @brief the program version */
#define CMDLINE_PARSER_VERSION VERSION
#endif

/** @brief Where the command line options are stored */
struct gengetopt_args_info
{
  const char *help_help; /**< @brief Print help and exit help description.  */
  const char *version_help; /**< @brief Print version and exit help description.  */
  int bins_arg;	/**< @brief Number of (equidistant) histogram bins.  */
  char * bins_orig;	/**< @brief Number of (equidistant) histogram bins original value given at command line.  */
  const char *bins_help; /**< @brief Number of (equidistant) histogram bins help description.  */
  #ifdef HAVE_LONG_LONG
  long long int checksteps_arg;	/**< @brief Number of Wang-Landau steps before histogram is checked for flatness (default=100000).  */
  #else
  long checksteps_arg;	/**< @brief Number of Wang-Landau steps before histogram is checked for flatness (default=100000).  */
  #endif
  char * checksteps_orig;	/**< @brief Number of Wang-Landau steps before histogram is checked for flatness original value given at command line.  */
  const char *checksteps_help; /**< @brief Number of Wang-Landau steps before histogram is checked for flatness help description.  */
  double elow_arg;	/**< @brief Lower limit of sampling window (currently n/a).  */
  char * elow_orig;	/**< @brief Lower limit of sampling window (currently n/a) original value given at command line.  */
  const char *elow_help; /**< @brief Lower limit of sampling window (currently n/a) help description.  */
  double ehigh_arg;	/**< @brief Upper limit of sampling window (currently n/a).  */
  char * ehigh_orig;	/**< @brief Upper limit of sampling window (currently n/a) original value given at command line.  */
  const char *ehigh_help; /**< @brief Upper limit of sampling window (currently n/a) help description.  */
  float flat_arg;	/**< @brief Flatness criterion for the histogram.  */
  char * flat_orig;	/**< @brief Flatness criterion for the histogram original value given at command line.  */
  const char *flat_help; /**< @brief Flatness criterion for the histogram help description.  */
  int info_flag;	/**< @brief Show settings (default=off).  */
  const char *info_help; /**< @brief Show settings help description.  */
  double max_arg;	/**< @brief Upper energy bound for sampling.  */
  char * max_orig;	/**< @brief Upper energy bound for sampling original value given at command line.  */
  const char *max_help; /**< @brief Upper energy bound for sampling help description.  */
  double mod_arg;	/**< @brief Final value of Wang-Landau modification factor.  */
  char * mod_orig;	/**< @brief Final value of Wang-Landau modification factor original value given at command line.  */
  const char *mod_help; /**< @brief Final value of Wang-Landau modification factor help description.  */
  int norm_arg;	/**< @brief Number of bins used for normalization.  */
  char * norm_orig;	/**< @brief Number of bins used for normalization original value given at command line.  */
  const char *norm_help; /**< @brief Number of bins used for normalization help description.  */
  double resolution_arg;	/**< @brief Sampling resolution (histogram bin width) (default='0.5').  */
  char * resolution_orig;	/**< @brief Sampling resolution (histogram bin width) original value given at command line.  */
  const char *resolution_help; /**< @brief Sampling resolution (histogram bin width) help description.  */
  #ifdef HAVE_LONG_LONG
  long long int steplimit_arg;	/**< @brief Maximum number of MC steps to perform (default=100000000).  */
  #else
  long steplimit_arg;	/**< @brief Maximum number of MC steps to perform (default=100000000).  */
  #endif
  char * steplimit_orig;	/**< @brief Maximum number of MC steps to perform original value given at command line.  */
  const char *steplimit_help; /**< @brief Maximum number of MC steps to perform help description.  */
  long seed_arg;	/**< @brief Seed for random number generation.  */
  char * seed_orig;	/**< @brief Seed for random number generation original value given at command line.  */
  const char *seed_help; /**< @brief Seed for random number generation help description.  */
  float Temp_arg;	/**< @brief Temperatur in Celsius.  */
  char * Temp_orig;	/**< @brief Temperatur in Celsius original value given at command line.  */
  const char *Temp_help; /**< @brief Temperatur in Celsius help description.  */
  int verbose_flag;	/**< @brief Verbose output (default=off).  */
  const char *verbose_help; /**< @brief Verbose output help description.  */
  int debug_flag;	/**< @brief Debugging output (default=off).  */
  const char *debug_help; /**< @brief Debugging output help description.  */
  
  unsigned int help_given ;	/**< @brief Whether help was given.  */
  unsigned int version_given ;	/**< @brief Whether version was given.  */
  unsigned int bins_given ;	/**< @brief Whether bins was given.  */
  unsigned int checksteps_given ;	/**< @brief Whether checksteps was given.  */
  unsigned int elow_given ;	/**< @brief Whether elow was given.  */
  unsigned int ehigh_given ;	/**< @brief Whether ehigh was given.  */
  unsigned int flat_given ;	/**< @brief Whether flat was given.  */
  unsigned int info_given ;	/**< @brief Whether info was given.  */
  unsigned int max_given ;	/**< @brief Whether max was given.  */
  unsigned int mod_given ;	/**< @brief Whether mod was given.  */
  unsigned int norm_given ;	/**< @brief Whether norm was given.  */
  unsigned int resolution_given ;	/**< @brief Whether resolution was given.  */
  unsigned int steplimit_given ;	/**< @brief Whether steplimit was given.  */
  unsigned int seed_given ;	/**< @brief Whether seed was given.  */
  unsigned int Temp_given ;	/**< @brief Whether Temp was given.  */
  unsigned int verbose_given ;	/**< @brief Whether verbose was given.  */
  unsigned int debug_given ;	/**< @brief Whether debug was given.  */

  char **inputs ; /**< @brief unamed options (options without names) */
  unsigned inputs_num ; /**< @brief unamed options number */
} ;

/** @brief The additional parameters to pass to parser functions */
struct cmdline_parser_params
{
  int override; /**< @brief whether to override possibly already present options (default 0) */
  int initialize; /**< @brief whether to initialize the option structure gengetopt_args_info (default 1) */
  int check_required; /**< @brief whether to check that all required options were provided (default 1) */
  int check_ambiguity; /**< @brief whether to check for options already specified in the option structure gengetopt_args_info (default 0) */
  int print_errors; /**< @brief whether getopt_long should print an error message for a bad option (default 1) */
} ;

/** @brief the purpose string of the program */
extern const char *gengetopt_args_info_purpose;
/** @brief the usage string of the program */
extern const char *gengetopt_args_info_usage;
/** @brief all the lines making the help output */
extern const char *gengetopt_args_info_help[];

/**
 * The command line parser
 * @param argc the number of command line options
 * @param argv the command line options
 * @param args_info the structure where option information will be stored
 * @return 0 if everything went fine, NON 0 if an error took place
 */
int cmdline_parser (int argc, char **argv,
  struct gengetopt_args_info *args_info);

/**
 * The command line parser (version with additional parameters - deprecated)
 * @param argc the number of command line options
 * @param argv the command line options
 * @param args_info the structure where option information will be stored
 * @param override whether to override possibly already present options
 * @param initialize whether to initialize the option structure my_args_info
 * @param check_required whether to check that all required options were provided
 * @return 0 if everything went fine, NON 0 if an error took place
 * @deprecated use cmdline_parser_ext() instead
 */
int cmdline_parser2 (int argc, char **argv,
  struct gengetopt_args_info *args_info,
  int override, int initialize, int check_required);

/**
 * The command line parser (version with additional parameters)
 * @param argc the number of command line options
 * @param argv the command line options
 * @param args_info the structure where option information will be stored
 * @param params additional parameters for the parser
 * @return 0 if everything went fine, NON 0 if an error took place
 */
int cmdline_parser_ext (int argc, char **argv,
  struct gengetopt_args_info *args_info,
  struct cmdline_parser_params *params);

/**
 * Save the contents of the option struct into an already open FILE stream.
 * @param outfile the stream where to dump options
 * @param args_info the option struct to dump
 * @return 0 if everything went fine, NON 0 if an error took place
 */
int cmdline_parser_dump(FILE *outfile,
  struct gengetopt_args_info *args_info);

/**
 * Save the contents of the option struct into a (text) file.
 * This file can be read by the config file parser (if generated by gengetopt)
 * @param filename the file where to save
 * @param args_info the option struct to save
 * @return 0 if everything went fine, NON 0 if an error took place
 */
int cmdline_parser_file_save(const char *filename,
  struct gengetopt_args_info *args_info);

/**
 * Print the help
 */
void cmdline_parser_print_help(void);
/**
 * Print the version
 */
void cmdline_parser_print_version(void);

/**
 * Initializes all the fields a cmdline_parser_params structure 
 * to their default values
 * @param params the structure to initialize
 */
void cmdline_parser_params_init(struct cmdline_parser_params *params);

/**
 * Allocates dynamically a cmdline_parser_params structure and initializes
 * all its fields to their default values
 * @return the created and initialized cmdline_parser_params structure
 */
struct cmdline_parser_params *cmdline_parser_params_create(void);

/**
 * Initializes the passed gengetopt_args_info structure's fields
 * (also set default values for options that have a default)
 * @param args_info the structure to initialize
 */
void cmdline_parser_init (struct gengetopt_args_info *args_info);
/**
 * Deallocates the string fields of the gengetopt_args_info structure
 * (but does not deallocate the structure itself)
 * @param args_info the structure to deallocate
 */
void cmdline_parser_free (struct gengetopt_args_info *args_info);

/**
 * Checks that all the required options were specified
 * @param args_info the structure to check
 * @param prog_name the name of the program that will be used to print
 *   possible errors
 * @return
 */
int cmdline_parser_required (struct gengetopt_args_info *args_info,
  const char *prog_name);


#ifdef __cplusplus
}
#endif /* __cplusplus */
#endif /* WL_CMDLINE_H */
