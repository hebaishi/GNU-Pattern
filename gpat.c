/*************************************************************************
*   Module: gpat.c                                                       *
*                                                                        *
*   gpat main program                                                    *
*                                                                        *
*   This file is part of the gpat 1.0 distribution                       *
*                                                                        *
*     Copyright (C) 2001 - Ying Xu                                       *
*                                                                        *
*  This program is free software; you can redistribute it and/or modify  *
*  it under the terms of the GNU General Public License as published by  *
*  the Free Software Foundation; either version 2 of the License, or     *
*  (at your option) any later version.                                   *
*                                                                        *
*  This program is distributed in the hope that it will be useful,       *
*  but WITHOUT ANY WARRANTY; without even the implied warranty of        *
*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
*  GNU General Public License for more details.                          *
*                                                                        *
*  You should have received a copy of the GNU General Public License     *
*  along with this program; if not, write to the Free Software           *
*  Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.             *
*************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <fcntl.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <sys/times.h>
#include <sys/types.h>
#include <pthread.h>
#include <unistd.h>
#include <getopt.h>
#include "ASM1.h"

const char *program_name;
int patt_tmp_uid = 1;		//id number used in seed set
int patt_fix_uid = 1;		//id number used in result set
int patt_dis_uid = 1;		//id number used in output
int seqn, spcn;			//number of sequences in file &  number of spectrum combinations
int verbose = 0, no_ext = 0, need_pos = 1, need_patt = 0;
struct phase1_jobs *p1_jobs = NULL;

//mutex locks for global variables
pthread_mutex_t p1_jobs_mutex = PTHREAD_MUTEX_INITIALIZER;
pthread_mutex_t patt_fix_uid_mutex = PTHREAD_MUTEX_INITIALIZER;
pthread_mutex_t fsw_mutex = PTHREAD_MUTEX_INITIALIZER;


void
print_usage (FILE * stream, int exit_code)
{
  fprintf (stream, "Usage: %s options [inputfile...]\n", program_name);
  fprintf (stream,
	   "-i input filename                                        input data from file.\n"
	   "-o output filename          (default my_pattern.pat)     write output to file.\n"
	   "-k min tokens in seed       (default 5)                  minimum number of tokens in seed.\n"
	   "-l win width of seed        (default 8)                  window width when scanning seed.\n"
	   "-s min supports for seed    (default 2)                  minimum supports needed to be a seed.\n"
	   "-J min tokens in pattern    (default 1)                  minimum number of tokens in pattern.\n"
	   "-K min supports for pattern (default 2)                  minimum supports needed to be a pattern.\n"
	   "-p min percentage covered   (not used)                   minimum percentage of sequences covered by pattern.\n"
	   "-w min pattern width        (default 0)                  minimum pattern length that should be gained.\n"
	   "-t type sequence type       (default 2)                  dna or rna or protein.\n"
	   "-m multi-threaded           (default 1)                  switch to multi-threaded state\n"
	   "-n store-pattern till end   (default 0)                  store pattern information till thread dies\n"
	   "-h --help                                                display this usage information.\n"
	   "-v --verbose                                             print verbose messages.\n"
	   "-x no extension/quick count                              no extension to achieve quick count(k must equal l)\n");
  exit (exit_code);
}


int
main (int argc, char *argv[])
{
  int i, j, k, aculen;
  int next_option;
  const char *short_options = "i:o:k:l:s:J:K:p:w:t:m:nhvx";
  const struct option long_options[] = {
    {"input", required_argument, NULL, 'i'},
    {"output", required_argument, NULL, 'o'},
    {"min_s_tokens", required_argument, NULL, 'k'},
    {"win_width", required_argument, NULL, 'l'},
    {"support4seed", required_argument, NULL, 's'},
    {"min_p_tokens", required_argument, NULL, 'J'},
    {"support4pattern", required_argument, NULL, 'K'},
    {"percentage", required_argument, NULL, 'p'},
    {"patt_size", required_argument, NULL, 'w'},
    {"type", required_argument, NULL, 't'},
    {"thread", required_argument, NULL, 'm'},
    {"not_store", no_argument, NULL, 'n'},
    {"help", no_argument, NULL, 'h'},
    {"verbose", no_argument, NULL, 'v'},
    {"quick_count", no_argument, NULL, 'x'},
    {NULL, 0, NULL, 0}
  };
  const char *output_filename = NULL;
  program_name = argv[0];
  paramptr pptr = (paramptr) malloc (sizeof (struct parameters));
  int *nr;
  t_sequence **seqs;

  /* default values */
  pptr->input_file = NULL;
  pptr->output_file = "my_pattern.pat";
  pptr->win_width = 8;
  pptr->min_s_tokens = 5;
  pptr->min_p_tokens = 1;
  pptr->support4seed = 2;
  pptr->support4pattern = 2;
  pptr->percentage = 0;
  pptr->patt_size = 0;
  pptr->type = 2;
  pptr->thread = 1;

  do
    {
      next_option =
	getopt_long (argc, argv, short_options, long_options, NULL);
      switch (next_option)
	{
	case 'i':
	  pptr->input_file = (char *) optarg;
	  break;
	case 'o':
	  pptr->output_file = (char *) optarg;
	  break;
	case 'k':
	  pptr->min_s_tokens = atoi (optarg);
	  break;
	case 'l':
	  pptr->win_width = atoi (optarg);
	  break;
	case 's':
	  pptr->support4seed = atoi (optarg);
	  break;
	case 'J':
	  pptr->min_p_tokens = atoi (optarg);
	  break;
	case 'K':
	  pptr->support4pattern = atoi (optarg);
	  break;
	case 't':
	  pptr->type = atoi (optarg);
	  break;
	case 'm':
	  pptr->thread = atoi (optarg);
	  if (pptr->thread > 256)
	    pptr->thread = 256;
	  if (pptr->thread == 0)
	    pptr->thread = 1;
	  break;
	case 'p':
	  pptr->percentage = atoi (optarg);
	  break;
	case 'w':
	  pptr->patt_size = atoi (optarg);
	  break;
	case 'n':
	  need_patt = 1;
	  break;
	case 'h':
	  print_usage (stdout, 0);
	  break;
	case 'v':
	  verbose = 1;
	  break;
	case 'x':
	  if (pptr->min_s_tokens != pptr->win_width)
	    print_usage (stderr, 1);
	  no_ext = 1;
	  need_pos = 0;
	  break;
	case '?':
	  print_usage (stderr, 1);
	case -1:
	  break;
	default:
	  abort ();
	}
    }
  while (next_option != -1);

  if (!pptr->input_file || !pptr->output_file
      || pptr->win_width <= 0 || pptr->min_s_tokens <= 0
      || pptr->min_s_tokens > pptr->win_width
      || pptr->min_p_tokens <= 0 || pptr->support4seed < 0
      || pptr->support4pattern < 0 || pptr->percentage < 0
      || pptr->patt_size < 0 || pptr->type < 0 || pptr->type > 2
      || pptr->thread < 0 || pptr->thread > 256)
    print_usage (stderr, 1);

  if (pptr->min_s_tokens == 1 && pptr->min_s_tokens != pptr->win_width)
    pptr->min_s_tokens += 1;
  pptr->support4seed < pptr->support4pattern ? pptr->support4pattern : pptr->support4seed;

  seqs = seqs_from_fasta_file (pptr->input_file, &seqn, pptr->type);

  //calculate the number of spectrum combinations
  spcn = spec_num (pptr->min_s_tokens - 1, pptr->win_width - 1);

  printf ("\n>>> Gpat, open source pattern scanner under x86-linux.\n");
  printf (">>> Copyright (C) 2004, and GNU GPL'd, by Ying Xu.\n");
  printf
    (">>> For more comments, write to bio_xy@hotmail.com/xuying@sibs.ac.cn\n\n");
  printf ("k = %d, l = %d, s = %d, K = %d\n", pptr->min_s_tokens,
	  pptr->win_width, pptr->support4seed, pptr->support4pattern);
  printf ("Job = %d activated.\n", spcn);

  if (pptr->thread > 1)
    printf ("Thread = %d running...\n", pptr->thread);

  //to measure the process time
  struct tms tmsstart, tmsend;
  clock_t start, end;
  start = times (&tmsstart);

  btreeptr spc_tree, curlf, nxtlf, leaf =
    (btreeptr) malloc (sizeof (struct btree));
  leaf->left = NULL;
  leaf->right = NULL;

  //generate spectrum tree
  spc_tree =
    g_tree (leaf, NULL, pptr, pptr->min_s_tokens, pptr->win_width, 1);
  curlf = leaf->right;
  nxtlf = curlf->right;

  pthread_t thread_id[pptr->thread];
  struct pthread_param thread_arg[pptr->thread];

  //produce spectrum from spectrum tree
  for (i = 0; i < spcn; i++)
    {
      spcptr_list spec = (spcptr_list) malloc (sizeof (struct spc_list));
      spec->data.id = 0;
      spec->data.spseg = 0;
      spec->previous = NULL;
      spec->next = NULL;
      get_spectrum (curlf, spcn, spec);
      delete_branch (curlf);
      curlf = nxtlf;
      if (curlf)
	nxtlf = curlf->right;
      struct phase1_jobs *job =
	(struct phase1_jobs *) malloc (sizeof (struct phase1_jobs));
      job->job = spec;
      job->next = NULL;
      insert_job_queue (job);	//add the job (spectrum) into the job queue
    }

  FILE *fsw;
  if (!(fsw = fopen (pptr->output_file, "w+")))
    {
      perror ("open");
      return 1;
    }

  for (i = 0; i < pptr->thread; i++)
    {
      thread_arg[i].seqs = seqs;
      thread_arg[i].fsw = fsw;
      thread_arg[i].pptr = pptr;
      pthread_create (&thread_id[i], NULL, &thread_find, &thread_arg[i]);
    }

  //wait all threads to finish or global variables will be lost if the main thread quits first
  for (i = 0; i < pptr->thread; i++)
    pthread_join (thread_id[i], NULL);

  free (leaf);
  for (i = 0; i < seqn; i++)
    {
      free (seqs[i]->seq);
      free (seqs[i]->filename);
      free (seqs[i]);
    }
  free (seqs);

  end = times (&tmsend);
  printf ("\nruntime report\n");
  pr_times (end - start, &tmsstart, &tmsend);

  //close (fdr);
  //fclose (fsr);
  //remove (tmp);
  //free (tmp);
  pptr->input_file = "";
  pptr->output_file = "";
  free (pptr);
  printf ("\nIT IS DONE...\n\n");

  return 0;
}
