/*************************************************************************
*   Module: ASM1.h                                                       *
*                                                                        *
*   gpat headfile of application specific module 1                       *
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

#ifndef ASM1_DEF
#define ASM1_DEF

/* special structs used in pattern discovery */

typedef enum
{ dna, rna, protein, unknown }
t_seq_type;

typedef struct t_sequence
{
  int index;
  char *filename;
  t_seq_type type;
  int length;
  int max_length;
  char *seq;
  struct t_sequence *next;
}
t_sequence, *t_sequence_ptr;

/* all parameters are here */
typedef struct parameters
{
  char *input_file;
  char *output_file;
  int win_width;		/* minimum window width */
  int min_s_tokens;		/* minimum tokens in seed */
  int min_p_tokens;		/* minimum tokens in pattern */
  int support4seed;		/* minimum supports for seed */
  int support4pattern;		/* minimum supports for pattern */
  int percentage;		/* minimum percentage in sequences */
  int patt_size;
  t_seq_type type;		/* 1=DNA 2=RNA 3=PROTEIN 4=TEXT */
  int thread;
}
parameters, *paramptr;

/* tree node for spectrum generation */
typedef struct btree_data
{
  int nlay;
  int nzero;
  int npos;
}
btree_data, *btreeptr_data;

/* tree structure for spectrum */
typedef struct btree
{
  struct btree_data data;
  struct btree *previous;
  struct btree *left;
  struct btree *right;
}
btree, *btreeptr;

/* divide spectrum into many smaller pieces */
typedef struct spc_node
{
  int id;
  unsigned int spseg;
}
spc_node, *spcptr_node;

/* use double linked list to realize the whole spectrum */
typedef struct spc_list
{
  struct spc_node data;
  struct spc_list *previous;
  struct spc_list *next;
}
spc_list, *spcptr_list;

/* sequence index structure */
typedef struct seq_idx
{
  int id;
  int left;
  int right;
  struct seq_idx *next;
}
seq_idx, *seqptr_idx;

/* one position node */
typedef struct position
{
  int seqnum;			/* in which sequence */
  int offset;			/* offset in the sequence */
}
position, *posptr;

/* list of positions */
typedef struct position_list
{
  struct position data;
  struct position_list *next;
}
position_list, *poslist;

/* pattern details */
typedef struct pattern
{
  int id;			/* unique number to identify the pattern */
  unsigned char *text;		/* string of tokens in the pattern */
  struct spc_list *spectrum;	/* double linked list version */
  int support;			/* supports of the pattern */
  int size;			/* the span of the pattern */
  struct position_list *location;	/* position_list of the pattern */
  struct position_list *rear;	/* point to the tail of position list */
}
pattern, *pattptr;

/* list of patterns */
typedef struct pattern_list
{
  struct pattern data;
  struct pattern_list *next;
}
pattern_list, *pattlist;

typedef struct phase1_jobs
{
  spcptr_list job;
  struct phase1_jobs *next;
}
phase1_jobs;

typedef struct pthread_param
{
  t_sequence **seqs;
  FILE *fsw;
  paramptr pptr;
}
pthread_param;


void pr_times (clock_t, struct tms *, struct tms *);

t_sequence *seq_init_sequence (int len, t_seq_type type, char *name);

static void seq_expand (t_sequence *s);

void seq_add_res (t_sequence *s, char t);

t_sequence **seqs_from_fasta_file (char *filename, int *seqn,
				   t_seq_type type);
/* count seq num from fasta file */
int count_seq (int fd);

/* get the size of alphabet */
int alph_num (int type);

char get_char (int k, int type);

/* the function for formatting long int number (not in used!) */
void formal (int *a);

/* the function for looking for the greatest common divider */
int gcd (int p, int q);

/* function for computing the num of spectrums */
unsigned int spec_num (int m, int n);

/* function for tree node insertion */
int insert_tree_node (btreeptr prev, btreeptr node, int tag);

/* function for spectrum tree generation */
btreeptr g_tree (btreeptr leaves, btreeptr prev, paramptr pptr, int m, int n,
		 int tag);

/* function for generating valid spectrums */
int get_spectrum (btreeptr leaves, int c, spcptr_list spec);

/* function for spc tree deletion */
int delete_spc_tree (btreeptr spc_tree);

/* calculate region for l-max judgement */
int getlx (paramptr pp, pattlist p);

/* calculate region for l-max judgement */
int getrx (paramptr pp, pattlist p, t_sequence **seqs);

/* get the gap num for the size of the gap array */
int gap_num (pattlist p);

/* get the token num for the size of the ngap array */
int ngap_num (pattlist p);

/* is compositional maximality? */
int not_cmax (int gp, pattlist p, t_sequence **seqs);

/* is left maximality? */
int not_lmax (int lx, pattlist p, t_sequence **seqs);

/* get gap position in spectrum */
void get_gap (int gp[], spcptr_list sp);

/* get token position in spectrum */
void get_ngap (int ngp[], spcptr_list sp);

/* count number of positions in pattern */
int count_pos (pattlist p);

/* is lc-maximal? */
int is_lcmaximal (paramptr pp, pattlist p, t_sequence **seqs);

/* plan to replace all of these with stl functions */

int is_equal_spc (spcptr_list spc1, spcptr_list spc2);

int insert_spc_node (spcptr_list spc, spcptr_list pspc);

int delete_spc_node (spcptr_list pspc);

void print_spectrum (FILE * fsw, spcptr_list spc);

int delete_spectrum (spcptr_list spc);

int insert_position (poslist pl, poslist rear, poslist p);

int delete_position (poslist pl, poslist p);

void print_position_list (FILE * fsw, poslist pl);

int delete_position_list (poslist pl);

int insert_pattern_rear (pattlist pl, pattlist r, pattlist p, int f);

int insert_pattern_after (pattlist p1, pattlist p2);

int delete_pattern (pattlist pl, pattlist p, int f);

void print_text (FILE * fsw, pattlist p);

void print_pattern_list (FILE * fsw, pattlist pl);

int delete_pattern_list (pattlist pl);

pattlist create_pattern_list (char *s);

/* determine the position of the given character in corrresponding alphebat array */
int get_ox (char x, char alph[]);

/* compute the pattern size from spectrum */
int get_size (spcptr_list spc);

/* clone spectrum */
spcptr_list clone_spec (spcptr_list spc);

/* just clone this node */
poslist clone_pos (poslist pl);

/* create a new pattern from the given spectrum */
pattlist mak_patt (char *d, spcptr_list spc, int support, int size,
		   poslist pos);

/* copy partial data from p2 to p1 */
void clone_patt (pattlist p1, pattlist p2);

/* clone all data from pl */
pattlist exactly_clone_patt (pattlist pl);

/* get the size of the blank region in pspc */
int get_blank (unsigned int spseg);

/* get the block of data that will be shifted out */
unsigned int get_flow_block (unsigned int spseg, int rx);

/* transport bits to the left */
void teleport (spcptr_list pspc, int rx, int n);

/* extend the pattern to a longer one and update the corresponding pattern information */
void update_patt (pattlist p, char c, int rx);

/* see whether there already exists this seed, if there is one, update seed information, if not, add it */
void check_and_update (char *d, spcptr_list spc, poslist pos, pattlist hs,
		       pattlist hr, paramptr pptr);

/* function for generating seed */
void generate_seed (t_sequence **seqs, paramptr pptr, spcptr_list spec, pattlist seed);

/* $$$ core function $$$ */
void find_patt (t_sequence **seqs, FILE * fsw, paramptr pptr,
		pattlist seed, pattlist result_set);

void insert_job_queue (struct phase1_jobs *job);

int process_phase1_jobs (struct phase1_jobs *next_job,
			 struct pthread_param *p);

void *thread_find (void *param);

#endif
