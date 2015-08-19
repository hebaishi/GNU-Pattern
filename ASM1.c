/*************************************************************************
*   Module: ASM1.c                                                       *
*                                                                        *
*   gpat application specific module 1                                   *
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

#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <fcntl.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <sys/times.h>
#include <sys/types.h>
#include <time.h>
#include <unistd.h>
#include <math.h>

#define ASM1_FUNC
#include "ASM1.h"

#define UINT_SIZE 32
#define MAXN 64
#define min(a,b) (a<=b)?a:b
#define max(a,b) (a>=b)?a:b


/* set array size to num of letters + 1, strlen(alph) will return the right num */
char NA_alphabet[] = "ACGTN\0";	//ACGTN
char EXTNA_alphabet[] = "ACGTN-BDHKMRSVWY\0";	//ACGTN-BDHKMRSVWY
char PROTEIN_alphabet[] = "ACGTN-BDHKMRSVWYLEFIOPQUXZ*\0";	//ACGTN-BDHKMRSVWYLEFIOPQUXZ*

extern int patt_tmp_uid;
extern int patt_fix_uid, patt_dis_uid;
extern int seqn, spcn, verbose, no_ext, need_pos, need_patt;
btreeptr btrear;
int job_count = 0;

extern struct phase1_jobs *p1_jobs;
extern pthread_mutex_t p1_jobs_mutex, patt_fix_uid_mutex, fsw_mutex;

void dummy()
{
  printf ("dummy\n");
}

/* assign an index number to every letters in the alphebat */
int
get_hash_key (char *d, int n, int alphn, int f)
{
  int i, j, m = 0, h = 0;
  for (i = 0; i < n; i++)
    {
      switch (d[i])
	{
	case 'A':
	  m = 0;
	  break;
	case 'C':
	  m = 1;
	  break;
	case 'G':
	  m = 2;
	  break;
	case 'T':
	  m = 3;
	  break;
	case 'N':
	  m = 4;
	  break;
	case '-':
	  m = 5;
	  break;
	case 'B':
	  m = 6;
	  break;
	case 'D':
	  m = 7;
	  break;
	case 'H':
	  m = 8;
	  break;
	case 'K':
	  m = 9;
	  break;
	case 'M':
	  m = 10;
	  break;
	case 'R':
	  m = 11;
	  break;
	case 'S':
	  m = 12;
	  break;
	case 'V':
	  m = 13;
	  break;
	case 'W':
	  m = 14;
	  break;
	case 'Y':
	  m = 15;
	  break;
	case 'L':
	  m = 16;
	  break;
	case 'E':
	  m = 17;
	  break;
	case 'F':
	  m = 18;
	  break;
	case 'I':
	  m = 19;
	  break;
	case 'O':
	  m = 20;
	  break;
	case 'P':
	  m = 21;
	  break;
	case 'Q':
	  m = 22;
	  break;
	case 'U':
	  m = 23;
	  break;
	case 'X':
	  m = 24;
	  break;
	case 'Z':
	  m = 25;
	  break;
	case '*':
	  m = 26;
	  break;
	}
      if (m + 1 > alphn)
	{
	  printf ("no such symbol (%c)\n", d[i]);
	  m = 0;
	  exit (0);
	}
      if (f == 0)
	break;
      h += m * pow (alphn, i);
    }
  if (f == 1)
    return h;
  else
    return m;
}

void
pr_times (clock_t real, struct tms *tmsstart, struct tms *tmsend)
{
  static long clktck = 0;
  if (clktck == 0)		/* fetch clock ticks per second first time */
    if ((clktck = sysconf (_SC_CLK_TCK)) < 0)
      printf ("sysconf error\n");
  fprintf (stderr, "real:  %7.2f\n", real / (double) clktck);
  fprintf (stderr, "user:  %7.2f\n",
	   (tmsend->tms_utime - tmsstart->tms_utime) / (double) clktck);
  fprintf (stderr, "sys:   %7.2f\n",
	   (tmsend->tms_stime - tmsstart->tms_stime) / (double) clktck);
  fprintf (stderr, "child user:  %7.2f\n",
	   (tmsend->tms_cutime - tmsstart->tms_cutime) / (double) clktck);
  fprintf (stderr, "child sys:   %7.2f\n",
	   (tmsend->tms_cstime - tmsstart->tms_cstime) / (double) clktck);
}

t_sequence *
seq_init_sequence (int len, t_seq_type type, char *name)
{
  t_sequence *temp;
  temp = (t_sequence *) calloc (1, sizeof (t_sequence));
  temp->filename = (char *) calloc (strlen (name) + 1, sizeof (char));
  strcpy (temp->filename, name);
  temp->length = 0;
  temp->max_length = len;
  temp->seq = (char *) calloc (len, sizeof (char));
  temp->type = type;
  return temp;
}

static void
seq_expand (t_sequence * s)
{
  s->max_length *= 2;
  if (!(s->seq = (char *) realloc (s->seq, sizeof (char) * s->max_length)))
    printf ("seq_expand_error\n");
}

void
seq_add_res (t_sequence * s, char t)
{
  if (s->length == s->max_length)
    seq_expand (s);
  s->seq[s->length] = t;
  s->length++;
}

t_sequence **
seqs_from_fasta_file (char *filename, int *seqn, t_seq_type type)
{
  FILE *fil;
  t_sequence **seqs;
  t_sequence *seq;
  char buf[100];
  char name[20];
  int line;
  int max_nr;
  int i;
  char *ok;

  fil = fopen (filename, "r");
  if (fil == NULL)
    {
      printf ("Could not open file %s\n", filename);
      return NULL;
    }
  max_nr = 10000;
  seqs = (t_sequence **) calloc (max_nr, sizeof (t_sequence *));

  line = 0;
  *seqn = 0;
  ok = buf;
  ok = fgets (buf, 100, fil);
  while (ok != NULL)
    {
      while (buf[0] != '>' && ok)
	{
	  ok = fgets (buf, 100, fil);
	  line++;
	}
      if (!ok)
	continue;
      sscanf (buf + 1, "%s", name);
      seq = seq_init_sequence (400, type, name);
      while ((ok = fgets (buf, 100, fil)) && (buf[0] != '>'))
	{
	  line++;
	  for (i = 0; i < (int) strlen (buf); i++)
	    if (isupper (buf[i]) || islower (buf[i]))
	      {
		if (islower (buf[i]))
		  buf[i] = toupper (buf[i]);
		seq_add_res (seq, buf[i]);
	      }
	}
      seqs[(*seqn)] = seq;
      (*seqn)++;
    }
  return seqs;
}

/* count seq num from fasta file */
int
count_seq (int fd)
{
  int i = 0, j;
  size_t bytes_read = 0;
  unsigned char buf[1024];
  lseek (fd, 0, SEEK_SET);
  do
    {
      bytes_read = read (fd, buf, sizeof (buf));
      for (j = 0; j < bytes_read; j++)
	if (buf[j] == '>')
	  i++;
    }
  while (bytes_read == sizeof (buf));
  return i;
}

/* get the size of alphabet */
int
alph_num (int type)
{
  switch (type)
    {
    case 0:
      return 5;
    case 1:
      return 17;
    case 2:
      return 27;
    default:
      abort ();
    }
}

char
get_char (int k, int type)
{
  switch (type)
    {
    case 0:
      return NA_alphabet[k];
    case 1:
      return EXTNA_alphabet[k];
    case 2:
      return PROTEIN_alphabet[k];
    default:
      abort ();
    }
}

void
formal (int *a)
{
  int p;
  for (p = 1; p < a[0] || a[p] >= 10; p++)
    {
      if (p >= a[0])
	a[p + 1] = 0;
      a[p + 1] += a[p] / 10;
      a[p] = a[p] % 10;
    }
  if (p > a[0])
    a[0] = p;
}

/* calculate greatest common divisor */
int
gcd (int p, int q)
{
  int r;
  while (q > 0)
    {
      r = p % q;
      p = q;
      q = r;
    }
  return p;
}

/* calculate the number of spectrum combinations */
unsigned int
spec_num (int m, int n)
{
  if (m == 0)
    return 1;
  int i, j, k, p = m, q = n, x;
  unsigned int pu = 1.0, pd = 1.0;
  if ((double) n / (double) m < 2.0)
    p = q - p;
  int d[MAXN], u[MAXN];
  for (k = 0, i = q; i >= q - p + 1; i--)
    u[++k] = i;
  u[0] = k;
  for (i = 1; i <= p; i++)
    d[i] = i;
  for (i = 1; i <= u[0]; i++)
    if (u[i] != 1)
      for (j = 1; j <= p; j++)
	if (d[j] != 1)
	  {
	    x = gcd (u[i], d[j]);
	    u[i] /= x;
	    d[j] /= x;
	  }
  for (i = 1; i <= u[0]; i++)
    if (u[i] != 1)
      pu *= u[i];
  return pu;
}

int
insert_tree_node (btreeptr prev, btreeptr node, int tag)
{
  if (prev && node)
    {
      if (tag == 0)
	prev->left = node;
      else
	prev->right = node;
      node->previous = prev;
      return 1;
    }
  else
    {
      printf ("insert tree node error\n");
      return 0;
    }
}

/* generate spectrum tree, heavy use of recursion */
btreeptr
g_tree (btreeptr leaf, btreeptr prev, paramptr pptr, int m, int n, int tag)
{
  int k, l;
  k = m - 1, l = n - 1;
  if (pptr->min_s_tokens < 1 || pptr->min_s_tokens > pptr->win_width)
    {
      printf ("spectrum structure error\n");
      return NULL;
    }
  //generating the first node, root whose ->previous points to NULL
  if (!prev)
    {
      btrear = leaf;
      btreeptr root = (btreeptr) malloc (sizeof (struct btree));
      //root->data.nlay equals 0
      root->data.nlay = 0;
      root->data.nzero = 0;
      root->data.npos = 1;
      root->previous = NULL;
      root->left = NULL;
      root->right = NULL;
      if (pptr->min_s_tokens == 1)
	{
	  //in this situation root is a leaf node
	  leaf->right = root;
	  btrear = root;
	  return root;
	}
      if (m == n)
	root->right = g_tree (btrear, root, pptr, k, l, 1);
      else if (k > 0)
	{
	  //root->data.npos always equals 1
	  root->left = g_tree (btrear, root, pptr, k, l, 0);
	  root->right = g_tree (btrear, root, pptr, k, l, 1);
	}
      return root;
    }
  else
    {
      btreeptr node = (btreeptr) malloc (sizeof (struct btree));
      node->previous = NULL;
      node->left = NULL;
      node->right = NULL;
      //left tree or right tree?
      if (!tag)
	{
	  node->data.nlay = pptr->win_width - n;
	  node->data.nzero = prev->data.nzero + 1;
	  node->data.npos = prev->data.npos;
	}
      else
	{
	  node->data.nlay = pptr->win_width - n;
	  node->data.nzero = 0;
	  node->data.npos = prev->data.npos + 1;
	}
      //insert the newly created node under the prev node
      insert_tree_node (prev, node, tag);
      //spectrum density constraint
      if (node->data.nzero == (pptr->win_width - pptr->min_s_tokens)
	  || m == n)
	{
	  //some problems here
	  if (node->data.npos == pptr->min_s_tokens)
	    {
	      leaf->right = node;
	      btrear = node;
	      return node;
	    }
	  node->right = g_tree (btrear, node, pptr, k, l, 1);
	  return node;
	}
      else
	{
	  if (node->data.npos == pptr->min_s_tokens)
	    {
	      leaf->right = node;
	      btrear = node;
	      return node;
	    }
	  else
	    {
	      if (pptr->win_width - node->data.nlay - 1 <=
		  pptr->min_s_tokens - node->data.npos)
		node->left = NULL;
	      else
		node->left = g_tree (btrear, node, pptr, k + 1, l, 0);
	      node->right = g_tree (btrear, node, pptr, k, l, 1);
	      return node;
	    }
	}
    }
}

/* produce spectrum from spectrum tree */
int
get_spectrum (btreeptr leaf, int count, spcptr_list spec)
{
  int i, j, k, m, p, q, t, x;
  k = leaf->data.nlay + 1;
  m = ceil ((double) k / (double) UINT_SIZE);
  spcptr_list node[m];
  for (j = 0; j < m; j++)
    {
      node[j] = (spcptr_list) malloc (sizeof (struct spc_list));
      node[j]->data.spseg = 0;
      node[j]->previous = NULL;
      node[j]->next = NULL;
      insert_spc_node (spec, node[j]);
    }
  node[0]->data.spseg |= 1;
  //point to the first node in the spectrum
  spcptr_list cur = node[0];
  //trace up from the leaf to the root
  btreeptr tmpn = leaf->previous;
  while (tmpn)
    {
      //distance from the node to the bottom plus 1
      t = k - tmpn->data.nlay;
      //target spc node to be processed
      p = ceil ((double) t / (double) UINT_SIZE);
      //the exponent corresponding to this position
      q = (t - 1) % UINT_SIZE;
      //look for the target spc node
      while (cur->data.id != p && cur != spec)
	cur = cur->next;
      if (cur == spec)
	printf ("no such spc node\n");
      if (!tmpn->data.nzero)
	{
	  x = pow (2, q);
	  //set the position q+1 to 1
	  cur->data.spseg |= x;
	}
      //go up one layer
      tmpn = tmpn->previous;
    }
  return 1;
}

int
delete_spc_tree (btreeptr spc_tree)
{
  btreeptr node = spc_tree;
  if (!node)
    return 1;
  if (node->left)
    {
      delete_spc_tree (node->left);
      if (node->previous)
	node->previous->right = NULL;
    }
  if (node->right)
    {
      delete_spc_tree (node->right);
      if (node->previous)
	node->previous->left = NULL;
    }
  node->previous = NULL;
  free (node);
  return 1;
}

int
delete_branch (btreeptr leaf)
{
  btreeptr cur = leaf, prev;
  cur->right = NULL;
  prev = cur;
  cur = cur->previous;
  if (!cur)
    delete_spc_tree (prev);
  else
    {
      while (cur && (!cur->left || !cur->right))
	{
	  prev = cur;
	  cur = cur->previous;
	}
      delete_spc_tree (prev);
    }
  return 1;
}

/* get the region size within which the pattern can extend to the left */
int
getlx (paramptr pp, pattlist p)
{
  spcptr_list cur = p->data.spectrum->previous;
  unsigned int s = cur->data.spseg;
  unsigned int z;
  int i = 0, j = 1, t, k = 0, l = UINT_SIZE, n;
  do
    {
      n = p->data.size - (j % l);
      //can't be greater than pow(2,32)-1
      z = pow (2, n % l);
      t = s & z;
      if (t)
	i++;
      if (pp->min_s_tokens == 1 || i == pp->min_s_tokens - 1)
	{
	  if (pp->min_s_tokens == 1)
	    j = 0;
	  k = pp->win_width - j;
	  poslist curpos;
	  curpos = p->data.location->next;
	  while (curpos)
	    {
	      k = min (k, curpos->data.offset);
	      if (k == 0)
		break;
	      curpos = curpos->next;
	    }
	  return k;
	}
      j++;
      if (n % l == 0 && cur != p->data.spectrum)
	{
	  cur = cur->previous;
	  s = cur->data.spseg;
	}
    }
  while (j < p->data.size);
}

/* get the region size within which the pattern can extend to the right */
int
getrx (paramptr pp, pattlist p, t_sequence **seqs)
{
  spcptr_list cur = p->data.spectrum->next;
  unsigned int s = cur->data.spseg;
  int i = 0, j = 0, t, m, n, k, l, x, y = UINT_SIZE;
  l = pp->min_s_tokens - 1;
  unsigned int z;
  do
    {
      z = pow (2, j % y);
      j++;
      if (j >= y && j % y == 0)
	{
	  cur = cur->next;
	  s = cur->data.spseg;
	}
      t = s & z;
      if (t)
	i++;
      if (pp->min_s_tokens == 1 || i == l)
	{
	  if (pp->min_s_tokens == 1)
	    j = 0;
	  n = pp->win_width - j;
	  poslist curpos;
	  curpos = p->data.location->next;
	  x =
	    seqs[curpos->data.seqnum - 1]->length -
	    curpos->data.offset - p->data.size;
	  while (curpos)
	    {
	      k =
		seqs[curpos->data.seqnum - 1]->length - curpos->data.offset -
		p->data.size;
	      if (n <= k)
		return n;
	      else
		{
		  m = max (x, k);
		  x = m;
		}
	      curpos = curpos->next;
	    }
	  return m;
	}
    }
  while (j < p->data.size);
}

int pop(unsigned x) {
  x = x - ((x >> 1) & 0x55555555);
  x = (x & 0x33333333) + ((x >> 2) & 0x33333333);
  x = (x + (x >> 4)) & 0x0F0F0F0F;
  x = x + (x >> 8);
  x = x + (x >> 16);
  return x & 0x0000003F;
}

/* get the number of gaps in pattern */
int
gap_num (pattlist p)
{
  spcptr_list cur = p->data.spectrum->next;
  unsigned int s = cur->data.spseg;
  int i, j = 0, k, f = 0, m, n = UINT_SIZE;
  while (cur->next != p->data.spectrum)
    {
      j += pop (~s);
      cur = cur->next;
      s = cur->data.spseg;
    }
  j += pop (~s);
  return j - get_blank (s);
}

/* get the token num */
int
ngap_num (pattlist p)
{
  return (p->data.size - gap_num (p));
}

/* is compositional maximality? */
int
is_cmax (int gp, pattlist p, t_sequence **seqs)
{
  char c, d;
  poslist l = p->data.location->next;
  int offset = l->data.offset + gp;
  c = toupper ((seqs[l->data.seqnum - 1]->seq)[offset]);
  l = l->next;
  while (l)
    {
      offset = l->data.offset + gp;
      d = toupper ((seqs[l->data.seqnum - 1]->seq)[offset]);
      if (d != c)
	return 1;
      l = l->next;
    }
  return 0;
}

/* is left maximality? */
int
is_lmax (int lx, pattlist p, t_sequence **seqs)
{
  char c, d;
  int offset, i, f;
  int n = count_pos (p);
  poslist l;
  if (lx == 0)
    return 1;
  for (i = 1; i <= lx; i++)
    {
      f = 1;
      l = p->data.location->next;
      offset = l->data.offset - i;
      c = toupper ((seqs[l->data.seqnum - 1]->seq)[offset]);
      l = l->next;
      while (l)
	{
	  offset = l->data.offset - i;
	  d = toupper ((seqs[l->data.seqnum - 1]->seq)[offset]);
	  if (d != c)
	    break;
	  else
	    f++;
	  l = l->next;
	}
      if (f == n)
	return 0;
    }
  return 1;
}

/* get gap position in spectrum */
void
get_gap (int gp[], spcptr_list spc)
{
  spcptr_list cur = spc->previous;
  unsigned int s = cur->data.spseg, h = cur->data.spseg;
  int i, j = 0, k, f = 0, g = 0, m, n = UINT_SIZE, blank = get_blank (h);
  for (i = n - blank - 1; i >= 0; i--)
    {
      m = pow (2, i);
      k = s & m;
      if (k != 0)
	f = 1;
      if (f && !k)
	gp[g++] = j;
      j++;
    }
  cur = cur->previous;
  s = cur->data.spseg;
  while (cur != spc)
    {
      for (i = n - 1; i >= 0; i--)
	{
	  m = pow (2, i);
	  k = s & m;
	  if (!k)
	    gp[g++] = j;
	  j++;
	}
      cur = cur->previous;
      s = cur->data.spseg;
    }
}

/* get token position information into ngp[] */
void
get_ngap (int ngp[], spcptr_list spc)
{
  spcptr_list cur = spc->previous;
  unsigned int s = cur->data.spseg, h = cur->data.spseg;
  int i, n = UINT_SIZE, k, f = 0, g = 0, m, j = 0, blank = get_blank (h);
  for (i = n - blank - 1; i >= 0; i--)
    {
      m = pow (2, i);
      k = cur->data.spseg & m;
      if (k != 0)
	{
	  //the left most position in the spectrum, begin j counting
	  f = 1;
	  ngp[g++] = j;
	}
      if (f == 1)
	j++;
    }
  cur = cur->previous;
  s = cur->data.spseg;
  while (cur != spc)
    {
      for (i = n - 1; i >= 0; i--)
	{
	  m = pow (2, i);
	  k = cur->data.spseg & m;
	  if (k != 0)
	    ngp[g++] = j;
	  j++;
	}
      cur = cur->previous;
      s = cur->data.spseg;
    }
}

/* count the number of positions in pattern */
int
count_pos (pattlist p)
{
  int n = 0;
  poslist cur = p->data.location->next;
  while (cur)
    {
      n++;
      cur = cur->next;
    }
  return n;
}

/* is lc-maximal? */
int
is_lcmaximal (paramptr pp, pattlist p, t_sequence **seqs)
{
  int i, j, k;
  int lx = getlx (pp, p);
  int gap = gap_num (p);
  if (gap)
    {
      int gp[gap];
      get_gap (gp, p->data.spectrum);
      for (i = 0; i < gap; i++)
	if (!is_cmax (gp[i], p, seqs))
	  return 0;
    }
  if (!is_lmax (lx, p, seqs))
    return 0;
  return 1;
}

/* see if the footprints are the same */
int
is_equal_spc (spcptr_list spc1, spcptr_list spc2)
{
  spcptr_list cur1, cur2;
  unsigned int s1, s2;
  cur1 = spc1->next;
  cur2 = spc2->next;
  if (spc1->data.spseg != spc2->data.spseg)
    return 0;
  while (cur1 != spc1 && cur2 != spc2)
    {
      if (cur1->data.spseg != cur2->data.spseg)
	return 0;
      cur1 = cur1->next;
      cur2 = cur2->next;
    }
  return 1;
}

/* insert spectrum node into the double link list */
int
insert_spc_node (spcptr_list spc, spcptr_list pspc)
{
  if (!spc->next)
    {
      spc->next = pspc;
      spc->previous = pspc;
      pspc->previous = spc;
      pspc->next = spc;
    }
  else
    {
      spc->previous->next = pspc;
      pspc->previous = spc->previous;
      spc->previous = pspc;
      pspc->next = spc;
    }
  pspc->data.id = pspc->previous->data.id + 1;
  spc->data.spseg++;
  return 1;
}

void
print_spectrum (FILE * fsw, spcptr_list spc)
{
  spcptr_list cur = spc->previous;
  while (cur != spc)
    {
      fprintf (fsw, "%08x", cur->data.spseg);
      cur = cur->previous;
    }
  fprintf (fsw, "\n");
}

int
delete_spc_node (spcptr_list pspc)
{
  pspc->previous->next = pspc->next;
  pspc->next->previous = pspc->previous;
  free (pspc);
  return 1;
}

int
delete_spectrum (spcptr_list spc)
{
  if (!spc)
    return 0;
  spcptr_list prev, cur;
  prev = spc->next;
  cur = prev->next;
  while (cur && cur != spc)
    {
      delete_spc_node (prev);
      prev = cur;
      cur = cur->next;
    }
  delete_spc_node (prev);
  delete_spc_node (cur);
  return 1;
}

int
insert_position (poslist pl, poslist rear, poslist p)
{
  if (!pl->next)
    pl->next = p;
  else
    rear->next = p;
  return 1;
}

void
print_position_list (FILE * fsw, poslist pl)
{
  poslist list = pl;
  if (!list)
    return;
  int i = 0;
  poslist p1;
  p1 = list->next;
  if (!p1)
    return;
  while (p1)
    {
      i++;
      fprintf (fsw, "[%d,%d]", p1->data.seqnum, p1->data.offset);
      p1 = p1->next;
    }
  fprintf (fsw, "\n");
}

int
delete_position_list (poslist pl)
{
  poslist list = pl, prev, cur;
  if (!list)
    return 0;
  prev = list;
  cur = list->next;
  while (cur)
    {
      prev->next = NULL;
      free (prev);
      prev = cur;
      cur = cur->next;
    }
  free (prev);
  return 1;
}

/* insert a pattern into pattern list */
int
insert_pattern_rear (pattlist pl, pattlist r, pattlist p, int f)
{
  //f=1 means pl is result set
  if (f)
    p->data.id = patt_fix_uid++;
  pattlist prev, cur;
  if (!r)
    r = pl;
  prev = r;
  r->next = p;
  return 1;
}

/* insert p2 after p1 */
int
insert_pattern_after (pattlist p1, pattlist p2)
{
  pattlist tmp, prev = p1, cur = p2;
  tmp = prev->next;
  prev->next = cur;
  cur->next = tmp;
  return 1;
}

/* delete pattern by id */
int
delete_pattern (pattlist pl, pattlist p, int f)
{
  int n = 0;
  pattlist cur, prev;
  cur = prev = pl;
  if (!p || !cur)
    return 0;
  if (!f)
    {
      pl->next = p->next;
      free (p->data.text);
      delete_spectrum (p->data.spectrum);
      delete_position_list (p->data.location);
      free (p);
      return 1;
    }
  else
    {
      while (cur)
	{
	  if (!strcmp (cur->data.text, p->data.text)
	      && is_equal_spc (cur->data.spectrum, p->data.spectrum))
	    {
	      prev->next = cur->next;
	      free (cur->data.text);
	      delete_spectrum (cur->data.spectrum);
	      delete_position_list (cur->data.location);
	      free (cur);
	      return 1;
	    }
	  prev = cur;
	  cur = cur->next;
	  n++;
	}
      printf ("No such pattern! Can't delete this pattern!\n");
      return 0;
    }
}

/* target the first token in 5 steps */
int
get_blank (unsigned int spseg)
{
  unsigned y;
  int n;
  n = 32;
  y = spseg >>16;  if (y != 0) {n = n -16;  spseg = y;}
  y = spseg >> 8;  if (y != 0) {n = n - 8;  spseg = y;}
  y = spseg >> 4;  if (y != 0) {n = n - 4;  spseg = y;}
  y = spseg >> 2;  if (y != 0) {n = n - 2;  spseg = y;}
  y = spseg >> 1;  if (y != 0) return n - 2;
  return n - spseg;
}

void
print_text (FILE * fsw, pattlist p)
{
  spcptr_list cur = p->data.spectrum->previous;
  unsigned int s = cur->data.spseg;
  int i, j = 0, m = UINT_SIZE, n;
  int k = get_blank (s);
  for (i = m - k - 1; i >= 0; i--)
    {
      n = pow (2, i);
      if (!(cur->data.spseg & n))
	fprintf (fsw, ".");
      else
	fprintf (fsw, "%c", p->data.text[j++]);
    }
  cur = cur->previous;
  s = cur->data.spseg;
  while (cur != p->data.spectrum)
    {
      for (i = m - 1; i >= 0; i--)
	{
	  n = pow (2, i);
	  if (!(cur->data.spseg & n))
	    fprintf (fsw, ".");
	  else
	    fprintf (fsw, "%c", p->data.text[j++]);
	}
      cur = cur->previous;
      s = cur->data.spseg;
    }
  fprintf (fsw, " ");
}

void
print_pattern_list (FILE * fsw, pattlist pl, paramptr pptr)
{
  pattlist p1;
  poslist pos;
  if (!pl)
    return;
  pos = pl->data.location;
  p1 = pl->next;
  while (p1 && (pptr->curr_pattern < pptr->max_patterns))
    {
      fprintf (fsw, "id %d support %d size %d ", patt_dis_uid++,
	       p1->data.support, p1->data.size);
      print_text (fsw, p1);
      print_spectrum (fsw, p1->data.spectrum);
      print_position_list (fsw, p1->data.location);
      p1 = p1->next;
      pptr->curr_pattern++;
    }
}

int
delete_pattern_list (pattlist pl)
{
  pattlist prev, cur, list = pl;
  if (!list)
    return 0;
  prev = list;
  cur = list->next;
  while (cur)
    {
      free (prev->data.text);
      delete_spectrum (prev->data.spectrum);
      delete_position_list (prev->data.location);
      free (prev);
      prev = cur;
      cur = cur->next;
    }
  free (prev->data.text);
  delete_spectrum (prev->data.spectrum);
  delete_position_list (prev->data.location);
  free (prev);
  return 1;
}

pattlist
create_pattern_list (char *s)
{
  pattlist pl = (pattlist) malloc (sizeof (struct pattern_list));
  pl->data.text =
    (unsigned char *) malloc ((strlen (s) + 1) * sizeof (unsigned char));
  strcpy (pl->data.text, s);
  pl->data.location = NULL;
  pl->data.rear = NULL;
  pl->data.spectrum = NULL;
  pl->next = NULL;
  return pl;
}

/* compute pattern size from spectrum */
int
get_size (spcptr_list spc)
{
  spcptr_list cur = spc->previous;
  unsigned int s = cur->data.spseg;
  int i, j, left, n = UINT_SIZE;
  left = 31 - get_blank (s);
  return left + 1 + n * (spc->data.spseg - 1);
}

spcptr_list
clone_spc_node (spcptr_list node)
{
  spcptr_list clone = (spcptr_list) malloc (sizeof (struct spc_list));
  clone->data.id = node->data.id;
  clone->data.spseg = node->data.spseg;
  clone->previous = node->previous;
  clone->next = node->previous;
  return clone;
}

spcptr_list
clone_spec (spcptr_list spc)
{
  spcptr_list clone, cur;
  int i, nodenum = spc->data.spseg;;
  clone = (spcptr_list) malloc (sizeof (struct spc_list));
  clone->data.id = spc->data.id;
  clone->data.spseg = 0;
  clone->previous = NULL;
  clone->next = NULL;
  cur = spc->next;
  while (cur != spc)
    {
      spcptr_list node = (spcptr_list) malloc (sizeof (struct spc_list));
      node->data.id = cur->data.id;
      node->data.spseg = cur->data.spseg;
      node->previous = NULL;
      node->next = NULL;
      insert_spc_node (clone, clone_spc_node (node));
      cur = cur->next;
      free (node);
    }
  return clone;
}

/* just clone this node */
poslist
clone_pos (poslist pl)
{
  poslist clone = (poslist) malloc (sizeof (struct position_list));
  clone->data.seqnum = pl->data.seqnum;
  clone->data.offset = pl->data.offset;
  clone->next = NULL;
  return clone;
}

void
clone_patt (pattlist tmp_patt, pattlist curpatt)
{
  int newlen = curpatt->data.size + 1;
  tmp_patt->data.id = patt_tmp_uid++;
  if (tmp_patt->data.text)
    free (tmp_patt->data.text);
  tmp_patt->data.text =
    (unsigned char *) malloc ((newlen + 1) * sizeof (unsigned char));
  strcpy (tmp_patt->data.text, curpatt->data.text);
  tmp_patt->data.text[newlen - 1] = '\0';
  tmp_patt->data.spectrum = clone_spec (curpatt->data.spectrum);
  tmp_patt->data.size = curpatt->data.size;
  //support and location have no need to be cloned here
}

pattlist
exactly_clone_patt (pattlist pl)
{
  poslist cur = pl->data.location->next;
  pattlist clone = (pattlist) malloc (sizeof (struct pattern_list));
  clone->data.id = pl->data.id;
  clone->data.text =
    (unsigned char *) malloc ((strlen (pl->data.text) + 1) *
			      sizeof (unsigned char));
  strcpy (clone->data.text, pl->data.text);
  clone->data.support = pl->data.support;
  clone->data.spectrum = clone_spec (pl->data.spectrum);
  clone->data.size = pl->data.size;
  clone->data.location = (poslist) malloc (sizeof (struct position_list));
  clone->data.location->next = NULL;
  clone->data.rear = cur;

  while (cur)
    {
      poslist pos = clone_pos (cur);
      insert_position (clone->data.location, clone->data.rear, pos);
      clone->data.rear = pos;
      cur = cur->next;
    }
  clone->next = NULL;
  return clone;
}

/* get the value of the block that will be left-shifted */
unsigned int
get_flow_block (unsigned int spseg, int rx)
{
  int i, j;
  unsigned int k, l = 0, n;
  for (i = UINT_SIZE - 1; i >= UINT_SIZE - rx - 1; i--)
    {
      k = pow (2, i);
      if (k & spseg)
	{
	  j = i - UINT_SIZE + rx;
	  n = pow (2, j);
	  l = l | n;
	}
    }
  return l;
}

/* update spectrum information */
void
teleport (spcptr_list pspc, int rx, int n)
{
  int i;
  unsigned int p, q, x, y;
  spcptr_list flo;
  for (i = n - 1; i >= 0; i--)
    {
      flo = pspc;
      p = pow (2, i);
      q =
	ceil ((double) ((pspc->data.id - 1) * UINT_SIZE + i + 1 + rx) /
	      (double) UINT_SIZE);
      while (flo->data.id != q)
	{
	  flo = flo->next;
	}
      // target offset in a pspc
      x = ((pspc->data.id - 1) * UINT_SIZE + i + rx) % UINT_SIZE;
      y = pow (2, x);
      if (p & pspc->data.spseg)
	flo->data.spseg |= y;
      else
	flo->data.spseg &= ~y;
    }
}

/* extend the pattern to a longer one and update the corresponding pattern information */
void
update_patt (pattlist p, char c, int rx)
{
  spcptr_list prev, cur, start, ostart;
  int i, j, m, n = UINT_SIZE, blank, x;
  unsigned k, l;
  char ec[2] = { c };
  strcat (p->data.text, ec);
  blank = get_blank (p->data.spectrum->previous->data.spseg);
  x = rx - blank;
  //calculate the number of nodes that should be added to the extended spectrum
  if (x <= 0)
    m = 0;
  else
    m = ceil ((double) (rx - blank) / (double) n);
  spcptr_list tmp[m];
  //the left-most node in the orignal spectrum
  ostart = p->data.spectrum->previous;
  for (i = 0; i < m; i++)
    {
      tmp[i] = (spcptr_list) malloc (sizeof (struct spc_list));
      tmp[i]->data.spseg = 0;
      tmp[i]->previous = NULL;
      tmp[i]->next = NULL;
      insert_spc_node (p->data.spectrum, clone_spc_node (tmp[i]));
      free (tmp[i]);
    }
  l = p->data.spectrum->data.spseg;
  //the left-most node in the extended spectrum
  start = prev = p->data.spectrum->previous;
  if (!(rx % n))
    {
      cur = ostart;
      while (cur != p->data.spectrum)
	{
	  prev->data.spseg = cur->data.spseg;
	  cur->data.spseg = 0;
	  prev = prev->previous;
	  cur = cur->previous;
	}
      p->data.spectrum->next->data.spseg = 1;
    }
  else if (rx < n && m == 0)
    {
      cur = prev->previous;
      // m = 0 means there is enough space for left shift
      prev->data.spseg <<= rx;
      while (cur != p->data.spectrum)
	{
	  k = get_flow_block (cur->data.spseg, rx);
	  prev->data.spseg |= k;
	  cur->data.spseg <<= rx;
	  prev = cur;
	  cur = cur->previous;
	}
      p->data.spectrum->next->data.spseg += 1;
    }
  else
    {
      cur = ostart;
      int mask;
      //to transport partial data from one node to the next
      teleport (cur, rx, n - get_blank (cur->data.spseg));
      //cur->data.spseg = (rx < n ? cur->data.spseg << rx: 0);

      if (rx < n)
	{
	  mask = pow (2, rx) - 1;
	  cur->data.spseg &= ~mask;
	}
      else
	cur->data.spseg = 0;

      cur = cur->previous;
      while (cur != p->data.spectrum)
	{
	  teleport (cur, rx, n);
	  //cur->data.spseg = (rx < n ? cur->data.spseg << rx: 0);

	  if (rx < n)
	    {
	      mask = pow (2, rx) - 1;
	      cur->data.spseg &= ~mask;
	    }
	  else
	    cur->data.spseg = 0;

	  prev = cur;
	  cur = cur->previous;
	}
      p->data.spectrum->next->data.spseg += 1;
    }
  p->data.size += rx;
}

/* create a new seed pattern */
pattlist
mak_patt (char *d, spcptr_list spc, int support, int size, poslist pos)
{
  pattlist pl = (pattlist) malloc (sizeof (struct pattern_list));
  pl->data.id = patt_tmp_uid++;
  pl->data.text =
    (unsigned char *) malloc ((strlen (d) + 1) * sizeof (unsigned char));
  strcpy (pl->data.text, d);
  pl->data.spectrum = clone_spec (spc);
  pl->data.support = support;
  pl->data.size = size;
  pl->data.location = (poslist) malloc (sizeof (struct position_list));
  pl->data.location->next = pos;
  pl->data.rear = pos;
  pl->next = NULL;
  return pl;
}

/* see whether there already exists this seed, if yes, update seed information, if no, create it */
void
check_and_update (char *d, spcptr_list spc, poslist pos, pattlist seed_sub,
		  pattlist rear_sub, paramptr pptr)
{
  pattlist p, cur;
  if (!seed_sub->next)
    rear_sub->next = NULL;
  cur = seed_sub->next;
  while (cur)
    {
      if (!memcmp (cur->data.text, d, pptr->min_s_tokens))
	{
	  cur->data.support++;
	  if (need_pos)
	    {
	      insert_position (cur->data.location, cur->data.rear, pos);
	      cur->data.rear = pos;
	    }
	  else
	    free (pos);
	  return;
	}
      cur = cur->next;
    }
  int size = get_size (spc);
  if (!need_pos && pos)
    {
      free (pos);
      pos = NULL;
    }
  p = mak_patt (d, spc, 1, size, pos);
  insert_pattern_rear (seed_sub, rear_sub->next, p, 0);
  rear_sub->next = p;
}

/* heuristic estimation of division number */
int
estimate_divs (paramptr p)
{
  int divs, l, alphn;
  l = p->min_s_tokens;
  //make sure that divs is not too big
  if (p->type == 0 && l > 5)
    {
      l = 5;
      if (no_ext)
	{
	  if (l == 6)
	    l = 6;
	  else
	    l = 7;
	}
    }
  if (p->type == 1 && l > 4)
    l = 4;
  if (p->type == 2 && l > 3)
    l = 3;
  alphn = alph_num (p->type);
  divs = pow (alphn, l);
  return divs;
}

/* generate seed pattern */
void
generate_seed (t_sequence **seqs, paramptr pptr, spcptr_list spec, pattlist seed)
{
  int i, j, k, l, m, size, offset, h;
  size_t bytes_read = 0;
  char c;
  char d[pptr->min_s_tokens + 1];
  int ngp[pptr->min_s_tokens];	/* indicate the token positions in the spectrum */
  /* fill the ngp array with the token positions in the spectrum */
  int divs = estimate_divs (pptr);
  /* the head array and the rear array */
  pattlist prefix_hash[divs], seed_sub;
  pattlist nexus_hash[divs], rear_sub;
  // n = the number of layers u wanna go
  int n = pptr->min_s_tokens;

  // change the number if the memory is not enough
  if (pptr->type == 0 && n > 5)
    {
      n = 5;
      if (no_ext)
	{
	  if (n == 6)
	    n = 6;
	  else
	    n = 7;
	}
    }
  if (pptr->type == 1 && n > 4)
    n = 4;
  if (pptr->type == 2 && n > 3)
    n = 3;
  for (i = 0; i < divs; i++)
    {
      prefix_hash[i] = create_pattern_list ("sub_seed");
      nexus_hash[i] = create_pattern_list ("sub_rear");
    }
  get_ngap (ngp, spec);

  /* for every sequence */
  for (j = 0; j < seqn; j++)
    {
      /* the sequence region that should be scanned */
      for (k = 0; k <= seqs[j]->length - get_size (spec); k++)
	{
	  poslist tmp_pos = (poslist) malloc (sizeof (struct position_list));
	  tmp_pos->data.seqnum = j + 1;
	  tmp_pos->data.offset = k;
	  tmp_pos->next = NULL;
	  for (l = 0; l < pptr->min_s_tokens; l++)
	    {
	      /* k offset in file, ngp offset in pattern */
	      offset = k + ngp[l];
	      d[l] = toupper ((seqs[j]->seq)[offset]);
	    }
	  d[l] = '\0';

	  h = get_hash_key (d, n, alph_num (pptr->type), 1);
	  seed_sub = prefix_hash[h];
	  rear_sub = nexus_hash[h];
	  /* check if the seed is already in the seed_set, if not, add it into the set, otherwise, update the support */
	  check_and_update (d, spec, tmp_pos, seed_sub, rear_sub, pptr);
	}
    }
  int sign = 0;
  pattlist nexus = NULL;
  for (i = 0; i < divs; i++)
    {
      if (prefix_hash[i]->next && sign == 0)
	{
	  seed->next = prefix_hash[i]->next;
	  nexus = nexus_hash[i]->next;
	  sign = 1;
	  continue;
	}
      if (nexus_hash[i]->next && sign == 1)
	{
	  nexus->next = prefix_hash[i]->next;
	  nexus = nexus_hash[i]->next;
	}
    }
  for (i = 0; i < divs; i++)
    {
      free (prefix_hash[i]->data.text);
      free (prefix_hash[i]);
      free (nexus_hash[i]->data.text);
      free (nexus_hash[i]);
    }
}

void go_into_result (paramptr pptr, FILE *fsw, pattlist rearpatt, pattlist rrear, pattlist result_set)
{
  if (rearpatt->data.support >= pptr->support4pattern
      && rearpatt->data.size >= pptr->patt_size
      && strlen (rearpatt->data.text) >= pptr->min_p_tokens)
    {
      if (!need_patt)
	{
	  if (pptr->curr_pattern < pptr->max_patterns) {
          pthread_mutex_lock (&fsw_mutex);
    	  fprintf (fsw, "id %d support %d size %d ", patt_dis_uid++, rearpatt->data.support, rearpatt->data.size);
    	  print_text (fsw, rearpatt);
    	  print_spectrum (fsw, rearpatt->data.spectrum);
    	  print_position_list (fsw, rearpatt->data.location);
    	  pptr->curr_pattern++;
          pthread_mutex_unlock (&fsw_mutex);
      }

      delete_spectrum (rearpatt->data.spectrum);
	  delete_position_list (rearpatt->data.location);
	  free (rearpatt->data.text);
	  free (rearpatt);
	  rearpatt = NULL;
	}
      else
	{
	  pthread_mutex_lock (&patt_fix_uid_mutex);
	  insert_pattern_rear (result_set, rrear, rearpatt, 1);
	  pthread_mutex_unlock (&patt_fix_uid_mutex);
	}
    }
  else
    {
      delete_spectrum (rearpatt->data.spectrum);
      delete_position_list (rearpatt->data.location);
      free (rearpatt->data.text);
      free (rearpatt);
      rearpatt = NULL;
    }
}

/* $$$ core function $$$ */
/* it feeds on seeds and try its best to extend these seeds */
void
find_patt (t_sequence **seqs, FILE * fsw, paramptr pptr,
	   pattlist seed, pattlist result_set)
{
  poslist curpos;
  char x, y;
  size_t bytes_read = 0;
  int rx, i, j, k, l, outrx, ox, result = 0, aln = alph_num (pptr->type);
  pattlist pl, curpatt, nexp, rearpatt, rrear = NULL, tmp_patt[aln];

  curpatt = seed->next;
  if (!curpatt)
    return;

  /* for every tmp pattern in seed */
  while (curpatt)
    {
      outrx = 0;
      /* calculate the rx for pattern extention */
      rx = getrx (pptr, curpatt, seqs);
      if (rx == 0)
	{
	  rearpatt = exactly_clone_patt (curpatt);
	  go_into_result (pptr, fsw, rearpatt, rrear, result_set);
	  if (rearpatt)
	    rrear = rearpatt;
	  nexp = curpatt->next;
	  delete_pattern (seed, curpatt, 1);
	}
      else
	{
	  for (i = 1; i <= rx; i++)
	    {
	      /* for every position of the current pattern, record the number of
	         the appearance of each character at the given position */
	      curpos = curpatt->data.location->next;

	      /* initialize tmp_patt array to produce new extended pattern */
	      for (l = 0; l < aln; l++)
		{
		  tmp_patt[l] =
		    (pattlist) malloc (sizeof (struct pattern_list));
		  tmp_patt[l]->data.id = tmp_patt[l]->data.size =
		    tmp_patt[l]->data.support = 0;
		  tmp_patt[l]->data.text = NULL;
		  tmp_patt[l]->data.spectrum = NULL;
		  tmp_patt[l]->data.location =
		    (poslist) malloc (sizeof (struct position_list));
		  tmp_patt[l]->data.rear = curpos;
		  tmp_patt[l]->data.location->next = NULL;
		  tmp_patt[l]->next = NULL;
		}

	      while (curpos)
		{
		  j =
		    curpos->data.offset +
		    curpatt->data.size + i - 1;
		  if (j <= seqs[curpos->data.seqnum - 1]->length - 1)
		    {
		      x = toupper ((seqs[curpos->data.seqnum - 1]->seq)[j]);
		      /* record the occurence of the letter at this position */
		      /* this is an very important function. be careful!!! */
		      ox = 0;
		      if (pptr->type == 0)
			ox = get_hash_key (&x, 1, alph_num (0), 0);
		      /* allocate proper size to tmp_patt array according to pptr->type */
		      /* can NOT insert curpos, we must clone a curpos , or the original curpos
		         will get lost from original poslist!!! */
		      else if (pptr->type == 1)
			ox = get_hash_key (&x, 1, alph_num (1), 0);
		      else if (pptr->type == 2)
			ox = get_hash_key (&x, 1, alph_num (2), 0);
		      tmp_patt[ox]->data.support++;
		      poslist rpos = clone_pos (curpos);
		      insert_position (tmp_patt[ox]->data.location,
				       tmp_patt[ox]->data.rear, rpos);
		      tmp_patt[ox]->data.rear = rpos;
		    }
		  //try next position of this pattern
		  curpos = curpos->next;
		}		//end of while curpos
	      /* keep those derived patterns which have support equal to or greater than pptr->support4pattern */
	      for (k = 0; k < aln; k++)
		{
		  pattlist wp = tmp_patt[k];
		  /* tmp_patt[k]->data.support must be less than curpatt->data.support in this case */
		  if (tmp_patt[k]->data.support >= pptr->support4pattern)
		    {
		      if (tmp_patt[k]->data.support == curpatt->data.support)
			{
			  /* fill the content of the temporary pattern with the information of the current pattern
			     and update some related fields */
			  clone_patt (tmp_patt[k], curpatt);
			  update_patt (tmp_patt[k], get_char (k, pptr->type),
				       i);
			  //put this temporary pattern after the current pattern
			  if (is_lcmaximal
			      (pptr, tmp_patt[k], seqs))
			    insert_pattern_after (curpatt,
						  exactly_clone_patt (tmp_patt
								      [k]));
			  nexp = curpatt->next;
			  //delete the curpatt from seed, we have better ones. to search the curpatt in the seed
			  //is very expensive. some algorithms should be adopted
			  delete_pattern (seed, curpatt, 1);
			  /* there wont be a second pattern at this position with such support.
			     jump out the loop */
			  outrx = 1;
			  for (j = k; j < aln; j++)
			    {
			      if (tmp_patt[j])
				{
				  delete_spectrum (tmp_patt[j]->data.
						   spectrum);
				  delete_position_list (tmp_patt[j]->data.
							location);
				  free (tmp_patt[j]->data.text);
				  free (tmp_patt[j]);
				}
			    }
			  break;
			}
		      //a candidate extended pattern, fill and update the content of this temporary pattern
		      clone_patt (tmp_patt[k], curpatt);
		      update_patt (tmp_patt[k], get_char (k, pptr->type), i);
		      //is this extended pattern lc-maximal?
		      if (is_lcmaximal
			  (pptr, tmp_patt[k], seqs))
			insert_pattern_after (curpatt,
					      exactly_clone_patt (tmp_patt
								  [k]));

		    }
		  //those temporary patterns that dont meet the criteria should be removed.
		  delete_spectrum (tmp_patt[k]->data.spectrum);
		  delete_position_list (tmp_patt[k]->data.location);
		  free (tmp_patt[k]->data.text);
		  free (tmp_patt[k]);
		  if (i == rx && k == aln - 1)
		    {
		      //the current pattern is a maximal pattern that meets the criteria. it shoud be added
		      //into the result-set and removed from the seed
		      rearpatt = exactly_clone_patt (curpatt);
		      go_into_result (pptr, fsw, rearpatt, rrear, result_set);
		      if (rearpatt)
			rrear = rearpatt;
		      nexp = curpatt->next;
		      //as mentioned above, this is an expensive operation
		      delete_pattern (seed, curpatt, 1);
		    }
		}		//end of for aln
	      if (outrx == 1)
		break;
	    }			//end of for rx
	}			//end of else
      curpatt = nexp;
    }				//end of while curpatt
}

void
insert_job_queue (struct phase1_jobs *job)
{
  pthread_mutex_lock (&p1_jobs_mutex);
  job->next = p1_jobs;
  p1_jobs = job;
  pthread_mutex_unlock (&p1_jobs_mutex);
}

/* package of tasks for each thread */
int
process_phase1_jobs (struct phase1_jobs *next_job, struct pthread_param *p)
{
  pattlist result_set = create_pattern_list ("result");
  pattlist seed = create_pattern_list ("seed");
  generate_seed (p->seqs, p->pptr, next_job->job, seed);
  delete_spectrum (next_job->job);
  //to lcmaximalize seed
  pattlist curpatt = seed->next;
  pattlist prev = seed;
  if (!no_ext)
    {
      while (curpatt)
	{
	  if (curpatt->data.support < p->pptr->support4seed
	      || !is_lcmaximal (p->pptr, curpatt, p->seqs))
	    {
	      delete_pattern (prev, curpatt, 0);
	      curpatt = prev->next;
	    }
	  else
	    {
	      prev = curpatt;
	      curpatt = curpatt->next;
	    }
	}
      find_patt (p->seqs, p->fsw, p->pptr, seed, result_set);

      free (seed->data.text);
      free (seed);

      if (need_patt)
	{
	  pthread_mutex_lock (&fsw_mutex);
	  print_pattern_list (p->fsw, result_set, p->pptr);
	  pthread_mutex_unlock (&fsw_mutex);
	}
    }
  else
    {
      while (curpatt)
	{
	  if (curpatt->data.support < p->pptr->support4seed)
	    {
	      delete_pattern (prev, curpatt, 0);
	      curpatt = prev->next;
	    }
	  else
	    {
	      prev = curpatt;
	      curpatt = curpatt->next;
	    }
	}
      pthread_mutex_lock (&fsw_mutex);
      print_pattern_list (p->fsw, seed, p->pptr);
      pthread_mutex_unlock (&fsw_mutex);

      delete_pattern_list (seed);
    }

  delete_pattern_list (result_set);
  return 1;
}

/* thread function */
void *
thread_find (void *param)
{
  struct pthread_param *p = (struct pthread_param *) param;
  while (1)
    {
      struct phase1_jobs *next_job;
      pthread_mutex_lock (&p1_jobs_mutex);
      if (!p1_jobs)
	next_job = NULL;
      else
	{
	  next_job = p1_jobs;
	  p1_jobs = p1_jobs->next;
	  job_count++;
	}
      pthread_mutex_unlock (&p1_jobs_mutex);
      if (!next_job)
	break;
      process_phase1_jobs (next_job, p);
      pthread_mutex_lock (&patt_fix_uid_mutex);
      if (!no_ext)
	{
	  if (need_patt)
	    printf ("Job No.: %d   Patt No.: %d\r", job_count,
		    patt_fix_uid - 1);
	  else
	    printf ("Job No.: %d   Patt No.: %d\r", job_count,
		    patt_dis_uid - 1);
	}
      fflush (stdout);
      pthread_mutex_unlock (&patt_fix_uid_mutex);

      free (next_job);
    }
  return NULL;
}
