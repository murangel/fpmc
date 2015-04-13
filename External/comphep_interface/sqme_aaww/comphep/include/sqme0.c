/*
   * Copyright (C) 2001-2003, CompHEP COllaboration
   *------------------------------------------------------
   * $Id: sqme0.c,v 1.4 2009/04/23 08:05:33 okepka Exp $
   *
   * $Log: sqme0.c,v $
   * Revision 1.4  2009/04/23 08:05:33  okepka
   * Fixing files in sqme_aaww
   *
   * Revision 1.2  2009/04/22 13:54:58  okepka
   * Namespaces in comphep, anomalous ZZ, change in energy cut to select survival p.
   *
   * Revision 1.1  2008/12/15 15:58:39  okepka
   * comphep
   *
   * Revision 1.1.1.1  2008-12-13 19:25:18  olda
   * Start
   *
   * Revision 1.7  2003/06/25 18:46:41  kryukov
   * *** empty log message ***
   *
   * Revision 1.6  2003/04/22 09:01:18  kryukov
   * Apply indent to improve readability of files
   *
   * Revision 1.5  2003/04/13 20:41:06  kryukov
   * Bug in comment
   *
   * Revision 1.4  2003/04/06 11:39:45  kryukov
   * Missprint in cvs tag
 */

#include<stdio.h>

int *calcCoef = NULL;
double *Q0 = NULL;
double *Q1 = NULL;
double *Q2 = NULL;
double *color_weights = NULL;

double computer_eps;
double Fmax;

static int **c_perm = NULL;
static double *cw_buff = NULL;
static int *particls = NULL;
static int cBasisPower;


static int 
indx_ (int k, int l)
{
  int i, j;
  if (k < l)
    {
      i = k;
      j = l;
    }
  else
    {
      i = l;
      j = k;
    }
  return i + (j * (j - 1)) / 2;
}


static void 
sprod_ (double *momenta)
{
  int k, i, j, ntot;

  ntot = nin_ + nout_;
  for (i = 0; i < ntot - 1; ++i)
    {
      for (j = i + 1; j < ntot; ++j)
	{
	  double *sum = DP + indx_ (i, j);
	  double *v1 = momenta + 4 * (i);
	  double *v2 = momenta + 4 * (j);
	  *sum = *v1 ** v2;
	  for (k = 1; k <= 3; k++)
	    (*sum) -= v1[k] * v2[k];
	}
    }
}				/* sprod_ */


double 
sqrMom (char *momnum, double *lv)
{
  char *ii;
  double s[4] =
  {0, 0, 0, 0};
  char nin_char;
  nin_char = nin_;

  ii = momnum;
  while (*ii)
    {
      int k;
      if (*ii > nin_char)
	for (k = 0; k < 4; k++)
	  s[k] -= lv[k + 4 * (*ii - 1)];
      else
	for (k = 0; k < 4; k++)
	  s[k] += lv[k + 4 * (*ii - 1)];
      ii++;
    }
  return (s[0] - s[1]) * (s[0] + s[1]) - s[2] * s[2] - s[3] * s[3];
}

static double 
simSqme (int nsub, double *momenta, int ntot, int level, int *err)
{
  int n, i, k;
  double ans, buff;

  if (level == ntot)
    {
      buff = smpl (nsub, momenta, err);
      calcCoef[nsub] = 0;
      return buff;
    }

  ans = simSqme (nsub, momenta, ntot, level + 1, err);
  n = 1;

  for (i = level + 1; i <= ntot; i++)
    {
      int *p = c_perm[indx_ (i - nin_ - 1, level - nin_ - 1)];

      if (particls[i] == particls[level])
	{
	  for (k = 0; k < 4; k++)
	    {
	      buff = momenta[k + 4 * level - 4];
	      momenta[k + 4 * level - 4] = momenta[k + 4 * i - 4];
	      momenta[k + 4 * i - 4] = buff;
	    }
	  if (color_weights && p)
	    {
	      int m;
	      memcpy (cw_buff, color_weights, sizeof (double) * cBasisPower);
	      for (m = 0; m < cBasisPower; m++)
		color_weights[m] = cw_buff[p[m]];
	    }

	  ans += simSqme (nsub, momenta, ntot, level + 1, err);
	  for (k = 0; k < 4; k++)
	    {
	      buff = momenta[k + 4 * level - 4];
	      momenta[k + 4 * level - 4] = momenta[k + 4 * i - 4];
	      momenta[k + 4 * i - 4] = buff;
	    }
	  if (color_weights && p)
	    {
	      int m;
	      memcpy (cw_buff, color_weights, sizeof (double) * cBasisPower);
	      for (m = 0; m < cBasisPower; m++)
		color_weights[p[m]] = cw_buff[m];
	    }
	  n++;
	}
    }
  return ans / n;
}


static void 
initSqme (int nsub)
{
  int nC, *cChains;
  int i, j;
  int ntot;
  char name_i[10], name_j[10];

  cStrings (nsub, &nC, &cBasisPower, &cChains);

  if (cw_buff)
    free (cw_buff);
  cw_buff = (double *) malloc (sizeof (double) * cBasisPower);

  for (i = 0; i < (nout_ * (nout_ - 1)) / 2; i++)
    {
      if (c_perm[i])
	{
	  free (c_perm[i]);
	  c_perm[i] = NULL;
	}
    }

  for (i = 0; i < nout_ - 1; i++)
    for (j = i + 1; j < nout_; j++)
      {
	pinf_ (nsub, i + nin_ + 1, name_i, NULL);
	pinf_ (nsub, j + nin_ + 1, name_j, NULL);
	if (!strcmp (name_i, name_j))
	  {
	    int k, l, l2;

	    int *cChains_ = (int*)malloc (2 * sizeof (int) * nC * cBasisPower);
	    memcpy (cChains_, cChains, 2 * sizeof (int) * nC * cBasisPower);

	    c_perm[indx_ (i, j)] = (int*)malloc (sizeof (int) * cBasisPower);

	    for (k = 0; k < 2 * nC * cBasisPower; k++)
	      if (cChains_[k] == i + nin_ + 1)
		cChains_[k] = j + nin_ + 1;
	      else if (cChains_[k] == j + nin_ + 1)
		cChains_[k] = i + nin_ + 1;

	    for (l = 0; l < cBasisPower; l++)
	      {
		int *c = cChains_ + 2 * nC * l;
		k = 0;
		while (k < 2 * nC - 2)
		  {
		    if (c[k] > c[k + 2])
		      {
			int c2 = c[k];
			int c3 = c[k + 1];
			c[k] = c[k + 2];
			c[k + 1] = c[k + 3];
			c[k + 2] = c2;
			c[k + 3] = c3;
			if (k > 0)
			  k -= 2;
			else
			  k += 2;
		      }
		    else
		      k += 2;
		  }

		for (l2 = 0; l2 < cBasisPower; l2++)
		  {
		    int *cc = cChains + 2 * nC * l2;
		    for (k = 0; k < 2 * nC; k++)
		      if (c[k] != cc[k])
			break;
		    if (k == 2 * nC)
		      {
			c_perm[indx_ (i, j)][l] = l2;
			break;
		      }
		  }
		if (l2 == cBasisPower)
		  printf ("Can not construct permutation\n");
	      }
	    free (cChains_);
	  }
      }

  ntot = nin_ + nout_;

  if (particls)
    free (particls);
  particls = (int *) malloc ((ntot + 1) * sizeof (int));

  for (i = 1; i <= ntot; i++)
    particls[i] = i;

  for (i = 1; i < ntot; i++)
    {
      if (particls[i] == i)
	{
	  pinf_ (nsub, i, name_i, NULL);
	  for (j = i + 1; j <= ntot; j++)
	    {
	      if (particls[j] == j)
		{
		  pinf_ (nsub, j, name_j, NULL);
		  if (strcmp (name_i, name_j) == 0)
		    particls[j] = i;
		}
	    }
	}
    }
}


double 
sqme_ (int nsub, double *momenta, int *err)
{
  char name[10];
  int i;
  double val;
  static double *amem = NULL;
  int recalc = 0;
  static int nsub0 = 0;

  if (nsub == 0)
    {
      nsub0 = 0;
      return 0;
    }
  Fmax = 0;

  if (nsub0 == 0)
    {
      double one = 1;
      double one_plus_eps;
      computer_eps = 1;

      do
	{
	  computer_eps = computer_eps / 2;
	  one_plus_eps = one + computer_eps;
	}
      while (one_plus_eps != one);
      computer_eps *= 2;

      c_perm = (int **)malloc (sizeof (int *) * (1 + indx_ (nout_, nout_ - 1)));
      for (i = 0; i < (nout_ * (nout_ - 1)) / 2; i++)
	c_perm[i] = NULL;
    }

  if (nsub0 != nsub)
    {
      initSqme (nsub);
      nsub0 = nsub;
    }

  if (color_weights)
    for (i = 0; i < cBasisPower; i++)
      color_weights[i] = 0;

  if (!calcCoef)
    {
      if (amem)
	free (amem);
      amem = (double *) malloc (sizeof (double) * (1 + nvar_));
      calcCoef = (int *) malloc (sizeof (int) * (nprc_ + 1));
      for (i = 1; i <= nvar_; i++)
	vinf_ (i, NULL, amem + i);
      recalc = 1;
    }
  else
    {
      for (i = 1; i <= nvar_; i++)
	{
	  vinf_ (i, name, &val);
	  if (strcmp ("GG", name) && val != amem[i])
	    {
	      amem[i] = val;
	      recalc = 1;
	    }
	}
    }
  if (recalc)
    {
      *err = calcFunc ();
      if (*err)
	{
	  *err = 1;
	  return 0;
	}
      for (i = 1; i <= nprc_; i++)
	calcCoef[i] = 1;
    }

  strfun_calc = 0;
  return simSqme (nsub, momenta, nin_ + nout_, nin_ + 1, err);
}
