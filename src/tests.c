#include <math.h>
#include <time.h>
#include <string.h>
#include <stdlib.h>
#include <R.h>
#include <Rmath.h>


/* ***********************************************************************
 *                                                                       *
 *                          Declaration of functions                     *
 *                                                                       *
 * ********************************************************************* */



/* Functions coming from the package ade4 */

void vecpermut (double *A, int *num, double *B);
double alea (void);
void aleapermutvec (double *a);
void trirapideintswap (int *v, int i, int j);
void trirapideint (int *x , int *num, int gauche, int droite);
void sqrvec (double *v1);
void getpermutation (int *numero, int repet);
void prodmatABC (double **a, double **b, double **c);
void prodmatAtAB (double **a, double **b);
void prodmatAtBC (double **a, double **b, double **c);
void prodmatAAtB (double **a, double **b);
void prodmatAtBrandomC (double **a, double **b, double **c, int*permut);
void taballoc (double ***tab, int l1, int c1);
void vecalloc (double **vec, int n);
void vecintalloc (int **vec, int n);
void freetab (double **tab);
void freevec (double *vec);
void freeintvec (int *vec);
void matcentrage (double **A, double *poili, char *typ);
void matmodiffc (double **tab, double *poili);
void matmodifcp (double **tab, double *poili);
void matmodifcs (double **tab, double *poili);
void matmodifcn (double **tab, double *poili);
void matmodifcm (double **tab, double *poili);
void DiagobgComp (int n0, double **w, double *d, int *rang);





/* Functions from the package adehabitat */
void mvtfreeman(int *in, int *jn, int *dir, int *np);
void getcontourc(double *grille, int *nlig, int *ncol, int *indicelig, 
		int *indicecol, int *lcont);
void lcontour(double *grille, int *nlig, int *ncol, int *lcont);
void levels(double *vec, double *lev, int *lvec);
void seqeticorr(double *grille, int *nlig, int *ncol);
void erodil(double *grille, int *nlig, int *ncol, int *ntour, int *oper);
void inout(double *x, double *y, double *xp, double *yp,
	   int *deds);
void inoutr(double *xr, double *yr, double *xpr, double *ypr,
	    int *dedsr, int *nxr, int *npr);
void rastpol(double *xp, double *yp, double *xg, double *yg,
	     double **carte);
void rastpolaire(double *xpr, double *ypr, double *xgr, double *ygr,
		 double *carter, int *nlg, int *ncg, int *nvp);
void regroufacasc(double **asce, double **ascs, int *np,
		  int *nlev);
void regroufacascr(double *ascer, double *ascsr, int *npr,
		   int *nlevr, int *nle, int *nce, int *nls, 
		   int *ncs);
void regrouascnum(double **ascent, double **ascso);
void regrouascnumr(double *ascentr, double *ascsor, 
		   double *nler, double *ncer,
		   double *nlsr, double *ncsr);
void ligpoly(double *x, double *y, double r, double *xp, double *yp);
void buflig(double **x, double r, double **carte, double *xg, double *yg);
void bufligr(double *xr, double *rr, double *carter, 
	     double *xgr, double *ygr, int *nlr, int *ncr, 
	     int *nlocr);





/*********************************************************************
 *********************************************************************
 *********                                                       *****
 *********               The sources of ADE-4                    *****
 *********               --------------------                    *****
 *********************************************************************
 *********************************************************************
 */



/**************************/
double alea (void)
{
    double w;
    GetRNGstate();
    w = unif_rand();
    PutRNGstate();
    return (w);
}

/*************************/
void aleapermutvec (double *a)
{
    /* Randomly permutes the elements of a vector a
       Manly p. 42 The vector is modified
       from Knuth 1981 p. 139 */
    int lig, i,j, k;
    double z;
    
    lig = a[0];
    for (i=1; i<=lig-1; i++) {
	j=lig-i+1;
	k = (int) (j*alea()+1);
	/* k = (int) (j*genrand()+1); */
	if (k>j) k=j;
	z = a[j];
	a[j]=a[k];
	a[k] = z;
    }
}


/*******************/	
void vecpermut (double *A, int *num, double *B)
{
/*---------------------------------------
 * A is a vector with n elements
 * B is a vector with n elements
 * num is a random permutation of the n first integers
 * B contains in output the permuted elements of A
 * ---------------------------------------*/
    
    int lig, lig1, lig2, i, k;
    
    lig = A[0];
    lig1 = B[0];
    lig2 = num[0];
    
    
    if ( (lig!=lig1) || (lig!=lig2) ) {
	/* err_message ("Illegal parameters (vecpermut)");
	   closelisting(); */
    }
    
    for (i=1; i<=lig; i++) {
	k=num[i];
	B[i] = A[k];
    }
}

/********* Centring accrding to row weights poili **********/	
void matcentrage (double **A, double *poili, char *typ)
{
    
    if (strcmp (typ,"nc") == 0) {
	return;
    } else if (strcmp (typ,"cm") == 0) {
	matmodifcm (A, poili);
	return;
    } else if (strcmp (typ,"cn") == 0) {
	matmodifcn (A, poili);
	return;
    } else if (strcmp (typ,"cp") == 0) {
	matmodifcp (A, poili);
	return;
    } else if (strcmp (typ,"cs") == 0) {
	matmodifcs (A, poili);
	return;
    } else if (strcmp (typ,"fc") == 0) {
	matmodiffc (A, poili);
	return;
    } else if (strcmp (typ,"fl") == 0) {
	matmodifcm (A, poili);
	return;
    }
}

/*********************/
void matmodifcm (double **tab, double *poili)
/*--------------------------------------------------
 * tab is a complete disjonctive table with n rows and m columns
 * poili is a vector with n components
 * The process returns tab centred by column
 * with weighting poili (sum=1)
 * centring type multple correspondances
 --------------------------------------------------*/
{
    double		poid;
    int 			i, j, l1, m1;
    double		*poimoda;
    double		x, z;
    
    l1 = tab[0][0];
    m1 = tab[1][0];
    vecalloc(&poimoda, m1);
    
    
    for (i=1;i<=l1;i++) {
	poid = poili[i];
	for (j=1;j<=m1;j++) {
	    poimoda[j] = poimoda[j] + tab[i][j] * poid;
	}
    }
    
    for (j=1;j<=m1;j++) {
	x = poimoda[j];
	if (x==0) {
	    for (i=1;i<=l1;i++) tab[i][j] = 0;
	} else {
	    
	    for (i=1;i<=l1;i++) {
		z = tab[i][j]/x - 1.0;
		tab[i][j] = z;
	    }
	}
    }
    freevec (poimoda);
}

/*********************************************************/
void matmodifcn (double **tab, double *poili)
/*--------------------------------------------------
 * tab is a table n rows and p columns
 * poili is a vector with n components
 * the function returns tab normed by column
 * with the weighting poili (sum=1)
 --------------------------------------------------*/
{
    double		poid, x, z, y, v2;
    int 			i, j, l1, c1;
    double		*moy, *var;
    
    l1 = tab[0][0];
    c1 = tab[1][0];
    
    vecalloc(&moy, c1);
    vecalloc(&var, c1);
    
    
/*--------------------------------------------------
 * centred and normed table
 --------------------------------------------------*/
    
    for (i=1;i<=l1;i++) {
	poid = poili[i];
	for (j=1;j<=c1;j++) {
	    moy[j] = moy[j] + tab[i][j] * poid;
	}
    }
    
    for (i=1;i<=l1;i++) {
	poid=poili[i];
	for (j=1;j<=c1;j++) {
	    x = tab[i][j] - moy[j];
	    var[j] = var[j] + poid * x * x;
	}
    }
    
    for (j=1;j<=c1;j++) {
	v2 = var[j];
	if (v2<=0) v2 = 1;
	v2 = sqrt(v2);
	var[j] = v2;
    }
    
    for (i=1;i<=c1;i++) {
	x = moy[i];
	y = var[i];
	for (j=1;j<=l1;j++) {
	    z = tab[j][i] - x;
	    z = z / y;
	    tab[j][i] = z;
	}
    }
    
    freevec(moy);
    freevec(var);
    
}

/*********************************************************/
void matmodifcs (double **tab, double *poili)
/*--------------------------------------------------
 * tab is a table n rows, p columns
 * poili is a vector with n components
 * The function returns tab standardised by column
 * for the weighting poili (sum=1)
 --------------------------------------------------*/
{
	double		x,poid, z, y, v2;
	int 			i, j, l1, c1;
	double		*var;
	
	l1 = tab[0][0];
	c1 = tab[1][0];
	vecalloc(&var, c1);
	

/*--------------------------------------------------
 * calculation of the standardised table
 --------------------------------------------------*/
	
	for (i=1;i<=l1;i++) {
	    poid=poili[i];
	    for (j=1;j<=c1;j++) {
		x = tab[i][j];
		var[j] = var[j] + poid * x * x;
	    }
	}
	
	for (j=1;j<=c1;j++) {
	    v2 = var[j];
	    if (v2<=0) v2 = 1;
	    v2 = sqrt(v2);
	    var[j] = v2;
	}
	
	for (i=1;i<=c1;i++) {
	    y = var[i];
	    for (j=1;j<=l1;j++) {
		z = tab[j][i];
		z = z / y;
		tab[j][i] = z;
	    }
	}
	freevec(var);
}


/**********/
void matmodifcp (double **tab, double *poili)
/*--------------------------------------------------
 * tab is a table with n rows and p colonnes
 * poili is a vector with n components
 * The function returns tab centred by column
 * for the weighting poili (sum=1)
 --------------------------------------------------*/
{
    double		poid;
    int 			i, j, l1, c1;
    double		*moy, x, z;
    
    l1 = tab[0][0];
    c1 = tab[1][0];
    vecalloc(&moy, c1);
    
    
/*--------------------------------------------------
 * Centred table
 --------------------------------------------------*/
    
    for (i=1;i<=l1;i++) {
	poid = poili[i];
	for (j=1;j<=c1;j++) {
	    moy[j] = moy[j] + tab[i][j] * poid;
	}
    }
    
    
    for (i=1;i<=c1;i++) {
	x = moy[i];
	for (j=1;j<=l1;j++) {
	    z = tab[j][i] - x;
	    tab[j][i] = z;
	}
    }
    freevec(moy);
}

/*********************/
void matmodiffc (double **tab, double *poili)
/*--------------------------------------------------
 * tab is a table with n rows and m columns
 * of number >=0
 * poili is a vector with n components
 * The function returns tab doubly centred
 * for the weighting poili (sum=1)
 * centring type simple correspondance analysis
 --------------------------------------------------*/
{
    double		poid;
    int 			i, j, l1, m1;
    double		*poimoda;
    double		x, z;
    
    l1 = tab[0][0];
    m1 = tab[1][0];
    vecalloc(&poimoda, m1);
    
    
    for (i=1;i<=l1;i++) {
	x = 0;
	for (j=1;j<=m1;j++) {
	    x = x + tab[i][j];
	}
	if (x!=0) {
	    for (j=1;j<=m1;j++) {
		tab[i][j] = tab[i][j]/x;
	    }
	}	
    }
    
    for (i=1;i<=l1;i++) {
	poid = poili[i];
	for (j=1;j<=m1;j++) {
	    poimoda[j] = poimoda[j] + tab[i][j] * poid;
	}
    }
    
    for (j=1;j<=m1;j++) {
	x = poimoda[j];
	if (x==0) {
	    /* err_message("column has a nul weight (matmodiffc)"); */
	}
	
	for (i=1;i<=l1;i++) {
	    z = tab[i][j]/x - 1.0;
	    tab[i][j] = z;
	}
    }
    freevec (poimoda);
}









/*****************/
void getpermutation (int *numero, int repet)
/*----------------------
 * affects a random permutation of the first n integers
 * in an integer vector of length n
 * First vecintalloc is needed
 * *numero is a vector of integer
 * repet is an integer which can take any arbitrary value
 * used in the seed of the pseudo-random number generation process
 * if it is increased in repeated calls (e.g. simulation), it is ensured that
 * two calls returns different results (seed=clock+repet)
 ------------------------*/
{
    int i, n;
    int *alea;
    
    n=numero[0];
    vecintalloc (&alea,n);
    
    /*-------------
     * numbering in numero
     -----------*/
    for (i=1;i<=n;i++) {
	numero[i]=i;
    }
    
    /*-------------
     * affects random numbers in alea
     ----------------*/
    for (i=1;i<=n;i++) {
 	GetRNGstate();
	alea[i] = ((int) (1e8)*(unif_rand()));
	PutRNGstate();
    }
    
    trirapideint (alea , numero, 1, n);
    freeintvec (alea);
}

/*****************************************/
/* Sorting: used in getpermutation */

void trirapideint (int *x , int *num, int gauche, int droite)
{
    int j, dernier, milieu, t;
    
    if ( (droite-gauche)<=0) return;
    
    milieu = (gauche+droite)/2;
    trirapideintswap (x, gauche, milieu);
    trirapideintswap (num, gauche, milieu);
    
    t=x[gauche];
    dernier=gauche;
    for (j = gauche+1; j<=droite; j++) {
	if (x[j] < t) {
	    dernier = dernier + 1;
	    trirapideintswap (x, dernier, j);	
	    trirapideintswap (num, dernier, j);
	}
    }
    trirapideintswap (x, gauche, dernier);
    trirapideintswap (num, gauche, dernier);
    
    trirapideint (x, num, gauche, dernier-1);
    trirapideint (x, num, dernier+1, droite);
    
}

/**************************************/
/* Sorting: used in trirapideint */

void trirapideintswap (int *v, int i, int j)
{
    int provi;
    
    provi=v[i];
    v[i]=v[j];
    v[j]=provi;
}

/***********************************************************************/
void sqrvec (double *v1)
/*--------------------------------------------------
 * Square root of the elements of a vector
 --------------------------------------------------*/
{
    int i, c1;
    double v2;
    
    c1 = v1[0];
    
    for (i=1;i<=c1;i++) {
	v2 = v1[i];
	/* if (v2 < 0.0) err_message("Error: Square root of negative number (sqrvec)"); */
	v2 = sqrt(v2);
	v1[i] = v2;
    }
}

/***********************************************************************/
void DiagobgComp (int n0, double **w, double *d, int *rang)
/*--------------------------------------------------
 * Eigenstructure of a matrix. See
 * T. FOUCART Analyse factorielle de tableaux multiples,
 * Masson, Paris 1984,185p., p. 62. D'apr?s VPROP et TRIDI,
 * de LEBART et coll.
 --------------------------------------------------*/
{
    double			*s;
    double			a, b, c, x, xp, q, bp, ab, ep, h, t, u , v;
    double			dble;
    int				ni, i, i2, j, k, jk, ijk, ij, l, ix, m, m1, isnou;
    
    vecalloc(&s, n0);
    a = 0.000000001;
    ni = 100;
    if (n0 == 1) {
	d[1] = w[1][1];
	w[1][1] = 1.0;
	*rang = 1;
	freevec (s);
	return;
    }
    
    for (i2=2;i2<=n0;i2++) {
	
	b=0.0;
	c=0.0;
	i=n0-i2+2;
	k=i-1;
	if (k < 2) goto Et1;
	for (l=1;l<=k;l++) {
	    c = c + fabs((double) w[i][l]);
	}
	if (c != 0.0) goto Et2;
	
    Et1:	s[i] = w[i][k];
	goto Etc;
	
    Et2:	for (l=1;l<=k;l++) {
	x = w[i][l] / c;
	w[i][l] = x;
	b = b + x * x;
    }
	xp = w[i][k];
	ix = 1;
	if (xp < 0.0) ix = -1;
		
/*		q = -sqrt(b) * ix; */
	dble = b;
	dble = -sqrt(dble);
	q = dble * ix;
	
	s[i] = c * q;
	b = b - xp * q;
	w[i][k] = xp - q;
	xp = 0;
	for (m=1;m<=k;m++) {
	    w[m][i] = w[i][m] / b / c;
	    q = 0;
	    for (l=1;l<=m;l++) {
		q = q + w[m][l] * w[i][l];
	    }
	    m1 = m + 1;
	    if (k < m1) goto Et3;
	    for (l=m1;l<=k;l++) {
		q = q + w[l][m] * w[i][l];
	    }
	    
	Et3:		s[m] = q / b;
	    xp = xp + s[m] * w[i][m];
	}
	bp = xp * 0.5 / b;
	for (m=1;m<=k;m++) {
	    xp = w[i][m];
	    q = s[m] - bp * xp;
	    s[m] = q;
	    for (l=1;l<=m;l++) {
		w[m][l] = w[m][l] - xp * s[l] - q * w[i][l];
	    }
	}
	for (l=1;l<=k;l++) {
	    w[i][l] = c * w[i][l];
	}
	
    Etc:	d[i] = b;
    } /* for (i2=2;i2<n0;i2++) */
    
    s[1] = 0.0;
    d[1] = 0.0;
    
    for (i=1;i<=n0;i++) {
	
	k = i - 1;
	if (d[i] == 0.0) goto Et4;
	for (m=1;m<=k;m++) {
	    q = 0.0;
	    for (l=1;l<=k;l++) {
		q = q + w[i][l] * w[l][m];
	    }
	    for (l=1;l<=k;l++) {
		w[l][m] = w[l][m] - q * w[l][i];
	    }
	}
	
    Et4:	d[i] = w[i][i];
	w[i][i] = 1.0;
	if (k < 1) goto Et5;
	for (m=1;m<=k;m++) {
	    w[i][m] = 0.0;
	    w[m][i] = 0.0;
	}
	
    Et5:;
    }
    
    for (i=2;i<=n0;i++) {
	s[i-1] = s[i];
    }
    s[n0] = 0.0;
    
    for (k=1;k<=n0;k++) {
	
	m = 0;
	
    Et6: 	for (j=k;j<=n0;j++) {
	if (j == n0) goto Et7;
	ab = fabs((double) s[j]);
	ep = a * (fabs((double) d[j]) + fabs((double) d[j+1]));
	if (ab < ep) goto Et7;
    }
	
    Et7: 	isnou = 1;
	h = d[k];
	if (j == k) goto Eta;
	if (m < ni) goto Etd;
	
	/* err_message("Error: can't compute matrix eigenvalues"); */
	
    Etd:	m = m + 1;
	q = (d[k+1]-h) * 0.5 / s[k];
	
/*		t = sqrt(q * q + 1.0); */
	dble = q * q + 1.0;
	dble = sqrt(dble);
	t = dble;
	
	if (q < 0.0) isnou = -1;
	q = d[j] - h + s[k] / (q + t * isnou);
	u = 1.0;
	v = 1.0;
	h = 0.0;
	jk = j-k;
	for (ijk=1;ijk<=jk;ijk++) {
	    i = j - ijk;
	    xp = u * s[i];
	    b = v * s[i];
	    if (fabs((double) xp) < fabs((double) q)) goto Et8;
	    u = xp / q;
	    
/*			t = sqrt(u * u + 1); */
	    dble = u * u + 1.0;
	    dble = sqrt(dble);
	    t = dble;
	    
	    s[i+1] = q * t;
	    v = 1 / t;
	    u = u * v;
	    goto Et9;
	    
	Et8:		v = q / xp;
	    
/*			t = sqrt(1 + v * v); */
	    dble = 1.0 + v * v;
	    dble = sqrt(dble);
	    t = dble;
	    
	    s[i+1] = t * xp;
	    u = 1 / t;
	    v = v * u;
	    
	Et9:
	    q = d[i+1] - h;
	    t = (d[i] - q) * u + 2.0 * v * b;
	    h = u * t;
	    d[i+1] = q + h;
	    q = v * t - b;
	    for (l=1;l<=n0;l++) {
		xp = w[l][i+1];
		w[l][i+1] = u * w[l][i] + v * xp;
		w[l][i] = v * w[l][i] - u * xp;
	    }
	}
	d[k] = d[k] - h;
	s[k] = q;
	s[j] = 0.0;
	
	goto Et6;
	
    Eta:;
    } /* for (k=1;k<=n0;k++) */
    
    for (ij=2;ij<=n0;ij++) {
	
	i = ij - 1;
	l = i;
	h = d[i];
	for (m=ij;m<=n0;m++) {
	    if (d[m] >= h) {
		l = m;
		h = d[m];
	    }
	}
	if (l == i) {
	    goto Etb;
	} else {
	    d[l] = d[i];
	    d[i] = h;
	}
	for (m=1;m<=n0;m++) {
	    h = w[m][i];
	    w[m][i] = w[m][l];
	    w[m][l] = h;
	}
	
    Etb:;
    } /* for (ij=2;ij<=n0;ij++) */
    
    /* final:; */
    *rang = 0;
    for (i=1;i<=n0;i++) {
	/*
	  if (d[i] / d[1] < 0.00001) d[i] = 0.0;
	  if (d[i] != 0.0) *rang = *rang + 1;
	*/
	if (d[i] > 0.0) *rang = *rang + 1;
    }
    freevec(s);
} /* DiagoCompbg */







/***********************************************************************/
void prodmatABC (double **a, double **b, double **c)
/*--------------------------------------------------
* Matrix product AB
--------------------------------------------------*/
{
    int j, k, i, lig, col, col2;
    double s;
    
    lig = a[0][0];
    col = a[1][0];
    
    col2 = b[1][0];
    
    for (i=1;i<=lig;i++) {
	for (k=1;k<=col2;k++) {
	    s = 0;
	    for (j=1;j<=col;j++) {
		s = s + a[i][j] * b[j][k];
	    }
	    c[i][k] = s;
	}		
    }
}

/***********************************************************************/
void prodmatAtAB (double **a, double **b)
/*--------------------------------------------------
* Matrix product AtA
--------------------------------------------------*/
{
    int j, k, i, lig, col;
    double s;
    
    lig = a[0][0];
    col = a[1][0];
    
    for (j=1;j<=col;j++) {
	for (k=j;k<=col;k++) {
	    s = 0;
	    for (i=1;i<=lig;i++) {
		s = s + a[i][k] * a[i][j];
	    }
	    b[j][k] = s;
	    b[k][j] = s;
	}		
    }
}

/***********************************************************************/
void prodmatAtBC (double **a, double **b, double **c)
/*--------------------------------------------------
 * Matrix product AtB
 --------------------------------------------------*/
{
    int j, k, i, lig, col, col2;
    double s;
    
    lig = a[0][0];
    col = a[1][0];
    
    col2 = b[1][0];
    
    for (j=1;j<=col;j++) {
	for (k=1;k<=col2;k++) {
	    s = 0;
	    for (i=1;i<=lig;i++) {
		s = s + a[i][j] * b[i][k];
	    }
	    c[j][k] = s;
	}		
    }
}


/***********************************************************************/
void prodmatAAtB (double **a, double **b)
/*--------------------------------------------------
 * Matrix product B = AAt
 --------------------------------------------------*/
{
    int j, k, i, lig, col;
    double s;
    
    lig = a[0][0];
    col = a[1][0];
    
    for (j=1;j<=lig;j++) {
	for (k=j;k<=lig;k++) {
	    s = 0;
	    for (i=1;i<=col;i++) {
		s = s + a[j][i] * a[k][i];
	    }
	    b[j][k] = s;
	    b[k][j] = s;
	}		
    }
}

/*******************/
void prodmatAtBrandomC (double **a, double **b, double **c, int*permut)
/*--------------------------------------------------
 * Produit matriciel AtB
 * les lignes de B sont permutees par la permutation permut
 --------------------------------------------------*/
{
    int j, k, i, i0, lig, col, col2;
    double s;
    
    lig = a[0][0];
    col = a[1][0];
    
    col2 = b[1][0];
    
    for (j=1;j<=col;j++) {
	for (k=1;k<=col2;k++) {
	    s = 0;
	    for (i=1;i<=lig;i++) {
		i0 = permut[i];
		s = s + a[i][j] * b[i0][k];
	    }
	    c[j][k] = s;
	}		
    }
}

/***********************************************************************/
void taballoc (double ***tab, int l1, int c1)
/*--------------------------------------------------
 * Dynamic Memory Allocation for a table (l1, c1)
 --------------------------------------------------*/
{
    int i, j;
    
    if ( (*tab = (double **) calloc(l1+1, sizeof(double *))) != 0) {
	for (i=0;i<=l1;i++) {
	    if ( (*(*tab+i)=(double *) calloc(c1+1, sizeof(double))) == 0 ) {
		return;
		for (j=0;j<i;j++) {
		    free(*(*tab+j));
		}
	    }
	}
    }
    
    **(*tab) = l1;
    **(*tab+1) = c1;
}

/***********************************************************************/
void vecalloc (double **vec, int n)
/*--------------------------------------------------
 * Memory Allocation for a vector of length n
 --------------------------------------------------*/
{
    if ( (*vec = (double *) calloc(n+1, sizeof(double))) != 0) {
	**vec = n;
	return;
    } else {
	return;
    }
}

/*****************/
void vecintalloc (int **vec, int n)
/*--------------------------------------------------
 * Memory allocation for an integer vector of length  n
 --------------------------------------------------*/
{
    if ( (*vec = (int *) calloc(n+1, sizeof(int))) != NULL) {
	**vec = n;
	return;
    } else {
	return;
    }
}

/***********************************************************************/
void freetab (double **tab)
/*--------------------------------------------------
 * Free memory for a table
 --------------------------------------------------*/
{
    int 	i, n;
    
    n = *(*(tab));
    for (i=0;i<=n;i++) {
	free((char *) *(tab+i) );
    }
    free((char *) tab);
}

/***********************************************************************/
void freevec (double *vec)
/*--------------------------------------------------
 * Free memory for a vector
 --------------------------------------------------*/
{
    free((char *) vec);	
}

/***********************************************************************/
void freeintvec (int *vec)
/*--------------------------------------------------
* Free memory for an integer  vector
--------------------------------------------------*/
{
    
    free((char *) vec);
    
}














/*********************************************************************
 *********************************************************************
 *********                                                       *****
 *********               The sources of adehabitatMA             *****
 *********               ---------------------------             *****
 *********************************************************************
 *********************************************************************
 */




/* ****************************************************************
   *                                                              *
   * mvtfreeman: arguments = indices of the rows and columns      *
   * (in et jn), of the freeman direction (dir), and we get       *
   * the indices of rows and columns (in the vector np) after the *
   * move.                                                        *
   *                                                              *
   **************************************************************** */

void mvtfreeman(int *in, int *jn, int *dir, int *np)
{
    int i,j;
    i=*in;
    j=*jn;
    
    if ((*dir == 0) | (*dir == 1) | (*dir == 7)) 
	i++;
    if ((*dir == 3) | (*dir == 4) | (*dir == 5)) 
	i--;
    if ((*dir == 1) | (*dir == 2) | (*dir == 3)) 
	j++;
    if ((*dir == 5) | (*dir == 6) | (*dir == 7)) 
	j--;
    
    np[1]=i;
    np[2]=j;
}


/* ****************************************************************
   *                                                              *
   * algorithm of contour monitoring (suivi de contour) to get    *
   * contour polygon                                              *
   *                                                              *
   **************************************************************** */


void getcontourc(double *grille, int *nlig, int *ncol, int *indicelig, 
		int *indicecol, int *lcont)
{
    /* Declaration of local variables */
    int i, j, k, nl, nc, *idlig, *idcol, *P0, *P1, fini, *np, dirprec, dir;
    int lidlig;
    double **x;
    
    /* Memory allocation */
    nl=*nlig;
    nc=*ncol;
    vecintalloc(&P0,2);
    vecintalloc(&P1,2);
    vecintalloc(&np,2);
    taballoc(&x, nl,nc);
    vecintalloc(&idlig, *lcont);
    vecintalloc(&idcol, *lcont);
    
    /* Copy from R -> C variables */
    k=0;
    for (i=1; i<=nl; i++) {
	for(j=1; j<=nc; j++) {
	    x[i][j] = grille[k];
	    k++;
	}
    }
    
    /* Search the indices of the rows and columns
       of the first cell with a non-missing value */
    k=0;
    i=0;
    j=1;
    
    while (k==0) {
	if (i != nl) {
	    i++;
	}
	else {
	    i=1;
	    j++;
	}
	k = (int) x[i][j];
    }
    
    /* When it is found, the algorithm begins */
    idlig[1] = i;
    idcol[1] = j;
    lidlig = 1;
    P0[1] = i;
    P0[2] = j;
    dir = 4;
    
    fini = 0;
    k = 0;
  
    while (fini==0) {
	
	/* finds the next direction */
	while (k==0) {
	    dir = (dir + 1)%8;
	    mvtfreeman(&i, &j, &dir, np);
	    dirprec = (dir + 5)%8;
	    k = (int) x[np[1]][np[2]];
	}
	/* once found, stores the new coordinate */
	if (lidlig == 1) {
	    P1[1] = np[1];
	    P1[2] = np[2];
	}
	else {
	    /* P0 is the first point of the contour and P1, the last
	       Is the contour closed? */
	    if ((i==P0[1])&&(j==P0[2])&&(np[1]==P1[1])&&(np[2]==P1[2])) 
		fini =1;
	}
	
	/* If it is not, then stores the result, and perform the move
	 in the found direction */
	if (fini==0) {
	    lidlig++;
	    idlig[lidlig] = np[1];
	    idcol[lidlig] = np[2];
	    i = np[1];
	    j = np[2];
	    mvtfreeman(&i, &j, &dirprec, np);
	    k = (int) x[np[1]][np[2]];
	    dir = dirprec;
	}
    }
    
    /* Copy from C -> R variables */
    for (i=1; i<=lidlig; i++) {
	indicelig[i-1]=idlig[i];
	indicecol[i-1]=idcol[i];
    }
    
    /* Free Memory */
    freeintvec(idlig);
    freeintvec(idcol);
    freeintvec(P0);
    freeintvec(P1);
    freeintvec(np);
    freetab(x);
}



/* ****************************************************************
   *                                                              *
   * Interactive version of getcontour for R                      *
   *                                                              *
   **************************************************************** */

void lcontour(double *grille, int *nlig, int *ncol, int *lcont)
{
    /* Declaration of local variables*/
    int i, j, k, nl, nc, *P0, *P1, fini, *np, dirprec, dir;
    int lidlig;
    double **x, **vois;
    
    /* Memory allocation */
    nl=*nlig;
    nc=*ncol;
    vecintalloc(&P0,2);
    vecintalloc(&P1,2);
    vecintalloc(&np,2);
    taballoc(&vois, 3,3);
    taballoc(&x, nl,nc);
    
    /* R objects -> C objects */
    k=0;
    for (i=1; i<=nl; i++) {
	for(j=1; j<=nc; j++) {
	    x[i][j] = grille[k];
	    k++;
	}
    }
    
    
    /* Search of the first cell for which the value is not NA */
    k=0;
    i=0;
    j=1;
    
    while (k==0) {
	if (i != nl) {
	    i++;
	}
	else {
	    i=1;
	    j++;
	}
	k = (int) x[i][j];
    }
    
    
    /* When found, performs the algorithm */
    lidlig = 1;
    P0[1] = i;
    P0[2] = j;
    dir = 4;
        
    fini = 0;
    k = 0;
    
    /* Same algorithm as in the previous function */
    
    while (fini==0) {
	while (k==0) {
	    dir = (dir + 1)%8;
	    mvtfreeman(&i, &j, &dir, np);
	    dirprec = (dir + 5)%8;
	    k = (int) x[np[1]][np[2]];
	}
	if (lidlig == 1) {
	    P1[1] = np[1];
	    P1[2] = np[2];
	}
	else {
	    if ((i==P0[1])&&(j==P0[2])&&(np[1]==P1[1])&&(np[2]==P1[2])) 
		fini = 1;
	}
	
	if (fini==0) {
	    lidlig++;
	    i = np[1];
	    j = np[2];
	    mvtfreeman(&i, &j, &dirprec, np);
	    k = (int) x[np[1]][np[2]];
	    dir = dirprec;
	}
    }
    
    /* and return to R */
    *lcont = lidlig;

    /* Free Memory */
    freeintvec(P0);
    freeintvec(P1);
    freeintvec(np);
    freetab(vois);
    freetab(x);
}




/* ****************************************************************
   *                                                              *
   *                 Gets the levels of a factor                  *
   *                                                              *
   **************************************************************** */


void levels(double *vec, double *lev, int *lvec)
{
    /* Declaration of local variables */
    int i,j,k,n, l;
    lev[1] = vec[1];
    k=1;
    n=*lvec;
  
    /* gets the levels */
    for (i=2; i<=n; i++) {
	l=0;
	for (j=1; j<=k; j++) {
	    if (fabs(vec[i] - lev[j]) < 0.000000001)
		l=1;
	}
	if (l==0) {
	    k++;
	    lev[k] = vec[i];
	}
    }
    *lvec = k;
}


/* ****************************************************************
   *                                                              *
   *             Sequential labelling of connex components        *
   *                                                              *
   **************************************************************** */


void seqeticorr(double *grille, int *nlig, int *ncol)
{
    /* Declaration of local variables */
    int i, j, k, l, m, n, o, nl, nc, pr, beta, nniv, eticour;
    double **x, *Tc, *prec, *tmp, *tmp1, *tmp2, *etcons, *lf;
    
    /* Memory allocation */
    nl=*nlig;
    nc=*ncol;
    taballoc(&x, nl, nc);
    vecalloc(&Tc, nl*nc);
    
    /* R objects -> C objects */
    k=0;
    for (i=1; i<=nl; i++) {
	for (j=1; j<=nc; j++) {
	    x[i][j]=grille[k];
	    k++;
	}
    }
    
    Tc[1]=1;
    eticour=1;
    
    for (j=2; j<=nc; j++) {
	for (i=2; i<=nl; i++) {
	    if (((int) x[i][j])!=0) {
		vecalloc(&prec, 4);
		prec[1] = x[i-1][j-1];
		prec[2] = x[i][j-1];
		prec[3] = x[i+1][j-1];
		prec[4] = x[i-1][j];
		
		k=0;
		for (l=1; l<=4; l++) {
		    if (((int) prec[l])!=0)
			k++;
		}
		
		/* k contains the number of non null predecessors */
		if (k!=0) {
		    vecalloc(&tmp, k); /* tmp contains the non null pred */
		    m=1;
		    for (l=1; l<=4; l++) {
			if (((int) prec[l])>0) {
			    tmp[m] = prec[l];
			    m++;
			}
		    }
		    
		    freevec(prec);
		    vecalloc(&prec, k);
		    for (l=1; l<=k; l++)
			prec[l] = tmp[l];
		    freevec(tmp);
		    /* Now, prec contains the non null preds */
		    
		    
		    
		    /* Number of levels of the factor prec */
		    vecalloc(&tmp1, 4);
		    m=k;
		    levels(prec, tmp1, &m);
		    /* m contains the number of levels */
		    vecalloc(&tmp2, m);
		    /* tmp2 contains the levels of prec
		       (equivalent of etiprec in R) */
		    for (l=1; l<=m; l++)
			tmp2[l]=tmp1[l];
		    freevec(tmp1);
		    
		    if (m == 1) {
			x[i][j] = tmp2[1];
		    } else {
			/* computation of the minimum 
			   level and storage in xij */
			x[i][j] = tmp2[1];
			for (l = 1; l <= m; l++) {
			    if (tmp2[l]<x[i][j])
				x[i][j] = tmp2[l];
			}
			
			/* etcons will contain the different labels of 
			   xij */
			vecalloc(&etcons, m-1);
			n=1;
			for (l=1; l<=m; l++) {
			    if (fabs(x[i][j] - tmp2[l]) > 0.000000001) {
				etcons[n] = tmp2[l];
				n++;
			    }
			}
			
			/* loop to fill the correspondence table */
			for (l=1; l<=(m-1); l++) {
			    pr = (int) etcons[l];
			    beta = pr;
			    while (((int) Tc[beta])!=beta) {
				o = (int) Tc[beta];
				Tc[beta] = Tc[(int) x[i][j]];
				beta = o;
			    }
			    Tc[beta] = Tc[(int) x[i][j]];
			}
			freevec(prec);
			freevec(tmp2);
			freevec(etcons);
		    }
		} else {
		    Tc[eticour] = eticour;
		    x[i][j]= eticour;
		    eticour++;
		}
	    }
	}
    }
    
    eticour--;
    
    /* Actualisation of the correspondence table */
    for (i=1; i<=eticour; i++) {
	j = i;
	while (((int) Tc[j])!=j)
	    j = (int) Tc[j];
	Tc[i] = j;
    }
    j=eticour;
    vecalloc(&tmp1, j);
    vecalloc(&tmp2, eticour);
    for (i=1; i<=eticour; i++) {
	tmp2[i]=Tc[i];
    }
    
    levels(tmp2, tmp1, &j);
    freevec(tmp2);
    
    vecalloc(&lf, j);
    for (i=1; i<=j; i++)
	lf[i]=tmp1[i];
    freevec(tmp1);
    nniv=j;
    
    /* Second pass */
    for (i=1; i<=nl; i++) {
	for (j=1; j<=nc; j++) {
	    if (fabs(x[i][j]) > 0.000000001) {
		x[i][j] = Tc[(int) x[i][j]];
	    }
	}
    }

    /* Last pass: levels varying from 1 to p */
    k = 1;
    for (j=1; j<=nniv; j++) {
      i = (int) lf[j];
      if (i != k) {
	for (l = 1; l <= nl; l++) {
	  for (m = 1; m <= nc; m++) {
	    if (((int) x[l][m]) == i)
	      x[l][m]=k;
	  }
	}
      }
      k++;
    }

    /* grid */
    k=0;
    for (i=1; i<=nl; i++) {
      for (j=1; j<=nc; j++) {
	grille[k]=x[i][j];
	k++;
      }
    }

    freetab(x);
    freevec(Tc);
    freevec(lf);
  }





/* ****************************************************************
   *                                                              *
   *   Morphological dilatation and erosion                       *
   *                                                              *
   **************************************************************** */

void erodil(double *grille, int *nlig, int *ncol, int *ntour, int *oper)
{
    /* declaration */
    int i,j,k,l,nl,nc, nt, etat0, etat1;
    double **x, **xm, *voisin;
    
    nl = *nlig;
    nc = *ncol;
    nt = *ntour;
    etat0 = 0;
    etat1 = 0;
    
    /* Memory alloocation */
    taballoc(&x,nl,nc);
    taballoc(&xm,nl,nc);
    vecalloc(&voisin, 9);
    
    /* R to C */
    k=0;
    for (i=1; i<=nl; i++) {
	for (j=1; j<=nc; j++) {
	    x[i][j]=grille[k];
	    k++;
	}
    }
    
    /* Morphology */
    for (k=1; k<=nt; k++) {
	for (i=2; i<= (nl-1); i++) {
	    for (j=2; j<= (nc-1); j++) {
		voisin[1] = x[i-1][j-1];
		voisin[2] = x[i-1][j];
		voisin[3] = x[i-1][j+1];
		voisin[4] = x[i][j-1];
		voisin[5] = x[i][j+1];
		voisin[6] = x[i+1][j-1];
		voisin[7] = x[i+1][j];
		voisin[8] = x[i+1][j+1];
		voisin[9] = x[i][j];
		for (l=1;l<=9; l++) {
		    if (((int) voisin[l])==0) {
			etat0 = etat0 + 1;
		    } else {
			etat1 = etat1 + 1;
		    }
		}
		if (*oper==1) {
		    if (etat1 > 0)
			xm[i][j] = 1;
		    if (etat1 == 0)
			xm[i][j] = 0;
		} else {
		    if (etat0 == 0)
			xm[i][j] =1;
		    if (etat0 > 0)
			xm[i][j] =0;
		}
		etat1 = 0;
		etat0 = 0;
	    }
	}
    
	
	for (i=1; i<=nl; i++) {
	    for (j=1; j<=nc; j++) {
		x[i][j]=xm[i][j];
	    }
	}
    }
    
    
    /* C to R */
    k=0;
    for (i=1; i<=nl; i++) {
	for (j=1; j<=nc; j++) {
	    grille[k]=xm[i][j];
	    k++;
	}
    }
    
    /* Memory */
    freetab(x);
    freetab(xm);
    freevec(voisin);
    
}










/********************************************************
 x and y are the coordinates of the points, xp and yp 
 the coordinates of the vertices of the polygon (first and
 last vertices should be the same). deds is a vector of the
 same length as x and y: deds take the value 1 if the point is
 in the polygon and 0 otherwise
********************************************************/

void inout(double *x, double *y, double *xp, double *yp,
	   int *deds)
{
    /* Declaration of variables */
    int i, j, n, wm, np;
    double *xpc, *ypc, sig, a, b, x0;
    
    /* Memory allocation */
    n = x[0];
    np = xp[0];
    
    vecalloc(&xpc, np);
    vecalloc(&ypc, np);
    
    for (i = 1; i <= n; i++) {
	deds[i] = 1;
    }
    
    for (j = 1; j <= n; j++) {
	
	/* Centring on the point */
	for (i = 1; i <= np; i++) {
	    xpc[i] = xp[i] - x[j];
	    ypc[i] = yp[i] - y[j];
	}
	
	/* Number of intersections with X axis, for X >0 */
	wm = 0;
	for (i = 1; i <= (np-1); i++) {
	    sig = ypc[i] * ypc[i+1];
	    if (sig < 0) {
		/* The slope and intercept */
		/* Case 1: The slope is not infinite */
		if (fabs(xpc[i+1] - xpc[i]) > 0.000000001)
		{
		    a = (ypc[i+1] - ypc[i]) / (xpc[i+1] - xpc[i]);
		    b = (ypc[i]- a * xpc[i]);
		    /* value of x for y = 0 */
		    /* makes sense only if a != 0 */
		    if ((fabs(ypc[i+1] - ypc[i]) > 0.000000001)) {
			x0 = - b / a;
			if (x0 >= 0)
			    wm = abs(wm - 1);
		    } 	    
		}
		/* Case 2: Infinite slope
		   verify on the right of the point, i.e. 
		   xi >0 */
		if ((fabs(xpc[i+1] - xpc[i]) < 0.000000001))
		{
		    if (xpc[i] >= 0)
			wm = abs(wm - 1);
		}
	    }
	}
	
	/* If even number: outside. Inside otherwise */
	if (wm == 0)
	    deds[j] = 0;
    }
    
    
    /* Free memory */
    freevec(xpc);
    freevec(ypc);
}




/* verification of inout with R */

void inoutr(double *xr, double *yr, double *xpr, double *ypr,
	    int *dedsr, int *nxr, int *npr)
{
    /* Declaration */
    int i, nx, np, *deds;
    double *x, *y, *xp, *yp;
  
    /* Memory allocation */
    nx = *nxr;
    np = *npr;
    vecalloc(&x, nx);
    vecalloc(&y, nx);
    vecalloc(&xp, np);
    vecalloc(&yp, np);
    vecintalloc(&deds, nx);
    
    /* R to C */
    for (i = 1; i <= nx; i++) {
	x[i] = xr[i-1];
	y[i] = yr[i-1];
    }
    
    for (i = 1; i <= np; i++) {
	xp[i] = xpr[i-1];
	yp[i] = ypr[i-1];
    }
    
    /* test of inout */
    inout(x, y, xp, yp, deds);
    
    /* C to R */
    for (i=1; i<=nx; i++) {
	dedsr[i-1] = deds[i];
    }
    
    /* Free memory */
    freevec(x);
    freevec(y);
    freevec(xp);
    freevec(yp);
    freeintvec(deds);
}




/***********************************************************
  Rasterization of a polygon: xp and yp are the coordinates
  of the polygon, xg and yg (not of the same length) are the
  coordinates of the rows and columns of the grid and 
  carte is a raster map.
*************************************************************/


void rastpol(double *xp, double *yp, double *xg, double *yg,
	     double **carte)
{
    /* Declaration of variables */
    int i, j, nl, nc, k, *deds;
    double *nxc, *nyc;
    
    /* Memory allocation */
    nl = xg[0];
    nc = yg[0];
    vecalloc(&nxc, nl*nc);
    vecalloc(&nyc, nl*nc);
    vecintalloc(&deds, nl*nc);
    
    /* Empties the map */
    for (i = 1; i <= nl; i++) {
	for (j = 1; j <= nc; j++) {
	    carte[i][j] = 0;
	}
    }
    
    /* Output of the coordinates of the pixels of the grid */
    k = 1;
    for (i = 1; i <= nl; i++) {
	for (j = 1; j <= nc; j++) {
	    nxc[k] = xg[i];
	    nyc[k] = yg[j];
	    k++;
	}
    }
    
    /* inout on these pixels */
    inout(nxc, nyc, xp, yp, deds);
    
    /* Fills the grid */
    k = 1;
    for (i = 1; i <= nl; i++) {
	for (j = 1; j <= nc; j++) {
	    carte[i][j] = (double) deds[k];
	    k++;
	}
    }
    
    /* Free memory */
    freevec(nxc);
    freevec(nyc);
    freeintvec(deds);
}



/* ****************************************************************
   *                                                              *
   *   Verification of rastpol with R                             *
   *                                                              *
   **************************************************************** */

void rastpolaire(double *xpr, double *ypr, double *xgr, double *ygr,
		 double *carter, int *nlg, int *ncg, int *nvp)
{
    /* Declaration */
    int i, j, k, nl, nc, nv;
    double *xp, *yp, *xg, *yg, **carte;
    
    /* Memory allocation */
    nl = *nlg;
    nc = *ncg;
    nv = *nvp;
    vecalloc(&xp, nv);
    vecalloc(&yp, nv);
    vecalloc(&xg, nl);
    vecalloc(&yg, nc);
    taballoc(&carte, nl, nc);
    
    /* R to C */
    for (i = 1; i <= nv; i++) {
	xp[i] = xpr[i-1];
	yp[i] = ypr[i-1];
    }
    
    for (i = 1; i <= nl; i++) {
	xg[i] = xgr[i-1];
    }
    
    for (i = 1; i <= nc; i++) {
	yg[i] = ygr[i-1];
    }
    
    k=0;
    for (i=1; i<=nl; i++) {
	for(j=1; j<=nc; j++) {
	    carte[i][j] = carter[k];
	    k++;
	}
    }
    
    /* call to rastpol */
    rastpol(xp, yp, xg, yg, carte);
    
    /* C to R */
    k=0;
    for (i=1; i<=nl; i++) {
	for (j=1; j<=nc; j++) {
	    carter[k] = carte[i][j];
      k++;
	}
    }
    
    /* Free memory */
    freevec(xp);
    freevec(yp);
    freevec(xg);
    freevec(yg);
    freetab(carte);
}





/* ********************************************************************
 *                                                                    *
 *            Diminish the resolution of a map                        *
 *                                                                    *
 * ******************************************************************** */



/* For factor maps */

void regroufacasc(double **asce, double **ascs, int *np,
		  int *nlev)
{
    /* declaration of variables */
    int i, j, k, l, m, dr, fr, dc, fc, nrs;
    int ncs, nl, *ll, max, vm, na, *vecmax;
    
    /* Memory allocation */
    nrs = ascs[0][0];
    ncs = ascs[1][0];
    nl = *nlev;
    vecintalloc(&ll, nl);
  
    /* loop to delete */
    for (i = 1; i <= nrs; i++) {
	for (j = 1; j <= ncs; j++) {
	    
	    /* extracts the corresponding subtable */
	    dr = (i-1)*(*np) + 1;
	    fr = i*(*np);
	    dc = (j-1)*(*np) + 1;
	    fc = j*(*np);
	    
	    /* empty ll */
	    for (m = 1; m <= nl; m++) {
		ll[m] = 0;
	    }
	    
	    /* One numbers the levels */
	    na = 0;
	    for (k = dr; k <= fr; k++) {
		for (l = dc; l <= fc; l++) {
		    if (fabs(asce[k][l] + 9999) > 0.000000001)
			ll[(int) asce[k][l]]++;
		    if (fabs(asce[k][l] + 9999) < 0.000000001)
			na++;
		}
	    }
	    
	    if (na != (*np)*(*np)) {
		/* One computes the maximum number */
		vm = ll[1];
		for (k = 2; k <= nl; k++) {
		    if (ll[k] >= vm) {
			vm = ll[k];
		    }
		}
		
		/* ... and the number OF max */
		max = 0;
		for (k = 1; k <= nl; k++) {
		    if (ll[k] == vm) {
			max++;
		    }
		}
		
		/* one identifies the levels for which the number is max */
		vecintalloc(&vecmax, max);
		l = 1;
		for (k = 1; k<=nl; k++) {
		    if (ll[k] == vm) {
			vecmax[1] = k;
		    }
		}
		
		/* Random sample of the levels in case of equality */
		if (max > 1) {
		    getpermutation(vecmax, i*j); /* random row */
		}
		ascs[i][j] = (double) vecmax[1];
		freeintvec(vecmax);
	    } else {
		ascs[i][j] = -9999;
	    }
	    
	}
    }
    /* free memory */
    freeintvec(ll);
}



/* Regroufacasc version for R */

void regroufacascr(double *ascer, double *ascsr, int *npr,
		   int *nlevr, int *nle, int *nce, int *nls, 
		   int *ncs)
{
    /* Declaration of the variables */
    int i,j,k, np, nlev;
    double **asce, **ascs;
    
    /* Memory Allocation */
    np = *npr;
    nlev = *nlevr;
    taballoc(&asce, *nle, *nce);
    taballoc(&ascs, *nls, *ncs);
    
    /* R to C */
    k =0;
    for (i = 1; i <= *nle; i++) {
	for (j = 1; j <= *nce; j++) {
	    asce[i][j] = ascer[k];
	    k++;
	}
    }
    
    /* function */
    regroufacasc(asce, ascs, &np, &nlev);
    
    k =0;
    for (i = 1; i <= *nls; i++) {
	for (j = 1; j <= *ncs; j++) {
	    ascsr[k] = ascs[i][j];
	    k++;
	}
    }
    
    /* Free memory */
    freetab(asce);
    freetab(ascs);
}




/* regrouascnum for numeric maps */

void regrouascnum(double **ascent, double **ascso)
{
    /* Declaration */
    int i, j, k, l, n, nle, nls, ncs, nreg;
    double moy, tmp;
    
    /* Definition of the variables */
    nle = ascent[0][0];
    nls = ascso[0][0];
    ncs = ascso[1][0];
    nreg = nle/nls;
    
    /* Computes the mean */
    for (i = 1; i <= nls; i++) {
	for (j = 1; j <= ncs; j++) {
	    moy = 0;
	    n = 0;
	    for (k = 1; k <= nreg; k++) {
		for (l = 1; l <= nreg; l++) {
		    tmp = ascent[((i - 1) * nreg) + k][((j - 1) * nreg) + l];
		    if (fabs(tmp + 9999) > 0.000000001) {
			moy = tmp + moy;
		    }
		    if (fabs(tmp + 9999) < 0.000000001) {
			n++;
		    }
		}
	    }
	    if (n == (nreg * nreg)) {
		ascso[i][j] = -9999;
	    } else {
		ascso[i][j] = moy / (((double) (nreg * nreg))- ((double) n));
		
	    }
	}
    }
}


/* Version for R */

void regrouascnumr(double *ascentr, double *ascsor, double *nler, double *ncer,
		   double *nlsr, double *ncsr)
{
    /* Declaration of variables */
    int i, j, k, nle, nce, nls, ncs;
    double **ascent, **ascso;
    
    /* Memory Allocation */
    nle = *nler;
    nce = *ncer;
    nls = *nlsr;
    ncs = *ncsr;
    
    taballoc(&ascent, nle, nce);
    taballoc(&ascso, nls, ncs);
  
    
    /* R to C */
    k = 0;
    for (i = 1; i <= nle; i++) {
	for (j = 1; j <= nce; j++) {
	    ascent[i][j] = ascentr[k];
	    k++;
	}
    }
    
    k = 0;
    for (i = 1; i <= nls; i++) {
	for (j = 1; j <= ncs; j++) {
	    ascso[i][j] = ascsor[k];
	    k++;
	}
    }
    
    /* procedure C */
    regrouascnum(ascent, ascso);
    
    /* C to R */
    k = 0;
    for (i = 1; i <= nls; i++) {
	for (j = 1; j <= ncs; j++) {
	    ascsor[k] = ascso[i][j];
	    k++;
	}
    }
    
    /* Free memory */
    freetab(ascso);
    freetab(ascent);
}










/* *********************************************************************
 *                                                                     *
 *                   Buffer on a line                                  *
 *                                                                     *
 ***********************************************************************/

/* given a line, ligpoly returns a buffer polygon containing the line */

void ligpoly(double *x, double *y, double r, double *xp, double *yp)
{
    /* Declaration */
    double x1, x2, y1, y2, xx, yy, alpha, beta, xim, xsm, yim, ysm, gamma;
    double xip, xsp, yip, ysp;
    
    x1 = x[1];
    x2 = x[2];
    y1 = y[1];
    y2 = y[2];
    xx = x2 - x1;
    yy = y2 - y1;
    
    alpha = atan(yy/xx);
    beta = alpha - (3.1415926/2);
    xim = x1 + r * (cos(beta));
    xsm = x2 + r * (cos(beta));
    yim = y1 + r * (sin(beta));
    ysm = y2 + r * (sin(beta));
    
    gamma = alpha + (3.1415926/2);
    xip = x1 + r * (cos(gamma));
    xsp = x2 + r * (cos(gamma));
    yip = y1 + r * (sin(gamma));
    ysp = y2 + r * (sin(gamma));
    
    xp[1] = xim;
    xp[2] = xsm;
    xp[3] = xsp;
    xp[4] = xip;
    xp[5] = xim;
    
    yp[1] = yim;
    yp[2] = ysm;
    yp[3] = ysp;
    yp[4] = yip;
    yp[5] = yim;
    
}


/* main function */

void buflig(double **x, double r, double **carte, double *xg, double *yg)
{
    /* Declaration */
    int i, j, k, nloc, nr, nc;
    double **x1, **x2, *xl, *yl, *xp, *yp, **cartebis;
    
    /* Memory allocation */
    nloc = x[0][0];
    k = 0;
    nr = carte[0][0];
    nc = carte[1][0];
    
    vecalloc(&xl, 2);
    vecalloc(&yl, 2);
    vecalloc(&xp, 5);
    vecalloc(&yp, 5);
    taballoc(&x1, nloc-1, 2);
    taballoc(&x2, nloc-1, 2);
    taballoc(&cartebis, nr, nc);
    
    /* Creates the two tables */
    for (i = 1; i <= nloc; i++) {
	if (i > 1) {
	    x2[i-1][1] = x[i][1];
	    x2[i-1][2] = x[i][2];
	}
	if (i < nloc) {
	    x1[i][1] = x[i][1];
	    x1[i][2] = x[i][2];
	}
    }
    
    /* Sets the map to 0 */
    for (i = 1; i <= nr; i++) {
	for (j = 1; j <= nc; j++) {
	    carte[i][j] = 0;
	}
    }
    
    
    /* Buffer around the line */
    for (i = 1; i <= (nloc-1); i++) {
	xl[1] = x1[i][1];
	xl[2] = x2[i][1];
	yl[1] = x1[i][2];
	yl[2] = x2[i][2];
	
	ligpoly(xl, yl, r, xp, yp);
	
	rastpol(xp, yp, xg, yg, cartebis);
	
	for (j = 1; j <= nr; j++) {
	    for (k = 1; k <= nc; k++) {
		carte[j][k] = cartebis[j][k] + carte[j][k];
	    }
	}
    }
    
    /* Free memory */
    freevec(xl);
    freevec(yl);
    freevec(xp);
    freevec(yp);
    freetab(x1);
    freetab(x2);
    freetab(cartebis);
    
}


/* For external call from xithin R */
void bufligr(double *xr, double *rr, double *carter, 
	     double *xgr, double *ygr, int *nlr, int *ncr, 
	     int *nlocr)
{
    /* Declaration */
    int i, j, k, nc, nl, nloc;
    double **x, r, **carte, *xg, *yg;
    
    /* Memory allocation */
    nc = *ncr;
    nl = *nlr;
    nloc = *nlocr;
    r = *rr;
    
    taballoc(&x, nloc, 2);
    taballoc(&carte, nl, nc);
    vecalloc(&xg, nl);
    vecalloc(&yg, nc);
    
    
    /* R to C */
    k = 0;
    for (i=1; i<= nl; i++) {
	for (j = 1; j<=nc; j++) {
	    carte[i][j]=carter[k];
	    k++;
	}
    }
    
    k = 0;
    for (i=1; i<= nloc; i++) {
	for (j = 1; j<=2; j++) {
	    x[i][j]=xr[k];
	    k++;
	}
    }
    
    for (i = 1; i <= nl; i++) {
	xg[i] = xgr[i-1];
    }
    
    for (i = 1; i <= nc; i++) {
	yg[i] = ygr[i-1];
    }
    
    /* Main function */
    buflig(x, r, carte, xg, yg);
    
    /* C to R */
    k = 0;
    for (i=1; i<= nl; i++) {
	for (j = 1; j<=nc; j++) {
	    carter[k]=carte[i][j];
	    k++;
	}
    }
    
    /* Free memory */
    freetab(x);
    freetab(carte);
    freevec(xg);
    freevec(yg);
    
}


/* Computes Euclidean distances from a map of class asc */

void distxy(double **xy1, double **xy2, double *di)
{
    /* Declaration */
    int i, j, n1, n2;
    double *dib, mi;
    
    /* Memory allocation */
    n1 = xy1[0][0];
    n2 = xy2[0][0];
    
    vecalloc(&dib, n2);
    
    /* Euclidean distances */
    for (i = 1; i <= n1; i++) {
	for (j = 1; j <= n2; j++) {
	    dib[j] = sqrt( ((xy1[i][1] - xy2[j][1]) * (xy1[i][1] - xy2[j][1])) + 
			   ((xy1[i][2] - xy2[j][2]) * (xy1[i][2] - xy2[j][2])));
	}
	mi = dib[1];
	for (j = 2; j <= n2; j++) {
	    if (mi > dib[j]) {
		mi = dib[j];
	    }
	}
	di[i] = mi;
    }
    
    /* Free memory */
    freevec(dib);
}


/* For external call from R */
void distxyr(double *xy1r, double *xy2r, int *n1r, 
	     int *n2r, double *dire)
{
    /* Declaration */
    int i, j, k, n1, n2;
    double **xy1, **xy2, *di;
    
    /* Memory allocation */
    n1 = *n1r;
    n2 = *n2r;
    
    taballoc(&xy1, n1, 2);
    taballoc(&xy2, n2, 2);
    vecalloc(&di, n1);
    
    /* R to C */
    k = 0;
    for (i = 1; i <= n1; i++) {
	for (j = 1; j <= 2; j++) {
	    xy1[i][j] = xy1r[k];
	    k++;
	}
    }
    
    k = 0;
    for (i = 1; i <= n2; i++) {
	for (j = 1; j <= 2; j++) {
	    xy2[i][j] = xy2r[k];
	    k++;
	}
    }
    
    /* The function */
    distxy(xy1, xy2, di);
    
    /* C to R */
    for (i = 1; i <= n1; i++) {
	dire[i-1] = di[i];
    }
    
    /* Free memory */
    freetab(xy1);
    freetab(xy2);
    freevec(di);
}





