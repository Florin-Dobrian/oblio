// Amd.cc — Approximate Minimum Degree ordering.
//
// Merged from SuiteSparse AMD 3.3.4 (Davis, Amestoy, Duff).
// Original license: BSD-3-clause.
// C sources merged into single C++ file, debug code stripped,
// malloc/free replaced with std::vector in entry point.

#include <cstdint>
#include <cstddef>
#include <cstdlib>
#include <cmath>
#include <cstdio>
#include <vector>
#include <algorithm>

namespace {  // anonymous namespace — all internal

using Int = int32_t;
using UInt = uint32_t;
constexpr Int Int_MAX_VAL = INT32_MAX;

#define EMPTY (-1)
#define FLIP(i) (-(i)-2)
#define UNFLIP(i) ((i < EMPTY) ? FLIP(i) : (i))
#define MAX(a,b) (((a) > (b)) ? (a) : (b))
#define MIN(a,b) (((a) < (b)) ? (a) : (b))
#define IMPLIES(p,q) (!(p) || (q))
#define TRUE (1)
#define FALSE (0)

#define AMD_CONTROL 5
#define AMD_INFO 20
#define AMD_DENSE 0
#define AMD_AGGRESSIVE 1
#define AMD_DEFAULT_DENSE 10.0
#define AMD_DEFAULT_AGGRESSIVE 1
#define AMD_STATUS 0
#define AMD_N 1
#define AMD_NZ 2
#define AMD_SYMMETRY 3
#define AMD_NZDIAG 4
#define AMD_NZ_A_PLUS_AT 5
#define AMD_NDENSE 6
#define AMD_MEMORY 7
#define AMD_NCMPA 8
#define AMD_LNZ 9
#define AMD_NDIV 10
#define AMD_NMULTSUBS_LDL 11
#define AMD_NMULTSUBS_LU 12
#define AMD_DMAX 13
#define AMD_OK 0
#define AMD_OUT_OF_MEMORY -1
#define AMD_INVALID -2
#define AMD_OK_BUT_JUMBLED 1

#ifndef SIZE_T_MAX
#define SIZE_T_MAX ((size_t)(-1))
#endif

// === amd_valid ===
static int AMD_valid
(
    /* inputs, not modified on output: */
    Int n_row,		/* A is n_row-by-n_col */
    Int n_col,
    const Int Ap [ ],	/* column pointers of A, of size n_col+1 */
    const Int Ai [ ]	/* row indices of A, of size nz = Ap [n_col] */
)
{
    Int nz, j, p1, p2, ilast, i, p ;
    int result = AMD_OK ;

    if (n_row < 0 || n_col < 0 || Ap == NULL || Ai == NULL)
    {
	return (AMD_INVALID) ;
    }
    nz = Ap [n_col] ;
    if (Ap [0] != 0 || nz < 0)
    {
	/* column pointers must start at Ap [0] = 0, and Ap [n] must be >= 0 */
	return (AMD_INVALID) ;
    }
    for (j = 0 ; j < n_col ; j++)
    {
	p1 = Ap [j] ;
	p2 = Ap [j+1] ;
	if (p1 > p2)
	{
	    /* column pointers must be ascending */
	    return (AMD_INVALID) ;
	}
	ilast = EMPTY ;
	for (p = p1 ; p < p2 ; p++)
	{
	    i = Ai [p] ;
	    if (i < 0 || i >= n_row)
	    {
		/* row index out of range */
		return (AMD_INVALID) ;
	    }
	    if (i <= ilast)
	    {
		/* row index unsorted, or duplicate entry present */
		result = AMD_OK_BUT_JUMBLED ;
	    }
	    ilast = i ;
	}
    }
    return (result) ;
}

// === amd_aat ===
static size_t AMD_aat	/* returns nz in A+A' */
(
    Int n,
    const Int Ap [ ],
    const Int Ai [ ],
    Int Len [ ],	/* Len [j]: length of column j of A+A', excl diagonal*/
    Int Tp [ ],		/* workspace of size n */
    double Info [ ]
)
{
    Int p1, p2, p, i, j, pj, pj2, k, nzdiag, nzboth, nz ;
    double sym ;
    size_t nzaat ;

#ifndef NDEBUG
    for (k = 0 ; k < n ; k++) Tp [k] = EMPTY ;
#endif

    if (Info != (double *) NULL)
    {
	/* clear the Info array, if it exists */
	for (i = 0 ; i < AMD_INFO ; i++)
	{
	    Info [i] = EMPTY ;
	}
	Info [AMD_STATUS] = AMD_OK ;
    }

    for (k = 0 ; k < n ; k++)
    {
	Len [k] = 0 ;
    }

    nzdiag = 0 ;
    nzboth = 0 ;
    nz = Ap [n] ;

    for (k = 0 ; k < n ; k++)
    {
	p1 = Ap [k] ;
	p2 = Ap [k+1] ;

	/* construct A+A' */
	for (p = p1 ; p < p2 ; )
	{
	    /* scan the upper triangular part of A */
	    j = Ai [p] ;
	    if (j < k)
	    {
		/* entry A (j,k) is in the strictly upper triangular part,
		 * add both A (j,k) and A (k,j) to the matrix A+A' */
		Len [j]++ ;
		Len [k]++ ;
		p++ ;
	    }
	    else if (j == k)
	    {
		/* skip the diagonal */
		p++ ;
		nzdiag++ ;
		break ;
	    }
	    else /* j > k */
	    {
		/* first entry below the diagonal */
		break ;
	    }
	    /* scan lower triangular part of A, in column j until reaching
	     * row k.  Start where last scan left off. */
	    pj2 = Ap [j+1] ;
	    for (pj = Tp [j] ; pj < pj2 ; )
	    {
		i = Ai [pj] ;
		if (i < k)
		{
		    /* A (i,j) is only in the lower part, not in upper.
		     * add both A (i,j) and A (j,i) to the matrix A+A' */
		    Len [i]++ ;
		    Len [j]++ ;
		    pj++ ;
		}
		else if (i == k)
		{
		    /* entry A (k,j) in lower part and A (j,k) in upper */
		    pj++ ;
		    nzboth++ ;
		    break ;
		}
		else /* i > k */
		{
		    /* consider this entry later, when k advances to i */
		    break ;
		}
	    }
	    Tp [j] = pj ;
	}
	/* Tp [k] points to the entry just below the diagonal in column k */
	Tp [k] = p ;
    }

    /* clean up, for remaining mismatched entries */
    for (j = 0 ; j < n ; j++)
    {
	for (pj = Tp [j] ; pj < Ap [j+1] ; pj++)
	{
	    i = Ai [pj] ;
	    /* A (i,j) is only in the lower part, not in upper.
	     * add both A (i,j) and A (j,i) to the matrix A+A' */
	    Len [i]++ ;
	    Len [j]++ ;
	}
    }

    /* --------------------------------------------------------------------- */
    /* compute the symmetry of the nonzero pattern of A */
    /* --------------------------------------------------------------------- */

    /* Given a matrix A, the symmetry of A is:
     *	B = tril (spones (A), -1) + triu (spones (A), 1) ;
     *  sym = nnz (B & B') / nnz (B) ;
     *  or 1 if nnz (B) is zero.
     */

    if (nz == nzdiag)
    {
	sym = 1 ;
    }
    else
    {
	sym = (2 * (double) nzboth) / ((double) (nz - nzdiag)) ;
    }

    nzaat = 0 ;
    for (k = 0 ; k < n ; k++)
    {
	nzaat += Len [k] ;
    }


    if (Info != (double *) NULL)
    {
	Info [AMD_STATUS] = AMD_OK ;
	Info [AMD_N] = n ;
	Info [AMD_NZ] = nz ;
	Info [AMD_SYMMETRY] = sym ;	    /* symmetry of pattern of A */
	Info [AMD_NZDIAG] = nzdiag ;	    /* nonzeros on diagonal of A */
	Info [AMD_NZ_A_PLUS_AT] = nzaat ;   /* nonzeros in A+A' */
    }

    return (nzaat) ;
}

// === amd_preprocess ===
static void AMD_preprocess
(
    Int n,		/* input matrix: A is n-by-n */
    const Int Ap [ ],	/* size n+1 */
    const Int Ai [ ],	/* size nz = Ap [n] */

    /* output matrix R: */
    Int Rp [ ],		/* size n+1 */
    Int Ri [ ],		/* size nz (or less, if duplicates present) */

    Int W [ ],		/* workspace of size n */
    Int Flag [ ]	/* workspace of size n */
)
{

    /* --------------------------------------------------------------------- */
    /* local variables */
    /* --------------------------------------------------------------------- */

    Int i, j, p, p2 ;


    /* --------------------------------------------------------------------- */
    /* count the entries in each row of A (excluding duplicates) */
    /* --------------------------------------------------------------------- */

    for (i = 0 ; i < n ; i++)
    {
	W [i] = 0 ;		/* # of nonzeros in row i (excl duplicates) */
	Flag [i] = EMPTY ;	/* Flag [i] = j if i appears in column j */
    }
    for (j = 0 ; j < n ; j++)
    {
	p2 = Ap [j+1] ;
	for (p = Ap [j] ; p < p2 ; p++)
	{
	    i = Ai [p] ;
	    if (Flag [i] != j)
	    {
		/* row index i has not yet appeared in column j */
		W [i]++ ;	    /* one more entry in row i */
		Flag [i] = j ;	    /* flag row index i as appearing in col j*/
	    }
	}
    }

    /* --------------------------------------------------------------------- */
    /* compute the row pointers for R */
    /* --------------------------------------------------------------------- */

    Rp [0] = 0 ;
    for (i = 0 ; i < n ; i++)
    {
	Rp [i+1] = Rp [i] + W [i] ;
    }
    for (i = 0 ; i < n ; i++)
    {
	W [i] = Rp [i] ;
	Flag [i] = EMPTY ;
    }

    /* --------------------------------------------------------------------- */
    /* construct the row form matrix R */
    /* --------------------------------------------------------------------- */

    /* R = row form of pattern of A */
    for (j = 0 ; j < n ; j++)
    {
	p2 = Ap [j+1] ;
	for (p = Ap [j] ; p < p2 ; p++)
	{
	    i = Ai [p] ;
	    if (Flag [i] != j)
	    {
		/* row index i has not yet appeared in column j */
		Ri [W [i]++] = j ;  /* put col j in row i */
		Flag [i] = j ;	    /* flag row index i as appearing in col j*/
	    }
	}
    }

#ifndef NDEBUG
    for (j = 0 ; j < n ; j++)
    {
    }
#endif
}

// === amd_post_tree ===
static Int AMD_post_tree
(
    Int root,			/* root of the tree */
    Int k,			/* start numbering at k */
    Int Child [ ],		/* input argument of size nn, undefined on
				 * output.  Child [i] is the head of a link
				 * list of all nodes that are children of node
				 * i in the tree. */
    const Int Sibling [ ],	/* input argument of size nn, not modified.
				 * If f is a node in the link list of the
				 * children of node i, then Sibling [f] is the
				 * next child of node i.
				 */
    Int Order [ ],		/* output order, of size nn.  Order [i] = k
				 * if node i is the kth node of the reordered
				 * tree. */
    Int Stack [ ]		/* workspace of size nn */
#ifndef NDEBUG
    , Int nn			/* nodes are in the range 0..nn-1. */
#endif
)
{
    Int f, head, h, i ;

#if 0
    /* --------------------------------------------------------------------- */
    /* recursive version (Stack [ ] is not used): */
    /* --------------------------------------------------------------------- */

    /* this is simple, but can cause stack overflow if nn is large */
    i = root ;
    for (f = Child [i] ; f != EMPTY ; f = Sibling [f])
    {
	k = AMD_post_tree (f, k, Child, Sibling, Order, Stack, nn) ;
    }
    Order [i] = k++ ;
    return (k) ;
#endif

    /* --------------------------------------------------------------------- */
    /* non-recursive version, using an explicit stack */
    /* --------------------------------------------------------------------- */

    /* push root on the stack */
    head = 0 ;
    Stack [0] = root ;

    while (head >= 0)
    {
	/* get head of stack */
	i = Stack [head] ;

	if (Child [i] != EMPTY)
	{
	    /* the children of i are not yet ordered */
	    /* push each child onto the stack in reverse order */
	    /* so that small ones at the head of the list get popped first */
	    /* and the biggest one at the end of the list gets popped last */
	    for (f = Child [i] ; f != EMPTY ; f = Sibling [f])
	    {
		head++ ;
	    }
	    h = head ;
	    for (f = Child [i] ; f != EMPTY ; f = Sibling [f])
	    {
		Stack [h--] = f ;
	    }

	    /* delete child list so that i gets ordered next time we see it */
	    Child [i] = EMPTY ;
	}
	else
	{
	    /* the children of i (if there were any) are already ordered */
	    /* remove i from the stack and order it.  Front i is kth front */
	    head-- ;
	    Order [i] = k++ ;
	}

#ifndef NDEBUG
	for (h = head ; h >= 0 ; h--)
	{
	    Int j = Stack [h] ;
	}
#endif

    }
    return (k) ;
}

// === amd_postorder ===
static void AMD_postorder
(
    /* inputs, not modified on output: */
    Int nn,		/* nodes are in the range 0..nn-1 */
    Int Parent [ ],	/* Parent [j] is the parent of j, or EMPTY if root */
    Int Nv [ ],		/* Nv [j] > 0 number of pivots represented by node j,
			 * or zero if j is not a node. */
    Int Fsize [ ],	/* Fsize [j]: size of node j */

    /* output, not defined on input: */
    Int Order [ ],	/* output post-order */

    /* workspaces of size nn: */
    Int Child [ ],
    Int Sibling [ ],
    Int Stack [ ]
)
{
    Int i, j, k, parent, frsize, f, fprev, maxfrsize, bigfprev, bigf, fnext ;

    for (j = 0 ; j < nn ; j++)
    {
	Child [j] = EMPTY ;
	Sibling [j] = EMPTY ;
    }

    /* --------------------------------------------------------------------- */
    /* place the children in link lists - bigger elements tend to be last */
    /* --------------------------------------------------------------------- */

    for (j = nn-1 ; j >= 0 ; j--)
    {
	if (Nv [j] > 0)
	{
	    /* this is an element */
	    parent = Parent [j] ;
	    if (parent != EMPTY)
	    {
		/* place the element in link list of the children its parent */
		/* bigger elements will tend to be at the end of the list */
		Sibling [j] = Child [parent] ;
		Child [parent] = j ;
	    }
	}
    }

#ifndef NDEBUG
    {
	Int nels, ff, nchild ;
	nels = 0 ;
	for (j = 0 ; j < nn ; j++)
	{
	    if (Nv [j] > 0)
	    {
		/* this is an element */
		/* dump the link list of children */
		nchild = 0 ;
		for (ff = Child [j] ; ff != EMPTY ; ff = Sibling [ff])
		{
		    nchild++ ;
		}
		parent = Parent [j] ;
		if (parent != EMPTY)
		{
		}
		nels++ ;
	    }
	}
    }
#endif

    /* --------------------------------------------------------------------- */
    /* place the largest child last in the list of children for each node */
    /* --------------------------------------------------------------------- */

    for (i = 0 ; i < nn ; i++)
    {
	if (Nv [i] > 0 && Child [i] != EMPTY)
	{

#ifndef NDEBUG
	    Int nchild ;
	    nchild = 0 ;
	    for (f = Child [i] ; f != EMPTY ; f = Sibling [f])
	    {
		nchild++ ;
	    }
#endif

	    /* find the biggest element in the child list */
	    fprev = EMPTY ;
	    maxfrsize = EMPTY ;
	    bigfprev = EMPTY ;
	    bigf = EMPTY ;
	    for (f = Child [i] ; f != EMPTY ; f = Sibling [f])
	    {
		frsize = Fsize [f] ;
		if (frsize >= maxfrsize)
		{
		    /* this is the biggest seen so far */
		    maxfrsize = frsize ;
		    bigfprev = fprev ;
		    bigf = f ;
		}
		fprev = f ;
	    }

	    fnext = Sibling [bigf] ;


	    if (fnext != EMPTY)
	    {
		/* if fnext is EMPTY then bigf is already at the end of list */

		if (bigfprev == EMPTY)
		{
		    /* delete bigf from the element of the list */
		    Child [i] = fnext ;
		}
		else
		{
		    /* delete bigf from the middle of the list */
		    Sibling [bigfprev] = fnext ;
		}

		/* put bigf at the end of the list */
		Sibling [bigf] = EMPTY ;
		Sibling [fprev] = bigf ;
	    }

#ifndef NDEBUG
	    for (f = Child [i] ; f != EMPTY ; f = Sibling [f])
	    {
		nchild-- ;
	    }
#endif

	}
    }

    /* --------------------------------------------------------------------- */
    /* postorder the assembly tree */
    /* --------------------------------------------------------------------- */

    for (i = 0 ; i < nn ; i++)
    {
	Order [i] = EMPTY ;
    }

    k = 0 ;

    for (i = 0 ; i < nn ; i++)
    {
	if (Parent [i] == EMPTY && Nv [i] > 0)
	{
	    k = AMD_post_tree (i, k, Child, Sibling, Order, Stack
#ifndef NDEBUG
		, nn
#endif
		) ;
	}
    }
}

// Forward declarations
static void AMD_2(Int, Int[], Int[], Int[], Int, Int, Int[], Int[], Int[], Int[], Int[], Int[], Int[], double[], double[]);

// === clear_flag (from amd_2) ===
static Int clear_flag(Int wflg, Int wbig, Int W[], Int n)
{
    Int x;
    if (wflg < 2 || wflg >= wbig)
    {
        for (x = 0; x < n; x++)
        {
            if (W[x] != 0) W[x] = 1;
        }
        wflg = 2;
    }
    return (wflg);
}

// === amd_1 ===
static void AMD_1
(
    Int n,		/* n > 0 */
    const Int Ap [ ],	/* input of size n+1, not modified */
    const Int Ai [ ],	/* input of size nz = Ap [n], not modified */
    Int P [ ],		/* size n output permutation */
    Int Pinv [ ],	/* size n output inverse permutation */
    Int Len [ ],	/* size n input, undefined on output */
    Int slen,		/* slen >= sum (Len [0..n-1]) + 7n,
			 * ideally slen = 1.2 * sum (Len) + 8n */
    Int S [ ],		/* size slen workspace */
    double Control [ ],	/* input array of size AMD_CONTROL */
    double Info [ ]	/* output array of size AMD_INFO */
)
{
    Int i, j, k, p, pfree, iwlen, pj, p1, p2, pj2, *Iw, *Pe, *Nv, *Head,
	*Elen, *Degree, *s, *W, *Sp, *Tp ;

    /* --------------------------------------------------------------------- */
    /* construct the matrix for AMD_2 */
    /* --------------------------------------------------------------------- */


    iwlen = slen - 6*n ;
    s = S ;
    Pe = s ;	    s += n ;
    Nv = s ;	    s += n ;
    Head = s ;	    s += n ;
    Elen = s ;	    s += n ;
    Degree = s ;    s += n ;
    W = s ;	    s += n ;
    Iw = s ;	    s += iwlen ;


    /* construct the pointers for A+A' */
    Sp = Nv ;			/* use Nv and W as workspace for Sp and Tp [ */
    Tp = W ;
    pfree = 0 ;
    for (j = 0 ; j < n ; j++)
    {
	Pe [j] = pfree ;
	Sp [j] = pfree ;
	pfree += Len [j] ;
    }

    /* Note that this restriction on iwlen is slightly more restrictive than
     * what is strictly required in AMD_2.  AMD_2 can operate with no elbow
     * room at all, but it will be very slow.  For better performance, at
     * least size-n elbow room is enforced. */

#ifndef NDEBUG
    for (p = 0 ; p < iwlen ; p++) Iw [p] = EMPTY ;
#endif

    for (k = 0 ; k < n ; k++)
    {
	p1 = Ap [k] ;
	p2 = Ap [k+1] ;

	/* construct A+A' */
	for (p = p1 ; p < p2 ; )
	{
	    /* scan the upper triangular part of A */
	    j = Ai [p] ;
	    if (j < k)
	    {
		/* entry A (j,k) in the strictly upper triangular part */
		Iw [Sp [j]++] = k ;
		Iw [Sp [k]++] = j ;
		p++ ;
	    }
	    else if (j == k)
	    {
		/* skip the diagonal */
		p++ ;
		break ;
	    }
	    else /* j > k */
	    {
		/* first entry below the diagonal */
		break ;
	    }
	    /* scan lower triangular part of A, in column j until reaching
	     * row k.  Start where last scan left off. */
	    pj2 = Ap [j+1] ;
	    for (pj = Tp [j] ; pj < pj2 ; )
	    {
		i = Ai [pj] ;
		if (i < k)
		{
		    /* A (i,j) is only in the lower part, not in upper */
		    Iw [Sp [i]++] = j ;
		    Iw [Sp [j]++] = i ;
		    pj++ ;
		}
		else if (i == k)
		{
		    /* entry A (k,j) in lower part and A (j,k) in upper */
		    pj++ ;
		    break ;
		}
		else /* i > k */
		{
		    /* consider this entry later, when k advances to i */
		    break ;
		}
	    }
	    Tp [j] = pj ;
	}
	Tp [k] = p ;
    }

    /* clean up, for remaining mismatched entries */
    for (j = 0 ; j < n ; j++)
    {
	for (pj = Tp [j] ; pj < Ap [j+1] ; pj++)
	{
	    i = Ai [pj] ;
	    /* A (i,j) is only in the lower part, not in upper */
	    Iw [Sp [i]++] = j ;
	    Iw [Sp [j]++] = i ;
	}
    }

#ifndef NDEBUG
#endif

    /* Tp and Sp no longer needed ] */

    /* --------------------------------------------------------------------- */
    /* order the matrix */
    /* --------------------------------------------------------------------- */

    AMD_2 (n, Pe, Iw, Len, iwlen, pfree,
	Nv, Pinv, P, Head, Elen, Degree, W, Control, Info) ;
}

// === amd_2 ===
static void AMD_2
(
    Int n,		/* A is n-by-n, where n > 0 */
    Int Pe [ ],		/* Pe [0..n-1]: index in Iw of row i on input */
    Int Iw [ ],		/* workspace of size iwlen. Iw [0..pfree-1]
			 * holds the matrix on input */
    Int Len [ ],	/* Len [0..n-1]: length for row/column i on input */
    Int iwlen,		/* length of Iw. iwlen >= pfree + n */
    Int pfree,		/* Iw [pfree ... iwlen-1] is empty on input */

    /* 7 size-n workspaces, not defined on input: */
    Int Nv [ ],		/* the size of each supernode on output */
    Int Next [ ],	/* the output inverse permutation */
    Int Last [ ],	/* the output permutation */
    Int Head [ ],
    Int Elen [ ],	/* the size columns of L for each supernode */
    Int Degree [ ],
    Int W [ ],

    /* control parameters and output statistics */
    double Control [ ],	/* array of size AMD_CONTROL */
    double Info [ ]	/* array of size AMD_INFO */
)
{

/*
 * Given a representation of the nonzero pattern of a symmetric matrix, A,
 * (excluding the diagonal) perform an approximate minimum (UMFPACK/MA38-style)
 * degree ordering to compute a pivot order such that the introduction of
 * nonzeros (fill-in) in the Cholesky factors A = LL' is kept low.  At each
 * step, the pivot selected is the one with the minimum UMFAPACK/MA38-style
 * upper-bound on the external degree.  This routine can optionally perform
 * aggresive absorption (as done by MC47B in the Harwell Subroutine
 * Library).
 *
 * The approximate degree algorithm implemented here is the symmetric analog of
 * the degree update algorithm in MA38 and UMFPACK (the Unsymmetric-pattern
 * MultiFrontal PACKage, both by Davis and Duff).  The routine is based on the
 * MA27 minimum degree ordering algorithm by Iain Duff and John Reid.
 *
 * This routine is a translation of the original AMDBAR and MC47B routines,
 * in Fortran, with the following modifications:
 *
 * (1) dense rows/columns are removed prior to ordering the matrix, and placed
 *	last in the output order.  The presence of a dense row/column can
 *	increase the ordering time by up to O(n^2), unless they are removed
 *	prior to ordering.
 *
 * (2) the minimum degree ordering is followed by a postordering (depth-first
 *	search) of the assembly tree.  Note that mass elimination (discussed
 *	below) combined with the approximate degree update can lead to the mass
 *	elimination of nodes with lower exact degree than the current pivot
 *	element.  No additional fill-in is caused in the representation of the
 *	Schur complement.  The mass-eliminated nodes merge with the current
 *	pivot element.  They are ordered prior to the current pivot element.
 *	Because they can have lower exact degree than the current element, the
 *	merger of two or more of these nodes in the current pivot element can
 *	lead to a single element that is not a "fundamental supernode".  The
 *	diagonal block can have zeros in it.  Thus, the assembly tree used here
 *	is not guaranteed to be the precise supernodal elemination tree (with
 *	"funadmental" supernodes), and the postordering performed by this
 *	routine is not guaranteed to be a precise postordering of the
 *	elimination tree.
 *
 * (3) input parameters are added, to control aggressive absorption and the
 *	detection of "dense" rows/columns of A.
 *
 * (4) additional statistical information is returned, such as the number of
 *	nonzeros in L, and the flop counts for subsequent LDL' and LU
 *	factorizations.  These are slight upper bounds, because of the mass
 *	elimination issue discussed above.
 *
 * (5) additional routines are added to interface this routine to MATLAB
 *	to provide a simple C-callable user-interface, to check inputs for
 *	errors, compute the symmetry of the pattern of A and the number of
 *	nonzeros in each row/column of A+A', to compute the pattern of A+A',
 *	to perform the assembly tree postordering, and to provide debugging
 *	ouput.  Many of these functions are also provided by the Fortran
 *	Harwell Subroutine Library routine MC47A.
 *
 * (6) both int32_t and int64_t versions are provided.  In the
 *      descriptions below an integer is int32_t or int64_t depending
 *      on which version is being used.

 **********************************************************************
 ***** CAUTION:  ARGUMENTS ARE NOT CHECKED FOR ERRORS ON INPUT.  ******
 **********************************************************************
 ** If you want error checking, a more versatile input format, and a **
 ** simpler user interface, use amd_order or amd_l_order instead.    **
 ** This routine is not meant to be user-callable.                   **
 **********************************************************************

 * ----------------------------------------------------------------------------
 * References:
 * ----------------------------------------------------------------------------
 *
 *  [1] Timothy A. Davis and Iain Duff, "An unsymmetric-pattern multifrontal
 *	method for sparse LU factorization", SIAM J. Matrix Analysis and
 *	Applications, vol. 18, no. 1, pp. 140-158.  Discusses UMFPACK / MA38,
 *	which first introduced the approximate minimum degree used by this
 *	routine.
 *
 *  [2] Patrick Amestoy, Timothy A. Davis, and Iain S. Duff, "An approximate
 *	minimum degree ordering algorithm," SIAM J. Matrix Analysis and
 *	Applications, vol. 17, no. 4, pp. 886-905, 1996.  Discusses AMDBAR and
 *	MC47B, which are the Fortran versions of this routine.
 *
 *  [3] Alan George and Joseph Liu, "The evolution of the minimum degree
 *	ordering algorithm," SIAM Review, vol. 31, no. 1, pp. 1-19, 1989.
 *	We list below the features mentioned in that paper that this code
 *	includes:
 *
 *	mass elimination:
 *	    Yes.  MA27 relied on supervariable detection for mass elimination.
 *
 *	indistinguishable nodes:
 *	    Yes (we call these "supervariables").  This was also in the MA27
 *	    code - although we modified the method of detecting them (the
 *	    previous hval was the true degree, which we no longer keep track
 *	    of).  A supervariable is a set of rows with identical nonzero
 *	    pattern.  All variables in a supervariable are eliminated together.
 *	    Each supervariable has as its numerical name that of one of its
 *	    variables (its principal variable).
 *
 *	quotient graph representation:
 *	    Yes.  We use the term "element" for the cliques formed during
 *	    elimination.  This was also in the MA27 code.  The algorithm can
 *	    operate in place, but it will work more efficiently if given some
 *	    "elbow room."
 *
 *	element absorption:
 *	    Yes.  This was also in the MA27 code.
 *
 *	external degree:
 *	    Yes.  The MA27 code was based on the true degree.
 *
 *	incomplete degree update and multiple elimination:
 *	    No.  This was not in MA27, either.  Our method of degree update
 *	    within MC47B is element-based, not variable-based.  It is thus
 *	    not well-suited for use with incomplete degree update or multiple
 *	    elimination.
 *
 * Authors, and Copyright (C) 2004 by:
 * Timothy A. Davis, Patrick Amestoy, Iain S. Duff, John K. Reid.
 *
 * Acknowledgements: This work (and the UMFPACK package) was supported by the
 * National Science Foundation (ASC-9111263, DMS-9223088, and CCR-0203270).
 * The UMFPACK/MA38 approximate degree update algorithm, the unsymmetric analog
 * which forms the basis of AMD, was developed while Tim Davis was supported by
 * CERFACS (Toulouse, France) in a post-doctoral position.  This C version, and
 * the etree postorder, were written while Tim Davis was on sabbatical at
 * Stanford University and Lawrence Berkeley National Laboratory.

 * ----------------------------------------------------------------------------
 * INPUT ARGUMENTS (unaltered):
 * ----------------------------------------------------------------------------

 * n:  The matrix order.  Restriction:  n >= 1.
 *
 * iwlen:  The size of the Iw array.  On input, the matrix is stored in
 *	Iw [0..pfree-1].  However, Iw [0..iwlen-1] should be slightly larger
 *	than what is required to hold the matrix, at least iwlen >= pfree + n.
 *	Otherwise, excessive compressions will take place.  The recommended
 *	value of iwlen is 1.2 * pfree + n, which is the value used in the
 *	user-callable interface to this routine (amd_order.c).  The algorithm
 *	will not run at all if iwlen < pfree.  Restriction: iwlen >= pfree + n.
 *	Note that this is slightly more restrictive than the actual minimum
 *	(iwlen >= pfree), but AMD_2 will be very slow with no elbow room.
 *	Thus, this routine enforces a bare minimum elbow room of size n.
 *
 * pfree: On input the tail end of the array, Iw [pfree..iwlen-1], is empty,
 *	and the matrix is stored in Iw [0..pfree-1].  During execution,
 *	additional data is placed in Iw, and pfree is modified so that
 *	Iw [pfree..iwlen-1] is always the unused part of Iw.
 *
 * Control:  A double array of size AMD_CONTROL containing input parameters
 *	that affect how the ordering is computed.  If NULL, then default
 *	settings are used.
 *
 *	Control [AMD_DENSE] is used to determine whether or not a given input
 *	row is "dense".  A row is "dense" if the number of entries in the row
 *	exceeds Control [AMD_DENSE] times sqrt (n), except that rows with 16 or
 *	fewer entries are never considered "dense".  To turn off the detection
 *	of dense rows, set Control [AMD_DENSE] to a negative number, or to a
 *	number larger than sqrt (n).  The default value of Control [AMD_DENSE]
 *	is AMD_DEFAULT_DENSE, which is defined in amd.h as 10.
 *
 *	Control [AMD_AGGRESSIVE] is used to determine whether or not aggressive
 *	absorption is to be performed.  If nonzero, then aggressive absorption
 *	is performed (this is the default).

 * ----------------------------------------------------------------------------
 * INPUT/OUPUT ARGUMENTS:
 * ----------------------------------------------------------------------------
 *
 * Pe:  An integer array of size n.  On input, Pe [i] is the index in Iw of
 *	the start of row i.  Pe [i] is ignored if row i has no off-diagonal
 *	entries.  Thus Pe [i] must be in the range 0 to pfree-1 for non-empty
 *	rows.
 *
 *	During execution, it is used for both supervariables and elements:
 *
 *	Principal supervariable i:  index into Iw of the description of
 *	    supervariable i.  A supervariable represents one or more rows of
 *	    the matrix with identical nonzero pattern.  In this case,
 *	    Pe [i] >= 0.
 *
 *	Non-principal supervariable i:  if i has been absorbed into another
 *	    supervariable j, then Pe [i] = FLIP (j), where FLIP (j) is defined
 *	    as (-(j)-2).  Row j has the same pattern as row i.  Note that j
 *	    might later be absorbed into another supervariable j2, in which
 *	    case Pe [i] is still FLIP (j), and Pe [j] = FLIP (j2) which is
 *	    < EMPTY, where EMPTY is defined as (-1) in amd_internal.h.
 *
 *	Unabsorbed element e:  the index into Iw of the description of element
 *	    e, if e has not yet been absorbed by a subsequent element.  Element
 *	    e is created when the supervariable of the same name is selected as
 *	    the pivot.  In this case, Pe [i] >= 0.
 *
 *	Absorbed element e:  if element e is absorbed into element e2, then
 *	    Pe [e] = FLIP (e2).  This occurs when the pattern of e (which we
 *	    refer to as Le) is found to be a subset of the pattern of e2 (that
 *	    is, Le2).  In this case, Pe [i] < EMPTY.  If element e is "null"
 *	    (it has no nonzeros outside its pivot block), then Pe [e] = EMPTY,
 *	    and e is the root of an assembly subtree (or the whole tree if
 *	    there is just one such root).
 *
 *	Dense variable i:  if i is "dense", then Pe [i] = EMPTY.
 *
 *	On output, Pe holds the assembly tree/forest, which implicitly
 *	represents a pivot order with identical fill-in as the actual order
 *	(via a depth-first search of the tree), as follows.  If Nv [i] > 0,
 *	then i represents a node in the assembly tree, and the parent of i is
 *	Pe [i], or EMPTY if i is a root.  If Nv [i] = 0, then (i, Pe [i])
 *	represents an edge in a subtree, the root of which is a node in the
 *	assembly tree.  Note that i refers to a row/column in the original
 *	matrix, not the permuted matrix.
 *
 * Info:  A double array of size AMD_INFO.  If present, (that is, not NULL),
 *	then statistics about the ordering are returned in the Info array.
 *	See amd.h for a description.

 * ----------------------------------------------------------------------------
 * INPUT/MODIFIED (undefined on output):
 * ----------------------------------------------------------------------------
 *
 * Len:  An integer array of size n.  On input, Len [i] holds the number of
 *	entries in row i of the matrix, excluding the diagonal.  The contents
 *	of Len are undefined on output.
 *
 * Iw:  An integer array of size iwlen.  On input, Iw [0..pfree-1] holds the
 *	description of each row i in the matrix.  The matrix must be symmetric,
 *	and both upper and lower triangular parts must be present.  The
 *	diagonal must not be present.  Row i is held as follows:
 *
 *	    Len [i]:  the length of the row i data structure in the Iw array.
 *	    Iw [Pe [i] ... Pe [i] + Len [i] - 1]:
 *		the list of column indices for nonzeros in row i (simple
 *		supervariables), excluding the diagonal.  All supervariables
 *		start with one row/column each (supervariable i is just row i).
 *		If Len [i] is zero on input, then Pe [i] is ignored on input.
 *
 *	    Note that the rows need not be in any particular order, and there
 *	    may be empty space between the rows.
 *
 *	During execution, the supervariable i experiences fill-in.  This is
 *	represented by placing in i a list of the elements that cause fill-in
 *	in supervariable i:
 *
 *	    Len [i]:  the length of supervariable i in the Iw array.
 *	    Iw [Pe [i] ... Pe [i] + Elen [i] - 1]:
 *		the list of elements that contain i.  This list is kept short
 *		by removing absorbed elements.
 *	    Iw [Pe [i] + Elen [i] ... Pe [i] + Len [i] - 1]:
 *		the list of supervariables in i.  This list is kept short by
 *		removing nonprincipal variables, and any entry j that is also
 *		contained in at least one of the elements (j in Le) in the list
 *		for i (e in row i).
 *
 *	When supervariable i is selected as pivot, we create an element e of
 *	the same name (e=i):
 *
 *	    Len [e]:  the length of element e in the Iw array.
 *	    Iw [Pe [e] ... Pe [e] + Len [e] - 1]:
 *		the list of supervariables in element e.
 *
 *	An element represents the fill-in that occurs when supervariable i is
 *	selected as pivot (which represents the selection of row i and all
 *	non-principal variables whose principal variable is i).  We use the
 *	term Le to denote the set of all supervariables in element e.  Absorbed
 *	supervariables and elements are pruned from these lists when
 *	computationally convenient.
 *
 *  CAUTION:  THE INPUT MATRIX IS OVERWRITTEN DURING COMPUTATION.
 *  The contents of Iw are undefined on output.

 * ----------------------------------------------------------------------------
 * OUTPUT (need not be set on input):
 * ----------------------------------------------------------------------------
 *
 * Nv:  An integer array of size n.  During execution, ABS (Nv [i]) is equal to
 *	the number of rows that are represented by the principal supervariable
 *	i.  If i is a nonprincipal or dense variable, then Nv [i] = 0.
 *	Initially, Nv [i] = 1 for all i.  Nv [i] < 0 signifies that i is a
 *	principal variable in the pattern Lme of the current pivot element me.
 *	After element me is constructed, Nv [i] is set back to a positive
 *	value.
 *
 *	On output, Nv [i] holds the number of pivots represented by super
 *	row/column i of the original matrix, or Nv [i] = 0 for non-principal
 *	rows/columns.  Note that i refers to a row/column in the original
 *	matrix, not the permuted matrix.
 *
 * Elen:  An integer array of size n.  See the description of Iw above.  At the
 *	start of execution, Elen [i] is set to zero for all rows i.  During
 *	execution, Elen [i] is the number of elements in the list for
 *	supervariable i.  When e becomes an element, Elen [e] = FLIP (esize) is
 *	set, where esize is the size of the element (the number of pivots, plus
 *	the number of nonpivotal entries).  Thus Elen [e] < EMPTY.
 *	Elen (i) = EMPTY set when variable i becomes nonprincipal.
 *
 *	For variables, Elen (i) >= EMPTY holds until just before the
 *	postordering and permutation vectors are computed.  For elements,
 *	Elen [e] < EMPTY holds.
 *
 *	On output, Elen [i] is the degree of the row/column in the Cholesky
 *	factorization of the permuted matrix, corresponding to the original row
 *	i, if i is a super row/column.  It is equal to EMPTY if i is
 *	non-principal.  Note that i refers to a row/column in the original
 *	matrix, not the permuted matrix.
 *
 *	Note that the contents of Elen on output differ from the Fortran
 *	version (Elen holds the inverse permutation in the Fortran version,
 *	which is instead returned in the Next array in this C version,
 *	described below).
 *
 * Last: In a degree list, Last [i] is the supervariable preceding i, or EMPTY
 *	if i is the head of the list.  In a hval bucket, Last [i] is the hval
 *	key for i.
 *
 *	Last [Head [hval]] is also used as the head of a hval bucket if
 *	Head [hval] contains a degree list (see the description of Head,
 *	below).
 *
 *	On output, Last [0..n-1] holds the permutation.  That is, if
 *	i = Last [k], then row i is the kth pivot row (where k ranges from 0 to
 *	n-1).  Row Last [k] of A is the kth row in the permuted matrix, PAP'.
 *
 * Next: Next [i] is the supervariable following i in a link list, or EMPTY if
 *	i is the last in the list.  Used for two kinds of lists:  degree lists
 *	and hval buckets (a supervariable can be in only one kind of list at a
 *	time).
 *
 *	On output Next [0..n-1] holds the inverse permutation. 	That is, if
 *	k = Next [i], then row i is the kth pivot row. Row i of A appears as
 *	the (Next[i])-th row in the permuted matrix, PAP'.
 *
 *	Note that the contents of Next on output differ from the Fortran
 *	version (Next is undefined on output in the Fortran version).

 * ----------------------------------------------------------------------------
 * LOCAL WORKSPACE (not input or output - used only during execution):
 * ----------------------------------------------------------------------------
 *
 * Degree:  An integer array of size n.  If i is a supervariable, then
 *	Degree [i] holds the current approximation of the external degree of
 *	row i (an upper bound).  The external degree is the number of nonzeros
 *	in row i, minus ABS (Nv [i]), the diagonal part.  The bound is equal to
 *	the exact external degree if Elen [i] is less than or equal to two.
 *
 *	We also use the term "external degree" for elements e to refer to
 *	|Le \ Lme|.  If e is an element, then Degree [e] is |Le|, which is the
 *	degree of the off-diagonal part of the element e (not including the
 *	diagonal part).
 *
 * Head:   An integer array of size n.  Head is used for degree lists.
 *	Head [deg] is the first supervariable in a degree list.  All
 *	supervariables i in a degree list Head [deg] have the same approximate
 *	degree, namely, deg = Degree [i].  If the list Head [deg] is empty then
 *	Head [deg] = EMPTY.
 *
 *	During supervariable detection Head [hval] also serves as a pointer to
 *	a hval bucket.  If Head [hval] >= 0, there is a degree list of degree
 *	hval.  The hval bucket head pointer is Last [Head [hval]].  If
 *	Head [hval] = EMPTY, then the degree list and hval bucket are both
 *	empty.  If Head [hval] < EMPTY, then the degree list is empty, and
 *	FLIP (Head [hval]) is the head of the hval bucket.  After supervariable
 *	detection is complete, all hval buckets are empty, and the
 *	(Last [Head [hval]] = EMPTY) condition is restored for the non-empty
 *	degree lists.
 *
 * W:  An integer array of size n.  The flag array W determines the status of
 *	elements and variables, and the external degree of elements.
 *
 *	for elements:
 *	    if W [e] = 0, then the element e is absorbed.
 *	    if W [e] >= wflg, then W [e] - wflg is the size of the set
 *		|Le \ Lme|, in terms of nonzeros (the sum of ABS (Nv [i]) for
 *		each principal variable i that is both in the pattern of
 *		element e and NOT in the pattern of the current pivot element,
 *		me).
 *	    if wflg > W [e] > 0, then e is not absorbed and has not yet been
 *		seen in the scan of the element lists in the computation of
 *		|Le\Lme| in Scan 1 below.
 *
 *	for variables:
 *	    during supervariable detection, if W [j] != wflg then j is
 *	    not in the pattern of variable i.
 *
 *	The W array is initialized by setting W [i] = 1 for all i, and by
 *	setting wflg = 2.  It is reinitialized if wflg becomes too large (to
 *	ensure that wflg+n does not cause integer overflow).

 * ----------------------------------------------------------------------------
 * LOCAL INTEGERS:
 * ----------------------------------------------------------------------------
 */

    Int deg, degme, dext, lemax, e, elenme, eln, i, ilast, inext, j,
	jlast, jnext, k, knt1, knt2, knt3, lenj, ln, me, mindeg, nel, nleft,
	nvi, nvj, nvpiv, slenme, wbig, we, wflg, wnvi, ok, ndense, ncmpa,
	dense, aggressive ;

    UInt hval ;	    /* unsigned, so that hval % n is well defined.*/

/*
 * deg:		the degree of a variable or element
 * degme:	size, |Lme|, of the current element, me (= Degree [me])
 * dext:	external degree, |Le \ Lme|, of some element e
 * lemax:	largest |Le| seen so far (called dmax in Fortran version)
 * e:		an element
 * elenme:	the length, Elen [me], of element list of pivotal variable
 * eln:		the length, Elen [...], of an element list
 * hval:	the computed value of the hval function
 * i:		a supervariable
 * ilast:	the entry in a link list preceding i
 * inext:	the entry in a link list following i
 * j:		a supervariable
 * jlast:	the entry in a link list preceding j
 * jnext:	the entry in a link list, or path, following j
 * k:		the pivot order of an element or variable
 * knt1:	loop counter used during element construction
 * knt2:	loop counter used during element construction
 * knt3:	loop counter used during compression
 * lenj:	Len [j]
 * ln:		length of a supervariable list
 * me:		current supervariable being eliminated, and the current
 *		    element created by eliminating that supervariable
 * mindeg:	current minimum degree
 * nel:		number of pivots selected so far
 * nleft:	n - nel, the number of nonpivotal rows/columns remaining
 * nvi:		the number of variables in a supervariable i (= Nv [i])
 * nvj:		the number of variables in a supervariable j (= Nv [j])
 * nvpiv:	number of pivots in current element
 * slenme:	number of variables in variable list of pivotal variable
 * wbig:	= (INT32_MAX - n) for the int32_t version, (INT64_MAX - n)
 *                  for the int64_t version.  wflg is not allowed to
 *                  be >= wbig.
 * we:		W [e]
 * wflg:	used for flagging the W array.  See description of Iw.
 * wnvi:	wflg - Nv [i]
 * x:		either a supervariable or an element
 *
 * ok:		true if supervariable j can be absorbed into i
 * ndense:	number of "dense" rows/columns
 * dense:	rows/columns with initial degree > dense are considered "dense"
 * aggressive:	true if aggressive absorption is being performed
 * ncmpa:	number of garbage collections

 * ----------------------------------------------------------------------------
 * LOCAL DOUBLES, used for statistical output only (except for alpha):
 * ----------------------------------------------------------------------------
 */

    double f, r, ndiv, s, nms_lu, nms_ldl, dmax, alpha, lnz, lnzme ;

/*
 * f:		nvpiv
 * r:		degme + nvpiv
 * ndiv:	number of divisions for LU or LDL' factorizations
 * s:		number of multiply-subtract pairs for LU factorization, for the
 *		    current element me
 * nms_lu	number of multiply-subtract pairs for LU factorization
 * nms_ldl	number of multiply-subtract pairs for LDL' factorization
 * dmax:	the largest number of entries in any column of L, including the
 *		    diagonal
 * alpha:	"dense" degree ratio
 * lnz:		the number of nonzeros in L (excluding the diagonal)
 * lnzme:	the number of nonzeros in L (excl. the diagonal) for the
 *		    current element me

 * ----------------------------------------------------------------------------
 * LOCAL "POINTERS" (indices into the Iw array)
 * ----------------------------------------------------------------------------
*/

    Int p, p1, p2, p3, p4, pdst, pend, pj, pme, pme1, pme2, pn, psrc ;

/*
 * Any parameter (Pe [...] or pfree) or local variable starting with "p" (for
 * Pointer) is an index into Iw, and all indices into Iw use variables starting
 * with "p."  The only exception to this rule is the iwlen input argument.
 *
 * p:           pointer into lots of things
 * p1:          Pe [i] for some variable i (start of element list)
 * p2:          Pe [i] + Elen [i] -  1 for some variable i
 * p3:          index of first supervariable in clean list
 * p4:		
 * pdst:        destination pointer, for compression
 * pend:        end of memory to compress
 * pj:          pointer into an element or variable
 * pme:         pointer into the current element (pme1...pme2)
 * pme1:        the current element, me, is stored in Iw [pme1...pme2]
 * pme2:        the end of the current element
 * pn:          pointer into a "clean" variable, also used to compress
 * psrc:        source pointer, for compression
*/

/* ========================================================================= */
/*  INITIALIZATIONS */
/* ========================================================================= */

    /* Note that this restriction on iwlen is slightly more restrictive than
     * what is actually required in AMD_2.  AMD_2 can operate with no elbow
     * room at all, but it will be slow.  For better performance, at least
     * size-n elbow room is enforced. */

    /* initialize output statistics */
    lnz = 0 ;
    ndiv = 0 ;
    nms_lu = 0 ;
    nms_ldl = 0 ;
    dmax = 1 ;
    me = EMPTY ;

    mindeg = 0 ;
    ncmpa = 0 ;
    nel = 0 ;
    lemax = 0 ;

    /* get control parameters */
    if (Control != (double *) NULL)
    {
	alpha = Control [AMD_DENSE] ;
	aggressive = (Control [AMD_AGGRESSIVE] != 0) ;
    }
    else
    {
	alpha = AMD_DEFAULT_DENSE ;
	aggressive = AMD_DEFAULT_AGGRESSIVE ;
    }
    /* Note: if alpha is NaN, this is undefined: */
    if (alpha < 0)
    {
	/* only remove completely dense rows/columns */
	dense = n-2 ;
    }
    else
    {
	dense = alpha * sqrt ((double) n) ;
    }
    dense = MAX (16, dense) ;
    dense = MIN (n,  dense) ;

    for (i = 0 ; i < n ; i++)
    {
	Last [i] = EMPTY ;
	Head [i] = EMPTY ;
	Next [i] = EMPTY ;
	/* if separate Hhead array is used for hval buckets: *
	Hhead [i] = EMPTY ;
	*/
	Nv [i] = 1 ;
	W [i] = 1 ;
	Elen [i] = 0 ;
	Degree [i] = Len [i] ;
    }

#ifndef NDEBUG
#endif

    /* initialize wflg */
    wbig = Int_MAX_VAL - n ;
    wflg = clear_flag (0, wbig, W, n) ;

    /* --------------------------------------------------------------------- */
    /* initialize degree lists and eliminate dense and empty rows */
    /* --------------------------------------------------------------------- */

    ndense = 0 ;

    for (i = 0 ; i < n ; i++)
    {
	deg = Degree [i] ;
	if (deg == 0)
	{

	    /* -------------------------------------------------------------
	     * we have a variable that can be eliminated at once because
	     * there is no off-diagonal non-zero in its row.  Note that
	     * Nv [i] = 1 for an empty variable i.  It is treated just
	     * the same as an eliminated element i.
	     * ------------------------------------------------------------- */

	    Elen [i] = FLIP (1) ;
	    nel++ ;
	    Pe [i] = EMPTY ;
	    W [i] = 0 ;

	}
	else if (deg > dense)
	{

	    /* -------------------------------------------------------------
	     * Dense variables are not treated as elements, but as unordered,
	     * non-principal variables that have no parent.  They do not take
	     * part in the postorder, since Nv [i] = 0.  Note that the Fortran
	     * version does not have this option.
	     * ------------------------------------------------------------- */

	    ndense++ ;
	    Nv [i] = 0 ;		/* do not postorder this node */
	    Elen [i] = EMPTY ;
	    nel++ ;
	    Pe [i] = EMPTY ;

	}
	else
	{

	    /* -------------------------------------------------------------
	     * place i in the degree list corresponding to its degree
	     * ------------------------------------------------------------- */

	    inext = Head [deg] ;
	    if (inext != EMPTY) Last [inext] = i ;
	    Next [i] = inext ;
	    Head [deg] = i ;

	}
    }

/* ========================================================================= */
/* WHILE (selecting pivots) DO */
/* ========================================================================= */

    while (nel < n)
    {

#ifndef NDEBUG
	{
	}
#endif

/* ========================================================================= */
/* GET PIVOT OF MINIMUM DEGREE */
/* ========================================================================= */

	/* ----------------------------------------------------------------- */
	/* find next supervariable for elimination */
	/* ----------------------------------------------------------------- */

	for (deg = mindeg ; deg < n ; deg++)
	{
	    me = Head [deg] ;
	    if (me != EMPTY) break ;
	}
	mindeg = deg ;

	/* ----------------------------------------------------------------- */
	/* remove chosen variable from link list */
	/* ----------------------------------------------------------------- */

	inext = Next [me] ;
	if (inext != EMPTY) Last [inext] = EMPTY ;
	Head [deg] = inext ;

	/* ----------------------------------------------------------------- */
	/* me represents the elimination of pivots nel to nel+Nv[me]-1. */
	/* place me itself as the first in this set. */
	/* ----------------------------------------------------------------- */

	elenme = Elen [me] ;
	nvpiv = Nv [me] ;
	nel += nvpiv ;

/* ========================================================================= */
/* CONSTRUCT NEW ELEMENT */
/* ========================================================================= */

	/* -----------------------------------------------------------------
	 * At this point, me is the pivotal supervariable.  It will be
	 * converted into the current element.  Scan list of the pivotal
	 * supervariable, me, setting tree pointers and constructing new list
	 * of supervariables for the new element, me.  p is a pointer to the
	 * current position in the old list.
	 * ----------------------------------------------------------------- */

	/* flag the variable "me" as being in Lme by negating Nv [me] */
	Nv [me] = -nvpiv ;
	degme = 0 ;

	if (elenme == 0)
	{

	    /* ------------------------------------------------------------- */
	    /* construct the new element in place */
	    /* ------------------------------------------------------------- */

	    pme1 = Pe [me] ;
	    pme2 = pme1 - 1 ;

	    for (p = pme1 ; p <= pme1 + Len [me] - 1 ; p++)
	    {
		i = Iw [p] ;
		nvi = Nv [i] ;
		if (nvi > 0)
		{

		    /* ----------------------------------------------------- */
		    /* i is a principal variable not yet placed in Lme. */
		    /* store i in new list */
		    /* ----------------------------------------------------- */

		    /* flag i as being in Lme by negating Nv [i] */
		    degme += nvi ;
		    Nv [i] = -nvi ;
		    Iw [++pme2] = i ;

		    /* ----------------------------------------------------- */
		    /* remove variable i from degree list. */
		    /* ----------------------------------------------------- */

		    ilast = Last [i] ;
		    inext = Next [i] ;
		    if (inext != EMPTY) Last [inext] = ilast ;
		    if (ilast != EMPTY)
		    {
			Next [ilast] = inext ;
		    }
		    else
		    {
			/* i is at the head of the degree list */
			Head [Degree [i]] = inext ;
		    }
		}
	    }
	}
	else
	{

	    /* ------------------------------------------------------------- */
	    /* construct the new element in empty space, Iw [pfree ...] */
	    /* ------------------------------------------------------------- */

	    p = Pe [me] ;
	    pme1 = pfree ;
	    slenme = Len [me] - elenme ;

	    for (knt1 = 1 ; knt1 <= elenme + 1 ; knt1++)
	    {

		if (knt1 > elenme)
		{
		    /* search the supervariables in me. */
		    e = me ;
		    pj = p ;
		    ln = slenme ;
		}
		else
		{
		    /* search the elements in me. */
		    e = Iw [p++] ;
		    pj = Pe [e] ;
		    ln = Len [e] ;
		}

		/* ---------------------------------------------------------
		 * search for different supervariables and add them to the
		 * new list, compressing when necessary. this loop is
		 * executed once for each element in the list and once for
		 * all the supervariables in the list.
		 * --------------------------------------------------------- */

		for (knt2 = 1 ; knt2 <= ln ; knt2++)
		{
		    i = Iw [pj++] ;
		    nvi = Nv [i] ;

		    if (nvi > 0)
		    {

			/* ------------------------------------------------- */
			/* compress Iw, if necessary */
			/* ------------------------------------------------- */

			if (pfree >= iwlen)
			{


			    /* prepare for compressing Iw by adjusting pointers
			     * and lengths so that the lists being searched in
			     * the inner and outer loops contain only the
			     * remaining entries. */

			    Pe [me] = p ;
			    Len [me] -= knt1 ;
			    /* check if nothing left of supervariable me */
			    if (Len [me] == 0) Pe [me] = EMPTY ;
			    Pe [e] = pj ;
			    Len [e] = ln - knt2 ;
			    /* nothing left of element e */
			    if (Len [e] == 0) Pe [e] = EMPTY ;

			    ncmpa++ ;	/* one more garbage collection */

			    /* store first entry of each object in Pe */
			    /* FLIP the first entry in each object */
			    for (j = 0 ; j < n ; j++)
			    {
				pn = Pe [j] ;
				if (pn >= 0)
				{
				    Pe [j] = Iw [pn] ;
				    Iw [pn] = FLIP (j) ;
				}
			    }

			    /* psrc/pdst point to source/destination */
			    psrc = 0 ;
			    pdst = 0 ;
			    pend = pme1 - 1 ;

			    while (psrc <= pend)
			    {
				/* search for next FLIP'd entry */
				j = FLIP (Iw [psrc++]) ;
				if (j >= 0)
				{
				    Iw [pdst] = Pe [j] ;
				    Pe [j] = pdst++ ;
				    lenj = Len [j] ;
				    /* copy from source to destination */
				    for (knt3 = 0 ; knt3 <= lenj - 2 ; knt3++)
				    {
					Iw [pdst++] = Iw [psrc++] ;
				    }
				}
			    }

			    /* move the new partially-constructed element */
			    p1 = pdst ;
			    for (psrc = pme1 ; psrc <= pfree-1 ; psrc++)
			    {
				Iw [pdst++] = Iw [psrc] ;
			    }
			    pme1 = p1 ;
			    pfree = pdst ;
			    pj = Pe [e] ;
			    p = Pe [me] ;

			}

			/* ------------------------------------------------- */
			/* i is a principal variable not yet placed in Lme */
			/* store i in new list */
			/* ------------------------------------------------- */

			/* flag i as being in Lme by negating Nv [i] */
			degme += nvi ;
			Nv [i] = -nvi ;
			Iw [pfree++] = i ;

			/* ------------------------------------------------- */
			/* remove variable i from degree link list */
			/* ------------------------------------------------- */

			ilast = Last [i] ;
			inext = Next [i] ;
			if (inext != EMPTY) Last [inext] = ilast ;
			if (ilast != EMPTY)
			{
			    Next [ilast] = inext ;
			}
			else
			{
			    /* i is at the head of the degree list */
			    Head [Degree [i]] = inext ;
			}
		    }
		}

		if (e != me)
		{
		    /* set tree pointer and flag to indicate element e is
		     * absorbed into new element me (the parent of e is me) */
		    Pe [e] = FLIP (me) ;
		    W [e] = 0 ;
		}
	    }

	    pme2 = pfree - 1 ;
	}

	/* ----------------------------------------------------------------- */
	/* me has now been converted into an element in Iw [pme1..pme2] */
	/* ----------------------------------------------------------------- */

	/* degme holds the external degree of new element */
	Degree [me] = degme ;
	Pe [me] = pme1 ;
	Len [me] = pme2 - pme1 + 1 ;

	Elen [me] = FLIP (nvpiv + degme) ;
	/* FLIP (Elen (me)) is now the degree of pivot (including
	 * diagonal part). */

#ifndef NDEBUG
#endif

	/* ----------------------------------------------------------------- */
	/* make sure that wflg is not too large. */
	/* ----------------------------------------------------------------- */

	/* With the current value of wflg, wflg+n must not cause integer
	 * overflow */

	wflg = clear_flag (wflg, wbig, W, n) ;

/* ========================================================================= */
/* COMPUTE (W [e] - wflg) = |Le\Lme| FOR ALL ELEMENTS */
/* ========================================================================= */

	/* -----------------------------------------------------------------
	 * Scan 1:  compute the external degrees of previous elements with
	 * respect to the current element.  That is:
	 *       (W [e] - wflg) = |Le \ Lme|
	 * for each element e that appears in any supervariable in Lme.  The
	 * notation Le refers to the pattern (list of supervariables) of a
	 * previous element e, where e is not yet absorbed, stored in
	 * Iw [Pe [e] + 1 ... Pe [e] + Len [e]].  The notation Lme
	 * refers to the pattern of the current element (stored in
	 * Iw [pme1..pme2]).   If aggressive absorption is enabled, and
	 * (W [e] - wflg) becomes zero, then the element e will be absorbed
	 * in Scan 2.
	 * ----------------------------------------------------------------- */

	for (pme = pme1 ; pme <= pme2 ; pme++)
	{
	    i = Iw [pme] ;
	    eln = Elen [i] ;
	    if (eln > 0)
	    {
		/* note that Nv [i] has been negated to denote i in Lme: */
		nvi = -Nv [i] ;
		wnvi = wflg - nvi ;
		for (p = Pe [i] ; p <= Pe [i] + eln - 1 ; p++)
		{
		    e = Iw [p] ;
		    we = W [e] ;
		    if (we >= wflg)
		    {
			/* unabsorbed element e has been seen in this loop */
			we -= nvi ;
		    }
		    else if (we != 0)
		    {
			/* e is an unabsorbed element */
			/* this is the first we have seen e in all of Scan 1 */
			we = Degree [e] + wnvi ;
		    }
		    W [e] = we ;
		}
	    }
	}

/* ========================================================================= */
/* DEGREE UPDATE AND ELEMENT ABSORPTION */
/* ========================================================================= */

	/* -----------------------------------------------------------------
	 * Scan 2:  for each i in Lme, sum up the degree of Lme (which is
	 * degme), plus the sum of the external degrees of each Le for the
	 * elements e appearing within i, plus the supervariables in i.
	 * Place i in hval list.
	 * ----------------------------------------------------------------- */

	for (pme = pme1 ; pme <= pme2 ; pme++)
	{
	    i = Iw [pme] ;
	    p1 = Pe [i] ;
	    p2 = p1 + Elen [i] - 1 ;
	    pn = p1 ;
	    hval = 0 ;
	    deg = 0 ;

	    /* ------------------------------------------------------------- */
	    /* scan the element list associated with supervariable i */
	    /* ------------------------------------------------------------- */

	    /* UMFPACK/MA38-style approximate degree: */
	    if (aggressive)
	    {
		for (p = p1 ; p <= p2 ; p++)
		{
		    e = Iw [p] ;
		    we = W [e] ;
		    if (we != 0)
		    {
			/* e is an unabsorbed element */
			/* dext = | Le \ Lme | */
			dext = we - wflg ;
			if (dext > 0)
			{
			    deg += dext ;
			    Iw [pn++] = e ;
			    hval += e ;
			}
			else
			{
			    /* external degree of e is zero, absorb e into me*/
			    Pe [e] = FLIP (me) ;
			    W [e] = 0 ;
			}
		    }
		}
	    }
	    else
	    {
		for (p = p1 ; p <= p2 ; p++)
		{
		    e = Iw [p] ;
		    we = W [e] ;
		    if (we != 0)
		    {
			/* e is an unabsorbed element */
			dext = we - wflg ;
			deg += dext ;
			Iw [pn++] = e ;
			hval += e ;
		    }
		}
	    }

	    /* count the number of elements in i (including me): */
	    Elen [i] = pn - p1 + 1 ;

	    /* ------------------------------------------------------------- */
	    /* scan the supervariables in the list associated with i */
	    /* ------------------------------------------------------------- */

	    /* The bulk of the AMD run time is typically spent in this loop,
	     * particularly if the matrix has many dense rows that are not
	     * removed prior to ordering. */
	    p3 = pn ;
	    p4 = p1 + Len [i] ;
	    for (p = p2 + 1 ; p < p4 ; p++)
	    {
		j = Iw [p] ;
		nvj = Nv [j] ;
		if (nvj > 0)
		{
		    /* j is unabsorbed, and not in Lme. */
		    /* add to degree and add to new list */
		    deg += nvj ;
		    Iw [pn++] = j ;
		    hval += j ;
		}
	    }

	    /* ------------------------------------------------------------- */
	    /* update the degree and check for mass elimination */
	    /* ------------------------------------------------------------- */

	    /* with aggressive absorption, deg==0 is identical to the
	     * Elen [i] == 1 && p3 == pn test, below. */

	    if (Elen [i] == 1 && p3 == pn)
	    {

		/* --------------------------------------------------------- */
		/* mass elimination */
		/* --------------------------------------------------------- */

		/* There is nothing left of this node except for an edge to
		 * the current pivot element.  Elen [i] is 1, and there are
		 * no variables adjacent to node i.  Absorb i into the
		 * current pivot element, me.  Note that if there are two or
		 * more mass eliminations, fillin due to mass elimination is
		 * possible within the nvpiv-by-nvpiv pivot block.  It is this
		 * step that causes AMD's analysis to be an upper bound.
		 *
		 * The reason is that the selected pivot has a lower
		 * approximate degree than the true degree of the two mass
		 * eliminated nodes.  There is no edge between the two mass
		 * eliminated nodes.  They are merged with the current pivot
		 * anyway.
		 *
		 * No fillin occurs in the Schur complement, in any case,
		 * and this effect does not decrease the quality of the
		 * ordering itself, just the quality of the nonzero and
		 * flop count analysis.  It also means that the post-ordering
		 * is not an exact elimination tree post-ordering. */

		Pe [i] = FLIP (me) ;
		nvi = -Nv [i] ;
		degme -= nvi ;
		nvpiv += nvi ;
		nel += nvi ;
		Nv [i] = 0 ;
		Elen [i] = EMPTY ;

	    }
	    else
	    {

		/* --------------------------------------------------------- */
		/* update the upper-bound degree of i */
		/* --------------------------------------------------------- */

		/* the following degree does not yet include the size
		 * of the current element, which is added later: */

		Degree [i] = MIN (Degree [i], deg) ;

		/* --------------------------------------------------------- */
		/* add me to the list for i */
		/* --------------------------------------------------------- */

		/* move first supervariable to end of list */
		Iw [pn] = Iw [p3] ;
		/* move first element to end of element part of list */
		Iw [p3] = Iw [p1] ;
		/* add new element, me, to front of list. */
		Iw [p1] = me ;
		/* store the new length of the list in Len [i] */
		Len [i] = pn - p1 + 1 ;

		/* --------------------------------------------------------- */
		/* place in hval bucket.  Save hval key of i in Last [i]. */
		/* --------------------------------------------------------- */

		/* NOTE: this can fail if hval is negative, because the ANSI C
		 * standard does not define a % b when a and/or b are negative.
		 * That's why hval is defined as an unsigned Int, to avoid this
		 * problem. */
		hval = hval % n ;

		/* if the Hhead array is not used: */
		j = Head [hval] ;
		if (j <= EMPTY)
		{
		    /* degree list is empty, hval head is FLIP (j) */
		    Next [i] = FLIP (j) ;
		    Head [hval] = FLIP (i) ;
		}
		else
		{
		    /* degree list is not empty, use Last [Head [hval]] as
		     * hval head. */
		    Next [i] = Last [j] ;
		    Last [j] = i ;
		}

		/* if a separate Hhead array is used: *
		Next [i] = Hhead [hval] ;
		Hhead [hval] = i ;
		*/

		Last [i] = hval ;
	    }
	}

	Degree [me] = degme ;

	/* ----------------------------------------------------------------- */
	/* Clear the counter array, W [...], by incrementing wflg. */
	/* ----------------------------------------------------------------- */

	/* make sure that wflg+n does not cause integer overflow */
	lemax =  MAX (lemax, degme) ;
	wflg += lemax ;
	wflg = clear_flag (wflg, wbig, W, n) ;
	/*  at this point, W [0..n-1] < wflg holds */

/* ========================================================================= */
/* SUPERVARIABLE DETECTION */
/* ========================================================================= */

	for (pme = pme1 ; pme <= pme2 ; pme++)
	{
	    i = Iw [pme] ;
	    if (Nv [i] < 0)
	    {
		/* i is a principal variable in Lme */

		/* ---------------------------------------------------------
		 * examine all hval buckets with 2 or more variables.  We do
		 * this by examing all unique hval keys for supervariables in
		 * the pattern Lme of the current element, me
		 * --------------------------------------------------------- */

		/* let i = head of hval bucket, and empty the hval bucket */
		hval = Last [i] ;

		/* if Hhead array is not used: */
		j = Head [hval] ;
		if (j == EMPTY)
		{
		    /* hval bucket and degree list are both empty */
		    i = EMPTY ;
		}
		else if (j < EMPTY)
		{
		    /* degree list is empty */
		    i = FLIP (j) ;
		    Head [hval] = EMPTY ;
		}
		else
		{
		    /* degree list is not empty, restore Last [j] of head j */
		    i = Last [j] ;
		    Last [j] = EMPTY ;
		}

		/* if separate Hhead array is used: *
		i = Hhead [hval] ;
		Hhead [hval] = EMPTY ;
		*/


		while (i != EMPTY && Next [i] != EMPTY)
		{

		    /* -----------------------------------------------------
		     * this bucket has one or more variables following i.
		     * scan all of them to see if i can absorb any entries
		     * that follow i in hval bucket.  Scatter i into w.
		     * ----------------------------------------------------- */

		    ln = Len [i] ;
		    eln = Elen [i] ;
		    /* do not flag the first element in the list (me) */
		    for (p = Pe [i] + 1 ; p <= Pe [i] + ln - 1 ; p++)
		    {
			W [Iw [p]] = wflg ;
		    }

		    /* ----------------------------------------------------- */
		    /* scan every other entry j following i in bucket */
		    /* ----------------------------------------------------- */

		    jlast = i ;
		    j = Next [i] ;

		    while (j != EMPTY)
		    {
			/* ------------------------------------------------- */
			/* check if j and i have identical nonzero pattern */
			/* ------------------------------------------------- */


			/* check if i and j have the same Len and Elen */
			ok = (Len [j] == ln) && (Elen [j] == eln) ;
			/* skip the first element in the list (me) */
			for (p = Pe [j] + 1 ; ok && p <= Pe [j] + ln - 1 ; p++)
			{
			    if (W [Iw [p]] != wflg) ok = 0 ;
			}
			if (ok)
			{
			    /* --------------------------------------------- */
			    /* found it!  j can be absorbed into i */
			    /* --------------------------------------------- */

			    Pe [j] = FLIP (i) ;
			    /* both Nv [i] and Nv [j] are negated since they */
			    /* are in Lme, and the absolute values of each */
			    /* are the number of variables in i and j: */
			    Nv [i] += Nv [j] ;
			    Nv [j] = 0 ;
			    Elen [j] = EMPTY ;
			    /* delete j from hval bucket */
			    j = Next [j] ;
			    Next [jlast] = j ;

			}
			else
			{
			    /* j cannot be absorbed into i */
			    jlast = j ;
			    j = Next [j] ;
			}
		    }

		    /* -----------------------------------------------------
		     * no more variables can be absorbed into i
		     * go to next i in bucket and clear flag array
		     * ----------------------------------------------------- */

		    wflg++ ;
		    i = Next [i] ;

		}
	    }
	}

/* ========================================================================= */
/* RESTORE DEGREE LISTS AND REMOVE NONPRINCIPAL SUPERVARIABLES FROM ELEMENT */
/* ========================================================================= */

	p = pme1 ;
	nleft = n - nel ;
	for (pme = pme1 ; pme <= pme2 ; pme++)
	{
	    i = Iw [pme] ;
	    nvi = -Nv [i] ;
	    if (nvi > 0)
	    {
		/* i is a principal variable in Lme */
		/* restore Nv [i] to signify that i is principal */
		Nv [i] = nvi ;

		/* --------------------------------------------------------- */
		/* compute the external degree (add size of current element) */
		/* --------------------------------------------------------- */

		deg = Degree [i] + degme - nvi ;
		deg = MIN (deg, nleft - nvi) ;

		/* --------------------------------------------------------- */
		/* place the supervariable at the head of the degree list */
		/* --------------------------------------------------------- */

		inext = Head [deg] ;
		if (inext != EMPTY) Last [inext] = i ;
		Next [i] = inext ;
		Last [i] = EMPTY ;
		Head [deg] = i ;

		/* --------------------------------------------------------- */
		/* save the new degree, and find the minimum degree */
		/* --------------------------------------------------------- */

		mindeg = MIN (mindeg, deg) ;
		Degree [i] = deg ;

		/* --------------------------------------------------------- */
		/* place the supervariable in the element pattern */
		/* --------------------------------------------------------- */

		Iw [p++] = i ;

	    }
	}

/* ========================================================================= */
/* FINALIZE THE NEW ELEMENT */
/* ========================================================================= */

	Nv [me] = nvpiv ;
	/* save the length of the list for the new element me */
	Len [me] = p - pme1 ;
	if (Len [me] == 0)
	{
	    /* there is nothing left of the current pivot element */
	    /* it is a root of the assembly tree */
	    Pe [me] = EMPTY ;
	    W [me] = 0 ;
	}
	if (elenme != 0)
	{
	    /* element was not constructed in place: deallocate part of */
	    /* it since newly nonprincipal variables may have been removed */
	    pfree = p ;
	}

	/* The new element has nvpiv pivots and the size of the contribution
	 * block for a multifrontal method is degme-by-degme, not including
	 * the "dense" rows/columns.  If the "dense" rows/columns are included,
	 * the frontal matrix is no larger than
	 * (degme+ndense)-by-(degme+ndense).
	 */

	if (Info != (double *) NULL)
	{
	    f = nvpiv ;
	    r = degme + ndense ;
	    dmax = MAX (dmax, f + r) ;

	    /* number of nonzeros in L (excluding the diagonal) */
	    lnzme = f*r + (f-1)*f/2 ;
	    lnz += lnzme ;

	    /* number of divide operations for LDL' and for LU */
	    ndiv += lnzme ;

	    /* number of multiply-subtract pairs for LU */
	    s = f*r*r + r*(f-1)*f + (f-1)*f*(2*f-1)/6 ;
	    nms_lu += s ;

	    /* number of multiply-subtract pairs for LDL' */
	    nms_ldl += (s + lnzme)/2 ;
	}

#ifndef NDEBUG
	for (pme = Pe [me] ; pme <= Pe [me] + Len [me] - 1 ; pme++)
	{
	}
#endif

    }

/* ========================================================================= */
/* DONE SELECTING PIVOTS */
/* ========================================================================= */

    if (Info != (double *) NULL)
    {

	/* count the work to factorize the ndense-by-ndense submatrix */
	f = ndense ;
	dmax = MAX (dmax, (double) ndense) ;

	/* number of nonzeros in L (excluding the diagonal) */
	lnzme = (f-1)*f/2 ;
	lnz += lnzme ;

	/* number of divide operations for LDL' and for LU */
	ndiv += lnzme ;

	/* number of multiply-subtract pairs for LU */
	s = (f-1)*f*(2*f-1)/6 ;
	nms_lu += s ;

	/* number of multiply-subtract pairs for LDL' */
	nms_ldl += (s + lnzme)/2 ;

	/* number of nz's in L (excl. diagonal) */
	Info [AMD_LNZ] = lnz ;

	/* number of divide ops for LU and LDL' */
	Info [AMD_NDIV] = ndiv ;

	/* number of multiply-subtract pairs for LDL' */
	Info [AMD_NMULTSUBS_LDL] = nms_ldl ;

	/* number of multiply-subtract pairs for LU */
	Info [AMD_NMULTSUBS_LU] = nms_lu ;

	/* number of "dense" rows/columns */
	Info [AMD_NDENSE] = ndense ;

	/* largest front is dmax-by-dmax */
	Info [AMD_DMAX] = dmax ;

	/* number of garbage collections in AMD */
	Info [AMD_NCMPA] = ncmpa ;

	/* successful ordering */
	Info [AMD_STATUS] = AMD_OK ;
    }

/* ========================================================================= */
/* POST-ORDERING */
/* ========================================================================= */

/* -------------------------------------------------------------------------
 * Variables at this point:
 *
 * Pe: holds the elimination tree.  The parent of j is FLIP (Pe [j]),
 *	or EMPTY if j is a root.  The tree holds both elements and
 *	non-principal (unordered) variables absorbed into them.
 *	Dense variables are non-principal and unordered.
 *
 * Elen: holds the size of each element, including the diagonal part.
 *	FLIP (Elen [e]) > 0 if e is an element.  For unordered
 *	variables i, Elen [i] is EMPTY.
 *
 * Nv: Nv [e] > 0 is the number of pivots represented by the element e.
 *	For unordered variables i, Nv [i] is zero.
 *
 * Contents no longer needed:
 *	W, Iw, Len, Degree, Head, Next, Last.
 *
 * The matrix itself has been destroyed.
 *
 * n: the size of the matrix.
 * No other scalars needed (pfree, iwlen, etc.)
 * ------------------------------------------------------------------------- */

    /* restore Pe */
    for (i = 0 ; i < n ; i++)
    {
	Pe [i] = FLIP (Pe [i]) ;
    }

    /* restore Elen, for output information, and for postordering */
    for (i = 0 ; i < n ; i++)
    {
	Elen [i] = FLIP (Elen [i]) ;
    }

/* Now the parent of j is Pe [j], or EMPTY if j is a root.  Elen [e] > 0
 * is the size of element e.  Elen [i] is EMPTY for unordered variable i. */

#ifndef NDEBUG
    for (i = 0 ; i < n ; i++)
    {
	if (Nv [i] > 0)
	{
	    /* this is an element */
	    e = i ;
	}
    }
    for (e = 0 ; e < n ; e++)
    {
	if (Nv [e] > 0)
	{
	}
    }
    for (i = 0 ; i < n ; i++)
    {
	Int cnt ;
	if (Nv [i] == 0)
	{
	    j = Pe [i] ;
	    cnt = 0 ;
	    if (j == EMPTY)
	    {
	    }
	    else
	    {
		while (Nv [j] == 0)
		{
		    j = Pe [j] ;
		    cnt++ ;
		    if (cnt > n) break ;
		}
		e = j ;
	    }
	}
    }
#endif

/* ========================================================================= */
/* compress the paths of the variables */
/* ========================================================================= */

    for (i = 0 ; i < n ; i++)
    {
	if (Nv [i] == 0)
	{

	    /* -------------------------------------------------------------
	     * i is an un-ordered row.  Traverse the tree from i until
	     * reaching an element, e.  The element, e, was the principal
	     * supervariable of i and all nodes in the path from i to when e
	     * was selected as pivot.
	     * ------------------------------------------------------------- */

	    j = Pe [i] ;
	    if (j == EMPTY)
	    {
		/* Skip a dense variable.  It has no parent. */
		continue ;
	    }

	    /* while (j is a variable) */
	    while (Nv [j] == 0)
	    {
		j = Pe [j] ;
	    }
	    /* got to an element e */
	    e = j ;

	    /* -------------------------------------------------------------
	     * traverse the path again from i to e, and compress the path
	     * (all nodes point to e).  Path compression allows this code to
	     * compute in O(n) time.
	     * ------------------------------------------------------------- */

	    j = i ;
	    /* while (j is a variable) */
	    while (Nv [j] == 0)
	    {
		jnext = Pe [j] ;
		Pe [j] = e ;
		j = jnext ;
	    }
	}
    }

/* ========================================================================= */
/* postorder the assembly tree */
/* ========================================================================= */

    AMD_postorder (n, Pe, Nv, Elen,
	W,			/* output order */
	Head, Next, Last) ;	/* workspace */

/* ========================================================================= */
/* compute output permutation and inverse permutation */
/* ========================================================================= */

    /* W [e] = k means that element e is the kth element in the new
     * order.  e is in the range 0 to n-1, and k is in the range 0 to
     * the number of elements.  Use Head for inverse order. */

    for (k = 0 ; k < n ; k++)
    {
	Head [k] = EMPTY ;
	Next [k] = EMPTY ;
    }
    for (e = 0 ; e < n ; e++)
    {
	k = W [e] ;
	if (k != EMPTY)
	{
	    Head [k] = e ;
	}
    }

    /* construct output inverse permutation in Next,
     * and permutation in Last */
    nel = 0 ;
    for (k = 0 ; k < n ; k++)
    {
	e = Head [k] ;
	if (e == EMPTY) break ;
	Next [e] = nel ;
	nel += Nv [e] ;
    }

    /* order non-principal variables (dense, & those merged into supervar's) */
    for (i = 0 ; i < n ; i++)
    {
	if (Nv [i] == 0)
	{
	    e = Pe [i] ;
	    if (e != EMPTY)
	    {
		/* This is an unordered variable that was merged
		 * into element e via supernode detection or mass
		 * elimination of i when e became the pivot element.
		 * Place i in order just before e. */
		Next [i] = Next [e] ;
		Next [e]++ ;
	    }
	    else
	    {
		/* This is a dense unordered variable, with no parent.
		 * Place it last in the output order. */
		Next [i] = nel++ ;
	    }
	}
    }

    for (i = 0 ; i < n ; i++)
    {
	k = Next [i] ;
	Last [k] = i ;
    }
}


// === amd_order (rewritten with std::vector) ===
static int AMD_order(Int n, const Int Ap[], const Int Ai[], Int P[],
                     double Control[], double Info[])
{
    Int i, status;
    bool info = (Info != nullptr);
    if (info) {
        for (i = 0; i < AMD_INFO; i++) Info[i] = EMPTY;
        Info[AMD_N] = n;
        Info[AMD_STATUS] = AMD_OK;
    }
    if (Ai == nullptr || Ap == nullptr || P == nullptr || n < 0) {
        if (info) Info[AMD_STATUS] = AMD_INVALID;
        return AMD_INVALID;
    }
    if (n == 0) return AMD_OK;

    Int nz = Ap[n];
    if (info) Info[AMD_NZ] = nz;
    if (nz < 0) {
        if (info) Info[AMD_STATUS] = AMD_INVALID;
        return AMD_INVALID;
    }
    if (((size_t)n) >= (size_t)Int_MAX_VAL / sizeof(Int)
     || ((size_t)nz) >= (size_t)Int_MAX_VAL / sizeof(Int)) {
        if (info) Info[AMD_STATUS] = AMD_OUT_OF_MEMORY;
        return AMD_OUT_OF_MEMORY;
    }

    status = AMD_valid(n, n, Ap, Ai);
    if (status == AMD_INVALID) {
        if (info) Info[AMD_STATUS] = AMD_INVALID;
        return AMD_INVALID;
    }

    std::vector<Int> Len(n), Pinv(n);
    const Int *Cp, *Ci;
    std::vector<Int> Rp, Ri;

    if (status == AMD_OK_BUT_JUMBLED) {
        Rp.resize(n + 1);
        Ri.resize(std::max(nz, (Int)1));
        AMD_preprocess(n, Ap, Ai, Rp.data(), Ri.data(), Len.data(), Pinv.data());
        Cp = Rp.data();
        Ci = Ri.data();
    } else {
        Cp = Ap;
        Ci = Ai;
    }

    size_t nzaat = AMD_aat(n, Cp, Ci, Len.data(), P, Info);

    size_t slen = nzaat + nzaat / 5 + 7 * (size_t)n;
    if (slen < nzaat || slen >= SIZE_T_MAX / sizeof(Int)) {
        if (info) Info[AMD_STATUS] = AMD_OUT_OF_MEMORY;
        return AMD_OUT_OF_MEMORY;
    }
    if (info) Info[AMD_MEMORY] = (double)((2*n + slen) * sizeof(Int));

    std::vector<Int> S(slen);
    AMD_1(n, Cp, Ci, P, Pinv.data(), Len.data(), (Int)slen, S.data(),
          Control, Info);

    if (info) Info[AMD_STATUS] = status;
    return status;
}

} // anonymous namespace

// Public C-linkage entry point
extern "C" int amd_order(int32_t n, const int32_t Ap[], const int32_t Ai[],
                          int32_t P[], double Control[], double Info[])
{
    return AMD_order(n, Ap, Ai, P, Control, Info);
}
