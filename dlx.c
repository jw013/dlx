/**
 * @author jw013
 * @date 2012
 * Implementation of Donald Knuth's
 * <a href="http://www-cs-faculty.stanford.edu/~uno/papers/dancing-color.ps.gz">
 * Dancing Links Algorithm</a> for solving Exact Cover (binary matrix
 * formulation, columns are the objects being covered by 1 elements in rows).
 *
 * All algorithms taken straight out of Knuth's DLX paper, translated fairly
 * literally into C.
 *
 * Summary of fundamental idea behind Knuth's DLX algorithm:
 * (1) Remove x from list:
 *	x->left->right = x->right;
 *	x->right->left = x->left;
 * (2) Restore x to its original position:
 *	x->left->right = x;
 *	x->right->left = x;
 */

#include "dlx.h"

/*
 * Private utility functions for manipulating node links
 * @{
 */

/** Remove node n from its left-right list */
static void
remove_lr(struct dlx_node *n)
{
	n->left->right = n->right;
	n->right->left = n->left;
}

/** Remove node n from its up-down list */
static void
remove_ud(struct dlx_node *n)
{
	n->up->down = n->down;
	n->down->up = n->up;
}

/** Restore node n to its left-right list */
static void
insert_lr(struct dlx_node *n)
{
	n->left->right = n->right->left = n;
}

/** Restore node n to its up-down list */
static void
insert_ud(struct dlx_node *n)
{
	n->up->down = n->down->up = n;
}

/** @return 1 if node has been removed from its up-down list, 0 otherwise */
static int
is_removed_ud(struct dlx_node *n)
{
	/*
	 * A node has been removed from its list if and only if both neighbors
	 * do not point to itself.  However, it is not possible for a node to
	 * be half in the list (unless the list is corrupted), so only one side
	 * is checked.
	 */
	return n->up->down != n;
}

/**
 * Insert node n into bottom of column c and update column node count.
 *
 * @param n	new node to insert.  If n is already part of the column, this
 * 		function will break the matrix horribly.
 */
static void
append_node_to_column(struct dlx_node *n, struct dlx_hnode *c)
{
	n->header = c;
	n->up	  = ((struct dlx_node *) c)->up;
	n->down	  = ((struct dlx_node *) c);
	insert_ud(n);
	c->node_count++;
}

/**
 * Cover column n.  The procedure is the same no matter which row is actually
 * used for covering.  First, the column header is removed from the header
 * list.  Then, the rest of the nodes in the column have their rows removed
 * from the matrix, by removing each row node not in column n from its column.
 *
 * @param n	header node for the column to cover
 */
static void
cover(struct dlx_node *h)
{
	struct dlx_node *i, *j;

	remove_lr(h);
	i = h;
	while ((i = i->down) != h) {	/* for each row, except header row */
		j = i;
		while ((j = j->right) != i) {	/* for each column node except i */
			remove_ud(j);
			j->header->node_count--;
		}
	}
}

/**
 * Reverse the procedure in cover(h), for backtracking.
 *
 * Must be called in exact reverse order as cover() to ensure matrix is
 * correctly restored to original state.
 *
 * @param n	header node for the column to uncover
 */
static void
uncover(struct dlx_node *h)
{
	struct dlx_node *i, *j;

	/* all loops MUST traverse in OPPOSITE order from cover() */
	i = h;
	while ((i = i->up) != h) {	/* for each row except header row */
		j = i;
		while ((j = j->left) != i) {	/* for each column except i */
			j->header->node_count++;
			insert_ud(j);
		}
	}
	insert_lr(h);
}

/**
 * Cover all columns of nodes in the same row as i, except the column of i
 * itself.
 *
 * The inverse function is uncover_other_columns, which ensures that
 * matrix links are restored in the correct order.
 */
static void
cover_other_columns(struct dlx_node *i)
{
	struct dlx_node *j = i;
	while ((j = j->right) != i)	/* for each column except i */
		cover((struct dlx_node *) j->header);
}

/**
 * Reverse the procedure in cover_other_columns(i), for backtracking.
 *
 * Must be called in exact reverse order as cover_other_columns() to ensure
 * matrix is correctly restored to original state.
 */
static void
uncover_other_columns(struct dlx_node *i)
{
	struct dlx_node *j = i;
	while ((j = j->left) != i)	/* for each column except i */
		uncover((struct dlx_node *) j->header);
}

/**
 * @return 	a column header with the smallest node count, or NULL if column
 * 		list is empty
 */
static struct dlx_hnode *
hnode_min_count(struct dlx_hnode *root)
{
	struct dlx_node *h = (struct dlx_node *) root;
	struct dlx_node *min = NULL;	/* return value */

	while ((h = h->right) != (struct dlx_node *) root) {
		if (min == NULL || ((struct dlx_hnode *) h)->node_count <
					((struct dlx_hnode *) min)->node_count)
			min = h;
	}
	return (struct dlx_hnode *) min;
}

/* @} */

/*
 * Functions for specifying exact cover problems by allowing user to pre-select
 * certain rows as part of the solution.
 * @{
 */

/**
 * Cover all columns that row r covers.
 *
 * This can be useful if you want to force a certain row to be included in the
 * solution (hence the name).
 *
 * @return 	0 on success, -1 if r has already been removed from the matrix and
 * 		cannot be selected.
 */
int
dlx_force_row(struct dlx_node *r)
{
	if (is_removed_ud(r))
		return -1;

	cover((struct dlx_node *) r->header);
	cover_other_columns(r);
	return 0;
}

/**
 * Inverse operation of dlx_force_row.
 *
 * Must be called in exact reverse order as calls to dlx_force_row for matrix
 * links to be restored properly.
 *
 * @return 0 on success, -1 if r is still in the matrix.
 */
int
dlx_unselect_row(struct dlx_node *r)
{
	if (!is_removed_ud(r))
		return -1;

	uncover_other_columns(r);
	uncover((struct dlx_node *) r->header);
	return 0;
}

/* @} */

/*
 * Functions to aid in initializing node links when constructing a DLX matrix.
 * @{
 */

/**
 * Set up the header nodes for a DLX matrix by connecting the array of n
 * pre-allocated column headers and the root node into a circularly linked
 * list (oriented left-right, as a row).
 *
 * The id members of the nodes are left un-touched.
 *
 * @param root		pointer to root node
 * @param headers	pre-allocated array of n header nodes
 * @param n		number of column headers, not including root node
 */
void
dlx_make_header_row(struct dlx_hnode *root, struct dlx_hnode *headers, size_t n)
{
	struct dlx_node *hi;
	size_t i;

	/* special case: zero columns needs to be handled separately */
	if (n < 1) {
		/* set up the root node */
		hi = (struct dlx_node *) root;
		hi->left	= hi;
		hi->right	= hi;
		hi->up		= NULL;
		hi->down	= NULL;
		hi->header	= NULL;
		root->node_count = 1;
		return;
	}

	/* set up the root node */
	hi = (struct dlx_node *) root;
	hi->left	= (struct dlx_node *) (headers + n - 1);
	hi->right	= (struct dlx_node *) headers;
	hi->up		= NULL;
	hi->down	= NULL;
	hi->header	= NULL;
	root->node_count = 1;

	/*
	 * set up the actual column headers:
	 * (*) left and right links point left and right
	 * (*) up and down links point to self
	 * (*) header points to self
	 * (*) initial node count is 0
	 * (*) id is not touched
	 */

	/* special case: a single header needs to be handled separately */
	if (n == 1) {
		hi = (struct dlx_node *) headers;
		hi->left	= (struct dlx_node *) root;
		hi->right	= (struct dlx_node *) root;
		hi->up		= hi;
		hi->down	= hi;
		hi->header	= (struct dlx_hnode *) hi;
		((struct dlx_hnode *) hi)->node_count = 1;
		return;
	}

	/*
	 * this code assumes at least two headers: the first header must be
	 * followed by another header and the last header must precededed by a
	 * header
	 */
	/* first column header */
	hi = (struct dlx_node *) headers;
	hi->left	= (struct dlx_node *) root;
	hi->right	= (struct dlx_node *) (headers + 1);
	hi->up		= hi;
	hi->down	= hi;
	hi->header	= (struct dlx_hnode *) hi;
	((struct dlx_hnode *) hi)->node_count = 1;

	/* from 2nd to 2nd to last column header */
	for (i = 1; i < n - 1; i++) {
		hi = (struct dlx_node *) (headers + i);
		hi->left  = (struct dlx_node *) (((struct dlx_hnode *) hi) - 1);
		hi->right = (struct dlx_node *) (((struct dlx_hnode *) hi) + 1);
		hi->up		= hi;
		hi->down	= hi;
		hi->header	= (struct dlx_hnode *) hi;
		((struct dlx_hnode *) hi)->node_count = 1;
	}

	/* last column header */
	hi = (struct dlx_node *) (headers + n - 1);
	hi->left  = (struct dlx_node *) (((struct dlx_hnode *) hi) - 1);
	hi->right = (struct dlx_node *) root;
	hi->up		= hi;
	hi->down	= hi;
	hi->header	= (struct dlx_hnode *) hi;
	((struct dlx_hnode *) hi)->node_count = 1;
}

/**
 * Set up the nodes for a DLX matrix by connecting the array of n pre-allocated
 * nodes into a circularly linked list (oriented left-right, as a row).
 *
 * @param nodes		pre-allocated array of n nodes
 * @param n		number of nodes in the row nodes[]
 */
void
dlx_make_row(struct dlx_node *nodes, void *row_id, size_t n)
{
	size_t i;
	struct dlx_node *ni = nodes;

	if (n < 1)
		return;

	/* special case: a single node row needs to be handled separately */
	if (n == 1) {
		ni->left	= ni;
		ni->right	= ni;
		ni->row_id	= row_id;
		return;
	}

	/*
	 * this code assumes at least two headers: the first header must be
	 * followed by another header and the last header must precededed by a
	 * header
	 */
	/* first node */
	ni->left	= ni + n - 1;
	ni->right	= ni + 1;
	ni->row_id	= row_id;

	/* from 2nd node to 2nd from last node */
	for (ni++, i = 1; i < n - 1; i++, ni++) {
		ni->left	= ni - 1;
		ni->right	= ni + 1;
		ni->row_id	= row_id;
	}

	/* last node */
	ni->left	= ni - 1;
	ni->right	= nodes;
	ni->row_id	= row_id;
}

/**
 * Add pre-initialized row to the DLX matrix by appending each row node in the
 * array to the corresponding column.
 *
 * @param nodes		row nodes to append
 * @param headers	parallel array to nodes of pointers to headers of
 * 			columns to which the row nodes will be appended.  Must
 * 			have the same number of elements as nodes and the
 * 			headers must occur in the same order.
 * @param n		number of nodes in row
 */
void
dlx_add_row(struct dlx_node *nodes, struct dlx_hnode **headers, size_t n)
{
	size_t i;
	for (i = 0; i < n; i++)
		append_node_to_column(&nodes[i], headers[i]);
}

/* @} */

/*
 * Variations on the core DLX algorithm for solving Exact Cover by D. Knuth.
 * @{
 */

/**
 * Exact cover DLX algorithm by D. Knuth, with some modification to allow for
 * skipping over a specified number of solutions.
 *
 * @param solution
 * 		Make sure to allocate enough space to contain the largest
 * 		possible cover solution to prevent buffer overflow.
 * @param root	pointer to root node of a valid DLX matrix structure.  The
 * 		matrix is modified by the function, but is restored to its
 * 		original state before the function returns.
 * @param k	must be 0; value is used internally.
 * @param pnsol	pointer to a value specifying the number of solutions to look
 * 		for.  The pointed to value will be decremented by the number of
 * 		solutions found, to a minimum of 0, and must be positive
 * 		initally.  The behavior of the function when *pnsol is zero is
 * 		undefined and most likely undesirable.
 * @return 	size of *pnsol'th solution, or 0 if not that many solutions
 * 		exist.  Note that the 0 case is ambiguous for an empty matrix,
 * 		which has a solution of size 0, but handling that situation is
 * 		left to the caller.
 */
size_t
dlx_exact_cover(struct dlx_srow *solution, struct dlx_hnode *root,
		size_t k, size_t *pnsol)
{
	/*
	 * Base cases: Recursion ends when either
	 *   * matrix is empty (success: entire matrix has been covered)
	 *   * empty column is found (failure: uncover-able column)
	 * Recursive/branching step:
	 *   * select column with fewest candidate rows
	 *   * select a row within that column for the solution
	 *   * recurse
	 */
	size_t n = 0;		/* return value, default 0 := no solution */
	struct dlx_node *i;	/* iterator pointer */
	struct dlx_node *col;	/* column to cover in this iteration */

	/* root->right == root means empty matrix */
	if (((struct dlx_node *) root)->right == (struct dlx_node *) root) {
		(*pnsol)--;
		return k;
	}

	col = (struct dlx_node *) hnode_min_count(root);
	cover(col);

	solution[k].id		= ((struct dlx_hnode *) col)->id;
	solution[k].n_choices	= ((struct dlx_hnode *) col)->node_count;

	/* try selecting each row in the column one at a time and recurse */
	i = col;
	while ((i = i->down) != col) {	/* for each row except header row */
		cover_other_columns(i);
		n = dlx_exact_cover(solution, root, k + 1, pnsol);
		uncover_other_columns(i);
		if (n > 0)	/* solution found */
			solution[k].row_node = i;
		if (*pnsol == 0)
			break;
	}

	/* restore node links and return */
	uncover(col);
	return n;
}

/* @} */

/* misc @{ */

/**
 * @param node
 * @return Pointer to row id of the row, or NULL of node is null.
 */
void *
dlx_row_id(struct dlx_node *node)
{
	return node ? node->row_id : NULL;
}

/* @} */
