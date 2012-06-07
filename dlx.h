/**
 * @author jw013
 * @date 2012
 * @brief Implementation of Donald Knuth's
 * <a href="http://en.wikipedia.org/wiki/Dancing_Links">"Dancing Links"</a>
 * algorithm.
 */

#ifndef DLX_H
#define DLX_H

#include <stddef.h>

/* ===============
 * Data Structures
 */

struct dlx_hnode;

/*
 * DLX node.
 * OPAQUE structure: do NOT access members.
 *
 * Each node can simultaneously be a member of two linked lists, one oriented
 * vertical (up-down) and the other oriented horizontally (left-right).  Nodes
 * also contain a pointer to the list header to allow access to the header in
 * O(1) time.
 */
struct dlx_node {
	void *row_id;
	struct dlx_node *left;
	struct dlx_node *right;
	struct dlx_node *up;
	struct dlx_node *down;
	struct dlx_hnode *header;
};

/*
 * List header node.
 * OPAQUE structure: do NOT access members.
 *
 * Special dlx_node that also tracks the number of nodes in its list
 */
struct dlx_hnode {
	struct dlx_node base;	/**< inherit fields from struct dlx_node */
	size_t node_count;	/**< number of nodes in the list */
	const void *id;		/**< pointer to unique list id */
};

/**
 * Convenient handle for holding pointers to the various data structures
 * necessary to store a DLX matrix.  While a simple root node pointer suffices
 * to solve exact cover on one, a few more pointers are needed for allocating
 * and tracking memory.  The row_off member behaves just like the typical
 * row_ptr in a compressed sparse row representation of a sparse matrix.
 *
 * All members are public.
 */
struct dlx_matrix {
	struct dlx_hnode root;
	struct dlx_hnode *headers;
	struct dlx_node *nodes;
	size_t *row_off;	/* row_off[i] ... row_off[i+1] = row i's nodes */
	size_t n_col;
	size_t n_row;		/* number of rows = length of row_off - 1 */
};

/**
 * Solution row.
 * Additional information about a selected solution row (you can think of the
 * s in srow as either solution or selected).
 *
 * If a row is included in the solution, it may be useful to know
 * (*) which column the row was selected to cover, a.k.a. its primary column
 *  	(though the row may also incidentally cover other columns, it was the
 *  	"best" choice for covering its primary column);
 * (*) the number of other possible row choices for that column.
 */
struct dlx_srow {
	struct dlx_node *row_node;	/**< row included in the solution */
	const void *cid;  		/**< id of "primary" column */
	size_t n_choices;		/**< number of possible choices of row */
};


/* ===
 * functions for INITIALIZING matrices by initializing node links
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
dlx_make_header_row(struct dlx_hnode *root, struct dlx_hnode *headers, size_t n);

/**
 * Set up the nodes for a DLX matrix by connecting the array of n pre-allocated
 * nodes into a circularly linked list (oriented left-right, as a row).
 *
 * @param nodes		pre-allocated array of n nodes
 * @param n		number of nodes in the row nodes[]
 */
void dlx_make_row(struct dlx_node *nodes, void *row_id, size_t n);

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
void dlx_add_row(struct dlx_node *nodes, struct dlx_hnode **columns, size_t n);

/* ===
 * functions for SPECIFYING exact cover problems by pre-selecting rows
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
int dlx_force_row(struct dlx_node *r);

/**
 * Inverse operation of dlx_force_row.
 *
 * Must be called in exact reverse order as calls to dlx_force_row for matrix
 * links to be restored properly.
 *
 * @return 0 on success, -1 if r is still in the matrix.
 */
int dlx_unselect_row(struct dlx_node *r);

/* ===
 * functions for SOLVING exact cover
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
		size_t k, size_t *pnsol);

/* ===
 * misc
 */

/**
 * @param node
 * @return Pointer to row id of the row, or NULL of node is null.
 */
void *dlx_row_id(struct dlx_node *node);

#endif
