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
	const void *id;	  		/**< id of "primary" column */
	size_t n_choices;		/**< number of possible choices of row */
};

/* functions for INITIALIZING matrices by setting up node links */
void
dlx_make_header_row(struct dlx_hnode *root, struct dlx_hnode *headers, size_t n);
void dlx_make_row(struct dlx_node *nodes, void *row_id, size_t n);
void dlx_add_row(struct dlx_node *nodes, struct dlx_hnode **columns, size_t n);

/* functions for SPECIFYING exact cover problems by pre-selecting rows */
int dlx_force_row(struct dlx_node *r);
int dlx_unselect_row(struct dlx_node *r);

/* functions for SOLVING exact cover */
size_t
dlx_exact_cover(struct dlx_srow *solution, struct dlx_hnode *root,
		size_t k, size_t *pnsol);

/* misc */
void *dlx_row_id(struct dlx_node *node);

#endif
