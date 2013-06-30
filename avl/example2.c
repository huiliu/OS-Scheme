#include "redblack.h"
#include <stdlib.h>
#include <stdio.h>

/*
 * This script demonstrates the worst case scenario - entering
 * the data already in sequence. This program enters 200 numbers
 * in reverse sequence (200 to 0) and then prints them out in
 * the usual order. This would kill a regular tree algorithm.
 *
 * This is the same as example1, except that the output is done
 * using rbwalk.
 */

void *xmalloc(unsigned n)
{
	void *p;
	p = malloc(n);
	if(p) return p;
	fprintf(stderr, "insufficient memory\n");
	exit(1);
}

int compare(const void *pa, const void *pb, const void *config)
{
    ITEM *a, *b;
    a = (ITEM *)pa; b = (ITEM *)pb;
	if(a->w < b->w) return -1;
	if(a->w > b->w) return 1;

/*
	if(*(int *)pa < *(int *)pb) return -1;
	if(*(int *)pa > *(int *)pb) return 1;
*/
	return 0;
}

int minleaf=-1;
int maxleaf=-1;

void walkact(const rbdata_t *p, const VISIT which, const int depth, void *msg)
{
	if (which == postorder || which == leaf) {
		printf("%s: %4d (depth=%2d)\n", (char *) msg, *(int *)p, depth);
	}

	if (which == leaf) {
		if (minleaf==-1 || depth < minleaf)
			minleaf=depth;
		if (maxleaf==-1 || depth > maxleaf)
			maxleaf=depth;
	}
}

int main(int argc, char** argv)
{
    ITEM *ptr, *result;
	struct rbtree *rb;
    ITEM c;
    FILE *f;
    char *dbName = "f_inc.0.u.test";

    c.w = atof(argv[1]);
    c.value = -1;

	if ((rb=rbinit(compare, NULL))==NULL)
	{
		fprintf(stderr, "insufficient memory\n");
		exit(1);
	}

    f = fopen(dbName, "r");
    ptr = (ITEM *)xmalloc(sizeof(ITEM));
    while (fscanf(f, "%lf", &ptr->w) != EOF) {
        fscanf(f, "%lf", &ptr->value);
        rbsearch(ptr, rb);
		ptr = (ITEM *)xmalloc(sizeof(ITEM));
    }

    result = (ITEM *)rbsearch(&c, rb);
    fprintf(stdout, "%.8lf\t%.8lf\n", result->w, result->value);

/*
	for (i = 50000; i > 0; i--)
	{
		ptr = (ITEM *)xmalloc(sizeof(ITEM));
		ptr->w = i;
		ptr->value = (double)i/7;
		val = rbsearch((void *)ptr, rb);
		if(val == NULL)
		{
			fprintf(stderr, "insufficient memory\n");
			exit(1);
		}
	}

	rbwalk(rb, walkact, "No");

	printf("Minimum branch length: %d\n", minleaf);
	printf("Maximum branch length: %d\n", maxleaf);
*/
	rbdestroy(rb);
	
	return 0;
}
