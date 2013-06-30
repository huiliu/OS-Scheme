#include "avl.h"
#include <string.h>
#include <stdlib.h>
#include <stdio.h>

int main(int argc, char** argv)
{
    AVL_DATA *item;
    struct avl_node *node;
    int (*comparison_func)(const void*, const void*, void*);
    struct libavl_allocator *allocator;
    struct avl_table *tree;
    void *param;

    allocator = calloc(sizeof(struct libavl_allocator),1);
    comparison_func = calloc(sizeof(int), 1);
    item = malloc(sizeof(AVL_DATA));
    item->w = 1.2;
    item->value =  3.1415;
    node = malloc(sizeof(struct avl_node));
    node->avl_data = item;

    tree = avl_create(comparison_func, NULL, NULL);
    avl_insert(tree, node);
    param = avl_find(tree, node);
    if (param != NULL)
        fprintf(stdout, "I found it!\n");

    return 0;
}
