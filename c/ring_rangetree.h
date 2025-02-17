#pragma once

#include <math.h>
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>

#include <motionlib/utils.h>

typedef struct point2 {
    union {
        motion_dtype x[2];
        struct {
            motion_dtype r;
            motion_dtype theta;
        };
    };
} point2;

typedef struct treenode {
    struct treenode* left;
    struct treenode* right;
    point2 val;
    size_t num_subpoints;
    point2 data_ptr[0];
} treenode;


static inline treenode* make_treenode(point2* val, size_t num_subpoints, point2* data) {
    treenode* t = (treenode*) malloc(sizeof(treenode) + num_subpoints * sizeof(point2));
    t->left = NULL;
    t->right = NULL;
    memcpy(&t->val, val, sizeof(point2));
    t->num_subpoints = num_subpoints;
    memcpy(t->data_ptr, data, num_subpoints * sizeof(point2));
    return t;
}

static inline void treenode_free_recursive(treenode* node) {
    if (node) {
        treenode_free_recursive(node->left);
        treenode_free_recursive(node->right);
        free(node);
    }
}


typedef struct cylindertree {
    size_t num_points;
    motion_dtype max;
    treenode* root;

    // "List" for making reports of points.
    size_t report_size;
    point2* report_scratch;
} cylindertree;


static inline void inplace_make_cylindertree(cylindertree* t, motion_dtype max) {
    t->max = max;
    t->root = NULL;
    t->num_points = 0;
    t->report_size = 0;
    t->report_scratch = NULL;
}

static inline void cylindertree_free(cylindertree* t) {
    treenode_free_recursive(t->root);
    free(t->report_scratch);
}

static inline int cmp_r(const point2* arg1, const point2* arg2) {
    return signum(arg1->r - arg2->r);
}

// NOTE: CLOBBERS point_data!
static inline treenode* construct_recurse(size_t n_points, point2* point_data) {
    if (n_points == 0) { return NULL; }
    if (n_points == 1) { return make_treenode(point_data, 1, point_data); }

    // "Random" pivot. lol. In our usecase theta is usually sorted so this is good but mmmm
    // Too lazy to implement out-of-place MoM+quickselect median find.
    size_t mid_idx = n_points / 2;
    motion_dtype mid_theta = point_data[mid_idx].theta;

    // Makes a copy of point_data. We're about to clobber it :)
    treenode* root = make_treenode(&point_data[mid_idx], n_points, point_data);

    // Overkill size. It's OK because this is not gonna be stored.
    point2* left_points = point_data;
    point2* right_points = (point2*) malloc(sizeof(point2) * n_points);

    size_t left_size = 0;
    size_t right_size = 0;
    for (size_t i = 0; i < n_points; ++i) {
        if (i == mid_idx) { continue; }
        if (point_data[i].theta < mid_theta) {
            left_points[left_size] = point_data[i];
            ++left_size;
        }
        else {
            right_points[right_size] = point_data[i];
            ++right_size;
        }
    }
    root->left = construct_recurse(left_size, left_points);
    root->right = construct_recurse(right_size, right_points);
    // left_points is held by whoever is above us in the tree.
    free(right_points);
    return root;
}

// NOTE: CLOBBERS point_data!
static inline void cylindertree_construct(cylindertree* t, size_t n_points, point2* point_data) {
    treenode_free_recursive(t->root);
    t->num_points = n_points;
    qsort(point_data, n_points, sizeof(point2), (int(*)(const void*, const void*)) &cmp_r);
    t->root = construct_recurse(n_points, point_data);
    t->report_size = 0;
    t->report_scratch = (point2*) realloc(t->report_scratch, n_points * sizeof(point2));
}

#define report_item(t, item) { \
    (t)->report_scratch[(t)->report_size] = (item); \
    (t)->report_size += 1; \
}

static inline void _cylindertree_traverse(cylindertree* t, treenode* root) {
    if (root == NULL) { return; }
    _cylindertree_traverse(t, root->left);
    report_item(t, root->val);
    _cylindertree_traverse(t, root->right);
}

static inline void cylindertree_traverse(cylindertree* t) {
    t->report_size = 0;
    _cylindertree_traverse(t, t->root);
}

// Following function copied from C++ stl.
static inline size_t first_greater_index(size_t n_points, point2* arr, motion_dtype x) {
    // NOTE: inclusive index.
    size_t search_window = n_points;
    size_t cur_idx = 0;
    while (search_window > 0) {
        size_t half_dist = search_window >> 1;
        size_t mid_idx = cur_idx + half_dist;
        if (arr[mid_idx].r > x) {
            /* We jumped too far.*/
            search_window = half_dist;
        }
        else {
            /* Jump to the item tested.
             * But we know the item tested doesn't satisfy the condition,
             * so we can inch forward by 1 more. */
            cur_idx = mid_idx + 1;
            search_window -= half_dist + 1;
        }
    }
    return cur_idx;
}
static inline size_t last_smaller_index(size_t n_points, point2* arr, motion_dtype x) {
    // NOTE: exclusive index.
    size_t search_window = n_points;
    size_t cur_idx = 0;
    while (search_window > 0) {
        size_t half_dist = search_window >> 1;
        size_t mid_idx = cur_idx + half_dist;
        if (arr[mid_idx].r < x) {
            /* Jump to the item tested.
             * But we know the item tested doesn't satisfy the condition,
             * so we can inch forward by 1 more. */
            cur_idx = mid_idx + 1;
            search_window -= half_dist + 1;
        }
        else {
            /* We jumped too far.*/
            search_window = half_dist;
        }
    }
    return cur_idx;
}

static inline void report_right(cylindertree* t, treenode* root, motion_dtype theta, motion_dtype r0, motion_dtype r1) {
    if (root == NULL) { return; }
    if (root->val.theta > theta) {
        report_right(t, root->left, theta, r0, r1);
        if (root->val.r > r0 && root->val.r < r1) {
            report_item(t, root->val);
        }
        if (root->right != NULL) {
            size_t left_idx = first_greater_index(root->right->num_subpoints, root->right->data_ptr, r0);
            size_t right_idx = last_smaller_index(root->right->num_subpoints, root->right->data_ptr, r1);
            if (left_idx < right_idx) {
                size_t n_copy = right_idx - left_idx;
                memcpy(&t->report_scratch[t->report_size], &root->right->data_ptr[left_idx], n_copy * sizeof(point2));
                t->report_size += n_copy;
            }
        }
    }
    else {
        report_right(t, root->right, theta, r0, r1);
    }
}

static inline void report_left(cylindertree* t, treenode* root, motion_dtype theta, motion_dtype r0, motion_dtype r1) {
    if (root == NULL) { return; }
    if (root->val.theta < theta) {
        if (root->left != NULL) {
            size_t left_idx = first_greater_index(root->left->num_subpoints, root->left->data_ptr, r0);
            size_t right_idx = last_smaller_index(root->left->num_subpoints, root->left->data_ptr, r1);
            if (left_idx < right_idx) {
                size_t n_copy = right_idx - left_idx;
                memcpy(&t->report_scratch[t->report_size], &root->left->data_ptr[left_idx], n_copy * sizeof(point2));
                t->report_size += n_copy;
            }
        }
        if (root->val.r > r0 && root->val.r < r1) {
            report_item(t, root->val);
        }
        report_left(t, root->right, theta, r0, r1);
    }
    else {
        report_left(t, root->left, theta, r0, r1);
    }
}

static inline void get_range(cylindertree* t, motion_dtype r_start, motion_dtype r_end, motion_dtype theta_start, motion_dtype theta_end) {
    // Return value is stored in t itself, under `report_size` and `report_scratch`.
    t->report_size = 0;
    
    // Always sweeping from start -> end.
    if (theta_start > theta_end) {
        // Actually easy case. We just report everything to the right of x_start and left of x_end.
        report_right(t, t->root, theta_start, r_start, r_end);
        report_left(t, t->root, theta_end, r_start, r_end);
        return;
    }

    treenode* node = t->root;
    while (node) {
        if (node->val.theta >= theta_start && node->val.theta <= theta_end) {
            // Bracketed the root. Start reporting!
            report_right(t, node->left, theta_start, r_start, r_end);
            if (node->val.r > r_start && node->val.r < r_end) {
                report_item(t, node->val);
            }
            report_left(t, node->right, theta_end, r_start, r_end);
            return;
        }
        if (node->val.theta < theta_start) {
            node = node->right;
        }
        else {
            node = node->left;
        }
    }
}

static inline bool cc_right(treenode* root, motion_dtype theta, motion_dtype r0, motion_dtype r1) {
    if (root == NULL) { return false; }
    if (root->val.theta > theta) {
        if (root->val.r > r0 && root->val.r < r1) { return true; }
        if (root->right != NULL) {
            size_t left_idx = first_greater_index(root->right->num_subpoints, root->right->data_ptr, r0);
            size_t right_idx = last_smaller_index(root->right->num_subpoints, root->right->data_ptr, r1);
            if (left_idx < right_idx) { return true; }
        }
        return cc_right(root->left, theta, r0, r1);
    }
    else {
        return cc_right(root->right, theta, r0, r1);
    }
}

static inline bool cc_left(treenode* root, motion_dtype theta, motion_dtype r0, motion_dtype r1) {
    if (root == NULL) { return false; }
    if (root->val.theta < theta) {
        if (root->val.r > r0 && root->val.r < r1) { return true; }
        if (root->left != NULL) {
            size_t left_idx = first_greater_index(root->left->num_subpoints, root->left->data_ptr, r0);
            size_t right_idx = last_smaller_index(root->left->num_subpoints, root->left->data_ptr, r1);
            if (left_idx < right_idx) { return true; }
        }
        return cc_left(root->right, theta, r0, r1);
    }
    else {
        return cc_left(root->left, theta, r0, r1);
    }
}

static inline bool collision_check(cylindertree* t, motion_dtype r_start, motion_dtype r_end, motion_dtype theta_start, motion_dtype theta_end) {
    // Always sweeping from start -> end.
    if (theta_start > theta_end) {
        return (
            cc_right(t->root, theta_start, r_start, r_end)
            || cc_left(t->root, theta_end, r_start, r_end)
        );
    }

    treenode* node = t->root;
    while (node) {
        if (node->val.theta >= theta_start && node->val.theta <= theta_end) {
            // Bracketed the root. Start reporting!
            if (node->val.r > r_start && node->val.r < r_end) {
                return true;
            }
            return (
                cc_right(node->left, theta_start, r_start, r_end)
                || cc_left(node->right, theta_end, r_start, r_end)
            );
        }
        if (node->val.theta < theta_start) {
            node = node->right;
        }
        else {
            node = node->left;
        }
    }
    return false;
}
