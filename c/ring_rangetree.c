#include <stdio.h>

#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include <structmember.h>

#include "ring_rangetree.h"

typedef struct {
    PyObject_HEAD
    /* Type-specific fields go here. */
    cylindertree tree;

} CylinderTreeObject;

static void
CylinderTree_dealloc(CylinderTreeObject* self) {
    cylindertree_free(&self->tree);
}

static PyObject*
CylinderTree_new(PyTypeObject* type, PyObject* args, PyObject* kwds) {
    static char* kwlist[] = {
        "theta_max", NULL
    };

    double theta_max;
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "d", kwlist, &theta_max)){
        return NULL;
    }
    if (theta_max <= 0) { PyErr_SetString(PyExc_TypeError, "theta_max must be greater than zero"); return NULL; }

    CylinderTreeObject* self = (CylinderTreeObject*) type->tp_alloc(type, 0);
    inplace_make_cylindertree(&self->tree, theta_max);
    return (PyObject*) self;
}

static int
CylinderTree_init(CylinderTreeObject* self, PyObject* args, PyObject* kwds) {
    return 0;
}

static PyObject*
CylinderTree_construct(CylinderTreeObject* self, PyObject* data) {
    int res = PyObject_CheckBuffer(data);
    if (res == 0) { PyErr_SetString(PyExc_TypeError, "Expect an input array of shape (N, 2)"); return NULL; }

    Py_buffer* buffer = malloc(sizeof(Py_buffer));
    res = PyObject_GetBuffer(data, buffer, PyBUF_RECORDS_RO);
    if (res == -1) { PyErr_SetString(PyExc_TypeError, "Could not create simple buffer from input data"); return NULL; }
    if (buffer->ndim != 2) { PyErr_SetString(PyExc_TypeError, "Expect an input array of shape (N, 2)"); return NULL; }
    if (buffer->shape[1] != 2) { PyErr_SetString(PyExc_TypeError, "Expect an input array of shape (N, 2)"); return NULL; }
    

    size_t element_size;
    if (strlen(buffer->format) > 1) {
        PyErr_SetString(PyExc_TypeError, "Unexpected buffer type (expected float or double)");
        PyBuffer_Release(buffer);
        return NULL;
    }
    if (buffer->format[0] == 'd') {
        element_size = 8;
    }
    else if (buffer->format[0] == 'f') {
        element_size = 4;
    }
    else {
        PyErr_SetString(PyExc_TypeError, "Unexpected buffer type (expected float or double)");
        PyBuffer_Release(buffer);
        return NULL;
    }

    size_t num_elements = buffer->shape[0] * buffer->shape[1];
    void* buf_contig = malloc(num_elements * buffer->itemsize);
    PyBuffer_ToContiguous(buf_contig, buffer, buffer->len, 'C');

    motion_dtype* point_data;
    if (element_size == sizeof(motion_dtype)) {
        point_data = buf_contig;
    }
    else {
        point_data = malloc(num_elements * sizeof(motion_dtype));
        if (element_size == 4) {
            float* buf_ptr = (float*) buf_contig;
            for (size_t i = 0; i < num_elements; ++i) {
                point_data[i] = buf_ptr[i];
            }
        }
        else if (element_size == 8) {
            double* buf_ptr = (double*) buf_contig;
            for (size_t i = 0; i < num_elements; ++i) {
                point_data[i] = buf_ptr[i];
            }
        }
        free(buf_contig);
    }

    // NOTE: point_data is clobbered after this call :)
    cylindertree_construct(&self->tree, buffer->shape[0], (point2*) point_data);

    free(point_data);
    PyBuffer_Release(buffer);
    free(buffer);
    Py_RETURN_NONE;
}

static PyObject* points_to_py(size_t n, point2* points) {
    PyObject* ret = PyList_New(n);
    for (size_t i = 0; i < n; ++i) {
        PyObject* point = vector_to_list(points[i].x, 2);
        PyList_SetItem(ret, i, point);
    }
    return ret;
}

int parse_range_args(CylinderTreeObject* self, PyObject* const* args, Py_ssize_t nargs,
        motion_dtype* r_bounds, motion_dtype* theta_bounds) {
#ifdef MOTION_DEBUG
    if (nargs != 2) {
        PyErr_SetString(PyExc_TypeError, "Invalid arguments (expected: (r_bounds, theta_bounds))");
        return 1;
    }
#endif
    PyObject* _r_bounds = args[0];
    Py_ssize_t n;
#ifdef MOTION_DEBUG
    n = PyObject_Length(_r_bounds);
    if (n < 0) {
        PyErr_SetString(PyExc_TypeError, "r_bounds has no length");
        return 1;
    }
    if (n != 2) {    // No Forgiveness
        PyErr_SetString(PyExc_TypeError, "r_bounds should be length 2 (min, max)");
        return 1;
    }
#endif

    PyObject* _theta_bounds = args[1];
#ifdef MOTION_DEBUG
    n = PyObject_Length(_theta_bounds);
    if (n < 0) {
        PyErr_SetString(PyExc_TypeError, "theta_bounds has no length");
        return 1;
    }
    if (n != 2) {    // No Forgiveness
        PyErr_SetString(PyExc_TypeError, "theta_bounds should be length 2 (min, max)");
        return 1;
    }
#endif
    if (list_to_vector_n(_r_bounds, r_bounds, 2) == NULL) {
        return 1;
    }
    if (list_to_vector_n(_theta_bounds, theta_bounds, 2) == NULL) {
        return 1;
    }

#ifdef MOTION_DEBUG
    if (r_bounds[0] > r_bounds[1]) {
        PyErr_SetString(PyExc_ValueError, "r_bounds[1] should be greater than r_bounds[0]");
        return 1;
    }
    if (theta_bounds[0] > self->tree.max) {
        PyErr_SetString(PyExc_ValueError, "theta_bounds[0] exceeds maximum range of ring coordinate");
        return 1;
    }
    if (theta_bounds[1] > self->tree.max) {
        PyErr_SetString(PyExc_ValueError, "theta_bounds[1] exceeds maximum range of ring coordinate");
        return 1;
    }
#endif
    return 0;
}

PyObject* CylinderTree_find_points(PyObject* _self, PyObject* const* args, Py_ssize_t nargs) {
    CylinderTreeObject* self = (CylinderTreeObject*) _self;
    motion_dtype r_bounds[2];
    motion_dtype theta_bounds[2];
    if (parse_range_args(self, args, nargs, r_bounds, theta_bounds)) {
        return NULL;
    }
    
    get_range(&self->tree, r_bounds[0], r_bounds[1], theta_bounds[0], theta_bounds[1]);

    return points_to_py(self->tree.report_size, self->tree.report_scratch);
}

PyObject* CylinderTree_contains_point(PyObject* _self, PyObject* const* args, Py_ssize_t nargs) {
    CylinderTreeObject* self = (CylinderTreeObject*) _self;
    motion_dtype r_bounds[2];
    motion_dtype theta_bounds[2];
    if (parse_range_args(self, args, nargs, r_bounds, theta_bounds)) {
        return NULL;
    }
    bool res = collision_check(&self->tree, r_bounds[0], r_bounds[1], theta_bounds[0], theta_bounds[1]);
    return PyBool_FromLong(res);
}

static PyObject*
CylinderTree_traverse(CylinderTreeObject* self, PyObject* _args) {
    cylindertree_traverse(&self->tree);
    return points_to_py(self->tree.report_size, self->tree.report_scratch);
}

static PyMethodDef CylinderTree_methods[] = {
    {"construct", (PyCFunction) CylinderTree_construct, METH_O,
            PyDoc_STR("Create a rangetree from the given points.")},
    {"get_range", (PyCFunction) CylinderTree_find_points, METH_FASTCALL,
            PyDoc_STR("Return all points contained in the bounding box.")},
    {"contains_point", (PyCFunction) CylinderTree_contains_point, METH_FASTCALL,
            PyDoc_STR("Check if the given bounding box contains any points.")},
    {"traverse", (PyCFunction) CylinderTree_traverse, METH_NOARGS,
            PyDoc_STR("Get all points in the tree in tree order.")},
    {NULL}  /* Sentinel */
};

static PyObject*
CylinderTree_get_num_points(CylinderTreeObject* self, void* _closure) {
    return PyFloat_FromDouble(self->tree.num_points);
}

static PyObject*
CylinderTree_get_all_points(CylinderTreeObject* self, void* _closure) {
    return points_to_py(self->tree.num_points, self->tree.root->data_ptr);
}

static PyMemberDef CylinderTree_members[] = {
    {NULL}  /* Sentinel */
};

static PyGetSetDef CylinderTree_getsetters[] = {
    {"num_points", (getter) CylinderTree_get_num_points, (setter) NULL,
     "Number of points total stored in this rangetree.", NULL},
    {"points", (getter) CylinderTree_get_all_points, (setter) NULL,
     "All points stored in this rangetree.", NULL},
    {NULL}  /* Sentinel */
};

static PyTypeObject CylinderTreeType = {
    PyVarObject_HEAD_INIT(NULL, 0)
    .tp_name = "ring_rangetree.CylinderTree",
    .tp_basicsize = sizeof(CylinderTreeObject),
    .tp_itemsize = 0,
    .tp_dealloc = (destructor) CylinderTree_dealloc,
    .tp_flags = Py_TPFLAGS_DEFAULT,
    .tp_doc = PyDoc_STR("Franka Driver, implemented in C++"),
    .tp_methods = CylinderTree_methods,
    .tp_members = CylinderTree_members,
    .tp_getset = CylinderTree_getsetters,
    .tp_init = (initproc) CylinderTree_init,
    .tp_new = (newfunc) CylinderTree_new,
};

static PyMethodDef rangetree_methods[] = {
    {NULL, NULL, 0, NULL}
};

static PyModuleDef rangetree_module = {
    PyModuleDef_HEAD_INIT,
    "ring_rangetree",
    NULL,   // Documentation
    -1,     /* size of per-interpreter state of the module,
                or -1 if the module keeps state in global variables. */
    rangetree_methods
};

PyMODINIT_FUNC
PyInit_ring_rangetree(void) {
    PyObject *m;
    if (PyType_Ready(&CylinderTreeType) < 0)
        return NULL;

    m = PyModule_Create(&rangetree_module);
    if (m == NULL)
        return NULL;

    Py_INCREF(&CylinderTreeType);
    if (PyModule_AddObject(m, "CylinderTree", (PyObject *) &CylinderTreeType) < 0) {
        Py_DECREF(&CylinderTreeType);
        Py_DECREF(m);
        return NULL;
    }

    return m;
}
