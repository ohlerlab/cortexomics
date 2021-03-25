import numpy as np
cimport numpy as np
cimport cython
from cpython cimport array
import array

cdef extern from "stdlib.h":
    double drand48()
    void srand48(long int seedval)

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.nonecheck(False)
def choose_col(np.ndarray[double, ndim=1] weights, Py_ssize_t length):
    cdef Py_ssize_t idx, i
    cdef double cs
    cdef double random
    random = drand48()
    cs = 0.0
    i = 0
    while cs < random and i < length:
        cs += weights[i]
        i += 1
    return i - 1

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.nonecheck(False)
def choose_cols(np.ndarray weights,int r, int c):
    
    # cdef array.array cols = array.array('i', [0]*r)
    cdef np.ndarray cols
    cols = np.zeros((r, 1), np.int64)
    
    for i in range(r):
        cols[i] = choose_col(weights[i, :], c)
    
    return cols

# trinds = find(readgraph.transpose())
#choose_colssp(trinds[1],trinds[0],trinds[2])
@cython.boundscheck(False)
@cython.wraparound(False)
@cython.nonecheck(False)
def choose_colssp(np.ndarray rows, np.ndarray cols, np.ndarray weights):
    
    cdef int N
    cdef int startind
    cdef int rind
    cdef int row
    N = rows.shape[0]
    # cdef int row;
    # cdef int col;
    # cdef int newrow;
    # cdef int newcol;
    cdef np.ndarray choices
    choices = np.zeros((rows.max()+1, 1), np.int64)
    startind = 0    
    for i in range(N+1):    
        # print('row'+str(i))
        row = rows[i-1]+1 if i==N else rows[i]
    #when we hit a new row, process the old one
        if((row!=rows[startind])):
            rind = choose_col(weights[startind:i],i-startind)
            # print('rind'+str(rind))
            rind += startind
            # print('trind'+str(rind))
            choices[rows[startind]]=cols[rind]
            # print('assind'+str(rows[startind]))
            startind = i

    return choices