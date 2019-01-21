#ifndef MATRIX_H
#define MATRIX_H

extern "C"
void matrixAddScalar(
    double* input, 
    const size_t width, 
    const size_t height,
    const double scalar
);

extern "C"
void matrixMultiplyColVector(
    const double* input,
    double* output,
    const double* scalarVector,
    const size_t width,
    const size_t height
);

extern "C"
void matrixDivideMatrix(
    const double* divided,
    const double* divisor,
    double* output,
    const size_t width,
    const size_t height
);

extern "C"
void matrixSubstractMatrix(
    const double* input,
    const double* substrahend,
    double* output,
    const size_t width,
    const size_t height
);

extern "C"
void matrixLog10(
    double* h_input, 
    const size_t width, 
    size_t height
);

extern "C"
void matrixTranspose(
    double *h_input,
    double *h_output,
    const size_t width,
    const size_t height
);

#endif