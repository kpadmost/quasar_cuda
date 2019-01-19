#ifndef FEFIT_H
#define FEFIT_H

extern "C"
void reduceFeChisq(
    double *h_chisq,
    const size_t *h_sizes,
    const size_t width
);

#endif