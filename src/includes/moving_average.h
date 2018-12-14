#ifndef MOVING_AVERAGE_H
#define MOVING_AVERAGE_H

extern "C"
void movingAverage( // centered moving average
    const double *h_input,
    double *h_output,
    const uint *cols,
    const size_t width,
    const size_t height,
    const uint window
);

#endif