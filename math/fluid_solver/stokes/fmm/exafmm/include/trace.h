#ifndef TRACE_H
#define TRACE_H

//! Structure for pthread based trace
struct Trace
{
    pthread_t thread;
    double    begin;
    double    end;
    int       color;
};


#endif
