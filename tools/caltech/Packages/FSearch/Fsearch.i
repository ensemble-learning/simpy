%module Fsearch
%include "typemaps.i"
%{
#include "Fsearch.h"
%}
void setResidual(int *, int, int, int, double *, int *);
int getResidual(int, double **, double *, int *, int, int);
