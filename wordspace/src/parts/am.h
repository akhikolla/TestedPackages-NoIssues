/*
 *  Compute different association measures from frequency signatures
 */

#ifndef wordspace_am_h
#define wordspace_am_h

typedef double (*am_func)(double f, double f1, double f2, double N, int sparse); 

extern int am_table_entries; /* number of AM function pointers in am_func table */
extern am_func am_table[];

double transform(double x, int method);

#endif /* wordspace_am_h */
