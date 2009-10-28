!! pre-processor macro for handling memory allocation
!! \author Magnus Hagdorn
!! \date October 2009

#ifdef HAVE_CONFIG_H
#include<config.inc>
#endif

#define GLIMMER_ALLOC1D(array,n1) allocate(array(n1),stat=merr); call glimmer_allocErr("(n1)","array",merr,__FILE__,__LINE__)
#define GLIMMER_ALLOC2D(array,n1,n2) allocate(array(n1,n2),stat=merr); call glimmer_allocErr("(n1,n2)","array",merr,__FILE__,__LINE__)
#define GLIMMER_ALLOC3D(array,n1,n2,n3) allocate(array(n1,n2,n3),stat=merr); call glimmer_allocErr("(n1,n2,n3)","array",merr,__FILE__,__LINE__)
#define GLIMMER_ALLOC4D(array,n1,n2,n3,n4) allocate(array(n1,n2,n3,n4),stat=merr); call glimmer_allocErr("(n1,n2,n3,n4)","array",merr,__FILE__,__LINE__)
#define GLIMMER_DEALLOC(array) deallocate(array,stat=merr); call glimmer_deallocErr("array",merr,__FILE__,__LINE__)
