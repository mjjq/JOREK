/* no suffix needed as we have one function */
#include <stdint.h>
#include <stdio.h>
#include <mpi.h>
#define INTSIZE32
#  include <murge.h>
#ifdef MURGE_INTERFACE_MAJOR_VERSION
#  if MURGE_INTERFACE_MAJOR_VERSION >= 1
#    ifdef MURGE_INTERFACE_MINOR_VERSION
#      if MURGE_INTERFACE_MINOR_VERSION >= 1
#         define MURGE_WITH_GETELEMENTLIST
#      endif
#    endif
#  endif
#endif

#ifdef MURGE_WITH_GETELEMENTLIST
#define SUFFIX(name) name

#define MURGE_UserData_t SUFFIX(MURGE_UserData_t)
#define MURGE_UserData_  SUFFIX(MURGE_UserData_)
#if (defined X_ARCHpower_ibm_aix)
#define FORTRAN_CALL(nom) nom
#else
#define FORTRAN_CALL(nom) nom ## _
#endif

typedef struct MURGE_UserData_ {
  int nVertexMax;
} MURGE_UserData_t;

#undef MURGE_user_data_t
#undef MURGE_user_data_
#define VERT_PER_ELEMENT(d) d->nVertexMax
int getVertices(int i, int * idx);
#define GET_VERTICES(i, idx, d) \
  do {                          \
    getVertices(i+1, idx);	\
  } while (0)

#  include "MURGE_GetLocalElementList.c"
#endif
