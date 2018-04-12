#include <stdio.h>
#include <stdint.h>

#if (defined X_ARCHpower_ibm_aix)
#define FORTRAN_CALL(nom) nom
#else
#define FORTRAN_CALL(nom) nom ## _
#endif

unsigned long memAllocGetCurrent ();

void FORTRAN_CALL(pastix_getmem)(int64_t * mem) {
#ifdef MEMORY_USAGE
  *mem = memAllocGetCurrent();
#endif
}
