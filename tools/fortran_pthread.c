#include <pthread.h>
#include <stdlib.h>
#include <stdio.h>

#define FORTRAN_NAME(nu,nl,pl,pc) \
void nu pl;                       \
void nl pl                        \
{ nu pc; }                        \
void nl##_ pl                     \
{ nu pc; }                        \
void nl##__ pl                    \
{ nu pc; }                        \
void nu pl

FORTRAN_NAME(FORTRAN_PTHREAD_CREATE,
	     fortran_pthread_create,
	     (long * thread_ptr, void *(*start_routine)(void *), void * arg, int * ret),
	     (thread_ptr, start_routine, arg, ret)) {
  pthread_t *thread = (pthread_t*)malloc(sizeof(pthread_t));
  *thread_ptr = (long)thread;
  *ret = pthread_create((pthread_t *)thread, 
		 NULL, 
		 start_routine,
		 arg);
}

FORTRAN_NAME(FORTRAN_PTHREAD_JOIN,
	     fortran_pthread_join,
	     (long * thread_ptr, int * value_ptr, int * ret),
	     (thread_ptr, value_ptr, ret)) {
  pthread_t * thread = (pthread_t*)(*thread_ptr);
  
  *ret = pthread_join((pthread_t)(*thread), (void**)value_ptr);
  free(thread);
}

FORTRAN_NAME(FORTRAN_PTHREAD_MUTEX_INIT,
	     fortran_pthread_mutex_init,
	     (long *mutex_ptr, int* ret),
	     (mutex_ptr, ret)) {
  pthread_mutex_t*mutex = (pthread_mutex_t*)malloc(sizeof(pthread_mutex_t));

  *mutex_ptr = (long)mutex;
  *ret = pthread_mutex_init((pthread_mutex_t *)mutex, NULL);
}

FORTRAN_NAME(FORTRAN_PTHREAD_MUTEX_DESTROY,
	     fortran_pthread_mutex_destroy,
	     (long *mutex_ptr, int * ret),
	     (mutex_ptr, ret)) {
  pthread_mutex_t * mutex = (pthread_mutex_t*)(*mutex_ptr);
  *ret = pthread_mutex_destroy((pthread_mutex_t *)mutex);
  free(mutex);
}

FORTRAN_NAME(FORTRAN_PTHREAD_MUTEX_LOCK,
	     fortran_pthread_mutex_lock,
	     (long *mutex_ptr, int * ret),
	     (mutex_ptr, ret)) {
  pthread_mutex_t * mutex = (pthread_mutex_t*)(*mutex_ptr);
  *ret=pthread_mutex_lock((pthread_mutex_t *)mutex);
}
FORTRAN_NAME(FORTRAN_PTHREAD_MUTEX_TRYLOCK,
	     fortran_pthread_mutex_trylock,
	     (long *mutex_ptr, int * ret),
	     (mutex_ptr, ret)) {
  pthread_mutex_t * mutex = (pthread_mutex_t*)(*mutex_ptr);
  *ret=pthread_mutex_trylock((pthread_mutex_t *)mutex);
}
FORTRAN_NAME(FORTRAN_PTHREAD_MUTEX_UNLOCK,
	     fortran_pthread_mutex_unlock,
	     (long *mutex_ptr, int * ret),
	     (mutex_ptr, ret)) {
  pthread_mutex_t * mutex = (pthread_mutex_t*)(*mutex_ptr);
  *ret=pthread_mutex_unlock((pthread_mutex_t *)mutex);
}

FORTRAN_NAME(FORTRAN_PTHREAD_COND_INIT,
	     fortran_pthread_cond_init,
	     (long *cond_ptr, int*ret),
	     (cond_ptr, ret)) {
  pthread_cond_t *cond = malloc(sizeof(pthread_cond_t));
  *cond_ptr = (long)(cond);
  *ret = pthread_cond_init(cond, NULL);
}

FORTRAN_NAME(FORTRAN_PTHREAD_COND_BROADCAST,
	     fortran_pthread_cond_broadcast,
	     (long *cond_ptr, int * ret),
	     (cond_ptr, ret)) {

  pthread_cond_t *cond = (pthread_cond_t*)cond_ptr;
  *ret = pthread_cond_broadcast(cond);
}

FORTRAN_NAME(FORTRAN_PTHREAD_COND_WAIT,
	     fortran_pthread_cond_wait,
	     (long *cond_ptr, pthread_mutex_t *mutex, int* ret),
	     (cond_ptr, mutex, ret)) {
  *ret = pthread_cond_wait((pthread_cond_t*)(cond_ptr), mutex);
}

FORTRAN_NAME(FORTRAN_PTHREAD_COND_DESTROY,
	     fortran_pthread_cond_destroy,
	     (long *cond_ptr, int * ret),
	     (cond_ptr, ret)) {
  *ret = pthread_cond_destroy((pthread_cond_t*)(cond_ptr));
  free((void*)(*cond_ptr));
}

typedef struct mypthread_barrier_ {
  int volatile    instance;         /*+ ID of the barrier                +*/
  int volatile    blocked_threads;  /*+ Number of threads in the barrier +*/
  pthread_mutex_t sync_lock;        /*+ mutex for the barrier            +*/
  pthread_cond_t  sync_cond;        /*+ cond for the barrier             +*/
} mypthread_barrier_t;


FORTRAN_NAME(INIT_BARRIER,
	     init_barrier,
	     (long * barrier_ptr),
	     (barrier_ptr)) {
  mypthread_barrier_t * barrier = (mypthread_barrier_t * )malloc(sizeof(mypthread_barrier_t));
  barrier->instance = 0;
  barrier->blocked_threads = 0;
  pthread_cond_init(&(barrier->sync_cond), NULL);
  pthread_mutex_init(&(barrier->sync_lock), NULL);
  *barrier_ptr = (long)barrier;
  
}
FORTRAN_NAME(DESTROY_BARRIER,
	     destroy_barrier,
	     (long * barrier_ptr), 
	     (barrier_ptr)) {
  mypthread_barrier_t * barrier = (mypthread_barrier_t * )(*barrier_ptr);
  pthread_cond_destroy(&(barrier->sync_cond));
  pthread_mutex_destroy(&(barrier->sync_lock));
  free(barrier);
}

FORTRAN_NAME(SYNCHRO_X_THREADS,
	     synchro_x_threads,
	     (int * nbthreads, long * barrier_ptr),
	     (nbthreads, barrier_ptr)) {
  int instance;							
  mypthread_barrier_t * barrier = (mypthread_barrier_t * )(*barrier_ptr);
  pthread_mutex_lock(&((barrier)->sync_lock));				
  instance = (barrier)->instance;					
  (barrier)->blocked_threads++;	
  if ((barrier)->blocked_threads == (*nbthreads))			
    {									
      (barrier)->blocked_threads = 0;					
      (barrier)->instance++;						
      pthread_cond_broadcast(&((barrier)->sync_cond));			
    }									
  while (instance == (barrier)->instance)				
    {									
      pthread_cond_wait(&((barrier)->sync_cond), &((barrier)->sync_lock));
    }									
  pthread_mutex_unlock(&((barrier)->sync_lock));			
}
