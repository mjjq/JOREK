#include <signal.h>
#include <string.h>
 
volatile sig_atomic_t sigtermed = 0;
 
void term(int signum)
{
  sigtermed = 1;
}
 
void set_trap_sigterm()
{
  struct sigaction action;
  memset(&action, 0, sizeof(struct sigaction));
  action.sa_handler = term;
  sigaction(SIGTERM, &action, NULL);
}

char sigterm_called()
{
  return (sigtermed == 1);
}
