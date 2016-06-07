#include <unistd.h>
void fgetpid(int *id)
{
  *id = (int)getpid();
}
