//
// element_rtree.cpp
//
// This is a direct port of the C version of the RTree test program.
//

#include "RTree.h"

extern "C" { // prevent name mangling
using namespace std;

typedef int ValueType;
typedef RTree<ValueType, double, 2, double> MyTree;
// Persistent tree
static MyTree ElementTree;

void PopulateTree(int n, double minx[], double miny[], double maxx[], double maxy[])
{
  int i;
  double min[2], max[2];
  ElementTree.RemoveAll();
  for(i=0; i<n; i++)
  {
    min[0] = minx[i];
    min[1] = miny[i];
    max[0] = maxx[i];
    max[1] = maxy[i];
    ElementTree.Insert(min, max, i+1); // store element number (1-based)
  }
}

// Return the number of elements in a rectangle
int NumElementsInRect(double minx, double miny, double maxx, double maxy)
{
  double min[2], max[2];
  min[0] = minx;
  min[1] = miny;
  max[0] = maxx;
  max[1] = maxy;
  return ElementTree.Search(min, max, NULL, NULL);
}

// Return element indices of elements contained within the rectangle in element_tree
// i_elm must be allocated by the caller to size at least nelm.
int ElementsInRect(double minx, double miny, double maxx, double maxy, int *ielm)
{
  double min[2], max[2];
  min[0] = minx;
  min[1] = miny;
  max[0] = maxx;
  max[1] = maxy;
  return ElementTree.Search(min, max, NULL, ielm);
}
}
