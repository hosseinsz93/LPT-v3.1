#include "variables.h"

void initlist(List *ilist) {
  ilist->head = PETSC_NULL;
  ilist->number = 0;
}

void insertnode(List *ilist, Location Node)
{
  node *tnew;
  node *current;
  current = ilist->head;

  PetscTruth Exist = PETSC_FALSE;
  while(current) {
    if (Node.i == current->Node.i &&
	Node.j == current->Node.j &&
	Node.k == current->Node.k) {
      Exist = PETSC_TRUE;
    }
    if (Exist) break;
    current = current->next;
  }
  if (!Exist) {
    PetscMalloc(sizeof(node), &tnew);
    tnew->next = ilist->head;
    tnew->Node = Node;
    ilist->head = tnew;
    ilist->number++;
  }
}

void destroy(List *ilist)
{
  node *current;
  while (ilist->head) {
    current = ilist->head->next;
    PetscFree(ilist->head);
    ilist->head = current;
  }
}
