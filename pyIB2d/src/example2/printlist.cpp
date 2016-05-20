#include "printlist.h"

using namespace std;

void printlist(list<int> &l){
    for(list<int>::const_iterator i = l.begin(); i != l.end(); i++)
    cout << *i << ' ';
    cout << endl;}
