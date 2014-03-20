#include <stdio.h>
#include <stl.h>
using namespace std;
string getHomePath()
{
	char *pPath;
  	pPath = getenv ("SEQBMAHOME");
  	if (pPath!=NULL)
		return pPath;
	else 
		return "";

}



