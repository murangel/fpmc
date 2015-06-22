#include <iostream>
#include "fostream.h"
#include <cstring>

extern"C" {
  void f77opn_(int&,char*,int);
  void f77out_(char*,int);
  void f77cls_();
  void f77rwd_();
}

#define f77opn f77opn_
#define f77out f77out_
#define f77cls f77cls_
#define f77rwd f77rwd_

fostream::fostream(int funit, const char* filename)
{
   if(filename){
      f77opn(funit,(char*)filename,strlen(filename));
   } else{
      const char* filename = "stdout.txt";
      f77opn(funit,(char*)filename,strlen(filename));
   }
}

fostream::~fostream()
{
  f77cls();
}

fostream& fostream::operator <<(const char ch)
{
   char str[2] = { ch,'\0'};
   f77out(str,strlen(str));
   return *this;
}

fostream& fostream::operator <<(const char* txt)
{
   char* mytxt = (char*)txt;
   f77out(mytxt,strlen(mytxt)); 
   return *this;
}

void fostream::rewind()
{
   f77rwd();
}
