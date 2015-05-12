//From somewhere in Google

#ifndef FOSTREAM_H
#define FOSTREAM_H
// declaration of a SIMPLE iostream-like class
// that is using FORTRAN to implement I/O
// This is useful when doing file I/O in a
// mixed C++/F77 program

//#include <stdlib.h>
class fostream {
public:
    //fostream(int funit=6, char* filename=NULL);
    fostream(int funit=6, const char* filename="");	
    fostream& operator <<(const char  ch);
    fostream& operator <<(const char* txt);
    void rewind();
    ~fostream();
};

#endif
