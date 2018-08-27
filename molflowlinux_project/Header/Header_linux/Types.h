#include <stdint.h>
typedef unsigned char BYTE; 
typedef uint64_t llong;
typedef int GLint;
typedef unsigned int GLuint;
typedef unsigned long DWORD;
typedef uint64_t sizt_t;
typedef float GLfloat;
//typedef sizt_t size_type;



class Error {

public:
	Error(const char *message);
	const char *GetMsg();

private:
	char msg[1024];

};
