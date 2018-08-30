//Hier muss noch eine Autor- und Lizenznachricht hin!

#ifndef MOLFLOWLINUX_SUB_HEADER_HEADER_LINUX_TYPES_H_
#define MOLFLOWLINUX_SUB_HEADER_HEADER_LINUX_TYPES_H_

#include <stdint.h>
typedef unsigned char BYTE; 
//typedef uint64_t llong;
typedef int GLint;
typedef unsigned int GLuint;
typedef unsigned long DWORD;
//typedef uint64_t sizt_t;
//typedef long long __int64;
typedef float GLfloat;
//typedef sizt_t size_type;
//typedef unsigned int size_t;



class Error {

public:
	Error(const char *message);
	const char *GetMsg();

private:
	char msg[1024];

};



#endif /* MOLFLOWLINUX_SUB_HEADER_HEADER_LINUX_TYPES_H_ */
