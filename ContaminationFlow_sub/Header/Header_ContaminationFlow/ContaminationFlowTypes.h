/*
Program:     ContaminationFlow
Description: Monte Carlo simulator for satellite contanimation studies
Authors:     Rudolf Schönmann / Hoai My Van
Copyright:   TU Munich
Forked from: Molflow (CERN) (https://cern.ch/molflow)

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

Full license text: https://www.gnu.org/licenses/old-licenses/gpl-2.0.en.html
*/

/*
 * User defined data types
 */

#ifndef MOLFLOWLINUX_SUB_HEADER_HEADER_LINUX_TYPES_H_
#define MOLFLOWLINUX_SUB_HEADER_HEADER_LINUX_TYPES_H_

#include <cstdint>
#define __int64 int64_t
#define  llong uint64_t

#define size_t uint64_t

typedef unsigned char BYTE; 
typedef int GLint;
typedef unsigned int GLuint;
typedef unsigned long DWORD;
//typedef uint64_t sizt_t;
//typedef long long __int64;
typedef float GLfloat;
//typedef sizt_t size_type;
//typedef unsigned int size_t;


/*
class Error {

public:
	Error(const char *message);
	const char *GetMsg();

private:
	char msg[1024];

};
*/


#endif /* MOLFLOWLINUX_SUB_HEADER_HEADER_LINUX_TYPES_H_ */
