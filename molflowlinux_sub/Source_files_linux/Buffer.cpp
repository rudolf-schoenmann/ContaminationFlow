/*
 * Buffer.cpp
 *
 *  Created on: 30.08.2018
 *      Author: schoenmann
 */

#include <iostream>     // std::ios, std::istream, std::cout
#include <fstream>      // std::filebuf

typedef unsigned char BYTE;

typedef struct {
	signed int size;
	BYTE *buffer;
}Databuff;

void importBuff(char *fileName, Databuff *buff)
{

	std::filebuf fb;
	if (fb.open(fileName, std::ios::in))
	{

		std::istream is(&fb);
		if (is) {
			is.seekg(0, is.end);
			signed int length = is.tellg();
			is.seekg(0, is.beg);

			char *temp = new char[length];
			buff->buffer = new BYTE[length];

			is.read(temp, length);
			buff->size = length;
			buff->buffer = (BYTE*)temp;
		}

		fb.close();
	}
	else{std::cout << "Could not open file to read data from." << std::endl;}
	std::cout << "Buffersize (read in): " << buff->size << std::endl;
}

void exportBuff(char *fileName, Databuff *buff)
{
	std::filebuf f;
	if (f.open(fileName, std::ios::out))
	{
		std::ostream os(&f);
		if (os)
		{
			signed int length = buff->size;
			os.write(reinterpret_cast<const char *>(buff->buffer), length);
		}
		f.close();
	}
	else {std::cout << "Could not open file to write data in." << std::endl;}
	std::cout << "Buffersize (write in file): " << buff->size << std::endl;
}




/* Test main function
 *
 * int main() {
	char fileimport[] = "C:/Users/van/Documents/buffer";
	char fileexport[] = "C:/Users/van/Documents/buffertest";
	Databuff buff;
	buff.buffer = NULL;

	importBuff(fileimport, &buff);
	exportBuff(fileexport, &buff);

	delete[] buff.buffer;
	return 0;
}*/


