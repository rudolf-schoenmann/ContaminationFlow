/*
 * Buffer.cpp
 *
 *  Created on: 30.08.2018
 *      Author: schoenmann
 */

#include <iostream>     // std::ios, std::istream, std::cout
#include <fstream>      // std::filebuf
#include "Buffer.h"



void importBuff(char *fileName, Databuff *databuffer)
{

	std::filebuf fb;
	databuffer->buff=NULL;
	if (fb.open(fileName, std::ios::in))
	{

		std::istream is(&fb);
		if (is) {
			is.seekg(0, is.end);
			signed int length = is.tellg();
			is.seekg(0, is.beg);

			char *temp = new char[length];
			databuffer->buff = new BYTE[length];

			is.read(temp, length);
			databuffer->size = length;
			databuffer->buff = (BYTE*)temp;
		}

		fb.close();
	}
	else{
		std::cout << "Could not open file to read data from." << std::endl;
	}
	std::cout << "Buffer '" << fileName <<"' imported. Buffersize (read in): " << databuffer->size << std::endl;
}




void exportBuff(char *fileName, Databuff *databuffer)
{
	std::filebuf f;
	if (f.open(fileName, std::ios::out))
	{
		std::ostream os(&f);
		if (os)
		{
			signed int length = databuffer->size;
			os.write(reinterpret_cast<const char *>(databuffer->buff), length);
		}
		f.close();
	}
	else {std::cout << "Could not open file to write data in." << std::endl;}
	std::cout << "Buffersize (write in file): " << databuffer->size << std::endl;
}

void importBuff(std::string fileName, Databuff *databuffer)
{

	std::filebuf fb;
	databuffer->buff=NULL;
	if (fb.open(fileName, std::ios::in))
	{

		std::istream is(&fb);
		if (is) {
			is.seekg(0, is.end);
			signed int length = is.tellg();
			is.seekg(0, is.beg);

			char *temp = new char[length];
			databuffer->buff = new BYTE[length];

			is.read(temp, length);
			databuffer->size = length;
			databuffer->buff = (BYTE*)temp;
		}

		fb.close();
	}
	else{
		std::cout << "Could not open file to read data from." << std::endl;
	}
	std::cout << "Buffer '" << fileName <<"' imported. Buffersize (read in): " << databuffer->size << std::endl;
}




void exportBuff(std::string fileName, Databuff *databuffer)
{
	std::filebuf f;
	if (f.open(fileName, std::ios::out))
	{
		std::ostream os(&f);
		if (os)
		{
			signed int length = databuffer->size;
			os.write(reinterpret_cast<const char *>(databuffer->buff), length);
		}
		f.close();
	}
	else {std::cout << "Could not open file to write data in." << std::endl;}
	std::cout << "Buffersize (write in file): " << databuffer->size << std::endl;
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


