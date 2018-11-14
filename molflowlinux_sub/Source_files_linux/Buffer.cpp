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
 * This file contains the functions used to import and export Buffer files
 */

#include <iostream>     // std::ios, std::istream, std::cout
#include <fstream>      // std::filebuf
#include "Buffer.h"

//Import/export of buffer files with filename given as char* or std::string
//-> no need for conversion

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


