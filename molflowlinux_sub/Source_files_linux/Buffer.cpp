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
#include "SimulationLinux.h"
#include <cstring>		//std::memcpy



// check if files writeable/readable
bool checkReadable(std::string fileName,int rank){
	std::filebuf fb;
	if (fb.open(fileName, std::ios::in))
	{
		return true;
	}
	std::ostringstream tmpstream (std::ostringstream::app);
	tmpstream <<home_to_tilde(fileName) <<" not readable" <<std::endl;
	if(rank==0)
		printStream(tmpstream.str());
	return false;

}

bool checkWriteable(std::string fileName, int rank){
	std::filebuf fb;
	if (fb.open(fileName, std::ios::app))
	{
		return true;
	}
	std::ostringstream tmpstream (std::ostringstream::app);
	tmpstream <<home_to_tilde(fileName) <<" not writable" <<std::endl;
	if(rank==0)
		printStream(tmpstream.str());
	return false;

}

//-----------------------------------------------------------
//Import/export of buffer files with filename given as char* or std::string

void importBuff(const char *fileName, Databuff *databuffer)
{
	std::filebuf fb;
	databuffer->buff=NULL;
	if (fb.open(fileName, std::ios::in))
	{
		std::istream is(&fb);
		if (is) {
			// Get length of file
			is.seekg(0, is.end);
			signed int length = is.tellg();
			is.seekg(0, is.beg);

			// Read file
			char *temp = new char[length];
			is.read(temp, length);

			// Set length and content of buffer
			databuffer->size = length;
			databuffer->buff = (BYTE*)temp;
		}
		fb.close();
	}
	else{
		std::cout << "Could not open " << home_to_tilde(fileName) <<" to read data from." << std::endl;
	}
	std::cout << "Buffer '" << home_to_tilde(fileName) <<"' imported. Buffersize (read in): " << databuffer->size << std::endl;
}




void exportBuff(const char *fileName, Databuff *databuffer)
{
	std::filebuf f;
	if (f.open(fileName, std::ios::out))
	{
		std::ostream os(&f);
		if (os)
		{
			// Write buffer to file
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
	const char *c =fileName.c_str(); // Convert std::string to const char*
	importBuff(c, databuffer);
}




void exportBuff(std::string fileName, Databuff *databuffer)
{
	const char *c =fileName.c_str(); // Convert std::string to const char*
	exportBuff(c, databuffer);
}
