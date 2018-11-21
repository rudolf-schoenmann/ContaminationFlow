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
 * Import and export of buffer file
 */

#include "MolflowLinuxTypes.h"

#ifndef MOLFLOWLINUX_SUB_SOURCE_FILES_LINUX_BUFFER_H_
#define MOLFLOWLINUX_SUB_SOURCE_FILES_LINUX_BUFFER_H_


typedef struct {
	signed int size;
	BYTE *buff;
}Databuff;


void importBuff(char *fileName, Databuff *databuffer);
void exportBuff(char *fileName, Databuff *databuffer);
void importBuff(std::string fileName, Databuff *databuffer);
void exportBuff(std::string fileName, Databuff *databuffer);

/*
// Tried Databuff class
class Databuff{

public:
	Databuff();
	~Databuff();

	unsigned int size;
	BYTE *buff;

	void importBuff(char *fileName);
	void exportBuff(char *fileName);
	void importBuff(std::string fileName);
	void exportBuff(std::string fileName);
};
*/

#endif /* MOLFLOWLINUX_SUB_SOURCE_FILES_LINUX_BUFFER_H_ */
