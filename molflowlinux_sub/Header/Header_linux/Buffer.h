/*
 * Buffer.h
 *
 *  Created on: 30.08.2018
 *      Author: schoenmann
 */

#ifndef MOLFLOWLINUX_SUB_SOURCE_FILES_LINUX_BUFFER_H_
#define MOLFLOWLINUX_SUB_SOURCE_FILES_LINUX_BUFFER_H_


typedef unsigned char BYTE;

typedef struct {
	signed int size;
	BYTE *buffer;
}Databuff;

void importBuff(char *fileName, Databuff *buff);
void exportBuff(char *fileName, Databuff *buff);


#endif /* MOLFLOWLINUX_SUB_SOURCE_FILES_LINUX_BUFFER_H_ */
