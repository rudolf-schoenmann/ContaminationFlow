/*
 * Databuffer.h
 *
 *  Created on: 20.08.2018
 *      Author: schoenmann
 */
//ier muss noch eine Autor- und Lizenznachricht hin!


#ifndef MOLFLOWLINUX_SUB_HEADER_HEADER_LINUX_DATABUFFER_H_
#define MOLFLOWLINUX_SUB_HEADER_HEADER_LINUX_DATABUFFER_H_

typedef struct{
	sizt_t size; //keep track of buffer size
	void *buff; // View handle
} Databuffer;



#endif /* MOLFLOWLINUX_SUB_HEADER_HEADER_LINUX_DATABUFFER_H_ */
