################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../molflowlinux_project/Source_Files_Linux/MolflowLinux.cpp \
../molflowlinux_project/Source_Files_Linux/emptyGeometry.cpp \
../molflowlinux_project/Source_Files_Linux/load.cpp \
../molflowlinux_project/Source_Files_Linux/loadConfig.cpp 

OBJS += \
./molflowlinux_project/Source_Files_Linux/MolflowLinux.o \
./molflowlinux_project/Source_Files_Linux/emptyGeometry.o \
./molflowlinux_project/Source_Files_Linux/load.o \
./molflowlinux_project/Source_Files_Linux/loadConfig.o 

CPP_DEPS += \
./molflowlinux_project/Source_Files_Linux/MolflowLinux.d \
./molflowlinux_project/Source_Files_Linux/emptyGeometry.d \
./molflowlinux_project/Source_Files_Linux/load.d \
./molflowlinux_project/Source_Files_Linux/loadConfig.d 


# Each subdirectory must supply rules for building sources it contributes
molflowlinux_project/Source_Files_Linux/%.o: ../molflowlinux_project/Source_Files_Linux/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


