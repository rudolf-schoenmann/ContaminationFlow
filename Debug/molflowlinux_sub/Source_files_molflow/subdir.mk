################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../molflowlinux_sub/Source_files_molflow/MolflowTypes.cpp \
../molflowlinux_sub/Source_files_molflow/Parameter.cpp \
../molflowlinux_sub/Source_files_molflow/SimulationControl.cpp \
../molflowlinux_sub/Source_files_molflow/SimulationMC.cpp \
../molflowlinux_sub/Source_files_molflow/molflowSub.cpp 

OBJS += \
./molflowlinux_sub/Source_files_molflow/MolflowTypes.o \
./molflowlinux_sub/Source_files_molflow/Parameter.o \
./molflowlinux_sub/Source_files_molflow/SimulationControl.o \
./molflowlinux_sub/Source_files_molflow/SimulationMC.o \
./molflowlinux_sub/Source_files_molflow/molflowSub.o 

CPP_DEPS += \
./molflowlinux_sub/Source_files_molflow/MolflowTypes.d \
./molflowlinux_sub/Source_files_molflow/Parameter.d \
./molflowlinux_sub/Source_files_molflow/SimulationControl.d \
./molflowlinux_sub/Source_files_molflow/SimulationMC.d \
./molflowlinux_sub/Source_files_molflow/molflowSub.d 


# Each subdirectory must supply rules for building sources it contributes
molflowlinux_sub/Source_files_molflow/%.o: ../molflowlinux_sub/Source_files_molflow/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	mpic++ -std=c++0x -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


