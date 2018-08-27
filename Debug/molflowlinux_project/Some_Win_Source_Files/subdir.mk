################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../molflowlinux_project/Some_Win_Source_Files/Geometry_shared.cpp \
../molflowlinux_project/Some_Win_Source_Files/GlobalSettings.cpp \
../molflowlinux_project/Some_Win_Source_Files/Interface.cpp \
../molflowlinux_project/Some_Win_Source_Files/MolFlow.cpp \
../molflowlinux_project/Some_Win_Source_Files/MolflowGeometry.cpp \
../molflowlinux_project/Some_Win_Source_Files/MolflowWorker.cpp \
../molflowlinux_project/Some_Win_Source_Files/Process.cpp \
../molflowlinux_project/Some_Win_Source_Files/Random.cpp \
../molflowlinux_project/Some_Win_Source_Files/ShMemory.cpp \
../molflowlinux_project/Some_Win_Source_Files/Vector.cpp \
../molflowlinux_project/Some_Win_Source_Files/Worker_shared.cpp \
../molflowlinux_project/Some_Win_Source_Files/test.cpp 

OBJS += \
./molflowlinux_project/Some_Win_Source_Files/Geometry_shared.o \
./molflowlinux_project/Some_Win_Source_Files/GlobalSettings.o \
./molflowlinux_project/Some_Win_Source_Files/Interface.o \
./molflowlinux_project/Some_Win_Source_Files/MolFlow.o \
./molflowlinux_project/Some_Win_Source_Files/MolflowGeometry.o \
./molflowlinux_project/Some_Win_Source_Files/MolflowWorker.o \
./molflowlinux_project/Some_Win_Source_Files/Process.o \
./molflowlinux_project/Some_Win_Source_Files/Random.o \
./molflowlinux_project/Some_Win_Source_Files/ShMemory.o \
./molflowlinux_project/Some_Win_Source_Files/Vector.o \
./molflowlinux_project/Some_Win_Source_Files/Worker_shared.o \
./molflowlinux_project/Some_Win_Source_Files/test.o 

CPP_DEPS += \
./molflowlinux_project/Some_Win_Source_Files/Geometry_shared.d \
./molflowlinux_project/Some_Win_Source_Files/GlobalSettings.d \
./molflowlinux_project/Some_Win_Source_Files/Interface.d \
./molflowlinux_project/Some_Win_Source_Files/MolFlow.d \
./molflowlinux_project/Some_Win_Source_Files/MolflowGeometry.d \
./molflowlinux_project/Some_Win_Source_Files/MolflowWorker.d \
./molflowlinux_project/Some_Win_Source_Files/Process.d \
./molflowlinux_project/Some_Win_Source_Files/Random.d \
./molflowlinux_project/Some_Win_Source_Files/ShMemory.d \
./molflowlinux_project/Some_Win_Source_Files/Vector.d \
./molflowlinux_project/Some_Win_Source_Files/Worker_shared.d \
./molflowlinux_project/Some_Win_Source_Files/test.d 


# Each subdirectory must supply rules for building sources it contributes
molflowlinux_project/Some_Win_Source_Files/%.o: ../molflowlinux_project/Some_Win_Source_Files/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


