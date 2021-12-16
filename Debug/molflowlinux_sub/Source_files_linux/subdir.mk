################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../molflowlinux_sub/Source_files_linux/Buffer.cpp \
../molflowlinux_sub/Source_files_linux/Iteration.cpp \
../molflowlinux_sub/Source_files_linux/SimulationCalc.cpp \
../molflowlinux_sub/Source_files_linux/SimulationLinux.cpp \
../molflowlinux_sub/Source_files_linux/UpdateMainProcess.cpp \
../molflowlinux_sub/Source_files_linux/UpdateSubProcess.cpp \
../molflowlinux_sub/Source_files_linux/worker.cpp 

OBJS += \
./molflowlinux_sub/Source_files_linux/Buffer.o \
./molflowlinux_sub/Source_files_linux/Iteration.o \
./molflowlinux_sub/Source_files_linux/SimulationCalc.o \
./molflowlinux_sub/Source_files_linux/SimulationLinux.o \
./molflowlinux_sub/Source_files_linux/UpdateMainProcess.o \
./molflowlinux_sub/Source_files_linux/UpdateSubProcess.o \
./molflowlinux_sub/Source_files_linux/worker.o 

CPP_DEPS += \
./molflowlinux_sub/Source_files_linux/Buffer.d \
./molflowlinux_sub/Source_files_linux/Iteration.d \
./molflowlinux_sub/Source_files_linux/SimulationCalc.d \
./molflowlinux_sub/Source_files_linux/SimulationLinux.d \
./molflowlinux_sub/Source_files_linux/UpdateMainProcess.d \
./molflowlinux_sub/Source_files_linux/UpdateSubProcess.d \
./molflowlinux_sub/Source_files_linux/worker.d 


# Each subdirectory must supply rules for building sources it contributes
molflowlinux_sub/Source_files_linux/%.o: ../molflowlinux_sub/Source_files_linux/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	mpic++ -std=c++0x -std=c++11 -D__GXX_EXPERIMENTAL_CXX0X__ -I"../include/include" -I".." -I"../molflowlinux_sub" -I"../molflowlinux_sub/Header/Header_shared_sources" -I"../molflowlinux_sub/Header" -I"../molflowlinux_sub/cereal/external/rapidjson" -I"../molflowlinux_sub/Header/Header_linux" -I"../molflowlinux_sub/Header/Header_molflow" -I"../molflowlinux_sub/TruncatedGaussian" -I"../molflowlinux_sub/boost" -I"../molflowlinux_sub/GLApp" -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


