################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../molflowlinux_sub/Source_files_molflow/MolflowTypes.cpp \
../molflowlinux_sub/Source_files_molflow/Parameter.cpp \
../molflowlinux_sub/Source_files_molflow/Simulation.cpp \
../molflowlinux_sub/Source_files_molflow/SimulationControl.cpp \
../molflowlinux_sub/Source_files_molflow/SimulationMC.cpp 

OBJS += \
./molflowlinux_sub/Source_files_molflow/MolflowTypes.o \
./molflowlinux_sub/Source_files_molflow/Parameter.o \
./molflowlinux_sub/Source_files_molflow/Simulation.o \
./molflowlinux_sub/Source_files_molflow/SimulationControl.o \
./molflowlinux_sub/Source_files_molflow/SimulationMC.o 

CPP_DEPS += \
./molflowlinux_sub/Source_files_molflow/MolflowTypes.d \
./molflowlinux_sub/Source_files_molflow/Parameter.d \
./molflowlinux_sub/Source_files_molflow/Simulation.d \
./molflowlinux_sub/Source_files_molflow/SimulationControl.d \
./molflowlinux_sub/Source_files_molflow/SimulationMC.d 


# Each subdirectory must supply rules for building sources it contributes
molflowlinux_sub/Source_files_molflow/%.o: ../molflowlinux_sub/Source_files_molflow/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	mpic++ -std=c++0x -std=c++11 -D__GXX_EXPERIMENTAL_CXX0X__ -I"../include/include" -I".." -I"../molflowlinux_sub" -I"../molflowlinux_sub/Header/Header_shared_sources" -I"../molflowlinux_sub/Header" -I"../molflowlinux_sub/cereal/external/rapidjson" -I"../molflowlinux_sub/Header/Header_linux" -I"../molflowlinux_sub/Header/Header_molflow" -I"../molflowlinux_sub/TruncatedGaussian" -I"../molflowlinux_sub/boost" -I"../molflowlinux_sub/GLApp" -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


