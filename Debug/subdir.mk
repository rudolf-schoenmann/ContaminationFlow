################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../molflowlinux_main.cpp 

OBJS += \
./molflowlinux_main.o 

CPP_DEPS += \
./molflowlinux_main.d 


# Each subdirectory must supply rules for building sources it contributes
%.o: ../%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	mpic++ -std=c++0x -I/usr/include/gsl -I"/home/van/MolflowLinux" -I"/home/van/MolflowLinux/molflowlinux_sub" -I"/home/van/MolflowLinux/molflowlinux_sub/include1" -I"/home/van/MolflowLinux/molflowlinux_sub/Header/Header_shared_sources" -I"/home/van/MolflowLinux/molflowlinux_sub/Header" -I"/home/van/MolflowLinux/molflowlinux_sub/cereal/external/rapidjson" -I"/home/van/MolflowLinux/molflowlinux_sub/Header/Header_linux" -I"/home/van/MolflowLinux/molflowlinux_sub/Header/Header_molflow" -I"/home/van/MolflowLinux/molflowlinux_sub/TruncatedGaussian" -I"/home/van/MolflowLinux/molflowlinux_sub/GLApp" -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


