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
	mpic++ -I/usr/include/gsl -I"/scratch/schoenmann/MolflowLinux" -I"/scratch/schoenmann/MolflowLinux/molflowlinux_sub" -I"/scratch/schoenmann/MolflowLinux/molflowlinux_sub/include1" -I"/scratch/schoenmann/MolflowLinux/molflowlinux_sub/Header/GLApp" -I"/scratch/schoenmann/MolflowLinux/molflowlinux_sub/Header/Header_shared_sources" -I"/scratch/schoenmann/MolflowLinux/molflowlinux_sub/Header/TruncatedGaussian" -I"/scratch/schoenmann/MolflowLinux/molflowlinux_sub/Header" -I"/scratch/schoenmann/MolflowLinux/molflowlinux_sub/cereal/external/rapidjson" -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


