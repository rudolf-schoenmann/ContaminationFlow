################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../molflowlinux_sub/include2/delayhlp.cpp 

OBJS += \
./molflowlinux_sub/include2/delayhlp.o 

CPP_DEPS += \
./molflowlinux_sub/include2/delayhlp.d 


# Each subdirectory must supply rules for building sources it contributes
molflowlinux_sub/include2/%.o: ../molflowlinux_sub/include2/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	mpic++ -I"/scratch/schoenmann/MolflowLinux" -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


