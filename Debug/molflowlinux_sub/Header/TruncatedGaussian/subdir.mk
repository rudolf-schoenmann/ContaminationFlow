################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../molflowlinux_sub/Header/TruncatedGaussian/main.cpp \
../molflowlinux_sub/Header/TruncatedGaussian/rtnorm.cpp 

OBJS += \
./molflowlinux_sub/Header/TruncatedGaussian/main.o \
./molflowlinux_sub/Header/TruncatedGaussian/rtnorm.o 

CPP_DEPS += \
./molflowlinux_sub/Header/TruncatedGaussian/main.d \
./molflowlinux_sub/Header/TruncatedGaussian/rtnorm.d 


# Each subdirectory must supply rules for building sources it contributes
molflowlinux_sub/Header/TruncatedGaussian/%.o: ../molflowlinux_sub/Header/TruncatedGaussian/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	mpic++ -std=c++0x -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


