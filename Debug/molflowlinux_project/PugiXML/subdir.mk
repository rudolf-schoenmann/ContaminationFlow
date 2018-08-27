################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../molflowlinux_project/PugiXML/pugixml.cpp 

OBJS += \
./molflowlinux_project/PugiXML/pugixml.o 

CPP_DEPS += \
./molflowlinux_project/PugiXML/pugixml.d 


# Each subdirectory must supply rules for building sources it contributes
molflowlinux_project/PugiXML/%.o: ../molflowlinux_project/PugiXML/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


