################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../ContaminationFlow_sub/boost/spirit/home/support/char_encoding/unicode/create_tables.cpp 

OBJS += \
./ContaminationFlow_sub/boost/spirit/home/support/char_encoding/unicode/create_tables.o 

CPP_DEPS += \
./ContaminationFlow_sub/boost/spirit/home/support/char_encoding/unicode/create_tables.d 


# Each subdirectory must supply rules for building sources it contributes
ContaminationFlow_sub/boost/spirit/home/support/char_encoding/unicode/%.o: ../ContaminationFlow_sub/boost/spirit/home/support/char_encoding/unicode/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	mpic++ -std=c++0x -std=c++11 -D__GXX_EXPERIMENTAL_CXX0X__ -I"../include/include" -I".." -I"../ContaminationFlow_sub" -I"../ContaminationFlow_sub/Header/Header_shared_sources" -I"../ContaminationFlow_sub/Header" -I"../ContaminationFlow_sub/cereal/external/rapidjson" -I"../ContaminationFlow_sub/Header/Header_ContaminationFlow" -I"../ContaminationFlow_sub/Header/Header_molflow" -I"../ContaminationFlow_sub/TruncatedGaussian" -I"../ContaminationFlow_sub/boost" -I"../ContaminationFlow_sub/GLApp" -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


