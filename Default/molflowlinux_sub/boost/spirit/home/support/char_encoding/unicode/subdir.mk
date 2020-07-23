################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../molflowlinux_sub/boost/spirit/home/support/char_encoding/unicode/create_tables.cpp 

OBJS += \
./molflowlinux_sub/boost/spirit/home/support/char_encoding/unicode/create_tables.o 

CPP_DEPS += \
./molflowlinux_sub/boost/spirit/home/support/char_encoding/unicode/create_tables.d 


# Each subdirectory must supply rules for building sources it contributes
molflowlinux_sub/boost/spirit/home/support/char_encoding/unicode/%.o: ../molflowlinux_sub/boost/spirit/home/support/char_encoding/unicode/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	mpic++ -std=c++0x -std=c++11 -D__GXX_EXPERIMENTAL_CXX0X__ -I"/scratch/schoenmann/MolflowLinux" -I"/scratch/schoenmann/MolflowLinux/include/include/gsl" -I"/scratch/schoenmann/MolflowLinux/molflowlinux_sub" -I"/scratch/schoenmann/MolflowLinux/molflowlinux_sub/Header/Header_shared_sources" -I"/scratch/schoenmann/MolflowLinux/molflowlinux_sub/Header" -I"/scratch/schoenmann/MolflowLinux/molflowlinux_sub/cereal/external/rapidjson" -I"/scratch/schoenmann/MolflowLinux/molflowlinux_sub/Header/Header_linux" -I"/scratch/schoenmann/MolflowLinux/molflowlinux_sub/Header/Header_molflow" -I"/scratch/schoenmann/MolflowLinux/molflowlinux_sub/TruncatedGaussian" -I"/scratch/schoenmann/MolflowLinux/molflowlinux_sub/boost" -I"/scratch/schoenmann/MolflowLinux/molflowlinux_sub/GLApp" -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


