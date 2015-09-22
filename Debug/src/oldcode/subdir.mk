################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/oldcode/gzstream.cpp \
../src/oldcode/ssw_cpp.cpp 

C_SRCS += \
../src/oldcode/ssw.c 

OBJS += \
./src/oldcode/gzstream.o \
./src/oldcode/ssw.o \
./src/oldcode/ssw_cpp.o 

C_DEPS += \
./src/oldcode/ssw.d 

CPP_DEPS += \
./src/oldcode/gzstream.d \
./src/oldcode/ssw_cpp.d 


# Each subdirectory must supply rules for building sources it contributes
src/oldcode/%.o: ../src/oldcode/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -D__GXX_EXPERIMENTAL_CXX0X__ -O3 -g3 -Wall -c -fmessage-length=0 -fopenmp -std=c++0x -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '

src/oldcode/%.o: ../src/oldcode/%.c
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C Compiler'
	gcc -O3 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


