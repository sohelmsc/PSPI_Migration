################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
O_SRCS += \
../src/Migration.o \
../src/SeismicImaging.o 

CPP_SRCS += \
../src/Migration.cpp \
../src/SeismicImaging.cpp 

OBJS += \
./src/Migration.o \
./src/SeismicImaging.o 

CPP_DEPS += \
./src/Migration.d \
./src/SeismicImaging.d 


# Each subdirectory must supply rules for building sources it contributes
src/%.o: ../src/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -O0 -g3 -Wall -Wextra -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


