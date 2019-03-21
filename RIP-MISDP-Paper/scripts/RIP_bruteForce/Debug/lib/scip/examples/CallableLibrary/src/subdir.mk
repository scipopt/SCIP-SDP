################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
C_SRCS += \
../lib/scip/examples/CallableLibrary/src/circle.c \
../lib/scip/examples/CallableLibrary/src/gastrans.c \
../lib/scip/examples/CallableLibrary/src/string.c 

OBJS += \
./lib/scip/examples/CallableLibrary/src/circle.o \
./lib/scip/examples/CallableLibrary/src/gastrans.o \
./lib/scip/examples/CallableLibrary/src/string.o 

C_DEPS += \
./lib/scip/examples/CallableLibrary/src/circle.d \
./lib/scip/examples/CallableLibrary/src/gastrans.d \
./lib/scip/examples/CallableLibrary/src/string.d 


# Each subdirectory must supply rules for building sources it contributes
lib/scip/examples/CallableLibrary/src/%.o: ../lib/scip/examples/CallableLibrary/src/%.c
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C Compiler'
	gcc -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


