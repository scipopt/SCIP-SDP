################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
C_SRCS += \
../lib/scip/unittests/src/unittest-setters/unittest-setters.c 

OBJS += \
./lib/scip/unittests/src/unittest-setters/unittest-setters.o 

C_DEPS += \
./lib/scip/unittests/src/unittest-setters/unittest-setters.d 


# Each subdirectory must supply rules for building sources it contributes
lib/scip/unittests/src/unittest-setters/%.o: ../lib/scip/unittests/src/unittest-setters/%.c
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C Compiler'
	gcc -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


