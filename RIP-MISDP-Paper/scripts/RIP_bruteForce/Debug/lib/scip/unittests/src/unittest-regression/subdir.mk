################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
C_SRCS += \
../lib/scip/unittests/src/unittest-regression/unittest-regression.c 

OBJS += \
./lib/scip/unittests/src/unittest-regression/unittest-regression.o 

C_DEPS += \
./lib/scip/unittests/src/unittest-regression/unittest-regression.d 


# Each subdirectory must supply rules for building sources it contributes
lib/scip/unittests/src/unittest-regression/%.o: ../lib/scip/unittests/src/unittest-regression/%.c
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C Compiler'
	gcc -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


