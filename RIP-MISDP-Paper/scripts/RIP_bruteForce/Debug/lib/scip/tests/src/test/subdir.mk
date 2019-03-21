################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
C_SRCS += \
../lib/scip/tests/src/test/stages.c 

OBJS += \
./lib/scip/tests/src/test/stages.o 

C_DEPS += \
./lib/scip/tests/src/test/stages.d 


# Each subdirectory must supply rules for building sources it contributes
lib/scip/tests/src/test/%.o: ../lib/scip/tests/src/test/%.c
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C Compiler'
	gcc -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


