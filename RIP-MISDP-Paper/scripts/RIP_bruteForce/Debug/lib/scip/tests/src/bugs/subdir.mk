################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
C_SRCS += \
../lib/scip/tests/src/bugs/depthlevel.c 

OBJS += \
./lib/scip/tests/src/bugs/depthlevel.o 

C_DEPS += \
./lib/scip/tests/src/bugs/depthlevel.d 


# Each subdirectory must supply rules for building sources it contributes
lib/scip/tests/src/bugs/%.o: ../lib/scip/tests/src/bugs/%.c
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C Compiler'
	gcc -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


