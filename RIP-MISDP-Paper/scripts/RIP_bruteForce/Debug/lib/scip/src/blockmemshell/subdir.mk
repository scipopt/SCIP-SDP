################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
C_SRCS += \
../lib/scip/src/blockmemshell/memory.c 

OBJS += \
./lib/scip/src/blockmemshell/memory.o 

C_DEPS += \
./lib/scip/src/blockmemshell/memory.d 


# Each subdirectory must supply rules for building sources it contributes
lib/scip/src/blockmemshell/%.o: ../lib/scip/src/blockmemshell/%.c
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C Compiler'
	gcc -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


