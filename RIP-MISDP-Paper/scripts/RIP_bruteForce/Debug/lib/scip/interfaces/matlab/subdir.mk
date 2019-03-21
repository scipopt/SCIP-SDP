################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
C_SRCS += \
../lib/scip/interfaces/matlab/matscip.c 

OBJS += \
./lib/scip/interfaces/matlab/matscip.o 

C_DEPS += \
./lib/scip/interfaces/matlab/matscip.d 


# Each subdirectory must supply rules for building sources it contributes
lib/scip/interfaces/matlab/%.o: ../lib/scip/interfaces/matlab/%.c
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C Compiler'
	gcc -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


