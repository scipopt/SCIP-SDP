################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
C_SRCS += \
../lib/scip/applications/STP/doc/xternal.c 

OBJS += \
./lib/scip/applications/STP/doc/xternal.o 

C_DEPS += \
./lib/scip/applications/STP/doc/xternal.d 


# Each subdirectory must supply rules for building sources it contributes
lib/scip/applications/STP/doc/%.o: ../lib/scip/applications/STP/doc/%.c
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C Compiler'
	gcc -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


