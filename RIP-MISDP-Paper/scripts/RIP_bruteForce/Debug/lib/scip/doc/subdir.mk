################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
C_SRCS += \
../lib/scip/doc/xternal.c 

OBJS += \
./lib/scip/doc/xternal.o 

C_DEPS += \
./lib/scip/doc/xternal.d 


# Each subdirectory must supply rules for building sources it contributes
lib/scip/doc/%.o: ../lib/scip/doc/%.c
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C Compiler'
	gcc -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


