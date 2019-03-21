################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
C_SRCS += \
../lib/scip/examples/Eventhdlr/doc/xternal_eventhdlr.c 

OBJS += \
./lib/scip/examples/Eventhdlr/doc/xternal_eventhdlr.o 

C_DEPS += \
./lib/scip/examples/Eventhdlr/doc/xternal_eventhdlr.d 


# Each subdirectory must supply rules for building sources it contributes
lib/scip/examples/Eventhdlr/doc/%.o: ../lib/scip/examples/Eventhdlr/doc/%.c
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C Compiler'
	gcc -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


