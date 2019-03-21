################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
C_SRCS += \
../lib/scip/examples/TSP/doc/xternal_tsp.c 

OBJS += \
./lib/scip/examples/TSP/doc/xternal_tsp.o 

C_DEPS += \
./lib/scip/examples/TSP/doc/xternal_tsp.d 


# Each subdirectory must supply rules for building sources it contributes
lib/scip/examples/TSP/doc/%.o: ../lib/scip/examples/TSP/doc/%.c
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C Compiler'
	gcc -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


