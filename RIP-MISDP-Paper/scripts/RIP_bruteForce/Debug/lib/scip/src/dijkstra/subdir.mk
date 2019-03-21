################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
C_SRCS += \
../lib/scip/src/dijkstra/dijkstra.c 

OBJS += \
./lib/scip/src/dijkstra/dijkstra.o 

C_DEPS += \
./lib/scip/src/dijkstra/dijkstra.d 


# Each subdirectory must supply rules for building sources it contributes
lib/scip/src/dijkstra/%.o: ../lib/scip/src/dijkstra/%.c
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C Compiler'
	gcc -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


