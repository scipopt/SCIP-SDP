################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
C_SRCS += \
../lib/scip/tests/src/cons/knapsack/solveknapsackapprox.c 

OBJS += \
./lib/scip/tests/src/cons/knapsack/solveknapsackapprox.o 

C_DEPS += \
./lib/scip/tests/src/cons/knapsack/solveknapsackapprox.d 


# Each subdirectory must supply rules for building sources it contributes
lib/scip/tests/src/cons/knapsack/%.o: ../lib/scip/tests/src/cons/knapsack/%.c
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C Compiler'
	gcc -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '

