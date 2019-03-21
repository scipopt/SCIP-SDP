################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
C_SRCS += \
../lib/scip/tests/src/cons/quadratic/gauge.c \
../lib/scip/tests/src/cons/quadratic/projection.c 

OBJS += \
./lib/scip/tests/src/cons/quadratic/gauge.o \
./lib/scip/tests/src/cons/quadratic/projection.o 

C_DEPS += \
./lib/scip/tests/src/cons/quadratic/gauge.d \
./lib/scip/tests/src/cons/quadratic/projection.d 


# Each subdirectory must supply rules for building sources it contributes
lib/scip/tests/src/cons/quadratic/%.o: ../lib/scip/tests/src/cons/quadratic/%.c
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C Compiler'
	gcc -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


