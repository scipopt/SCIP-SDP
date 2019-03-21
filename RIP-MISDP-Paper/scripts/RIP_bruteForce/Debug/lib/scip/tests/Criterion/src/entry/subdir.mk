################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
C_SRCS += \
../lib/scip/tests/Criterion/src/entry/entry.c \
../lib/scip/tests/Criterion/src/entry/options.c \
../lib/scip/tests/Criterion/src/entry/params.c 

OBJS += \
./lib/scip/tests/Criterion/src/entry/entry.o \
./lib/scip/tests/Criterion/src/entry/options.o \
./lib/scip/tests/Criterion/src/entry/params.o 

C_DEPS += \
./lib/scip/tests/Criterion/src/entry/entry.d \
./lib/scip/tests/Criterion/src/entry/options.d \
./lib/scip/tests/Criterion/src/entry/params.d 


# Each subdirectory must supply rules for building sources it contributes
lib/scip/tests/Criterion/src/entry/%.o: ../lib/scip/tests/Criterion/src/entry/%.c
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C Compiler'
	gcc -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


