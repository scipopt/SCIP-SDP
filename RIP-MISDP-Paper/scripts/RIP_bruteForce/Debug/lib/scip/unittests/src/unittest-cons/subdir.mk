################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
C_SRCS += \
../lib/scip/unittests/src/unittest-cons/cons_unittest.c \
../lib/scip/unittests/src/unittest-cons/unittest-cons.c 

OBJS += \
./lib/scip/unittests/src/unittest-cons/cons_unittest.o \
./lib/scip/unittests/src/unittest-cons/unittest-cons.o 

C_DEPS += \
./lib/scip/unittests/src/unittest-cons/cons_unittest.d \
./lib/scip/unittests/src/unittest-cons/unittest-cons.d 


# Each subdirectory must supply rules for building sources it contributes
lib/scip/unittests/src/unittest-cons/%.o: ../lib/scip/unittests/src/unittest-cons/%.c
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C Compiler'
	gcc -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


