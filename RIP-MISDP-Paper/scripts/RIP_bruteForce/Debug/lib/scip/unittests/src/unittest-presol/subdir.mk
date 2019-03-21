################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
C_SRCS += \
../lib/scip/unittests/src/unittest-presol/presol_unittest.c \
../lib/scip/unittests/src/unittest-presol/unittest-presol.c 

OBJS += \
./lib/scip/unittests/src/unittest-presol/presol_unittest.o \
./lib/scip/unittests/src/unittest-presol/unittest-presol.o 

C_DEPS += \
./lib/scip/unittests/src/unittest-presol/presol_unittest.d \
./lib/scip/unittests/src/unittest-presol/unittest-presol.d 


# Each subdirectory must supply rules for building sources it contributes
lib/scip/unittests/src/unittest-presol/%.o: ../lib/scip/unittests/src/unittest-presol/%.c
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C Compiler'
	gcc -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


