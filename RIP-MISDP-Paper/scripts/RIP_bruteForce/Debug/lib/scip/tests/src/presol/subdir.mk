################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
C_SRCS += \
../lib/scip/tests/src/presol/presol.c \
../lib/scip/tests/src/presol/qpkktref.c 

OBJS += \
./lib/scip/tests/src/presol/presol.o \
./lib/scip/tests/src/presol/qpkktref.o 

C_DEPS += \
./lib/scip/tests/src/presol/presol.d \
./lib/scip/tests/src/presol/qpkktref.d 


# Each subdirectory must supply rules for building sources it contributes
lib/scip/tests/src/presol/%.o: ../lib/scip/tests/src/presol/%.c
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C Compiler'
	gcc -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


