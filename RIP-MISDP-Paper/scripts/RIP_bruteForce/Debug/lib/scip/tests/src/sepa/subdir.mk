################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
C_SRCS += \
../lib/scip/tests/src/sepa/convexproj.c \
../lib/scip/tests/src/sepa/gauge.c 

OBJS += \
./lib/scip/tests/src/sepa/convexproj.o \
./lib/scip/tests/src/sepa/gauge.o 

C_DEPS += \
./lib/scip/tests/src/sepa/convexproj.d \
./lib/scip/tests/src/sepa/gauge.d 


# Each subdirectory must supply rules for building sources it contributes
lib/scip/tests/src/sepa/%.o: ../lib/scip/tests/src/sepa/%.c
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C Compiler'
	gcc -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


