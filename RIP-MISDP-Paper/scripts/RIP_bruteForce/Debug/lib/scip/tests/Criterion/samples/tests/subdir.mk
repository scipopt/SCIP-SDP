################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
C_SRCS += \
../lib/scip/tests/Criterion/samples/tests/exit.c \
../lib/scip/tests/Criterion/samples/tests/failmessages.c \
../lib/scip/tests/Criterion/samples/tests/long-messages.c \
../lib/scip/tests/Criterion/samples/tests/other-crashes.c \
../lib/scip/tests/Criterion/samples/tests/theories_regression.c \
../lib/scip/tests/Criterion/samples/tests/with-time.c 

OBJS += \
./lib/scip/tests/Criterion/samples/tests/exit.o \
./lib/scip/tests/Criterion/samples/tests/failmessages.o \
./lib/scip/tests/Criterion/samples/tests/long-messages.o \
./lib/scip/tests/Criterion/samples/tests/other-crashes.o \
./lib/scip/tests/Criterion/samples/tests/theories_regression.o \
./lib/scip/tests/Criterion/samples/tests/with-time.o 

C_DEPS += \
./lib/scip/tests/Criterion/samples/tests/exit.d \
./lib/scip/tests/Criterion/samples/tests/failmessages.d \
./lib/scip/tests/Criterion/samples/tests/long-messages.d \
./lib/scip/tests/Criterion/samples/tests/other-crashes.d \
./lib/scip/tests/Criterion/samples/tests/theories_regression.d \
./lib/scip/tests/Criterion/samples/tests/with-time.d 


# Each subdirectory must supply rules for building sources it contributes
lib/scip/tests/Criterion/samples/tests/%.o: ../lib/scip/tests/Criterion/samples/tests/%.c
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C Compiler'
	gcc -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


