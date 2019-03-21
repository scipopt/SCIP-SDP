################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
C_SRCS += \
../lib/scip/tests/Criterion/samples/asserts.c \
../lib/scip/tests/Criterion/samples/description.c \
../lib/scip/tests/Criterion/samples/fixtures.c \
../lib/scip/tests/Criterion/samples/more-suites.c \
../lib/scip/tests/Criterion/samples/parameterized.c \
../lib/scip/tests/Criterion/samples/redirect.c \
../lib/scip/tests/Criterion/samples/report.c \
../lib/scip/tests/Criterion/samples/signal.c \
../lib/scip/tests/Criterion/samples/simple.c \
../lib/scip/tests/Criterion/samples/suites.c \
../lib/scip/tests/Criterion/samples/theories.c \
../lib/scip/tests/Criterion/samples/timeout.c 

OBJS += \
./lib/scip/tests/Criterion/samples/asserts.o \
./lib/scip/tests/Criterion/samples/description.o \
./lib/scip/tests/Criterion/samples/fixtures.o \
./lib/scip/tests/Criterion/samples/more-suites.o \
./lib/scip/tests/Criterion/samples/parameterized.o \
./lib/scip/tests/Criterion/samples/redirect.o \
./lib/scip/tests/Criterion/samples/report.o \
./lib/scip/tests/Criterion/samples/signal.o \
./lib/scip/tests/Criterion/samples/simple.o \
./lib/scip/tests/Criterion/samples/suites.o \
./lib/scip/tests/Criterion/samples/theories.o \
./lib/scip/tests/Criterion/samples/timeout.o 

C_DEPS += \
./lib/scip/tests/Criterion/samples/asserts.d \
./lib/scip/tests/Criterion/samples/description.d \
./lib/scip/tests/Criterion/samples/fixtures.d \
./lib/scip/tests/Criterion/samples/more-suites.d \
./lib/scip/tests/Criterion/samples/parameterized.d \
./lib/scip/tests/Criterion/samples/redirect.d \
./lib/scip/tests/Criterion/samples/report.d \
./lib/scip/tests/Criterion/samples/signal.d \
./lib/scip/tests/Criterion/samples/simple.d \
./lib/scip/tests/Criterion/samples/suites.d \
./lib/scip/tests/Criterion/samples/theories.d \
./lib/scip/tests/Criterion/samples/timeout.d 


# Each subdirectory must supply rules for building sources it contributes
lib/scip/tests/Criterion/samples/%.o: ../lib/scip/tests/Criterion/samples/%.c
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C Compiler'
	gcc -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


