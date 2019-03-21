################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
C_SRCS += \
../lib/scip/tests/Criterion/src/core/abort.c \
../lib/scip/tests/Criterion/src/core/ordered-set.c \
../lib/scip/tests/Criterion/src/core/report.c \
../lib/scip/tests/Criterion/src/core/runner.c \
../lib/scip/tests/Criterion/src/core/runner_coroutine.c \
../lib/scip/tests/Criterion/src/core/stats.c \
../lib/scip/tests/Criterion/src/core/test.c \
../lib/scip/tests/Criterion/src/core/theories.c \
../lib/scip/tests/Criterion/src/core/worker.c 

OBJS += \
./lib/scip/tests/Criterion/src/core/abort.o \
./lib/scip/tests/Criterion/src/core/ordered-set.o \
./lib/scip/tests/Criterion/src/core/report.o \
./lib/scip/tests/Criterion/src/core/runner.o \
./lib/scip/tests/Criterion/src/core/runner_coroutine.o \
./lib/scip/tests/Criterion/src/core/stats.o \
./lib/scip/tests/Criterion/src/core/test.o \
./lib/scip/tests/Criterion/src/core/theories.o \
./lib/scip/tests/Criterion/src/core/worker.o 

C_DEPS += \
./lib/scip/tests/Criterion/src/core/abort.d \
./lib/scip/tests/Criterion/src/core/ordered-set.d \
./lib/scip/tests/Criterion/src/core/report.d \
./lib/scip/tests/Criterion/src/core/runner.d \
./lib/scip/tests/Criterion/src/core/runner_coroutine.d \
./lib/scip/tests/Criterion/src/core/stats.d \
./lib/scip/tests/Criterion/src/core/test.d \
./lib/scip/tests/Criterion/src/core/theories.d \
./lib/scip/tests/Criterion/src/core/worker.d 


# Each subdirectory must supply rules for building sources it contributes
lib/scip/tests/Criterion/src/core/%.o: ../lib/scip/tests/Criterion/src/core/%.c
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C Compiler'
	gcc -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


