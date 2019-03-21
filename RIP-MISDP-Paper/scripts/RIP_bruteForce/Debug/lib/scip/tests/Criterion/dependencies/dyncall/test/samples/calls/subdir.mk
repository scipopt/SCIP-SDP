################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
C_SRCS += \
../lib/scip/tests/Criterion/dependencies/dyncall/test/samples/calls/d16.c \
../lib/scip/tests/Criterion/dependencies/dyncall/test/samples/calls/d20.c \
../lib/scip/tests/Criterion/dependencies/dyncall/test/samples/calls/d40.c \
../lib/scip/tests/Criterion/dependencies/dyncall/test/samples/calls/f16.c \
../lib/scip/tests/Criterion/dependencies/dyncall/test/samples/calls/f20.c \
../lib/scip/tests/Criterion/dependencies/dyncall/test/samples/calls/f40.c \
../lib/scip/tests/Criterion/dependencies/dyncall/test/samples/calls/fd40.c \
../lib/scip/tests/Criterion/dependencies/dyncall/test/samples/calls/i3.c \
../lib/scip/tests/Criterion/dependencies/dyncall/test/samples/calls/i7.c \
../lib/scip/tests/Criterion/dependencies/dyncall/test/samples/calls/id40.c \
../lib/scip/tests/Criterion/dependencies/dyncall/test/samples/calls/l16.c \
../lib/scip/tests/Criterion/dependencies/dyncall/test/samples/calls/many.c 

OBJS += \
./lib/scip/tests/Criterion/dependencies/dyncall/test/samples/calls/d16.o \
./lib/scip/tests/Criterion/dependencies/dyncall/test/samples/calls/d20.o \
./lib/scip/tests/Criterion/dependencies/dyncall/test/samples/calls/d40.o \
./lib/scip/tests/Criterion/dependencies/dyncall/test/samples/calls/f16.o \
./lib/scip/tests/Criterion/dependencies/dyncall/test/samples/calls/f20.o \
./lib/scip/tests/Criterion/dependencies/dyncall/test/samples/calls/f40.o \
./lib/scip/tests/Criterion/dependencies/dyncall/test/samples/calls/fd40.o \
./lib/scip/tests/Criterion/dependencies/dyncall/test/samples/calls/i3.o \
./lib/scip/tests/Criterion/dependencies/dyncall/test/samples/calls/i7.o \
./lib/scip/tests/Criterion/dependencies/dyncall/test/samples/calls/id40.o \
./lib/scip/tests/Criterion/dependencies/dyncall/test/samples/calls/l16.o \
./lib/scip/tests/Criterion/dependencies/dyncall/test/samples/calls/many.o 

C_DEPS += \
./lib/scip/tests/Criterion/dependencies/dyncall/test/samples/calls/d16.d \
./lib/scip/tests/Criterion/dependencies/dyncall/test/samples/calls/d20.d \
./lib/scip/tests/Criterion/dependencies/dyncall/test/samples/calls/d40.d \
./lib/scip/tests/Criterion/dependencies/dyncall/test/samples/calls/f16.d \
./lib/scip/tests/Criterion/dependencies/dyncall/test/samples/calls/f20.d \
./lib/scip/tests/Criterion/dependencies/dyncall/test/samples/calls/f40.d \
./lib/scip/tests/Criterion/dependencies/dyncall/test/samples/calls/fd40.d \
./lib/scip/tests/Criterion/dependencies/dyncall/test/samples/calls/i3.d \
./lib/scip/tests/Criterion/dependencies/dyncall/test/samples/calls/i7.d \
./lib/scip/tests/Criterion/dependencies/dyncall/test/samples/calls/id40.d \
./lib/scip/tests/Criterion/dependencies/dyncall/test/samples/calls/l16.d \
./lib/scip/tests/Criterion/dependencies/dyncall/test/samples/calls/many.d 


# Each subdirectory must supply rules for building sources it contributes
lib/scip/tests/Criterion/dependencies/dyncall/test/samples/calls/%.o: ../lib/scip/tests/Criterion/dependencies/dyncall/test/samples/calls/%.c
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C Compiler'
	gcc -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


