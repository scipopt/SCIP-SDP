################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
C_SRCS += \
../lib/scip/tests/Criterion/src/io/asprintf.c \
../lib/scip/tests/Criterion/src/io/event.c \
../lib/scip/tests/Criterion/src/io/file.c \
../lib/scip/tests/Criterion/src/io/json.c \
../lib/scip/tests/Criterion/src/io/output.c \
../lib/scip/tests/Criterion/src/io/redirect.c \
../lib/scip/tests/Criterion/src/io/tap.c \
../lib/scip/tests/Criterion/src/io/xml.c 

OBJS += \
./lib/scip/tests/Criterion/src/io/asprintf.o \
./lib/scip/tests/Criterion/src/io/event.o \
./lib/scip/tests/Criterion/src/io/file.o \
./lib/scip/tests/Criterion/src/io/json.o \
./lib/scip/tests/Criterion/src/io/output.o \
./lib/scip/tests/Criterion/src/io/redirect.o \
./lib/scip/tests/Criterion/src/io/tap.o \
./lib/scip/tests/Criterion/src/io/xml.o 

C_DEPS += \
./lib/scip/tests/Criterion/src/io/asprintf.d \
./lib/scip/tests/Criterion/src/io/event.d \
./lib/scip/tests/Criterion/src/io/file.d \
./lib/scip/tests/Criterion/src/io/json.d \
./lib/scip/tests/Criterion/src/io/output.d \
./lib/scip/tests/Criterion/src/io/redirect.d \
./lib/scip/tests/Criterion/src/io/tap.d \
./lib/scip/tests/Criterion/src/io/xml.d 


# Each subdirectory must supply rules for building sources it contributes
lib/scip/tests/Criterion/src/io/%.o: ../lib/scip/tests/Criterion/src/io/%.c
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C Compiler'
	gcc -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


