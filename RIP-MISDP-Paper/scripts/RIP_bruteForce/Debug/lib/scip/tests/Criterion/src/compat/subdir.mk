################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
C_SRCS += \
../lib/scip/tests/Criterion/src/compat/alloc.c \
../lib/scip/tests/Criterion/src/compat/basename.c \
../lib/scip/tests/Criterion/src/compat/mockfile.c \
../lib/scip/tests/Criterion/src/compat/pipe.c \
../lib/scip/tests/Criterion/src/compat/process.c \
../lib/scip/tests/Criterion/src/compat/processor.c \
../lib/scip/tests/Criterion/src/compat/section.c \
../lib/scip/tests/Criterion/src/compat/time.c 

OBJS += \
./lib/scip/tests/Criterion/src/compat/alloc.o \
./lib/scip/tests/Criterion/src/compat/basename.o \
./lib/scip/tests/Criterion/src/compat/mockfile.o \
./lib/scip/tests/Criterion/src/compat/pipe.o \
./lib/scip/tests/Criterion/src/compat/process.o \
./lib/scip/tests/Criterion/src/compat/processor.o \
./lib/scip/tests/Criterion/src/compat/section.o \
./lib/scip/tests/Criterion/src/compat/time.o 

C_DEPS += \
./lib/scip/tests/Criterion/src/compat/alloc.d \
./lib/scip/tests/Criterion/src/compat/basename.d \
./lib/scip/tests/Criterion/src/compat/mockfile.d \
./lib/scip/tests/Criterion/src/compat/pipe.d \
./lib/scip/tests/Criterion/src/compat/process.d \
./lib/scip/tests/Criterion/src/compat/processor.d \
./lib/scip/tests/Criterion/src/compat/section.d \
./lib/scip/tests/Criterion/src/compat/time.d 


# Each subdirectory must supply rules for building sources it contributes
lib/scip/tests/Criterion/src/compat/%.o: ../lib/scip/tests/Criterion/src/compat/%.c
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C Compiler'
	gcc -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


