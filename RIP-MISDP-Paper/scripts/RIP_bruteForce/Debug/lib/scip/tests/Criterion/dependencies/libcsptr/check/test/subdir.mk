################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
C_SRCS += \
../lib/scip/tests/Criterion/dependencies/libcsptr/check/test/array.c \
../lib/scip/tests/Criterion/dependencies/libcsptr/check/test/misc.c \
../lib/scip/tests/Criterion/dependencies/libcsptr/check/test/scalar.c \
../lib/scip/tests/Criterion/dependencies/libcsptr/check/test/shared.c \
../lib/scip/tests/Criterion/dependencies/libcsptr/check/test/test.c 

OBJS += \
./lib/scip/tests/Criterion/dependencies/libcsptr/check/test/array.o \
./lib/scip/tests/Criterion/dependencies/libcsptr/check/test/misc.o \
./lib/scip/tests/Criterion/dependencies/libcsptr/check/test/scalar.o \
./lib/scip/tests/Criterion/dependencies/libcsptr/check/test/shared.o \
./lib/scip/tests/Criterion/dependencies/libcsptr/check/test/test.o 

C_DEPS += \
./lib/scip/tests/Criterion/dependencies/libcsptr/check/test/array.d \
./lib/scip/tests/Criterion/dependencies/libcsptr/check/test/misc.d \
./lib/scip/tests/Criterion/dependencies/libcsptr/check/test/scalar.d \
./lib/scip/tests/Criterion/dependencies/libcsptr/check/test/shared.d \
./lib/scip/tests/Criterion/dependencies/libcsptr/check/test/test.d 


# Each subdirectory must supply rules for building sources it contributes
lib/scip/tests/Criterion/dependencies/libcsptr/check/test/%.o: ../lib/scip/tests/Criterion/dependencies/libcsptr/check/test/%.c
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C Compiler'
	gcc -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


