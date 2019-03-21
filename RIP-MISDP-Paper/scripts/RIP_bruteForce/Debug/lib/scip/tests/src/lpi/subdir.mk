################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
C_SRCS += \
../lib/scip/tests/src/lpi/bases.c \
../lib/scip/tests/src/lpi/boundchg.c \
../lib/scip/tests/src/lpi/change.c \
../lib/scip/tests/src/lpi/matrix.c \
../lib/scip/tests/src/lpi/solve.c 

OBJS += \
./lib/scip/tests/src/lpi/bases.o \
./lib/scip/tests/src/lpi/boundchg.o \
./lib/scip/tests/src/lpi/change.o \
./lib/scip/tests/src/lpi/matrix.o \
./lib/scip/tests/src/lpi/solve.o 

C_DEPS += \
./lib/scip/tests/src/lpi/bases.d \
./lib/scip/tests/src/lpi/boundchg.d \
./lib/scip/tests/src/lpi/change.d \
./lib/scip/tests/src/lpi/matrix.d \
./lib/scip/tests/src/lpi/solve.d 


# Each subdirectory must supply rules for building sources it contributes
lib/scip/tests/src/lpi/%.o: ../lib/scip/tests/src/lpi/%.c
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C Compiler'
	gcc -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


