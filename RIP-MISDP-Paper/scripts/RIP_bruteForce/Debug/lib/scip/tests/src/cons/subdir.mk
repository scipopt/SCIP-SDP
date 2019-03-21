################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
C_SRCS += \
../lib/scip/tests/src/cons/cons.c \
../lib/scip/tests/src/cons/enforelaximplemented.c \
../lib/scip/tests/src/cons/initlp.c 

OBJS += \
./lib/scip/tests/src/cons/cons.o \
./lib/scip/tests/src/cons/enforelaximplemented.o \
./lib/scip/tests/src/cons/initlp.o 

C_DEPS += \
./lib/scip/tests/src/cons/cons.d \
./lib/scip/tests/src/cons/enforelaximplemented.d \
./lib/scip/tests/src/cons/initlp.d 


# Each subdirectory must supply rules for building sources it contributes
lib/scip/tests/src/cons/%.o: ../lib/scip/tests/src/cons/%.c
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C Compiler'
	gcc -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


