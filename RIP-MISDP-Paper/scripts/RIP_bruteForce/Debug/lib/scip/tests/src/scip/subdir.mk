################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
C_SRCS += \
../lib/scip/tests/src/scip/bilinenvelope.c \
../lib/scip/tests/src/scip/intervalarith.c \
../lib/scip/tests/src/scip/probingobj.c \
../lib/scip/tests/src/scip/setters.c 

OBJS += \
./lib/scip/tests/src/scip/bilinenvelope.o \
./lib/scip/tests/src/scip/intervalarith.o \
./lib/scip/tests/src/scip/probingobj.o \
./lib/scip/tests/src/scip/setters.o 

C_DEPS += \
./lib/scip/tests/src/scip/bilinenvelope.d \
./lib/scip/tests/src/scip/intervalarith.d \
./lib/scip/tests/src/scip/probingobj.d \
./lib/scip/tests/src/scip/setters.d 


# Each subdirectory must supply rules for building sources it contributes
lib/scip/tests/src/scip/%.o: ../lib/scip/tests/src/scip/%.c
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C Compiler'
	gcc -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


