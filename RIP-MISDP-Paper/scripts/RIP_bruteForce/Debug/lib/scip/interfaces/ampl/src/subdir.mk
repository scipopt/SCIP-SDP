################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
C_SRCS += \
../lib/scip/interfaces/ampl/src/cmain.c \
../lib/scip/interfaces/ampl/src/reader_nl.c 

OBJS += \
./lib/scip/interfaces/ampl/src/cmain.o \
./lib/scip/interfaces/ampl/src/reader_nl.o 

C_DEPS += \
./lib/scip/interfaces/ampl/src/cmain.d \
./lib/scip/interfaces/ampl/src/reader_nl.d 


# Each subdirectory must supply rules for building sources it contributes
lib/scip/interfaces/ampl/src/%.o: ../lib/scip/interfaces/ampl/src/%.c
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C Compiler'
	gcc -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


