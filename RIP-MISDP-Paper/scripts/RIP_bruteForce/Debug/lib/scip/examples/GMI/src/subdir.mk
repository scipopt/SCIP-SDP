################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
C_SRCS += \
../lib/scip/examples/GMI/src/cmain.c \
../lib/scip/examples/GMI/src/sepa_gmi.c 

OBJS += \
./lib/scip/examples/GMI/src/cmain.o \
./lib/scip/examples/GMI/src/sepa_gmi.o 

C_DEPS += \
./lib/scip/examples/GMI/src/cmain.d \
./lib/scip/examples/GMI/src/sepa_gmi.d 


# Each subdirectory must supply rules for building sources it contributes
lib/scip/examples/GMI/src/%.o: ../lib/scip/examples/GMI/src/%.c
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C Compiler'
	gcc -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


