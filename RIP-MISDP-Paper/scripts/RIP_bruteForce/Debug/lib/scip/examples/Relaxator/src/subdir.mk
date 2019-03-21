################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
C_SRCS += \
../lib/scip/examples/Relaxator/src/cmain.c \
../lib/scip/examples/Relaxator/src/relax_lp.c \
../lib/scip/examples/Relaxator/src/relax_nlp.c 

OBJS += \
./lib/scip/examples/Relaxator/src/cmain.o \
./lib/scip/examples/Relaxator/src/relax_lp.o \
./lib/scip/examples/Relaxator/src/relax_nlp.o 

C_DEPS += \
./lib/scip/examples/Relaxator/src/cmain.d \
./lib/scip/examples/Relaxator/src/relax_lp.d \
./lib/scip/examples/Relaxator/src/relax_nlp.d 


# Each subdirectory must supply rules for building sources it contributes
lib/scip/examples/Relaxator/src/%.o: ../lib/scip/examples/Relaxator/src/%.c
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C Compiler'
	gcc -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


