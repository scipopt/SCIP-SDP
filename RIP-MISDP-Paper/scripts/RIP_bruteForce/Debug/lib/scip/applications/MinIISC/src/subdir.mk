################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
C_SRCS += \
../lib/scip/applications/MinIISC/src/benders.c \
../lib/scip/applications/MinIISC/src/miniisc.c \
../lib/scip/applications/MinIISC/src/readargs.c 

OBJS += \
./lib/scip/applications/MinIISC/src/benders.o \
./lib/scip/applications/MinIISC/src/miniisc.o \
./lib/scip/applications/MinIISC/src/readargs.o 

C_DEPS += \
./lib/scip/applications/MinIISC/src/benders.d \
./lib/scip/applications/MinIISC/src/miniisc.d \
./lib/scip/applications/MinIISC/src/readargs.d 


# Each subdirectory must supply rules for building sources it contributes
lib/scip/applications/MinIISC/src/%.o: ../lib/scip/applications/MinIISC/src/%.c
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C Compiler'
	gcc -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


