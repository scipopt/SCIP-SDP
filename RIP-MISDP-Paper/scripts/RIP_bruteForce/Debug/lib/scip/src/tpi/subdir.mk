################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
C_SRCS += \
../lib/scip/src/tpi/tpi_none.c \
../lib/scip/src/tpi/tpi_openmp.c \
../lib/scip/src/tpi/tpi_tnycthrd.c 

OBJS += \
./lib/scip/src/tpi/tpi_none.o \
./lib/scip/src/tpi/tpi_openmp.o \
./lib/scip/src/tpi/tpi_tnycthrd.o 

C_DEPS += \
./lib/scip/src/tpi/tpi_none.d \
./lib/scip/src/tpi/tpi_openmp.d \
./lib/scip/src/tpi/tpi_tnycthrd.d 


# Each subdirectory must supply rules for building sources it contributes
lib/scip/src/tpi/%.o: ../lib/scip/src/tpi/%.c
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C Compiler'
	gcc -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


