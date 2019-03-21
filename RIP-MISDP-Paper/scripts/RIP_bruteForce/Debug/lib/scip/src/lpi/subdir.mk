################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
C_SRCS += \
../lib/scip/src/lpi/lpi_cpx.c \
../lib/scip/src/lpi/lpi_grb.c \
../lib/scip/src/lpi/lpi_msk.c \
../lib/scip/src/lpi/lpi_none.c \
../lib/scip/src/lpi/lpi_qso.c \
../lib/scip/src/lpi/lpi_xprs.c 

OBJS += \
./lib/scip/src/lpi/lpi_cpx.o \
./lib/scip/src/lpi/lpi_grb.o \
./lib/scip/src/lpi/lpi_msk.o \
./lib/scip/src/lpi/lpi_none.o \
./lib/scip/src/lpi/lpi_qso.o \
./lib/scip/src/lpi/lpi_xprs.o 

C_DEPS += \
./lib/scip/src/lpi/lpi_cpx.d \
./lib/scip/src/lpi/lpi_grb.d \
./lib/scip/src/lpi/lpi_msk.d \
./lib/scip/src/lpi/lpi_none.d \
./lib/scip/src/lpi/lpi_qso.d \
./lib/scip/src/lpi/lpi_xprs.d 


# Each subdirectory must supply rules for building sources it contributes
lib/scip/src/lpi/%.o: ../lib/scip/src/lpi/%.c
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C Compiler'
	gcc -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


