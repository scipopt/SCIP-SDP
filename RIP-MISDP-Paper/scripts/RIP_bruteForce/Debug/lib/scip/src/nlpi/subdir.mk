################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
C_SRCS += \
../lib/scip/src/nlpi/expr.c \
../lib/scip/src/nlpi/exprinterpret_none.c \
../lib/scip/src/nlpi/nlpi.c \
../lib/scip/src/nlpi/nlpi_ipopt_dummy.c \
../lib/scip/src/nlpi/nlpi_xyz.c \
../lib/scip/src/nlpi/nlpioracle.c 

OBJS += \
./lib/scip/src/nlpi/expr.o \
./lib/scip/src/nlpi/exprinterpret_none.o \
./lib/scip/src/nlpi/nlpi.o \
./lib/scip/src/nlpi/nlpi_ipopt_dummy.o \
./lib/scip/src/nlpi/nlpi_xyz.o \
./lib/scip/src/nlpi/nlpioracle.o 

C_DEPS += \
./lib/scip/src/nlpi/expr.d \
./lib/scip/src/nlpi/exprinterpret_none.d \
./lib/scip/src/nlpi/nlpi.d \
./lib/scip/src/nlpi/nlpi_ipopt_dummy.d \
./lib/scip/src/nlpi/nlpi_xyz.d \
./lib/scip/src/nlpi/nlpioracle.d 


# Each subdirectory must supply rules for building sources it contributes
lib/scip/src/nlpi/%.o: ../lib/scip/src/nlpi/%.c
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C Compiler'
	gcc -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


