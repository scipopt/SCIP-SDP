################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
C_SRCS += \
../lib/scip/tests/src/misc/dbldblarith.c \
../lib/scip/tests/src/misc/normaldistribution.c \
../lib/scip/tests/src/misc/rbtreetest.c \
../lib/scip/tests/src/misc/regression.c \
../lib/scip/tests/src/misc/select.c 

OBJS += \
./lib/scip/tests/src/misc/dbldblarith.o \
./lib/scip/tests/src/misc/normaldistribution.o \
./lib/scip/tests/src/misc/rbtreetest.o \
./lib/scip/tests/src/misc/regression.o \
./lib/scip/tests/src/misc/select.o 

C_DEPS += \
./lib/scip/tests/src/misc/dbldblarith.d \
./lib/scip/tests/src/misc/normaldistribution.d \
./lib/scip/tests/src/misc/rbtreetest.d \
./lib/scip/tests/src/misc/regression.d \
./lib/scip/tests/src/misc/select.d 


# Each subdirectory must supply rules for building sources it contributes
lib/scip/tests/src/misc/%.o: ../lib/scip/tests/src/misc/%.c
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C Compiler'
	gcc -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


