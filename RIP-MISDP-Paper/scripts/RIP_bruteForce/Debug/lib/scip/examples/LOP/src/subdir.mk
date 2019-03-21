################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
C_SRCS += \
../lib/scip/examples/LOP/src/cmain.c \
../lib/scip/examples/LOP/src/cons_lop.c \
../lib/scip/examples/LOP/src/genRandomLOPInstance.c \
../lib/scip/examples/LOP/src/reader_lop.c 

OBJS += \
./lib/scip/examples/LOP/src/cmain.o \
./lib/scip/examples/LOP/src/cons_lop.o \
./lib/scip/examples/LOP/src/genRandomLOPInstance.o \
./lib/scip/examples/LOP/src/reader_lop.o 

C_DEPS += \
./lib/scip/examples/LOP/src/cmain.d \
./lib/scip/examples/LOP/src/cons_lop.d \
./lib/scip/examples/LOP/src/genRandomLOPInstance.d \
./lib/scip/examples/LOP/src/reader_lop.d 


# Each subdirectory must supply rules for building sources it contributes
lib/scip/examples/LOP/src/%.o: ../lib/scip/examples/LOP/src/%.c
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C Compiler'
	gcc -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


