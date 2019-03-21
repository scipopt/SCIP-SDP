################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
C_SRCS += \
../lib/scip/examples/Binpacking/src/branch_ryanfoster.c \
../lib/scip/examples/Binpacking/src/cmain.c \
../lib/scip/examples/Binpacking/src/cons_samediff.c \
../lib/scip/examples/Binpacking/src/pricer_binpacking.c \
../lib/scip/examples/Binpacking/src/probdata_binpacking.c \
../lib/scip/examples/Binpacking/src/reader_bpa.c \
../lib/scip/examples/Binpacking/src/vardata_binpacking.c 

OBJS += \
./lib/scip/examples/Binpacking/src/branch_ryanfoster.o \
./lib/scip/examples/Binpacking/src/cmain.o \
./lib/scip/examples/Binpacking/src/cons_samediff.o \
./lib/scip/examples/Binpacking/src/pricer_binpacking.o \
./lib/scip/examples/Binpacking/src/probdata_binpacking.o \
./lib/scip/examples/Binpacking/src/reader_bpa.o \
./lib/scip/examples/Binpacking/src/vardata_binpacking.o 

C_DEPS += \
./lib/scip/examples/Binpacking/src/branch_ryanfoster.d \
./lib/scip/examples/Binpacking/src/cmain.d \
./lib/scip/examples/Binpacking/src/cons_samediff.d \
./lib/scip/examples/Binpacking/src/pricer_binpacking.d \
./lib/scip/examples/Binpacking/src/probdata_binpacking.d \
./lib/scip/examples/Binpacking/src/reader_bpa.d \
./lib/scip/examples/Binpacking/src/vardata_binpacking.d 


# Each subdirectory must supply rules for building sources it contributes
lib/scip/examples/Binpacking/src/%.o: ../lib/scip/examples/Binpacking/src/%.c
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C Compiler'
	gcc -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


