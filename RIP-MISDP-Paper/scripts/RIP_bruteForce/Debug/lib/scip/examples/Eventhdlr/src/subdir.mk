################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
C_SRCS += \
../lib/scip/examples/Eventhdlr/src/cmain.c \
../lib/scip/examples/Eventhdlr/src/event_bestsol.c \
../lib/scip/examples/Eventhdlr/src/event_boundwriting.c 

OBJS += \
./lib/scip/examples/Eventhdlr/src/cmain.o \
./lib/scip/examples/Eventhdlr/src/event_bestsol.o \
./lib/scip/examples/Eventhdlr/src/event_boundwriting.o 

C_DEPS += \
./lib/scip/examples/Eventhdlr/src/cmain.d \
./lib/scip/examples/Eventhdlr/src/event_bestsol.d \
./lib/scip/examples/Eventhdlr/src/event_boundwriting.d 


# Each subdirectory must supply rules for building sources it contributes
lib/scip/examples/Eventhdlr/src/%.o: ../lib/scip/examples/Eventhdlr/src/%.c
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C Compiler'
	gcc -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


