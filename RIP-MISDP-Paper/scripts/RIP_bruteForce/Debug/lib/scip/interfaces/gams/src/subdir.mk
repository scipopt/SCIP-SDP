################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
C_SRCS += \
../lib/scip/interfaces/gams/src/GamsSolveTrace.c \
../lib/scip/interfaces/gams/src/event_solvetrace.c \
../lib/scip/interfaces/gams/src/reader_gmo.c 

OBJS += \
./lib/scip/interfaces/gams/src/GamsSolveTrace.o \
./lib/scip/interfaces/gams/src/event_solvetrace.o \
./lib/scip/interfaces/gams/src/reader_gmo.o 

C_DEPS += \
./lib/scip/interfaces/gams/src/GamsSolveTrace.d \
./lib/scip/interfaces/gams/src/event_solvetrace.d \
./lib/scip/interfaces/gams/src/reader_gmo.d 


# Each subdirectory must supply rules for building sources it contributes
lib/scip/interfaces/gams/src/%.o: ../lib/scip/interfaces/gams/src/%.c
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C Compiler'
	gcc -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


