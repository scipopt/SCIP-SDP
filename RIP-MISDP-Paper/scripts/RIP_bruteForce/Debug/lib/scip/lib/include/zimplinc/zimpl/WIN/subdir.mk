################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
C_SRCS += \
../lib/scip/lib/include/zimplinc/zimpl/WIN/getopt.c 

OBJS += \
./lib/scip/lib/include/zimplinc/zimpl/WIN/getopt.o 

C_DEPS += \
./lib/scip/lib/include/zimplinc/zimpl/WIN/getopt.d 


# Each subdirectory must supply rules for building sources it contributes
lib/scip/lib/include/zimplinc/zimpl/WIN/%.o: ../lib/scip/lib/include/zimplinc/zimpl/WIN/%.c
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C Compiler'
	gcc -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


