################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
C_SRCS += \
../lib/scip/tests/src/cons/nonlinear/getCoeffsAndConstantFromLinearExpr.c \
../lib/scip/tests/src/cons/nonlinear/issue1326.c 

OBJS += \
./lib/scip/tests/src/cons/nonlinear/getCoeffsAndConstantFromLinearExpr.o \
./lib/scip/tests/src/cons/nonlinear/issue1326.o 

C_DEPS += \
./lib/scip/tests/src/cons/nonlinear/getCoeffsAndConstantFromLinearExpr.d \
./lib/scip/tests/src/cons/nonlinear/issue1326.d 


# Each subdirectory must supply rules for building sources it contributes
lib/scip/tests/src/cons/nonlinear/%.o: ../lib/scip/tests/src/cons/nonlinear/%.c
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C Compiler'
	gcc -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


