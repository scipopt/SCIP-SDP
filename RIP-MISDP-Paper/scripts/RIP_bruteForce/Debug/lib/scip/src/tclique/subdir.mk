################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
C_SRCS += \
../lib/scip/src/tclique/tclique_branch.c \
../lib/scip/src/tclique/tclique_coloring.c \
../lib/scip/src/tclique/tclique_graph.c 

OBJS += \
./lib/scip/src/tclique/tclique_branch.o \
./lib/scip/src/tclique/tclique_coloring.o \
./lib/scip/src/tclique/tclique_graph.o 

C_DEPS += \
./lib/scip/src/tclique/tclique_branch.d \
./lib/scip/src/tclique/tclique_coloring.d \
./lib/scip/src/tclique/tclique_graph.d 


# Each subdirectory must supply rules for building sources it contributes
lib/scip/src/tclique/%.o: ../lib/scip/src/tclique/%.c
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C Compiler'
	gcc -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


