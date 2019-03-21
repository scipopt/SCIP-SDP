################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
C_SRCS += \
../lib/scip/applications/Coloring/src/branch_coloring.c \
../lib/scip/applications/Coloring/src/branch_strongcoloring.c \
../lib/scip/applications/Coloring/src/coloringplugins.c \
../lib/scip/applications/Coloring/src/cons_storeGraph.c \
../lib/scip/applications/Coloring/src/heur_init.c \
../lib/scip/applications/Coloring/src/main.c \
../lib/scip/applications/Coloring/src/pricer_coloring.c \
../lib/scip/applications/Coloring/src/probdata_coloring.c \
../lib/scip/applications/Coloring/src/reader_col.c \
../lib/scip/applications/Coloring/src/reader_csol.c 

OBJS += \
./lib/scip/applications/Coloring/src/branch_coloring.o \
./lib/scip/applications/Coloring/src/branch_strongcoloring.o \
./lib/scip/applications/Coloring/src/coloringplugins.o \
./lib/scip/applications/Coloring/src/cons_storeGraph.o \
./lib/scip/applications/Coloring/src/heur_init.o \
./lib/scip/applications/Coloring/src/main.o \
./lib/scip/applications/Coloring/src/pricer_coloring.o \
./lib/scip/applications/Coloring/src/probdata_coloring.o \
./lib/scip/applications/Coloring/src/reader_col.o \
./lib/scip/applications/Coloring/src/reader_csol.o 

C_DEPS += \
./lib/scip/applications/Coloring/src/branch_coloring.d \
./lib/scip/applications/Coloring/src/branch_strongcoloring.d \
./lib/scip/applications/Coloring/src/coloringplugins.d \
./lib/scip/applications/Coloring/src/cons_storeGraph.d \
./lib/scip/applications/Coloring/src/heur_init.d \
./lib/scip/applications/Coloring/src/main.d \
./lib/scip/applications/Coloring/src/pricer_coloring.d \
./lib/scip/applications/Coloring/src/probdata_coloring.d \
./lib/scip/applications/Coloring/src/reader_col.d \
./lib/scip/applications/Coloring/src/reader_csol.d 


# Each subdirectory must supply rules for building sources it contributes
lib/scip/applications/Coloring/src/%.o: ../lib/scip/applications/Coloring/src/%.c
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C Compiler'
	gcc -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


