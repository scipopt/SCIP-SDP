################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
C_SRCS += \
../lib/scip/applications/STP/src/branch_stp.c \
../lib/scip/applications/STP/src/cmain.c \
../lib/scip/applications/STP/src/cons_stp.c \
../lib/scip/applications/STP/src/dialog_stp.c \
../lib/scip/applications/STP/src/event_bestsol.c \
../lib/scip/applications/STP/src/grphbase.c \
../lib/scip/applications/STP/src/grphload.c \
../lib/scip/applications/STP/src/grphmcut.c \
../lib/scip/applications/STP/src/grphpath.c \
../lib/scip/applications/STP/src/grphsave.c \
../lib/scip/applications/STP/src/heur_ascendprune.c \
../lib/scip/applications/STP/src/heur_local.c \
../lib/scip/applications/STP/src/heur_prune.c \
../lib/scip/applications/STP/src/heur_rec.c \
../lib/scip/applications/STP/src/heur_slackprune.c \
../lib/scip/applications/STP/src/heur_tm.c \
../lib/scip/applications/STP/src/misc_stp.c \
../lib/scip/applications/STP/src/pricer_stp.c \
../lib/scip/applications/STP/src/probdata_stp.c \
../lib/scip/applications/STP/src/prop_stp.c \
../lib/scip/applications/STP/src/reader_stp.c \
../lib/scip/applications/STP/src/reduce.c \
../lib/scip/applications/STP/src/reduce_alt.c \
../lib/scip/applications/STP/src/reduce_bnd.c \
../lib/scip/applications/STP/src/reduce_simple.c \
../lib/scip/applications/STP/src/validate.c 

OBJS += \
./lib/scip/applications/STP/src/branch_stp.o \
./lib/scip/applications/STP/src/cmain.o \
./lib/scip/applications/STP/src/cons_stp.o \
./lib/scip/applications/STP/src/dialog_stp.o \
./lib/scip/applications/STP/src/event_bestsol.o \
./lib/scip/applications/STP/src/grphbase.o \
./lib/scip/applications/STP/src/grphload.o \
./lib/scip/applications/STP/src/grphmcut.o \
./lib/scip/applications/STP/src/grphpath.o \
./lib/scip/applications/STP/src/grphsave.o \
./lib/scip/applications/STP/src/heur_ascendprune.o \
./lib/scip/applications/STP/src/heur_local.o \
./lib/scip/applications/STP/src/heur_prune.o \
./lib/scip/applications/STP/src/heur_rec.o \
./lib/scip/applications/STP/src/heur_slackprune.o \
./lib/scip/applications/STP/src/heur_tm.o \
./lib/scip/applications/STP/src/misc_stp.o \
./lib/scip/applications/STP/src/pricer_stp.o \
./lib/scip/applications/STP/src/probdata_stp.o \
./lib/scip/applications/STP/src/prop_stp.o \
./lib/scip/applications/STP/src/reader_stp.o \
./lib/scip/applications/STP/src/reduce.o \
./lib/scip/applications/STP/src/reduce_alt.o \
./lib/scip/applications/STP/src/reduce_bnd.o \
./lib/scip/applications/STP/src/reduce_simple.o \
./lib/scip/applications/STP/src/validate.o 

C_DEPS += \
./lib/scip/applications/STP/src/branch_stp.d \
./lib/scip/applications/STP/src/cmain.d \
./lib/scip/applications/STP/src/cons_stp.d \
./lib/scip/applications/STP/src/dialog_stp.d \
./lib/scip/applications/STP/src/event_bestsol.d \
./lib/scip/applications/STP/src/grphbase.d \
./lib/scip/applications/STP/src/grphload.d \
./lib/scip/applications/STP/src/grphmcut.d \
./lib/scip/applications/STP/src/grphpath.d \
./lib/scip/applications/STP/src/grphsave.d \
./lib/scip/applications/STP/src/heur_ascendprune.d \
./lib/scip/applications/STP/src/heur_local.d \
./lib/scip/applications/STP/src/heur_prune.d \
./lib/scip/applications/STP/src/heur_rec.d \
./lib/scip/applications/STP/src/heur_slackprune.d \
./lib/scip/applications/STP/src/heur_tm.d \
./lib/scip/applications/STP/src/misc_stp.d \
./lib/scip/applications/STP/src/pricer_stp.d \
./lib/scip/applications/STP/src/probdata_stp.d \
./lib/scip/applications/STP/src/prop_stp.d \
./lib/scip/applications/STP/src/reader_stp.d \
./lib/scip/applications/STP/src/reduce.d \
./lib/scip/applications/STP/src/reduce_alt.d \
./lib/scip/applications/STP/src/reduce_bnd.d \
./lib/scip/applications/STP/src/reduce_simple.d \
./lib/scip/applications/STP/src/validate.d 


# Each subdirectory must supply rules for building sources it contributes
lib/scip/applications/STP/src/%.o: ../lib/scip/applications/STP/src/%.c
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C Compiler'
	gcc -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


