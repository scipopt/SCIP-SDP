################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
C_SRCS += \
../lib/scip/tests/Criterion/dependencies/klib/bgzf.c \
../lib/scip/tests/Criterion/dependencies/klib/kexpr.c \
../lib/scip/tests/Criterion/dependencies/klib/khmm.c \
../lib/scip/tests/Criterion/dependencies/klib/kmath.c \
../lib/scip/tests/Criterion/dependencies/klib/knetfile.c \
../lib/scip/tests/Criterion/dependencies/klib/knhx.c \
../lib/scip/tests/Criterion/dependencies/klib/kopen.c \
../lib/scip/tests/Criterion/dependencies/klib/ksa.c \
../lib/scip/tests/Criterion/dependencies/klib/kson.c \
../lib/scip/tests/Criterion/dependencies/klib/kstring.c \
../lib/scip/tests/Criterion/dependencies/klib/ksw.c \
../lib/scip/tests/Criterion/dependencies/klib/kthread.c \
../lib/scip/tests/Criterion/dependencies/klib/kurl.c 

OBJS += \
./lib/scip/tests/Criterion/dependencies/klib/bgzf.o \
./lib/scip/tests/Criterion/dependencies/klib/kexpr.o \
./lib/scip/tests/Criterion/dependencies/klib/khmm.o \
./lib/scip/tests/Criterion/dependencies/klib/kmath.o \
./lib/scip/tests/Criterion/dependencies/klib/knetfile.o \
./lib/scip/tests/Criterion/dependencies/klib/knhx.o \
./lib/scip/tests/Criterion/dependencies/klib/kopen.o \
./lib/scip/tests/Criterion/dependencies/klib/ksa.o \
./lib/scip/tests/Criterion/dependencies/klib/kson.o \
./lib/scip/tests/Criterion/dependencies/klib/kstring.o \
./lib/scip/tests/Criterion/dependencies/klib/ksw.o \
./lib/scip/tests/Criterion/dependencies/klib/kthread.o \
./lib/scip/tests/Criterion/dependencies/klib/kurl.o 

C_DEPS += \
./lib/scip/tests/Criterion/dependencies/klib/bgzf.d \
./lib/scip/tests/Criterion/dependencies/klib/kexpr.d \
./lib/scip/tests/Criterion/dependencies/klib/khmm.d \
./lib/scip/tests/Criterion/dependencies/klib/kmath.d \
./lib/scip/tests/Criterion/dependencies/klib/knetfile.d \
./lib/scip/tests/Criterion/dependencies/klib/knhx.d \
./lib/scip/tests/Criterion/dependencies/klib/kopen.d \
./lib/scip/tests/Criterion/dependencies/klib/ksa.d \
./lib/scip/tests/Criterion/dependencies/klib/kson.d \
./lib/scip/tests/Criterion/dependencies/klib/kstring.d \
./lib/scip/tests/Criterion/dependencies/klib/ksw.d \
./lib/scip/tests/Criterion/dependencies/klib/kthread.d \
./lib/scip/tests/Criterion/dependencies/klib/kurl.d 


# Each subdirectory must supply rules for building sources it contributes
lib/scip/tests/Criterion/dependencies/klib/%.o: ../lib/scip/tests/Criterion/dependencies/klib/%.c
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C Compiler'
	gcc -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


