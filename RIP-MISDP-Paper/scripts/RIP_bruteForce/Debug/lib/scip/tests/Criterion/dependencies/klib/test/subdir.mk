################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
C_SRCS += \
../lib/scip/tests/Criterion/dependencies/klib/test/kbit_test.c \
../lib/scip/tests/Criterion/dependencies/klib/test/kbtree_test.c \
../lib/scip/tests/Criterion/dependencies/klib/test/kgraph_test.c \
../lib/scip/tests/Criterion/dependencies/klib/test/khash_keith.c \
../lib/scip/tests/Criterion/dependencies/klib/test/khash_keith2.c \
../lib/scip/tests/Criterion/dependencies/klib/test/khash_test.c \
../lib/scip/tests/Criterion/dependencies/klib/test/klist_test.c \
../lib/scip/tests/Criterion/dependencies/klib/test/kmin_test.c \
../lib/scip/tests/Criterion/dependencies/klib/test/kseq_bench.c \
../lib/scip/tests/Criterion/dependencies/klib/test/kseq_bench2.c \
../lib/scip/tests/Criterion/dependencies/klib/test/kseq_test.c \
../lib/scip/tests/Criterion/dependencies/klib/test/ksort_test.c \
../lib/scip/tests/Criterion/dependencies/klib/test/kstring_bench.c \
../lib/scip/tests/Criterion/dependencies/klib/test/kstring_bench2.c \
../lib/scip/tests/Criterion/dependencies/klib/test/kstring_test.c \
../lib/scip/tests/Criterion/dependencies/klib/test/kthread_test.c \
../lib/scip/tests/Criterion/dependencies/klib/test/kthread_test2.c 

OBJS += \
./lib/scip/tests/Criterion/dependencies/klib/test/kbit_test.o \
./lib/scip/tests/Criterion/dependencies/klib/test/kbtree_test.o \
./lib/scip/tests/Criterion/dependencies/klib/test/kgraph_test.o \
./lib/scip/tests/Criterion/dependencies/klib/test/khash_keith.o \
./lib/scip/tests/Criterion/dependencies/klib/test/khash_keith2.o \
./lib/scip/tests/Criterion/dependencies/klib/test/khash_test.o \
./lib/scip/tests/Criterion/dependencies/klib/test/klist_test.o \
./lib/scip/tests/Criterion/dependencies/klib/test/kmin_test.o \
./lib/scip/tests/Criterion/dependencies/klib/test/kseq_bench.o \
./lib/scip/tests/Criterion/dependencies/klib/test/kseq_bench2.o \
./lib/scip/tests/Criterion/dependencies/klib/test/kseq_test.o \
./lib/scip/tests/Criterion/dependencies/klib/test/ksort_test.o \
./lib/scip/tests/Criterion/dependencies/klib/test/kstring_bench.o \
./lib/scip/tests/Criterion/dependencies/klib/test/kstring_bench2.o \
./lib/scip/tests/Criterion/dependencies/klib/test/kstring_test.o \
./lib/scip/tests/Criterion/dependencies/klib/test/kthread_test.o \
./lib/scip/tests/Criterion/dependencies/klib/test/kthread_test2.o 

C_DEPS += \
./lib/scip/tests/Criterion/dependencies/klib/test/kbit_test.d \
./lib/scip/tests/Criterion/dependencies/klib/test/kbtree_test.d \
./lib/scip/tests/Criterion/dependencies/klib/test/kgraph_test.d \
./lib/scip/tests/Criterion/dependencies/klib/test/khash_keith.d \
./lib/scip/tests/Criterion/dependencies/klib/test/khash_keith2.d \
./lib/scip/tests/Criterion/dependencies/klib/test/khash_test.d \
./lib/scip/tests/Criterion/dependencies/klib/test/klist_test.d \
./lib/scip/tests/Criterion/dependencies/klib/test/kmin_test.d \
./lib/scip/tests/Criterion/dependencies/klib/test/kseq_bench.d \
./lib/scip/tests/Criterion/dependencies/klib/test/kseq_bench2.d \
./lib/scip/tests/Criterion/dependencies/klib/test/kseq_test.d \
./lib/scip/tests/Criterion/dependencies/klib/test/ksort_test.d \
./lib/scip/tests/Criterion/dependencies/klib/test/kstring_bench.d \
./lib/scip/tests/Criterion/dependencies/klib/test/kstring_bench2.d \
./lib/scip/tests/Criterion/dependencies/klib/test/kstring_test.d \
./lib/scip/tests/Criterion/dependencies/klib/test/kthread_test.d \
./lib/scip/tests/Criterion/dependencies/klib/test/kthread_test2.d 


# Each subdirectory must supply rules for building sources it contributes
lib/scip/tests/Criterion/dependencies/klib/test/%.o: ../lib/scip/tests/Criterion/dependencies/klib/test/%.c
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C Compiler'
	gcc -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


