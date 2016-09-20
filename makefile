
RM := rm -rf


CPP_SRCS += \
src/KmerLight.cpp 

OBJS += \
./KmerLight.o 



# Each subdirectory must supply rules for building sources it contributes
%.o: ./src/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -O3 -g3 -Wall -c -fmessage-length=0 -std=c++11 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '




LIBS := -lpthread -lz


# Add inputs and outputs from these tool invocations to the build variables 

# All Target
all: kmerlight

# Tool invocations
kmerlight: $(OBJS) $(USER_OBJS)
	@echo 'Building target: $@'
	@echo 'Invoking: GCC C++ Linker'
	g++ -L./libs -o "kmerlight" $(OBJS) $(USER_OBJS) $(LIBS)
	@echo 'Finished building target: $@'
	@echo ' '

# Other Targets
clean:
	-$(RM) $(CC_DEPS)$(C++_DEPS)$(EXECUTABLES)$(C_UPPER_DEPS)$(CXX_DEPS)$(OBJS)$(CPP_DEPS)$(C_DEPS) kmerlight
	-@echo ' '

.PHONY: all clean dependents
.SECONDARY:
