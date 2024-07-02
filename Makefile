
ifdef DEBUG
CHPL_FLAGS=-g
else
CHPL_FLAGS=--fast
endif

all: pi_cpu pi_gpu heat_cpu heat_gpu

pi_cpu: pi.chpl
	CHPL_LOCALE_MODEL=flat chpl $(CHPL_FLAGS) $< -o $@

pi_gpu: pi.chpl
	CHPL_LOCALE_MODEL=gpu CHPL_GPU=nvidia chpl $(CHPL_FLAGS) $< -o $@

heat_cpu: heat.chpl
	CHPL_LOCALE_MODEL=flat chpl $(CHPL_FLAGS) $< -o $@

heat_gpu: heat.chpl
	CHPL_LOCALE_MODEL=gpu CHPL_GPU=nvidia chpl $(CHPL_FLAGS) $< -o $@

.PHONY: clean

clean:
	rm pi_cpu pi_gpu heat_cpu heat_gpu
