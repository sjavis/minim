BENCHMARKS = pf-interface pf-volume pf-pillar-collapse

.PHONY: all clean $(BENCHMARKS)

# Use 'failed' file to determine success so that all tests are performed even if one fails
all: $(BENCHMARKS) clean
	@[ -f failed ] && rm failed && exit 1 || exit 0

$(BENCHMARKS): %:
	@echo "Running $*"
	@$(MAKE) -C $* > /dev/null && (echo "\033[0;32mSUCCESS\033[0m") || (echo "\033[0;31mFAILED\033[0m" && touch failed)

clean:
	-@for dir in $(BENCHMARKS); do $(MAKE) -C $$dir clean > /dev/null; done
