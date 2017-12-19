SUBDIRS = burgers stvenant godunov

.PHONY: all clean

all clean:
	@for dir in $(SUBDIRS); do \
	$(MAKE) -C $$dir -f makefile $@; \
	done
