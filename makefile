SUBDIRS = burgers stvenant

.PHONY: all clean

all clean:
	for dir in $(SUBDIRS); do \
	$(MAKE) -C $$dir -f Makefile $@; \
	done
