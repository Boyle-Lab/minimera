.PHONY: all

all: build/minimera build/minimera.1

# Build -----------------------------------------------------------------------
lisps := $(shell ffind '\.(asd|lisp)$$')

build/minimera: $(lisps) Makefile
	mkdir -p build/
	sbcl --disable-debugger --load "src/build-binary.lisp"

build/minimera.1: $(lisps) Makefile
	mkdir -p build/
	sbcl --disable-debugger --load "src/build-manual.lisp" --quit

clean:
	rm build/minimera*
