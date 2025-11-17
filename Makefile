.PHONY: all

all: build/minimera build/minimera.1

# Build -----------------------------------------------------------------------
lisps := $(shell ffind '\.(asd|lisp)$$')

build/asdf-manifest: Makefile minimera.asd
	mkdir -p build/
	sbcl --disable-debugger --quit --eval '(ql:write-asdf-manifest-file "build/asdf-manifest")'

# build/minimera: $(lisps) Makefile
# 	mkdir -p build/
# 	sbcl --disable-debugger --load "src/build-binary.lisp"

build/minimera: $(lisps) Makefile build/asdf-manifest
	mkdir -p build/
	buildapp \
		--load-system 'minimera' \
		--entry 'minimera:toplevel' \
		--manifest-file 'build/asdf-manifest' \
		--compress-core \
		--output 'build/minimera'

build/minimera.1: $(lisps) Makefile
	mkdir -p build/
	sbcl --disable-debugger --load "src/build-manual.lisp" --quit

clean:
	rm build/minimera*
	rm build/asdf-manifest
