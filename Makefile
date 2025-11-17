.PHONY: all release

all: build/minimera build/minimera.1

# Build -----------------------------------------------------------------------
lisps := $(shell ffind '\.(asd|lisp)$$')

build/asdf-manifest: Makefile minimera.asd
	mkdir -p build/
	sbcl --disable-debugger --quit --eval '(ql:write-asdf-manifest-file "build/asdf-manifest")'

build/minimera: $(lisps) Makefile build/asdf-manifest
	mkdir -p build/

	buildapp \
		--load-system 'minimera' \
		--eval '(setf minimera::*version* "'"$(shell git describe --dirty)"'")' \
		--entry 'minimera:toplevel' \
		--manifest-file 'build/asdf-manifest' \
		--compress-core \
		--output 'build/minimera'

build/minimera.1: $(lisps) Makefile
	mkdir -p build/
	sbcl --disable-debugger --load "build-manual.lisp" --quit

build/minimera.sif: build/minimera build/minimera.1 contrib/minimera.def
	singularity build --fakeroot build/minimera.sif contrib/minimera.def


release: build/minimera build/minimera.1 build/minimera.sif build-release.sh Makefile
	./build-release.sh

clean:
	rm -r ./build/
