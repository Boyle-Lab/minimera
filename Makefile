.PHONY: all release binary clean contrib test

binary: build/minimera

all: build/minimera build/minimera.1 build/minimera.sif DOCUMENTATION.markdown

contrib: build/fastq-stats build/minimera.fish

# Config ----------------------------------------------------------------------
lisps := $(shell ffind '\.(asd|lisp)$$')

# Main Build ------------------------------------------------------------------
build/asdf-manifest: Makefile minimera.asd
	mkdir -p build/
	sbcl --disable-debugger --quit --eval '(ql:write-asdf-manifest-file "build/asdf-manifest")'

build/minimera: $(lisps) Makefile build/asdf-manifest build-binary.sh
	mkdir -p build/
	./build-binary.sh

build/minimera.1: $(lisps) Makefile
	mkdir -p build/
	sbcl --disable-debugger --load "build-manual.lisp" --quit

build/minimera.sif: build/minimera build/minimera.1 contrib/minimera.def
	mkdir -p build/
	singularity build --force --fakeroot build/minimera.sif contrib/minimera.def

DOCUMENTATION.markdown: Makefile src/ui.lisp build-manual.lisp
	sbcl --disable-debugger --load "build-manual.lisp" --quit

# Releases --------------------------------------------------------------------
release: build/minimera build/minimera.1 build/minimera.fish build/minimera.sif build-release.sh Makefile
	mkdir -p build/
	./build-release.sh

# Contrib ---------------------------------------------------------------------
build/minimera.fish: $(lisps) Makefile
	mkdir -p build/
	sbcl --disable-debugger --load "build-manual.lisp" --quit

build/fastq-stats: contrib/fastq-stats.lisp Makefile
	mkdir -p build/
	sbcl --disable-debugger --load "contrib/fastq-stats.lisp" --eval "(fastq-stats:build)" --quit

# Test ------------------------------------------------------------------------
test:
	cd tests && cram minimera.t

# Clean -----------------------------------------------------------------------
clean:
	rm -r ./build/
