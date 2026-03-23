# Minimera

```
┏┳┓ ┳ ┏┓╻ ┳
┃┃┃ ┃ ┃┗┫ ┃
╹ ╹ ┻ ╹ ╹ ┻
┏┳┓┏━╸┏━┓┏━┓
┃┃┃┣╸ ┣┳┛┣━┫
╹ ╹┗━╸╹┗╸╹ ╹
```

* **License:** [GPLv3 or later](https://www.gnu.org/licenses/gpl-3.0.html)
* **Git:** <https://github.com/Boyle-Lab/minimera/>

Check the [man page][] for information about the command-line `minimera`
interface.

[man page]: https://github.com/Boyle-Lab/minimera/blob/main/DOCUMENTATION.markdown

## Installation

We provide [releases][] of Minimera as binaries (Linux AMD64 only) and
Singularity containers.

[releases]: https://github.com/Boyle-Lab/minimera/releases

## Building

If you just want to *use* Minimera you don't need to build it from scratch.  You
can grab the latest binary or Singularity container from the releases page and
don't need to worry about any of this.

If you want to build Minimera from scratch (e.g. to modify it), here are some
rough notes to get you started.  You should be at least a little familiar with
Common Lisp or it's going to be pretty confusing, sorry.  You'll need at least:

* SBCL
* Quicklisp
* Buildapp

You'll also need to clone a few projects that aren't in Quicklisp into your
Quickload `local-projects` directory (I'll try to eventually get them into
Quicklisp some day, sorry for the fiddliness in the mean time):

* cl-losh
* conserve
* faster

Make sure you can `(ql:quickload :minimera)` successfully.  Then you should be
able to run `make` to:

1. Generate the ASDF manifest.
2. Build a `minimera` binary without Quicklisp by using `buildapp` with that manifest.
3. Build a `minimera.1` man page.

## Tests

There's a (very small) [Cram](https://bitheap.org/cram/) test suite, mostly so
I can to avoid breaking things inadvertently before a release.  Once you've got
cram installed `make test` to run it.
