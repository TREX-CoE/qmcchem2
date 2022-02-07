#!/bin/bash

export srcdir="."
git submodule init
git submodule update
autoreconf -f -i -Wall --no-recursive
