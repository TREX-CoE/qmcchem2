#!/bin/bash

export srcdir="."
git submodule init
git submodule update
autoreconf -i -Wall --no-recursive
