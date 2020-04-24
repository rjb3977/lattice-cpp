#!/bin/bash

if [[ ! -d build ]]; then
    meson setup --wipe --buildtype=release build
fi

meson install -C build
