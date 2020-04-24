#!/bin/bash

meson setup --wipe --buildtype=release build
meson install -C build
