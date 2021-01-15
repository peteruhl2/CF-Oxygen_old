#!/bin/bash
# Make a gif from a bunch of .png files in the current folder

convert -delay 100 -loop 0 *.png thing.gif
