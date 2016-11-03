#!/bin/bash

ffmpeg -framerate 1 -i tiled_domain_%04d.png -vf scale=800:-1 -c:v libx264 -r 30 -pix_fmt yuv420p tiled_domain.mp4
