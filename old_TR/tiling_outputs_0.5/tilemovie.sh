#!/bin/bash
rm video_concat_list.txt
ffmpeg -framerate 2 -i tiled_domain_%04d.png -vf scale=800:-1 -c:v libx264 -r 30 -pix_fmt yuv420p tiled_domain.mp4
echo "file 'tiled_domain.mp4'" >> video_concat_list.txt
ffmpeg -framerate 2 -pattern_type glob -i 'tiled_domain_ncd*.png' -vf scale=800:-1 -c:v libx264 -r 30 -pix_fmt yuv420p virtual_tiled_domain.mp4
echo "file 'virtual_tiled_domain.mp4'" >> video_concat_list.txt
ffmpeg -f concat -i video_concat_list.txt -c copy full-virtual_tiled_domain.mp4
rm video_concat_list.txt
