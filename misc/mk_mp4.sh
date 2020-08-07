#!/bin/bash

## Make video with ffmpeg
# First make a list of the files and put in a file

rm joint.txt

im="../_plots/div_movie_frames_1h/*$1*"

for if in $(ls $im)
do	
	echo "file '$if'" >> joint.txt
done 

# Make video from the files in the list
#ffmpeg -r 04 -f concat -i joint.txt -c:v libx264 sss_96.mp4
#ffmpeg -r 08 -f concat -i joint.txt -c:v libx264 -vf scale=1080:-2 -pix_fmt yuv420p test.mp4 
ffmpeg -r 06 -f concat -safe 0 -i joint.txt -c:v libx264 -vf scale=2160:-2 -pix_fmt yuv420p $1.mp4
