#!/bin/bash

#https://superuser.com/questions/135117/how-to-convert-video-to-images
# if the command would create more than 10,000 frames change the %04d to %05d


ffmpeg -i IMG_1810.mp4 output_%04d.jpg