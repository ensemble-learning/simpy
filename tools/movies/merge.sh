#mencoder -really-quiet -ovc lavc -lavcopts vcodec=mjpeg -mf fps=${FPS} -vf scale=${videoX}:${videoY} -o $output_video_file_name video_*.avi

FPS=24
videoX=640
videoY=480
output_video_file_name="ml.mp4"

mencoder -really-quiet -ovc lavc -lavcopts vcodec=mjpeg -mf fps=${FPS} -vf scale=${videoX}:${videoY} -o $output_video_file_name 0*.mp4 1*.mp4

mencoder -really-quiet -ovc lavc -lavcopts vcodec=mjpeg -mf fps=24 -vf scale=640:480 -o ml.mp4 0*.mp4 1*.mp4
